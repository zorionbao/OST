"""Ancestry-constrained partial transport and its exact family quotient."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from scipy.spatial.distance import cdist

SEED = 20261002
PROGRAMMES = ("Replication", "Redox", "Matrix")
REJECTION_COSTS = (0.25, 0.5, 0.9, 1.5, 3.0)
DUPLICATE_COPIES = (1, 2, 5, 20, 100, 500, 2000)


def _inputs(cost, families, capacities=None):
    cost = np.asarray(cost, dtype=float)
    families = np.asarray(families)
    if cost.ndim != 2 or min(cost.shape) < 1 or not np.isfinite(cost).all():
        raise ValueError("Cost must be a finite nonempty patient-by-model matrix")
    if families.ndim != 1 or len(families) != cost.shape[1] or pd.isna(families).any():
        raise ValueError("Supply one complete family identifier per model column")
    unique, codes = np.unique(families, return_inverse=True)
    capacity = np.repeat(1.0 / len(unique), len(unique)) if capacities is None else np.asarray(capacities, dtype=float)
    if capacity.shape != (len(unique),) or not np.isfinite(capacity).all() or (capacity < 0).any():
        raise ValueError("Supply a nonnegative capacity for each sorted family identifier")
    return cost, families, unique, codes, capacity


def solve_cost(
    cost: np.ndarray,
    families: Sequence[str],
    kappa: float = 0.9,
    capacities: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Solve partial transport with equal patient mass and shared family capacities."""
    cost, families, unique, codes, capacities = _inputs(cost, families, capacities)
    if not np.isfinite(kappa) or kappa < 0:
        raise ValueError("Unmatched-mass cost must be finite and nonnegative")
    n, m = cost.shape
    nv = n * m + n
    rr, cc = np.repeat(np.arange(n), m), np.arange(n * m)
    eq = sparse.coo_matrix(
        (np.ones(n * m + n), (np.r_[rr, np.arange(n)], np.r_[cc, n * m + np.arange(n)])),
        shape=(n, nv),
    ).tocsr()
    ub = sparse.coo_matrix((np.ones(n * m), (np.tile(codes, n), cc)), shape=(len(unique), nv)).tocsr()
    result = linprog(
        np.r_[cost.ravel(), np.repeat(kappa, n)],
        A_ub=ub,
        b_ub=capacities,
        A_eq=eq,
        b_eq=np.repeat(1.0 / n, n),
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.x[: n * m].reshape(n, m), result.x[n * m :], float(result.fun), unique


def quotient_cost(
    cost: np.ndarray,
    families: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Select the cheapest descendant for each patient and ancestry family."""
    cost, families, unique, _, _ = _inputs(cost, families)
    qcost = np.empty((len(cost), len(unique)))
    chosen = np.empty_like(qcost, dtype=int)
    for j, family in enumerate(unique):
        descendants = np.where(families == family)[0]
        arg = np.argmin(cost[:, descendants], axis=1)
        chosen[:, j] = descendants[arg]
        qcost[:, j] = cost[np.arange(len(cost)), chosen[:, j]]
    return qcost, chosen, unique


def solve_quotient(
    cost: np.ndarray,
    families: Sequence[str],
    kappa: float = 0.9,
    capacities: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, np.ndarray]:
    """Return family coupling, descendant lift, unmatched mass and optimal objective."""
    cost, families, _, _, _ = _inputs(cost, families, capacities)
    qcost, chosen, unique = quotient_cost(cost, families)
    coupling, unmatched, objective, _ = solve_cost(qcost, unique, kappa, capacities)
    lifted = np.zeros_like(cost)
    for j in range(len(unique)):
        lifted[np.arange(len(cost)), chosen[:, j]] += coupling[:, j]
    return coupling, lifted, unmatched, objective, qcost, unique


def patient_model_cost(
    patient_scores: pd.DataFrame,
    model_scores: pd.DataFrame,
    programmes: Sequence[str] = PROGRAMMES,
) -> np.ndarray:
    """Standardize state scores using patients and compute mean squared distance."""
    if len(patient_scores) < 2 or model_scores.empty or not programmes:
        raise ValueError("At least two patients, one model and one programme are required")
    patient_columns = ["CRMP_" + name if "CRMP_" + name in patient_scores else name for name in programmes]
    model_columns = ["CRMP_" + name if "CRMP_" + name in model_scores else name for name in programmes]
    patient = patient_scores[patient_columns].astype(float).copy()
    models = model_scores[model_columns].astype(float).copy()
    patient.columns = models.columns = list(programmes)
    if not np.isfinite(patient.to_numpy()).all() or not np.isfinite(models.to_numpy()).all():
        raise ValueError("State scores must be finite")
    mean, sd = patient.mean(), patient.std()
    if (sd <= 0).any():
        raise ValueError("Every patient programme requires a positive standard deviation")
    x, y = ((patient - mean) / sd).to_numpy(), ((models - mean) / sd).to_numpy()
    return cdist(x, y, "sqeuclidean") / len(programmes)


def duplication_benchmark(
    cost: np.ndarray,
    families: Sequence[str],
    kappa: float = 0.9,
    capacities: np.ndarray | None = None,
    copies: Sequence[int] = DUPLICATE_COPIES,
    repeats: int = 3,
) -> pd.DataFrame:
    """Measure transport under identical-record expansion of the most allocated model."""
    cost, families, unique, _, capacities = _inputs(cost, families, capacities)
    if not copies or any(count < 1 or int(count) != count for count in copies) or repeats < 1:
        raise ValueError("Duplicate counts and timing repeats must be positive integers")
    full, _, objective, _ = solve_cost(cost, families, kappa, capacities)
    selected = int(np.argmax(full.sum(axis=0)))
    duplicate_family = families[selected]
    records = []
    for count in copies:
        expanded = np.column_stack([cost] + [cost[:, [selected]]] * (int(count) - 1))
        expanded_families = np.r_[families, np.repeat(duplicate_family, int(count) - 1)]
        timings = {"Full sparse": [], "Ancestry quotient": []}
        values = {}
        for _ in range(repeats):
            for method in timings:
                start = time.perf_counter()
                if method == "Full sparse":
                    coupling, unmatched, obj, _ = solve_cost(expanded, expanded_families, kappa, capacities)
                    mass = float(coupling[:, expanded_families == duplicate_family].sum())
                else:
                    family_coupling, coupling, unmatched, obj, _, ordered = solve_quotient(expanded, expanded_families, kappa, capacities)
                    mass = float(family_coupling[:, ordered == duplicate_family].sum())
                timings[method].append(time.perf_counter() - start)
                values[method] = (obj, mass, float(coupling.sum()), float(np.max(np.abs(coupling.sum(axis=1) + unmatched - 1.0 / len(cost)))))
        kernel = np.exp(-expanded)
        denominator = kernel.sum(axis=1, keepdims=True)
        kernel = np.divide(kernel, denominator, out=np.zeros_like(kernel), where=denominator != 0)
        kernel_mass = float(kernel[:, expanded_families == duplicate_family].sum(axis=1).mean())
        for method, timings_for_method in timings.items():
            obj, mass, total, residual = values[method]
            ncols = expanded.shape[1] if method == "Full sparse" else len(unique)
            records.append({"copies": int(count), "models": expanded.shape[1], "families": len(unique),
                            "patients": len(cost), "method": method,
                            "median_seconds": float(np.median(timings_for_method)),
                            "min_seconds": min(timings_for_method), "max_seconds": max(timings_for_method),
                            "lp_variables": len(cost) * (ncols + 1),
                            "lp_nonzero_constraints": 2 * len(cost) * ncols + len(cost),
                            "objective": obj, "objective_difference_from_original": abs(obj - objective),
                            "duplicated_family_mass": mass, "total_allocated_mass": total,
                            "row_mass_residual": residual, "kernel_family_mass": kernel_mass})
    return pd.DataFrame(records)


def randomized_equivalence(seed: int = SEED, cases: int = 100) -> pd.DataFrame:
    """Compare independent full and quotient solutions on general-cost problems."""
    if cases < 1:
        raise ValueError("At least one randomized case is required")
    rng = np.random.default_rng(seed)
    records = []
    for replicate in range(cases):
        n, nf = int(rng.integers(3, 41)), int(rng.integers(2, 16))
        counts = rng.integers(1, 7, nf)
        families = np.repeat(np.arange(nf), counts)
        cost = rng.gamma(shape=2.0, scale=0.5, size=(n, len(families)))
        capacities = rng.dirichlet(np.repeat(2.0, nf))
        kappa = float(rng.uniform(0.1, 2.0))
        _, _, full_objective, _ = solve_cost(cost, families, kappa, capacities)
        family, lifted, unmatched, quotient_objective, _, _ = solve_quotient(cost, families, kappa, capacities)
        records.append({"replicate": replicate, "patients": n, "families": nf, "descendants": len(families),
                        "kappa": kappa, "objective_difference": abs(full_objective - quotient_objective),
                        "lift_objective_difference": abs(float((lifted * cost).sum() + kappa * unmatched.sum()) - quotient_objective),
                        "row_residual": float(np.max(np.abs(lifted.sum(axis=1) + unmatched - 1.0 / n))),
                        "family_capacity_violation": float(np.maximum(family.sum(axis=0) - capacities, 0).max())})
    results = pd.DataFrame(records)
    if results[["objective_difference", "lift_objective_difference", "row_residual", "family_capacity_violation"]].to_numpy().max() >= 1e-8:
        raise ArithmeticError("Full and quotient solution tolerance exceeded")
    return results


def run_transport(
    cost: np.ndarray,
    families: Sequence[str],
    patient_ids: Sequence[str] | None = None,
    model_ids: Sequence[str] | None = None,
    kappa: float = 0.9,
    capacities: np.ndarray | None = None,
    rejection_costs: Sequence[float] = REJECTION_COSTS,
    benchmark: bool = False,
    copies: Sequence[int] = DUPLICATE_COPIES,
    repeats: int = 3,
    seed: int = SEED,
) -> dict[str, pd.DataFrame | dict]:
    """Compute full allocation, quotient lift, sensitivity and optional benchmarks."""
    cost, families, unique, _, capacities = _inputs(cost, families, capacities)
    n, m = cost.shape
    patient_ids = list(range(n)) if patient_ids is None else list(patient_ids)
    model_ids = list(range(m)) if model_ids is None else list(model_ids)
    if len(patient_ids) != n or len(model_ids) != m:
        raise ValueError("Patient and model identifiers must align with the cost matrix")
    full, unmatched, objective, _ = solve_cost(cost, families, kappa, capacities)
    family, lifted, quotient_unmatched, quotient_objective, qcost, ordered = solve_quotient(cost, families, kappa, capacities)
    if abs(objective - quotient_objective) >= 1e-10:
        raise ArithmeticError("Full and quotient objective tolerance exceeded")
    rows = []
    for penalty in rejection_costs:
        coupling, u, _, _ = solve_cost(cost, families, penalty, capacities)
        rows.append({"rejection_cost": penalty, "represented_mass": float(coupling.sum()),
                     "mean_matching_cost": float((coupling * cost).sum() / max(coupling.sum(), 1e-12)),
                     "families_used": int(sum(coupling[:, families == f].sum() > 1e-8 for f in unique))})
    outputs = {
        "patient_model_cost": pd.DataFrame(cost, index=patient_ids, columns=model_ids),
        "full_patient_model_transport": pd.DataFrame(full, index=patient_ids, columns=model_ids),
        "patient_family_transport": pd.DataFrame(family, index=patient_ids, columns=ordered),
        "lifted_patient_model_transport": pd.DataFrame(lifted, index=patient_ids, columns=model_ids),
        "quotient_cost": pd.DataFrame(qcost, index=patient_ids, columns=ordered),
        "patient_representation": pd.DataFrame({"sample": patient_ids, "represented_fraction": family.sum(axis=1) * n,
                                                "unmatched_fraction": quotient_unmatched * n}),
        "family_capacities": pd.DataFrame({"family": ordered, "capacity": capacities}),
        "transport_sensitivity": pd.DataFrame(rows),
    }
    summary = {"seed": int(seed), "patients": n, "models": m, "families": len(unique),
               "unmatched_mass_cost": float(kappa), "full_objective": objective,
               "quotient_objective": quotient_objective, "represented_mass": float(family.sum()),
               "objective_difference": abs(objective - quotient_objective),
               "lift_objective_difference": abs(float((lifted * cost).sum() + kappa * quotient_unmatched.sum()) - quotient_objective),
               "row_mass_residual": float(np.max(np.abs(lifted.sum(axis=1) + quotient_unmatched - 1.0 / n))),
               "family_capacity_violation": float(np.maximum(family.sum(axis=0) - capacities, 0).max()),
               "benchmark": bool(benchmark)}
    if benchmark:
        outputs["duplication_benchmark"] = duplication_benchmark(cost, families, kappa, capacities, copies, repeats)
        outputs["randomized_equivalence"] = randomized_equivalence(seed)
        summary["benchmark_repeats"] = int(repeats)
        summary["duplicate_copies"] = list(copies)
        summary["randomized_cases"] = len(outputs["randomized_equivalence"])
    outputs["run_summary"] = summary
    return outputs


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--programmes", nargs="+", default=list(PROGRAMMES))
    parser.add_argument("--kappa", type=float, default=0.9)
    parser.add_argument("--rejection-costs", type=float, nargs="+", default=list(REJECTION_COSTS))
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--copies", type=int, nargs="+", default=list(DUPLICATE_COPIES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args(argv)
    cost_path = args.input_dir / "cost.npz"
    if cost_path.exists():
        with np.load(cost_path, allow_pickle=False) as arrays:
            cost = arrays["cost"]
            patient_ids = arrays["patient_ids"].astype(str).tolist() if "patient_ids" in arrays else None
        models = pd.read_csv(args.input_dir / "families.csv", dtype=str)
        if "family" not in models:
            raise ValueError("families.csv requires a family column")
        families = models.family.to_numpy()
        model_ids = models.model_id.tolist() if "model_id" in models else None
    else:
        patients = pd.read_csv(args.input_dir / "patient_scores.csv", index_col=0)
        models = pd.read_csv(args.input_dir / "model_scores.csv", index_col=0)
        if "family" not in models:
            raise ValueError("model_scores.csv requires a family column")
        cost = patient_model_cost(patients, models, args.programmes)
        families = models.family.to_numpy()
        patient_ids, model_ids = patients.index.tolist(), models.index.tolist()
    capacity_path = args.input_dir / "capacities.csv"
    capacities = None
    if capacity_path.exists():
        capacity = pd.read_csv(capacity_path, dtype={"family": str}).set_index("family")
        capacities = capacity.loc[np.unique(families).astype(str), "capacity"].to_numpy(dtype=float)
    outputs = run_transport(cost, families, patient_ids, model_ids, args.kappa, capacities,
                            args.rejection_costs, args.benchmark, args.copies, args.repeats, args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    matrices = {"patient_model_cost", "full_patient_model_transport", "patient_family_transport",
                "lifted_patient_model_transport", "quotient_cost"}
    for name, value in outputs.items():
        if isinstance(value, pd.DataFrame):
            value.to_csv(args.output_dir / (name + ".csv"), index=name in matrices,
                         index_label="sample" if name in matrices else None)
        else:
            (args.output_dir / "run_summary.json").write_text(json.dumps(value, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
