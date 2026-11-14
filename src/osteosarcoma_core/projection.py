"""Rank-mixture projection with donor-balanced pure-state constraints."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.linalg import null_space, solve
from scipy.stats import rankdata, spearmanr, wilcoxon

SEED = 20261002
PROGRAMMES = ("Replication", "Redox", "Matrix")
RNA_FRACTIONS = (0.2, 0.4, 0.6, 0.8, 1.0)
BACKGROUNDS = ("Myeloid", "T cell", "B cell", "Fibroblast", "Endothelial", "Osteoclast")
GRID = tuple(itertools.product((0.0, 0.25, 1.0, 4.0), (0.01, 0.1, 1.0)))
METHODS = {
    "CRMP fixed": {"balance": False, "constrained": False, "fixed": (1.0, 0.1)},
    "CRMP nested": {"balance": False, "constrained": False},
    "Balanced CRMP": {"balance": True, "constrained": False},
    "Pure-state CRMP": {"balance": True, "constrained": True},
    "Pure-state no drift": {
        "balance": True,
        "constrained": True,
        "grid": ((0.0, 0.01), (0.0, 0.1), (0.0, 1.0)),
    },
}


def balanced_weights(frame: pd.DataFrame, balance: bool = True) -> np.ndarray:
    """Assign equal total weight to each donor and each available background."""
    if frame.empty:
        raise ValueError("At least one mixture record is required")
    if not {"donor", "reference"}.issubset(frame.columns):
        raise ValueError("Mixture metadata require donor and reference columns")
    if frame[["donor", "reference"]].isna().any().any():
        raise ValueError("Donor and reference identifiers must be complete")
    if not balance:
        return np.repeat(1.0 / len(frame), len(frame))
    n = frame.groupby(["donor", "reference"], observed=True).donor.transform("size").to_numpy()
    k = frame.groupby("donor", observed=True).reference.transform("nunique").to_numpy()
    return 1.0 / (frame.donor.nunique() * k * n)


def fit_projection(
    X: np.ndarray,
    Y: np.ndarray,
    X0: np.ndarray,
    frame: pd.DataFrame,
    W: np.ndarray,
    lam: float = 1.0,
    alpha: float = 0.1,
    balance: bool = False,
    constrained: bool = False,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit anchored CRMP coefficients and an intercept."""
    X, Y, X0, W = (np.asarray(a, dtype=float) for a in (X, Y, X0, W))
    if X.ndim != 2 or Y.ndim != 2 or X0.shape != X.shape:
        raise ValueError("X and X0 must share a matrix shape; Y must be a matrix")
    if len(frame) != len(X) or len(Y) != len(X) or W.shape != (X.shape[1], Y.shape[1]):
        raise ValueError("Projection arrays and metadata have inconsistent dimensions")
    if X.shape[1] == 0 or Y.shape[1] == 0 or not all(np.isfinite(a).all() for a in (X, Y, X0, W)):
        raise ValueError("Projection arrays must be finite and have nonempty feature and target axes")
    if lam < 0 or alpha <= 0:
        raise ValueError("The drift penalty must be nonnegative and the coefficient penalty positive")
    w = balanced_weights(frame, balance)
    D = X - X0
    if not constrained:
        xm, ym = w @ X, w @ Y
        xc, yc = X - xm, Y - ym
        B = solve(
            xc.T @ (w[:, None] * xc) + lam * D.T @ (w[:, None] * D) + alpha * np.eye(X.shape[1]),
            xc.T @ (w[:, None] * yc) + alpha * W,
            assume_a="pos",
        )
        return B, ym - xm @ B, float("nan")
    first = np.unique(frame.donor.to_numpy(), return_index=True)[1]
    pure_aug = np.column_stack([X0[first], np.ones(len(first))])
    N = null_space(pure_aug, rcond=1e-12)
    baseline = np.vstack([W, np.zeros(W.shape[1])])
    if N.shape[1] == 0:
        return W.copy(), np.zeros(W.shape[1]), 0.0
    Z = np.column_stack([X, np.ones(len(X))]) @ N
    V = np.column_stack([D, np.zeros(len(D))]) @ N
    penalty = N[:-1].T @ N[:-1]
    h = solve(
        Z.T @ (w[:, None] * Z) + lam * V.T @ (w[:, None] * V) + alpha * penalty,
        Z.T @ (w[:, None] * (Y - X @ W)) - lam * V.T @ (w[:, None] * (D @ W)),
        assume_a="pos",
    )
    corrected = baseline + N @ h
    error = float(np.max(np.abs(pure_aug @ corrected - pure_aug @ baseline)))
    if error > 1e-10:
        raise ArithmeticError("Pure-state equality tolerance exceeded")
    return corrected[:-1], corrected[-1], error


def rank_features(
    expression: pd.DataFrame,
    gene_universe: Sequence[str],
    feature_genes: Sequence[str],
) -> np.ndarray:
    """Rank genes within each expression profile and select feature columns."""
    genes, features = list(gene_universe), list(feature_genes)
    if not genes or len(set(genes)) != len(genes) or len(set(features)) != len(features):
        raise ValueError("Gene and feature identifiers must be unique, with a nonempty gene universe")
    if not expression.index.is_unique:
        raise ValueError("Expression gene identifiers must be unique")
    lookup = {gene: i for i, gene in enumerate(genes)}
    if any(gene not in lookup for gene in features):
        raise ValueError("Every feature gene must belong to the gene universe")
    values = expression.reindex(genes).fillna(0).to_numpy(dtype=float).T
    if not np.isfinite(values).all():
        raise ValueError("Expression profiles must contain finite values")
    return (rankdata(values, axis=1, method="average") / len(genes))[:, [lookup[g] for g in features]]


def programme_baseline(
    feature_genes: Sequence[str],
    programme_genes: Mapping[str, Sequence[str]],
    programmes: Sequence[str] = PROGRAMMES,
) -> np.ndarray:
    """Build programme-average coefficients over the feature genes."""
    features = list(feature_genes)
    if len(set(features)) != len(features):
        raise ValueError("Feature identifiers must be unique")
    lookup = {gene: i for i, gene in enumerate(features)}
    W = np.zeros((len(features), len(programmes)))
    for j, programme in enumerate(programmes):
        matched = [gene for gene in programme_genes[programme] if gene in lookup]
        if not matched:
            raise ValueError(f"Programme {programme} has no feature genes")
        W[[lookup[gene] for gene in matched], j] = 1.0 / len(matched)
    return W


def build_mixtures(
    pseudobulks: Mapping[tuple[str, str], pd.DataFrame],
    gene_universe: Sequence[str],
    feature_genes: Sequence[str],
    W: np.ndarray,
    rna_fractions: Sequence[float] = RNA_FRACTIONS,
    backgrounds: Sequence[str] = BACKGROUNDS,
) -> tuple[dict[str, np.ndarray], pd.DataFrame, pd.DataFrame]:
    """Construct within-donor RNA mixtures from gene-by-compartment pseudobulks."""
    genes = list(gene_universe)
    W = np.asarray(W, dtype=float)
    if W.ndim != 2 or W.shape[0] != len(feature_genes):
        raise ValueError("W must have one row per feature gene")
    if not rna_fractions or any(not 0 <= p <= 1 for p in rna_fractions):
        raise ValueError("RNA fractions must be a nonempty sequence between zero and one")
    records, pure_records, xs, ys, x0s, pure_x, pure_y = [], [], [], [], [], [], []
    for (cohort, donor), expression in sorted(pseudobulks.items()):
        if "Osteogenic" not in expression.columns:
            continue
        values = expression.reindex(genes).fillna(0).astype(float)
        if not np.isfinite(values.to_numpy()).all() or (values.to_numpy() < 0).any():
            raise ValueError("Pseudobulk expression must be finite and nonnegative")
        available = [name for name in backgrounds if name in values.columns]
        if not available:
            continue
        selected = values[["Osteogenic", *available]]
        totals = selected.sum(axis=0)
        if (totals <= 0).any():
            raise ValueError("Every selected compartment must have a positive expression total")
        normalized = selected / totals
        pure = normalized[["Osteogenic"]]
        x0 = rank_features(pure, genes, feature_genes)[0]
        y = x0 @ W
        pure_records.append({"cohort": str(cohort), "donor": str(donor)})
        pure_x.append(x0)
        pure_y.append(y)
        for background in available:
            for fraction in rna_fractions:
                mixture = fraction * normalized.Osteogenic + (1.0 - fraction) * normalized[background]
                xs.append(rank_features(mixture.to_frame("mixture"), genes, feature_genes)[0])
                ys.append(y)
                x0s.append(x0)
                records.append({"cohort": str(cohort), "donor": str(donor), "reference": background, "rna_fraction": float(fraction)})
    if not records:
        raise ValueError("No eligible donor/background mixtures were found")
    arrays = {"X": np.asarray(xs), "Y": np.asarray(ys), "X0": np.asarray(x0s), "W": W,
              "PX": np.asarray(pure_x), "PY": np.asarray(pure_y)}
    return arrays, pd.DataFrame(records), pd.DataFrame(pure_records)


def _validate_inputs(arrays, mixtures, pure_profiles, programmes, discovery_cohort):
    required = {"X", "Y", "X0", "W", "PX", "PY"}
    if not required.issubset(arrays):
        raise ValueError(f"Projection arrays require {', '.join(sorted(required))}")
    a = {key: np.asarray(arrays[key], dtype=float) for key in required}
    X, Y, X0, W, PX, PY = (a[key] for key in ("X", "Y", "X0", "W", "PX", "PY"))
    if any(value.ndim != 2 or not np.isfinite(value).all() for value in a.values()):
        raise ValueError("Every projection array must be a finite matrix")
    if X0.shape != X.shape or Y.shape != (len(X), len(programmes)) or W.shape != (X.shape[1], len(programmes)):
        raise ValueError("Mixture arrays have inconsistent dimensions")
    if PX.shape != (len(pure_profiles), X.shape[1]) or PY.shape != (len(pure_profiles), len(programmes)) or len(mixtures) != len(X):
        raise ValueError("Pure-profile arrays and metadata have inconsistent dimensions")
    if not {"cohort", "donor", "reference", "rna_fraction"}.issubset(mixtures.columns):
        raise ValueError("Mixture metadata require cohort, donor, reference and rna_fraction")
    if not {"cohort", "donor"}.issubset(pure_profiles.columns):
        raise ValueError("Pure-profile metadata require cohort and donor")
    meta, pm = mixtures.copy().reset_index(drop=True), pure_profiles.copy().reset_index(drop=True)
    if meta[["cohort", "donor", "reference"]].isna().any().any() or pm[["cohort", "donor"]].isna().any().any():
        raise ValueError("Cohort, donor and background identifiers must be complete")
    for name in ("cohort", "donor"):
        meta[name], pm[name] = meta[name].astype(str), pm[name].astype(str)
    if pm.duplicated(["cohort", "donor"]).any():
        raise ValueError("Supply one pure profile per cohort and donor")
    mixture_units = set(map(tuple, meta[["cohort", "donor"]].to_numpy()))
    if mixture_units != set(map(tuple, pm[["cohort", "donor"]].to_numpy())):
        raise ValueError("Mixture and pure-profile donor sets must agree")
    train = meta.cohort.eq(discovery_cohort).to_numpy()
    if meta.loc[train, "donor"].nunique() < 3 or not (~train).any():
        raise ValueError("Grouped nested validation requires at least three discovery donors and an independent cohort")
    if meta.donor.eq("EXTERNAL").any():
        raise ValueError("EXTERNAL is reserved for the independent evaluation fold")
    return a, meta, pm, train


def run_projection(
    arrays: Mapping[str, np.ndarray],
    mixtures: pd.DataFrame,
    pure_profiles: pd.DataFrame,
    programmes: Sequence[str] = PROGRAMMES,
    discovery_cohort: str = "GSE152048",
    seed: int = SEED,
    grid: Sequence[tuple[float, float]] = GRID,
    bootstrap_resamples: int = 10000,
) -> dict[str, pd.DataFrame | dict]:
    """Run donor-grouped nested selection and independent-cohort evaluation."""
    if bootstrap_resamples < 1 or not grid:
        raise ValueError("A positive bootstrap count and nonempty parameter grid are required")
    programmes = tuple(programmes)
    a, meta, pm, train = _validate_inputs(arrays, mixtures, pure_profiles, programmes, discovery_cohort)
    X, Y, X0, W, PX, PY = (a[key] for key in ("X", "Y", "X0", "W", "PX", "PY"))
    discovery = sorted(meta.loc[train, "donor"].unique())
    external = sorted(meta.loc[~train, "donor"].unique())
    predictions, pure_predictions, selections, surfaces, assignments, coefficients = [], [], [], [], [], []
    cache = {}

    def fitted(donors, method, la, al):
        key = (tuple(sorted(donors)), method, la, al)
        if key not in cache:
            t = train & meta.donor.isin(donors).to_numpy()
            params = METHODS[method]
            cache[key] = fit_projection(X[t], Y[t], X0[t], meta.loc[t], W, la, al, params["balance"], params["constrained"])
        return cache[key]

    def append_predictions(method, test, pred, test_pure, pp):
        for i, idx in enumerate(np.where(test)[0]):
            for j, programme in enumerate(programmes):
                predictions.append({**meta.iloc[idx].to_dict(), "programme": programme, "method": method,
                                    "truth": Y[idx, j], "prediction": pred[i, j]})
        for i, idx in enumerate(np.where(test_pure)[0]):
            for j, programme in enumerate(programmes):
                pure_predictions.append({**pm.iloc[idx].to_dict(), "programme": programme, "method": method,
                                         "truth": PY[idx, j], "prediction": pp[i, j]})

    for heldout in discovery + ["EXTERNAL"]:
        eligible = [donor for donor in discovery if donor != heldout]
        tested = external if heldout == "EXTERNAL" else [heldout]
        test = meta.donor.isin(tested).to_numpy() & (~train if heldout == "EXTERNAL" else train)
        pure_cohort = pm.cohort.ne(discovery_cohort) if heldout == "EXTERNAL" else pm.cohort.eq(discovery_cohort)
        test_pure = pm.donor.isin(tested).to_numpy() & pure_cohort.to_numpy()
        for method, params in METHODS.items():
            if "fixed" in params:
                la, al = params["fixed"]
                selected_loss = float("nan")
            else:
                candidates = params.get("grid", grid)
                losses = []
                for la, al in candidates:
                    donor_losses = []
                    for inner in eligible:
                        inner_train = [donor for donor in eligible if donor != inner]
                        inner_test = train & meta.donor.eq(inner).to_numpy()
                        B, b, _ = fitted(inner_train, method, la, al)
                        row_risk = np.mean((X[inner_test] @ B + b - Y[inner_test]) ** 2, axis=1)
                        donor_losses.append(float(balanced_weights(meta.loc[inner_test], True) @ row_risk))
                    loss = float(np.mean(donor_losses))
                    losses.append(loss)
                    surfaces.append({"outer_fold": heldout, "method": method, "lambda": la, "alpha": al,
                                     "inner_donor_mse": loss, "inner_donors": len(eligible)})
                choice = int(np.argmin(losses))
                la, al = candidates[choice]
                selected_loss = losses[choice]
            B, b, constraint_error = fitted(eligible, method, la, al)
            selections.append({"outer_fold": heldout, "method": method, "lambda": la, "alpha": al,
                               "inner_donor_mse": selected_loss, "training_donors": len(eligible),
                               "test_donors": len(tested), "max_training_pure_error": constraint_error})
            assignments.append({"outer_fold": heldout, "method": method,
                                "training_donors": "|".join(eligible), "heldout_donors": "|".join(tested),
                                "overlap": len(set(eligible) & set(tested)) if heldout != "EXTERNAL" else 0})
            append_predictions(method, test, X[test] @ B + b, test_pure, PX[test_pure] @ B + b)
            if heldout == "EXTERNAL":
                coefficients.append({"method": method, "lambda": float(la), "alpha": float(al),
                                     "coefficients": B.tolist(), "intercept": b.tolist()})
        tr = train & meta.donor.isin(eligible).to_numpy()
        raw, pure_raw = X @ W, PX @ W
        scalar_pred, scalar_pure = np.empty((test.sum(), len(programmes))), np.empty((test_pure.sum(), len(programmes)))
        for j in range(len(programmes)):
            reg = np.linalg.lstsq(np.column_stack([np.ones(tr.sum()), raw[tr, j]]), Y[tr, j], rcond=None)[0]
            scalar_pred[:, j] = reg[0] + reg[1] * raw[test, j]
            scalar_pure[:, j] = reg[0] + reg[1] * pure_raw[test_pure, j]
        append_predictions("Pathway ranks", test, raw[test], test_pure, pure_raw[test_pure])
        append_predictions("Scalar calibration", test, scalar_pred, test_pure, scalar_pure)
        append_predictions("Training mean", test, np.tile(Y[tr].mean(axis=0), (test.sum(), 1)),
                           test_pure, np.tile(Y[tr].mean(axis=0), (test_pure.sum(), 1)))

    pred, pp = pd.DataFrame(predictions), pd.DataFrame(pure_predictions)
    pred["squared_error"] = (pred.prediction - pred.truth) ** 2
    pred["absolute_error"] = np.abs(pred.prediction - pred.truth)
    metrics = pred.groupby(["cohort", "donor", "programme", "method"]).squared_error.mean().pow(0.5).rename("rmse").reset_index()
    background = pred.groupby(["cohort", "donor", "reference", "programme", "method"]).squared_error.mean().pow(0.5).rename("rmse").reset_index()
    purity = pred.groupby(["cohort", "donor", "rna_fraction", "programme", "method"])[["absolute_error", "squared_error"]].mean().reset_index()
    calibration = []
    for (cohort, programme, method), d in pp.groupby(["cohort", "programme", "method"]):
        slope = np.linalg.lstsq(np.column_stack([np.ones(len(d)), d.truth]), d.prediction, rcond=None)[0]
        calibration.append({"cohort": cohort, "programme": programme, "method": method, "donors": len(d),
                            "intercept": slope[0], "slope": slope[1],
                            "sd_ratio": d.prediction.std() / d.truth.std() if d.truth.std() > 0 else float("nan"),
                            "rmse": np.sqrt(np.mean((d.prediction - d.truth) ** 2)),
                            "spearman": spearmanr(d.truth, d.prediction).statistic if d.prediction.std() > 1e-12 and d.truth.std() > 0 else float("nan")})
    rng = np.random.default_rng(seed)
    summaries = []
    for (cohort, programme), d in metrics.groupby(["cohort", "programme"]):
        wide = d.pivot(index="donor", columns="method", values="rmse")
        draws = rng.integers(0, len(wide), (bootstrap_resamples, len(wide)))
        for method in wide:
            raw, z = wide["Pathway ranks"].to_numpy(), wide[method].to_numpy()
            denominator = raw[draws].mean(axis=1)
            boot = np.divide(z[draws].mean(axis=1), denominator, out=np.full(len(draws), np.nan), where=denominator != 0)
            boot = 1.0 - boot
            delta = raw - z
            summaries.append({"cohort": cohort, "programme": programme, "method": method, "donors": len(wide),
                              "mean_rmse": z.mean(), "relative_reduction": 1.0 - z.mean() / raw.mean() if raw.mean() else float("nan"),
                              "ci_low": np.quantile(boot, 0.025), "ci_high": np.quantile(boot, 0.975),
                              "paired_wilcoxon_p": wilcoxon(delta).pvalue if np.any(delta) else 1.0,
                              "donors_improved": int((delta > 0).sum())})
    summary = {"seed": int(seed), "bootstrap_replicates": int(bootstrap_resamples), "features": X.shape[1],
               "mixtures": len(X), "programmes": list(programmes), "discovery_cohort": discovery_cohort,
               "discovery_donors": discovery, "independent_donors": external, "grid": list(grid),
               "fit_count": len(cache), "frozen_models": coefficients}
    return {"heldout_predictions": pred, "donor_metrics": metrics, "background_metrics": background,
            "purity_metrics": purity, "pure_predictions": pp, "pure_calibration": pd.DataFrame(calibration),
            "fold_parameters": pd.DataFrame(selections), "inner_parameter_surface": pd.DataFrame(surfaces),
            "fold_assignments": pd.DataFrame(assignments),
            "mixture_design": meta.groupby(["cohort", "donor", "reference"]).size().rename("mixture_rows").reset_index(),
            "summary": pd.DataFrame(summaries), "run_summary": summary}


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--discovery-cohort", default="GSE152048")
    parser.add_argument("--programmes", nargs="+", default=list(PROGRAMMES))
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--bootstrap-resamples", type=int, default=10000)
    args = parser.parse_args(argv)
    with np.load(args.input_dir / "projection_arrays.npz", allow_pickle=False) as source:
        arrays = {key: source[key] for key in ("X", "Y", "X0", "W", "PX", "PY")}
    mixtures = pd.read_csv(args.input_dir / "mixtures.csv", dtype={"cohort": str, "donor": str, "reference": str})
    pure = pd.read_csv(args.input_dir / "pure_profiles.csv", dtype={"cohort": str, "donor": str})
    outputs = run_projection(arrays, mixtures, pure, args.programmes, args.discovery_cohort, args.seed,
                             bootstrap_resamples=args.bootstrap_resamples)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, value in outputs.items():
        if isinstance(value, pd.DataFrame):
            suffix = ".csv.gz" if name == "heldout_predictions" else ".csv"
            value.to_csv(args.output_dir / (name + suffix), index=False)
        else:
            (args.output_dir / "run_summary.json").write_text(json.dumps(value, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
