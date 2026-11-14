"""Paired Context Harmonic Transfer and whole-model nested validation."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import eigh
from scipy.spatial.distance import cdist
from scipy.stats import rankdata, wilcoxon
from statsmodels.stats.multitest import multipletests


LAMBDA_GRID = (0.0, 0.25, 1.0, 4.0)
ETA_GRID = (0.0, 1.0, 4.0)
ALPHA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
GRID = tuple(itertools.product(LAMBDA_GRID, ETA_GRID, ALPHA_GRID))
REFERENCE_METHODS = ("Global anchor", "Paired midpoint")
ESTIMATORS = (
    "No change", "Target mean", "Paired mean", "Unpaired mean",
    "Nearest pair", "Global anchor", "Paired midpoint", "PCHT",
)
ABLATIONS = ("PCHT without graph", "PCHT without context")


@dataclass(frozen=True)
class ContextBenchmark:
    """Ranked screens, paired conditions and functional-class graph."""

    ranks: pd.DataFrame
    metadata: pd.DataFrame
    library: pd.DataFrame
    pairs: pd.DataFrame
    genes: tuple[str, ...]
    classes: tuple[str, ...]
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    source: np.ndarray
    target: np.ndarray
    shifts: np.ndarray


def _boolean_column(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.strip().str.lower()
    if not values.isin(["true", "false", "1", "0"]).all():
        raise ValueError("Library membership values must be true/false or 1/0.")
    return values.isin(["true", "1"])


def build_benchmark(
    effects: pd.DataFrame, metadata: pd.DataFrame, library: pd.DataFrame
) -> ContextBenchmark:
    """Construct paired ranks, weighting every biological model equally."""
    effects = effects.copy()
    metadata = metadata.copy()
    library = library.copy()
    required = {"model", "growth_format", "serum_status"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Screen metadata requires {sorted(required)}.")
    if "all_test_genes" not in library:
        raise ValueError("The library requires all_test_genes.")
    for name, table in [("effects", effects), ("metadata", metadata), ("library", library)]:
        if not table.index.is_unique or table.index.hasnans:
            raise ValueError(f"The {name} index must contain unique nonmissing identifiers.")
    if not set(metadata.index).issubset(effects.index):
        raise ValueError("Every metadata screen must occur in the effects matrix.")
    if metadata[list(required)].isna().any().any():
        raise ValueError("Model and condition identifiers must be complete.")
    values = effects.to_numpy(dtype=float)
    if values.shape[1] == 0 or not np.isfinite(values).all():
        raise ValueError("Effects must have finite values and at least one gene.")
    ranks = pd.DataFrame(rankdata(values, axis=1) / values.shape[1],
                         index=effects.index, columns=effects.columns)
    selected = library.index[_boolean_column(library["all_test_genes"])]
    genes = tuple(selected.intersection(effects.columns, sort=False))
    if not genes:
        raise ValueError("The test library and effect matrix have no shared genes.")
    excluded_classes = {"all_test_genes", "Other", "Actin Nucleation", "Gap Junction"}
    classes = tuple(c for c in library.columns if c not in excluded_classes)
    if not classes:
        raise ValueError("The library requires functional-class membership columns.")
    z = np.column_stack([_boolean_column(library.loc[list(genes), c]).to_numpy(float)
                         for c in classes])
    adjacency = (z / np.maximum(z.sum(axis=0), 1)) @ z.T
    np.fill_diagonal(adjacency, 0)
    laplacian = np.diag(adjacency.sum(axis=1)) - adjacency
    eigenvalues, eigenvectors = eigh(laplacian)
    records = []
    for model, screens in metadata.groupby("model", sort=True):
        for source, target in itertools.combinations(screens.index, 2):
            a, b = screens.loc[source], screens.loc[target]
            if a.serum_status == b.serum_status and a.growth_format != b.growth_format:
                if {a.growth_format, b.growth_format} != {"Plastic", "Dome"}:
                    raise ValueError("Format contrasts require Plastic and Dome labels.")
                source, target = ((source, target) if a.growth_format == "Plastic"
                                  else (target, source))
                contrast = "Format"
            elif a.growth_format == b.growth_format and a.serum_status != b.serum_status:
                if {a.serum_status, b.serum_status} != {"Serum", "Serum-free"}:
                    raise ValueError("Medium contrasts require Serum and Serum-free labels.")
                source, target = ((source, target) if a.serum_status == "Serum"
                                  else (target, source))
                contrast = "Medium"
            else:
                continue
            records.append({"model": model, "source": source, "target": target,
                            "contrast": contrast,
                            "source_label": metadata.loc[source, "growth_format"] + " / "
                                            + metadata.loc[source, "serum_status"],
                            "target_label": metadata.loc[target, "growth_format"] + " / "
                                            + metadata.loc[target, "serum_status"]})
    pairs = pd.DataFrame(records)
    if pairs.empty or pairs.model.nunique() < 3:
        raise ValueError("Nested validation requires paired contrasts from at least three models.")
    pairs.index.name = "pair_id"
    source = ranks.loc[pairs.source, list(genes)].to_numpy()
    target = ranks.loc[pairs.target, list(genes)].to_numpy()
    return ContextBenchmark(ranks, metadata, library, pairs, genes, classes,
                            eigenvalues, eigenvectors, source, target, target - source)


def model_splits(pairs: pd.DataFrame, indices=None):
    """Yield train/test pair indices with all conditions of one model removed."""
    indices = pairs.index.to_numpy() if indices is None else np.asarray(indices, dtype=int)
    for model in sorted(pairs.loc[indices, "model"].unique()):
        heldout = pairs.loc[indices, "model"].eq(model).to_numpy()
        yield model, indices[~heldout], indices[heldout]


def _contrast_indices(data: ContextBenchmark, indices, contrast: str) -> np.ndarray:
    indices = np.asarray(indices, dtype=int)
    result = indices[data.pairs.loc[indices, "contrast"].eq(contrast).to_numpy()]
    if result.size == 0:
        raise ValueError(f"No training pairs are available for {contrast}.")
    return result


def _model_mean(data: ContextBenchmark, values: np.ndarray, indices) -> np.ndarray:
    return pd.DataFrame(values[indices], index=data.pairs.loc[indices, "model"]).groupby(
        level=0, sort=True).mean().mean().to_numpy()


def fit_shift(data: ContextBenchmark, indices, contrast: str,
              shrinkage: float, graph_penalty: float) -> np.ndarray:
    """Solve the graph-regularized, model-balanced paired-shift objective."""
    indices = _contrast_indices(data, indices, contrast)
    mean = _model_mean(data, data.shifts, indices)
    return data.eigenvectors @ ((data.eigenvectors.T @ mean) /
                               (1 + shrinkage + graph_penalty * data.eigenvalues))


def target_anchor(data: ContextBenchmark, indices, contrast: str) -> np.ndarray:
    return _model_mean(data, data.target, _contrast_indices(data, indices, contrast))


def reference_anchor(data: ContextBenchmark, indices, contrast: str, method: str) -> np.ndarray:
    if method == "Global anchor":
        models = data.pairs.loc[indices, "model"].unique()
        screens = data.metadata.index[data.metadata.model.isin(models)]
        return data.ranks.loc[screens, list(data.genes)].set_axis(
            data.metadata.loc[screens, "model"], axis=0).groupby(level=0).mean().mean().to_numpy()
    if method == "Paired midpoint":
        ix = _contrast_indices(data, indices, contrast)
        return _model_mean(data, (data.source + data.target) / 2, ix)
    raise ValueError(f"Unknown reference method: {method}.")


def inner_model_loss(data: ContextBenchmark, indices, parameters, method="PCHT") -> float:
    """Average each inner model's gene-pair squared errors, then average models."""
    model_losses = []
    for _, train, test in model_splits(data.pairs, indices):
        errors = []
        for i in test:
            contrast = data.pairs.loc[i, "contrast"]
            if method != "Global anchor" and not data.pairs.loc[train, "contrast"].eq(contrast).any():
                continue
            if method == "PCHT":
                shrinkage, graph_penalty, alpha = parameters
                shift = fit_shift(data, train, contrast, shrinkage, graph_penalty)
                anchor = target_anchor(data, train, contrast)
                prediction = alpha * (data.source[i] + shift) + (1 - alpha) * anchor
            else:
                alpha = float(parameters)
                prediction = alpha * data.source[i] + (1 - alpha) * reference_anchor(
                    data, train, contrast, method)
            errors.extend((np.clip(prediction, 0, 1) - data.target[i]) ** 2)
        if errors:
            model_losses.append(np.mean(errors))
    if not model_losses:
        raise ValueError("The inner split has no evaluable models.")
    return float(np.mean(model_losses))


def nested_predictions(data: ContextBenchmark) -> dict[str, pd.DataFrame]:
    """Evaluate eight estimators and two component ablations in whole-model folds."""
    predictions, fold_rows, reference_rows, sensitivity = [], [], [], []
    for model, train, test in model_splits(data.pairs):
        losses = [inner_model_loss(data, train, params) for params in GRID]
        selected = int(np.argmin(losses))
        shrinkage, graph_penalty, alpha = GRID[selected]
        common = {"heldout_model": model, "training_models": data.pairs.loc[train, "model"].nunique(),
                  "heldout_pairs": len(test)}
        fold_rows.append({**common, "lambda": shrinkage, "eta": graph_penalty,
                          "alpha": alpha, "inner_mse": losses[selected]})
        selected_alpha = {}
        for method in REFERENCE_METHODS:
            candidate_losses = [inner_model_loss(data, train, aa, method) for aa in ALPHA_GRID]
            choice = int(np.argmin(candidate_losses))
            selected_alpha[method] = ALPHA_GRID[choice]
            reference_rows.append({**common, "method": method, "alpha": ALPHA_GRID[choice],
                                   "inner_mse": candidate_losses[choice]})
        for i in test:
            contrast = data.pairs.loc[i, "contrast"]
            ix = _contrast_indices(data, train, contrast)
            delta = fit_shift(data, train, contrast, shrinkage, graph_penalty)
            mean = fit_shift(data, train, contrast, 0, 0)
            nearest = ix[int(np.argmin(cdist(data.source[[i]], data.source[ix])[0]))]
            screens = data.metadata[data.metadata.model != model]
            column = "growth_format" if contrast == "Format" else "serum_status"
            low, high = ("Plastic", "Dome") if contrast == "Format" else ("Serum", "Serum-free")
            unpaired = (data.ranks.loc[screens.index[screens[column] == high], list(data.genes)].mean().to_numpy()
                        - data.ranks.loc[screens.index[screens[column] == low], list(data.genes)].mean().to_numpy())
            anchor = target_anchor(data, train, contrast)
            methods = {
                "PCHT": alpha * (data.source[i] + delta) + (1 - alpha) * anchor,
                "Paired mean": data.source[i] + mean,
                "Unpaired mean": data.source[i] + unpaired,
                "Nearest pair": data.source[i] + data.shifts[nearest],
                "No change": data.source[i], "Target mean": anchor,
                "PCHT without graph": alpha * (data.source[i] + fit_shift(
                    data, train, contrast, shrinkage, 0)) + (1 - alpha) * anchor,
                "PCHT without context": alpha * data.source[i] + (1 - alpha) * anchor,
            }
            for method, aa in selected_alpha.items():
                methods[method] = aa * data.source[i] + (1 - aa) * reference_anchor(
                    data, train, contrast, method)
            for name, prediction in methods.items():
                for j, gene in enumerate(data.genes):
                    predictions.append({"pair_id": i, "model": model, "contrast": contrast,
                                        "gene": gene, "method": name, "source": data.source[i, j],
                                        "observed": data.target[i, j],
                                        "predicted": np.clip(prediction, 0, 1)[j]})
            for la, et, al in GRID:
                prediction = np.clip(al * (data.source[i] + fit_shift(data, train, contrast, la, et))
                                     + (1 - al) * anchor, 0, 1)
                sensitivity.append({"model": model, "pair_id": i, "contrast": contrast,
                                    "lambda": la, "eta": et, "alpha": al,
                                    "mse": float(np.mean((prediction - data.target[i]) ** 2))})
    pred = pd.DataFrame(predictions)
    pred["squared_error"] = (pred.predicted - pred.observed) ** 2
    metrics = pred.groupby(["model", "contrast", "method"]).squared_error.mean().pow(
        0.5).rename("rmse").reset_index()
    return {"heldout_predictions": pred, "nested_fold_parameters": pd.DataFrame(fold_rows),
            "baseline_fold_parameters": pd.DataFrame(reference_rows), "model_validation": metrics,
            "parameter_sensitivity": pd.DataFrame(sensitivity)}


def _signed_rank(values) -> float:
    values = np.asarray(values)
    return float(wilcoxon(values).pvalue) if np.any(values != 0) else 1.0


def validation_summaries(predictions: pd.DataFrame, resamples=10000,
                         seed=20260923, comparison_seed=20260926):
    """Calculate paired-model percentile intervals and signed-rank comparisons."""
    predictions = predictions.copy()
    predictions["squared_error"] = (predictions.predicted - predictions.observed) ** 2
    tables = {}
    for contrast in ["All", "Format", "Medium"]:
        subset = predictions if contrast == "All" else predictions[predictions.contrast == contrast]
        if not subset.empty:
            tables[contrast] = subset.groupby(["model", "method"]).squared_error.mean().pow(0.5).unstack()
    rng = np.random.default_rng(seed)
    rows = []
    for references_only in [False, True]:
        for contrast, errors in tables.items():
            methods = list(REFERENCE_METHODS) if references_only else [
                m for m in errors.columns if m not in REFERENCE_METHODS]
            for method in methods:
                draws = rng.integers(0, len(errors), (resamples, len(errors)))
                boot = 1 - errors[method].to_numpy()[draws].mean(axis=1) / errors[
                    "No change"].to_numpy()[draws].mean(axis=1)
                rows.append({"contrast": contrast, "method": method, "models": len(errors),
                             "mean_rmse": errors[method].mean(),
                             "relative_reduction": 1 - errors[method].mean() / errors["No change"].mean(),
                             "ci_low": np.quantile(boot, .025), "ci_high": np.quantile(boot, .975),
                             "paired_p": _signed_rank(errors["No change"] - errors[method])})
    comparison_rng = np.random.default_rng(comparison_seed)
    comparisons = []
    for contrast, errors in tables.items():
        for method in REFERENCE_METHODS:
            difference = errors[method].to_numpy() - errors.PCHT.to_numpy()
            draws = comparison_rng.integers(0, len(errors), (resamples, len(errors)))
            boot = difference[draws].mean(axis=1)
            comparisons.append({"contrast": contrast, "method": method, "reference": "PCHT",
                                "models": len(errors), "mean_rmse_difference": difference.mean(),
                                "ci_low": np.quantile(boot, .025), "ci_high": np.quantile(boot, .975),
                                "paired_p": _signed_rank(difference),
                                "models_lower_error": int((difference < 0).sum())})
    return pd.DataFrame(rows), pd.DataFrame(comparisons)


def paired_effects(data: ContextBenchmark):
    """Summarize model-balanced paired gene effects and functional-class shifts."""
    effects, pathways = [], []
    for contrast in ["Format", "Medium"]:
        ix = data.pairs.index[data.pairs.contrast.eq(contrast)].to_numpy()
        if len(ix) == 0:
            continue
        balanced = pd.DataFrame(data.shifts[ix], index=data.pairs.loc[ix, "model"],
                                columns=data.genes).groupby(level=0).mean()
        for gene in data.genes:
            values = balanced[gene]
            effects.append({"contrast": contrast, "gene": gene, "models": len(values),
                            "mean_rank_shift": values.mean(), "median_rank_shift": values.median(),
                            "sd": values.std(), "negative_fraction": (values < 0).mean(),
                            "p": _signed_rank(values)})
        for functional_class in data.classes:
            gene_indices = np.flatnonzero(_boolean_column(
                data.library.loc[list(data.genes), functional_class]).to_numpy())
            if gene_indices.size == 0:
                continue
            for model in data.pairs.loc[ix, "model"].unique():
                model_indices = ix[data.pairs.loc[ix, "model"].eq(model).to_numpy()]
                pathways.append({"contrast": contrast, "model": model,
                                 "pathway": functional_class, "genes": len(gene_indices),
                                 "rank_shift": data.shifts[np.ix_(model_indices, gene_indices)].mean()})
    effects = pd.DataFrame(effects)
    effects["q"] = effects.groupby("contrast").p.transform(
        lambda p: multipletests(p, method="fdr_bh")[1])
    return effects, pd.DataFrame(pathways)


def analyze_context_transfer(effects, metadata, library, resamples=10000):
    data = build_benchmark(effects, metadata, library)
    results = nested_predictions(data)
    results["validation_summary"], results["baseline_comparisons"] = validation_summaries(
        results["heldout_predictions"], resamples=resamples)
    results["paired_gene_effects"], results["pathway_shifts"] = paired_effects(data)
    results["screen_rank_matrix"] = data.ranks
    results["paired_contexts"] = data.pairs
    results["paired_rank_shifts"] = pd.DataFrame(data.shifts, index=data.pairs.index, columns=data.genes)
    return data, results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    effects = pd.read_csv(args.input_dir / "screen_effects.csv").set_index("screen_id")
    metadata = pd.read_csv(args.input_dir / "screen_metadata.csv").set_index("screen_id")
    library = pd.read_csv(args.input_dir / "library.csv").set_index("gene")
    data, results = analyze_context_transfer(effects, metadata, library)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    indexed = {"screen_rank_matrix", "paired_contexts", "paired_rank_shifts"}
    for name, table in results.items():
        suffix = ".csv.gz" if name == "heldout_predictions" else ".csv"
        output = table.drop(columns="squared_error") if name == "heldout_predictions" else table
        output.to_csv(args.output_dir / (name + suffix), index=name in indexed)
    design = {"name": "Paired Context Harmonic Transfer", "abbreviation": "PCHT",
              "screens": len(data.ranks), "library_genes": data.ranks.shape[1],
              "test_genes": len(data.genes), "models": data.metadata.model.nunique(),
              "paired_contrasts": len(data.pairs),
              "rank_direction": "Lower percentiles indicate stronger depletion.",
              "lambda_grid": list(LAMBDA_GRID), "eta_grid": list(ETA_GRID),
              "alpha_grid": list(ALPHA_GRID), "seed": 20260923,
              "baseline_comparison_seed": 20260926, "bootstrap_resamples": 10000,
              "estimators": list(ESTIMATORS), "ablations": list(ABLATIONS),
              "outer_split": "All conditions of one model are held out together.",
              "inner_split": "All conditions of each training model are held out together."}
    (args.output_dir / "design.json").write_text(json.dumps(design, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
