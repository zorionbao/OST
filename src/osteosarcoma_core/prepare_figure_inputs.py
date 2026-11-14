"""Prepare publication panel tables from externally supplied analysis results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.cbook import boxplot_stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import gaussian_kde
from sklearn.decomposition import PCA

PROGRAMMES = ["Replication", "Redox", "Matrix"]
PROJECTION_METHODS = ["Pathway ranks", "CRMP fixed", "Pure-state CRMP"]
TRANSFER_METHODS = ["No change", "Target mean", "Paired mean", "Unpaired mean", "Nearest pair", "Global anchor", "Paired midpoint", "PCHT"]


def prepare(input_dir, output_dir):
    source, destination = Path(input_dir), Path(output_dir)
    if source.resolve() == destination.resolve():
        raise ValueError("Input and output directories must differ.")
    destination.mkdir(parents=True, exist_ok=True)
    records = []

    def read(name, indexed=False):
        aliases = {"library_composition": ["library"], "cell_model_scores": ["model_scores"],
                   "screen_metadata": ["metadata"], "patient_scores": ["patient_programme_scores"]}
        candidates = [source / f"{base}.{suffix}" for base in [name, *aliases.get(name, [])] for suffix in ["csv", "csv.gz"]]
        path = next((candidate for candidate in candidates if candidate.exists()), candidates[0])
        return pd.read_csv(path, index_col=0 if indexed else None)

    def emit(panel, data, *, indexed=False, tables=None, metadata=None):
        datasets = {"data": (data, indexed)}
        datasets.update(tables or {})
        descriptors = {}
        for key, (table, use_index) in datasets.items():
            filename = f"{panel}.csv" if key == "data" else f"{panel}_{key}.csv"
            table.to_csv(destination / filename, index=use_index, float_format="%.17g")
            descriptors[key] = {"file": filename, "indexed": use_index,
                                "columns": [str(column) for column in table.columns], "rows": len(table)}
        settings = {"panel": panel, "datasets": descriptors, **(metadata or {})}
        (destination / f"{panel}.json").write_text(json.dumps(settings, indent=2, allow_nan=False), encoding="utf-8")
        records.append({"panel": panel, "rows": len(data), "datasets": list(descriptors)})

    pseudobulk = read("pseudobulk_programmes")
    means = pseudobulk.groupby("compartment")[PROGRAMMES].mean()
    rows = means.rename_axis("compartment").reset_index().melt("compartment", var_name="programme", value_name="value")
    low, high = float(rows.value.min()), float(rows.value.max())
    rows["area"] = 15 + 90 * (rows.value - low) / (high - low if high > low else 1)
    emit("figure_1a", rows, metadata={"compartments": means.index.tolist(), "programmes": PROGRAMMES, "vmin": low, "vmax": high})
    purity = read("purity_metrics")
    selection = purity[(purity.cohort == "GSE162454") & (purity.programme == "Matrix") & purity.method.isin(PROJECTION_METHODS)]
    curves = selection.groupby(["method", "rna_fraction"]).absolute_error.mean().rename("value").reset_index()
    emit("figure_1b", curves, metadata={"method_order": PROJECTION_METHODS})
    donor = read("donor_metrics")
    improvements, densities, quantiles = [], [], []
    for programme in PROGRAMMES:
        table = donor[(donor.cohort == "GSE162454") & (donor.programme == programme)].pivot(index="donor", columns="method", values="rmse")
        values = 100 * (1 - table["Pure-state CRMP"] / table["Pathway ranks"])
        improvements.extend({"donor": str(identifier), "programme": programme, "value": float(value)} for identifier, value in values.items())
        grid = np.linspace(values.min() - .05 * np.ptp(values), values.max() + .05 * np.ptp(values), 180)
        if values.std() > 1e-12:
            density = gaussian_kde(values)(grid)
            densities.extend({"programme": programme, "x": float(x), "height": float(y / density.max() * .33)} for x, y in zip(grid, density))
        q1, median, q3 = np.quantile(values, [.25, .5, .75])
        quantiles.append({"programme": programme, "q1": float(q1), "median": float(median), "q3": float(q3)})
    emit("figure_1c", pd.DataFrame(improvements), tables={"density": (pd.DataFrame(densities, columns=["programme", "x", "height"]), False), "quantiles": (pd.DataFrame(quantiles), False)}, metadata={"programme_order": PROGRAMMES})
    pure = read("pure_predictions")
    external_pure = pure[pure.cohort == "GSE162454"]
    calibration_points = external_pure[(external_pure.programme == "Matrix") & external_pure.method.isin(["CRMP fixed", "Pure-state CRMP"])]
    emit("figure_1d", calibration_points[["donor", "method", "truth", "prediction"]], metadata={"identity_min": float(calibration_points.truth.min() - .025), "identity_max": float(calibration_points.truth.max() + .025)})
    errors = external_pure[external_pure.method == "Pure-state CRMP"].copy()
    errors["value"] = errors.prediction - errors.truth
    boxes = []
    for programme in PROGRAMMES:
        stats = boxplot_stats(errors[errors.programme == programme].value.to_numpy(), whis=1.5)[0]
        boxes.append({"programme": programme, **{key: float(stats[key]) for key in ["med", "q1", "q3", "whislo", "whishi"]}})
    emit("figure_1e", errors[["donor", "programme", "value"]], tables={"boxes": (pd.DataFrame(boxes), False)}, metadata={"programme_order": PROGRAMMES})
    summary = read("crmp_summary")
    calibration = read("pure_calibration")
    mixture_metrics = summary[(summary.cohort == "GSE162454") & (summary.programme == "Matrix")].set_index("method")
    pure_metrics = calibration[(calibration.cohort == "GSE162454") & (calibration.programme == "Matrix")].set_index("method")
    values = np.array([[mixture_metrics.loc[method, "mean_rmse"], pure_metrics.loc[method, "rmse"], pure_metrics.loc[method, "sd_ratio"]] for method in PROJECTION_METHODS])
    scale = np.array([values[:, 0].max(), values[:, 1].max(), 1.])
    metrics = [{"method": method, "metric": metric, "value": float(value)} for method, row in zip(PROJECTION_METHODS, values / scale) for metric, value in zip(["Mixture RMSE", "Pure RMSE", "Score SD"], row)]
    emit("figure_1f", pd.DataFrame(metrics), metadata={"method_order": PROJECTION_METHODS, "metric_order": ["Mixture RMSE", "Pure RMSE", "Score SD"], "normalization": scale.tolist()})
    intervals = summary[(summary.cohort == "GSE162454") & (summary.method == "Pure-state CRMP")].set_index("programme").loc[PROGRAMMES].reset_index()
    for column in ["relative_reduction", "ci_low", "ci_high"]:
        intervals[column] *= 100
    emit("figure_1g", intervals[["programme", "relative_reduction", "ci_low", "ci_high"]], metadata={"programme_order": PROGRAMMES})
    surface = read("inner_parameter_surface")
    surface = surface[(surface.method == "Pure-state CRMP") & (surface.outer_fold != "Final")]
    surface = surface.groupby(["lambda", "alpha"]).inner_donor_mse.mean().rename("value").reset_index()
    surface["x"] = np.log10(surface.alpha)
    surface["value"] *= 1000
    emit("figure_1h", surface[["x", "lambda", "value"]], metadata={"x_transform": "log10 alpha", "value_unit": "MSE × 1000", "contour_levels": 7})

    screen = read("context_specific_screen", indexed=True)
    screen = screen[screen.named_gene.astype(bool)].copy()
    screen["category"] = np.select([screen.lung_hit & screen.culture_hit, screen.lung_hit & ~screen.culture_hit], ["Both contexts", "Lung threshold only"], "Other genes")
    points = screen.rename_axis("gene").reset_index().rename(columns={"log2FC_culture": "x", "log2FC_lung": "y"})
    emit("figure_2a", points[["gene", "x", "y", "category"]])
    curves = []
    for column, category in [("log2FC_lung", "Lung"), ("log2FC_culture", "Culture")]:
        values = np.sort(screen[column].dropna().to_numpy())
        curves.extend({"category": category, "x": float(x), "y": float(y)} for x, y in zip(values, np.arange(1, len(values) + 1) / len(values)))
    emit("figure_2b", pd.DataFrame(curves))
    effects = read("osteosarcoma_gene_effects", indexed=True)
    selected_genes = [gene for gene in screen[screen.lung_hit & ~screen.culture_hit].sort_values("log2FC_lung").index if gene in effects.columns]
    dependency = effects[selected_genes].rename_axis("model").reset_index().melt("model", var_name="gene", value_name="value")
    emit("figure_2c", dependency, metadata={"gene_order": selected_genes})
    patients, models = read("patient_scores", indexed=True), read("cell_model_scores", indexed=True)
    states = pd.concat([pd.DataFrame({"identifier": patients.index, "x": patients.CRMP_Matrix, "y": patients.CRMP_Replication, "category": "Patient tumours"}), pd.DataFrame({"identifier": models.index, "x": models.Matrix, "y": models.Replication, "category": "Cell models"})], ignore_index=True)
    emit("figure_2d", states)
    allocation = read("patient_family_transport", indexed=True)
    displayed = allocation.loc[allocation.sum(axis=1) > 1e-12, allocation.sum(axis=0) > 1e-12]
    displayed = displayed.loc[displayed.sum(axis=1).sort_values(ascending=False).index] * len(allocation)
    family_labels = {"CCLFPEDS0019T": "CCLF\nPEDS\n0019T", "NCCOS1X2C1": "NCCOS1\nX2C1", "HSOS1": "HSOS\n1"}
    emit("figure_2e", displayed, indexed=True, metadata={"patient_order": displayed.index.astype(str).tolist(), "family_order": displayed.columns.tolist(), "family_display_labels": [family_labels.get(name, name) for name in displayed.columns], "patient_input_mass": 1 / len(allocation)})
    sensitivity = read("transport_sensitivity")
    sensitivity["represented_percent"] = sensitivity.represented_mass * 100
    emit("figure_2f", sensitivity[["rejection_cost", "represented_percent"]])
    benchmark = read("duplication_benchmark")
    invariant = benchmark[benchmark.method == "Ancestry quotient"]
    emit("figure_2g", invariant[["copies", "duplicated_family_mass", "kernel_family_mass"]])
    timings = benchmark[benchmark.method.isin(["Full sparse", "Ancestry quotient"])].copy()
    timings["milliseconds"] = timings.median_seconds * 1000
    emit("figure_2h", timings[["method", "models", "milliseconds"]])

    ranks = read("screen_rank_matrix", indexed=True)
    library = read("library_composition", indexed=True)
    library = library.rename(columns={"all_test_genes": "All test genes"})
    genes = library.index[library["All test genes"].astype(bool)].intersection(ranks.columns)
    rank_matrix = ranks[genes]
    order = leaves_list(linkage(rank_matrix.T, method="average", metric="euclidean"))
    ordered = rank_matrix.iloc[:, order]
    emit("figure_3a", ordered, indexed=True, metadata={"screen_order": ordered.index.astype(str).tolist(), "gene_order": ordered.columns.tolist(), "clustering": "Average-linkage Euclidean gene profiles"})
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(rank_matrix)
    metadata = read("screen_metadata", indexed=True)
    metadata = metadata.rename(columns={"growth_format": "Growth Format"})
    coordinates = pd.DataFrame({"screen": rank_matrix.index, "x": coordinates[:, 0], "y": coordinates[:, 1], "format": metadata.loc[rank_matrix.index, "Growth Format"].to_numpy()})
    emit("figure_3b", coordinates, metadata={"feature_order": rank_matrix.columns.tolist(), "screen_order": rank_matrix.index.astype(str).tolist(), "explained_variance_ratio": pca.explained_variance_ratio_.tolist(), "pca_mean": pca.mean_.tolist(), "pca_components": pca.components_.tolist()})
    paired = read("paired_gene_effects").pivot(index="gene", columns="contrast", values="mean_rank_shift")
    emit("figure_3c", paired.rename(columns={"Format": "x", "Medium": "y"}).rename_axis("gene").reset_index())
    predictions = read("heldout_predictions")
    predictions["squared_error"] = (predictions.predicted - predictions.observed) ** 2
    model_rmse = predictions.groupby(["model", "method"]).squared_error.mean().pow(.5).unstack()[TRANSFER_METHODS]
    emit("figure_3d", model_rmse.rename_axis("model").reset_index().melt("model", var_name="method", value_name="rmse"), metadata={"model_order": model_rmse.index.tolist(), "method_order": TRANSFER_METHODS})
    emit("figure_3e", model_rmse, indexed=True, metadata={"model_order": model_rmse.index.tolist(), "method_order": TRANSFER_METHODS})
    pcht_predictions = predictions[predictions.method == "PCHT"]
    emit("figure_3f", pcht_predictions[["pair_id", "model", "gene", "observed", "predicted"]], metadata={"hexbin_gridsize": 26})
    pathway = read("pathway_shifts").groupby(["pathway", "contrast"]).rank_shift.mean().unstack()
    emit("figure_3g", pathway, indexed=True, metadata={"pathway_order": pathway.index.tolist(), "contrast_order": pathway.columns.tolist()})
    atlas = read("osteosarcoma_context_gene_atlas", indexed=True).dropna(subset=["OS_Chronos_mean", "MeanDifference"])
    atlas = atlas.rename_axis("gene").reset_index().rename(columns={"OS_Chronos_mean": "x", "MeanDifference": "y", "Format": "color"})
    emit("figure_3h", atlas[["gene", "x", "y", "color"]], metadata={"color_limit": float(atlas.color.abs().max())})

    kp, pt = read("kp_spots"), read("pt_spots")
    categories = ["KP", "Tgfbr2", "Jak2", "periphery", "Unassigned"]
    spatial_annotations = pd.DataFrame({"x": kp.array_x, "y": -kp.array_y, "category": kp.perturbation})
    emit("figure_4a", spatial_annotations, metadata={"category_order": categories, "vertical_coordinate": "Negated array_y"})
    signature = "Fibroblast / TGFB1"
    for panel, frame, column in [("figure_4b", kp, signature), ("figure_4c", pt, signature), ("figure_4d", pt, "CRMP_Matrix")]:
        data = pd.DataFrame({"x": frame.array_x, "y": -frame.array_y, "value": frame[column]})
        emit(panel, data, metadata={"vmin": float(data.value.min()), "vmax": float(data.value.max()), "score": column, "vertical_coordinate": "Negated array_y"})
    inference = read("signature_inference").pivot(index="signature", columns="programme", values="median_rho")
    inference.index = inference.index.str.replace("Lung_tumor", "Tumour").str.replace(" / ", " ")
    inference = inference[["CRMP_Matrix", "CRMP_Redox", "Hypoxia"]].rename(columns={"CRMP_Matrix": "Matrix", "CRMP_Redox": "Redox"})
    emit("figure_4e", inference, indexed=True, metadata={"signature_order": inference.index.tolist(), "programme_order": inference.columns.tolist(), "color_limit": max(.5, float(inference.abs().max().max()))})
    contrasts = read("perturb_section_contrasts")
    contrasts = contrasts[contrasts.signature == signature]
    emit("figure_4f", contrasts[["sample", "perturbation", "difference"]], metadata={"perturbation_order": sorted(contrasts.perturbation.unique())})
    coupling = read("spatial_coupling")
    coupling.index = [str(value) for value in range(len(coupling))]
    emit("figure_4g", coupling, indexed=True, metadata={"primary_order": coupling.index.tolist(), "relapse_order": coupling.columns.tolist(), "vmax": float(coupling.to_numpy().max()), "mass_sum": float(coupling.to_numpy().sum())})
    changes = read("paired_spatial_transport").set_index("patient")[["CRMP_Replication", "CRMP_Redox", "CRMP_Matrix", "Hypoxia", "Immune", "Endothelial"]]
    changes.columns = ["Replication", "Redox", "Matrix", "Hypoxia", "Immune", "Endothelial"]
    emit("figure_4h", changes, indexed=True, metadata={"patient_order": changes.index.astype(str).tolist(), "programme_order": changes.columns.tolist(), "color_limit": max(.16, float(changes.abs().max().max()))})
    manifest = {"panel_count": len(records), "panels": records}
    (destination / "panel_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.input_dir, args.output_dir)))


if __name__ == "__main__":
    main()
