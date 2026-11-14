"""Gene-disjoint spatial response scores and mass-preserving metacell transport."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import logsumexp
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, wilcoxon
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


PRIMARY_PROGRAMMES = ("CRMP_Matrix", "CRMP_Redox", "Hypoxia")
TRANSPORT_FEATURES = ("CRMP_Replication", "CRMP_Redox", "CRMP_Matrix", "Hypoxia", "Immune", "Endothelial")
FIBROBLAST_MARKERS = ("DCN", "LUM", "COL3A1", "PDGFRA", "DPT", "CFD")
PRIMARY_COVARIATES = ("counts", "Immune", "Osteogenic")


def disjoint_signatures(signatures: pd.DataFrame, feature_genes, programme_gene_sets) -> pd.DataFrame:
    """Exclude representation features and all programme genes by human symbol."""
    required = {"cell_context", "treatment", "gene", "avg_fc"}
    if not required.issubset(signatures.columns):
        raise ValueError(f"Signatures require {sorted(required)}.")
    excluded = set(feature_genes)
    for genes in programme_gene_sets.values():
        excluded.update(genes)
    result = signatures[signatures.cell_context.isin(["Lung_tumor", "Fibroblast"])].copy()
    result["signature"] = result.cell_context + " / " + result.treatment
    result = result[~result.gene.isin(excluded)].copy()
    if result.empty or not np.isfinite(result.avg_fc.to_numpy(float)).all():
        raise ValueError("Gene-disjoint signatures require finite fold-change weights.")
    return result


def _validate_expression(expression, genes):
    expression = sparse.csr_matrix(expression, dtype=float)
    if expression.shape[0] != len(genes):
        raise ValueError("Expression rows must match gene symbols in row order.")
    if not np.isfinite(expression.data).all() or (expression.data < 0).any():
        raise ValueError("Expression must contain finite nonnegative counts.")
    return expression


def signature_scores(expression, genes, signatures: pd.DataFrame, species="human"):
    """Compute signed L1-weighted log1p counts-per-10,000 response scores."""
    expression = _validate_expression(expression, genes)
    if species not in {"human", "mouse"}:
        raise ValueError("Species must be human or mouse.")
    column = "gene" if species == "human" else "mouse_gene"
    if column not in signatures:
        raise ValueError(f"Signatures require {column} for {species} expression.")
    gene_index = {gene: i for i, gene in enumerate(genes)}
    weights = np.zeros((len(genes), signatures.signature.nunique()))
    names, coverage = [], []
    for j, (name, definitions) in enumerate(signatures.groupby("signature", sort=True)):
        values = definitions[definitions[column].notna()].groupby(column).avg_fc.mean()
        values = values[values.index.isin(gene_index)]
        if len(values):
            denominator = np.abs(values).sum()
            if denominator == 0:
                raise ValueError(f"Matched signature weights have zero absolute sum: {name}.")
            weights[[gene_index[g] for g in values.index], j] = values / denominator
        names.append(name)
        coverage.append({"signature": name, "species": species, "genes": len(values)})
    counts = np.asarray(expression.sum(axis=0)).ravel()
    normalized = (expression @ sparse.diags(10000 / np.maximum(counts, 1))).tocsr()
    normalized.data = np.log1p(normalized.data)
    scores = pd.DataFrame(np.asarray(normalized.T @ weights), columns=names)
    return scores, pd.DataFrame(coverage)


def marker_score(expression, genes, markers):
    """Average log1p counts-per-10,000 over a fixed marker set."""
    expression = _validate_expression(expression, genes)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    missing = set(markers) - set(gene_index)
    if missing:
        raise ValueError(f"Marker expression is missing: {sorted(missing)}.")
    if not markers:
        raise ValueError("The marker set must contain at least one gene.")
    counts = np.asarray(expression.sum(axis=0)).ravel()
    values = expression[[gene_index[g] for g in markers], :].toarray()
    return np.log1p(values * 10000 / np.maximum(counts, 1)).mean(axis=0)


def ols_residuals(outcomes, covariates):
    """Residualize every outcome separately on a common covariate matrix."""
    outcomes = np.asarray(outcomes, dtype=float)
    covariates = np.asarray(covariates, dtype=float)
    if outcomes.shape[0] != covariates.shape[0]:
        raise ValueError("Outcomes and covariates require the same rows.")
    if not np.isfinite(outcomes).all() or not np.isfinite(covariates).all():
        raise ValueError("Outcomes and covariates must be finite.")
    return outcomes - covariates @ np.linalg.lstsq(covariates, outcomes, rcond=None)[0]


def section_correlations(spots: pd.DataFrame, signature_names,
                         programmes=PRIMARY_PROGRAMMES, composition=False) -> pd.DataFrame:
    """Calculate within-section Spearman correlations of OLS residual scores."""
    rows = []
    signature_names = list(signature_names)
    programmes = list(programmes)
    required = {"sample", "patient", *PRIMARY_COVARIATES, *signature_names, *programmes}
    if composition:
        required.update(["fibroblast_score", "Endothelial"])
    if not required.issubset(spots.columns):
        raise ValueError(f"Spatial scores require {sorted(required - set(spots.columns))}.")
    for sample, section in spots.groupby("sample", sort=True):
        if section.patient.nunique() != 1:
            raise ValueError("Each section must belong to one patient.")
        primary = np.column_stack([np.ones(len(section)), np.log1p(section.counts),
                                   section.Immune, section.Osteogenic])
        adjustments = {"Primary adjustment": primary}
        if composition:
            adjustments.update({
                "Plus fibroblast": np.column_stack([primary, section.fibroblast_score]),
                "Plus fibroblast and endothelial": np.column_stack(
                    [primary, section.fibroblast_score, section.Endothelial]),
            })
        outcomes = section[signature_names + programmes].to_numpy(float)
        for adjustment, covariates in adjustments.items():
            residuals = ols_residuals(outcomes, covariates)
            for j, signature in enumerate(signature_names):
                for k, programme in enumerate(programmes):
                    rho = float(spearmanr(residuals[:, j], residuals[:, len(signature_names) + k]).statistic)
                    record = {"sample": sample, "patient": section.patient.iloc[0],
                              "signature": signature, "programme": programme, "rho": rho}
                    if composition:
                        record["adjustment"] = adjustment
                    rows.append(record)
    return pd.DataFrame(rows)


def patient_correlations(section_results: pd.DataFrame) -> pd.DataFrame:
    """Give each patient one mean correlation across available sections."""
    columns = ["patient", "signature", "programme"]
    if "adjustment" in section_results:
        columns.append("adjustment")
    return section_results.groupby(columns, sort=True).rho.mean().reset_index()


def correlation_inference(patient_results: pd.DataFrame, seed=20260923,
                          resamples=5000) -> pd.DataFrame:
    """Bootstrap patient medians, test signed ranks and adjust complete families."""
    grouped = ["signature", "programme"]
    composition = "adjustment" in patient_results
    if composition:
        grouped.insert(0, "adjustment")
    rng = np.random.default_rng(seed)
    records = []
    for key, patients in patient_results.groupby(grouped, sort=True):
        values = patients.rho.to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError("Patient correlations must be finite for inference.")
        draws = rng.integers(0, len(values), (resamples, len(values)))
        intervals = np.quantile(np.median(values[draws], axis=1), [.025, .975])
        record = dict(zip(grouped, key))
        record.update({"patients": len(values), "median_rho": float(np.median(values)),
                       "ci_low": float(intervals[0]), "ci_high": float(intervals[1]),
                       "p": float(wilcoxon(values).pvalue) if np.any(values != 0) else 1.0})
        if composition:
            record["positive_patients"] = int((values > 0).sum())
        records.append(record)
    inference = pd.DataFrame(records)
    if composition:
        inference["q"] = inference.groupby("adjustment").p.transform(
            lambda p: multipletests(p, method="fdr_bh")[1])
        inference["FDR_family_size"] = inference.groupby("adjustment").p.transform("size")
    else:
        inference["q"] = multipletests(inference.p, method="fdr_bh")[1]
    return inference


def composition_sensitivity(expression, genes, spots, signatures, programme_gene_sets,
                            feature_genes=()):
    """Add fixed fibroblast and endothelial covariates to section associations."""
    fibroblast = [g for g in FIBROBLAST_MARKERS if g not in programme_gene_sets["Matrix"]]
    endothelial = list(programme_gene_sets["Endothelial"])
    response_genes = set(signatures.gene)
    if set(fibroblast + endothelial) & response_genes:
        raise ValueError("Composition markers and response genes must be disjoint.")
    result = spots.copy()
    result["fibroblast_score"] = marker_score(expression, genes, fibroblast)
    correlations = section_correlations(result, sorted(signatures.signature.unique()), composition=True)
    overlap = {"fibroblast_markers": fibroblast,
               "excluded_matrix_markers": sorted(set(FIBROBLAST_MARKERS) & set(programme_gene_sets["Matrix"])),
               "endothelial_markers": endothelial,
               "fibroblast_CRMP_overlap": sorted(set(fibroblast) & set(feature_genes)),
               "endothelial_CRMP_overlap": sorted(set(endothelial) & set(feature_genes)),
               "response_signature_overlap": []}
    return result, correlations, overlap


def perturbation_cluster_contrasts(spots: pd.DataFrame, signature_names):
    """Compare section-wise means of lesion-cluster medians with KP clusters."""
    if not {"sample", "phenotypes"}.issubset(spots.columns):
        raise ValueError("Perturbation spots require sample and phenotypes.")
    spots = spots.copy()
    spots["perturbation"] = spots.phenotypes.str.extract(
        r"^(Tgfbr2|Ifngr2|Jak2|KP|periphery)", expand=False).fillna("Unassigned")
    lesions = spots.dropna(subset=["phenotypes"]).groupby(
        ["sample", "phenotypes", "perturbation"], sort=True)[list(signature_names)].median().reset_index()
    rows = []
    for sample, section in lesions.groupby("sample", sort=True):
        for perturbation in ["Tgfbr2", "Ifngr2", "Jak2"]:
            treated = section[section.perturbation == perturbation]
            control = section[section.perturbation == "KP"]
            if treated.empty or control.empty:
                continue
            for signature in signature_names:
                rows.append({"sample": sample, "perturbation": perturbation, "signature": signature,
                             "perturbed_clusters": len(treated), "KP_reference_clusters": len(control),
                             "difference": treated[signature].mean() - control[signature].mean()})
    return lesions, pd.DataFrame(rows)


def log_sinkhorn(cost, source_mass, target_mass, epsilon=.1,
                 max_iterations=2000, tolerance=1e-11):
    """Balance an entropic transport plan using alternating log-space updates."""
    cost = np.asarray(cost, dtype=float)
    source_mass, target_mass = np.asarray(source_mass, dtype=float), np.asarray(target_mass, dtype=float)
    if cost.shape != (len(source_mass), len(target_mass)) or not np.isfinite(cost).all():
        raise ValueError("Cost dimensions must match finite source and target masses.")
    if (source_mass <= 0).any() or (target_mass <= 0).any() or epsilon <= 0:
        raise ValueError("Masses and epsilon must be positive.")
    if not np.isclose(source_mass.sum(), target_mass.sum(), atol=1e-12, rtol=0):
        raise ValueError("Source and target masses must have the same total.")
    log_kernel = -cost / epsilon
    u, v = np.zeros(len(source_mass)), np.zeros(len(target_mass))
    log_a, log_b = np.log(source_mass), np.log(target_mass)
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        previous = u.copy()
        u = log_a - logsumexp(log_kernel + v[None, :], axis=1)
        v = log_b - logsumexp(log_kernel + u[:, None], axis=0)
        if np.max(np.abs(u - previous)) < tolerance:
            break
    plan = np.exp(log_kernel + u[:, None] + v[None, :])
    error = max(np.max(np.abs(plan.sum(axis=1) - source_mass)),
                np.max(np.abs(plan.sum(axis=0) - target_mass)))
    return plan, {"marginal_error": float(error), "iterations": iterations}


def spatial_metacells(spots: pd.DataFrame, features=TRANSPORT_FEATURES,
                      clusters=16, seed=20260923) -> pd.DataFrame:
    """Partition section coordinates and retain each metacell's spot mass."""
    columns = list(features) + ["array_x", "array_y"]
    if len(spots) < clusters or not np.isfinite(spots[columns].to_numpy(float)).all():
        raise ValueError("Metacells require finite scores and at least 16 spots.")
    if np.unique(spots[["array_x", "array_y"]].to_numpy(), axis=0).shape[0] < clusters:
        raise ValueError("Metacells require at least 16 distinct coordinates.")
    with threadpool_limits(limits=1, user_api="blas"), threadpool_limits(limits=2, user_api="openmp"):
        labels = KMeans(clusters, random_state=seed, n_init=10).fit(
            spots[["array_x", "array_y"]]).labels_
    section = spots.copy()
    section["metacell"] = labels
    nodes = section.groupby("metacell", sort=True)[columns].mean()
    nodes["mass"] = section.groupby("metacell").size() / len(section)
    return nodes


def metacell_transport(source: pd.DataFrame, target: pd.DataFrame,
                       features=TRANSPORT_FEATURES):
    """Transport all spot mass and report paired spot-mean programme differences."""
    features = list(features)
    a, b = spatial_metacells(source, features), spatial_metacells(target, features)
    standard_deviation = pd.concat([source[features], target[features]]).std().replace(0, 1)
    cost = cdist(a[features] / standard_deviation, b[features] / standard_deviation,
                 "sqeuclidean") / len(features)
    cost = cost / max(float(np.median(cost)), 1e-9)
    plan, diagnostics = log_sinkhorn(cost, a.mass.to_numpy(), b.mass.to_numpy())
    summary = {"primary_spots": len(source), "relapse_spots": len(target),
               "transport_cost": float((plan * cost).sum()), **diagnostics,
               **{feature: float(target[feature].mean() - source[feature].mean()) for feature in features}}
    return a, b, plan, summary


def paired_spatial_transport(spots: pd.DataFrame):
    """Evaluate every primary/relapse section pair within each patient."""
    records, nodes, plans = [], [], {}
    for patient, patient_spots in spots.groupby("patient", sort=True):
        primary = patient_spots[patient_spots.group == "Primary"]
        relapse = patient_spots[patient_spots.group == "Relapse"]
        for primary_sample, source in primary.groupby("sample", sort=True):
            for relapse_sample, target in relapse.groupby("sample", sort=True):
                a, b, plan, summary = metacell_transport(source, target)
                for sample, group, table in [(primary_sample, "Primary", a), (relapse_sample, "Relapse", b)]:
                    table = table.copy()
                    table["patient"], table["sample"], table["group"] = patient, sample, group
                    nodes.append(table.reset_index())
                plans[(primary_sample, relapse_sample)] = plan
                records.append({"patient": patient, "primary": primary_sample, "relapse": relapse_sample, **summary})
    metacells = (pd.concat(nodes, ignore_index=True).drop_duplicates(["sample", "metacell"])
                 if nodes else pd.DataFrame())
    return pd.DataFrame(records), metacells, plans


def _read_section(directory: Path, row):
    expression = sparse.load_npz(directory / row.expression_file).tocsr()
    genes = pd.read_csv(directory / row.genes_file)["gene"].tolist()
    metadata = pd.read_csv(directory / row.metadata_file, float_precision="round_trip")
    if "spot_id" not in metadata or not metadata.spot_id.is_unique:
        raise ValueError("Spot metadata requires unique spot_id values in expression-column order.")
    if expression.shape[1] != len(metadata):
        raise ValueError("Spot metadata rows must match expression columns.")
    expression = _validate_expression(expression, genes)
    computed_counts = np.asarray(expression.sum(axis=0)).ravel()
    if "counts" in metadata and not np.array_equal(computed_counts, metadata.counts.to_numpy()):
        raise ValueError("Metadata counts must equal expression-column sums.")
    metadata["counts"] = computed_counts
    metadata["sample"] = row.sample
    return expression, genes, metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    directory = args.input_dir
    signatures = pd.read_csv(directory / "signatures.csv")
    feature_genes = pd.read_csv(directory / "model_feature_genes.csv")["gene"].tolist()
    programme_gene_sets = json.loads((directory / "programme_gene_sets.json").read_text(encoding="utf-8"))
    signatures = disjoint_signatures(signatures, feature_genes, programme_gene_sets)
    sections = pd.read_csv(directory / "sections.csv")
    required = {"sample", "species", "analysis", "expression_file", "genes_file", "metadata_file"}
    if not required.issubset(sections.columns) or not sections["sample"].is_unique:
        raise ValueError(f"Section manifest requires unique samples and {sorted(required)}.")
    if not sections.analysis.isin(["osteosarcoma", "perturbation"]).all():
        raise ValueError("Section analysis must be osteosarcoma or perturbation.")
    names = sorted(signatures.signature.unique())
    primary_sections, composition_sections, os_means, os_spots = [], [], [], []
    perturb_spots, covers, marker_coverage, composition_covariates = [], [], [], []
    outputs, overlap = {}, None
    for row in sections.sort_values("sample").itertuples(index=False):
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(row.sample)):
            raise ValueError("Sample names require letters, digits, underscores, periods or hyphens.")
        expression, genes, metadata = _read_section(directory, row)
        scores, coverage = signature_scores(expression, genes, signatures, row.species)
        spots = metadata.copy()
        spots[names] = scores.to_numpy()
        covers.append(coverage)
        if row.analysis == "osteosarcoma":
            if row.species != "human":
                raise ValueError("Osteosarcoma sections require human expression.")
            primary_sections.append(section_correlations(spots, names))
            adjusted, correlations, overlap = composition_sensitivity(
                expression, genes, spots, signatures, programme_gene_sets, feature_genes)
            composition_sections.append(correlations)
            marker_coverage.append({"sample": row.sample, "patient": spots.patient.iloc[0],
                                    "spots": len(spots), "fibroblast_markers_available": len(overlap["fibroblast_markers"]),
                                    "fibroblast_markers": ";".join(overlap["fibroblast_markers"])})
            composition_covariates.append(adjusted[["spot_id", "sample", "patient", "fibroblast_score", "Endothelial"]].rename(
                columns={"Endothelial": "endothelial_score"}))
            os_means.append(spots.groupby(["sample", "patient", "group"])[names].mean().reset_index())
            os_spots.append(spots)
            outputs[row.sample + "_signature_spots"] = spots
        else:
            perturb_spots.append(spots)
            outputs[row.sample + "_perturb_spots"] = spots
    outputs["disjoint_signature_gene_definitions"] = signatures
    outputs["signature_gene_coverage"] = pd.concat(covers).drop_duplicates()
    plans = {}
    if primary_sections:
        correlations = pd.concat(primary_sections, ignore_index=True)
        patients = patient_correlations(correlations)
        inference = correlation_inference(patients)
        composition = pd.concat(composition_sections, ignore_index=True)
        composition_patients = patient_correlations(composition)
        outputs.update({"OS_section_signature_means": pd.concat(os_means, ignore_index=True),
                        "OS_section_signature_correlations": correlations,
                        "OS_patient_signature_correlations": patients,
                        "OS_signature_inference": inference,
                        "composition_section_correlations": composition,
                        "composition_patient_correlations": composition_patients,
                        "composition_signature_inference": correlation_inference(composition_patients, seed=20260926),
                        "composition_marker_coverage": pd.DataFrame(marker_coverage),
                        "composition_spot_covariates": pd.concat(composition_covariates, ignore_index=True)})
        transport, metacells, plans = paired_spatial_transport(pd.concat(os_spots, ignore_index=True))
        outputs["paired_spatial_transport"], outputs["spatial_metacells"] = transport, metacells
    if perturb_spots:
        lesions, contrasts = perturbation_cluster_contrasts(pd.concat(perturb_spots, ignore_index=True), names)
        outputs["perturb_cluster_signature_scores"], outputs["perturb_section_contrasts"] = lesions, contrasts
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in outputs.items():
        table.to_csv(args.output_dir / (name + ".csv"), index=False)
    for (primary, relapse), plan in plans.items():
        pd.DataFrame(plan).to_csv(args.output_dir / f"transport_{primary}_{relapse}.csv", index=False)
    if overlap is not None:
        (args.output_dir / "composition_gene_overlap.json").write_text(json.dumps(overlap, indent=2), encoding="utf-8")
    summary = {"signature_count": len(names), "OS_sections": len(primary_sections),
               "OS_patients": int(patients.patient.nunique()) if primary_sections else 0,
               "signature_associations": len(inference) if primary_sections else 0,
               "FDR_associations": int((inference.q < .05).sum()) if primary_sections else 0,
               "paired_spatial_comparisons": len(plans),
               "max_transport_marginal_error": float(transport.marginal_error.max())
                   if primary_sections and len(transport) else None,
               "response_seed": 20260923, "composition_seed": 20260926,
               "bootstrap_resamples": 5000, "metacells_per_section": 16,
               "kmeans_n_init": 10, "openmp_threads": 2, "blas_threads": 1,
               "sinkhorn_epsilon": .1, "sinkhorn_max_iterations": 2000,
               "sinkhorn_tolerance": 1e-11}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
