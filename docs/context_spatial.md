# Context transfer and spatial response analysis

The two modules expose numerical functions and explicit command-line entry points. Importing either module reads or writes no files. Inputs are supplied separately from this code repository. All commands use `--input-dir` and `--output-dir`.

## Paired Context Harmonic Transfer

```sh
python -m osteosarcoma_core.context_transfer --input-dir inputs/context --output-dir outputs/context
```

The input directory contains three CSV files:

| File | Required structure |
|---|---|
| `screen_effects.csv` | One row per screen. First identifier column: `screen_id`; remaining columns: gene symbols with finite numerical effects. Include the complete screened library for rank calculation. |
| `screen_metadata.csv` | `screen_id`, `model`, `growth_format`, `serum_status`. Screens retain their input order. `growth_format` uses `Plastic` or `Dome`; `serum_status` uses `Serum` or `Serum-free`. |
| `library.csv` | `gene`, `all_test_genes`, and one membership column per functional class. Membership values are `true`/`false` or `1`/`0`. The library row order defines the test-gene order. |

Average-tie ranks are computed within each screen over every input gene and divided by the complete library size. Test genes are selected afterwards. A format pair changes Plastic to Dome while holding serum status fixed. A medium pair changes Serum to Serum-free while holding growth format fixed. Pairs are formed within `model`.

Functional classes define adjacency `A = (Z / column_sums(Z)) Z.T`, with a zero diagonal; the denominator is bounded below by one. The graph includes library membership columns except `all_test_genes`, `Other`, `Actin Nucleation`, and `Gap Junction`. The Laplacian is `L = diag(A.sum(1)) - A`. The model-balanced mean paired shift is projected through eigenvalues of `L` with denominator `1 + lambda + eta * eigenvalue`.

Outer folds hold out every pair and screen of one model. Inner folds repeat whole-model removal among eligible training models. Each inner-model loss is the mean squared error across its eligible genes and pairs; these model losses are averaged equally. The fixed search order is the Cartesian product of `lambda = [0, 0.25, 1, 4]`, `eta = [0, 1, 4]`, and `alpha = [0, 0.25, 0.5, 0.75, 1]`. The first minimum determines the selected parameters. Predictions are clipped to `[0, 1]`.

The eight estimators are No change, Target mean, Paired mean, Unpaired mean, Nearest pair, Global anchor, Paired midpoint, and PCHT. Global anchor averages every eligible training screen within model and then across models. Paired midpoint averages paired midpoints within model and contrast, then across models. Both independently select `alpha` through inner whole-model folds. The two ablations retain the outer fold's selected parameters: PCHT without graph sets `eta = 0`; PCHT without context sets the explicit shift to zero.

For validation summaries, all eligible gene-pair squared errors are averaged within a model before taking the square root. Reported averages then weight models equally. Paired-model percentile intervals use 10,000 resamples with seed `20260923`; reference-anchor comparisons use seed `20260926`. The implementation preserves the random draw sequence across estimators and contrast groups. Two-sided Wilcoxon signed-rank tests compare paired model errors. Gene effects first average contrasts within model; Benjamini-Hochberg adjustment is performed separately for Format and Medium. Functional-class outputs average shifts across the selected class's genes and pairs within model.

Outputs include `screen_rank_matrix.csv`, `paired_contexts.csv`, `paired_rank_shifts.csv`, `heldout_predictions.csv.gz`, `nested_fold_parameters.csv`, `baseline_fold_parameters.csv`, `model_validation.csv`, `parameter_sensitivity.csv`, `validation_summary.csv`, `baseline_comparisons.csv`, `paired_gene_effects.csv`, `pathway_shifts.csv`, and `design.json`.

## Spatial response and transport

```sh
python -m osteosarcoma_core.spatial_context --input-dir inputs/spatial --output-dir outputs/spatial
```

The input directory contains:

| File | Required structure |
|---|---|
| `sections.csv` | `sample`, `species`, `analysis`, `expression_file`, `genes_file`, `metadata_file`. Each sample occurs once. `species` is `human` or `mouse`; `analysis` is `osteosarcoma` or `perturbation`. The three filename columns are relative to the input directory. |
| Expression NPZ per section | A SciPy sparse CSR matrix saved by `scipy.sparse.save_npz`: gene rows, spot columns, nonnegative counts. |
| Gene CSV per section | A `gene` column preserving expression-row order. Duplicate symbols use their final occurrence for indexed scoring; all expression rows contribute to spot count totals. |
| Spot metadata CSV per section | A unique `spot_id` column in expression-column order; all expression columns have one metadata row. A supplied `counts` column must equal the complete expression-column sum. Numeric metadata uses round-trip precision to preserve supplied coordinates. |
| `signatures.csv` | `cell_context`, `treatment`, `gene`, `avg_fc`; also `mouse_gene` for mouse sections. `avg_fc` contains signed response weights. Missing orthologues are empty cells. |
| `model_feature_genes.csv` | A `gene` column listing the fitted representation's feature genes. |
| `programme_gene_sets.json` | A JSON object mapping programme names to arrays of human gene symbols. Include all programme definitions, including `Matrix` and `Endothelial`. |

Osteosarcoma metadata additionally contains `patient`, `group`, `array_x`, `array_y`, `CRMP_Replication`, `CRMP_Redox`, `CRMP_Matrix`, `Hypoxia`, `Immune`, `Osteogenic`, and `Endothelial`. `group` uses `Primary` or `Relapse` for paired transport. These programme values are supplied from the frozen representation and programme-scoring analysis. At least 16 spots with distinct coordinates are required in sections participating in paired transport. Each section belongs to one patient.

Perturbation metadata additionally contains `phenotypes` and the desired spatial coordinates. All raw spots remain in expression-column order. Unannotated spots have an empty `phenotypes` cell. Lesion labels start with `Tgfbr2`, `Ifngr2`, `Jak2`, `KP`, or `periphery`.

Signatures use Lung_tumor and Fibroblast contexts. Human response genes overlapping either the model features or any programme gene set are removed before species mapping. Duplicate mapped symbols are averaged within signature. Matched signed weights are divided by their sum of absolute values. Expression is normalized to 10,000 counts per spot and transformed by `log1p`; each response score is its weighted sum. Coverage is recorded for every signature and species.

The primary association separately OLS-residualizes response and programme scores on an intercept, `log1p(counts)`, Immune, and Osteogenic, then calculates Spearman correlation within section. Section correlations are averaged within patient. Inference reports the median patient correlation, a 5,000-resample patient-bootstrap percentile interval, a two-sided signed-rank P value, and Benjamini-Hochberg Q values across the complete signature-programme family. The response seed is `20260923`.

Composition sensitivity adds the mean `log1p` counts-per-10,000 score of fixed fibroblast markers DCN, LUM, PDGFRA, DPT, and CFD; COL3A1 is excluded through the Matrix programme definition. The third adjustment also adds the supplied Endothelial programme score. Composition markers must be disjoint from all response genes. Patient means and bootstrap inference are recalculated, using seed `20260926` and separate Benjamini-Hochberg families per adjustment. Marker coverage and feature overlaps are written with the results.

Perturbation comparisons first calculate each lesion cluster's median response score. Within section, the mean of the perturbed clusters is compared with the mean of KP reference clusters for Tgfbr2, Ifngr2, and Jak2.

Paired spatial transport uses 16 coordinate-defined KMeans metacells per section, seed `20260923`, and `n_init = 10`. Fitting runs with two OpenMP threads and one BLAS thread. Node scores and positions are spot means; node mass is its fraction of all section spots. Cost is squared Euclidean distance over CRMP_Replication, CRMP_Redox, CRMP_Matrix, Hypoxia, Immune, and Endothelial, scaled by the pooled spot standard deviations, divided by six and by its median. Log-space Sinkhorn uses `epsilon = 0.1`, at most 2,000 iterations, and an update tolerance of `1e-11`. Every Primary/Relapse section pair within a patient is evaluated. Paired programme differences are target minus source spot means. The output records transport cost, marginal error, spot counts, metacell summaries and the full transport matrix.

Spatial outputs use CSV and JSON, including per-section signature scores, section and patient correlations, primary and composition inference, perturbation cluster scores and contrasts, signature coverage, `spatial_metacells.csv`, `paired_spatial_transport.csv`, and `transport_<primary>_<relapse>.csv`. The `summary.json` records numerical settings and counts.

## Verification

```sh
python -m unittest discover -s tests -p test_context_spatial.py
```

The tests use synthetic arrays to check full-library ranks, outer and inner whole-model exclusion, model-balanced loss, the ten estimator outputs, gene-disjoint weighted scores, patient averaging, Sinkhorn marginals, metacell mass and paired spot-mean differences.
