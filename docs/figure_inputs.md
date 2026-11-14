# Main-figure input specification

The main figures have eight panels each. Every panel has its own drawing module, `figures/figure_1a.py` through `figures/figure_4h.py`. Each module exports `draw(ax, inputs)` and accepts `--data-dir`, `--output-dir` and optional `--font-dir`. `compose_figures.py` assembles those same drawings into four 170 × 185 mm figures. Typography uses Liberation Sans, colors use the Okabe–Ito palette, vector fonts are embedded, and PNG export uses 350 dpi.

## Preparing external figure inputs

Supply analysis results to the input preparation command:

```text
python -m osteosarcoma_core.prepare_figure_inputs --input-dir analysis_results --output-dir panel_inputs
python figures/figure_1a.py --data-dir panel_inputs --output-dir figures_output
python compose_figures.py --data-dir panel_inputs --output-dir figures_output
```

The preparer creates a CSV and JSON descriptor for each panel, plus auxiliary CSVs for the raincloud and box summaries. Each descriptor identifies the datasets, column names, matrix indexing, labels, ordering and display scales. `panel_manifest.json` indexes all 32 panels. These files are written to the supplied output directory. The drawing modules load these prepared tables directly. Statistical inference, parameter selection, PCA, clustering, RMSE calculation and donor summaries belong to analysis and input preparation.

The input directory contains the following CSV basenames. `heldout_predictions.csv.gz` is also accepted. Indexed tables have their identifier in the first CSV column.

| Analysis component | Basenames | Required fields |
| --- | --- | --- |
| Programme representation | `pseudobulk_programmes`, `patient_scores`, `cell_model_scores` | Compartment table: `cohort, donor, compartment, Replication, Redox, Matrix`; patient indexed table: `CRMP_Matrix, CRMP_Replication`; model indexed table: `Matrix, Replication` |
| Pure-state projection | `purity_metrics`, `donor_metrics`, `pure_predictions`, `pure_calibration`, `crmp_summary`, `inner_parameter_surface` | `cohort, programme, method` plus respectively `rna_fraction, absolute_error`; `donor, rmse`; `donor, truth, prediction`; `rmse, sd_ratio`; `mean_rmse, relative_reduction, ci_low, ci_high`; `outer_fold, lambda, alpha, inner_donor_mse` |
| Dependency and ancestry | `context_specific_screen`, `osteosarcoma_gene_effects`, `patient_family_transport`, `transport_sensitivity`, `duplication_benchmark` | Gene-indexed screen: `named_gene, lung_hit, culture_hit, log2FC_culture, log2FC_lung`; model-indexed gene-effect matrix; patient-indexed family coupling matrix; `rejection_cost, represented_mass`; `method, copies, models, median_seconds, duplicated_family_mass, kernel_family_mass` |
| Context benchmark | `screen_rank_matrix`, `library_composition`, `screen_metadata`, `paired_gene_effects`, `heldout_predictions`, `pathway_shifts` | Screen-indexed fitness-rank matrix; gene-indexed library with `All test genes`; screen-indexed metadata with `Growth Format`; `gene, contrast, mean_rank_shift`; `pair_id, model, gene, method, observed, predicted`; `pathway, contrast, model, rank_shift` |
| Context atlas | `osteosarcoma_context_gene_atlas` | Gene-indexed table: `OS_Chronos_mean, MeanDifference, Format` |
| Tissue spatial context | `kp_spots`, `pt_spots`, `signature_inference`, `perturb_section_contrasts`, `spatial_coupling`, `paired_spatial_transport` | Spot tables: `array_x, array_y, Fibroblast / TGFB1`, plus `perturbation` for KP and `CRMP_Matrix` for PT; `signature, programme, median_rho`; `sample, perturbation, signature, difference`; numeric 16 × 16 coupling CSV; `patient, CRMP_Replication, CRMP_Redox, CRMP_Matrix, Hypoxia, Immune, Endothelial` |

Aliases accepted by the preparer are `library.csv` for `library_composition.csv`, `model_scores.csv` for `cell_model_scores.csv`, `metadata.csv` for `screen_metadata.csv`, and `patient_programme_scores.csv` for `patient_scores.csv`. Library field `all_test_genes` and screen field `growth_format` are accepted alongside the publication field names. Export indexed tables with their index as the first column. Context-model output supplies predictions, paired effects, pathway shifts and screen ranks; the source screen metadata and library annotation provide culture labels and test-gene eligibility for the corresponding display.

## Panel-ready tables

Every row below lists the exact CSV payloads and JSON descriptor loaded by its drawing module. `datasets` in the descriptor records each file, its columns, record count and indexed status. An indexed matrix uses its first CSV column as the row identifier. Tables retain all prepared rows and columns. Cartesian axes increase rightwards and upwards; image-style heatmaps put the first row at the top and the first column on the left. Spatial maps use equal aspect and the supplied `x = array_x`, `y = -array_y` coordinates.

Colors use blue `#0072B2`, vermilion `#D55E00`, green `#009E73`, purple `#CC79A7`, light blue `#56B4E9` and black `#000000`. In the tables, **S** means the sequential white–light-blue–blue map; **D** means the signed blue–white–vermilion map. Programme colors are Replication, blue; Redox, green; Matrix, vermilion. Projection method colors are Pathway ranks, black; CRMP fixed, vermilion; Pure-state CRMP, blue. Categorical colors have no numeric color limits.

### Figure 1

| Panel | Prepared filenames | Fields and units | Ordering and coordinate mapping | Color encoding and limits |
| --- | --- | --- | --- | --- |
| 1a | `figure_1a.csv`; `figure_1a.json` | `compartment, programme, value, area`; `value`: mean dimensionless programme rank; `area`: marker area in pt² | x follows descriptor `programmes`; y follows `compartments`, first at top | S with descriptor `vmin = min(value)`, `vmax = max(value)`; `area = 15 + 90 × (value − vmin)/(vmax − vmin)`, denominator 1 for a constant range |
| 1b | `figure_1b.csv`; `figure_1b.json` | `method, rna_fraction, value`; x: RNA fraction; y: mean absolute rank-score error | Methods follow `method_order`; each curve sorts `rna_fraction` ascending | Projection method colors |
| 1c | `figure_1c.csv`; `figure_1c_density.csv`; `figure_1c_quantiles.csv`; `figure_1c.json` | Main: `donor, programme, value` (% RMSE reduction); density: `programme, x, height` (x in %, height in categorical-axis units); quantiles: `programme, q1, median, q3` (%) | `programme_order` maps bottom to top; points retain donor order within programme; density x increases left to right; density height is normalized to maximum 0.33 | Programme colors |
| 1d | `figure_1d.csv`; `figure_1d.json` | `donor, method, truth, prediction`; x/y: dimensionless pure matrix rank | Source row order within fixed CRMP and Pure-state CRMP groups; x = `truth`, y = `prediction`; identity segment uses `identity_min = min(truth) − 0.025`, `identity_max = max(truth) + 0.025` | Fixed CRMP: vermilion squares; Pure-state CRMP: blue circles; identity: black dotted line |
| 1e | `figure_1e.csv`; `figure_1e_boxes.csv`; `figure_1e.json` | Main: `donor, programme, value`; boxes: `programme, med, q1, q3, whislo, whishi`; all values: signed dimensionless prediction-minus-truth error | x follows `programme_order` left to right; donor points retain source order; box summaries and 1.5-IQR whisker endpoints are supplied in the auxiliary file | Programme colors; black medians and whiskers |
| 1f | `figure_1f.csv`; `figure_1f.json` | `method, metric, value`; dimensionless within-metric relative values | x follows `metric_order`: Mixture RMSE, Pure RMSE, Score SD; curves follow `method_order`; `normalization` supplies the two displayed RMSE maxima and SD divisor 1 | Projection method colors |
| 1g | `figure_1g.csv`; `figure_1g.json` | `programme, relative_reduction, ci_low, ci_high`; estimate and interval endpoints in % | y follows `programme_order` bottom to top; x = reduction; horizontal intervals use supplied endpoints | Programme colors |
| 1h | `figure_1h.csv`; `figure_1h.json` | `x, lambda, value`; x = log10 coefficient penalty; y = dimensionless drift penalty; `value` = mean inner MSE × 1,000 | Grid is pivoted with ascending `lambda` rows and ascending x columns; contour coordinates increase rightwards/upwards; descriptor `x_transform`, `value_unit` and `contour_levels` preserve transforms | S; `contourf` uses `levels = contour_levels` (7); Matplotlib generates level boundaries from the supplied grid range, and first/last boundaries set the color limits; evaluated grid points are black |

### Figure 2

| Panel | Prepared filenames | Fields and units | Ordering and coordinate mapping | Color encoding and limits |
| --- | --- | --- | --- | --- |
| 2a | `figure_2a.csv`; `figure_2a.json` | `gene, x, y, category`; x = culture log2 fold change; y = lung log2 fold change | Source gene order within categories; drawing order: Both contexts, Lung threshold only, Other genes | Both contexts: blue; Lung threshold only: vermilion; Other genes: green |
| 2b | `figure_2b.csv`; `figure_2b.json` | `category, x, y`; x: depletion log2 fold change; y: cumulative gene fraction | Lung then Culture; each group's x is ascending and y = rank/N | Lung: vermilion; Culture: blue |
| 2c | `figure_2c.csv`; `figure_2c.json` | `model, gene, value`; x: dimensionless Chronos gene effect | y follows `gene_order` top to bottom, selected by ascending lung log2 fold change; model records retain source order within each gene | Blue points |
| 2d | `figure_2d.csv`; `figure_2d.json` | `identifier, x, y, category`; x/y: dimensionless matrix/replication programme scores | Patient rows then cell-model rows, preserving each source order; x = matrix, y = replication | Patient tumours: blue circles; Cell models: vermilion triangles |
| 2e | `figure_2e.csv`; `figure_2e.json` | Indexed patient × family matrix; cell values are fractions of each patient's input mass | `patient_order` top to bottom is descending allocated mass; `family_order` left to right retains positive-allocation source columns; `family_display_labels` supplies full label wrapping; `patient_input_mass = 1/N` | S, fixed limits 0–1 |
| 2f | `figure_2f.csv`; `figure_2f.json` | `rejection_cost, represented_percent`; x: dimensionless unmatched-mass cost; y: represented patient mass (%) | Ascending `rejection_cost` left to right | Blue line and circles |
| 2g | `figure_2g.csv`; `figure_2g.json` | `copies, duplicated_family_mass, kernel_family_mass`; x: record count; y: allocated mass fraction | Ascending copies; logarithmic x; linear y | Ancestry quotient: blue circles; Row kernel: vermilion squares |
| 2h | `figure_2h.csv`; `figure_2h.json` | `method, models, milliseconds`; x: total record count; y: median solve time (ms) | Full sparse then Ancestry quotient; ascending model count within method; logarithmic x and y | Full sparse: vermilion squares; Ancestry quotient: blue circles |

### Figure 3

| Panel | Prepared filenames | Fields and units | Ordering and coordinate mapping | Color encoding and limits |
| --- | --- | --- | --- | --- |
| 3a | `figure_3a.csv`; `figure_3a.json` | Indexed screen × test-gene matrix; cells: dimensionless within-screen fitness rank | `screen_order` top to bottom; `gene_order` left to right is average-linkage Euclidean order from complete eligible profiles; `clustering` records this rule | S, fixed limits 0–1 |
| 3b | `figure_3b.csv`; `figure_3b.json` | `screen, x, y, format`; x/y: dimensionless fitness-rank PCA coordinates | Rows follow `screen_order`; x = fitted PC1, y = PC2; `feature_order`, `pca_mean`, `pca_components`, `explained_variance_ratio` preserve projection and axis orientation | Plastic: blue circles; Dome: vermilion triangles |
| 3c | `figure_3c.csv`; `figure_3c.json` | `gene, x, y`; dimensionless rank shifts: x = Dome − Plastic; y = Serum-free − Serum | Pivot-generated ascending gene order; conventional Cartesian directions | Blue points |
| 3d | `figure_3d.csv`; `figure_3d.json` | `model, method, rmse`; y: dimensionless fitness-rank RMSE | Traces follow `model_order`; x follows `method_order` and is labelled M1–M8; the RMSE combines all eligible measurements of each model/method | Blue model traces; PCHT highlighted vermilion; Global anchor highlighted green |
| 3e | `figure_3e.csv`; `figure_3e.json` | Indexed model × method matrix; cells: dimensionless fitness-rank RMSE | `model_order` top to bottom; `method_order` left to right, displayed as M1–M8; model tick labels remove the `ACH-` prefix while descriptor identities retain it | S; lower limit `min(0.1, minimum matrix RMSE)`, upper limit `max(0.4, maximum matrix RMSE)` |
| 3f | `figure_3f.csv`; `figure_3f.json` | `pair_id, model, gene, observed, predicted`; x/y: dimensionless fitness ranks | PCHT rows retain source prediction order; x = observed, y = predicted; descriptor `hexbin_gridsize = 26` | S applied to measurements per nonempty hexagon; `mincnt = 1`; automatic limits are minimum/maximum displayed bin counts |
| 3g | `figure_3g.csv`; `figure_3g.json` | Indexed pathway × contrast matrix; cells: dimensionless mean paired rank shift | `pathway_order` top to bottom; `contrast_order` left to right; all pathway ticks are shown; PI3K–MTOR label retains a complete single-line display | D; symmetric limits ±`max(0.12, maximum absolute matrix value)` |
| 3h | `figure_3h.csv`; `figure_3h.json` | `gene, x, y, color`; x: osteosarcoma mean Chronos; y: organoid-minus-adherent Chronos; color: dimensionless paired format rank shift | Eligible atlas rows retain source order; x/y use conventional Cartesian directions | D; limits ±descriptor `color_limit`, calculated as maximum absolute `color` |

M1–M8 correspond to No change, Target mean, Paired mean, Unpaired mean, Nearest pair, Global anchor, Paired midpoint and PCHT. Their full names are retained in both `method_order` descriptors.

### Figure 4

| Panel | Prepared filenames | Fields and units | Ordering and coordinate mapping | Color encoding and limits |
| --- | --- | --- | --- | --- |
| 4a | `figure_4a.csv`; `figure_4a.json` | `x, y, category`; coordinates: spatial array units | Source spot order within descriptor `category_order`; x = array_x, y = −array_y; equal aspect; `vertical_coordinate` records direction | `category_order` KP, Tgfbr2, Jak2, periphery, Unassigned maps to blue, vermilion, green, purple, light blue |
| 4b | `figure_4b.csv`; `figure_4b.json` | `x, y, value`; coordinates: spatial array units; value: dimensionless weighted normalized-expression TGFB1-response score in KP_1 | Source spot order; x = array_x, y = −array_y; equal aspect; descriptor `score = Fibroblast / TGFB1` and `vertical_coordinate` | S; descriptor `vmin = min(value)`, `vmax = max(value)` for this panel |
| 4c | `figure_4c.csv`; `figure_4c.json` | `x, y, value`; coordinates: spatial array units; value: dimensionless weighted normalized-expression TGFB1-response score in PT1 | Source spot order; x = array_x, y = −array_y; equal aspect; descriptor `score = Fibroblast / TGFB1` and `vertical_coordinate` | S; descriptor `vmin = min(value)`, `vmax = max(value)` for this panel |
| 4d | `figure_4d.csv`; `figure_4d.json` | `x, y, value`; coordinates: spatial array units; value: dimensionless fixed CRMP matrix rank in PT1 | Source spot order; x = array_x, y = −array_y; equal aspect; descriptor `score = CRMP_Matrix` and `vertical_coordinate` | S; descriptor `vmin = min(value)`, `vmax = max(value)` for this panel |
| 4e | `figure_4e.csv`; `figure_4e.json` | Indexed signature × programme matrix; cells: dimensionless patient-median partial Spearman correlations | `signature_order` top to bottom; `programme_order` Matrix, Redox, Hypoxia left to right; signature labels replace Lung_tumor with Tumour and remove the slash separator | D; limits ±descriptor `color_limit = max(0.5, maximum absolute matrix value)` |
| 4f | `figure_4f.csv`; `figure_4f.json` | `sample, perturbation, difference`; y: dimensionless TGFB1-response cluster-minus-KP score difference | x follows descriptor `perturbation_order` (ascending category names); sections retain source order within each category; positive differences increase upwards | Tgfbr2: vermilion; Jak2: green; Ifngr2: purple |
| 4g | `figure_4g.csv`; `figure_4g.json` | Indexed primary × relapse metacell matrix; cells: transported spot-mass fractions | `primary_order` top to bottom; `relapse_order` left to right; row/column identities are the supplied metacell order; `mass_sum` records total coupling mass | S; lower limit 0, upper limit descriptor `vmax = maximum coupling entry` |
| 4h | `figure_4h.csv`; `figure_4h.json` | Indexed patient × programme matrix; cells: dimensionless original-score relapse-minus-primary changes | `patient_order` top to bottom; `programme_order` Replication, Redox, Matrix, Hypoxia, Immune, Endothelial left to right | D; limits ±descriptor `color_limit = max(0.16, maximum absolute matrix value)` |

The full model RMSE in 3d/e is computed as the square root of the mean squared error over all eligible measurements of each model and method. The context-specific RMSE summaries are separate quantities. The gene order and PCA rotation in 3a/b are generated together from the same eligible rank matrix. Spatial vertical coordinates are negated array y; tissue maps use equal aspect. The ancestry heatmap uses patient input mass 1/N, while the primary–relapse heatmap retains spot-mass fractions with total mass one.

Single-panel exports display the complete labels. In the compact composition, dense heatmaps retain every matrix cell and use spaced tick labels with complete strings. The JSON descriptors and figure legends retain the complete ordered identities, including the TGFB1 signatures and all metacells. All nine pathways are labelled in 3g. Legend files are `figure_1_legend.md` through `figure_4_legend.md`. Each export also writes a layout report containing its geometry, typography and any text extending outside the canvas.
