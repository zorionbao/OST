# Projection and ancestry transport

These modules implement Counterfactual Rank Mixture Projection (CRMP), Pure-state CRMP, and ancestry-quotient partial transport. Numerical interfaces accept arrays and tables. File loading occurs only through the command-line entry points.

## Projection inputs

The input directory contains three files:

| File | Required contents |
| --- | --- |
| `projection_arrays.npz` | Numeric arrays `X`, `Y`, `X0`, `W`, `PX`, and `PY` |
| `mixtures.csv` | One row per mixture, with `cohort`, `donor`, `reference`, and `rna_fraction` |
| `pure_profiles.csv` | One row per pure profile, with `cohort` and `donor` |

`X` and `X0` have shape `(mixtures, features)`. `Y` has shape `(mixtures, programmes)`. `W` has shape `(features, programmes)`. `PX` and `PY` have shapes `(pure_profiles, features)` and `(pure_profiles, programmes)`. Array rows follow their corresponding CSV row order. Programme columns default to Replication, Redox, and Matrix in that order. Supply one pure profile per cohort and donor. The discovery cohort defaults to `GSE152048`; every other cohort supplies independent evaluation. At least three discovery donors and one independent donor are required for grouped nested validation.

```bash
python -m osteosarcoma_core.projection --input-dir inputs/projection --output-dir outputs/projection
```

Use `--discovery-cohort` for a different discovery label and `--programmes` for alternative target column names. `--seed` defaults to `20261002`; `--bootstrap-resamples` defaults to `10000`.

`rank_features(expression, gene_universe, feature_genes)` accepts a gene-by-profile DataFrame and returns profile-by-feature average ranks divided by the gene-universe size. Missing genes receive zero expression before ranking. `programme_baseline` creates feature coefficients by averaging matched genes within each programme. `build_mixtures` accepts a mapping from `(cohort, donor)` to gene-by-compartment pseudobulks. It normalizes each available compartment by its total expression, pairs Osteogenic with each available background within the donor, and evaluates RNA fractions 0.2, 0.4, 0.6, 0.8, and 1.0. It returns the six input arrays and the two metadata tables.

## Projection estimates and outputs

`fit_projection` returns feature coefficients, an intercept, and the maximum training pure-state equality residual. Balanced fitting assigns each row weight `1 / (D * K_d * M_dk)`, where `D` is the donor count, `K_d` the donor's available background count, and `M_dk` the background's record count. Pure-state coefficients equal baseline coefficients plus a correction in the null space of the augmented pure-feature matrix. A full-rank pure constraint returns the baseline coefficients and a zero intercept.

`run_projection` preserves five fitted methods: CRMP fixed, CRMP nested, Balanced CRMP, Pure-state CRMP, and Pure-state no drift. Comparators include pathway ranks, scalar calibration, and the training mean. Fixed CRMP uses drift penalty 1 and coefficient penalty 0.1. Nested fits search drift penalties 0, 0.25, 1, and 4, and coefficient penalties 0.01, 0.1, and 1. The no-drift fit fixes the drift penalty at zero. Every outer discovery fold holds out one complete donor; every inner fold also holds out one complete donor. Inner risks average programmes and RNA records within backgrounds, then backgrounds within donors, then donors. Independent-cohort targets enter evaluation after selection.

The command writes:

- `heldout_predictions.csv.gz`, `pure_predictions.csv`: individual predictions and targets.
- `donor_metrics.csv`, `background_metrics.csv`, `purity_metrics.csv`: grouped prediction errors.
- `pure_calibration.csv`: slope, intercept, standard-deviation ratio, RMSE, and Spearman correlation.
- `fold_parameters.csv`, `inner_parameter_surface.csv`, `fold_assignments.csv`: selected penalties, inner losses, and donor membership.
- `mixture_design.csv`, `summary.csv`: mixture counts and donor-level error reductions, paired bootstrap intervals, and Wilcoxon tests.
- `run_summary.json`: programmes, seeds, grid, donor identifiers, and fitted coefficient matrices.

Evaluation takes the square root after averaging squared error within each donor, and summaries weight donors equally. Pure calibration uses one profile per donor. Bootstrap intervals resample complete donors 10000 times.

## Transport inputs

Choose either score tables or a supplied cost matrix. The input mode is determined by the presence of `cost.npz`.

| Mode | Required files and columns |
| --- | --- |
| State scores | `patient_scores.csv` and `model_scores.csv`; first column is the patient or model identifier; programme columns are Replication, Redox, Matrix or their `CRMP_`-prefixed forms; models also require `family` |
| Cost matrix | `cost.npz` with numeric array `cost`; `families.csv` with one `family` per cost column and optional `model_id`; the NPZ may include `patient_ids` as a string array |

An optional `capacities.csv` contains `family` and `capacity`. Each family receives capacity `1 / number_of_families` when this file is absent. State-score mode uses `CRMP_`-prefixed columns when both prefixed and unprefixed scores are present. It standardizes patient and model coordinates using patient means and sample standard deviations, then averages squared Euclidean distances over the selected programmes.

```bash
python -m osteosarcoma_core.transport --input-dir inputs/transport --output-dir outputs/transport
python -m osteosarcoma_core.transport --input-dir inputs/transport --output-dir outputs/transport --benchmark
```

`--kappa` defaults to 0.9. `--rejection-costs` defaults to 0.25, 0.5, 0.9, 1.5, and 3. `--programmes` changes the selected score columns. The duplication benchmark runs only with `--benchmark`, using total copies 1, 2, 5, 20, 100, 500, and 2000 and three timed solves per configuration. `--copies`, `--repeats`, and `--seed` set the benchmark configuration. The default random seed is `20261002`.

## Transport estimates and outputs

The full linear program assigns patient mass `1 / number_of_patients`, conserves each patient's allocated plus unmatched mass, and constrains total family mass. Descendants share family capacity. The objective adds transport costs and `kappa * unmatched_mass`. The quotient replaces each patient-family cost with its minimum descendant cost. A quotient solution lifts to the first minimum-cost descendant in input order. HiGHS solves both formulations independently.

`solve_cost` returns the patient-model coupling, unmatched mass, objective, and sorted families. `solve_quotient` returns the family coupling, descendant lift, unmatched mass, objective, quotient cost, and sorted families. Capacity arrays follow that sorted family order.

The command writes `patient_model_cost.csv`, `full_patient_model_transport.csv`, `patient_family_transport.csv`, `lifted_patient_model_transport.csv`, `quotient_cost.csv`, `patient_representation.csv`, `family_capacities.csv`, `transport_sensitivity.csv`, and `run_summary.json`. Matrix rows retain patient identifiers. The optional benchmark adds `duplication_benchmark.csv` and `randomized_equivalence.csv`. Timing includes sparse constraint construction and optimization; quotient timing also includes reduction and lifting. The common precomputed patient-model cost is outside the timer. Randomized equivalence evaluates 100 general-cost problems with varying dimensions, descendant counts, capacities, and unmatched-mass costs.

## Tests

```bash
python -m unittest discover -s tests -p test_projection_transport.py
```

Tests construct arrays during execution and check pure-state equality, complete donor/background record replication, grouped target separation, full-rank constraints, rank-mixture construction, full/quotient equality, lift feasibility, tie ordering, zero capacity, patient standardization, and 100 randomized transport cases.
