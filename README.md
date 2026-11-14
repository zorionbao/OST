# Osteosarcoma context framework

The software implements pure-state-preserving projection, exact ancestry-family transport, paired context transfer, and spatial response analysis.

## Installation

Use Python **3.12.10** and the pinned numerical dependencies. The release has been tested on Windows. File interfaces use Python's cross-platform path handling.

```shell
python -m venv .venv
```

Activate the environment with `.venv\Scripts\Activate.ps1` in PowerShell or `source .venv/bin/activate` in a POSIX shell, then install:

```shell
python -m pip install .
```

Install **Liberation Sans** for figure rendering, or pass `--font-dir` with a directory containing its TTF files. Figures use the Okabe–Ito palette, embedded PDF fonts, editable SVG text, and 350 dpi PNG previews.

## Core analyses

Supply external CSV/NPZ inputs as described in the method documentation. Analyses write to the requested output directory.

```shell
python -m osteosarcoma_core.projection --input-dir inputs/projection --output-dir outputs/projection
python -m osteosarcoma_core.transport --input-dir inputs/transport --output-dir outputs/transport
python -m osteosarcoma_core.context_transfer --input-dir inputs/context --output-dir outputs/context
python -m osteosarcoma_core.spatial_context --input-dir inputs/spatial --output-dir outputs/spatial
```

- [Projection and transport inputs, methods, and outputs](docs/projection_transport.md)
- [Context transfer and spatial analysis inputs, methods, and outputs](docs/context_spatial.md)

The numerical modules can also be imported and used as functions. Importing a module performs no data access or output creation. Donor and model holdouts, fitting penalties, family capacities, and spatial analysis settings are specified by the respective modules.

## Main figures

Prepare panel inputs from the external analysis tables, then render individual panels or assembled figures:

```shell
python -m osteosarcoma_core.prepare_figure_inputs --input-dir inputs/figure_sources --output-dir outputs/panel_inputs
python figures/figure_1a.py --data-dir outputs/panel_inputs --output-dir outputs/panels
python compose_figures.py --data-dir outputs/panel_inputs --output-dir outputs/main_figures
```

Scripts are named `figure_1a.py` through `figure_4h.py`. Their `draw(ax, inputs)` functions render prepared inputs. Statistical inference, rank transformations, PCA, clustering, and model error aggregation are performed during analysis or input preparation. Full heatmap identities and display order accompany each prepared matrix.

- [Panel inputs and rendering metadata](docs/figure_inputs.md)
- [Figure 1 legend](docs/figure_1_legend.md)
- [Figure 2 legend](docs/figure_2_legend.md)
- [Figure 3 legend](docs/figure_3_legend.md)
- [Figure 4 legend](docs/figure_4_legend.md)

## Tests

```shell
python -m unittest discover -s tests -v
```

Tests construct synthetic arrays during execution and examine projection equality, record-replication invariance, optimal transport, grouped holdouts, gene-disjoint scoring, patient aggregation, and coupling marginals.

## External inputs and sources

Research data, fitted coefficients, analysis results, and figure inputs are supplied externally. Input schemas are documented above. Numerical analyses create their output tables from those inputs; the repository distributes source code and documentation.

Relevant public resources are listed in [External data sources](docs/data_sources.md). Source accessions identify the expression, spatial, dependency-screen, and experimental-context resources used by the study. Provider-specific data access and licensing terms accompany those resources.

## Software citation and archiving

The release version is **1.0.0**. Citation metadata are in `CITATION.cff`.

## License

The software is distributed under the [MIT License](LICENSE).
