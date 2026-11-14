# External data sources

The methods accept externally supplied, schema-compatible numerical inputs. Public source resources relevant to the four main figures are listed below.

| Resource | Identifier or release | Study input |
|---|---|---|
| GEO | [GSE152048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152048) | Discovery osteosarcoma single-cell expression and donor pseudobulks |
| GEO | [GSE162454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162454) | Independent donor expression and pseudobulks |
| GDC | [TARGET-OS](https://portal.gdc.cancer.gov/projects/TARGET-OS) | Patient tumour expression for fixed programme scores |
| DepMap | [Public data portal](https://depmap.org/portal/data_page/) — 24Q4 | Osteosarcoma Chronos effects and model annotations |
| Figshare | [29472362, version 1](https://doi.org/10.6084/m9.figshare.29472362.v1) | Minipool fitness screens, test-gene library, culture metadata, and paired contexts |
| Figshare | [28339340, version 1](https://doi.org/10.6084/m9.figshare.28339340.v1) | Published organoid dependency and expression resources |
| GEO | [GSE193460](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193460) | Spatial CRISPR reference sections, annotations, and measured positions |
| GEO | [GSE293065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293065) | Osteosarcoma spatial sections |

Prepared inputs retain complete eligible observations, stable donor/model identifiers, programme definitions, gene order, and spatial section identifiers. Projection inputs include the programme baseline weights and rank-feature specification; spatial inputs include signature weights and the feature/programme exclusion sets. The corresponding method pages define the complete input interfaces.

Derived figure tables use the formats in [Panel inputs](figure_inputs.md). Figures read the supplied panel values and metadata. The release citation identifies the software; data accessions and the study's derived-data archive identify the numerical inputs.
