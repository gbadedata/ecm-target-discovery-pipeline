# FAIR Data Principles Alignment

This project follows the FAIR principles (Findable, Accessible, Interoperable, Reusable) as far as practical for a computational analysis of public data.

## Findable

- All data sources are documented in `metadata/data_sources.md` with URLs and citations.
- The sample manifest (`metadata/sample_manifest.csv`) provides a structured index of all samples.
- The repository is publicly accessible on GitHub with a descriptive README.

## Accessible

- All input data is from public repositories (CPTAC, GDC, MSigDB).
- No authentication is required beyond standard data portal access.
- Analysis outputs are generated in open formats (CSV, TSV, PDF, PNG).

## Interoperable

- Gene identifiers use standard Ensembl IDs and HGNC symbols.
- Clinical metadata follows CPTAC's published data dictionary.
- File formats are standard and tool-agnostic (CSV, Parquet, PNG, PDF).

## Reusable

- All analysis parameters are externalised in `params.yaml`.
- Environment files (`environment.yml`, `install_r_packages.R`) enable environment recreation.
- The decision log documents why each choice was made.
- Limitations are documented explicitly.
- The workflow can be adapted to other CPTAC cohorts or similar datasets by changing `params.yaml`.
