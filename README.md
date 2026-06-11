# ECM Target Discovery Pipeline

A reproducible multi-omics bioinformatics workflow for investigating extracellular matrix-associated molecular signatures in pancreatic ductal adenocarcinoma, integrating bulk transcriptomics, quantitative proteomics, and clinical metadata to identify candidate targets relevant to fibrosis and oncology research.

## The Problem

Drug discovery teams working in fibrosis and oncology need reproducible ways to move from complex molecular datasets to biological insight. In pancreatic cancer, the extracellular matrix is not incidental to the disease -- it defines it. The desmoplastic stroma in PDAC influences immune exclusion, drug delivery, tumour invasion, and patient survival. Understanding which ECM-associated molecular signatures drive these processes is a central question in translational research.

This project uses public multi-omics data from the CPTAC PDAC cohort to computationally investigate that question.

## What This Project Does

1. **Ingests and validates** matched transcriptomic, proteomic, and clinical data from the CPTAC PDAC cohort
2. **Performs differential expression analysis** (DESeq2) to identify genes associated with tumour biology
3. **Scores ECM signatures** per sample using ssGSEA across curated matrisome and stromal gene sets
4. **Runs pathway enrichment** (fgsea, clusterProfiler) to identify ECM, fibrosis, and TME pathway activity
5. **Analyses proteomics** data (limma) to identify differentially abundant ECM-associated proteins
6. **Integrates RNA and protein** evidence to assess concordance and identify multi-omics-supported signals
7. **Associates ECM signatures with clinical outcomes** using Kaplan-Meier survival analysis and Cox regression
8. **Applies dimensionality reduction, clustering, and machine learning** to identify molecular subgroups
9. **Prioritises candidate ECM targets** using a transparent multi-evidence scoring framework

## Architecture

```
CPTAC PDAC Data
      |
      v
[Data Ingestion & Validation]     Python (Pydantic, pandas)
      |
      v
[Transcriptomics Analysis]        R (DESeq2, GSVA, clusterProfiler)
      |
      v
[Proteomics Analysis]             R (limma) + Python
      |
      v
[Multi-Omics Integration]         Python (pandas, scipy)
      |
      v
[Clinical & Survival Analysis]    R (survival, survminer)
      |
      v
[ML & Clustering]                 Python (scikit-learn, UMAP)
      |
      v
[Target Prioritisation]           Python
      |
      v
[Scientific Report]               Quarto
```

## Dataset

**CPTAC Pancreatic Ductal Adenocarcinoma (PDAC)**

| Layer | Type | Source |
|---|---|---|
| Transcriptomics | Bulk RNA-seq | CPTAC / GDC |
| Proteomics | TMT quantitative | CPTAC Data Portal |
| Clinical | Demographics, staging, survival | CPTAC / GDC |

Citation: Cao L, Huang C, Cui Zhou D, et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell. 2021;184(19):5031-5052.e26.

## Repository Structure

```
ecm-target-discovery-pipeline/
├── README.md
├── params.yaml                  # All configurable parameters
├── environment.yml              # Python conda environment
├── install_r_packages.R         # R package installation
├── main.nf                      # Nextflow workflow
├── nextflow.config
│
├── docs/                        # Project documentation
│   ├── project-charter.md
│   ├── decision-log.md
│   ├── architecture.md
│   ├── reproducibility.md
│   ├── fair-principles.md
│   └── limitations.md
│
├── metadata/                    # Data provenance and dictionaries
│   ├── data_sources.md
│   ├── sample_manifest.csv
│   ├── clinical_dictionary.csv
│   └── checksums.txt
│
├── src/
│   ├── r/                       # R analysis scripts
│   └── python/                  # Python validation, integration, ML
│
├── workflows/modules/           # Nextflow modules
├── notebooks/                   # Exploratory analysis notebooks
├── reports/                     # Quarto/RMarkdown reports
├── figures/                     # Generated figures
├── results/                     # Analysis outputs
├── evidence/                    # Logs, screenshots, evidence index
└── tests/                       # Test suite
```

## Quickstart

### Prerequisites

- R >= 4.4
- Python 3.12 (via conda)
- Nextflow >= 24.0 (optional)

### Setup

```bash
git clone https://github.com/gbadedata/ecm-target-discovery-pipeline.git
cd ecm-target-discovery-pipeline

# Python
conda env create -f environment.yml
conda activate ecm-pipeline

# R
Rscript install_r_packages.R
```

### Run

```bash
# Full workflow
nextflow run main.nf -params-file params.yaml

# Or step by step (see docs/reproducibility.md)
```

## Key Decisions

| Decision | Choice | Reasoning |
|---|---|---|
| Dataset | CPTAC PDAC | Matched RNA + protein + clinical; strongest ECM narrative |
| DE method | DESeq2 | Gold standard for bulk RNA-seq differential expression |
| Signature scoring | ssGSEA (GSVA) | Per-sample continuous scores; standard in translational oncology |
| Proteomics DA | limma | Robust for TMT data with missing values |
| Survival | Kaplan-Meier + Cox | Standard clinical association framework |
| ML | Logistic regression, random forest, elastic net | Interpretable models suitable for biological feature ranking |
| Workflow | Nextflow DSL2 | Standard bioinformatics workflow manager; R + Python in one pipeline |

Full decision log: [docs/decision-log.md](docs/decision-log.md)

## Documentation

- [Project Charter](docs/project-charter.md)
- [Architecture](docs/architecture.md)
- [Decision Log](docs/decision-log.md)
- [Data Sources](metadata/data_sources.md)
- [Reproducibility](docs/reproducibility.md)
- [FAIR Principles](docs/fair-principles.md)
- [Limitations](docs/limitations.md)

## Tools and Technologies

**R:** DESeq2, edgeR, limma, GSVA, fgsea, clusterProfiler, msigdbr, survival, survminer, pheatmap, ggplot2, tidyverse

**Python:** pandas, numpy, scipy, scikit-learn, lifelines, statsmodels, seaborn, matplotlib, Pydantic, pytest, structlog

**Workflow:** Nextflow DSL2, conda, GitHub Actions

## Status

| Phase | Description | Status |
|---|---|---|
| 0 | Repository scaffold and environments | Complete |
| 1 | Dataset selection and data source documentation | In progress |
| 2 | Data ingestion and validation | Pending |
| 3 | Clinical metadata cleaning | Pending |
| 4 | Transcriptomics analysis | Pending |
| 5 | ECM signature scoring | Pending |
| 6 | Pathway enrichment | Pending |
| 7 | Proteomics analysis | Pending |
| 8 | Multi-omics integration | Pending |
| 9 | Clinical and survival analysis | Pending |
| 10 | ML and clustering | Pending |
| 11 | Target prioritisation | Pending |
| 12 | Nextflow orchestration | Pending |
| 13 | Testing and CI | Pending |
| 14 | Scientific reporting | Pending |

## Author

**Oluwagbade Odimayo**
Bioinformatics Data Engineer | MSc Applied Data Science

- [GitHub](https://github.com/gbadedata)
- [LinkedIn](https://www.linkedin.com/in/oluwagbade-odimayo-)
