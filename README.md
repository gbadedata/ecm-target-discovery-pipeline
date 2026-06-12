# ECM Target Discovery Pipeline

A reproducible multi-omics bioinformatics workflow investigating extracellular matrix-associated molecular signatures in pancreatic ductal adenocarcinoma (PDAC), integrating bulk transcriptomics, quantitative proteomics, and clinical metadata to identify candidate targets relevant to fibrosis and oncology research.

## The Problem

In pancreatic cancer, the extracellular matrix is not incidental to the disease -- it defines it. The desmoplastic stroma in PDAC influences immune exclusion, drug delivery, tumour invasion, and patient survival. Drug discovery teams need reproducible computational workflows to move from complex multi-omics data to evidence-based candidate targets. This pipeline provides that workflow.

## Key Results

| Analysis | Finding |
|---|---|
| Differential expression | 2,866 significant genes (1,385 up, 1,481 down); 306 ECM genes significantly altered |
| ECM signature validation | Computational ECM score correlates with pathologist stromal fraction (Spearman rho=0.39, p=1.4e-6) |
| Pathway enrichment | 15/31 ECM pathways significantly enriched; all upregulated in tumour |
| Proteomics | ECM remodelling confirmed at protein level; TGF-beta axis components identified |
| RNA-protein concordance | 111 concordant ECM targets; ECM rho=0.60 (stronger than genome average of 0.47) |
| Survival | ECM score does not predict survival (p=0.83), consistent with published PDAC literature |
| ML classification | Random Forest AUC=0.87 for ECM-high vs ECM-low prediction |
| Target prioritisation | Top 30 candidates all Tier 1 (multi-omics concordant) with transparent scoring |

### Top 10 ECM Target Candidates

| Rank | Gene | RNA log2FC | Protein log2FC | Score | Biological Role |
|---|---|---|---|---|---|
| 1 | SERPINB5 | 5.90 | 1.74 | 0.975 | Serine protease inhibitor |
| 2 | S100A11 | 1.56 | 1.33 | 0.967 | Calcium-binding, ECM remodelling |
| 3 | IL1RN | 2.49 | 1.53 | 0.963 | IL-1 receptor antagonist |
| 4 | INHBA | 1.75 | 1.64 | 0.962 | TGF-beta superfamily (fibrosis driver) |
| 5 | COL11A1 | 2.79 | 1.45 | 0.959 | Collagen XI (known PDAC stromal marker) |
| 6 | S100P | 4.93 | 2.15 | 0.955 | Known PDAC diagnostic biomarker |
| 7 | CTHRC1 | 1.29 | 1.66 | 0.951 | Promotes collagen deposition |
| 8 | COL12A1 | 1.55 | 1.43 | 0.951 | FACIT collagen, stromal component |
| 9 | THBS2 | 1.16 | 1.75 | 0.944 | Thrombospondin-2, ECM glycoprotein |
| 10 | LAMA3 | 2.59 | 1.37 | 0.943 | Laminin alpha-3, basement membrane |

Notable druggable targets: **LOXL2** (rank 13, collagen crosslinker, clinical trial target), **PLAU** (rank 17, urokinase plasminogen activator), **FN1** (rank 24, fibronectin).

## Architecture

```
CPTAC PDAC Data (140 tumour, 21 normal, matched RNA + protein + clinical)
      |
      v
[Data Validation]                Python (Pydantic, pandas)
      |
      v
[Clinical Cleaning]              Python (semicolon parsing, survival mapping)
      |
      +---------------------------+
      |                           |
      v                           v
[Differential Expression]    [ECM Signature Scoring]
R (limma)                    R (GSVA, ssGSEA)
      |                           |
      v                           v
[Pathway Enrichment]         [Survival Analysis]
R (fgsea, msigdbr)           R (survival, survminer, Cox PH)
      |
      v
[Proteomics Analysis]            R (limma, MinProb imputation)
      |
      v
[Multi-Omics Integration]       Python (RNA-protein concordance)
      |
      v
[ML & Clustering]                Python (scikit-learn, PCA, KMeans, RF)
      |
      v
[Target Prioritisation]          Python (multi-evidence weighted scoring)
```

## Dataset

**CPTAC Pancreatic Ductal Adenocarcinoma (PDAC)**

| Layer | Samples | Features | Source |
|---|---|---|---|
| Bulk RNA-seq (RSEM-UQ log2) | 140 tumour, 21 normal | 28,057 genes | LinkedOmics |
| TMT Proteomics | 140 tumour, 75 normal | 11,662 proteins | LinkedOmics |
| Clinical metadata | 140 patients | 39 variables | LinkedOmics |

Citation: Cao L, Huang C, et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell. 2021;184(19):5031-5052.e26.

## Evidence and Reproducibility

Every phase produces structured JSON evidence in `evidence/logs/`. All analysis is reproducible by running the scripts in `src/r/` and `src/python/`.

| Phase | Evidence | Reproducibility Command |
|---|---|---|
| Data validation | `data_validation_report.json` | `python -m src.python.validate_data` |
| Clinical cleaning | `clinical_cleaning_report.json` | `python -m src.python.clean_clinical` |
| Differential expression | `de_analysis_report.json` | `Rscript src/r/01_differential_expression.R` |
| ECM scoring | `ecm_scoring_report.json` | `Rscript src/r/02_ecm_signature_scoring.R` |
| Pathway enrichment | `enrichment_report.json` | `Rscript src/r/03_pathway_enrichment.R` |
| Proteomics | `proteomics_report.json` | `Rscript src/r/04_proteomics_analysis.R` |
| Survival | `survival_report.json` | `Rscript src/r/05_survival_analysis.R` |
| Multi-omics | `multiomics_report.json` | `python -m src.python.multiomics_integration` |
| ML/clustering | `ml_report.json` | `python -m src.python.ml_clustering` |
| Prioritisation | `prioritisation_report.json` | `python -m src.python.target_prioritisation` |

## Testing

```
25 tests passing
CI: lint (ruff) + pytest + R syntax check
```

## Quickstart

```bash
git clone https://github.com/gbadedata/ecm-target-discovery-pipeline.git
cd ecm-target-discovery-pipeline

# Python
python3 -m venv .venv
source .venv/bin/activate
pip install pyyaml pandas numpy scipy scikit-learn statsmodels lifelines matplotlib seaborn pydantic pytest ruff structlog gseapy

# R packages
Rscript install_r_packages.R

# Download data (see metadata/data_sources.md)
# Then run each phase in order, or use Nextflow:
nextflow run main.nf -params-file params.yaml
```

## Key Decisions

| Decision | Choice | Reasoning |
|---|---|---|
| Dataset | CPTAC PDAC | Matched RNA + protein + clinical; strongest ECM narrative |
| DE method | limma (not DESeq2) | Data is pre-normalised RSEM-UQ log2; limma is appropriate for continuous data |
| Signature scoring | ssGSEA (GSVA) | Per-sample continuous scores; standard in translational oncology |
| Proteomics DA | limma + MinProb imputation | Robust for TMT data with 24.6% missingness |
| Survival | Kaplan-Meier + Cox PH | Standard clinical association framework |
| Target scoring | Weighted multi-evidence composite | Transparent, reproducible, interpretable |
| Workflow | Nextflow DSL2 | Standard bioinformatics workflow manager |

Full decision log: [docs/decision-log.md](docs/decision-log.md)

## Challenges and Solutions

| Challenge | Solution |
|---|---|
| Vital status labels ("Deceased"/"Living" not "Dead"/"Alive") | Flexible mapping handling both conventions |
| Tissue composition as semicolon-separated multi-reviewer estimates | Parser averages multiple pathologist estimates, converts to fractions |
| 24.6% proteomics missingness | MinProb imputation with documented assumptions |
| ECM score does not predict PDAC survival | Reported as finding consistent with published literature; not hidden |
| R package BH/C++ compatibility (fgsea) | Set CXX11STD to gnu++14 in Makevars |

## Documentation

- [Project Charter](docs/project-charter.md)
- [Architecture](docs/architecture.md)
- [Decision Log](docs/decision-log.md)
- [Data Sources](metadata/data_sources.md)
- [Clinical Dictionary](metadata/clinical_dictionary.csv)
- [Reproducibility](docs/reproducibility.md)
- [FAIR Principles](docs/fair-principles.md)
- [Limitations](docs/limitations.md)

## Tools

**R:** limma, GSVA, fgsea, msigdbr, survival, survminer, tidyverse, pheatmap, ggrepel, ggplot2

**Python:** pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, Pydantic, pytest, structlog

**Workflow:** Nextflow DSL2, GitHub Actions CI

## Author

**Oluwagbade Odimayo**
Bioinformatics Data Engineer | MSc Applied Data Science

- [GitHub](https://github.com/gbadedata)
- [LinkedIn](https://www.linkedin.com/in/oluwagbade-odimayo-)
