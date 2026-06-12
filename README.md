# ECM Target Discovery Pipeline

A reproducible multi-omics bioinformatics workflow for investigating extracellular matrix-associated molecular signatures in pancreatic ductal adenocarcinoma. The pipeline integrates bulk RNA-seq transcriptomics, TMT quantitative proteomics, and clinical metadata from the CPTAC PDAC cohort to identify, validate, and prioritise candidate ECM-associated targets relevant to fibrosis and oncology research.

140 tumour samples. 28,057 genes. 11,662 proteins. 39 clinical variables. One question: which ECM molecules are driving stromal remodelling in pancreatic cancer, and can we identify them with multi-layer evidence?

The answer is a validated, transparent, 30-target shortlist supported by transcriptomic significance, proteomic confirmation, multi-omics concordance, pathway enrichment, and machine learning feature ranking. Every target is Tier 1: concordant across both RNA and protein.

---

## The Problem

Pancreatic ductal adenocarcinoma has a five-year survival rate of approximately 11%. A defining pathological feature is the desmoplastic reaction: cancer-associated fibroblasts deposit massive quantities of extracellular matrix proteins, creating a dense stromal barrier that blocks immune cell infiltration, impedes chemotherapy delivery, and promotes tumour invasion and metastasis. Most of what a pathologist sees under the microscope in a PDAC tumour is not cancer cells. It is stroma.

Drug discovery teams working in fibrosis and oncology need reproducible computational workflows to identify which ECM-associated molecules drive this stromal programme. The challenge is scale and integration. A single cohort generates measurements across tens of thousands of genes, thousands of proteins, and dozens of clinical variables. No spreadsheet can handle that. No single data layer is sufficient. A gene that is overexpressed at the RNA level but absent at the protein level is a weaker candidate than one confirmed in both.

This pipeline solves that problem. It takes matched multi-omics data, applies rigorous statistical and computational methods across six evidence layers, and produces a ranked candidate target table with transparent scoring and documented limitations.

---

## Key Results

| Analysis | Finding |
|---|---|
| Differential expression (limma) | 2,866 genes significantly altered in tumour vs normal (FDR<0.05, \|log2FC\|>1). 1,385 upregulated, 1,481 downregulated. 306 ECM genes significant |
| ECM signature validation | ssGSEA Core Matrisome score correlates with pathologist-estimated stromal fraction: Spearman rho=0.39, p=1.4e-6 (n=140) |
| Pathway enrichment (fgsea) | 15 of 31 ECM-related pathways significantly enriched (FDR<0.05). All upregulated in tumour. Hallmark EMT enriched at NES=2.09, p=1.85e-10 |
| Proteomics differential abundance | ECM remodelling confirmed at protein level. TGF-beta axis components (LTBP1, INHBA) identified. 24.6% missingness handled via MinProb imputation |
| RNA-protein concordance | 111 concordant ECM targets (85 up, 26 down). ECM-specific rho=0.60, stronger than genome-wide average of 0.47 |
| Survival analysis | ECM score does not significantly predict overall survival (log-rank p=0.83). Documented as finding consistent with published PDAC literature |
| Machine learning | Random Forest CV AUC=0.87 for ECM-high vs ECM-low classification. L1 logistic regression selected 26 sparse features |
| Target prioritisation | 30-target shortlist. All Tier 1 (multi-omics concordant). Transparent 6-layer weighted composite scoring |

### Top 10 Candidate ECM Targets

| Rank | Gene | RNA log2FC | Protein log2FC | Composite Score | Biological Role |
|---|---|---|---|---|---|
| 1 | SERPINB5 | 5.90 | 1.74 | 0.975 | Serine protease inhibitor, known PDAC marker |
| 2 | S100A11 | 1.56 | 1.33 | 0.967 | Calcium-binding protein, ECM remodelling signalling |
| 3 | IL1RN | 2.49 | 1.53 | 0.963 | IL-1 receptor antagonist, inflammatory ECM regulation |
| 4 | INHBA | 1.75 | 1.64 | 0.962 | TGF-beta superfamily member, master fibrosis driver |
| 5 | COL11A1 | 2.79 | 1.45 | 0.959 | Collagen XI alpha-1, published PDAC stromal biomarker |
| 6 | S100P | 4.93 | 2.15 | 0.955 | Calcium-binding protein, established PDAC diagnostic biomarker |
| 7 | CTHRC1 | 1.29 | 1.66 | 0.951 | Collagen triple helix repeat protein, promotes collagen deposition |
| 8 | COL12A1 | 1.55 | 1.43 | 0.951 | FACIT collagen, stromal scaffold component |
| 9 | THBS2 | 1.16 | 1.75 | 0.944 | Thrombospondin-2, ECM glycoprotein, anti-angiogenic |
| 10 | LAMA3 | 2.59 | 1.37 | 0.943 | Laminin alpha-3, basement membrane remodelling |

Notable druggable targets in the top 30: **LOXL2** (rank 13, lysyl oxidase-like 2, collagen crosslinker tested in clinical trials with simtuzumab), **PLAU** (rank 17, urokinase plasminogen activator, ECM degradation for invasion), **FN1** (rank 24, fibronectin, the most abundant ECM glycoprotein).

Downregulated targets of interest: **LGALS2** (rank 21, galectin-2), **CXCL12** (rank 29, stromal cell-derived factor 1), **SERPINI2** (rank 27) represent potentially protective factors lost during tumour ECM remodelling.

---

## Architecture

```
CPTAC PDAC Data
(140 tumour + 21 normal, matched RNA-seq + TMT proteomics + clinical metadata)
      |
      v
[Phase 1: Data Validation]              Python (Pydantic, pandas, SHA-256 checksums)
      |                                  5/5 files validated, structure and completeness checked
      v
[Phase 3: Clinical Cleaning]            Python (semicolon parsing, survival mapping)
      |                                  140 patients, 42 columns, stromal high/low groups
      |
      +------------------------------------+
      |                                    |
      v                                    v
[Phase 4: Differential Expression]   [Phase 5: ECM Signature Scoring]
R (limma, lmFit + eBayes)           R (GSVA, ssGSEA, 27 gene sets)
2,866 DE genes, 306 ECM              ECM score validated vs stromal
      |                              fraction (rho=0.39, p=1.4e-6)
      |                                    |
      v                                    v
[Phase 6: Pathway Enrichment]       [Phase 9: Survival Analysis]
R (fgsea, msigdbr)                  R (survival, survminer, Cox PH)
15/31 ECM pathways enriched          p=0.83 (documented negative result)
      |
      v
[Phase 7: Proteomics Analysis]          R (limma, MinProb imputation)
11,662 proteins, 24.6% missingness      TGF-beta axis confirmed at protein level
      |
      v
[Phase 8: Multi-Omics Integration]      Python (scipy, pandas)
8,861 matched genes, rho=0.47 global    111 concordant ECM targets (rho=0.60)
      |
      v
[Phase 10: ML and Clustering]           Python (scikit-learn)
PCA, KMeans (k=3), RF (AUC=0.87)       26 L1-selected features
      |
      v
[Phase 11: Target Prioritisation]       Python (6-layer weighted scoring)
1,070 ECM genes scored                   Top 30 shortlist, all Tier 1
```

### Language Responsibilities

| Language | Responsibility | Justification |
|---|---|---|
| R | Differential expression, ssGSEA scoring, pathway enrichment, proteomics DA, survival analysis | Gold-standard Bioconductor packages (limma, GSVA, fgsea, survminer) that are the published standard in translational oncology |
| Python | Data validation, clinical cleaning, multi-omics integration, ML, target prioritisation, testing | Stronger ecosystem for validation (Pydantic), ML (scikit-learn), data integration (pandas, scipy), and CI (pytest, ruff) |
| Nextflow | Workflow orchestration | Connects R and Python processes with dependency management and execution reporting |

---

## Dataset

**CPTAC Pancreatic Ductal Adenocarcinoma (PDAC)**

The Clinical Proteomic Tumor Analysis Consortium PDAC cohort is one of the few public datasets providing matched bulk RNA-seq transcriptomics, TMT-based quantitative proteomics, and clinical metadata for the same patients. This enables genuine multi-omics integration rather than correlating unrelated cohorts.

| Data Layer | Type | Samples | Features | Source |
|---|---|---|---|---|
| Transcriptomics | Bulk RNA-seq, RSEM upper-quartile normalised, log2 transformed | 140 tumour, 21 normal | 28,057 genes | LinkedOmics / CPTAC |
| Proteomics | TMT quantitative, median-normalised log2 abundance | 140 tumour, 75 normal | 11,662 proteins | LinkedOmics / CPTAC |
| Clinical | Demographics, staging, tissue composition, survival | 140 patients | 39 variables (42 after cleaning) | LinkedOmics / CPTAC |

**Citation:** Cao L, Huang C, Cui Zhou D, et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell. 2021;184(19):5031-5052.e26.

All data is publicly available and de-identified. No institutional review board approval is required for secondary analysis of public CPTAC data.

### Why PDAC for ECM Research

Pancreatic cancer was selected because it has the most ECM-dense stroma of any common solid tumour. The desmoplastic reaction, driven by cancer-associated fibroblasts, is a defining pathological feature and one of the most active areas in translational oncology research. The CPTAC PDAC dataset includes pathologist-estimated tissue composition fractions (neoplastic cellularity, stromal fraction, acinar fraction, inflammation fraction, and others), providing a rare opportunity to validate computational ECM scores against histological ground truth.

---

## Methodology

### Differential Expression (Phase 4)

Tumour (n=140) and normal adjacent tissue (n=21) expression profiles were compared using limma (lmFit + eBayes). The choice of limma over DESeq2 was deliberate: the CPTAC data is pre-normalised RSEM upper-quartile log2-transformed continuous expression values, not raw integer counts. DESeq2 requires raw counts and models count-based distributions. limma is the statistically appropriate method for already-normalised continuous expression data and is widely used in CPTAC publications.

Low-expression genes (log2 expression below 1 in fewer than 10% of samples) were filtered, reducing 28,057 genes to 24,133. Significance thresholds: FDR < 0.05, |log2 fold change| > 1.0.

ECM gene annotation used the Naba Matrisome collection from MSigDB (1,029 matrisome genes) plus the Hallmark Epithelial-Mesenchymal Transition gene set (200 genes).

### ECM Signature Scoring (Phase 5)

Per-sample ECM activity scores were computed using single-sample Gene Set Enrichment Analysis (ssGSEA) via the GSVA R package. 27 gene sets were scored, including Naba Matrisome subdivisions (Core Matrisome, Collagens, ECM Glycoproteins, ECM Regulators, Proteoglycans, Secreted Factors, Basement Membranes), Hallmark EMT, and Reactome ECM pathways (Extracellular Matrix Organization, Collagen Formation, Collagen Degradation, Collagen Biosynthesis and Modifying Enzymes).

The primary ECM score used the Core Matrisome gene set (268 genes present in the data). Samples were split at the median into ECM-high (n=70) and ECM-low (n=70) groups.

**Validation:** The computational ECM score was correlated against the pathologist-estimated stromal fraction using Spearman rank correlation. Result: rho=0.39, p=1.4e-6 (n=140). The molecular signature captures the same biological signal visible under the microscope. A moderate (not perfect) correlation is expected because gene expression of ECM components and histological stromal proportion measure related but distinct biological layers.

### Pathway Enrichment (Phase 6)

Gene Set Enrichment Analysis was performed using fgsea with 10,000 permutations. Genes were ranked by sign(log2FC) multiplied by -log10(p-value), which captures both direction and significance. Gene set collections included MSigDB Hallmark (50 sets), Naba Matrisome core (10 sets), and Reactome ECM-related pathways (33 sets filtered by ECM/collagen/integrin/matrix keywords). Size filter: 15-500 genes per set.

15 of 31 ECM-related pathways reached significance (FDR < 0.05). All were upregulated in tumour, indicating coordinated activation of the entire ECM programme rather than isolated gene-level changes. The top Hallmark pathways (E2F Targets, G2M Checkpoint, EMT, KRAS Signaling Up, Glycolysis, Hypoxia) form a coherent picture of PDAC biology: proliferation, invasion, metabolic reprogramming, and stromal remodelling.

### Proteomics Analysis (Phase 7)

TMT proteomics data (11,662 proteins, 140 tumour and 75 normal samples) was analysed using limma for differential abundance. Proteomics data had 24.6% missingness, typical for TMT mass spectrometry where low-abundance proteins fall below the detection threshold.

Missing values were imputed using the MinProb method: each missing value was replaced with a random draw from a normal distribution centred at the 1st percentile of observed values in that sample, with standard deviation equal to 30% of the observed standard deviation. This reflects the assumption that missing proteins are present at low abundance rather than truly absent.

Proteins detected in fewer than 50% of samples in either group were removed. Significance thresholds: FDR < 0.05, |log2FC| > 0.5 (lower than RNA because protein fold changes are typically compressed relative to RNA).

### Multi-Omics Integration (Phase 8)

RNA differential expression results and proteomics differential abundance results were merged on gene symbol. 8,861 gene-protein pairs were matched. For each pair, concordance was assessed: concordant means significant in both data layers with the same direction of change.

Global RNA-protein correlation: Spearman rho=0.47 (p approximately 0). ECM-specific correlation: rho=0.60 (p=8.54e-54). The stronger ECM correlation reflects the biology: structural ECM proteins (collagens, laminins, fibronectin) are produced in bulk with minimal post-transcriptional regulation, making RNA a reliable predictor of protein abundance for this gene class.

111 ECM genes are concordant (85 upregulated, 26 downregulated). 24 are discordant (significant in both but opposite directions). The 20:1 concordant-to-discordant ratio indicates strong multi-omics agreement.

### Survival Analysis (Phase 9)

Kaplan-Meier survival curves were generated for ECM-high vs ECM-low groups. Log-rank test: p=0.83. Median survival: ECM-low 19.7 months, ECM-high 19.8 months. Cox proportional hazards regression (univariate): HR=0.002, p=0.12. Multivariate (adjusting for age and stage): ECM score p=0.08.

**This is a negative result and is reported as such.** The ECM score does not significantly predict overall survival in this PDAC cohort. This is consistent with published literature: the relationship between stromal density and pancreatic cancer outcomes is genuinely debated, with studies reporting protective, harmful, and null associations depending on cohort, stroma definition, and analytical approach. PDAC is uniformly aggressive (median follow-up 13.1 months in this cohort) and molecular stratification of survival is challenging when the overall prognosis is poor.

A significant sex difference in ECM scores was observed (Wilcoxon p=0.026) and a borderline multivariate association (p=0.08) that may reach significance in larger cohorts. These are documented as hypothesis-generating observations.

### ML and Clustering (Phase 10)

Principal Component Analysis on 2,866 DE genes showed PC1 captures 23% of variance. K-means clustering identified 3 optimal clusters (silhouette score 0.59), with a borderline association to ECM groups (chi-squared p=0.056).

Random Forest classification (500 trees, 5-fold stratified CV) achieved AUC=0.87 for predicting ECM-high vs ECM-low from gene expression. L1-penalised logistic regression achieved AUC=0.86 and selected 26 sparse features. Top features include FGF10, SLIT3, SRPX, COL25A1 (ECM-associated) alongside non-ECM genes, indicating that ECM status is learnable from expression data and the signal is concentrated in specific genes.

### Target Prioritisation (Phase 11)

All 1,070 ECM genes present in the RNA dataset were scored using a transparent weighted composite across six evidence layers:

| Evidence Layer | Weight | Rationale |
|---|---|---|
| RNA significance (lower padj = better) | 20% | Statistical confidence of differential expression |
| RNA effect size (higher \|log2FC\| = better) | 10% | Magnitude of transcriptomic change |
| Protein significance (lower padj = better) | 20% | Independent proteomic confirmation |
| Protein effect size (higher \|log2FC\| = better) | 10% | Magnitude of protein-level change |
| Multi-omics concordance (binary) | 25% | Highest weight because concordance is the strongest single evidence type |
| ML feature importance (Random Forest) | 15% | Data-driven importance for ECM stratification |

Evidence tiers:
- **Tier 1: Multi-omics concordant** (111 genes) -- significant in both RNA and protein, same direction
- **Tier 2: RNA + Protein** (10 genes) -- significant in both but not concordant
- **Tier 3: RNA only** (949 genes) -- significant in RNA, no protein confirmation

All top 30 candidates are Tier 1.

---

## Challenges and Solutions

### Vital status labels do not match common conventions

The CPTAC clinical metadata encodes vital status as "Deceased" and "Living" rather than the more common "Dead" and "Alive". The survival mapping function handles both conventions to prevent silent data loss when computing binary event indicators.

### Tissue composition stored as semicolon-separated multi-reviewer estimates

Pathologist-estimated tissue fractions (stromal fraction, neoplastic cellularity, and others) contain multiple reviewer estimates separated by semicolons. For example, a Stromal_fraction value of "75;55;53" represents three independent pathologist estimates of 75%, 55%, and 53%. The cleaning pipeline parses these values, averages them, and converts percentages to 0-1 fractions. This is documented in the clinical cleaning module and the clinical data dictionary.

### Proteomics data has 24.6% missingness

TMT mass spectrometry has inherent detection limits. Low-abundance proteins produce missing values that are not randomly distributed: they are biased toward proteins present at low concentration. MinProb imputation addresses this by replacing missing values with small values drawn from the low end of the observed distribution, reflecting the assumption that undetected proteins are present at sub-threshold abundance rather than truly absent. The imputation method and its assumptions are documented in the evidence report.

### ECM score does not predict survival

The initial hypothesis was that ECM-high tumours would show worse outcomes. The data does not support this (log-rank p=0.83). Rather than omitting this result, it is reported transparently and contextualised against the published literature. A borderline multivariate association (p=0.08) suggests the relationship may reach significance in larger cohorts or with refined ECM scoring approaches.

### R package compilation failures

The BH (Boost Headers) R package requires C++14, but R 4.3.3 defaults to C++11 for packages requesting it. The fgsea package failed to compile until `CXX11STD = -std=gnu++14` was added to `~/.R/Makevars`. The fs package required `libuv1-dev` and gdtools required `libcairo2-dev`. These system dependencies are documented in the reproducibility guide.

### Pre-normalised expression data requires limma, not DESeq2

The CPTAC RNA-seq data from LinkedOmics is RSEM upper-quartile normalised and log2 transformed. DESeq2 requires raw integer counts and models count-based negative binomial distributions. Using DESeq2 on pre-normalised continuous data would be statistically incorrect. limma with lmFit and eBayes is the appropriate method for this data type and is the standard approach in CPTAC publications.

---

## Evidence and Reproducibility

Every analytical phase produces a structured JSON evidence report in `evidence/logs/`. Each report includes a timestamp, method description, software versions, parameter values, and key numerical results. The evidence is generated by running the analysis scripts, not by manual documentation.

| Phase | Evidence File | Reproducibility Command |
|---|---|---|
| 1. Data validation | `data_validation_report.json` | `python -m src.python.validate_data` |
| 3. Clinical cleaning | `clinical_cleaning_report.json` | `python -m src.python.clean_clinical` |
| 4. Differential expression | `de_analysis_report.json` | `Rscript src/r/01_differential_expression.R` |
| 5. ECM scoring | `ecm_scoring_report.json` | `Rscript src/r/02_ecm_signature_scoring.R` |
| 6. Pathway enrichment | `enrichment_report.json` | `Rscript src/r/03_pathway_enrichment.R` |
| 7. Proteomics | `proteomics_report.json` | `Rscript src/r/04_proteomics_analysis.R` |
| 9. Survival | `survival_report.json` | `Rscript src/r/05_survival_analysis.R` |
| 8. Multi-omics | `multiomics_report.json` | `python -m src.python.multiomics_integration` |
| 10. ML/clustering | `ml_report.json` | `python -m src.python.ml_clustering` |
| 11. Prioritisation | `prioritisation_report.json` | `python -m src.python.target_prioritisation` |

Data files are not committed to the repository. Raw data is downloaded from LinkedOmics following instructions in `metadata/data_sources.md`. File integrity is verified against SHA-256 checksums in `metadata/checksums.txt`.

---

## Testing

```
25 tests passing
```

Tests cover configuration loading, data validation (checksums, matrix structure, clinical structure), clinical cleaning (survival mapping, semicolon parsing, stromal groups, data dictionary), PCA dimensionality reduction, and target prioritisation scoring (ascending/descending, ties, NaN handling, output range).

CI pipeline runs on every push: ruff lint, pytest, and R script syntax checking.

---

## Repository Structure

```
ecm-target-discovery-pipeline/
├── README.md
├── params.yaml                          # All configurable parameters
├── environment.yml                      # Python conda environment
├── install_r_packages.R                 # R package installation
├── pyproject.toml                       # Python project config (pytest, ruff)
├── main.nf                              # Nextflow workflow
├── nextflow.config                      # Nextflow configuration
│
├── .github/workflows/ci.yml            # GitHub Actions CI
│
├── src/
│   ├── r/
│   │   ├── 01_differential_expression.R # limma DE analysis
│   │   ├── 02_ecm_signature_scoring.R   # ssGSEA via GSVA
│   │   ├── 03_pathway_enrichment.R      # fgsea GSEA
│   │   ├── 04_proteomics_analysis.R     # Proteomics DA
│   │   └── 05_survival_analysis.R       # KM + Cox PH
│   └── python/
│       ├── config.py                    # Parameter loading
│       ├── validate_data.py             # Data validation
│       ├── clean_clinical.py            # Clinical metadata cleaning
│       ├── multiomics_integration.py    # RNA-protein concordance
│       ├── ml_clustering.py             # PCA, KMeans, RF, L1 LR
│       └── target_prioritisation.py     # Multi-evidence scoring
│
├── tests/                               # 25 passing tests
│   ├── test_config.py
│   ├── test_validate_data.py
│   ├── test_clean_clinical.py
│   ├── test_ml_clustering.py
│   └── test_target_prioritisation.py
│
├── docs/
│   ├── project-charter.md
│   ├── decision-log.md
│   ├── architecture.md
│   ├── reproducibility.md
│   ├── fair-principles.md
│   └── limitations.md
│
├── metadata/
│   ├── data_sources.md                  # Download instructions and citations
│   ├── clinical_dictionary.csv          # 42-column data dictionary
│   └── checksums.txt                    # SHA-256 file integrity
│
├── evidence/logs/                       # 10 structured JSON evidence reports
├── figures/                             # Generated visualisations
├── results/                             # Analysis outputs (gitignored)
└── resources/gene_sets/                 # Curated gene set resources
```

---

## Quickstart

### Prerequisites

- R >= 4.3 with Bioconductor packages
- Python 3.12
- Nextflow >= 24.0 (optional, for full workflow orchestration)

### Setup

```bash
git clone https://github.com/gbadedata/ecm-target-discovery-pipeline.git
cd ecm-target-discovery-pipeline

# Python environment
python3 -m venv .venv
source .venv/bin/activate
pip install pyyaml pandas numpy scipy scikit-learn statsmodels lifelines \
    matplotlib seaborn plotly pydantic pytest pytest-cov ruff structlog \
    tqdm gseapy pyarrow openpyxl

# R packages (15-30 minutes first time)
Rscript install_r_packages.R

# Download data
cd data/raw
wget -O rnaseq_tumor.cct "https://linkedomics.org/data_download/CPTAC-PDAC/mRNA_RSEM_UQ_log2_Tumor.cct"
wget -O rnaseq_normal.cct "https://linkedomics.org/data_download/CPTAC-PDAC/mRNA_RSEM_UQ_log2_Normal.cct"
wget -O proteomics_tumor.cct "https://linkedomics.org/data_download/CPTAC-PDAC/proteomics_gene_level_MD_abundance_tumor.cct"
wget -O proteomics_normal.cct "https://linkedomics.org/data_download/CPTAC-PDAC/proteomics_gene_level_MD_abundance_normal.cct"
wget -O clinical.tsv "https://linkedomics.org/data_download/CPTAC-PDAC/clinical_table_140.tsv"
cd ../..

# Validate data
python -m src.python.validate_data

# Run full pipeline
nextflow run main.nf -params-file params.yaml

# Or step by step
python -m src.python.clean_clinical
Rscript src/r/01_differential_expression.R
Rscript src/r/02_ecm_signature_scoring.R
Rscript src/r/03_pathway_enrichment.R
Rscript src/r/04_proteomics_analysis.R
Rscript src/r/05_survival_analysis.R
python -m src.python.multiomics_integration
python -m src.python.ml_clustering
python -m src.python.target_prioritisation
```

---

## Key Decisions

| Decision | Choice | Reasoning |
|---|---|---|
| Dataset | CPTAC PDAC | Only public dataset with matched RNA + protein + clinical for the same patients in a cancer type where ECM is central to pathology |
| DE method | limma (not DESeq2) | Data is pre-normalised RSEM-UQ log2 continuous values. DESeq2 requires raw counts. Using it would be statistically incorrect |
| Signature scoring | ssGSEA via GSVA | Per-sample continuous scores enabling downstream stratification and correlation. Standard in translational oncology |
| Proteomics imputation | MinProb | Biologically motivated: missing proteins are likely low-abundance, not absent. Standard approach for TMT data |
| Proteomics DA | limma | Robust linear modelling appropriate for log2-transformed TMT abundance data |
| Survival | Kaplan-Meier + Cox PH | Standard clinical association framework. Negative result reported transparently |
| Target scoring weights | Concordance 25%, significance 20% each, ML 15%, effect size 10% each | Multi-omics concordance is the strongest single evidence type and receives the highest weight |
| Workflow | Nextflow DSL2 | Standard bioinformatics workflow manager connecting R and Python processes |

Full decision log: [docs/decision-log.md](docs/decision-log.md)

---

## Limitations

- **Single cohort.** All results are from one CPTAC PDAC cohort (n=140). Findings are hypothesis-generating and require external validation.
- **Bulk transcriptomics.** Bulk RNA-seq averages signal across all cell types. ECM scores may reflect stromal cell content rather than tumour-intrinsic biology.
- **Proteomics coverage.** TMT proteomics detected 11,662 proteins out of approximately 20,000 protein-coding genes. Absence from the dataset does not mean absence from tissue.
- **Imputation assumptions.** MinProb imputation assumes missing proteins are low-abundance. This is generally valid for TMT data but introduces noise for proteins that are genuinely absent in specific samples.
- **No wet-lab validation.** All findings are computational. Candidate targets are prioritised, not validated.
- **Association, not causation.** Differential expression and clinical associations do not establish causal mechanisms.
- **Survival power.** With 135 patients and 76 events, statistical power for survival subgroup analyses is limited.

Full limitations documentation: [docs/limitations.md](docs/limitations.md)

---

## Numbers

| Metric | Value |
|---|---|
| Tumour samples | 140 |
| Normal samples | 21 (RNA), 75 (protein) |
| Genes analysed | 28,057 (24,133 after filtering) |
| Proteins analysed | 11,662 (8,876 after filtering) |
| Proteomics missingness | 24.6% |
| DE significant genes | 2,866 |
| ECM genes significant (RNA) | 306 |
| ECM pathways enriched | 15/31 |
| Multi-omics concordant ECM targets | 111 |
| RNA-protein Spearman rho (ECM) | 0.60 |
| ECM score vs stromal fraction rho | 0.39 (p=1.4e-6) |
| Survival log-rank p | 0.83 (not significant) |
| RF classification AUC | 0.87 |
| Final target shortlist | 30 (all Tier 1) |
| Tests passing | 25 |
| Evidence reports | 10 |
| CI pipeline | Green (ruff + pytest + R syntax) |

---

## Tools and Technologies

**R:** limma, GSVA, fgsea, msigdbr, survival, survminer, tidyverse, pheatmap, ggrepel, ggplot2

**Python:** pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, Pydantic, pytest, ruff, structlog

**Workflow and Infrastructure:** Nextflow DSL2, GitHub Actions CI, conda/venv

---

## Documentation

- [Project Charter](docs/project-charter.md)
- [Architecture](docs/architecture.md)
- [Decision Log](docs/decision-log.md)
- [Data Sources](metadata/data_sources.md)
- [Clinical Dictionary](metadata/clinical_dictionary.csv)
- [Reproducibility Guide](docs/reproducibility.md)
- [FAIR Principles](docs/fair-principles.md)
- [Limitations](docs/limitations.md)

---

## Author

**Oluwagbade Odimayo**
Bioinformatics Data Engineer | MSc Applied Data Science

- [GitHub](https://github.com/gbadedata)
- [LinkedIn](https://www.linkedin.com/in/oluwagbade-odimayo-)
