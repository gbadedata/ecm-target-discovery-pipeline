# Decision Log

This document records key decisions made during the project, the reasoning behind them, and the alternatives considered.

## DL-001: Dataset Selection -- CPTAC PDAC

**Date:** 2026-06-11
**Decision:** Use the CPTAC Pancreatic Ductal Adenocarcinoma (PDAC) cohort as the flagship dataset.
**Reasoning:**
- PDAC has the most ECM-dense stroma of any common solid tumour. The desmoplastic reaction, driven by cancer-associated fibroblasts, is a defining pathological feature.
- CPTAC provides matched bulk RNA-seq, TMT proteomics, and clinical metadata for the same patients. This enables genuine multi-omics integration rather than correlating unrelated cohorts.
- The ECM/fibrosis narrative is strongest in PDAC: stromal remodelling influences drug delivery, immune exclusion, and patient outcomes.
- CPTAC data is publicly available, well-documented, and widely cited in translational research.
**Alternatives considered:**
- CPTAC-CCRCC: strong TME component but ECM remodelling is less central to the disease biology than in PDAC.
- TCGA-PAAD: RNA-seq and clinical metadata available but no matched proteomics.
- IPF lung fibrosis datasets: strong fibrosis relevance but limited to transcriptomics only.

## DL-002: R for Transcriptomics, Python for Integration and ML

**Date:** 2026-06-11
**Decision:** Use R (DESeq2, clusterProfiler, survival, GSVA) for transcriptomics, pathway analysis, and survival modelling. Use Python (pandas, scikit-learn, lifelines) for data validation, multi-omics integration, machine learning, and automation.
**Reasoning:**
- DESeq2 and clusterProfiler are the gold-standard tools in their domains. Reimplementing their statistical methods in Python would be scientifically weaker.
- Python is stronger for data validation (Pydantic), ML pipelines (scikit-learn), and orchestration.
- Bioinformatics roles expect fluency in both languages. Using each where it is strongest demonstrates practical judgment.
**Alternatives considered:**
- Python-only with pyDESeq2: possible but less established, weaker pathway analysis ecosystem.
- R-only: possible for analysis but weaker for data engineering, validation, and ML.

## DL-003: ssGSEA for ECM Signature Scoring

**Date:** 2026-06-11
**Decision:** Use single-sample Gene Set Enrichment Analysis (ssGSEA) via GSVA for computing per-sample ECM signature scores.
**Reasoning:**
- ssGSEA produces a continuous score per sample per gene set, enabling downstream stratification, survival analysis, and correlation with clinical variables.
- It does not require predefined group comparisons, making it suitable for unsupervised exploration.
- GSVA is the standard R implementation and is widely used in translational oncology publications.
**Alternatives considered:**
- Mean z-score of gene set members: simpler but does not account for gene-gene correlation or rank-based enrichment.
- AUCell: designed for single-cell data, not optimal for bulk RNA-seq.

## DL-004: Quarantine-Not-Delete for Data Validation

**Date:** 2026-06-11
**Decision:** Samples or records failing validation are flagged and documented, not silently removed.
**Reasoning:**
- Consistent with the approach used in previous portfolio projects (BioSeq, Biomarker Concordance Pipeline).
- Transparent data handling is a requirement for reproducibility and aligns with FAIR principles.
- Removed samples must be accounted for in the methods section of any scientific report.

## DL-005: Nextflow for Workflow Orchestration

**Date:** 2026-06-11
**Decision:** Use Nextflow to orchestrate the end-to-end analysis workflow.
**Reasoning:**
- Nextflow is the standard workflow manager in bioinformatics (alongside Snakemake).
- It handles R and Python processes in the same workflow, manages dependencies per process via containers or conda, and provides execution reports.
- Demonstrates continuity with the Biomarker Concordance Pipeline (also Nextflow DSL2).
**Alternatives considered:**
- Snakemake: equally valid but Nextflow is more common in clinical bioinformatics environments.
- Make/shell scripts: not reproducible across environments.
