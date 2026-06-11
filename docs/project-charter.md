# Project Charter

**ECM Target Discovery Pipeline: A Reproducible Multi-Omics Bioinformatics Workflow for Fibrosis and Oncology Research**

## Core Scientific Question

Which ECM-associated transcriptomic and proteomic signatures are associated with disease biology, tumour microenvironment remodelling, and clinical outcomes in pancreatic ductal adenocarcinoma?

## Dataset

CPTAC PDAC cohort: matched bulk RNA-seq, TMT proteomics, and clinical metadata.

## Scope

1. Data ingestion and validation
2. Bulk transcriptomics analysis (DESeq2)
3. ECM and fibrosis gene signature scoring (ssGSEA via GSVA)
4. Pathway enrichment analysis (fgsea, clusterProfiler)
5. Quantitative proteomics analysis (limma)
6. Multi-omics integration (RNA-protein correlation)
7. Clinical association and survival analysis
8. Dimensionality reduction, clustering, and ML stratification
9. Candidate ECM target prioritisation
10. Nextflow workflow orchestration
11. Testing, CI, and evidence capture
12. Scientific reporting and documentation

## Out of Scope

Wet-lab validation, clinical diagnostic claims, Kubernetes, streaming pipelines, production web applications, private patient data.

## Quality Bar

Real public data, reproducible analysis, correct statistics, conservative interpretation, documented limitations, evidence capture, full test suite.

## Guiding Principle

Every tool must earn its place. The centre of gravity is biological interpretation and translational relevance, not infrastructure.

The full charter is maintained separately and informed the design of this repository.
