# Data Sources

## Primary Dataset: CPTAC Pancreatic Ductal Adenocarcinoma (PDAC)

### Overview

The Clinical Proteomic Tumor Analysis Consortium (CPTAC) PDAC cohort provides matched multi-omics data for pancreatic ductal adenocarcinoma patients. This is one of the few public datasets offering both bulk RNA-seq transcriptomics and TMT-based quantitative proteomics for the same patient samples, alongside clinical metadata.

### Data Layers

| Layer | Type | Source | Format |
|---|---|---|---|
| Transcriptomics | Bulk RNA-seq (gene-level counts or FPKM/TPM) | CPTAC Data Portal / GDC | TSV/CSV |
| Proteomics | TMT quantitative proteomics (log2 ratios) | CPTAC Data Portal / LinkedOmics | TSV/CSV |
| Clinical | Patient demographics, tumour stage, survival, histology | CPTAC Data Portal / GDC | TSV/CSV |

### Access

- **CPTAC Data Portal:** https://proteomics.cancer.gov/data-portal
- **GDC (Genomic Data Commons):** https://portal.gdc.cancer.gov/
- **LinkedOmics:** http://www.linkedomics.org/
- **cBioPortal (CPTAC-3 PDAC):** https://www.cbioportal.org/

### Citation

Cao L, Huang C, Cui Zhou D, et al. Proteogenomic characterization of pancreatic ductal adenocarcinoma. Cell. 2021;184(19):5031-5052.e26.

### Ethical and Legal Status

All data is publicly available and de-identified. No institutional review board approval is required for secondary analysis of public CPTAC data. The project uses no private patient data.

## Gene Set Resources

| Resource | Description | URL |
|---|---|---|
| MSigDB Hallmark | Curated gene sets including EMT, inflammatory response | https://www.gsea-msigdb.org/ |
| Naba Matrisome | Comprehensive ECM gene annotation (core matrisome + matrisome-associated) | https://matrisome.org/ |
| Reactome | Pathway gene sets for ECM organisation, collagen formation, degradation | https://reactome.org/ |
| GO Biological Process | Gene Ontology terms for ECM, fibrosis, collagen | https://geneontology.org/ |

## Data Not Committed to Git

Raw data files are not stored in this repository. The `data/` directory is in `.gitignore`. To reproduce the analysis:

1. Download data from the sources listed above
2. Place files in `data/raw/` following the structure in `metadata/sample_manifest.csv`
3. Run the data validation step to confirm file integrity against `metadata/checksums.txt`
