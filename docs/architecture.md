# Architecture

## System Overview

The ECM Target Discovery Pipeline processes public multi-omics data through a structured analytical workflow to identify ECM-associated molecular signatures relevant to fibrosis and oncology research.

## Data Flow

```
CPTAC PDAC Data Portal
        |
        v
[1. Data Ingestion & Validation]  -- Python (Pydantic, pandas)
        |
        v
[2. Transcriptomics Analysis]     -- R (DESeq2, clusterProfiler, GSVA)
   - Normalisation
   - Differential expression
   - ECM signature scoring
   - Pathway enrichment
        |
        v
[3. Proteomics Analysis]          -- R (limma) + Python (pandas)
   - Missingness assessment
   - Filtering and normalisation
   - Differential abundance
   - ECM protein annotation
        |
        v
[4. Multi-Omics Integration]      -- Python (pandas, scipy)
   - RNA-protein correlation
   - Concordance assessment
   - Integrated feature matrix
        |
        v
[5. Clinical Association]         -- R (survival, survminer) + Python (lifelines)
   - ECM score vs clinical features
   - Survival analysis (KM, Cox)
   - Subgroup characterisation
        |
        v
[6. ML & Clustering]              -- Python (scikit-learn)
   - Dimensionality reduction (PCA, UMAP)
   - Consensus clustering
   - Feature importance ranking
        |
        v
[7. Target Prioritisation]        -- Python
   - Multi-evidence scoring
   - Candidate target table
   - Biological interpretation
        |
        v
[8. Reporting]                    -- Quarto/RMarkdown
   - Technical report
   - Scientist-facing summary
   - Figures and tables
```

## Language Responsibilities

| Language | Responsibility | Justification |
|---|---|---|
| R | Differential expression, pathway enrichment, ssGSEA scoring, survival analysis, biological visualisation | Gold-standard Bioconductor packages (DESeq2, clusterProfiler, GSVA, survminer) |
| Python | Data validation, proteomics QC, multi-omics integration, ML, workflow automation, testing | Stronger ecosystem for validation (Pydantic), ML (scikit-learn), and orchestration |
| Nextflow | Workflow orchestration | Connects R and Python processes, manages dependencies, produces execution reports |

## Key Design Principles

1. **Separation of concerns.** Each analytical step is a self-contained script or module with defined inputs and outputs.
2. **Reproducibility.** Random seeds, package versions, parameter files, and environment definitions are tracked.
3. **Evidence-first.** Every claim in the final report is traceable to a specific analysis output.
4. **Conservative interpretation.** The project identifies associations and candidates, not validated targets. Limitations are documented.
5. **FAIR alignment.** Data sources are documented, metadata is structured, analysis is reproducible, and outputs are reusable.
