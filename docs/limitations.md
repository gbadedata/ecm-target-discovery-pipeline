# Limitations

This document records known limitations, assumptions, and caveats. It is updated as the project progresses.

## Dataset Limitations

- **Single cohort.** The analysis uses one CPTAC PDAC cohort. Findings are hypothesis-generating, not externally validated.
- **Sample size.** CPTAC PDAC includes approximately 140 tumour samples. Statistical power for subgroup analyses is limited.
- **Bulk transcriptomics.** Bulk RNA-seq averages signal across all cell types in the tissue. ECM signatures may reflect stromal content rather than tumour-intrinsic biology. This is a known limitation of bulk deconvolution and is acknowledged, not corrected.

## Analytical Limitations

- **No wet-lab validation.** All findings are computational. Candidate targets are prioritised, not validated.
- **Proteomics coverage.** TMT proteomics detects fewer proteins than the genome encodes genes. Absence of a protein from the dataset does not mean absence from the tissue.
- **Missing data.** Proteomics data has inherent missingness. Imputation introduces assumptions. The imputation method and its limitations are documented.
- **Multiple testing.** Multiple comparisons are corrected using Benjamini-Hochberg FDR. Residual false discovery risk remains at the stated threshold.

## Interpretation Limitations

- **Association, not causation.** Differential expression and clinical association do not establish causal relationships.
- **ECM score stratification.** Splitting samples into ECM-high and ECM-low groups by median is a common but arbitrary threshold. Sensitivity to threshold choice is assessed.
- **Survival analysis.** Kaplan-Meier and Cox regression identify associations with outcomes. Confounders (stage, grade, treatment) are adjusted for where metadata permits, but residual confounding is possible.

## Scope Limitations

- This project does not claim to discover a novel therapeutic target.
- This project does not provide clinical diagnostic recommendations.
- This project does not process raw sequencing reads (FASTQ). It uses pre-processed expression matrices.
