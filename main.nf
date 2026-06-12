#!/usr/bin/env nextflow

/*
 * ECM Target Discovery Pipeline
 *
 * A reproducible multi-omics bioinformatics workflow for investigating
 * ECM-associated molecular signatures in pancreatic ductal adenocarcinoma.
 *
 * Author: Oluwagbade Odimayo
 * GitHub: github.com/gbadedata/ecm-target-discovery-pipeline
 */

nextflow.enable.dsl = 2

log.info """
    =============================================
      ECM Target Discovery Pipeline  v${workflow.manifest.version}
      Dataset: CPTAC PDAC
    =============================================
    Output directory : ${params.outdir}
    Figures directory: ${params.figdir}
    """.stripIndent()

// -----------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------

process VALIDATE_DATA {
    label 'python_process'

    output:
    path "evidence/logs/data_validation_report.json", emit: report

    script:
    """
    cd ${projectDir}
    python -m src.python.validate_data
    """
}

process CLEAN_CLINICAL {
    label 'python_process'

    output:
    path "data/processed/clinical_clean.csv", emit: clinical

    script:
    """
    cd ${projectDir}
    python -m src.python.clean_clinical
    """
}

process DIFFERENTIAL_EXPRESSION {
    label 'r_process'

    input:
    path clinical

    output:
    path "results/de_results_all.csv", emit: de_all
    path "results/de_results_significant.csv", emit: de_sig
    path "results/de_results_ecm_genes.csv", emit: de_ecm

    script:
    """
    cd ${projectDir}
    Rscript src/r/01_differential_expression.R
    """
}

process ECM_SCORING {
    label 'r_process'

    input:
    path clinical

    output:
    path "results/ecm_scores.csv", emit: scores

    script:
    """
    cd ${projectDir}
    Rscript src/r/02_ecm_signature_scoring.R
    """
}

process PATHWAY_ENRICHMENT {
    label 'r_process'

    input:
    path de_all

    output:
    path "results/fgsea_results_all.csv", emit: fgsea

    script:
    """
    cd ${projectDir}
    Rscript src/r/03_pathway_enrichment.R
    """
}

process PROTEOMICS_ANALYSIS {
    label 'r_process'

    output:
    path "results/proteomics_da_all.csv", emit: prot_da

    script:
    """
    cd ${projectDir}
    Rscript src/r/04_proteomics_analysis.R
    """
}

process SURVIVAL_ANALYSIS {
    label 'r_process'

    input:
    path scores

    output:
    path "results/survival_results.csv", emit: survival

    script:
    """
    cd ${projectDir}
    Rscript src/r/05_survival_analysis.R
    """
}

process MULTIOMICS_INTEGRATION {
    label 'python_process'

    input:
    path de_all
    path prot_da

    output:
    path "results/multiomics_integrated.csv", emit: integrated
    path "results/multiomics_concordant_ecm.csv", emit: concordant

    script:
    """
    cd ${projectDir}
    python -m src.python.multiomics_integration
    """
}

process ML_CLUSTERING {
    label 'python_process'

    input:
    path de_sig
    path scores

    output:
    path "results/ml_feature_importance.csv", emit: features

    script:
    """
    cd ${projectDir}
    python -m src.python.ml_clustering
    """
}

process TARGET_PRIORITISATION {
    label 'python_process'

    input:
    path de_all
    path prot_da
    path integrated
    path features

    output:
    path "results/target_prioritisation_table.csv", emit: targets
    path "results/target_shortlist_top30.csv", emit: shortlist

    script:
    """
    cd ${projectDir}
    python -m src.python.target_prioritisation
    """
}

// -----------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------

workflow {
    // Phase 1: Validation
    VALIDATE_DATA()

    // Phase 3: Clinical cleaning
    CLEAN_CLINICAL()

    // Phase 4-5: Transcriptomics (parallel)
    DIFFERENTIAL_EXPRESSION(CLEAN_CLINICAL.out.clinical)
    ECM_SCORING(CLEAN_CLINICAL.out.clinical)

    // Phase 6: Pathway enrichment (depends on DE)
    PATHWAY_ENRICHMENT(DIFFERENTIAL_EXPRESSION.out.de_all)

    // Phase 7: Proteomics (independent)
    PROTEOMICS_ANALYSIS()

    // Phase 9: Survival (depends on ECM scores)
    SURVIVAL_ANALYSIS(ECM_SCORING.out.scores)

    // Phase 8: Multi-omics integration (depends on DE + proteomics)
    MULTIOMICS_INTEGRATION(
        DIFFERENTIAL_EXPRESSION.out.de_all,
        PROTEOMICS_ANALYSIS.out.prot_da
    )

    // Phase 10: ML (depends on DE significant + ECM scores)
    ML_CLUSTERING(
        DIFFERENTIAL_EXPRESSION.out.de_sig,
        ECM_SCORING.out.scores
    )

    // Phase 11: Target prioritisation (depends on all upstream)
    TARGET_PRIORITISATION(
        DIFFERENTIAL_EXPRESSION.out.de_all,
        PROTEOMICS_ANALYSIS.out.prot_da,
        MULTIOMICS_INTEGRATION.out.integrated,
        ML_CLUSTERING.out.features
    )
}
