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
    ╔══════════════════════════════════════════════╗
    ║   ECM Target Discovery Pipeline  v${workflow.manifest.version}     ║
    ║   Dataset: CPTAC PDAC                        ║
    ╚══════════════════════════════════════════════╝
    Output directory : ${params.outdir}
    Figures directory: ${params.figdir}
    """.stripIndent()

// Workflow modules will be added as each analysis phase is built.
// Each phase produces a self-contained Nextflow process.

workflow {
    log.info "Pipeline scaffold ready. Analysis modules will be added per phase."
}
