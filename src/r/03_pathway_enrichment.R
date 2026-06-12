#!/usr/bin/env Rscript
# =============================================================================
# Phase 6: Pathway Enrichment Analysis (fgsea)
#
# Tests whether ECM, fibrosis, and TME-related pathways are enriched
# among differentially expressed genes using Gene Set Enrichment Analysis.
#
# Inputs:
#   results/de_results_all.csv
#
# Outputs:
#   results/fgsea_results_all.csv
#   results/fgsea_results_significant.csv
#   results/fgsea_results_ecm_focused.csv
#   figures/enrichment_barplot.png
#   figures/enrichment_dotplot_ecm.png
#   figures/gsea_waterfall.png
#   evidence/logs/enrichment_report.json
# =============================================================================

suppressPackageStartupMessages({
    library(fgsea)
    library(msigdbr)
    library(tidyverse)
    library(jsonlite)
})

set.seed(42)

cat("=== Phase 6: Pathway Enrichment Analysis ===\n\n")

# ---- 1. Load DE results and build ranked list --------------------------------

cat("Loading DE results...\n")
de <- read.csv("results/de_results_all.csv")
cat(sprintf("  Genes: %d\n", nrow(de)))

# Create ranked gene list: sign(log2FC) * -log10(p-value)
# This ranks genes by both magnitude and significance
de <- de %>%
    filter(!is.na(p_value) & !is.na(log2fc)) %>%
    mutate(rank_metric = sign(log2fc) * -log10(p_value))

ranked_genes <- setNames(de$rank_metric, de$gene)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

cat(sprintf("  Ranked genes: %d\n", length(ranked_genes)))
cat(sprintf("  Top ranked (up): %s (%.1f)\n",
            names(ranked_genes)[1], ranked_genes[1]))
cat(sprintf("  Bottom ranked (down): %s (%.1f)\n",
            names(ranked_genes)[length(ranked_genes)],
            ranked_genes[length(ranked_genes)]))

# ---- 2. Build gene set collections ------------------------------------------

cat("\nBuilding gene set collections...\n")

# Hallmark
hallmark <- msigdbr(species = "Homo sapiens") %>%
    filter(grepl("^HALLMARK_", gs_name))
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)

# Naba Matrisome (core sets only, not cancer-specific subsets)
naba <- msigdbr(species = "Homo sapiens") %>%
    filter(gs_name %in% c(
        "NABA_CORE_MATRISOME", "NABA_MATRISOME",
        "NABA_MATRISOME_ASSOCIATED",
        "NABA_COLLAGENS", "NABA_ECM_GLYCOPROTEINS",
        "NABA_ECM_REGULATORS", "NABA_ECM_AFFILIATED",
        "NABA_PROTEOGLYCANS", "NABA_SECRETED_FACTORS",
        "NABA_BASEMENT_MEMBRANES"
    ))
naba_sets <- split(naba$gene_symbol, naba$gs_name)

# Reactome ECM and fibrosis-related
reactome <- msigdbr(species = "Homo sapiens") %>%
    filter(grepl("^REACTOME_", gs_name)) %>%
    filter(grepl(
        "EXTRACELLULAR|COLLAGEN|INTEGRIN|LAMININ|MATRIX|ELASTIC|PROTEOGLYCAN|FIBR",
        gs_name
    ))
reactome_sets <- split(reactome$gene_symbol, reactome$gs_name)

# Combine
all_sets <- c(hallmark_sets, naba_sets, reactome_sets)

# Filter to sets with 15-500 genes
all_sets <- all_sets[sapply(all_sets, length) >= 15 &
                     sapply(all_sets, length) <= 500]

cat(sprintf("  Hallmark sets: %d\n", length(hallmark_sets)))
cat(sprintf("  Naba Matrisome sets: %d\n", length(naba_sets)))
cat(sprintf("  Reactome ECM sets: %d\n", length(reactome_sets)))
cat(sprintf("  Total after size filter: %d\n", length(all_sets)))

# ---- 3. Run fgsea ------------------------------------------------------------

cat("\nRunning fgsea...\n")

fgsea_results <- fgsea(
    pathways = all_sets,
    stats = ranked_genes,
    minSize = 15,
    maxSize = 500,
    nPermSimple = 10000
)

fgsea_results <- fgsea_results %>%
    arrange(padj) %>%
    mutate(
        direction = ifelse(NES > 0, "Up in Tumour", "Down in Tumour"),
        significant = padj < 0.05,
        # Clean pathway names for display
        pathway_short = gsub("^HALLMARK_|^NABA_|^REACTOME_", "", pathway),
        pathway_short = gsub("_", " ", pathway_short),
        pathway_short = str_to_title(pathway_short)
    )

n_sig <- sum(fgsea_results$significant)
n_up <- sum(fgsea_results$significant & fgsea_results$NES > 0)
n_down <- sum(fgsea_results$significant & fgsea_results$NES < 0)

cat(sprintf("  Tested: %d gene sets\n", nrow(fgsea_results)))
cat(sprintf("  Significant (FDR<0.05): %d (%d up, %d down)\n",
            n_sig, n_up, n_down))

# ---- 4. ECM-focused results -------------------------------------------------

cat("\nECM-focused pathway results...\n")

ecm_keywords <- c("NABA_", "MATRISOME", "COLLAGEN", "EXTRACELLULAR",
                   "MATRIX", "LAMININ", "INTEGRIN", "EMT", "ELASTIC",
                   "PROTEOGLYCAN", "BASEMENT")

ecm_results <- fgsea_results %>%
    filter(grepl(paste(ecm_keywords, collapse = "|"), pathway,
                 ignore.case = TRUE))

cat(sprintf("  ECM-related pathways: %d\n", nrow(ecm_results)))
cat(sprintf("  ECM significant: %d\n",
            sum(ecm_results$significant)))

cat("\n--- ECM Pathway Results ---\n")
ecm_display <- ecm_results %>%
    select(pathway_short, NES, padj, size, direction) %>%
    mutate(
        NES = round(NES, 3),
        padj = format(padj, digits = 3, scientific = TRUE)
    )
print(as.data.frame(ecm_display), row.names = FALSE)

# ---- 5. Top Hallmark results ------------------------------------------------

cat("\n--- Top 15 Hallmark Pathways ---\n")

hallmark_results <- fgsea_results %>%
    filter(grepl("^HALLMARK_", pathway)) %>%
    arrange(padj) %>%
    head(15)

hallmark_display <- hallmark_results %>%
    select(pathway_short, NES, padj, direction) %>%
    mutate(
        NES = round(NES, 3),
        padj = format(padj, digits = 3, scientific = TRUE)
    )
print(as.data.frame(hallmark_display), row.names = FALSE)

# ---- 6. Save results --------------------------------------------------------

cat("\nSaving results...\n")

# Convert leadingEdge list column to string for CSV
fgsea_out <- fgsea_results %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    select(pathway, pval, padj, log2err, ES, NES, size,
           leadingEdge, direction, significant, pathway_short)

write.csv(fgsea_out, "results/fgsea_results_all.csv", row.names = FALSE)

sig_out <- fgsea_out %>% filter(significant)
write.csv(sig_out, "results/fgsea_results_significant.csv",
          row.names = FALSE)

ecm_out <- fgsea_out %>%
    filter(grepl(paste(ecm_keywords, collapse = "|"), pathway,
                 ignore.case = TRUE))
write.csv(ecm_out, "results/fgsea_results_ecm_focused.csv",
          row.names = FALSE)

cat(sprintf("  Saved: results/fgsea_results_all.csv (%d rows)\n",
            nrow(fgsea_out)))
cat(sprintf("  Saved: results/fgsea_results_significant.csv (%d rows)\n",
            nrow(sig_out)))
cat(sprintf("  Saved: results/fgsea_results_ecm_focused.csv (%d rows)\n",
            nrow(ecm_out)))

# ---- 7. Enrichment bar plot (top pathways) -----------------------------------

cat("\nGenerating figures...\n")

top_pathways <- fgsea_results %>%
    filter(significant) %>%
    arrange(desc(abs(NES))) %>%
    head(25)

p_bar <- ggplot(top_pathways,
                aes(x = reorder(pathway_short, NES), y = NES,
                    fill = direction)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(
        values = c("Up in Tumour" = "#E74C3C",
                    "Down in Tumour" = "#3498DB"),
        name = "Direction"
    ) +
    labs(
        title = "Top 25 Enriched Pathways (PDAC Tumour vs Normal)",
        subtitle = "fgsea, FDR < 0.05, ranked by |NES|",
        x = NULL,
        y = "Normalised Enrichment Score (NES)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

ggsave("figures/enrichment_barplot.png", p_bar,
       width = 10, height = 8, dpi = 300)
cat("  Saved: figures/enrichment_barplot.png\n")

# ---- 8. ECM dot plot ---------------------------------------------------------

ecm_sig <- ecm_results %>% filter(significant | padj < 0.1)

if (nrow(ecm_sig) > 0) {
    p_ecm_dot <- ggplot(
        ecm_sig,
        aes(x = NES, y = reorder(pathway_short, NES),
            size = size, colour = -log10(padj))
    ) +
        geom_point(alpha = 0.8) +
        scale_colour_gradient(low = "grey60", high = "#8E44AD",
                              name = "-log10(FDR)") +
        scale_size_continuous(range = c(3, 10), name = "Gene Set Size") +
        geom_vline(xintercept = 0, linetype = "dashed",
                   colour = "grey50") +
        labs(
            title = "ECM-Related Pathway Enrichment in PDAC",
            subtitle = "fgsea, Naba Matrisome + Reactome ECM pathways",
            x = "Normalised Enrichment Score (NES)",
            y = NULL
        ) +
        theme_minimal(base_size = 11)

    ggsave("figures/enrichment_dotplot_ecm.png", p_ecm_dot,
           width = 10, height = 7, dpi = 300)
    cat("  Saved: figures/enrichment_dotplot_ecm.png\n")
}

# ---- 9. Waterfall plot of all significant pathways ---------------------------

sig_pathways <- fgsea_results %>%
    filter(significant) %>%
    arrange(NES)

p_waterfall <- ggplot(sig_pathways,
                      aes(x = reorder(pathway_short, NES), y = NES,
                          fill = direction)) +
    geom_col(alpha = 0.75, width = 0.7) +
    coord_flip() +
    scale_fill_manual(
        values = c("Up in Tumour" = "#E74C3C",
                    "Down in Tumour" = "#3498DB"),
        name = "Direction"
    ) +
    labs(
        title = sprintf(
            "All Significant Pathways (n=%d, FDR<0.05)",
            nrow(sig_pathways)
        ),
        x = NULL,
        y = "Normalised Enrichment Score (NES)"
    ) +
    theme_minimal(base_size = 8) +
    theme(legend.position = "bottom")

ggsave("figures/gsea_waterfall.png", p_waterfall,
       width = 10, height = max(6, nrow(sig_pathways) * 0.2),
       dpi = 300)
cat("  Saved: figures/gsea_waterfall.png\n")

# ---- 10. Evidence report -----------------------------------------------------

report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    analysis = "pathway_enrichment",
    method = "fgsea",
    ranking_metric = "sign(log2FC) * -log10(pvalue)",
    gene_sets = list(
        hallmark = length(hallmark_sets),
        naba = length(naba_sets),
        reactome_ecm = length(reactome_sets),
        total_tested = nrow(fgsea_results)
    ),
    results = list(
        significant_total = n_sig,
        up_in_tumour = n_up,
        down_in_tumour = n_down,
        ecm_pathways_tested = nrow(ecm_results),
        ecm_pathways_significant = sum(ecm_results$significant)
    ),
    top_ecm_pathways = ecm_results %>%
        filter(significant) %>%
        arrange(padj) %>%
        head(10) %>%
        select(pathway, NES, padj) %>%
        mutate(NES = round(NES, 3)) %>%
        as.list(),
    top_hallmark_pathways = hallmark_results %>%
        head(10) %>%
        select(pathway, NES, padj) %>%
        mutate(NES = round(NES, 3)) %>%
        as.list(),
    r_version = R.version.string,
    fgsea_version = as.character(packageVersion("fgsea"))
)

write_json(report, "evidence/logs/enrichment_report.json",
           pretty = TRUE, auto_unbox = TRUE)
cat("  Evidence: evidence/logs/enrichment_report.json\n")

cat("\n=== Phase 6 Complete ===\n")
