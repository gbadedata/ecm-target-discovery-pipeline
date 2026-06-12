#!/usr/bin/env Rscript
# =============================================================================
# Phase 4: Differential Expression Analysis (limma)
#
# Compares CPTAC PDAC tumour (n=140) vs normal adjacent tissue (n=21)
# using limma on pre-normalised RSEM-UQ log2 expression data.
#
# Inputs:
#   data/raw/rnaseq_tumor.cct
#   data/raw/rnaseq_normal.cct
#
# Outputs:
#   results/de_results_all.csv           - Full DE table (all genes)
#   results/de_results_significant.csv   - FDR < 0.05, |log2FC| > 1
#   results/de_results_ecm_genes.csv     - ECM-associated DE genes
#   figures/volcano_plot.png             - Volcano plot
#   figures/volcano_plot_ecm.png         - Volcano plot with ECM genes highlighted
#   evidence/logs/de_analysis_report.json - Structured evidence
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(tidyverse)
    library(pheatmap)
    library(ggrepel)
    library(msigdbr)
    library(jsonlite)
})

set.seed(42)

cat("=== Phase 4: Differential Expression Analysis ===\n\n")

# ---- 1. Load data -----------------------------------------------------------

cat("Loading expression data...\n")

tumor <- read.delim("data/raw/rnaseq_tumor.cct", row.names = 1, check.names = FALSE)
normal <- read.delim("data/raw/rnaseq_normal.cct", row.names = 1, check.names = FALSE)

cat(sprintf("  Tumour: %d genes x %d samples\n", nrow(tumor), ncol(tumor)))
cat(sprintf("  Normal: %d genes x %d samples\n", nrow(normal), ncol(normal)))

# Verify same gene set
stopifnot(all(rownames(tumor) == rownames(normal)))

# Combine into single matrix
expr <- cbind(tumor, normal)
cat(sprintf("  Combined: %d genes x %d samples\n", nrow(expr), ncol(expr)))

# ---- 2. Create design matrix ------------------------------------------------

cat("\nCreating design matrix...\n")

group <- factor(
    c(rep("Tumour", ncol(tumor)), rep("Normal", ncol(normal))),
    levels = c("Normal", "Tumour")
)

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "TumourVsNormal")

cat(sprintf("  Groups: Tumour=%d, Normal=%d\n",
            sum(group == "Tumour"), sum(group == "Normal")))

# ---- 3. Filter low-expression genes ----------------------------------------

cat("\nFiltering low-expression genes...\n")

# Keep genes expressed (log2 > 1) in at least 10% of samples
min_samples <- ceiling(0.1 * ncol(expr))
keep <- rowSums(expr > 1) >= min_samples
expr_filtered <- expr[keep, ]

cat(sprintf("  Before filtering: %d genes\n", nrow(expr)))
cat(sprintf("  After filtering:  %d genes\n", nrow(expr_filtered)))
cat(sprintf("  Removed:          %d genes (%.1f%%)\n",
            nrow(expr) - nrow(expr_filtered),
            (nrow(expr) - nrow(expr_filtered)) / nrow(expr) * 100))

# ---- 4. limma analysis ------------------------------------------------------

cat("\nRunning limma...\n")

# Data is already log2-transformed RSEM-UQ, so use lmFit directly (no voom)
fit <- lmFit(as.matrix(expr_filtered), design)
fit <- eBayes(fit)

# Extract results for Tumour vs Normal contrast
results <- topTable(fit, coef = "TumourVsNormal", number = Inf, sort.by = "none")

# Clean up column names
results$gene <- rownames(results)
results <- results %>%
    rename(
        log2fc = logFC,
        avg_expr = AveExpr,
        t_statistic = t,
        p_value = P.Value,
        padj = adj.P.Val,
        b_statistic = B
    ) %>%
    select(gene, log2fc, avg_expr, t_statistic, p_value, padj, b_statistic) %>%
    arrange(padj)

cat(sprintf("  Tested: %d genes\n", nrow(results)))

# ---- 5. Significance thresholds --------------------------------------------

cat("\nApplying significance thresholds...\n")

padj_thresh <- 0.05
lfc_thresh <- 1.0

results <- results %>%
    mutate(
        significant = padj < padj_thresh & abs(log2fc) > lfc_thresh,
        direction = case_when(
            significant & log2fc > 0 ~ "Up in Tumour",
            significant & log2fc < 0 ~ "Down in Tumour",
            TRUE ~ "Not significant"
        )
    )

n_sig <- sum(results$significant, na.rm = TRUE)
n_up <- sum(results$direction == "Up in Tumour", na.rm = TRUE)
n_down <- sum(results$direction == "Down in Tumour", na.rm = TRUE)

cat(sprintf("  Significant (FDR<%.2f, |log2FC|>%.1f): %d\n",
            padj_thresh, lfc_thresh, n_sig))
cat(sprintf("  Up in Tumour:   %d\n", n_up))
cat(sprintf("  Down in Tumour: %d\n", n_down))

# ---- 6. ECM gene annotation ------------------------------------------------

cat("\nAnnotating ECM-associated genes...\n")

# Get Naba Matrisome gene sets from MSigDB
matrisome_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>%
    filter(grepl("NABA_", gs_name))

matrisome_genes <- unique(matrisome_sets$gene_symbol)

# Also get Hallmark EMT
emt_genes <- msigdbr(species = "Homo sapiens", category = "H") %>%
    filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>%
    pull(gene_symbol) %>%
    unique()

# Combine ECM-related gene lists
ecm_genes <- unique(c(matrisome_genes, emt_genes))

results <- results %>%
    mutate(
        is_ecm = gene %in% ecm_genes,
        ecm_and_significant = significant & is_ecm
    )

n_ecm_total <- sum(results$is_ecm)
n_ecm_sig <- sum(results$ecm_and_significant, na.rm = TRUE)
n_ecm_up <- sum(results$ecm_and_significant & results$direction == "Up in Tumour",
                 na.rm = TRUE)
n_ecm_down <- sum(results$ecm_and_significant & results$direction == "Down in Tumour",
                   na.rm = TRUE)

cat(sprintf("  ECM genes in dataset: %d\n", n_ecm_total))
cat(sprintf("  ECM genes significant: %d (%d up, %d down)\n",
            n_ecm_sig, n_ecm_up, n_ecm_down))

# ---- 7. Save results -------------------------------------------------------

cat("\nSaving results...\n")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("evidence/logs", showWarnings = FALSE, recursive = TRUE)

write.csv(results, "results/de_results_all.csv", row.names = FALSE)

sig_results <- results %>% filter(significant)
write.csv(sig_results, "results/de_results_significant.csv", row.names = FALSE)

ecm_results <- results %>% filter(is_ecm) %>% arrange(padj)
write.csv(ecm_results, "results/de_results_ecm_genes.csv", row.names = FALSE)

cat(sprintf("  Saved: results/de_results_all.csv (%d rows)\n", nrow(results)))
cat(sprintf("  Saved: results/de_results_significant.csv (%d rows)\n", nrow(sig_results)))
cat(sprintf("  Saved: results/de_results_ecm_genes.csv (%d rows)\n", nrow(ecm_results)))

# ---- 8. Volcano plot --------------------------------------------------------

cat("\nGenerating volcano plots...\n")

# Standard volcano plot
p_volcano <- ggplot(results, aes(x = log2fc, y = -log10(padj))) +
    geom_point(aes(colour = direction), size = 0.5, alpha = 0.5) +
    scale_colour_manual(
        values = c("Up in Tumour" = "#E74C3C",
                    "Down in Tumour" = "#3498DB",
                    "Not significant" = "grey70"),
        name = "Direction"
    ) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    labs(
        title = "Differential Expression: PDAC Tumour vs Normal",
        subtitle = sprintf(
            "%d significant (FDR<%.2f, |log2FC|>%.1f): %d up, %d down",
            n_sig, padj_thresh, lfc_thresh, n_up, n_down
        ),
        x = "log2(Fold Change)",
        y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

ggsave("figures/volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("  Saved: figures/volcano_plot.png\n")

# ECM-highlighted volcano plot
top_ecm <- results %>%
    filter(ecm_and_significant) %>%
    arrange(padj) %>%
    head(20)

p_ecm_volcano <- ggplot(results, aes(x = log2fc, y = -log10(padj))) +
    geom_point(
        data = results %>% filter(!is_ecm),
        colour = "grey80", size = 0.3, alpha = 0.4
    ) +
    geom_point(
        data = results %>% filter(is_ecm & !significant),
        colour = "grey50", size = 0.5, alpha = 0.6
    ) +
    geom_point(
        data = results %>% filter(ecm_and_significant),
        aes(colour = direction), size = 1.5, alpha = 0.8
    ) +
    scale_colour_manual(
        values = c("Up in Tumour" = "#E74C3C",
                    "Down in Tumour" = "#3498DB"),
        name = "Direction"
    ) +
    geom_text_repel(
        data = top_ecm, aes(label = gene),
        size = 2.5, max.overlaps = 25, segment.size = 0.2
    ) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    labs(
        title = "ECM-Associated Genes in PDAC Tumour vs Normal",
        subtitle = sprintf(
            "%d ECM genes significant: %d up, %d down (Matrisome + Hallmark EMT)",
            n_ecm_sig, n_ecm_up, n_ecm_down
        ),
        x = "log2(Fold Change)",
        y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

ggsave("figures/volcano_plot_ecm.png", p_ecm_volcano, width = 8, height = 6, dpi = 300)
cat("  Saved: figures/volcano_plot_ecm.png\n")

# ---- 9. Top ECM genes table ------------------------------------------------

cat("\n--- Top 15 ECM Genes (by adjusted p-value) ---\n")

top15 <- results %>%
    filter(ecm_and_significant) %>%
    arrange(padj) %>%
    head(15) %>%
    select(gene, log2fc, padj, direction)

print(as.data.frame(top15), row.names = FALSE)

# ---- 10. Evidence report ---------------------------------------------------

report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    analysis = "differential_expression",
    method = "limma_lmFit_eBayes",
    comparison = "Tumour_vs_Normal",
    samples = list(
        tumour = ncol(tumor),
        normal = ncol(normal),
        total = ncol(expr)
    ),
    filtering = list(
        genes_before = nrow(expr),
        genes_after = nrow(expr_filtered),
        genes_removed = nrow(expr) - nrow(expr_filtered),
        criterion = "log2 > 1 in >= 10% of samples"
    ),
    thresholds = list(
        padj = padj_thresh,
        log2fc = lfc_thresh
    ),
    results = list(
        genes_tested = nrow(results),
        significant_total = n_sig,
        up_in_tumour = n_up,
        down_in_tumour = n_down
    ),
    ecm_analysis = list(
        ecm_genes_in_dataset = n_ecm_total,
        ecm_significant = n_ecm_sig,
        ecm_up = n_ecm_up,
        ecm_down = n_ecm_down,
        gene_set_sources = c("MSigDB_NABA_MATRISOME", "HALLMARK_EMT"),
        top_ecm_genes = top15$gene
    ),
    r_version = R.version.string,
    limma_version = as.character(packageVersion("limma"))
)

write_json(report, "evidence/logs/de_analysis_report.json", pretty = TRUE, auto_unbox = TRUE)
cat("\n  Evidence: evidence/logs/de_analysis_report.json\n")

cat("\n=== Phase 4 Complete ===\n")
