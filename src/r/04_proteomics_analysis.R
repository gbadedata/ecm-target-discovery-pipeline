#!/usr/bin/env Rscript
# =============================================================================
# Phase 7: Proteomics Differential Abundance Analysis (limma)
#
# Analyses TMT proteomics data from CPTAC PDAC to identify differentially
# abundant proteins between tumour and normal tissue, with focus on
# ECM-associated proteins.
#
# Inputs:
#   data/raw/proteomics_tumor.cct
#   data/raw/proteomics_normal.cct
#
# Outputs:
#   results/proteomics_da_all.csv
#   results/proteomics_da_significant.csv
#   results/proteomics_da_ecm.csv
#   results/proteomics_missingness.csv
#   figures/proteomics_volcano.png
#   figures/proteomics_volcano_ecm.png
#   figures/proteomics_missingness.png
#   evidence/logs/proteomics_report.json
# =============================================================================

suppressPackageStartupMessages({
    library(limma)
    library(tidyverse)
    library(msigdbr)
    library(ggrepel)
    library(jsonlite)
})

set.seed(42)

cat("=== Phase 7: Proteomics Differential Abundance Analysis ===\n\n")

# ---- 1. Load data -----------------------------------------------------------

cat("Loading proteomics data...\n")

prot_tumor <- read.delim("data/raw/proteomics_tumor.cct",
                          row.names = 1, check.names = FALSE)
prot_normal <- read.delim("data/raw/proteomics_normal.cct",
                           row.names = 1, check.names = FALSE)

cat(sprintf("  Tumour: %d proteins x %d samples\n",
            nrow(prot_tumor), ncol(prot_tumor)))
cat(sprintf("  Normal: %d proteins x %d samples\n",
            nrow(prot_normal), ncol(prot_normal)))

# ---- 2. Missingness assessment ----------------------------------------------

cat("\nAssessing missingness...\n")

# Tumour missingness per protein
tumor_missing <- data.frame(
    protein = rownames(prot_tumor),
    tumor_missing_n = rowSums(is.na(prot_tumor)),
    tumor_missing_pct = rowMeans(is.na(prot_tumor)) * 100,
    tumor_detected_n = rowSums(!is.na(prot_tumor)),
    stringsAsFactors = FALSE
)

# Normal missingness per protein
normal_missing <- data.frame(
    protein = rownames(prot_normal),
    normal_missing_n = rowSums(is.na(prot_normal)),
    normal_missing_pct = rowMeans(is.na(prot_normal)) * 100,
    normal_detected_n = rowSums(!is.na(prot_normal)),
    stringsAsFactors = FALSE
)

missingness <- merge(tumor_missing, normal_missing, by = "protein")

overall_tumor_missing <- mean(is.na(as.matrix(prot_tumor))) * 100
overall_normal_missing <- mean(is.na(as.matrix(prot_normal))) * 100

cat(sprintf("  Overall tumour missingness: %.1f%%\n", overall_tumor_missing))
cat(sprintf("  Overall normal missingness: %.1f%%\n", overall_normal_missing))

# Proteins with complete data
complete_both <- sum(
    missingness$tumor_missing_n == 0 & missingness$normal_missing_n == 0
)
cat(sprintf("  Proteins with zero missing (both): %d\n", complete_both))

write.csv(missingness, "results/proteomics_missingness.csv",
          row.names = FALSE)

# Missingness histogram
p_missing <- ggplot(missingness,
                    aes(x = tumor_missing_pct)) +
    geom_histogram(bins = 50, fill = "#3498DB", alpha = 0.7,
                   colour = "white") +
    labs(
        title = "Protein Missingness Distribution (Tumour Samples)",
        subtitle = sprintf(
            "n=%d proteins, overall missingness=%.1f%%",
            nrow(prot_tumor), overall_tumor_missing
        ),
        x = "% Missing Across 140 Tumour Samples",
        y = "Number of Proteins"
    ) +
    theme_minimal(base_size = 12)

ggsave("figures/proteomics_missingness.png", p_missing,
       width = 8, height = 5, dpi = 300)
cat("  Saved: figures/proteomics_missingness.png\n")

# ---- 3. Filter proteins by detection rate -----------------------------------

cat("\nFiltering proteins...\n")

# Keep proteins detected in >= 50% of samples in BOTH groups
min_detect_tumor <- ceiling(0.5 * ncol(prot_tumor))
min_detect_normal <- ceiling(0.5 * ncol(prot_normal))

keep_tumor <- rowSums(!is.na(prot_tumor)) >= min_detect_tumor
keep_normal <- rowSums(!is.na(prot_normal)) >= min_detect_normal

# Find common samples between tumour and normal for paired analysis
common_samples <- intersect(colnames(prot_tumor), colnames(prot_normal))
cat(sprintf("  Common samples (tumour & normal): %d\n", length(common_samples)))

# For DA analysis: use all samples, keep proteins detected in either group
keep <- keep_tumor | keep_normal

prot_tumor_filt <- prot_tumor[keep, ]
prot_normal_filt <- prot_normal[keep, ]

cat(sprintf("  Before filtering: %d proteins\n", nrow(prot_tumor)))
cat(sprintf("  After filtering:  %d proteins\n", nrow(prot_tumor_filt)))
cat(sprintf("  Removed:          %d proteins\n",
            nrow(prot_tumor) - nrow(prot_tumor_filt)))

# ---- 4. Combine and impute --------------------------------------------------

cat("\nPreparing expression matrix...\n")

# For tumour vs normal comparison, use samples present in each group
# Normal samples may have different IDs than tumour samples
expr_combined <- cbind(prot_tumor_filt, prot_normal_filt)

group <- factor(
    c(rep("Tumour", ncol(prot_tumor_filt)),
      rep("Normal", ncol(prot_normal_filt))),
    levels = c("Normal", "Tumour")
)

cat(sprintf("  Combined: %d proteins x %d samples\n",
            nrow(expr_combined), ncol(expr_combined)))
cat(sprintf("  Missing values in combined: %.1f%%\n",
            mean(is.na(as.matrix(expr_combined))) * 100))

# Minimum probability imputation for remaining NAs
# Replace NA with a small value drawn from the low end of observed values
impute_minprob <- function(mat, q = 0.01) {
    mat_imputed <- mat
    for (j in seq_len(ncol(mat))) {
        col_vals <- mat[, j]
        na_idx <- is.na(col_vals)
        if (any(na_idx)) {
            observed <- col_vals[!na_idx]
            if (length(observed) > 0) {
                imp_mean <- quantile(observed, probs = q, na.rm = TRUE)
                imp_sd <- sd(observed, na.rm = TRUE) * 0.3
                n_impute <- sum(na_idx)
                mat_imputed[na_idx, j] <- rnorm(n_impute,
                                                  mean = imp_mean,
                                                  sd = max(imp_sd, 0.1))
            }
        }
    }
    return(mat_imputed)
}

expr_imputed <- impute_minprob(as.matrix(expr_combined))

cat(sprintf("  After imputation: %.1f%% missing\n",
            mean(is.na(expr_imputed)) * 100))

# ---- 5. limma differential abundance ----------------------------------------

cat("\nRunning limma for differential abundance...\n")

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "TumourVsNormal")

fit <- lmFit(expr_imputed, design)
fit <- eBayes(fit)

da_results <- topTable(fit, coef = "TumourVsNormal",
                        number = Inf, sort.by = "none")

da_results$protein <- rownames(da_results)
da_results <- da_results %>%
    rename(
        log2fc = logFC,
        avg_abundance = AveExpr,
        t_statistic = t,
        p_value = P.Value,
        padj = adj.P.Val,
        b_statistic = B
    ) %>%
    select(protein, log2fc, avg_abundance, t_statistic,
           p_value, padj, b_statistic) %>%
    arrange(padj)

# Add detection info
da_results <- da_results %>%
    left_join(
        missingness %>% select(protein, tumor_missing_pct, normal_missing_pct),
        by = "protein"
    )

# Significance
padj_thresh <- 0.05
lfc_thresh <- 0.5  # Lower threshold for proteomics

da_results <- da_results %>%
    mutate(
        significant = padj < padj_thresh & abs(log2fc) > lfc_thresh,
        direction = case_when(
            significant & log2fc > 0 ~ "Up in Tumour",
            significant & log2fc < 0 ~ "Down in Tumour",
            TRUE ~ "Not significant"
        )
    )

n_tested <- nrow(da_results)
n_sig <- sum(da_results$significant, na.rm = TRUE)
n_up <- sum(da_results$direction == "Up in Tumour", na.rm = TRUE)
n_down <- sum(da_results$direction == "Down in Tumour", na.rm = TRUE)

cat(sprintf("  Tested: %d proteins\n", n_tested))
cat(sprintf("  Significant (FDR<%.2f, |log2FC|>%.1f): %d\n",
            padj_thresh, lfc_thresh, n_sig))
cat(sprintf("  Up in Tumour:   %d\n", n_up))
cat(sprintf("  Down in Tumour: %d\n", n_down))

# ---- 6. ECM protein annotation ----------------------------------------------

cat("\nAnnotating ECM-associated proteins...\n")

naba <- msigdbr(species = "Homo sapiens") %>%
    filter(grepl("NABA_", gs_name))
ecm_genes <- unique(naba$gene_symbol)

emt <- msigdbr(species = "Homo sapiens") %>%
    filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
ecm_genes <- unique(c(ecm_genes, emt$gene_symbol))

da_results <- da_results %>%
    mutate(
        is_ecm = protein %in% ecm_genes,
        ecm_and_significant = significant & is_ecm
    )

n_ecm <- sum(da_results$is_ecm)
n_ecm_sig <- sum(da_results$ecm_and_significant, na.rm = TRUE)
n_ecm_up <- sum(da_results$ecm_and_significant &
                da_results$direction == "Up in Tumour", na.rm = TRUE)
n_ecm_down <- sum(da_results$ecm_and_significant &
                  da_results$direction == "Down in Tumour", na.rm = TRUE)

cat(sprintf("  ECM proteins in dataset: %d\n", n_ecm))
cat(sprintf("  ECM significant: %d (%d up, %d down)\n",
            n_ecm_sig, n_ecm_up, n_ecm_down))

# ---- 7. Save results --------------------------------------------------------

cat("\nSaving results...\n")

write.csv(da_results, "results/proteomics_da_all.csv", row.names = FALSE)

sig_results <- da_results %>% filter(significant)
write.csv(sig_results, "results/proteomics_da_significant.csv",
          row.names = FALSE)

ecm_prot <- da_results %>% filter(is_ecm) %>% arrange(padj)
write.csv(ecm_prot, "results/proteomics_da_ecm.csv", row.names = FALSE)

cat(sprintf("  Saved: results/proteomics_da_all.csv (%d rows)\n",
            nrow(da_results)))
cat(sprintf("  Saved: results/proteomics_da_significant.csv (%d rows)\n",
            nrow(sig_results)))
cat(sprintf("  Saved: results/proteomics_da_ecm.csv (%d rows)\n",
            nrow(ecm_prot)))

# ---- 8. Volcano plots -------------------------------------------------------

cat("\nGenerating volcano plots...\n")

p_volcano <- ggplot(da_results,
                    aes(x = log2fc, y = -log10(padj))) +
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
        title = "Proteomics: PDAC Tumour vs Normal",
        subtitle = sprintf(
            "%d significant (FDR<%.2f, |log2FC|>%.1f): %d up, %d down",
            n_sig, padj_thresh, lfc_thresh, n_up, n_down
        ),
        x = "log2(Fold Change)",
        y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

ggsave("figures/proteomics_volcano.png", p_volcano,
       width = 8, height = 6, dpi = 300)
cat("  Saved: figures/proteomics_volcano.png\n")

# ECM-highlighted volcano
top_ecm_prot <- da_results %>%
    filter(ecm_and_significant) %>%
    arrange(padj) %>%
    head(20)

p_ecm_volcano <- ggplot(da_results,
                        aes(x = log2fc, y = -log10(padj))) +
    geom_point(
        data = da_results %>% filter(!is_ecm),
        colour = "grey80", size = 0.3, alpha = 0.4
    ) +
    geom_point(
        data = da_results %>% filter(is_ecm & !significant),
        colour = "grey50", size = 0.5, alpha = 0.6
    ) +
    geom_point(
        data = da_results %>% filter(ecm_and_significant),
        aes(colour = direction), size = 1.5, alpha = 0.8
    ) +
    scale_colour_manual(
        values = c("Up in Tumour" = "#E74C3C",
                    "Down in Tumour" = "#3498DB"),
        name = "Direction"
    ) +
    geom_text_repel(
        data = top_ecm_prot, aes(label = protein),
        size = 2.5, max.overlaps = 25, segment.size = 0.2
    ) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    labs(
        title = "ECM Proteins in PDAC Tumour vs Normal",
        subtitle = sprintf(
            "%d ECM proteins significant: %d up, %d down",
            n_ecm_sig, n_ecm_up, n_ecm_down
        ),
        x = "log2(Fold Change)",
        y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

ggsave("figures/proteomics_volcano_ecm.png", p_ecm_volcano,
       width = 8, height = 6, dpi = 300)
cat("  Saved: figures/proteomics_volcano_ecm.png\n")

# ---- 9. Top ECM proteins table ----------------------------------------------

cat("\n--- Top 15 ECM Proteins (by adjusted p-value) ---\n")

top15_prot <- da_results %>%
    filter(ecm_and_significant) %>%
    arrange(padj) %>%
    head(15) %>%
    select(protein, log2fc, padj, direction)

print(as.data.frame(top15_prot), row.names = FALSE)

# ---- 10. Evidence report ----------------------------------------------------

report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    analysis = "proteomics_differential_abundance",
    method = "limma_with_minprob_imputation",
    comparison = "Tumour_vs_Normal",
    samples = list(
        tumour = ncol(prot_tumor),
        normal = ncol(prot_normal)
    ),
    missingness = list(
        tumour_overall_pct = round(overall_tumor_missing, 2),
        normal_overall_pct = round(overall_normal_missing, 2),
        proteins_complete_both = complete_both,
        imputation_method = "MinProb (q=0.01, sd_factor=0.3)"
    ),
    filtering = list(
        proteins_before = nrow(prot_tumor),
        proteins_after = nrow(prot_tumor_filt),
        criterion = "Detected in >= 50% of samples in either group"
    ),
    thresholds = list(
        padj = padj_thresh,
        log2fc = lfc_thresh
    ),
    results = list(
        proteins_tested = n_tested,
        significant_total = n_sig,
        up_in_tumour = n_up,
        down_in_tumour = n_down
    ),
    ecm_analysis = list(
        ecm_proteins_in_dataset = n_ecm,
        ecm_significant = n_ecm_sig,
        ecm_up = n_ecm_up,
        ecm_down = n_ecm_down,
        top_ecm_proteins = top15_prot$protein
    ),
    r_version = R.version.string,
    limma_version = as.character(packageVersion("limma"))
)

write_json(report, "evidence/logs/proteomics_report.json",
           pretty = TRUE, auto_unbox = TRUE)
cat("\n  Evidence: evidence/logs/proteomics_report.json\n")

cat("\n=== Phase 7 Complete ===\n")
