#!/usr/bin/env Rscript
# =============================================================================
# Phase 5: ECM Signature Scoring (ssGSEA via GSVA)
#
# Computes per-sample ECM activity scores using ssGSEA across curated
# gene sets, then validates against pathologist-estimated stromal fraction.
#
# Inputs:
#   data/raw/rnaseq_tumor.cct
#   data/processed/clinical_clean.csv
#
# Outputs:
#   results/ecm_scores.csv
#   figures/ecm_score_vs_stromal_fraction.png
#   figures/ecm_score_distribution.png
#   figures/ecm_heatmap_top_genes.png
#   evidence/logs/ecm_scoring_report.json
# =============================================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(msigdbr)
    library(tidyverse)
    library(pheatmap)
    library(jsonlite)
})

set.seed(42)

cat("=== Phase 5: ECM Signature Scoring ===\n\n")

# ---- 1. Load data -----------------------------------------------------------

cat("Loading expression data...\n")
tumor <- read.delim("data/raw/rnaseq_tumor.cct",
                     row.names = 1, check.names = FALSE)
cat(sprintf("  Tumour: %d genes x %d samples\n", nrow(tumor), ncol(tumor)))

cat("Loading clinical data...\n")
clinical <- read.csv("data/processed/clinical_clean.csv")
cat(sprintf("  Clinical: %d patients\n", nrow(clinical)))

# ---- 2. Build ECM gene sets -------------------------------------------------

cat("\nBuilding ECM gene sets from MSigDB...\n")

# Naba Matrisome collections
naba <- msigdbr(species = "Homo sapiens") %>%
    filter(grepl("NABA_", gs_name))

naba_sets <- split(naba$gene_symbol, naba$gs_name)

# Hallmark EMT
emt <- msigdbr(species = "Homo sapiens") %>%
    filter(gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
emt_set <- list(HALLMARK_EMT = unique(emt$gene_symbol))

# Reactome ECM
reactome_ecm <- msigdbr(species = "Homo sapiens") %>%
    filter(gs_name %in% c(
        "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
        "REACTOME_COLLAGEN_FORMATION",
        "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
        "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES"
    ))
reactome_sets <- split(reactome_ecm$gene_symbol, reactome_ecm$gs_name)

# Combine all gene sets
all_gene_sets <- c(naba_sets, emt_set, reactome_sets)

# Report gene set sizes
cat(sprintf("  Gene sets loaded: %d\n", length(all_gene_sets)))
for (name in names(all_gene_sets)) {
    n_in_data <- sum(all_gene_sets[[name]] %in% rownames(tumor))
    cat(sprintf("    %s: %d genes (%d in data)\n",
                name, length(all_gene_sets[[name]]), n_in_data))
}

# Filter to sets with at least 10 genes in the data
gene_sets_filtered <- all_gene_sets[
    sapply(all_gene_sets, function(gs) sum(gs %in% rownames(tumor)) >= 10)
]
cat(sprintf("  Gene sets after filtering (>=10 genes in data): %d\n",
            length(gene_sets_filtered)))

# ---- 3. Run ssGSEA ----------------------------------------------------------

cat("\nRunning ssGSEA via GSVA (this may take a minute)...\n")

expr_matrix <- as.matrix(tumor)

gsva_param <- ssgseaParam(
    exprData = expr_matrix,
    geneSets = gene_sets_filtered,
    normalize = TRUE
)
ssgsea_scores <- gsva(gsva_param, verbose = FALSE)

cat(sprintf("  ssGSEA complete: %d gene sets x %d samples\n",
            nrow(ssgsea_scores), ncol(ssgsea_scores)))

# ---- 4. Create composite ECM score -----------------------------------------

cat("\nComputing composite ECM score...\n")

# Use core matrisome set as the primary ECM score
core_matrisome_name <- grep("NABA_CORE_MATRISOME$", rownames(ssgsea_scores),
                            value = TRUE)

if (length(core_matrisome_name) == 1) {
    ecm_primary <- ssgsea_scores[core_matrisome_name, ]
    cat(sprintf("  Primary ECM score: %s\n", core_matrisome_name))
} else {
    # Fallback: average of all NABA sets
    naba_rows <- grep("NABA_", rownames(ssgsea_scores))
    ecm_primary <- colMeans(ssgsea_scores[naba_rows, , drop = FALSE])
    cat("  Primary ECM score: average of all NABA sets\n")
}

# Build scores data frame
scores_df <- data.frame(
    case_id = colnames(ssgsea_scores),
    ecm_score = as.numeric(ecm_primary),
    stringsAsFactors = FALSE
)

# Add all individual gene set scores
for (i in seq_len(nrow(ssgsea_scores))) {
    col_name <- gsub("[^A-Za-z0-9_]", "_", rownames(ssgsea_scores)[i])
    col_name <- tolower(col_name)
    scores_df[[col_name]] <- as.numeric(ssgsea_scores[i, ])
}

# ---- 5. Merge with clinical data -------------------------------------------

cat("\nMerging with clinical data...\n")

scores_clinical <- scores_df %>%
    left_join(clinical, by = "case_id")

matched <- sum(!is.na(scores_clinical$fraction_stromal))
cat(sprintf("  Matched with stromal fraction: %d / %d\n",
            matched, nrow(scores_clinical)))

# Create ECM high/low groups (median split)
ecm_median <- median(scores_clinical$ecm_score)
scores_clinical <- scores_clinical %>%
    mutate(
        ecm_group = ifelse(ecm_score >= ecm_median, "ECM_high", "ECM_low"),
        ecm_group = factor(ecm_group, levels = c("ECM_low", "ECM_high"))
    )

cat(sprintf("  ECM median score: %.4f\n", ecm_median))
cat(sprintf("  ECM_high: %d, ECM_low: %d\n",
            sum(scores_clinical$ecm_group == "ECM_high"),
            sum(scores_clinical$ecm_group == "ECM_low")))

# ---- 6. Validate against stromal fraction -----------------------------------

cat("\nValidating ECM score vs pathologist stromal fraction...\n")

valid_pairs <- scores_clinical %>%
    filter(!is.na(fraction_stromal) & !is.na(ecm_score))

if (nrow(valid_pairs) >= 10) {
    cor_test <- cor.test(valid_pairs$ecm_score, valid_pairs$fraction_stromal,
                         method = "spearman")
    rho <- cor_test$estimate
    p_val <- cor_test$p.value

    cat(sprintf("  Spearman rho: %.4f\n", rho))
    cat(sprintf("  P-value: %.2e\n", p_val))
    cat(sprintf("  Interpretation: %s\n",
                ifelse(rho > 0.3 & p_val < 0.05,
                       "VALIDATED - significant positive correlation",
                       ifelse(rho > 0 & p_val < 0.05,
                              "Weak but significant positive correlation",
                              "No significant correlation"))))

    # Scatter plot: ECM score vs stromal fraction
    p_validation <- ggplot(valid_pairs,
                           aes(x = fraction_stromal, y = ecm_score)) +
        geom_point(aes(colour = ecm_group), size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, colour = "grey30",
                    linewidth = 0.8, linetype = "dashed") +
        scale_colour_manual(
            values = c("ECM_low" = "#3498DB", "ECM_high" = "#E74C3C"),
            name = "ECM Group"
        ) +
        labs(
            title = "ECM Score vs Pathologist Stromal Fraction",
            subtitle = sprintf(
                "Spearman rho=%.3f, p=%.2e (n=%d)",
                rho, p_val, nrow(valid_pairs)
            ),
            x = "Pathologist Stromal Fraction",
            y = "ssGSEA ECM Score (Core Matrisome)"
        ) +
        theme_minimal(base_size = 12)

    ggsave("figures/ecm_score_vs_stromal_fraction.png",
           p_validation, width = 8, height = 6, dpi = 300)
    cat("  Saved: figures/ecm_score_vs_stromal_fraction.png\n")
} else {
    rho <- NA
    p_val <- NA
    cat("  Insufficient paired data for correlation\n")
}

# ---- 7. ECM score distribution plot -----------------------------------------

p_dist <- ggplot(scores_clinical, aes(x = ecm_score, fill = ecm_group)) +
    geom_histogram(bins = 30, alpha = 0.7, colour = "white") +
    scale_fill_manual(
        values = c("ECM_low" = "#3498DB", "ECM_high" = "#E74C3C"),
        name = "ECM Group"
    ) +
    geom_vline(xintercept = ecm_median, linetype = "dashed", colour = "grey30") +
    labs(
        title = "Distribution of ECM Scores (CPTAC PDAC, n=140)",
        subtitle = "ssGSEA Core Matrisome score, median split",
        x = "ECM Score",
        y = "Count"
    ) +
    theme_minimal(base_size = 12)

ggsave("figures/ecm_score_distribution.png", p_dist,
       width = 8, height = 5, dpi = 300)
cat("  Saved: figures/ecm_score_distribution.png\n")

# ---- 8. Heatmap of top ECM genes --------------------------------------------

cat("\nGenerating ECM gene heatmap...\n")

# Get core matrisome genes
core_genes <- gene_sets_filtered[[core_matrisome_name]]
core_in_data <- core_genes[core_genes %in% rownames(tumor)]

# Select top variable core matrisome genes
gene_vars <- apply(tumor[core_in_data, ], 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]

# Order samples by ECM score
sample_order <- scores_clinical %>%
    arrange(ecm_score) %>%
    pull(case_id)

# Annotation bar
annot_df <- scores_clinical %>%
    select(case_id, ecm_group) %>%
    column_to_rownames("case_id")

annot_colours <- list(ecm_group = c(ECM_low = "#3498DB", ECM_high = "#E74C3C"))

# Scale expression for heatmap
heatmap_data <- t(scale(t(as.matrix(tumor[top_var_genes, sample_order]))))
heatmap_data[heatmap_data > 3] <- 3
heatmap_data[heatmap_data < -3] <- -3

png("figures/ecm_heatmap_top_genes.png", width = 12, height = 8,
    units = "in", res = 300)
pheatmap(
    heatmap_data,
    annotation_col = annot_df[sample_order, , drop = FALSE],
    annotation_colors = annot_colours,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    fontsize_row = 6,
    main = "Top 50 Variable Core Matrisome Genes\n(samples ordered by ECM score)"
)
dev.off()
cat("  Saved: figures/ecm_heatmap_top_genes.png\n")

# ---- 9. Save results -------------------------------------------------------

cat("\nSaving results...\n")

# Save scores with clinical data
output_cols <- c("case_id", "ecm_score", "ecm_group",
                 grep("^naba_|^hallmark_|^reactome_", names(scores_clinical),
                      value = TRUE))
# Also include key clinical columns if present
clinical_cols <- c("fraction_stromal", "stromal_group",
                   "survival_days", "survival_event", "age", "sex",
                   "stage_overall", "cellularity_neoplastic")
output_cols <- c(output_cols,
                 intersect(clinical_cols, names(scores_clinical)))

write.csv(scores_clinical[, output_cols, drop = FALSE],
          "results/ecm_scores.csv", row.names = FALSE)
cat(sprintf("  Saved: results/ecm_scores.csv (%d rows)\n",
            nrow(scores_clinical)))

# ---- 10. Evidence report ---------------------------------------------------

report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    analysis = "ecm_signature_scoring",
    method = "ssGSEA_via_GSVA",
    primary_score = core_matrisome_name,
    gene_sets = list(
        total_loaded = length(all_gene_sets),
        after_filtering = length(gene_sets_filtered),
        names = names(gene_sets_filtered)
    ),
    ecm_score_summary = list(
        mean = round(mean(scores_clinical$ecm_score), 4),
        median = round(ecm_median, 4),
        sd = round(sd(scores_clinical$ecm_score), 4),
        min = round(min(scores_clinical$ecm_score), 4),
        max = round(max(scores_clinical$ecm_score), 4)
    ),
    ecm_groups = list(
        ECM_high = sum(scores_clinical$ecm_group == "ECM_high"),
        ECM_low = sum(scores_clinical$ecm_group == "ECM_low")
    ),
    validation = list(
        spearman_rho = ifelse(is.na(rho), NA, round(as.numeric(rho), 4)),
        p_value = ifelse(is.na(p_val), NA, p_val),
        n_paired = nrow(valid_pairs),
        interpretation = ifelse(
            !is.na(rho) && rho > 0.3 && p_val < 0.05,
            "Validated: significant positive correlation with stromal fraction",
            "See report for details"
        )
    ),
    r_version = R.version.string,
    gsva_version = as.character(packageVersion("GSVA"))
)

write_json(report, "evidence/logs/ecm_scoring_report.json",
           pretty = TRUE, auto_unbox = TRUE)
cat("  Evidence: evidence/logs/ecm_scoring_report.json\n")

cat("\n=== Phase 5 Complete ===\n")
