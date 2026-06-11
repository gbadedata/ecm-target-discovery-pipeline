#!/usr/bin/env Rscript
# install_r_packages.R
# Installs all R packages required for the ECM Target Discovery Pipeline.
# Run once after setting up R: Rscript install_r_packages.R

cat("Installing R packages for ECM Target Discovery Pipeline\n")
cat("This may take 15-30 minutes on first run.\n\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
    # Data handling
    "tidyverse",
    "data.table",
    "readxl",
    "writexl",
    "jsonlite",
    "yaml",

    # Statistics
    "survival",
    "survminer",
    "broom",
    "ggfortify",

    # Visualisation
    "ggplot2",
    "pheatmap",
    "ggrepel",
    "patchwork",
    "RColorBrewer",
    "viridis",
    "scales",

    # ML and clustering
    "factoextra",
    "cluster",
    "umap",
    "Rtsne",
    "glmnet",
    "randomForest",
    "caret",

    # Reporting
    "rmarkdown",
    "knitr",

    # Utilities
    "here",
    "fs",
    "janitor",
    "testthat",
    "optparse"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  Installing %s...\n", pkg))
        install.packages(pkg, quiet = TRUE)
    } else {
        cat(sprintf("  %s already installed\n", pkg))
    }
}

# Bioconductor packages
bioc_packages <- c(
    # Core transcriptomics
    "DESeq2",
    "edgeR",
    "limma",

    # Annotation
    "org.Hs.eg.db",
    "AnnotationDbi",
    "biomaRt",

    # Pathway analysis
    "clusterProfiler",
    "enrichplot",
    "msigdbr",
    "GSVA",
    "fgsea",

    # Data structures
    "SummarizedExperiment",
    "GenomicRanges",
    "S4Vectors",
    "BiocGenerics",

    # Immune and TME estimation
    "ESTIMATE",
    "MCPcounter"
)

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  Installing %s...\n", pkg))
        tryCatch(
            BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE),
            error = function(e) {
                cat(sprintf("  WARNING: Failed to install %s: %s\n", pkg, e$message))
            }
        )
    } else {
        cat(sprintf("  %s already installed\n", pkg))
    }
}

# Verification
cat("\n--- Package verification ---\n")
all_packages <- c(cran_packages, bioc_packages)
installed <- vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
missing <- all_packages[!installed]

if (length(missing) == 0) {
    cat("All packages installed successfully.\n")
} else {
    cat(sprintf("Missing packages (%d): %s\n", length(missing), paste(missing, collapse = ", ")))
    cat("You may need to install system dependencies and retry.\n")
}

cat("\nR session info:\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("BiocManager version: %s\n", as.character(packageVersion("BiocManager"))))
cat("Done.\n")
