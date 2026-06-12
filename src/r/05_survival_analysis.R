#!/usr/bin/env Rscript
# =============================================================================
# Phase 9: Clinical Association and Survival Analysis
#
# Tests whether ECM-high patients have worse clinical outcomes using
# Kaplan-Meier survival analysis and Cox proportional hazards regression.
#
# Inputs:
#   results/ecm_scores.csv (from Phase 5, includes clinical data)
#
# Outputs:
#   results/survival_results.csv
#   figures/km_ecm_groups.png
#   figures/km_stromal_groups.png
#   figures/ecm_clinical_associations.png
#   evidence/logs/survival_report.json
# =============================================================================

suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(tidyverse)
    library(jsonlite)
})

set.seed(42)

cat("=== Phase 9: Clinical Association and Survival Analysis ===\n\n")

# ---- 1. Load data -----------------------------------------------------------

cat("Loading ECM scores with clinical data...\n")
scores <- read.csv("results/ecm_scores.csv")
cat(sprintf("  Patients: %d\n", nrow(scores)))
cat(sprintf("  ECM_high: %d, ECM_low: %d\n",
            sum(scores$ecm_group == "ECM_high", na.rm = TRUE),
            sum(scores$ecm_group == "ECM_low", na.rm = TRUE)))

# Prepare survival data
surv_data <- scores %>%
    filter(!is.na(survival_days) & !is.na(survival_event)) %>%
    mutate(
        survival_months = survival_days / 30.44,
        ecm_group = factor(ecm_group, levels = c("ECM_low", "ECM_high"))
    )

cat(sprintf("  Patients with survival data: %d\n", nrow(surv_data)))
cat(sprintf("  Events (deaths): %d\n", sum(surv_data$survival_event)))
cat(sprintf("  Censored: %d\n", sum(surv_data$survival_event == 0)))
cat(sprintf("  Median follow-up: %.1f months\n",
            median(surv_data$survival_months)))

# ---- 2. Kaplan-Meier: ECM groups -------------------------------------------

cat("\n--- Kaplan-Meier: ECM Score Groups ---\n")

surv_obj <- Surv(surv_data$survival_months, surv_data$survival_event)

km_ecm <- survfit(surv_obj ~ ecm_group, data = surv_data)

# Log-rank test
logrank_ecm <- survdiff(surv_obj ~ ecm_group, data = surv_data)
logrank_p <- 1 - pchisq(logrank_ecm$chisq, df = 1)

# Median survival per group
km_summary <- summary(km_ecm)$table
cat("  Median survival (months):\n")
print(km_summary[, c("records", "events", "median",
                      "0.95LCL", "0.95UCL")])
cat(sprintf("\n  Log-rank p-value: %.4f\n", logrank_p))
cat(sprintf("  Interpretation: %s\n",
            ifelse(logrank_p < 0.05,
                   "SIGNIFICANT - ECM groups differ in survival",
                   "Not significant at p<0.05")))

# KM plot
p_km_ecm <- ggsurvplot(
    km_ecm,
    data = surv_data,
    pval = TRUE,
    pval.method = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    palette = c("#3498DB", "#E74C3C"),
    legend.labs = c("ECM-low", "ECM-high"),
    legend.title = "ECM Group",
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = "Overall Survival by ECM Score Group (CPTAC PDAC)",
    ggtheme = theme_minimal(base_size = 12)
)

ggsave("figures/km_ecm_groups.png",
       plot = print(p_km_ecm),
       width = 8, height = 7, dpi = 300)
cat("  Saved: figures/km_ecm_groups.png\n")

# ---- 3. Kaplan-Meier: Stromal groups ---------------------------------------

cat("\n--- Kaplan-Meier: Stromal Fraction Groups ---\n")

if ("stromal_group" %in% names(surv_data)) {
    surv_stromal <- surv_data %>%
        filter(!is.na(stromal_group)) %>%
        mutate(stromal_group = factor(stromal_group,
                                       levels = c("stromal_low", "stromal_high")))

    if (nrow(surv_stromal) >= 20) {
        km_stromal <- survfit(
            Surv(survival_months, survival_event) ~ stromal_group,
            data = surv_stromal
        )

        logrank_stromal <- survdiff(
            Surv(survival_months, survival_event) ~ stromal_group,
            data = surv_stromal
        )
        stromal_p <- 1 - pchisq(logrank_stromal$chisq, df = 1)

        cat(sprintf("  Log-rank p-value: %.4f\n", stromal_p))

        p_km_stromal <- ggsurvplot(
            km_stromal,
            data = surv_stromal,
            pval = TRUE,
            pval.method = TRUE,
            risk.table = TRUE,
            risk.table.col = "strata",
            palette = c("#3498DB", "#E74C3C"),
            legend.labs = c("Stromal-low", "Stromal-high"),
            legend.title = "Stromal Group",
            xlab = "Time (months)",
            ylab = "Overall Survival Probability",
            title = "Overall Survival by Pathologist Stromal Fraction (CPTAC PDAC)",
            ggtheme = theme_minimal(base_size = 12)
        )

        ggsave("figures/km_stromal_groups.png",
               plot = print(p_km_stromal),
               width = 8, height = 7, dpi = 300)
        cat("  Saved: figures/km_stromal_groups.png\n")
    } else {
        stromal_p <- NA
        cat("  Insufficient data for stromal group analysis\n")
    }
} else {
    stromal_p <- NA
    cat("  Stromal group not available\n")
}

# ---- 4. Cox Proportional Hazards -------------------------------------------

cat("\n--- Cox Proportional Hazards Regression ---\n")

# Univariate: ECM score (continuous)
cox_ecm <- coxph(surv_obj ~ ecm_score, data = surv_data)
cox_ecm_summary <- summary(cox_ecm)

cat("\n  Univariate: ECM score (continuous)\n")
cat(sprintf("    HR: %.3f (95%% CI: %.3f-%.3f)\n",
            exp(coef(cox_ecm)),
            exp(confint(cox_ecm))[1],
            exp(confint(cox_ecm))[2]))
cat(sprintf("    P-value: %.4f\n", cox_ecm_summary$coefficients[, "Pr(>|z|)"]))

# Multivariate: adjust for age and stage
if ("age" %in% names(surv_data) & "stage_overall" %in% names(surv_data)) {
    surv_multi <- surv_data %>%
        filter(!is.na(age) & !is.na(stage_overall))

    if (nrow(surv_multi) >= 30) {
        cox_multi <- coxph(
            Surv(survival_months, survival_event) ~
                ecm_score + age + stage_overall,
            data = surv_multi
        )
        cox_multi_summary <- summary(cox_multi)

        cat("\n  Multivariate: ECM score + age + stage\n")
        cat(sprintf("    N: %d\n", nrow(surv_multi)))

        cox_table <- as.data.frame(cox_multi_summary$coefficients)
        cox_table$HR <- exp(cox_table$coef)
        cox_ci <- exp(confint(cox_multi))
        cox_table$CI_lower <- cox_ci[, 1]
        cox_table$CI_upper <- cox_ci[, 2]

        for (var in rownames(cox_table)) {
            cat(sprintf("    %s: HR=%.3f (%.3f-%.3f), p=%.4f\n",
                        var,
                        cox_table[var, "HR"],
                        cox_table[var, "CI_lower"],
                        cox_table[var, "CI_upper"],
                        cox_table[var, "Pr(>|z|)"]))
        }
    } else {
        cox_multi_summary <- NULL
        cat("  Insufficient data for multivariate model\n")
    }
} else {
    cox_multi_summary <- NULL
    cat("  Age or stage not available for multivariate model\n")
}

# ---- 5. Clinical associations -----------------------------------------------

cat("\n--- ECM Score vs Clinical Features ---\n")

# ECM score by stage
if ("stage_overall" %in% names(surv_data)) {
    stage_summary <- surv_data %>%
        filter(!is.na(stage_overall)) %>%
        group_by(stage_overall) %>%
        summarise(
            n = n(),
            ecm_mean = round(mean(ecm_score), 4),
            ecm_sd = round(sd(ecm_score), 4),
            .groups = "drop"
        )
    cat("\n  ECM score by stage:\n")
    print(as.data.frame(stage_summary), row.names = FALSE)
}

# ECM score by sex
if ("sex" %in% names(surv_data)) {
    sex_test <- wilcox.test(ecm_score ~ sex, data = surv_data)
    cat(sprintf("\n  ECM score by sex: Wilcoxon p=%.4f\n", sex_test$p.value))
}

# ECM score vs age correlation
if ("age" %in% names(surv_data)) {
    age_cor <- cor.test(surv_data$ecm_score, surv_data$age,
                        method = "spearman")
    cat(sprintf("  ECM score vs age: rho=%.3f, p=%.4f\n",
                age_cor$estimate, age_cor$p.value))
}

# ---- 6. Save results -------------------------------------------------------

cat("\nSaving results...\n")

surv_results <- surv_data %>%
    select(case_id, ecm_score, ecm_group,
           survival_days, survival_months, survival_event,
           any_of(c("age", "sex", "stage_overall",
                     "fraction_stromal", "stromal_group")))

write.csv(surv_results, "results/survival_results.csv", row.names = FALSE)
cat(sprintf("  Saved: results/survival_results.csv (%d rows)\n",
            nrow(surv_results)))

# ---- 7. Evidence report ----------------------------------------------------

report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    analysis = "survival_analysis",
    patients_with_survival = nrow(surv_data),
    events = sum(surv_data$survival_event),
    censored = sum(surv_data$survival_event == 0),
    median_followup_months = round(median(surv_data$survival_months), 1),
    kaplan_meier_ecm = list(
        logrank_p = round(logrank_p, 6),
        significant = logrank_p < 0.05
    ),
    kaplan_meier_stromal = list(
        logrank_p = ifelse(is.na(stromal_p), NA, round(stromal_p, 6))
    ),
    cox_univariate_ecm = list(
        hr = round(as.numeric(exp(coef(cox_ecm))), 3),
        ci_lower = round(as.numeric(exp(confint(cox_ecm))[1]), 3),
        ci_upper = round(as.numeric(exp(confint(cox_ecm))[2]), 3),
        p_value = round(
            as.numeric(cox_ecm_summary$coefficients[, "Pr(>|z|)"]), 6
        )
    ),
    r_version = R.version.string,
    survival_version = as.character(packageVersion("survival"))
)

write_json(report, "evidence/logs/survival_report.json",
           pretty = TRUE, auto_unbox = TRUE)
cat("  Evidence: evidence/logs/survival_report.json\n")

cat("\n=== Phase 9 Complete ===\n")
