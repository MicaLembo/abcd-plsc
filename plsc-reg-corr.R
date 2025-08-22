library("lmerTest")
library("readr")

data_dir <- "/Users/chloehampson/Desktop/projects/abcd-plsc/derivatives/none-reduced/regression/" # Make sure to leave the slash at the end
networks <- c("cgc-dt", "dt-dla", "dt-dt", "dt-vs", "vs-vs")
roi <- "rsfc"

# Level-1 Predictors
categorical_vars <- c("demo_sex_v2", "demo_prnt_gender_id_v2", "demo_origin_v2", "mri_info_manufacturer")
numerical_vars <- c("interview_age", "demo_prnt_age_v2", "demo_prnt_ed_v2_2yr_l", "demo_prtnr_ed_v2_2yr_l", "demo_comb_income_v2", "rsfmri_meanmotion")
phyhealth_vars <- c(
    "mctq_sdweek_calc",
    "mctq_msfsc_calc",
    "resp_wheeze_yn_y",
    "resp_pmcough_yn_y",
    "resp_diagnosis_yn_y",
    "resp_bronch_yn_y",
    "blood_pressure_sys_mean",
    "blood_pressure_dia_mean",
    "physical_activity1_y",
    "cbcl_scr_syn_internal_t",
    "cbcl_scr_syn_external_t"
)

phyhealth_cats <- c("resp_wheeze_yn_y", "resp_pmcough_yn_y", "resp_diagnosis_yn_y", "resp_bronch_yn_y")

for (network in networks) {
    data_path <- paste0(data_dir, "phyhealth_", network, "_data.csv")
    data <- read.table(file = data_path, sep = ",", header = TRUE)

    for (phyhealth_var in phyhealth_vars) {
        all_columns <- c(roi, categorical_vars, numerical_vars, phyhealth_var, "site_id_l", "rel_family_id")
        sub_data <- data[, all_columns]
        sub_data <- na.omit(sub_data)

        # Convert categorical variables to factors
        for (var in categorical_vars) {
            sub_data[[var]] <- factor(sub_data[[var]])
        }

        # Scale continuous predictors
        for (var in numerical_vars) {
            sub_data[[var]] <- scale(sub_data[[var]], center = TRUE, scale = TRUE)
        }

        if (phyhealth_var %in% phyhealth_cats) {
            sub_data[[phyhealth_var]] <- factor(sub_data[[phyhealth_var]])
        } else {
            sub_data[[phyhealth_var]] <- scale(sub_data[[phyhealth_var]],
                center = TRUE, scale = TRUE
            )
        }

        # Fixed effects formula
        fixed_effects <- paste(c(numerical_vars, categorical_vars), collapse = " + ")

        # Full model
        equation_lme <- paste(roi, "~", phyhealth_var, "+", fixed_effects, "+ (1|site_id_l/rel_family_id)") # "+ (1|site_id_l/rel_family_id)")

        # Run the full model
        model <- lmer(as.formula(equation_lme), data = sub_data)
        model_summary <- summary(model)
        print(model_summary)

        # Create formula for plotting
        plot_formula <- as.formula(paste(roi, "~ fitted(.) |", phyhealth_var))

        # Generate and print plot
        p <- plot(model, plot_formula, abline = c(0, 1))
        print(p)

        # Save results
        model_table <- as.data.frame(coef(summary(model)))
        out_file <- paste0(data_dir, "phyhealth_", network, "_", phyhealth_var, "_table.csv")
        write.csv(model_table, file = out_file, row.names = TRUE)


        # ===== Forest plot from existing LMM CSV outputs =====
# Uses files like: phyhealth_<network>_<predictor>_table.csv
# X-axis = beta (fixed effect) with 95% CI; Y-axis = predictor (with level if factor)

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)

# ---- match your current script ----
data_dir <- "/Users/chloehampson/Desktop/projects/abcd-plsc/derivatives/none-reduced/regression/"
networks <- c("cgc-dt", "dt-dla", "dt-dt", "dt-vs", "vs-vs")
phyhealth_vars <- c(
  "mctq_sdweek_calc","mctq_msfsc_calc",
  "resp_wheeze_yn_y","resp_pmcough_yn_y","resp_diagnosis_yn_y","resp_bronch_yn_y",
  "blood_pressure_sys_mean","blood_pressure_dia_mean",
  "physical_activity1_y",
  "cbcl_scr_syn_internal_t","cbcl_scr_syn_external_t"
)

# ---- helper: read one per-model CSV and pull predictor rows ----
read_model_table <- function(network, var) {
  f <- file.path(data_dir, paste0("phyhealth_", network, "_", var, "_table.csv"))
  if (!file.exists(f)) return(NULL)

  # try readr first, fall back to base if needed to capture rownames
  tab <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tab) || !any(names(tab) %in% c("Estimate","Std. Error"))) {
    tab <- tryCatch(as.data.frame(read.csv(f, check.names = FALSE)), error = function(e) NULL)
    if (is.null(tab)) return(NULL)
  }

  # recover term names (row names may be stored in X1 / blank / 'term')
  if ("term" %in% names(tab)) {
    tab$term <- tab$term
  } else if ("X1" %in% names(tab)) {
    tab$term <- tab$X1
  } else if ("" %in% names(tab)) {
    tab$term <- tab[[""]]
  } else {
    # final fallback: re-read with row.names=1
    tmp <- tryCatch(read.csv(f, row.names = 1, check.names = FALSE), error = function(e) NULL)
    if (is.null(tmp)) return(NULL)
    tab <- as.data.frame(tmp)
    tab$term <- rownames(tab)
  }

  est_col <- names(tab)[names(tab) %in% c("Estimate","estimate")][1]
  se_col  <- names(tab)[names(tab) %in% c("Std. Error","Std.Error","SE","std.error")][1]
  p_col   <- names(tab)[names(tab) %in% c("Pr(>|t|)","Pr(>|z|)","p.value","pvalue","p")][1]

  if (is.na(est_col) || is.na(se_col)) return(NULL)

  # keep predictor rows (handles factor levels like varLevel)
  tab <- tab %>% filter(str_starts(term, var))
  if (nrow(tab) == 0) return(NULL)

  z <- 1.96
  tibble(
    network = network,
    predictor = var,
    term = tab$term,
    beta = as.numeric(tab[[est_col]]),
    se   = as.numeric(tab[[se_col]]),
    ci_low  = beta - z*se,
    ci_high = beta + z*se,
    p_value = if (!is.na(p_col)) suppressWarnings(as.numeric(tab[[p_col]])) else NA_real_,
    label = ifelse(term == var, var, str_replace(term, paste0("^", var), paste0(var, " = ")))
  )
}

# ---- gather all model outputs ----
results <- map_dfr(networks, \(nw) map_dfr(phyhealth_vars, \(v) read_model_table(nw, v)))

# (optional) save combined table for transparency
write_csv(results, file.path(data_dir, "forestplot_all_results_from_existing_models.csv"))

# ---- keep only significant for poster (CI excludes 0 or p<.05 if present) ----
results_sig <- results %>%
  filter((ci_low > 0 | ci_high < 0) | (!is.na(p_value) & p_value < 0.05))

write_csv(results_sig, file.path(data_dir, "forestplot_significant_results_from_existing_models.csv"))

# ---- order and plot ----
if (nrow(results_sig) == 0) {
  message("No significant effects with current criteria.")
} else {
  plot_df <- results_sig %>%
    group_by(network) %>%
    arrange(beta, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(label = factor(label, levels = unique(label)))

  p <- ggplot(plot_df, aes(x = beta, y = label)) +
    geom_point(size = 2.8) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ network, scales = "free_y") +
    labs(
      title = "Significant LMM Associations with rsfc",
      subtitle = "β (fixed effect) with 95% CI; random intercepts for site/family; adjusted for covariates",
      x = "β (95% CI)",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  print(p)
    }
}
