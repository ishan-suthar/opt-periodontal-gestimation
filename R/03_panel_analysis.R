#==============================================================================
# 03_panel_analysis.R
# Repeated-measures models for BOP, PD4, PD5, CAL2 and PLI.
#
# Model: generalized least squares (nlme::gls) with an unstructured covariance
# across Baseline, Visit 3 and Visit 5 (corSymm + varIdent; the R equivalent of
# SAS REPEATED / TYPE=UN), fitted by REML. Type III (marginal) F tests on the
# residual denominator df; gls does not implement Kenward-Roger.
#
# Runs:  (1) primary, all available visits       -> Table 2, Table 4, Figures 1-3
#        (2) complete cases (both follow-ups)     -> Table S1, Table 4
#        (3) pregnancy-wide covariates            -> Table S3
#
# Outputs: output/Panel_AllEstimates.csv, output/Panel_LSMeans.csv
#==============================================================================

opt <- readRDS(out_path("opt_analysis.rds"))

# ---- Wide to long ------------------------------------------------------------
keep <- c("PID", "Treat", "Group", "Race3", "Edu3", "Age", "BMI",
          "HTN_base", "DM_base", "INF_base", "TOB_base",
          "HTN_preg", "DM_preg", "INF_preg", "TOB_preg", "complete")
mk <- function(v, pre) {
  out <- opt[, c(keep, paste0(pre, OUTCOMES))]
  names(out)[(length(keep) + 1):ncol(out)] <- OUTCOMES
  out$Visit <- v
  out
}
long <- bind_rows(mk(1L, "BL"), mk(2L, "V3"), mk(3L, "V5")) %>% arrange(PID, Visit)
long$Visit <- factor(long$Visit, levels = 1:3, labels = c("Baseline", "Visit3", "Visit5"))
long$PID   <- factor(long$PID)

# ---- Model ---------------------------------------------------------------------
fit_panel_model <- function(outcome_var, data, covset = c("base", "preg")) {
  covset <- match.arg(covset)
  sfx <- if (covset == "base") "_base" else "_preg"
  cov_terms <- paste0(c("HTN", "DM", "INF", "TOB"), sfx, collapse = " + ")
  f <- as.formula(paste0(
    "Y ~ Group * Visit + Age + I(Age^2) + BMI + I(BMI^2) + Race3 + Edu3 + ", cov_terms))
  data$Y   <- data[[outcome_var]]
  used     <- c("Y", "Age", "BMI", "Race3", "Edu3", paste0(c("HTN", "DM", "INF", "TOB"), sfx))
  data_fit <- data[stats::complete.cases(data[, used]), , drop = FALSE]
  model <- tryCatch(
    gls(f, data = data_fit, method = "REML",
        correlation = corSymm(form = ~ as.numeric(Visit) | PID),
        weights     = varIdent(form = ~ 1 | Visit),
        control     = glsControl(maxIter = 200, msMaxIter = 200)),
    error = function(e) { cat("  ERROR", outcome_var, ":", conditionMessage(e), "\n"); NULL })
  list(model = model, data = data_fit)
}

# Contrast weights on the emmeans grid (Group within Visit):
#   1 C-Baseline, 2 T-Baseline, 3 C-Visit3, 4 T-Visit3, 5 C-Visit5, 6 T-Visit5
# All contrasts are control minus treatment (positive = less disease with treatment).
contr_list <- list(
  BL_level_C_minus_T      = c( 1, -1, 0,  0, 0,  0),
  V3_level_C_minus_T      = c( 0,  0, 1, -1, 0,  0),
  V5_level_C_minus_T      = c( 0,  0, 0,  0, 1, -1),
  Change_BLtoV5_C_minus_T = c(-1,  1, 0,  0, 1, -1)   # same scale as g-estimation psi
)

lsm_store <- list()
run_all <- function(data, covset, label) {
  cat("\n==============================================================\n")
  cat("MODELS:", label, "\n")
  cat("==============================================================\n")
  res <- list()
  for (o in OUTCOMES) {
    fit <- fit_panel_model(o, data, covset)
    if (is.null(fit$model)) next
    m <- fit$model; fd <- fit$data

    aov_tab <- anova(m, type = "marginal")
    int_row <- aov_tab[grepl("Group:Visit", rownames(aov_tab)), ]
    den_df  <- suppressWarnings(as.numeric(gsub("\\D", "", attr(aov_tab, "label"))))

    emm  <- emmeans(m, ~ Group * Visit, data = fd, mode = "df.error")
    grid <- as.data.frame(summary(emm))
    stopifnot(identical(as.character(grid$Group), c("C", "T", "C", "T", "C", "T")))
    lsm_store[[paste(label, o)]] <<- data.frame(Outcome = o, Covariates = label,
      grid[, c("Group", "Visit", "emmean", "lower.CL", "upper.CL")], row.names = NULL)

    cs <- summary(contrast(emm, contr_list, adjust = MULTIPLICITY), infer = c(TRUE, TRUE))
    get <- function(name, col) cs[[col]][cs$contrast == name]

    cat(sprintf("%-5s n=%d obs=%d  F(%d,%s)=%.2f  V5 diff=%.2f [%.2f, %.2f]\n",
                o, length(unique(fd$PID)), nobs(m), int_row[["numDF"]], den_df,
                int_row[["F-value"]], get("V5_level_C_minus_T", "estimate"),
                get("V5_level_C_minus_T", "lower.CL"), get("V5_level_C_minus_T", "upper.CL")))

    res[[o]] <- data.frame(
      Outcome = o, Covariates = label,
      N_part = length(unique(fd$PID)), N_obs = nobs(m),
      F_GxV = int_row[["F-value"]], df_num = int_row[["numDF"]], df_den = den_df,
      p_GxV = int_row[["p-value"]],
      LSM_C_BL = grid$emmean[1], LSM_T_BL = grid$emmean[2],
      LSM_C_V3 = grid$emmean[3], LSM_T_V3 = grid$emmean[4],
      LSM_C_V5 = grid$emmean[5], LSM_T_V5 = grid$emmean[6],
      BL_diff = get("BL_level_C_minus_T", "estimate"), BL_p = get("BL_level_C_minus_T", "p.value"),
      V3_diff = get("V3_level_C_minus_T", "estimate"), V3_p = get("V3_level_C_minus_T", "p.value"),
      V5_diff = get("V5_level_C_minus_T", "estimate"),
      V5_lcl  = get("V5_level_C_minus_T", "lower.CL"), V5_ucl = get("V5_level_C_minus_T", "upper.CL"),
      V5_p    = get("V5_level_C_minus_T", "p.value"),
      Change_diff = get("Change_BLtoV5_C_minus_T", "estimate"),
      Change_lcl  = get("Change_BLtoV5_C_minus_T", "lower.CL"),
      Change_ucl  = get("Change_BLtoV5_C_minus_T", "upper.CL"),
      row.names = NULL)
  }
  do.call(rbind, res)
}

primary       <- run_all(long, "base", "Primary")
complete_case <- run_all(long %>% filter(complete == 1), "base", "Complete case")
sens          <- if (SENSITIVITY_PREG_COVARS) run_all(long, "preg", "Pregnancy-wide covariates") else NULL

panel_all <- bind_rows(primary, complete_case, sens)
write.csv(panel_all, out_path("Panel_AllEstimates.csv"), row.names = FALSE)
write.csv(do.call(rbind, lsm_store), out_path("Panel_LSMeans.csv"), row.names = FALSE)
cat("\n")
