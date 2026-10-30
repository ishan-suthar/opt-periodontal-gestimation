#==============================================================================
# 04_gestimation.R
# G-estimation of a structural nested mean model for change from Baseline to
# Visit 5, with stabilized inverse probability of censoring weights (IPCW).
#
#   Model:      E[Y^a - Y^0 | A = a, L] = psi * a
#   Censoring:  C = 1 if the Visit 5 examination was missed
#   Weights:    sw = P(C = 0 | A) / P(C = 0 | A, L); censored participants get 0
#   Estimator:  sum sw (A - p(L)) (Y - psi A) = 0, so
#               psi = sum sw (A - p(L)) Y / sum sw (A - p(L)) A
#               with p(L) from an sw-weighted logistic regression of A on L
#               among uncensored participants.
#
# psi is an intention-to-treat effect of randomized assignment. Negative values
# mean greater improvement in the early-treatment arm.
#
# Specifications:
#   base        enrollment covariates + Baseline outcome value     (PRIMARY, Table 3)
#   unweighted  treatment-only censoring models, so every sw = 1   (Table 3)
#   preg        pregnancy-wide covariates                          (Table S3)
#   base_nobl   enrollment covariates without the Baseline outcome (not reported)
#
# Uncertainty: 1,000 participant-level bootstrap resamples, BC-a intervals.
# Outputs: output/Gest_Bootstrap.csv, Gest_PointEstimates.csv,
#          Gest_WeightDiagnostics.csv, Gest_BootReplicates.csv
#==============================================================================

opt <- readRDS(out_path("opt_analysis.rds"))

model_covs <- c("Age", "BMI", "Race3", "Edu3", "HTN_base", "DM_base", "INF_base", "TOB_base")
optg <- opt[stats::complete.cases(opt[, model_covs]), ]
cat("G-estimation analytic sample:", nrow(optg), "| censored at Visit 5:", sum(optg$cens), "\n\n")

outcomes_g <- lapply(setNames(OUTCOMES, OUTCOMES),
                     function(o) list(diff = paste0("Diff", o), bl = paste0("BL", o)))

cov_string <- function(spec, blvar) {
  switch(spec,
    base      = paste0("Age + I(Age^2) + BMI + I(BMI^2) + Race3 + Edu3 + ",
                       "HTN_base + DM_base + INF_base + TOB_base + ", blvar),
    base_nobl = "Age + I(Age^2) + BMI + I(BMI^2) + Race3 + Edu3 + HTN_base + DM_base + INF_base + TOB_base",
    preg      = "Age + I(Age^2) + BMI + I(BMI^2) + Race3 + Edu3 + HTN_preg + DM_preg + INF_preg + TOB_preg",
    stop("unknown spec"))
}

est_psi <- function(data, diffvar, blvar, spec = "base", return_weights = FALSE) {
  d <- data
  d$Y   <- d[[diffvar]]
  d$unc <- 1 - d$cens

  # Stabilized censoring weights
  denom_formula <- if (spec == "unweighted") unc ~ Treat else
    as.formula(paste("unc ~ Treat +", cov_string(spec, blvar)))
  denom_fit <- tryCatch(suppressWarnings(glm(denom_formula, data = d, family = binomial())),
                        error = function(e) NULL)
  numer_fit <- tryCatch(suppressWarnings(glm(unc ~ Treat, data = d, family = binomial())),
                        error = function(e) NULL)
  if (is.null(denom_fit) || is.null(numer_fit)) return(if (return_weights) NULL else NA_real_)
  d$pd <- predict(denom_fit, newdata = d, type = "response")
  d$pn <- predict(numer_fit, newdata = d, type = "response")

  d_unc <- d[d$cens == 0, , drop = FALSE]
  d_unc$sw <- d_unc$pn / d_unc$pd
  d_unc <- d_unc[is.finite(d_unc$sw) & d_unc$sw > 0 & !is.na(d_unc$Y), , drop = FALSE]
  if (nrow(d_unc) == 0) return(if (return_weights) NULL else NA_real_)

  # sw-weighted propensity model p(L) = P(A = 1 | L) among the uncensored
  pa_formula <- as.formula(paste("Treat ~",
                  if (spec == "unweighted") cov_string("base", blvar) else cov_string(spec, blvar)))
  pa_fit <- tryCatch(suppressWarnings(glm(pa_formula, data = d_unc,
                                          family = binomial(), weights = sw)),
                     error = function(e) NULL)
  if (is.null(pa_fit)) return(if (return_weights) NULL else NA_real_)
  d_unc$pA <- predict(pa_fit, newdata = d_unc, type = "response")
  if (return_weights) return(d_unc)

  # Closed-form solution of the estimating equation
  A <- d_unc$Treat; Y <- d_unc$Y; sw <- d_unc$sw; dlt <- A - d_unc$pA
  den <- sum(sw * A * dlt)
  if (!is.finite(den) || abs(den) < 1e-6) return(NA_real_)
  psi <- sum(sw * Y * dlt) / den
  if (!is.finite(psi) || abs(psi) > 200) return(NA_real_)   # guard against numerical blow-up
  psi
}

# ---- Point estimates for every specification ----------------------------------
specs <- c("base", "unweighted", "preg", "base_nobl")
point <- expand.grid(Outcome = OUTCOMES, Spec = specs, stringsAsFactors = FALSE)
point$psi <- mapply(function(o, s) est_psi(optg, outcomes_g[[o]]$diff, outcomes_g[[o]]$bl, s),
                    point$Outcome, point$Spec)
cat("Point estimates (psi):\n")
print(reshape(point, idvar = "Outcome", timevar = "Spec", direction = "wide"),
      row.names = FALSE, digits = 4)

# ---- Weight diagnostics (reported in Methods) ---------------------------------
wt_diag <- do.call(rbind, lapply(OUTCOMES, function(o) {
  w <- est_psi(optg, outcomes_g[[o]]$diff, outcomes_g[[o]]$bl, "base", return_weights = TRUE)
  data.frame(Outcome = o, n_uncensored = nrow(w), mean = mean(w$sw), sd = sd(w$sw),
             min = min(w$sw), max = max(w$sw), row.names = NULL)
}))
cat("\nStabilized weights (primary specification):\n"); print(wt_diag, row.names = FALSE, digits = 3)
w_un <- est_psi(optg, outcomes_g$BOP$diff, outcomes_g$BOP$bl, "unweighted", return_weights = TRUE)
stopifnot(isTRUE(all.equal(unname(w_un$sw), rep(1, nrow(w_un)))))   # unweighted really is sw = 1

# ---- Bootstrap ------------------------------------------------------------------
cat("\nBootstrap:", R_BOOT, "resamples per outcome and specification\n")
set.seed(SEED)
res <- list(); reps <- list()
for (s in c("base", "unweighted")) {
  for (o in OUTCOMES) {
    stat <- local({ s_ <- s; o_ <- o
      function(data, i) est_psi(data[i, , drop = FALSE], outcomes_g[[o_]]$diff,
                                outcomes_g[[o_]]$bl, s_) })
    bo  <- boot(optg, statistic = stat, R = R_BOOT)
    ci  <- tryCatch(boot.ci(bo, type = "bca"),  error = function(e) NULL)
    pct <- tryCatch(boot.ci(bo, type = "perc"), error = function(e) NULL)
    reps[[paste(s, o)]] <- data.frame(Spec = s, Outcome = o, replicate = seq_len(R_BOOT),
                                      psi = as.numeric(bo$t))
    res[[paste(s, o)]] <- data.frame(
      Spec = s, Outcome = o, psi = bo$t0, boot_se = sd(bo$t, na.rm = TRUE),
      bca_lcl  = if (!is.null(ci))  ci$bca[4]      else NA,
      bca_ucl  = if (!is.null(ci))  ci$bca[5]      else NA,
      perc_lcl = if (!is.null(pct)) pct$percent[4] else NA,
      perc_ucl = if (!is.null(pct)) pct$percent[5] else NA,
      failed_replicates = sum(!is.finite(bo$t)), row.names = NULL)
    cat(sprintf("  %-10s %-5s psi=%8.3f  BC-a [%8.3f, %8.3f]  failed=%d\n", s, o, bo$t0,
                res[[paste(s, o)]]$bca_lcl, res[[paste(s, o)]]$bca_ucl,
                res[[paste(s, o)]]$failed_replicates))
  }
}

write.csv(do.call(rbind, res),  out_path("Gest_Bootstrap.csv"),         row.names = FALSE)
write.csv(point,                out_path("Gest_PointEstimates.csv"),     row.names = FALSE)
write.csv(wt_diag,              out_path("Gest_WeightDiagnostics.csv"),  row.names = FALSE)
write.csv(do.call(rbind, reps), out_path("Gest_BootReplicates.csv"),     row.names = FALSE)
cat("\n")
