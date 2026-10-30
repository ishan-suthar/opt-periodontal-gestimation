#==============================================================================
# 02_descriptives.R
# Table 1 (baseline characteristics), Supplementary Table S2 (unadjusted
# means by visit) and the follow-up counts used in the Results and in
# Supplementary Figure S1.
#
# Tests: Welch t-tests for continuous variables; Pearson chi-square tests
# without continuity correction for categorical variables.
#==============================================================================

opt <- readRDS(out_path("opt_analysis.rds"))
g   <- opt$Group

fmt_mean <- function(x, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f)"),
                                        mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_n    <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
p_t      <- function(v) t.test(v ~ g)$p.value
p_chi    <- function(v) chisq.test(table(v, g), correct = FALSE)$p.value

rows <- list()
add <- function(label, c_val, t_val, p = NA) {
  rows[[length(rows) + 1]] <<- data.frame(Characteristic = label, Control = c_val,
                                          Treatment = t_val, p = p)
}
add("n (randomized)", sum(g == "C"), sum(g == "T"))
add("Age, years, mean (SD)", fmt_mean(opt$Age[g == "C"]), fmt_mean(opt$Age[g == "T"]), p_t(opt$Age))
add("BMI, kg/m2, mean (SD)", fmt_mean(opt$BMI[g == "C"]), fmt_mean(opt$BMI[g == "T"]), p_t(opt$BMI))
add("Race/ethnicity", "", "", p_chi(opt$Race3))
for (lv in c("Black", "White", "Other"))
  add(paste0("  ", lv), fmt_n(opt$Race3[g == "C"] == lv), fmt_n(opt$Race3[g == "T"] == lv))
add("Education", "", "", p_chi(opt$Edu3))
for (lv in c("1", "2", "3"))
  add(paste0("  ", c("1" = "<8 years", "2" = "8-12 years", "3" = ">12 years")[lv]),
      fmt_n(opt$Edu3[g == "C"] == lv), fmt_n(opt$Edu3[g == "T"] == lv))
comorb <- c(HTN_base = "Chronic hypertension", DM_base = "Diabetes at enrollment",
            INF_base = "Bacterial vaginosis at Baseline", TOB_base = "Tobacco use at Baseline")
for (v in names(comorb))
  add(comorb[[v]], fmt_n(opt[[v]][g == "C"] == 1), fmt_n(opt[[v]][g == "T"] == 1), p_chi(opt[[v]]))
add("Attended both follow-up visits", fmt_n(opt$complete[g == "C"] == 1),
    fmt_n(opt$complete[g == "T"] == 1), p_chi(opt$complete))
add("Missed at least one follow-up visit", fmt_n(opt$complete[g == "C"] == 0),
    fmt_n(opt$complete[g == "T"] == 0))
for (o in OUTCOMES) {
  v <- opt[[paste0("BL", o)]]; d <- if (o == "PlI") 2 else 1
  add(paste("Baseline", OUTCOME_LABELS[o], "mean (SD)"),
      fmt_mean(v[g == "C"], d), fmt_mean(v[g == "T"], d), p_t(v))
}
table1 <- do.call(rbind, rows)
table1$p <- ifelse(is.na(table1$p), "", sprintf("%.2f", table1$p))
write.csv(table1, out_path("Table1_Baseline.csv"), row.names = FALSE)
cat("Table 1:\n"); print(table1, row.names = FALSE)

# ---- Supplementary Table S2: unadjusted means by visit -------------------------
visits <- c(Baseline = "BL", "Visit 3" = "V3", "Visit 5" = "V5")
s2 <- do.call(rbind, lapply(OUTCOMES, function(o) do.call(rbind, lapply(names(visits), function(vn) {
  y <- opt[[paste0(visits[[vn]], o)]]
  data.frame(Outcome = o, Visit = vn,
             Control_N = sum(!is.na(y[g == "C"])),   Control_Mean_SD   = fmt_mean(y[g == "C"], 2),
             Treatment_N = sum(!is.na(y[g == "T"])), Treatment_Mean_SD = fmt_mean(y[g == "T"], 2),
             Raw_C_minus_T = round(mean(y[g == "C"], na.rm = TRUE) - mean(y[g == "T"], na.rm = TRUE), 2))
}))))
write.csv(s2, out_path("TableS2_Unadjusted.csv"), row.names = FALSE)

# ---- Follow-up counts (Results; Supplementary Figure S1) -----------------------
fu <- opt %>% group_by(Group) %>% summarise(
  n = n(), both = sum(complete), missV3_retV5 = sum(obsV3 == 0 & obsV5 == 1),
  V3_missV5 = sum(obsV3 == 1 & obsV5 == 0), neither = sum(obsV3 == 0 & obsV5 == 0),
  incomplete = sum(complete == 0), examined_V5 = sum(obsV5), .groups = "drop")
chi_fu <- chisq.test(table(opt$Group, opt$complete == 0), correct = FALSE)
fu$chisq <- unname(chi_fu$statistic); fu$chisq_p <- chi_fu$p.value
write.csv(fu, out_path("FollowUp_Counts.csv"), row.names = FALSE)
cat("\nFollow-up by arm:\n"); print(as.data.frame(fu))

# Baseline severity by follow-up status (Results)
sev <- opt %>% group_by(Incomplete = complete == 0) %>%
  summarise(across(c(BLBOP, BLCAL2), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
cat("\nBaseline BOP and CAL2 by follow-up status:\n"); print(as.data.frame(sev))
cat("\n")
