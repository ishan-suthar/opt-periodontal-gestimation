#==============================================================================
# 01_data_preparation.R
# Imports the public-use OPT data and builds every analysis variable once,
# so the descriptive tables, the repeated-measures models and the
# g-estimation all use identical definitions.
#
# Output: output/opt_analysis.rds  (one row per randomized participant, n = 823)
#==============================================================================

yn <- function(x) tolower(trimws(as.character(x))) %in% c("y", "yes", "1")

opt_raw <- read_excel(DATA_FILE)
cat("Imported", nrow(opt_raw), "rows,", ncol(opt_raw), "columns.\n\n")

opt <- opt_raw %>%
  mutate(
    Group = factor(ifelse(yn(Group) | Group == "T", "T", "C"), levels = c("C", "T")),
    Treat = as.integer(Group == "T"),

    # ---- Covariates recorded at enrollment (primary analyses) ----------------
    HTN_base = as.integer(yn(Hypertension)),                     # chronic hypertension
    DM_base  = as.integer(yn(Diabetes)),                         # diabetes known at enrollment
    INF_base = as.integer(!is.na(BLBacvag) & BLBacvag == 1),     # bacterial vaginosis at Baseline
    TOB_base = as.integer(!is.na(BLCigDay) & BLCigDay > 0),      # tobacco use at Baseline

    # ---- Pregnancy-wide composites (sensitivity analysis only) ---------------
    # These include conditions recorded after randomization (preeclampsia,
    # gestational diabetes, infections and smoking during pregnancy).
    HTN_preg = as.integer(yn(Hypertension) | yn(Preeclamp)),
    DM_preg  = as.integer(yn(Diabetes)     | yn(Gestdiab)),
    INF_preg = as.integer(yn(Gonorrhea) | yn(Chlamydia) | yn(StrepB) | yn(UTI)),
    TOB_preg = as.integer(yn(UseTob)),

    # ---- Race: built from the binary race indicators, Black first ------------
    # Participants reporting more than one race are assigned hierarchically
    # (Black, then White); all others, including those reporting no race
    # indicator, are Other. Reference level = White.
    Race3 = factor(ifelse(yn(Black), "Black", ifelse(yn(White), "White", "Other")),
                   levels = c("White", "Black", "Other")),

    # Education: 1 = <8 years, 2 = 8-12 years (reference), 3 = >12 years
    Edu3 = relevel(factor(Edu), ref = "2"),

    # ---- Change from Baseline to Visit 5 (g-estimation outcomes) -------------
    DiffBOP  = V5BOP  - BLBOP,
    DiffPD4  = V5PD4  - BLPD4,
    DiffPD5  = V5PD5  - BLPD5,
    DiffCAL2 = V5CAL2 - BLCAL2,
    DiffPlI  = V5PlI  - BLPlI,

    # ---- Follow-up (all five outcomes are co-missing within an examination) --
    obsV3    = as.integer(!is.na(V3BOP)),
    obsV5    = as.integer(!is.na(V5BOP)),
    complete = as.integer(obsV3 == 1 & obsV5 == 1),
    cens     = as.integer(obsV5 == 0)            # censored = Visit 5 missed
  )

# ---- Audit printed for the Methods section ------------------------------------
race_ind <- cbind(yn(opt_raw$Black), yn(opt_raw$White), yn(opt_raw$NatAm), yn(opt_raw$Asian))
cat("Race: participants with more than one indicator:", sum(rowSums(race_ind) > 1), "\n")
cat("Race: participants with no indicator (Other):   ", sum(rowSums(race_ind) == 0), "\n")
cat("Race by arm:\n"); print(table(opt$Race3, opt$Group))

model_covs <- c("Age", "BMI", "Race3", "Edu3", "HTN_base", "DM_base", "INF_base", "TOB_base")
cat("\nMissing values in model covariates:\n")
print(sapply(opt[, model_covs], function(x) sum(is.na(x))))
cat("Participants with complete covariates (analytic sample):",
    sum(complete.cases(opt[, model_covs])), "\n")

cat("\nFollow-up patterns:\n")
cat("  Attended both follow-up examinations:", sum(opt$complete), "\n")
cat("  Missed Visit 3, returned at Visit 5 :", sum(opt$obsV3 == 0 & opt$obsV5 == 1), "\n")
cat("  Attended Visit 3, missed Visit 5    :", sum(opt$obsV3 == 1 & opt$obsV5 == 0), "\n")
cat("  Missed both follow-up examinations  :", sum(opt$obsV3 == 0 & opt$obsV5 == 0), "\n")
cat("  Censored (Visit 5 missed)           :", sum(opt$cens), "\n\n")

saveRDS(opt, out_path("opt_analysis.rds"))
