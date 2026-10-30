#==============================================================================
# 00_config.R
# Run options, packages and file paths shared by every script.
#==============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("readxl", "dplyr", "nlme", "emmeans", "boot", "ggplot2")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(nlme)
  library(emmeans); library(boot); library(ggplot2)
})

# ---- Paths -------------------------------------------------------------------
# Download the public-use OPT data (https://doi.org/10.13020/D6KW2D) and save
# the main merged participant file as data/opt_main_merged_new.xlsx.
DATA_FILE  <- file.path("data", "opt_main_merged_new.xlsx")
OUT_DIR    <- "output"
FIG_DIR    <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- function(name) file.path(OUT_DIR, name)

# ---- Analysis options ----------------------------------------------------------
SENSITIVITY_PREG_COVARS <- TRUE    # pregnancy-wide covariate sensitivity (Table S3)
MULTIPLICITY            <- "none"  # "none" or "sidak" across visit contrasts
TEST_MODE               <- FALSE   # TRUE = 25 bootstrap replicates, for a quick check
R_BOOT                  <- if (TEST_MODE) 25 else 1000
SEED                    <- 20241

OUTCOMES       <- c("BOP", "PD4", "PD5", "CAL2", "PlI")
OUTCOME_LABELS <- c(BOP = "BOP (%)", PD4 = "PD4 (%)", PD5 = "PD5 (%)",
                    CAL2 = "CAL2 (%)", PlI = "PLI (Silness-L\u00f6e units)")

if (!file.exists(DATA_FILE)) {
  stop("Data file not found: ", DATA_FILE,
       "\nDownload the public-use OPT data from https://doi.org/10.13020/D6KW2D",
       " and save the main merged file as ", DATA_FILE)
}
