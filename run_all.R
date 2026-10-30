#==============================================================================
# run_all.R
# Reproduces every table and figure in:
#   Suthar I, et al. Evaluating Periodontal Treatment During Pregnancy with
#   Repeated-Measures Mixed Models and G-Estimation to Address Incomplete
#   Follow-up. (Submitted to the Journal of Clinical Periodontology.)
#
# Secondary analysis of the Obstetrics and Periodontal Therapy (OPT) trial,
# ClinicalTrials.gov NCT00066131. Data: https://doi.org/10.13020/D6KW2D
#
# Usage: set the working directory to this folder, place the data file in
# data/ (see README), then run  source("run_all.R")
#==============================================================================

source("R/00_config.R")
source("R/01_data_preparation.R")   # analysis variables      -> output/opt_analysis.rds
source("R/02_descriptives.R")       # Table 1, Table S2, follow-up counts
source("R/03_panel_analysis.R")     # Tables 2 and S1, panel columns of Tables 4 and S3
source("R/04_gestimation.R")        # Table 3, g-estimation columns of Tables 4 and S3
source("R/05_figures.R")            # Figures 1-3, S1, S2

writeLines(capture.output(sessionInfo()), out_path("sessionInfo.txt"))
cat("All outputs written to", normalizePath(OUT_DIR), "\n")
