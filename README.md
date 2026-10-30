# Periodontal treatment during pregnancy: repeated-measures and g-estimation analysis of the OPT trial

Analysis code for:

> Suthar I, Chen B, Hasturk H, Rokonuzzaman SM, Toma NT, Suen S-C, Merchant AT. *Evaluating Periodontal Treatment During Pregnancy with Repeated-Measures Mixed Models and G-Estimation to Address Incomplete Follow-up.* Submitted to the Journal of Clinical Periodontology.

This is a secondary analysis of the Obstetrics and Periodontal Therapy (OPT) trial (ClinicalTrials.gov NCT00066131). Running `run_all.R` reproduces every table and figure in the manuscript from the public-use data.

## Data

The data are **not included** in this repository. They are the public-use OPT dataset:

> Hodges, J. S., & Michalowicz, B. S. (2009). *Public-use data from the Obstetrics and Periodontal Therapy (OPT) Study.* Data Repository for the University of Minnesota. https://doi.org/10.13020/D6KW2D

Download it and save the main merged participant file as `data/opt_main_merged_new.xlsx`.

## How to run

```r
setwd("path/to/this/repository")
source("run_all.R")
```

All results are written to `output/`, and figures to `output/figures/` as PNG (300 dpi) and PDF. The full run, with 1,000 bootstrap replicates for each outcome and specification, takes several minutes. For a quick check, set `TEST_MODE <- TRUE` in `R/00_config.R`, which uses 25 replicates.

Requirements: R (4.1 or later) and the packages `readxl`, `dplyr`, `nlme`, `emmeans`, `boot` and `ggplot2`. Any missing packages are installed on the first run. The exact versions used are saved to `output/sessionInfo.txt`.

## Scripts

| Script | What it does |
|---|---|
| `run_all.R` | Runs the scripts below in order |
| `R/00_config.R` | Packages, file paths, run options, random seed |
| `R/01_data_preparation.R` | Imports the data and builds every analysis variable once (covariates, race, change scores, follow-up and censoring indicators) |
| `R/02_descriptives.R` | Baseline characteristics, unadjusted means by visit, follow-up counts |
| `R/03_panel_analysis.R` | Repeated-measures models: primary, complete-case, and pregnancy-wide covariate sensitivity |
| `R/04_gestimation.R` | G-estimation with stabilized IPCW, weight diagnostics, bootstrap BC-a intervals |
| `R/05_figures.R` | All figures, drawn from the output files of scripts 02-04 |

## Where each manuscript result comes from

| Manuscript | Output file |
|---|---|
| Table 1 | `Table1_Baseline.csv` |
| Table 2 | `Panel_AllEstimates.csv` (Primary rows), `Panel_LSMeans.csv` |
| Table 3 | `Gest_Bootstrap.csv` |
| Table 4 | `Panel_AllEstimates.csv` (Primary and Complete case rows), `Gest_Bootstrap.csv` |
| Supplementary Table S1 | `Panel_AllEstimates.csv` (Complete case rows) |
| Supplementary Table S2 | `TableS2_Unadjusted.csv` |
| Supplementary Table S3 | `Panel_AllEstimates.csv` (Pregnancy-wide covariates rows), `Gest_PointEstimates.csv` (spec `preg`) |
| Stabilized weights (Methods) | `Gest_WeightDiagnostics.csv` |
| Follow-up counts (Results) | `FollowUp_Counts.csv` |
| Figures 1-3, S1, S2 | `figures/` |

## Methods in brief

- **Covariates.** All models adjust for covariates recorded at enrollment: age, BMI (each with a squared term), race, education, chronic hypertension, diabetes known at enrollment, bacterial vaginosis at Baseline and tobacco use at Baseline. Race is built from the trial's binary race indicators, Black first, then White; all others are Other.
- **Repeated-measures models.** Generalized least squares (`nlme::gls`) with an unstructured covariance across Baseline, Visit 3 and Visit 5, fitted by REML. Type III F tests use the residual degrees of freedom.
- **G-estimation.** A structural nested mean model for change from Baseline to Visit 5. Participants who missed Visit 5 are censored and reweighted with stabilized inverse probability of censoring weights. Uncertainty comes from participant-level bootstrap resampling with BC-a intervals.

## License

MIT. See `LICENSE`.

## Contact

Ishan Suthar, Department of Epidemiology & Biostatistics, Arnold School of Public Health, University of South Carolina. isuthar@email.sc.edu
