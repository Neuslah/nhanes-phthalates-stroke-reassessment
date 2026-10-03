options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 01_static_package_tests.R <package_root>")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)

required_files <- c(
  "README.md", "LICENSE", "config/nhanes-components.csv",
  "config/translated-cache-manifest.csv", "config/variable-dictionary.csv", "environment/package-versions.csv",
  "environment/session-info.txt", "R/00_run_deterministic_pipeline.R",
  "tests/03_inherited_mixture_output_tests.R",
  sprintf("R/%s", c(
    "01_download_or_import_nhanes.R", "02_build_analytic_dataset.R", "03_primary_models.R",
    "04_primary_fdr.R", "05_sex_interaction.R", "06_common_cycle_sensitivity.R",
    "07_rcs_explicit_3k.R", "08a_grouped_derive_current_weights.R",
    "08b_grouped_models_current_weights.R", "09_era_analyses.R",
    "10_male_subgroup_global_interactions.R", "11_wqs_full_sample.R",
    "12_qgcomp.R", "13_wqs_repeated_holdout.R"
  ))
)
missing_files <- required_files[!file.exists(file.path(root, required_files))]
if (length(missing_files)) stop("Missing package files: ", paste(missing_files, collapse = ", "))

components <- read.csv(file.path(root, "config", "nhanes-components.csv"), check.names = FALSE)
if (nrow(components) != 66L || anyDuplicated(components$Component) || any(!nzchar(components$Source_SHA))) {
  stop("Component manifest assertion failed")
}
translated <- read.csv(file.path(root, "config", "translated-cache-manifest.csv"), check.names = FALSE)
if (nrow(translated) != 66L || anyDuplicated(translated$Component) || any(!nzchar(translated$Translated_RDS_SHA))) {
  stop("Translated-cache manifest assertion failed")
}
dictionary <- read.csv(file.path(root, "config", "variable-dictionary.csv"), check.names = FALSE)
required_variables <- c(
  "SEQN", "cycle", "stroke", "sex", "RIDAGEYR", "race_eth_f", "edu_f", "marital_f",
  "pir", "bmi", "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "WTMEC_comb", "URXUCR",
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln", "MCPP_ln",
  "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln", "ln_sigma_LMW", "ln_sigma_HMW",
  "ln_sigma_DEHP"
)
missing_variables <- setdiff(required_variables, dictionary$Variable)
if (length(missing_variables)) stop("Variable dictionary missing: ", paste(missing_variables, collapse = ", "))

r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
r_text <- paste(unlist(lapply(r_files, readLines, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
forbidden_patterns <- c(
  "[A-Za-z]:[/\\\\]", "setwd\\s*\\(", "MainRegression_All", "OR_All", "P_All",
  paste0("P", "A-[0-9]"), paste0("S", "C-[0-9]"), paste0("A", "S-[0-9]"),
  paste0("PROJECT", "_STATUS"), paste0("改稿", "修订蓝图")
)
hits <- vapply(forbidden_patterns, function(pattern) grepl(pattern, r_text, perl = TRUE), logical(1))
if (any(hits)) stop("Forbidden code pattern(s): ", paste(forbidden_patterns[hits], collapse = ", "))

all_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
participant_data <- all_files[grepl("\\.(rds|rdata|xpt|sas7bdat)$", all_files, ignore.case = TRUE)]
if (length(participant_data)) stop("Participant-level data present in package")
cat("STATIC_PACKAGE_TESTS=PASS\n")
