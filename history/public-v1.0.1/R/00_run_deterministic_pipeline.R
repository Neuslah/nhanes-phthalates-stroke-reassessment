# Run all deterministic analyses and lightweight mixture validation in an explicit workspace.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript 00_run_deterministic_pipeline.R",
    "<package_root> <official_source_cache> <translated_cache_or_dash> <workspace_root> <output_root>"
  ))
}
package_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
source_cache <- gsub("\\\\", "/", args[[2]])
translated_cache <- gsub("\\\\", "/", args[[3]])
workspace_root <- gsub("\\\\", "/", args[[4]])
output_root <- gsub("\\\\", "/", args[[5]])
dir.create(workspace_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) stop("Rscript executable not found")

run_step <- function(label, script, step_args) {
  log_path <- file.path(output_root, paste0(label, ".log"))
  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(file.path(package_root, "R", script)), vapply(step_args, shQuote, character(1))),
    stdout = log_path, stderr = log_path
  )
  if (!identical(status, 0L)) stop("Pipeline step failed: ", label, "; see ", log_path)
  cat(label, "=PASS\n", sep = "")
}

dataset_out <- file.path(output_root, "dataset")
primary_out <- file.path(output_root, "primary")
secondary_out <- file.path(output_root, "secondary")
rcs_out <- file.path(output_root, "rcs")
grouped_out <- file.path(output_root, "grouped")
era_out <- file.path(output_root, "era")
subgroup_out <- file.path(output_root, "subgroup")
mixture_out <- file.path(output_root, "mixture")
invisible(vapply(c(dataset_out, primary_out, secondary_out, rcs_out, grouped_out, era_out,
                   subgroup_out, mixture_out), dir.create, logical(1), recursive = TRUE, showWarnings = FALSE))

run_step("01_source_import", "01_download_or_import_nhanes.R", c(
  package_root, source_cache, translated_cache, workspace_root
))
run_step("02_dataset_build", "02_build_analytic_dataset.R", c(workspace_root, dataset_out))
analytic_rds <- file.path(workspace_root, "data_derived", "NPS_DATA_MAIN_analytic_V5.rds")
run_step("03_primary", "03_primary_models.R", c(analytic_rds, primary_out))
run_step("04_primary_fdr", "04_primary_fdr.R", c(
  file.path(primary_out, "NPS_OUT_TBL_MainRegression_Female_V5.csv"),
  file.path(primary_out, "NPS_OUT_TBL_MainRegression_Male_V5.csv"),
  file.path(primary_out, "primary_fdr_20test.csv")
))
run_step("05_sex_interaction", "05_sex_interaction.R", c(analytic_rds, file.path(secondary_out, "sex_interaction.csv")))
run_step("06_common_cycle", "06_common_cycle_sensitivity.R", c(
  analytic_rds, file.path(secondary_out, "common_cycle_models.csv"),
  file.path(secondary_out, "common_cycle_fdr.csv")
))
run_step("07_rcs", "07_rcs_explicit_3k.R", c(
  analytic_rds, file.path(rcs_out, "rcs_combined.csv"), file.path(rcs_out, "rcs_female.csv"),
  file.path(rcs_out, "rcs_male.csv"), file.path(rcs_out, "rcs_knots_qc.csv"),
  file.path(rcs_out, "rcs_warnings.csv")
))
grouped_rds <- file.path(workspace_root, "data_derived", "grouped_current_weights.rds")
run_step("08a_grouped_derive", "08a_grouped_derive_current_weights.R", c(
  analytic_rds, grouped_rds, file.path(grouped_out, "grouped_weight_qc.csv")
))
run_step("08b_grouped_models", "08b_grouped_models_current_weights.R", c(
  analytic_rds, grouped_rds, file.path(grouped_out, "grouped_results.csv"),
  file.path(grouped_out, "grouped_model_qc.csv")
))
run_step("09_era", "09_era_analyses.R", c(
  analytic_rds, file.path(era_out, "era_stratified.csv"), file.path(era_out, "era_interaction.csv")
))
run_step("10_subgroup", "10_male_subgroup_global_interactions.R", c(
  analytic_rds, file.path(subgroup_out, "male_subgroup_global_interactions.csv")
))
wqs_out <- file.path(mixture_out, "wqs")
qgcomp_out <- file.path(mixture_out, "qgcomp")
rh_out <- file.path(mixture_out, "repeated_holdout")
run_step("11_wqs_validate", "11_wqs_full_sample.R", c(analytic_rds, wqs_out, "validate"))
mixture_rds <- file.path(wqs_out, "mixture_analysis_sample_2005_2018.rds")
run_step("12_qgcomp_validate", "12_qgcomp.R", c(mixture_rds, qgcomp_out, "validate"))
run_step("13_repeated_holdout_validate", "13_wqs_repeated_holdout.R", c(mixture_rds, wqs_out, rh_out, "validate"))
cat("DETERMINISTIC_PIPELINE_STATUS=PASS\n")
