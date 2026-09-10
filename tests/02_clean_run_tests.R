options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 02_clean_run_tests.R <workspace_root> <output_root>")
workspace <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)

anchors <- read.csv(file.path(output, "dataset", "phase_a_anchor_results.csv"), check.names = FALSE)
if (nrow(anchors) != 17L || !all(anchors$pass)) stop("Participant anchors failed")
dat <- readRDS(file.path(workspace, "data_derived", "NPS_DATA_MAIN_analytic_V5.rds"))
if (nrow(dat) != 9129L || sum(dat$stroke) != 375L ||
    sum(dat$sex == "Female") != 4611L || sum(dat$sex == "Male") != 4518L) {
  stop("Final sample anchors failed")
}
if (anyNA(dat$WTMEC_comb) || any(dat$WTMEC_comb <= 0) || !identical(dat$WTMEC_comb, dat$WTMEC2YR / 8)) {
  stop("Survey-weight assertion failed")
}
female <- read.csv(file.path(output, "primary", "NPS_OUT_TBL_MainRegression_Female_V5.csv"), check.names = FALSE)
male <- read.csv(file.path(output, "primary", "NPS_OUT_TBL_MainRegression_Male_V5.csv"), check.names = FALSE)
if (nrow(female) != 30L || nrow(male) != 30L || any(female$Model != rep(paste0("Model", 1:3), 10L))) {
  stop("Primary output row/schema assertion failed")
}
fdr <- read.csv(file.path(output, "primary", "primary_fdr_20test.csv"), check.names = FALSE)
if (nrow(fdr) != 20L || any(!is.finite(fdr$q_BH))) stop("Primary FDR family assertion failed")
common_models <- read.csv(file.path(output, "secondary", "common_cycle_models.csv"), check.names = FALSE)
common_fdr <- read.csv(file.path(output, "secondary", "common_cycle_fdr.csv"), check.names = FALSE)
if (nrow(common_models) != 20L || nrow(common_fdr) != 40L ||
    sum(common_fdr$view == "original_order") != 20L ||
    any(common_fdr$m != 20L)) stop("Common-cycle family assertion failed")
grouped <- read.csv(file.path(output, "grouped", "grouped_results.csv"), check.names = FALSE)
if (nrow(grouped) != 9L || !identical(grouped$n, c(9129L,4518L,4611L,8311L,4132L,4179L,9129L,4518L,4611L))) {
  stop("Grouped output assertion failed")
}
rcs <- read.csv(file.path(output, "rcs", "rcs_combined.csv"), check.names = FALSE)
knots <- read.csv(file.path(output, "rcs", "rcs_knots_qc.csv"), check.names = FALSE)
if (nrow(rcs) != 20L || nrow(knots) != 20L || any(knots$Number_of_knots != 3L)) stop("RCS assertion failed")
interaction <- read.csv(file.path(output, "secondary", "sex_interaction.csv"), check.names = FALSE)
era_stratified <- read.csv(file.path(output, "era", "era_stratified.csv"), check.names = FALSE)
era_interaction <- read.csv(file.path(output, "era", "era_interaction.csv"), check.names = FALSE)
subgroup <- read.csv(file.path(output, "subgroup", "male_subgroup_global_interactions.csv"), check.names = FALSE)
if (!nrow(interaction) || !nrow(era_stratified) || !nrow(era_interaction) || !nrow(subgroup)) {
  stop("Secondary deterministic output assertion failed")
}
wqs_flow <- read.csv(file.path(output, "mixture", "wqs", "mixture_sample_flow.csv"), check.names = FALSE)
if (!identical(wqs_flow$n[wqs_flow$stage %in% c("Female", "Male")], c(4179L, 4132L))) stop("Mixture anchors failed")
partition_hashes <- read.csv(file.path(output, "mixture", "repeated_holdout", "partition_hashes.csv"), check.names = FALSE)
if (nrow(partition_hashes) != 2L || !all(partition_hashes$object_identity_pass)) {
  stop("Repeated-holdout partition identity failed")
}
cat("CLEAN_RUN_TESTS=PASS\n")

# BEGIN S3 TESTS: append-only; all existing assertions above remain intact.
s3_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(s3_script_arg) != 1L) stop("Cannot locate S3 reference relative to test script")
s3_package_root <- dirname(dirname(normalizePath(
  sub("^--file=", "", s3_script_arg), winslash = "/", mustWork = TRUE)))
s3 <- read.csv(file.path(output, "dataset", "s3_four_variable_summary.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
s3_ref <- read.csv(file.path(s3_package_root, "config", "s3-four-variable-reference.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
s3_schema <- c("variable", "group", "group_n", "positive_n", "nonmissing_n",
               "missing_n", "prevalence", "SMD")
s3_vars <- c("hypertension", "diabetes", "hyperlipidemia", "smoking")
s3_groups <- c("analytic", "excluded_creatinine_unavailable")
for (z in list(s3, s3_ref)) {
  if (!identical(names(z), s3_schema) || nrow(z) != 8L ||
      !setequal(z$variable, s3_vars) || !setequal(z$group, s3_groups) ||
      anyDuplicated(paste(z$variable, z$group)) || anyNA(z)) {
    stop("S3 schema, group, variable, or privacy assertion failed")
  }
}
s3_key <- paste(s3_ref$variable, s3_ref$group)
s3 <- s3[match(s3_key, paste(s3$variable, s3$group)), , drop = FALSE]
for (nm in c("group_n", "positive_n", "nonmissing_n", "missing_n")) {
  if (any(s3[[nm]] != s3_ref[[nm]])) stop("S3 count mismatch: ", nm)
}
for (nm in c("prevalence", "SMD")) {
  if (any(!is.finite(s3[[nm]])) || any(abs(s3[[nm]] - s3_ref[[nm]]) > 1e-12)) {
    stop("S3 full-precision mismatch: ", nm)
  }
}
if (any(s3$group_n != s3$nonmissing_n + s3$missing_n) ||
    any(abs(s3$prevalence - s3$positive_n / s3$nonmissing_n) > 1e-12)) {
  stop("S3 denominator assertion failed")
}
cat("S3_FOUR_VARIABLE_TESTS=PASS\n")
# END S3 TESTS
