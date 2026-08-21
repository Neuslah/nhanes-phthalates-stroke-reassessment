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
