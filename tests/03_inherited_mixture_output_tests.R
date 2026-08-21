options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 03_inherited_mixture_output_tests.R <full_sample_results_dir> <repeated_holdout_summary_dir>")
}
full_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
rh_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
read_table <- function(directory, name) read.csv(file.path(directory, name), check.names = FALSE, stringsAsFactors = FALSE)
require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(label, " schema missing: ", paste(missing, collapse = ", "))
}

wqs_effects <- read_table(full_dir, "01_WQS_ALL_DIRECTIONS_EFFECTS.csv")
wqs_weights <- read_table(full_dir, "02_WQS_ALL_COMPONENT_WEIGHTS.csv")
wqs_diagnostics <- read_table(full_dir, "03_WQS_MODEL_DIAGNOSTICS.csv")
require_columns(wqs_effects, c("sex", "direction", "OR", "P_value", "n", "events", "q", "B", "seed", "validation", "signal", "sampling_design"), "WQS effects")
require_columns(wqs_weights, c("sex", "direction", "component", "weight", "weight_sum", "rank"), "WQS weights")
if (nrow(wqs_effects) != 4L || nrow(wqs_weights) != 40L || nrow(wqs_diagnostics) != 4L ||
    !identical(wqs_effects$n, c(4179L, 4179L, 4132L, 4132L)) ||
    !identical(wqs_effects$events, c(175L, 175L, 164L, 164L)) ||
    any(wqs_effects$q != 4L) || any(wqs_effects$B != 1000L) || any(wqs_effects$seed != 2026L) ||
    any(wqs_effects$validation != 0) || any(wqs_effects$signal != "t2") ||
    any(wqs_effects$sampling_design != "none") || any(abs(wqs_weights$weight_sum - 1) > 1e-8)) {
  stop("WQS aggregate/schema anchor failure")
}

qg_nb <- read_table(full_dir, "04_QGCOMP_NOBOOT_EFFECTS.csv")
qg_weights <- read_table(full_dir, "05_QGCOMP_NOBOOT_COMPONENT_WEIGHTS.csv")
qg_boot <- read_table(full_dir, "06_QGCOMP_BOOT_EFFECTS.csv")
qg_diag <- read_table(full_dir, "07_QGCOMP_MODEL_DIAGNOSTICS.csv")
require_columns(qg_nb, c("sex", "psi", "OR", "P_value", "n", "events", "q", "seed", "sampling_design"), "qgcomp noboot")
require_columns(qg_boot, c("sex", "psi", "OR", "P_value", "n", "events", "q", "B_requested", "B_actual", "seed", "sampling_design"), "qgcomp boot")
if (nrow(qg_nb) != 2L || nrow(qg_weights) != 20L || nrow(qg_boot) != 2L || nrow(qg_diag) != 4L ||
    !identical(qg_nb$n, c(4179L, 4132L)) || !identical(qg_nb$events, c(175L, 164L)) ||
    any(qg_nb$q != 4L) || any(qg_nb$seed != 2026L) || any(qg_nb$sampling_design != "none") ||
    any(qg_boot$q != 4L) || any(qg_boot$B_requested != 1000L) || any(qg_boot$B_actual != 1000L) ||
    any(qg_boot$seed != 2026L) || any(qg_boot$sampling_design != "none")) {
  stop("qgcomp aggregate/schema anchor failure")
}

rh_effects <- read_table(rh_dir, "01_WQS_RH_EFFECT_SUMMARY.csv")
rh_weights <- read_table(rh_dir, "02_WQS_RH_WEIGHT_STABILITY.csv")
rh_top <- read_table(rh_dir, "03_WQS_RH_TOP_COMPONENT_FREQUENCY.csv")
rh_class <- read_table(rh_dir, "04_WQS_RH_CLASSIFICATION.csv")
require_columns(rh_effects, c("sex", "direction", "attempted_holdouts", "successful_holdouts", "OR_median", "OR_p2_5", "OR_p97_5", "proportion_prespecified_direction"), "Repeated-holdout effects")
require_columns(rh_class, c("sex", "direction", "successful_holdouts", "classification"), "Repeated-holdout classification")
if (nrow(rh_effects) != 4L || nrow(rh_weights) != 40L || nrow(rh_top) != 40L || nrow(rh_class) != 4L ||
    any(rh_effects$attempted_holdouts != 100L) ||
    !identical(rh_effects$successful_holdouts, c(100L, 99L, 97L, 100L)) ||
    !identical(rh_class$successful_holdouts, c(100L, 99L, 97L, 100L))) {
  stop("Repeated-holdout aggregate/schema anchor failure")
}
cat("INHERITED_MIXTURE_OUTPUT_TESTS=PASS\n")
