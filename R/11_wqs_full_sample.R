# Public release-candidate implementation of the registered unweighted WQS analysis.
# The validation mode constructs and checks the frozen mixture sample without fitting models.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[3]] %in% c("validate", "run")) {
  stop("Usage: Rscript 11_wqs_full_sample.R <analytic_rds> <output_directory> <validate|run>")
}
input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- gsub("\\\\", "/", args[[2]])
mode <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_label <- "UNWEIGHTED EXPLORATORY MIXTURE ANALYSIS"
mixture_vars <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
covars <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f", "hyperlipidemia_f"
)
eligible_cycles <- c(
  "2005_2006", "2007_2008", "2009_2010", "2011_2012",
  "2013_2014", "2015_2016", "2017_2018"
)
required_fields <- c("SEQN", "cycle", mixture_vars, "stroke", "sex", covars)

analytic_main <- readRDS(input_rds)
missing_fields <- setdiff(required_fields, names(analytic_main))
if (length(missing_fields)) stop("Required fields missing: ", paste(missing_fields, collapse = ", "))
eligible <- analytic_main[analytic_main$cycle %in% eligible_cycles, required_fields, drop = FALSE]
complete_flag <- complete.cases(eligible[, c(mixture_vars, "stroke", "sex", covars), drop = FALSE])
frozen <- eligible[complete_flag, , drop = FALSE]
frozen$stroke <- as.integer(frozen$stroke)

actual_counts <- c(
  total_n = nrow(frozen), total_events = sum(frozen$stroke),
  female_n = sum(frozen$sex == "Female"),
  female_events = sum(frozen$stroke[frozen$sex == "Female"]),
  male_n = sum(frozen$sex == "Male"),
  male_events = sum(frozen$stroke[frozen$sex == "Male"])
)
expected_counts <- c(
  total_n = 8311L, total_events = 339L, female_n = 4179L,
  female_events = 175L, male_n = 4132L, male_events = 164L
)
if (!identical(actual_counts, expected_counts)) {
  print(rbind(expected = expected_counts, actual = actual_counts))
  stop("Mixture sample count assertion failed")
}
if (anyNA(frozen[, c(mixture_vars, "stroke", "sex", covars), drop = FALSE])) {
  stop("Mixture sample contains missing model fields")
}
if (anyDuplicated(frozen$SEQN) || !identical(sort(unique(frozen$stroke)), c(0L, 1L))) {
  stop("Mixture sample participant/outcome assertion failed")
}

sample_flow <- data.frame(
  stage = c("analytic_sample", "eligible_cycles_2005_2018", "complete_case", "Female", "Male"),
  n = c(nrow(analytic_main), nrow(eligible), nrow(frozen), expected_counts[["female_n"]], expected_counts[["male_n"]]),
  events = c(sum(analytic_main$stroke), sum(eligible$stroke), sum(frozen$stroke),
             expected_counts[["female_events"]], expected_counts[["male_events"]]),
  stringsAsFactors = FALSE
)
write.csv(sample_flow, file.path(output_dir, "mixture_sample_flow.csv"), row.names = FALSE)
specification <- data.frame(
  item = c("analysis_role", "sampling_design", "q", "B", "seed", "validation", "signal", "models"),
  value = c(analysis_label, "none", "4", "1000", "2026", "0", "t2", "Female/Male x Positive/Negative"),
  stringsAsFactors = FALSE
)
write.csv(specification, file.path(output_dir, "wqs_specification.csv"), row.names = FALSE)
frozen_path <- file.path(output_dir, "mixture_analysis_sample_2005_2018.rds")
saveRDS(frozen, frozen_path, version = 3)

if (mode == "validate") {
  cat("WQS_VALIDATION_STATUS=PASS\n")
  cat("MIXTURE_ROWS=8311\n")
  quit(status = 0, save = "no")
}

required_packages <- c("gWQS", "future")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(gWQS))
future::plan("sequential")

model_specs <- data.frame(
  sex = c("Female", "Female", "Male", "Male"),
  direction = c("Positive", "Negative", "Positive", "Negative"),
  b1_pos = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
effects <- list()
weights <- list()
diagnostics <- list()
for (i in seq_len(nrow(model_specs))) {
  spec <- model_specs[i, ]
  df_sex <- frozen[frozen$sex == spec$sex, , drop = FALSE]
  fml <- as.formula(paste("stroke ~", paste(c("wqs", covars), collapse = " + ")))
  observed_warnings <- character()
  set.seed(2026)
  elapsed <- system.time({
    fit <- withCallingHandlers(
      gwqs(
        formula = fml, mix_name = mixture_vars, data = df_sex,
        q = 4, validation = 0, b = 1000, b1_pos = spec$b1_pos,
        family = binomial(link = "logit"), seed = 2026, signal = "t2",
        plan_strategy = "sequential"
      ),
      warning = function(w) {
        observed_warnings <<- c(observed_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })[["elapsed"]]
  coefficient_table <- summary(fit)$coefficients
  beta <- unname(coefficient_table["wqs", "Estimate"])
  se <- unname(coefficient_table["wqs", "Std. Error"])
  p_value <- unname(coefficient_table["wqs", "Pr(>|z|)"])
  w <- fit$final_weights
  if (!all(c("mix_name", "mean_weight") %in% names(w))) stop("Unexpected WQS weight structure")
  weight_sum <- sum(w$mean_weight)
  weight_qc <- all(is.finite(w$mean_weight)) && all(w$mean_weight >= 0) && abs(weight_sum - 1) <= 1e-8
  if (!weight_qc) stop("WQS weight assertion failed")
  effects[[i]] <- data.frame(
    analysis_label = analysis_label, method = "WQS", sex = spec$sex,
    direction = spec$direction, coefficient = beta, standard_error = se,
    OR = exp(beta), CI_lower = exp(beta - 1.96 * se), CI_upper = exp(beta + 1.96 * se),
    P_value = p_value, n = nrow(df_sex), events = sum(df_sex$stroke),
    q = 4, B = 1000, seed = 2026, validation = 0, signal = "t2",
    family = "binomial(logit)", sampling_design = "none",
    elapsed_seconds = unname(elapsed), stringsAsFactors = FALSE
  )
  weight_frame <- data.frame(
    analysis_label = analysis_label, method = "WQS", sex = spec$sex,
    direction = spec$direction, component = w$mix_name, weight = w$mean_weight,
    weight_sum = weight_sum, stringsAsFactors = FALSE
  )
  weight_frame$rank <- rank(-weight_frame$weight, ties.method = "first")
  weights[[i]] <- weight_frame
  diagnostics[[i]] <- data.frame(
    sex = spec$sex, direction = spec$direction, status = "SUCCESS",
    n = nrow(df_sex), events = sum(df_sex$stroke), weight_qc = weight_qc,
    warnings = paste(unique(observed_warnings), collapse = " | "), stringsAsFactors = FALSE
  )
  saveRDS(fit, file.path(output_dir, paste0("WQS_", spec$sex, "_", spec$direction, ".rds")), version = 3)
}
effects <- do.call(rbind, effects)
weights <- do.call(rbind, weights)
diagnostics <- do.call(rbind, diagnostics)
if (nrow(effects) != 4L || nrow(weights) != 40L || !all(diagnostics$weight_qc)) stop("WQS completion assertion failed")
write.csv(effects, file.path(output_dir, "01_WQS_ALL_DIRECTIONS_EFFECTS.csv"), row.names = FALSE)
write.csv(weights, file.path(output_dir, "02_WQS_ALL_COMPONENT_WEIGHTS.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "03_WQS_MODEL_DIAGNOSTICS.csv"), row.names = FALSE)
cat("WQS_RUN_STATUS=PASS\n")
