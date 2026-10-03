# Public release-candidate implementation of the registered repeated-holdout WQS analysis.
# Validation mode recreates the frozen partition registry but does not fit WQS models.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !args[[4]] %in% c("validate", "run")) {
  stop(paste(
    "Usage: Rscript 13_wqs_repeated_holdout.R",
    "<mixture_rds> <full_wqs_output_directory> <output_directory> <validate|run>"
  ))
}
mixture_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
full_wqs_dir <- gsub("\\\\", "/", args[[2]])
output_dir <- gsub("\\\\", "/", args[[3]])
mode <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_label <- "UNWEIGHTED EXPLORATORY MIXTURE ANALYSIS"
partition_seed <- 2026L
model_seed <- 2026L
rh_n <- 100L
b_n <- 100L
validation_fraction <- 0.60
equal_weight_threshold <- 0.10
weight_sum_tolerance <- 1e-8
mixture_vars <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
covars <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f", "hyperlipidemia_f"
)

d <- readRDS(mixture_rds)
required <- c("SEQN", mixture_vars, covars, "stroke", "sex")
missing_required <- setdiff(required, names(d))
if (length(missing_required)) stop("Missing mixture fields: ", paste(missing_required, collapse = ", "))
d$source_row_id <- seq_len(nrow(d))
actual_counts <- c(
  total_n = nrow(d), total_events = sum(d$stroke),
  female_n = sum(d$sex == "Female"), female_events = sum(d$stroke[d$sex == "Female"]),
  male_n = sum(d$sex == "Male"), male_events = sum(d$stroke[d$sex == "Male"])
)
expected_counts <- c(total_n = 8311L, total_events = 339L, female_n = 4179L,
                     female_events = 175L, male_n = 4132L, male_events = 164L)
if (!identical(actual_counts, expected_counts) || anyNA(d[, required, drop = FALSE])) {
  stop("Repeated-holdout sample assertion failed")
}
if (!requireNamespace("digest", quietly = TRUE)) stop("Missing required package: digest")
hash_file <- function(path) toupper(digest::digest(path, algo = "sha256", file = TRUE))
write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

partition_dir <- file.path(output_dir, "partitions")
dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
partition_paths <- c(
  Female = file.path(partition_dir, "female_validation_rows.rds"),
  Male = file.path(partition_dir, "male_validation_rows.rds")
)
canonical_partition_file_sha <- c(
  Female = "CC85B54CAB24442EFF9622562AD5CE9EC513B01E46754B9404F46AEDC2F27802",
  Male = "777EDBB81B7FEC8859D702DF9F24294E3851C5ACDB5F32981FDA0C03C2CE7F84"
)
expected_partition_object_digest <- c(
  Female = "2B9FDA53B1C6ABAD5D0EA0D8AB2C87FE198FAE47C69AC159DD10DE95F00C5AF0",
  Male = "B56BAFE5512CA0E9316A8F667697DD53F3F5F65D4F755707EAC816127B42F5AA"
)
partition_counts <- list()
partition_hashes <- list()
for (sex_value in c("Female", "Male")) {
  ds <- d[d$sex == sex_value, , drop = FALSE]
  n_val <- as.integer(round(validation_fraction * nrow(ds)))
  set.seed(partition_seed)
  validation_rows <- lapply(seq_len(rh_n), function(i) {
    z <- rep(FALSE, nrow(ds))
    z[sample.int(nrow(ds), size = n_val, replace = FALSE)] <- TRUE
    z
  })
  saveRDS(validation_rows, partition_paths[[sex_value]], version = 3)
  observed_file_sha <- hash_file(partition_paths[[sex_value]])
  observed_object_digest <- toupper(digest::digest(validation_rows, algo = "sha256", serialize = TRUE))
  if (!identical(observed_object_digest, expected_partition_object_digest[[sex_value]])) {
    stop("Partition object identity mismatch for ", sex_value, ": ", observed_object_digest)
  }
  partition_hashes[[sex_value]] <- data.frame(
    sex = sex_value, partition_seed = partition_seed, rh = rh_n,
    validation_fraction = validation_fraction, validation_n = n_val,
    training_n = nrow(ds) - n_val, generated_file_sha256 = observed_file_sha,
    canonical_file_sha256 = canonical_partition_file_sha[[sex_value]],
    object_digest_sha256 = observed_object_digest,
    expected_object_digest_sha256 = expected_partition_object_digest[[sex_value]],
    object_identity_pass = TRUE,
    stringsAsFactors = FALSE
  )
  partition_counts[[sex_value]] <- do.call(rbind, lapply(seq_len(rh_n), function(i) {
    validation <- validation_rows[[i]]
    training <- !validation
    data.frame(
      sex = sex_value, holdout_id = i,
      training_n = sum(training), training_events = sum(ds$stroke[training]),
      validation_n = sum(validation), validation_events = sum(ds$stroke[validation]),
      stringsAsFactors = FALSE
    )
  }))
}
partition_counts <- do.call(rbind, partition_counts)
partition_hashes <- do.call(rbind, partition_hashes)
if (nrow(partition_counts) != 200L || any(partition_counts$training_events <= 0L) ||
    any(partition_counts$validation_events <= 0L)) stop("Partition event assertion failed")
write_csv(partition_hashes, file.path(output_dir, "partition_hashes.csv"))
write_csv(partition_counts, file.path(output_dir, "partition_event_counts.csv"))
specification <- data.frame(
  item = c("analysis_role", "sampling_design", "rh", "b", "q", "partition_seed",
           "model_seed", "training_fraction", "validation_fraction", "signal", "rs", "b_constr"),
  value = c(analysis_label, "none", "100", "100", "4", "2026", "2026", "0.40",
            "0.60", "t2", "FALSE", "FALSE"),
  stringsAsFactors = FALSE
)
write_csv(specification, file.path(output_dir, "repeated_holdout_specification.csv"))
if (mode == "validate") {
  cat("REPEATED_HOLDOUT_VALIDATION_STATUS=PASS\n")
  cat("PARTITION_OBJECT_IDENTITIES=IDENTICAL\n")
  quit(status = 0, save = "no")
}

required_packages <- c("gWQS", "future")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(gWQS))
future::plan("sequential")

model_specs <- data.frame(
  key = c("Female_Positive", "Female_Negative", "Male_Positive", "Male_Negative"),
  sex = c("Female", "Female", "Male", "Male"),
  direction = c("Positive", "Negative", "Positive", "Negative"),
  b1_pos = c(TRUE, FALSE, TRUE, FALSE), stringsAsFactors = FALSE
)
capture_model <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)), class = "controlled_wqs_error")
  )
  list(value = value, warnings = unique(warnings))
}
map_child_to_holdout <- function(child_validation, validation_rows) {
  hits <- which(vapply(validation_rows, function(v) {
    length(v) == length(child_validation) && all(v == child_validation)
  }, logical(1)))
  if (length(hits) != 1L) NA_integer_ else hits
}

effect_parts <- list()
weight_parts <- list()
for (i in seq_len(nrow(model_specs))) {
  spec <- model_specs[i, ]
  ds <- d[d$sex == spec$sex, , drop = FALSE]
  validation_rows <- readRDS(partition_paths[[spec$sex]])
  fml <- as.formula(paste("stroke ~", paste(c("wqs", covars), collapse = " + ")))
  set.seed(model_seed)
  captured <- capture_model(gwqs(
    formula = fml, data = ds, mix_name = mixture_vars,
    rh = rh_n, b = b_n, b1_pos = spec$b1_pos,
    b_constr = FALSE, q = 4, validation = validation_fraction,
    validation_rows = validation_rows, family = binomial(link = "logit"),
    signal = "t2", rs = FALSE, seed = model_seed,
    plan_strategy = "sequential"
  ))
  fit <- captured$value
  if (inherits(fit, "controlled_wqs_error")) stop("Repeated-holdout model failed: ", spec$key, ": ", fit$message)
  saveRDS(fit, file.path(output_dir, paste0("WQS_RH_", spec$key, ".rds")), version = 3)
  map_ids <- vapply(fit$gwqslist, function(x) map_child_to_holdout(x$validation_rows, validation_rows), integer(1))
  if (anyNA(map_ids) || anyDuplicated(map_ids)) stop("Model-to-partition mapping failure for ", spec$key)
  counts <- partition_counts[partition_counts$sex == spec$sex, ]
  effects <- vector("list", rh_n)
  weights <- vector("list", rh_n)
  for (holdout_id in seq_len(rh_n)) {
    child_position <- match(holdout_id, map_ids)
    count_row <- counts[counts$holdout_id == holdout_id, ]
    if (is.na(child_position)) {
      effects[[holdout_id]] <- data.frame(
        analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
        holdout_id = holdout_id, beta = NA_real_, standard_error = NA_real_, OR = NA_real_,
        lower_CI = NA_real_, upper_CI = NA_real_, P = NA_real_,
        training_n = count_row$training_n, training_events = count_row$training_events,
        validation_n = count_row$validation_n, validation_events = count_row$validation_events,
        success = FALSE, stringsAsFactors = FALSE
      )
      weights[[holdout_id]] <- data.frame(
        analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
        holdout_id = holdout_id, metabolite = mixture_vars, weight = NA_real_,
        rank = NA_real_, weight_sum = NA_real_, weight_qc = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      child <- fit$gwqslist[[child_position]]
      coefficient_table <- coef(summary(child$fit))
      beta <- unname(coefficient_table["wqs", "Estimate"])
      se <- unname(coefficient_table["wqs", "Std. Error"])
      p_value <- unname(coefficient_table["wqs", grep("Pr", colnames(coefficient_table), value = TRUE)[1]])
      effects[[holdout_id]] <- data.frame(
        analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
        holdout_id = holdout_id, beta = beta, standard_error = se, OR = exp(beta),
        lower_CI = exp(beta - 1.96 * se), upper_CI = exp(beta + 1.96 * se), P = p_value,
        training_n = count_row$training_n, training_events = count_row$training_events,
        validation_n = count_row$validation_n, validation_events = count_row$validation_events,
        success = TRUE, stringsAsFactors = FALSE
      )
      w <- child$final_weights[match(mixture_vars, child$final_weights$mix_name), ]
      weight_sum <- sum(w$mean_weight)
      weight_qc <- all(is.finite(w$mean_weight)) && all(w$mean_weight >= 0) &&
        abs(weight_sum - 1) <= weight_sum_tolerance
      weights[[holdout_id]] <- data.frame(
        analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
        holdout_id = holdout_id, metabolite = w$mix_name, weight = w$mean_weight,
        rank = rank(-w$mean_weight, ties.method = "first"), weight_sum = weight_sum,
        weight_qc = weight_qc, stringsAsFactors = FALSE
      )
    }
  }
  effect_parts[[i]] <- do.call(rbind, effects)
  weight_parts[[i]] <- do.call(rbind, weights)
}
effects <- do.call(rbind, effect_parts)
weights <- do.call(rbind, weight_parts)
write_csv(effects, file.path(output_dir, "01_WQS_RH_EFFECTS_LONG.csv"))
write_csv(weights, file.path(output_dir, "02_WQS_RH_WEIGHTS_LONG.csv"))

quantile7 <- function(x, p) unname(quantile(x[is.finite(x)], p, type = 7, names = FALSE))
effect_summary <- list()
weight_summary <- list()
for (i in seq_len(nrow(model_specs))) {
  spec <- model_specs[i, ]
  e <- effects[effects$sex == spec$sex & effects$direction == spec$direction & effects$success, ]
  effect_summary[[i]] <- data.frame(
    analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
    attempted_holdouts = rh_n, successful_holdouts = nrow(e), failed_holdouts = rh_n - nrow(e),
    beta_median = median(e$beta), beta_p2_5 = quantile7(e$beta, 0.025), beta_p97_5 = quantile7(e$beta, 0.975),
    OR_median = median(e$OR), OR_p2_5 = quantile7(e$OR, 0.025), OR_p97_5 = quantile7(e$OR, 0.975),
    beta_mean = mean(e$beta), beta_SD = sd(e$beta),
    proportion_OR_lt_1 = mean(e$OR < 1), proportion_OR_gt_1 = mean(e$OR > 1),
    proportion_prespecified_direction = if (spec$b1_pos) mean(e$OR > 1) else mean(e$OR < 1),
    proportion_P_lt_0_05_descriptive = mean(e$P < 0.05), stringsAsFactors = FALSE
  )
  w <- weights[weights$sex == spec$sex & weights$direction == spec$direction & is.finite(weights$weight), ]
  for (metabolite in mixture_vars) {
    wm <- w[w$metabolite == metabolite, ]
    weight_summary[[length(weight_summary) + 1L]] <- data.frame(
      analysis_label = analysis_label, sex = spec$sex, direction = spec$direction,
      metabolite = metabolite, successful_weight_rows = nrow(wm),
      mean_weight = mean(wm$weight), median_weight = median(wm$weight), SD_weight = sd(wm$weight),
      weight_p2_5 = quantile7(wm$weight, 0.025), weight_p97_5 = quantile7(wm$weight, 0.975),
      above_0_10_count = sum(wm$weight > equal_weight_threshold),
      above_0_10_proportion = mean(wm$weight > equal_weight_threshold),
      top1_count = sum(wm$rank == 1), top1_proportion = mean(wm$rank == 1),
      top3_count = sum(wm$rank <= 3), top3_proportion = mean(wm$rank <= 3),
      rank_median = median(wm$rank), rank_min = min(wm$rank), rank_max = max(wm$rank),
      stringsAsFactors = FALSE
    )
  }
}
effect_summary <- do.call(rbind, effect_summary)
weight_summary <- do.call(rbind, weight_summary)
write_csv(effect_summary, file.path(output_dir, "03_WQS_RH_EFFECT_SUMMARY.csv"))
write_csv(weight_summary, file.path(output_dir, "04_WQS_RH_WEIGHT_STABILITY.csv"))
write_csv(weight_summary[, c("analysis_label", "sex", "direction", "metabolite", "top1_count",
                             "top1_proportion", "top3_count", "top3_proportion", "rank_median",
                             "rank_min", "rank_max")],
          file.path(output_dir, "05_WQS_RH_TOP_COMPONENT_FREQUENCY.csv"))

full_effect_path <- file.path(full_wqs_dir, "01_WQS_ALL_DIRECTIONS_EFFECTS.csv")
if (!file.exists(full_effect_path)) stop("Full-sample WQS effects file not found: ", full_effect_path)
full_effects <- read.csv(full_effect_path, check.names = FALSE, stringsAsFactors = FALSE)
full_effects <- full_effects[match(model_specs$key, paste(full_effects$sex, full_effects$direction, sep = "_")), ]
comparison <- data.frame(
  analysis_label = analysis_label, sex = model_specs$sex, direction = model_specs$direction,
  full_sample_beta = log(full_effects$OR), full_sample_OR = full_effects$OR,
  full_sample_P = full_effects$P_value, rh_median_beta = effect_summary$beta_median,
  rh_mean_beta = effect_summary$beta_mean,
  beta_absolute_difference_median = abs(effect_summary$beta_median - log(full_effects$OR)),
  rh_median_OR = effect_summary$OR_median,
  OR_absolute_difference = abs(effect_summary$OR_median - full_effects$OR),
  direction_match = sign(effect_summary$beta_median) == sign(log(full_effects$OR)),
  stringsAsFactors = FALSE
)
write_csv(comparison, file.path(output_dir, "06_FULL_SAMPLE_VS_RH_EFFECTS.csv"))
if (nrow(effects) != 400L || any(table(effects$sex, effects$direction) != 100L) ||
    any(!weights$weight_qc[is.finite(weights$weight)])) stop("Repeated-holdout completion assertion failed")
cat("REPEATED_HOLDOUT_RUN_STATUS=PASS\n")
