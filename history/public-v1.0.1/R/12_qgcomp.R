# Public release-candidate implementation of the registered unweighted qgcomp analysis.
options(stringsAsFactors = FALSE, digits = 17, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[[3]] %in% c("validate", "run")) {
  stop("Usage: Rscript 12_qgcomp.R <mixture_rds> <output_directory> <validate|run>")
}
mixture_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
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
frozen <- readRDS(mixture_rds)
required <- c("SEQN", mixture_vars, covars, "stroke", "sex")
missing_required <- setdiff(required, names(frozen))
if (length(missing_required)) stop("Missing mixture fields: ", paste(missing_required, collapse = ", "))
actual_counts <- c(
  total_n = nrow(frozen), total_events = sum(frozen$stroke),
  female_n = sum(frozen$sex == "Female"), female_events = sum(frozen$stroke[frozen$sex == "Female"]),
  male_n = sum(frozen$sex == "Male"), male_events = sum(frozen$stroke[frozen$sex == "Male"])
)
expected_counts <- c(total_n = 8311L, total_events = 339L, female_n = 4179L,
                     female_events = 175L, male_n = 4132L, male_events = 164L)
if (!identical(actual_counts, expected_counts) || anyNA(frozen[, required, drop = FALSE])) {
  stop("qgcomp frozen-sample assertion failed")
}
specification <- data.frame(
  item = c("analysis_role", "sampling_design", "q", "noboot_seed", "boot_B", "boot_seed", "sex_models"),
  value = c(analysis_label, "none", "4", "2026", "1000", "2026", "Female and Male"),
  stringsAsFactors = FALSE
)
write.csv(specification, file.path(output_dir, "qgcomp_specification.csv"), row.names = FALSE)
if (mode == "validate") {
  cat("QGCOMP_VALIDATION_STATUS=PASS\n")
  quit(status = 0, save = "no")
}

if (!requireNamespace("qgcomp", quietly = TRUE)) stop("Missing required package: qgcomp")
suppressPackageStartupMessages(library(qgcomp))

extract_nobs <- function(object, fallback) {
  candidates <- list(object$fit, tryCatch(object$fit$fit, error = function(e) NULL))
  for (candidate in candidates) {
    if (is.null(candidate)) next
    n <- tryCatch(stats::nobs(candidate), error = function(e) NA_real_)
    if (length(n) == 1L && is.finite(n)) return(as.integer(n))
  }
  as.integer(fallback)
}
extract_convergence <- function(object) {
  candidates <- list(
    tryCatch(object$fit$converged, error = function(e) NULL),
    tryCatch(object$fit$fit$converged, error = function(e) NULL)
  )
  value <- Filter(Negate(is.null), candidates)
  if (!length(value)) return(NA)
  isTRUE(value[[1]])
}
expand_weights <- function(weight_vector) {
  out <- setNames(rep(0, length(mixture_vars)), mixture_vars)
  if (!length(weight_vector)) return(out)
  if (is.null(names(weight_vector))) stop("qgcomp weights are not named")
  unknown <- setdiff(names(weight_vector), mixture_vars)
  if (length(unknown)) stop("Unknown qgcomp weight names: ", paste(unknown, collapse = ", "))
  out[names(weight_vector)] <- as.numeric(weight_vector)
  out
}

noboot_effects <- list()
noboot_weights <- list()
boot_effects <- list()
diagnostics <- list()
for (i in seq_along(c("Female", "Male"))) {
  sex_value <- c("Female", "Male")[[i]]
  df_sex <- frozen[frozen$sex == sex_value, , drop = FALSE]
  model_n <- nrow(df_sex)
  events <- sum(df_sex$stroke)
  fml <- as.formula(paste("stroke ~", paste(c(mixture_vars, covars), collapse = " + ")))

  set.seed(2026)
  fit_nb <- qgcomp.noboot(
    f = fml, data = df_sex, expnms = mixture_vars,
    family = binomial(), q = 4
  )
  psi_nb <- as.numeric(fit_nb$psi[[1]])
  se_nb <- sqrt(as.numeric(fit_nb$var.psi[[1]]))
  ci_nb <- psi_nb + c(-1, 1) * 1.96 * se_nb
  pos <- expand_weights(fit_nb$pos.weights)
  neg <- expand_weights(fit_nb$neg.weights)
  pos_sum <- sum(pos)
  neg_sum <- sum(neg)
  weight_qc <- (sum(pos > 0) == 0L || abs(pos_sum - 1) <= 1e-8) &&
    (sum(neg > 0) == 0L || abs(neg_sum - 1) <= 1e-8) &&
    all(is.finite(c(pos, neg))) && all(c(pos, neg) >= 0)
  if (!weight_qc) stop("qgcomp noboot weight assertion failed")
  noboot_effects[[i]] <- data.frame(
    analysis_label = analysis_label, method = "qgcomp.noboot", sex = sex_value,
    psi = psi_nb, standard_error = se_nb, OR = exp(psi_nb),
    CI_lower = exp(ci_nb[[1]]), CI_upper = exp(ci_nb[[2]]),
    P_value = 2 * pnorm(-abs(psi_nb / se_nb)), n = model_n,
    model_object_n = extract_nobs(fit_nb, model_n), events = events,
    q = 4, seed = 2026, family = "binomial(logit)", CI_method = "Wald",
    sampling_design = "none", stringsAsFactors = FALSE
  )
  noboot_weights[[i]] <- data.frame(
    analysis_label = analysis_label, method = "qgcomp.noboot", sex = sex_value,
    component = mixture_vars, positive_weight = unname(pos), negative_weight = unname(neg),
    net_weight = unname(pos - neg), positive_weight_sum = pos_sum,
    negative_weight_sum = neg_sum, stringsAsFactors = FALSE
  )
  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    method = "QGCOMP_NOBOOT", sex = sex_value, status = "SUCCESS", n = model_n,
    events = events, q = 4, B_requested = 0, B_actual = 0, seed = 2026,
    convergence = extract_convergence(fit_nb), weight_qc = weight_qc,
    stringsAsFactors = FALSE
  )
  saveRDS(fit_nb, file.path(output_dir, paste0("QGCOMP_", sex_value, "_NoBoot.rds")), version = 3)

  set.seed(2026)
  fit_boot <- qgcomp.boot(
    f = fml, data = df_sex, expnms = mixture_vars,
    family = binomial(), q = 4, B = 1000, seed = 2026
  )
  psi_boot <- as.numeric(fit_boot$psi[[1]])
  se_boot <- sqrt(as.numeric(fit_boot$var.psi[[1]]))
  ci_boot <- psi_boot + c(-1, 1) * 1.96 * se_boot
  boot_effects[[i]] <- data.frame(
    analysis_label = analysis_label, method = "qgcomp.boot", sex = sex_value,
    psi = psi_boot, standard_error = se_boot, OR = exp(psi_boot),
    CI_lower = exp(ci_boot[[1]]), CI_upper = exp(ci_boot[[2]]),
    P_value = 2 * pnorm(-abs(psi_boot / se_boot)), n = model_n,
    model_object_n = extract_nobs(fit_boot, model_n), events = events,
    q = 4, B_requested = 1000, B_actual = 1000, seed = 2026,
    family = "binomial(logit)", CI_method = "bootstrap variance normal approximation",
    sampling_design = "none", stringsAsFactors = FALSE
  )
  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    method = "QGCOMP_BOOT", sex = sex_value, status = "SUCCESS", n = model_n,
    events = events, q = 4, B_requested = 1000, B_actual = 1000, seed = 2026,
    convergence = extract_convergence(fit_boot), weight_qc = NA,
    stringsAsFactors = FALSE
  )
  saveRDS(fit_boot, file.path(output_dir, paste0("QGCOMP_", sex_value, "_Boot.rds")), version = 3)
}
noboot_effects <- do.call(rbind, noboot_effects)
noboot_weights <- do.call(rbind, noboot_weights)
boot_effects <- do.call(rbind, boot_effects)
diagnostics <- do.call(rbind, diagnostics)
if (nrow(noboot_effects) != 2L || nrow(noboot_weights) != 20L || nrow(boot_effects) != 2L) {
  stop("qgcomp completion assertion failed")
}
write.csv(noboot_effects, file.path(output_dir, "04_QGCOMP_NOBOOT_EFFECTS.csv"), row.names = FALSE)
write.csv(noboot_weights, file.path(output_dir, "05_QGCOMP_NOBOOT_COMPONENT_WEIGHTS.csv"), row.names = FALSE)
write.csv(boot_effects, file.path(output_dir, "06_QGCOMP_BOOT_EFFECTS.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(output_dir, "07_QGCOMP_MODEL_DIAGNOSTICS.csv"), row.names = FALSE)
cat("QGCOMP_RUN_STATUS=PASS\n")
