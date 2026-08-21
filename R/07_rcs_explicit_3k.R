options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "Usage: Rscript 07_rcs_explicit_3k.R",
    "<analytic_rds> <combined_csv> <female_csv> <male_csv>",
    "<knots_qc_csv> <warnings_csv>"
  ))
}

suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(rms))

input_rds <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
combined_csv <- args[2]
female_csv <- args[3]
male_csv <- args[4]
knots_qc_csv <- args[5]
warnings_csv <- args[6]
for (path in c(combined_csv, female_csv, male_csv, knots_qc_csv, warnings_csv)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

dat <- readRDS(input_rds)
if (!is.data.frame(dat) || nrow(dat) != 9129L) stop("Current V5 identity/row-count assertion failed")

exposures <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
covariates <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)
required <- unique(c(
  "stroke", "sex", "SDMVPSU", "SDMVSTRA", "WTMEC_comb",
  exposures, covariates
))
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))
if (anyNA(dat$WTMEC_comb) || any(dat$WTMEC_comb <= 0)) stop("Invalid current V5 survey weights")

design_all <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_comb,
  data = dat, nest = TRUE
)
designs <- list(
  Female = subset(design_all, sex == "Female"),
  Male = subset(design_all, sex == "Male")
)
data_groups <- list(
  Female = dat[dat$sex == "Female", , drop = FALSE],
  Male = dat[dat$sex == "Male", , drop = FALSE]
)

capture_warnings <- function(expr) {
  observed <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      observed <<- c(observed, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(observed))
}

is_default_knot_warning <- function(messages) {
  grepl("defaulting to[[:space:]]+[0-9]+[[:space:]]+knots", messages, ignore.case = TRUE)
}

all_results <- list()
all_knots_qc <- list()
all_warnings <- list()

run_one <- function(sex_name, design, data_df, exposure) {
  model_vars <- unique(c(
    "stroke", exposure, covariates,
    "SDMVPSU", "SDMVSTRA", "WTMEC_comb"
  ))
  model_complete <- complete.cases(data_df[, model_vars, drop = FALSE]) &
    is.finite(data_df[[exposure]]) &
    is.finite(data_df$WTMEC_comb) & data_df$WTMEC_comb > 0
  model_n <- sum(model_complete)
  if (model_n < 1L) stop("No model-complete observations: ", sex_name, " / ", exposure)

  exposure_values <- data_df[[exposure]][model_complete]
  knots <- as.numeric(quantile(
    exposure_values,
    probs = c(0.10, 0.50, 0.90),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
  if (length(knots) != 3L || any(!is.finite(knots)) || any(diff(knots) <= 0)) {
    stop("Invalid/non-increasing explicit knots: ", sex_name, " / ", exposure)
  }

  basis_capture <- capture_warnings(rms::rcs(exposure_values, parms = knots))
  if (any(is_default_knot_warning(basis_capture$warnings))) {
    stop("Default-knot warning during basis QC: ", sex_name, " / ", exposure)
  }
  basis <- basis_capture$value
  basis_knots <- as.numeric(attr(basis, "parms"))
  basis_nonlinear <- as.logical(attr(basis, "nonlinear"))
  if (ncol(basis) != 2L) stop("Explicit 3-knot basis did not produce two columns: ", sex_name, " / ", exposure)
  if (!isTRUE(all.equal(basis_knots, knots, tolerance = 0))) {
    stop("Basis knot locations differ from supplied knots: ", sex_name, " / ", exposure)
  }
  if (!identical(basis_nonlinear, c(FALSE, TRUE))) {
    stop("Linear/nonlinear basis mapping is ambiguous: ", sex_name, " / ", exposure)
  }

  knot_text <- paste(format(knots, digits = 17, scientific = FALSE, trim = TRUE), collapse = ",")
  spline_call <- paste0("rms::rcs(", exposure, ", parms=c(", knot_text, "))")
  formula_text <- paste(
    "stroke ~", spline_call, "+",
    paste(covariates, collapse = " + ")
  )

  fit_capture <- capture_warnings(
    svyglm(
      as.formula(formula_text), design = design,
      family = quasibinomial(link = "logit")
    )
  )
  if (any(is_default_knot_warning(fit_capture$warnings))) {
    stop("Default-knot warning during model fit: ", sex_name, " / ", exposure)
  }
  fit <- fit_capture$value
  if (as.integer(nobs(fit)) != as.integer(model_n)) {
    stop(
      "Model N differs from knot-sample N: ", sex_name, " / ", exposure,
      "; knot N=", model_n, "; model N=", nobs(fit)
    )
  }

  coefficient_names <- names(coef(fit))
  spline_terms <- grep("rcs\\(", coefficient_names, value = TRUE)
  spline_terms <- spline_terms[grepl(exposure, spline_terms, fixed = TRUE)]
  if (length(spline_terms) != 2L) {
    stop(
      "Explicit 3-knot coefficient structure is ambiguous: ", sex_name, " / ", exposure,
      "; terms=", paste(spline_terms, collapse = " | ")
    )
  }
  linear_term <- spline_terms[1]
  nonlinear_terms <- spline_terms[-1]
  if (length(nonlinear_terms) != 1L) stop("Expected exactly one nonlinear spline term")

  covariance <- vcov(fit)
  beta_nonlinear <- coef(fit)[nonlinear_terms]
  covariance_nonlinear <- covariance[nonlinear_terms, nonlinear_terms, drop = FALSE]
  beta_all <- coef(fit)[spline_terms]
  covariance_all <- covariance[spline_terms, spline_terms, drop = FALSE]
  if (any(!is.finite(beta_all)) || any(!is.finite(covariance_all))) {
    stop("Non-finite RCS coefficient/covariance: ", sex_name, " / ", exposure)
  }

  chi2_nonlinear <- as.numeric(
    t(beta_nonlinear) %*% solve(covariance_nonlinear) %*% beta_nonlinear
  )
  chi2_overall <- as.numeric(t(beta_all) %*% solve(covariance_all) %*% beta_all)
  p_nonlinear <- pchisq(chi2_nonlinear, df = length(nonlinear_terms), lower.tail = FALSE)
  p_overall <- pchisq(chi2_overall, df = length(spline_terms), lower.tail = FALSE)
  if (!is.finite(p_nonlinear) || !is.finite(p_overall)) {
    stop("Non-finite RCS P value: ", sex_name, " / ", exposure)
  }

  observed_warnings <- unique(c(basis_capture$warnings, fit_capture$warnings))
  if (length(observed_warnings)) {
    all_warnings[[length(all_warnings) + 1L]] <<- data.frame(
      Sex = sex_name,
      Metabolite = exposure,
      Warning = observed_warnings,
      stringsAsFactors = FALSE
    )
  }

  result <- data.frame(
    Sex = sex_name,
    Metabolite = exposure,
    N = as.integer(model_n),
    P_overall = p_overall,
    P_nonlinear = p_nonlinear,
    DF_overall = length(spline_terms),
    DF_nonlinear = length(nonlinear_terms),
    stringsAsFactors = FALSE
  )

  knot_qc <- data.frame(
    Sex = sex_name,
    Metabolite = exposure,
    N = as.integer(model_n),
    Knot10 = knots[1],
    Knot50 = knots[2],
    Knot90 = knots[3],
    Number_of_knots = 3L,
    Basis_columns = ncol(basis),
    Basis_nonlinear_map = paste(basis_nonlinear, collapse = " | "),
    Linear_component = linear_term,
    Nonlinear_components = paste(nonlinear_terms, collapse = " | "),
    Implementation = "rms::rcs(exposure, parms=c(empirical p10,p50,p90)); ordinary unweighted quantile(type=7) in sex-specific model-complete sample",
    Formula = formula_text,
    stringsAsFactors = FALSE
  )

  list(result = result, knot_qc = knot_qc)
}

for (sex_name in c("Female", "Male")) {
  for (exposure in exposures) {
    fitted <- run_one(
      sex_name = sex_name,
      design = designs[[sex_name]],
      data_df = data_groups[[sex_name]],
      exposure = exposure
    )
    all_results[[length(all_results) + 1L]] <- fitted$result
    all_knots_qc[[length(all_knots_qc) + 1L]] <- fitted$knot_qc
  }
}

results <- do.call(rbind, all_results)
knots_qc <- do.call(rbind, all_knots_qc)
rownames(results) <- NULL
rownames(knots_qc) <- NULL

if (nrow(results) != 20L) stop("RCS 20-model completion assertion failed")
if (nrow(knots_qc) != 20L || any(knots_qc$Number_of_knots != 3L)) stop("RCS knot-QC completion assertion failed")
if (anyNA(results) || anyNA(knots_qc[, c("Knot10", "Knot50", "Knot90")])) stop("RCS output completeness assertion failed")
if (any(knots_qc$Basis_columns != 2L) || any(knots_qc$Basis_nonlinear_map != "FALSE | TRUE")) {
  stop("RCS basis mapping assertion failed")
}

warnings_df <- if (length(all_warnings)) {
  do.call(rbind, all_warnings)
} else {
  data.frame(Sex = character(), Metabolite = character(), Warning = character(), stringsAsFactors = FALSE)
}
if (nrow(warnings_df) && any(is_default_knot_warning(warnings_df$Warning))) {
  stop("Default-knot warning detected in accumulated warnings")
}

write.csv(results, combined_csv, row.names = FALSE)
write.csv(results[results$Sex == "Female", , drop = FALSE], female_csv, row.names = FALSE)
write.csv(results[results$Sex == "Male", , drop = FALSE], male_csv, row.names = FALSE)
write.csv(knots_qc, knots_qc_csv, row.names = FALSE)
write.csv(warnings_df, warnings_csv, row.names = FALSE)

cat("RCS_STATUS=PASS\n")
cat("RCS_MODELS_COMPLETED=", nrow(results), "\n", sep = "")
cat("RCS_EXPLICIT_THREE_KNOTS=", sum(knots_qc$Number_of_knots == 3L), "\n", sep = "")
cat("RCS_DEFAULT_KNOT_WARNINGS=0\n")
cat("RCS_NONLINEAR_P_LT_0_05=", sum(results$P_nonlinear < 0.05), "\n", sep = "")
cat("RCS_OVERALL_P_LT_0_05=", sum(results$P_overall < 0.05), "\n", sep = "")
