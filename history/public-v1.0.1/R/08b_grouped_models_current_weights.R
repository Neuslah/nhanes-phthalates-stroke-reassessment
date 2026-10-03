options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 08b_grouped_models_current_weights.R <analytic_rds> <grouped_rds> <result_csv> <model_qc_csv>")
}

suppressPackageStartupMessages(library(survey))

current_v5_rds <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
candidate_v6_rds <- args[2]
result_csv <- args[3]
model_qc_csv <- args[4]
dir.create(dirname(result_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_qc_csv), recursive = TRUE, showWarnings = FALSE)

v5 <- readRDS(current_v5_rds)
dat <- readRDS(candidate_v6_rds)
if (!identical(dat$SEQN, v5$SEQN)) stop("Candidate V6 participant/order mismatch")
if (!identical(dat$WTMEC2YR, v5$WTMEC2YR)) stop("Candidate V6 WTMEC2YR differs from current V5")
if (!identical(dat$WTMEC_comb, v5$WTMEC_comb)) stop("Candidate V6 WTMEC_comb differs from current V5")
if (!identical(dat$WTMEC_comb, dat$WTMEC2YR / 8)) stop("Candidate V6 weight divisor assertion failed")

required <- c(
  "stroke", "sex", "SDMVPSU", "SDMVSTRA", "WTMEC_comb",
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f", "ln_sigma_LMW", "ln_sigma_HMW", "ln_sigma_DEHP"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

design_full <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_comb,
  data = dat, nest = TRUE
)
design_female <- subset(design_full, sex == "Female")
design_male <- subset(design_full, sex == "Male")
design_hmw <- subset(design_full, !is.na(ln_sigma_HMW))
design_hmw_female <- subset(design_hmw, sex == "Female")
design_hmw_male <- subset(design_hmw, sex == "Male")

covariates <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)

run_one <- function(design, exposure, index, stratum) {
  formula_text <- paste("stroke ~", exposure, "+", paste(covariates, collapse = " + "))
  fit <- svyglm(
    as.formula(formula_text), design = design,
    family = quasibinomial(link = "logit")
  )
  beta <- unname(coef(fit)[exposure])
  se <- sqrt(unname(diag(vcov(fit))[exposure]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) stop("Invalid grouped model coefficient/SE: ", index, " / ", stratum)
  n_model <- as.integer(nobs(fit))
  events <- as.integer(sum(fit$y == 1, na.rm = TRUE))
  result <- data.frame(
    Label = paste(index, ifelse(stratum == "Overall", "Total", stratum), sep = "_"),
    OR = round(exp(beta), 4),
    CI_L = round(exp(beta - 1.96 * se), 4),
    CI_U = round(exp(beta + 1.96 * se), 4),
    P = round(2 * pnorm(-abs(beta / se)), 6),
    n = n_model,
    stringsAsFactors = FALSE
  )
  qc <- data.frame(
    Index = index,
    Stratum = stratum,
    Exposure = exposure,
    N = n_model,
    Events = events,
    Beta = beta,
    SE = se,
    Formula = formula_text,
    Weight = "WTMEC_comb",
    IDs = "SDMVPSU",
    Strata = "SDMVSTRA",
    Nest = TRUE,
    stringsAsFactors = FALSE
  )
  list(result = result, qc = qc)
}

specifications <- list(
  list(design = design_full, exposure = "ln_sigma_LMW", index = "LMW", stratum = "Overall"),
  list(design = design_male, exposure = "ln_sigma_LMW", index = "LMW", stratum = "Male"),
  list(design = design_female, exposure = "ln_sigma_LMW", index = "LMW", stratum = "Female"),
  list(design = design_hmw, exposure = "ln_sigma_HMW", index = "HMW", stratum = "Overall"),
  list(design = design_hmw_male, exposure = "ln_sigma_HMW", index = "HMW", stratum = "Male"),
  list(design = design_hmw_female, exposure = "ln_sigma_HMW", index = "HMW", stratum = "Female"),
  list(design = design_full, exposure = "ln_sigma_DEHP", index = "DEHP", stratum = "Overall"),
  list(design = design_male, exposure = "ln_sigma_DEHP", index = "DEHP", stratum = "Male"),
  list(design = design_female, exposure = "ln_sigma_DEHP", index = "DEHP", stratum = "Female")
)

fits <- lapply(specifications, function(spec) {
  run_one(spec$design, spec$exposure, spec$index, spec$stratum)
})
results <- do.call(rbind, lapply(fits, `[[`, "result"))
model_qc <- do.call(rbind, lapply(fits, `[[`, "qc"))
rownames(results) <- NULL
rownames(model_qc) <- NULL

expected_labels <- c(
  "LMW_Total", "LMW_Male", "LMW_Female",
  "HMW_Total", "HMW_Male", "HMW_Female",
  "DEHP_Total", "DEHP_Male", "DEHP_Female"
)
expected_n <- c(9129L, 4518L, 4611L, 8311L, 4132L, 4179L, 9129L, 4518L, 4611L)
if (!identical(results$Label, expected_labels)) stop("Grouped result label/order assertion failed")
if (!identical(results$n, expected_n)) stop("Grouped model N assertion failed: ", paste(results$n, collapse = ","))
if (anyNA(results) || any(!is.finite(results$OR)) || any(!is.finite(results$P))) stop("Grouped output completeness assertion failed")

write.csv(results, result_csv, row.names = FALSE)
write.csv(model_qc, model_qc_csv, row.names = FALSE)

cat("GROUPED_STEP2_STATUS=PASS\n")
cat("GROUPED_MODELS_COMPLETED=", nrow(results), "\n", sep = "")
cat("GROUPED_P_LT_0_05=", sum(results$P < 0.05), "\n", sep = "")
