options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript 06_common_cycle_sensitivity.R <analytic_rds> <model_csv> <fdr_csv>")

suppressPackageStartupMessages(library(survey))
input_rds <- args[1]
model_csv <- args[2]
fdr_csv <- args[3]
if (!file.exists(input_rds)) stop("Analytic input not found: ", input_rds)
dir.create(dirname(model_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(fdr_csv), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_rds)
exposures <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
covariates <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)
target_cycles <- c(
  "2005_2006", "2007_2008", "2009_2010", "2011_2012",
  "2013_2014", "2015_2016", "2017_2018"
)
required <- unique(c(
  "stroke", "sex", "cycle", "SDMVPSU", "SDMVSTRA", "WTMEC_comb",
  "WTMEC2YR", exposures, covariates
))
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

target <- dat$cycle %in% target_cycles
anchors <- c(
  n = sum(target), stroke = sum(target & dat$stroke == 1),
  female = sum(target & dat$sex == "Female"),
  female_stroke = sum(target & dat$sex == "Female" & dat$stroke == 1),
  male = sum(target & dat$sex == "Male"),
  male_stroke = sum(target & dat$sex == "Male" & dat$stroke == 1)
)
expected <- c(n = 8311, stroke = 339, female = 4179, female_stroke = 175, male = 4132, male_stroke = 164)
if (!identical(as.numeric(anchors), as.numeric(expected))) stop("Common-cycle sample anchors failed")
if (!all(complete.cases(dat[target, c("stroke", exposures, covariates, "SDMVPSU", "SDMVSTRA", "WTMEC_comb")]))) {
  stop("Unexpected missingness in common-cycle model variables")
}

dat$WTMEC_2005_2018 <- dat$WTMEC_comb * 8 / 7
if (max(abs(dat$WTMEC_2005_2018 - dat$WTMEC2YR / 7), na.rm = TRUE) >= 1e-12) {
  stop("Seven-cycle weight algebra check failed")
}
options(survey.lonely.psu = "fail")
design_all <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_2005_2018,
  data = dat, nest = TRUE
)
design_target <- subset(design_all, cycle %in% target_cycles)
designs <- list(
  Female = subset(design_target, sex == "Female"),
  Male = subset(design_target, sex == "Male")
)

run_one <- function(design, sex_value, exposure, order_index) {
  f <- as.formula(paste("stroke ~", exposure, "+", paste(covariates, collapse = " + ")))
  fit <- svyglm(f, design = design, family = quasibinomial(link = "logit"))
  beta <- unname(coef(fit)[exposure])
  se <- sqrt(unname(vcov(fit)[exposure, exposure]))
  model_frame <- model.frame(fit)
  data.frame(
    analysis_id = sprintf("OPT02-%s-%02d-%s", toupper(substr(sex_value, 1, 1)), match(exposure, exposures), exposure),
    sex = sex_value,
    exposure = exposure,
    cycle_window = "2005-2018",
    sample_n = nrow(model_frame),
    stroke_events = sum(model_frame$stroke == 1),
    unweighted_n = nrow(model_frame),
    weighted_population_total = sum(weights(design, type = "sampling")),
    beta = beta,
    standard_error = se,
    OR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    P_value = 2 * pnorm(-abs(beta / se)),
    survey_df = degf(design),
    convergence = isTRUE(fit$converged),
    aliased = is.na(beta) || is.na(se),
    original_row_order = order_index,
    stringsAsFactors = FALSE
  )
}

rows <- list()
index <- 0L
for (sex_value in c("Female", "Male")) {
  for (exposure in exposures) {
    index <- index + 1L
    rows[[index]] <- run_one(designs[[sex_value]], sex_value, exposure, index)
  }
}
result <- do.call(rbind, rows)
if (nrow(result) != 20L || !all(result$convergence) || any(result$aliased)) stop("Common-cycle model assertion failed")
result$OPT02_q_value <- p.adjust(result$P_value, method = "BH")

p <- result$P_value
m <- length(p)
ord <- order(p, seq_along(p))
rank_values <- integer(m)
rank_values[ord] <- seq_len(m)
raw_sorted <- p[ord] * m / seq_len(m)
manual_sorted <- pmin(1, rev(cummin(rev(raw_sorted))))
manual_q <- numeric(m)
manual_q[ord] <- manual_sorted
if (max(abs(manual_q - result$OPT02_q_value)) > 1e-12) stop("Manual and R BH values differ")

fdr <- data.frame(
  view = "original_order",
  original_row_order = result$original_row_order,
  analysis_id = result$analysis_id,
  sex = result$sex,
  exposure = result$exposure,
  P_value = p,
  rank = rank_values,
  m = m,
  raw_BH_value = p * m / rank_values,
  monotonic_adjusted_q = manual_q,
  R_p_adjust_q = result$OPT02_q_value,
  q_difference = manual_q - result$OPT02_q_value,
  stringsAsFactors = FALSE
)
fdr <- rbind(fdr, transform(fdr[order(fdr$rank), ], view = "sorted_by_P"))

model_out <- result[, c(
  "analysis_id", "sex", "exposure", "cycle_window", "sample_n", "stroke_events",
  "unweighted_n", "weighted_population_total", "beta", "standard_error", "OR",
  "CI_lower", "CI_upper", "P_value", "OPT02_q_value", "survey_df",
  "convergence", "aliased"
)]
write.csv(model_out, model_csv, row.names = FALSE)
write.csv(fdr, fdr_csv, row.names = FALSE)
cat("COMMON_CYCLE_STATUS=PASS\n")

