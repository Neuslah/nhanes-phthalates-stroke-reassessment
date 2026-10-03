options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 10_male_subgroup_global_interactions.R <analytic_rds> <output_csv>")

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
input_rds <- args[1]
output_csv <- args[2]
if (!file.exists(input_rds)) stop("Analytic input not found: ", input_rds)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_rds)
required <- c(
  "stroke", "sex", "RIDAGEYR", "race_eth_f", "edu_f", "marital_f",
  "pir", "bmi", "smoking_f", "alcohol_f", "hypertension_f",
  "diabetes_f", "hyperlipidemia_f", "SDMVPSU", "SDMVSTRA",
  "WTMEC_comb", "MCNP_ln", "MCOP_ln"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

dat$age_group_f <- factor(ifelse(dat$RIDAGEYR < 60, "<60", ">=60"), levels = c("<60", ">=60"))
dat$race_subgrp <- factor(
  ifelse(dat$race_eth_f == "Non-Hispanic White", "NH White",
  ifelse(dat$race_eth_f == "Non-Hispanic Black", "NH Black",
  ifelse(dat$race_eth_f %in% c("Mexican American", "Other Hispanic"), "Hispanic", "Other"))),
  levels = c("NH White", "NH Black", "Hispanic", "Other")
)
dat$bmi_group_f <- factor(
  ifelse(dat$bmi < 25, "<25", ifelse(dat$bmi < 30, "25-30", ifelse(!is.na(dat$bmi), ">=30", NA))),
  levels = c("<25", "25-30", ">=30")
)
male <- dat[dat$sex == "Male" & !is.na(dat$stroke) & dat$RIDAGEYR >= 20, ]
base_covariates <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)
exposures <- c("MCNP_ln", "MCOP_ln")
modifier_specs <- list(
  list(Modifier = "Age", variable = "age_group_f", covariates = c(base_covariates, "age_group_f"), expected_df = 1L),
  list(Modifier = "BMI", variable = "bmi_group_f", covariates = c(setdiff(base_covariates, "bmi"), "bmi_group_f"), expected_df = 2L),
  list(Modifier = "Race/ethnicity", variable = "race_subgrp", covariates = c(setdiff(base_covariates, "race_eth_f"), "race_subgrp"), expected_df = 3L)
)
design <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_comb,
  data = male, nest = TRUE
)
formula_text <- function(formula) paste(deparse(formula, width.cutoff = 500L), collapse = " ")

rows <- list()
index <- 0L
for (exposure in exposures) {
  for (spec in modifier_specs) {
    index <- index + 1L
    modifier <- spec$variable
    covariates_no_modifier <- setdiff(spec$covariates, modifier)
    f <- as.formula(paste(
      "stroke ~", paste(c(paste0(exposure, " * ", modifier), covariates_no_modifier), collapse = " + ")
    ))
    model_variables <- unique(all.vars(f))
    complete_index <- complete.cases(male[, model_variables, drop = FALSE])
    model_n <- sum(complete_index)
    model_events <- sum(male$stroke[complete_index] == 1, na.rm = TRUE)
    fit <- svyglm(f, design = design, family = quasibinomial(link = "logit"))
    coefficient_names <- rownames(summary(fit)$coefficients)
    interaction_names <- coefficient_names[
      grepl(":", coefficient_names, fixed = TRUE) &
      grepl(exposure, coefficient_names, fixed = TRUE) &
      grepl(modifier, coefficient_names, fixed = TRUE)
    ]
    if (length(interaction_names) != spec$expected_df) stop("Unexpected interaction coefficient count for ", exposure, " / ", spec$Modifier)
    term_test <- survey::regTermTest(fit, as.formula(paste("~", exposure, ":", modifier)), method = "Wald")
    numerator_df <- as.numeric(term_test$df)
    denominator_df <- as.numeric(term_test$ddf)
    test_statistic <- if (!is.null(term_test$Ftest)) as.numeric(term_test$Ftest) else as.numeric(term_test$chisq)
    global_p <- as.numeric(term_test$p)
    if (numerator_df != spec$expected_df || is.na(global_p)) stop("Incomplete global Wald test for ", exposure, " / ", spec$Modifier)
    rows[[index]] <- data.frame(
      Exposure = exposure,
      Modifier = spec$Modifier,
      Test_method = "survey::regTermTest (Wald)",
      Numerator_df = numerator_df,
      Denominator_df = denominator_df,
      Test_statistic = test_statistic,
      Global_P_interaction = global_p,
      model_formula = formula_text(f),
      model_n = model_n,
      stroke_events = model_events,
      interaction_coefficient_count = length(interaction_names),
      expected_interaction_df = spec$expected_df,
      status = "COMPLETED",
      warnings = "",
      stringsAsFactors = FALSE
    )
  }
}
result <- do.call(rbind, rows)
if (nrow(result) != 6L || !all(result$model_n == 4132L) || !all(result$stroke_events == 164L)) {
  stop("Male subgroup global-interaction assertion failed")
}
write.csv(result, output_csv, row.names = FALSE)
cat("MALE_SUBGROUP_GLOBAL_INTERACTION_STATUS=PASS\n")

