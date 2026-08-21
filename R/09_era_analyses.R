options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript 09_era_analyses.R <analytic_rds> <stratified_csv> <interaction_csv>")

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")
input_rds <- args[1]
stratified_csv <- args[2]
interaction_csv <- args[3]
if (!file.exists(input_rds)) stop("Analytic input not found: ", input_rds)
dir.create(dirname(stratified_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(interaction_csv), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_rds)
early_cycles <- c("2003_2004", "2005_2006", "2007_2008", "2009_2010")
late_cycles <- c("2011_2012", "2013_2014", "2015_2016", "2017_2018")
dat$era <- ifelse(dat$cycle %in% early_cycles, 0, ifelse(dat$cycle %in% late_cycles, 1, NA))
early <- dat[dat$cycle %in% early_cycles, ]
late <- dat[dat$cycle %in% late_cycles, ]
early$WTMEC_period <- early$WTMEC2YR / 4
late$WTMEC_period <- late$WTMEC2YR / 4

exposures <- c(
  "MBP_ln", "MBzP_ln", "MECPP_ln", "MEHHP_ln", "MEOHP_ln",
  "MCPP_ln", "MEP_ln", "MiBP_ln", "MCNP_ln", "MCOP_ln"
)
labels <- sub("_ln$", "", exposures)
covariates <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)
required <- unique(c(
  "stroke", "sex", "cycle", "SDMVPSU", "SDMVSTRA", "WTMEC_comb", "WTMEC2YR",
  exposures, covariates
))
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

run_stratified <- function(data, exposure) {
  keep_covariates <- covariates[vapply(covariates, function(variable) {
    values <- data[[variable]]
    if (is.factor(values) || is.character(values)) length(unique(na.omit(values))) >= 2L else sd(values, na.rm = TRUE) > 0
  }, logical(1))]
  design <- svydesign(
    ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_period,
    data = data, nest = TRUE
  )
  f <- as.formula(paste("stroke ~", exposure, "+", paste(keep_covariates, collapse = " + ")))
  fit <- svyglm(f, design = design, family = quasibinomial(link = "logit"))
  coefficient_table <- summary(fit)$coefficients
  beta <- unname(coefficient_table[exposure, "Estimate"])
  se <- unname(coefficient_table[exposure, "Std. Error"])
  p <- unname(coefficient_table[exposure, "Pr(>|t|)"])
  c(OR = round(exp(beta), 3), lower = round(exp(beta - 1.96 * se), 3), upper = round(exp(beta + 1.96 * se), 3), P = round(p, 4))
}

stratified_rows <- list()
index <- 0L
for (sex_value in c("Female", "Male")) {
  early_sex <- early[early$sex == sex_value, ]
  late_sex <- late[late$sex == sex_value, ]
  for (j in seq_along(exposures)) {
    index <- index + 1L
    early_result <- run_stratified(early_sex, exposures[j])
    late_result <- run_stratified(late_sex, exposures[j])
    stratified_rows[[index]] <- data.frame(
      Metabolite = labels[j], Sex = sex_value,
      Early_OR = early_result["OR"],
      Early_CI = sprintf("%.3f (%.3f-%.3f)", early_result["OR"], early_result["lower"], early_result["upper"]),
      Early_P = early_result["P"],
      Late_OR = late_result["OR"],
      Late_CI = sprintf("%.3f (%.3f-%.3f)", late_result["OR"], late_result["lower"], late_result["upper"]),
      Late_P = late_result["P"],
      stringsAsFactors = FALSE
    )
  }
}
stratified <- do.call(rbind, stratified_rows)
if (nrow(stratified) != 20L) stop("Era-stratified row-count assertion failed")
write.csv(stratified, stratified_csv, row.names = FALSE)

male <- dat[dat$sex == "Male", ]
male_design <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_comb,
  data = male, nest = TRUE
)
run_interaction <- function(exposure) {
  f <- as.formula(paste("stroke ~", exposure, "* era +", paste(covariates, collapse = " + ")))
  fit <- svyglm(f, design = male_design, family = quasibinomial())
  term <- paste0(exposure, ":era")
  coefficient <- summary(fit)$coefficients[term, ]
  complete_variables <- c("stroke", "era", exposure, covariates, "SDMVPSU", "SDMVSTRA", "WTMEC_comb")
  complete_data <- male[complete.cases(male[, complete_variables]), ]
  data.frame(
    exposure = sub("_ln$", "", exposure), term = term,
    beta_interaction = unname(coefficient["Estimate"]),
    SE_interaction = unname(coefficient["Std. Error"]),
    OR_interaction = exp(unname(coefficient["Estimate"])),
    CI_low = exp(unname(coefficient["Estimate"]) - 1.96 * unname(coefficient["Std. Error"])),
    CI_high = exp(unname(coefficient["Estimate"]) + 1.96 * unname(coefficient["Std. Error"])),
    P_interaction = unname(coefficient["Pr(>|t|)"]),
    n_complete = nrow(complete_data),
    stroke_early = sum(complete_data$stroke[complete_data$era == 0]),
    stroke_late = sum(complete_data$stroke[complete_data$era == 1]),
    male_stroke_early_all = sum(male$stroke[male$era == 0], na.rm = TRUE),
    male_stroke_late_all = sum(male$stroke[male$era == 1], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
interaction <- rbind(run_interaction("MCNP_ln"), run_interaction("MCOP_ln"))
if (nrow(interaction) != 2L || !all(interaction$n_complete == 4132L)) stop("Era-interaction assertion failed")
write.csv(interaction, interaction_csv, row.names = FALSE)
cat("ERA_ANALYSES_STATUS=PASS\n")

