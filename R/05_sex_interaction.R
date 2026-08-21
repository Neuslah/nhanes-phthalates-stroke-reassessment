options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 05_sex_interaction.R <analytic_rds> <output_csv>")

suppressPackageStartupMessages(library(survey))
options(survey.lonely.psu = "adjust")

input_rds <- args[1]
output_csv <- args[2]
if (!file.exists(input_rds)) stop("Analytic input not found: ", input_rds)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_rds)
base_covars <- c(
  "RIDAGEYR", "race_eth_f", "edu_f", "marital_f", "pir", "bmi",
  "smoking_f", "alcohol_f", "hypertension_f", "diabetes_f",
  "hyperlipidemia_f"
)
exposures <- c("MBzP_ln", "MCNP_ln", "MCOP_ln")
required <- unique(c(
  "stroke", "sex", "SDMVPSU", "SDMVSTRA", "WTMEC_comb",
  base_covars, exposures
))
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

design <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_comb,
  data = dat, nest = TRUE
)

run_one <- function(exposure) {
  f <- as.formula(paste(
    "stroke ~", exposure, "+ sex +", paste0(exposure, ":sex"), "+",
    paste(base_covars, collapse = " + ")
  ))
  fit <- svyglm(f, design = design, family = quasibinomial(link = "logit"))
  coefficient_table <- summary(fit)$coefficients
  interaction_row <- grep(paste0(exposure, ":"), rownames(coefficient_table))
  if (!length(interaction_row)) interaction_row <- grep(paste0(":", exposure), rownames(coefficient_table))
  if (length(interaction_row) != 1L) stop("Expected exactly one sex-interaction coefficient for ", exposure)
  data.frame(
    Metabolite = exposure,
    Subgroup = "Sex (all)",
    P_interaction = round(coefficient_table[interaction_row, "Pr(>|t|)"], 4),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, lapply(exposures, run_one))
if (nrow(result) != 3L || !identical(result$Metabolite, exposures)) stop("Sex-interaction output assertion failed")
write.csv(result, output_csv, row.names = FALSE)
cat("SEX_INTERACTION_STATUS=PASS\n")

