options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 08a_grouped_derive_current_weights.R <analytic_rds> <grouped_rds> <weight_qc_csv>")
}

input_rds <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
candidate_rds <- args[2]
weight_qc_csv <- args[3]
dir.create(dirname(candidate_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(weight_qc_csv), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_rds)
if (!is.data.frame(dat)) stop("Current V5 input is not a data.frame")
if (nrow(dat) != 9129L) stop("Current V5 row-count assertion failed: ", nrow(dat))

required <- c(
  "SEQN", "cycle", "WTMEC2YR", "WTMEC_comb",
  "MEP_ln", "MBP_ln", "MiBP_ln", "MBzP_ln", "MEHHP_ln",
  "MEOHP_ln", "MECPP_ln", "MCPP_ln", "MCNP_ln", "MCOP_ln"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) {
  stop("Missing required variables: ", paste(missing_required, collapse = ", "))
}

expected_new <- c(
  "MEP_cr", "MBP_cr", "MiBP_cr", "MBzP_cr", "MEHHP_cr",
  "MEOHP_cr", "MECPP_cr", "MCPP_cr", "MCNP_cr", "MCOP_cr",
  "sigma_LMW", "sigma_HMW", "sigma_DEHP",
  "ln_sigma_LMW", "ln_sigma_HMW", "ln_sigma_DEHP"
)
already_present <- intersect(expected_new, names(dat))
if (length(already_present)) {
  stop("Current V5 already contains grouped-derived variables: ", paste(already_present, collapse = ", "))
}

original_names <- names(dat)
original_seqn <- dat$SEQN
original_wtmec2yr <- dat$WTMEC2YR
original_wtmec_comb <- dat$WTMEC_comb

if (length(unique(dat$cycle)) != 8L) stop("Expected eight survey cycles")
if (anyNA(original_wtmec2yr) || any(original_wtmec2yr <= 0)) stop("WTMEC2YR contains missing/non-positive values")
if (anyNA(original_wtmec_comb) || any(original_wtmec_comb <= 0)) stop("WTMEC_comb contains missing/non-positive values")
if (!identical(original_wtmec_comb, original_wtmec2yr / 8)) {
  stop("Current V5 WTMEC_comb is not exactly WTMEC2YR / 8")
}

MW <- list(
  MEP = 194.18, MBP = 222.24, MiBP = 222.24, MBzP = 312.36,
  MEHHP = 294.35, MEOHP = 292.33, MECPP = 308.33,
  MCPP = 252.22, MCNP = 336.38, MCOP = 322.35
)

for (metabolite in names(MW)) {
  dat[[paste0(metabolite, "_cr")]] <- exp(dat[[paste0(metabolite, "_ln")]])
}

dat$sigma_LMW <- dat$MEP_cr / MW$MEP +
  dat$MBP_cr / MW$MBP +
  dat$MiBP_cr / MW$MiBP

dat$sigma_HMW <- dat$MBzP_cr / MW$MBzP +
  dat$MEHHP_cr / MW$MEHHP +
  dat$MEOHP_cr / MW$MEOHP +
  dat$MECPP_cr / MW$MECPP +
  dat$MCNP_cr / MW$MCNP +
  dat$MCOP_cr / MW$MCOP +
  dat$MCPP_cr / MW$MCPP

dat$sigma_DEHP <- dat$MEHHP_cr / MW$MEHHP +
  dat$MEOHP_cr / MW$MEOHP +
  dat$MECPP_cr / MW$MECPP

dat$ln_sigma_LMW <- log(dat$sigma_LMW)
dat$ln_sigma_HMW <- log(dat$sigma_HMW)
dat$ln_sigma_DEHP <- log(dat$sigma_DEHP)

actual_new <- setdiff(names(dat), original_names)
if (!identical(actual_new, expected_new)) {
  stop(
    "Derived-variable assertion failed. Expected: ", paste(expected_new, collapse = ", "),
    "; actual: ", paste(actual_new, collapse = ", ")
  )
}
if (!identical(dat$SEQN, original_seqn)) stop("SEQN order changed during Step1")
if (!identical(dat$WTMEC2YR, original_wtmec2yr)) stop("Step1 changed WTMEC2YR")
if (!identical(dat$WTMEC_comb, original_wtmec_comb)) stop("Step1 changed WTMEC_comb")

for (variable in c("sigma_LMW", "sigma_HMW", "sigma_DEHP")) {
  observed <- dat[[variable]]
  if (!all(observed[!is.na(observed)] > 0)) stop(variable, " contains non-positive values")
}

observed_hmw_cycles <- sort(unique(as.character(dat$cycle[!is.na(dat$ln_sigma_HMW)])))
expected_hmw_cycles <- sort(setdiff(unique(as.character(dat$cycle)), "2003_2004"))
if (!identical(observed_hmw_cycles, expected_hmw_cycles)) {
  stop("HMW cycle-availability assertion failed")
}

saveRDS(dat, candidate_rds)
candidate <- readRDS(candidate_rds)
if (!identical(candidate$SEQN, original_seqn)) stop("Persisted candidate SEQN order changed")
if (!identical(candidate$WTMEC2YR, original_wtmec2yr)) stop("Persisted candidate WTMEC2YR changed")
if (!identical(candidate$WTMEC_comb, original_wtmec_comb)) stop("Persisted candidate WTMEC_comb changed")
if (!identical(setdiff(names(candidate), original_names), expected_new)) stop("Persisted candidate derived-variable assertion failed")

qc <- data.frame(
  Check = c(
    "Current_V5_rows_9129",
    "Eight_cycles",
    "WTMEC_comb_equals_WTMEC2YR_div_8",
    "WTMEC2YR_preserved_in_memory",
    "WTMEC_comb_preserved_in_memory",
    "Exactly_16_grouped_variables_added",
    "HMW_available_only_2005_2018",
    "WTMEC2YR_preserved_after_save",
    "WTMEC_comb_preserved_after_save"
  ),
  Passed = TRUE,
  Detail = c(
    as.character(nrow(dat)),
    paste(sort(unique(as.character(dat$cycle))), collapse = " | "),
    "Exact identical() comparison",
    "Exact identical() comparison",
    "Exact identical() comparison",
    paste(expected_new, collapse = " | "),
    paste(observed_hmw_cycles, collapse = " | "),
    "Exact identical() comparison after read-back",
    "Exact identical() comparison after read-back"
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, weight_qc_csv, row.names = FALSE)

cat("GROUPED_STEP1_STATUS=PASS\n")
cat("INPUT_ROWS=", nrow(dat), "\n", sep = "")
cat("DERIVED_VARIABLES_ADDED=", length(expected_new), "\n", sep = "")
cat("WEIGHT_COLUMNS_CHANGED=0\n")
