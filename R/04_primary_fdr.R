# Stage 5-B release candidate: primary 20-test Benjamini-Hochberg correction.
# Calculation and purely mechanical canonical-schema normalization are consolidated here.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 04_primary_fdr.R <female_primary_csv> <male_primary_csv> <output_csv>")
}
female_file <- gsub("\\\\", "/", args[[1]])
male_file <- gsub("\\\\", "/", args[[2]])
output_file <- gsub("\\\\", "/", args[[3]])
if (!requireNamespace("digest", quietly = TRUE)) stop("Package digest is required")

read_m3 <- function(path, sex_label) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  keep <- which(x$Model == "Model3")
  if (length(keep) != 10L) stop("Expected 10 Model 3 rows in ", basename(path))
  data.frame(
    sex = sex_label,
    metabolite = sub("_ln$", "", x$Metabolite[keep]),
    model = "Model3",
    p_exact = as.numeric(x[["P value"]][keep]),
    input_source = basename(path),
    input_row = keep + 1L,
    input_source_sha256 = toupper(digest::digest(path, algo = "sha256", file = TRUE)),
    stringsAsFactors = FALSE
  )
}

d <- rbind(read_m3(female_file, "Female"), read_m3(male_file, "Male"))
if (nrow(d) != 20L || any(!is.finite(d$p_exact))) {
  stop("Primary FDR family must contain 20 finite Model 3 P values")
}
d$q_BH <- p.adjust(d$p_exact, method = "BH")
d$q_BH_display <- sprintf("%.3f", d$q_BH)
d$calculation_method <- "p.adjust across 20 Model 3 P values, method=BH"
d$survived_FDR_0_05 <- d$q_BH < 0.05
d <- d[, c(
  "sex", "metabolite", "model", "p_exact", "q_BH", "q_BH_display",
  "input_source", "input_row", "input_source_sha256",
  "calculation_method", "survived_FDR_0_05"
)]

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(d, output_file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
cat("FDR rows:", nrow(d), "\n")
cat("q < 0.05:", sum(d$survived_FDR_0_05), "\n")
