options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run_tests.R <package_root> <workspace_root> <output_root>")
rscript <- file.path(R.home("bin"), "Rscript.exe")
run <- function(script, script_args) {
  status <- system2(rscript, c("--vanilla", shQuote(script), vapply(script_args, shQuote, character(1))))
  if (!identical(status, 0L)) stop("Test failed: ", basename(script))
}
run(file.path(args[[1]], "tests", "01_static_package_tests.R"), args[[1]])
run(file.path(args[[1]], "tests", "02_clean_run_tests.R"), args[2:3])
cat("PACKAGE_TEST_SUITE=PASS\n")
