# Stage 5-B release candidate: acquire official NHANES XPT files and build translated workspace RDS files.
# This script contains no project-specific absolute path. Paths are supplied explicitly at runtime.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 01_download_or_import_nhanes.R",
    "<release_candidate_root> <raw_source_cache> <translated_cache_or_dash> <reproduction_workspace>"
  ))
}

package_root <- gsub("\\\\", "/", args[[1]])
source_cache <- gsub("\\\\", "/", args[[2]])
translated_cache <- if (identical(args[[3]], "-")) NA_character_ else gsub("\\\\", "/", args[[3]])
workspace_root <- gsub("\\\\", "/", args[[4]])
manifest_path <- file.path(package_root, "config", "nhanes-components.csv")
translated_manifest_path <- file.path(package_root, "config", "translated-cache-manifest.csv")
workspace_raw <- file.path(workspace_root, "data_raw")

required_packages <- c("foreign", "digest")
if (is.na(translated_cache)) required_packages <- c(required_packages, "nhanesA")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
if (!file.exists(manifest_path)) {
  stop("Component manifest not found: ", manifest_path)
}
if (!is.na(translated_cache) && !file.exists(translated_manifest_path)) {
  stop("Translated-cache manifest not found: ", translated_manifest_path)
}

dir.create(source_cache, recursive = TRUE, showWarnings = FALSE)
dir.create(workspace_raw, recursive = TRUE, showWarnings = FALSE)

components <- read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE,
  na.strings = character()
)
required_columns <- c(
  "Cycle", "Domain", "Component", "Official_source",
  "Local_cache_filename", "Source_SHA", "Role"
)
if (!identical(names(components), required_columns)) {
  stop("Unexpected component-manifest schema")
}
if (nrow(components) != 66L || anyDuplicated(components$Component)) {
  stop("Expected 66 unique NHANES components")
}
translated_manifest <- if (!is.na(translated_cache)) {
  x <- read.csv(translated_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  expected_names <- c("Cycle", "Component", "Workspace_RDS", "Translated_RDS_SHA", "Rows", "Columns")
  if (!identical(names(x), expected_names) || nrow(x) != 66L || anyDuplicated(x$Component)) {
    stop("Unexpected translated-cache manifest")
  }
  x
} else NULL

audit_rows <- vector("list", nrow(components))
for (i in seq_len(nrow(components))) {
  row <- components[i, , drop = FALSE]
  cache_file <- do.call(
    file.path,
    as.list(c(source_cache, strsplit(row$Local_cache_filename, "/", fixed = TRUE)[[1]]))
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(cache_file)) {
    message(sprintf("[%02d/%02d] Downloading %s", i, nrow(components), row$Component))
    download.file(row$Official_source, cache_file, mode = "wb", quiet = TRUE)
  } else {
    message(sprintf("[%02d/%02d] Using cached %s", i, nrow(components), row$Component))
  }
  if (!file.exists(cache_file) || file.info(cache_file)$size <= 0) {
    stop("Missing or empty source file: ", row$Component)
  }

  source_sha <- toupper(digest::digest(cache_file, algo = "sha256", file = TRUE))
  expected_sha <- toupper(row$Source_SHA)
  if (!nzchar(expected_sha) || !identical(source_sha, expected_sha)) {
    stop(
      "Source identity mismatch for ", row$Component,
      "; expected ", expected_sha, "; observed ", source_sha
    )
  }

  raw <- foreign::read.xport(cache_file)
  if (is.null(raw) || !"SEQN" %in% names(raw)) {
    stop("Invalid XPT structure for ", row$Component)
  }
  if (!is.na(translated_cache)) {
    translated_row <- translated_manifest[translated_manifest$Component == row$Component, , drop = FALSE]
    if (nrow(translated_row) != 1L || translated_row$Cycle != row$Cycle) {
      stop("Translated-cache manifest lookup failed for ", row$Component)
    }
    translated_file <- do.call(
      file.path,
      as.list(c(translated_cache, strsplit(translated_row$Workspace_RDS, "/", fixed = TRUE)[[1]]))
    )
    if (!file.exists(translated_file)) stop("Translated cache file missing: ", row$Component)
    translated_sha <- toupper(digest::digest(translated_file, algo = "sha256", file = TRUE))
    if (!identical(translated_sha, toupper(translated_row$Translated_RDS_SHA))) {
      stop("Translated cache identity mismatch for ", row$Component)
    }
    translated <- readRDS(translated_file)
    if (nrow(translated) != translated_row$Rows || ncol(translated) != translated_row$Columns) {
      stop("Translated cache dimension mismatch for ", row$Component)
    }
  } else {
    translated_sha <- NA_character_
    translated <- suppressWarnings(suppressMessages(
      nhanesA::nhanesTranslate(
        row$Component,
        colnames = names(raw)[names(raw) != "SEQN"],
        data = raw
      )
    ))
  }
  if (is.null(translated) || nrow(translated) != nrow(raw)) {
    stop("Translation failed for ", row$Component)
  }

  cycle_dir <- file.path(workspace_raw, row$Cycle)
  dir.create(cycle_dir, recursive = TRUE, showWarnings = FALSE)
  rds_file <- file.path(cycle_dir, paste0(row$Component, ".rds"))
  saveRDS(translated, rds_file, version = 3)

  audit_rows[[i]] <- data.frame(
    Cycle = row$Cycle,
    Component = row$Component,
    Source_SHA = source_sha,
    Translated_Cache_SHA = translated_sha,
    Source_bytes = unname(file.info(cache_file)$size),
    Rows = nrow(translated),
    Columns = ncol(translated),
    Workspace_RDS = file.path(row$Cycle, paste0(row$Component, ".rds")),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit_rows)
write.csv(
  audit,
  file.path(workspace_root, "source-import-audit.csv"),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
message("Official-source acquisition/import complete: ", nrow(audit), " components")
