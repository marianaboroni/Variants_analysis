#!/usr/bin/env Rscript
# DEPRECATED. Kept only for backward compatibility. Use the single CLI instead:
#   inst/exec/tumoronly run --config <config.yml>
# This wrapper contains no logic; it delegates to run_tumor_only().

warning("scripts/run_tumor_only_filter.R is deprecated. Use: tumoronly run --config <config.yml>",
        call. = FALSE, immediate. = TRUE)

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[[1]] == "--config") args[[2]] else NULL
if (is.null(config_path)) stop("Usage: Rscript scripts/run_tumor_only_filter.R --config config/example.yml")

r_dir <- file.path(getwd(), "R")
for (f in list.files(r_dir, pattern = "[.]R$", full.names = TRUE)) sys.source(f, envir = globalenv())
run_tumor_only(config_path)
