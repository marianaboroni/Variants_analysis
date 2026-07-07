#!/usr/bin/env Rscript
# DEPRECATED and REMOVED as a distinct pipeline. This script was a byte-for-byte
# duplicate of run_tumor_only_filter.R. Use the single CLI:
#   inst/exec/tumoronly run --config <config.yml>

stop(paste0(
  "scripts/run_automated_variant_analysis.R has been removed (it duplicated ",
  "run_tumor_only_filter.R).\nUse: tumoronly run --config <config.yml>"), call. = FALSE)
