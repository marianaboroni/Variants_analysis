#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
if (length(args) >= 2 && args[[1]] == "--config") {
  config_path <- args[[2]]
}
if (is.null(config_path)) {
  stop("Usage: Rscript scripts/run_automated_variant_analysis.R --config config/example_config.yml")
}

source("R/io.R")
source("R/metadata.R")
source("R/features.R")
source("R/hard_filters.R")
source("R/scoring.R")
source("R/validation.R")
source("R/reference.R")
source("R/driver_classification.R")
source("R/ml_filtering.R")
source("R/clonality.R")
source("R/tmb.R")
source("R/qc_metrics.R")
source("R/reporting.R")
source("R/visualization.R")
source("R/workflow.R")

run_variant_analysis(config_path)
