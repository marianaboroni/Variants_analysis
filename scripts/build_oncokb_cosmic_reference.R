#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- NULL
if (length(args) >= 2 && args[[1]] == "--config") {
  config_path <- args[[2]]
}
if (is.null(config_path)) {
  stop("Usage: Rscript scripts/build_oncokb_cosmic_reference.R --config config/reference_build_example.yml")
}

source("R/io.R")
source("R/validation.R")
source("R/reference.R")
source("R/driver_classification.R")

cfg <- read_config(config_path)
outdir <- cfg$output$dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

oncokb_input <- cfg_get(cfg, c("reference_inputs", "oncokb"), NULL)
cosmic_input <- cfg_get(cfg, c("reference_inputs", "cosmic"), NULL)
driver_genes_input <- cfg_get(cfg, c("reference_inputs", "driver_genes"), NULL)
hotspots_input <- cfg_get(cfg, c("reference_inputs", "hotspots"), NULL)

if (!is.null(oncokb_input) && !is.na(oncokb_input)) {
  message("Standardizing OncoKB reference: ", oncokb_input)
  standardize_oncokb_reference(
    oncokb_input,
    file.path(outdir, "oncokb_reference.tsv"),
    cfg_get(cfg, c("reference_inputs", "oncokb_delimiter"), "\t")
  )
}

if (!is.null(cosmic_input) && !is.na(cosmic_input)) {
  message("Standardizing COSMIC reference: ", cosmic_input)
  standardize_cosmic_reference(
    cosmic_input,
    file.path(outdir, "cosmic_reference.tsv"),
    cfg_get(cfg, c("reference_inputs", "cosmic_delimiter"), "\t")
  )
}

if (!is.null(driver_genes_input) && !is.na(driver_genes_input)) {
  message("Standardizing driver gene reference: ", driver_genes_input)
  ref <- read_variants(driver_genes_input, cfg_get(cfg, c("reference_inputs", "driver_genes_delimiter"), "\t"))
  write_tsv(standardize_driver_gene_reference(ref), file.path(outdir, "driver_genes_reference.tsv"))
}

if (!is.null(hotspots_input) && !is.na(hotspots_input)) {
  message("Standardizing hotspot reference: ", hotspots_input)
  ref <- read_variants(hotspots_input, cfg_get(cfg, c("reference_inputs", "hotspots_delimiter"), "\t"))
  write_tsv(standardize_hotspot_reference(ref), file.path(outdir, "hotspots_reference.tsv"))
}

message("Done: ", normalizePath(outdir, mustWork = FALSE))
