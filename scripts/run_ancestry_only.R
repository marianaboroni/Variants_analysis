#!/usr/bin/env Rscript
# Run ONLY the (experimental) genetic-ancestry module on one or more VCFs,
# without the somatic pipeline (filters, COSMIC, ABraOM, report). Used to
# compare ancestry inputs side by side, e.g.:
#   - caller VCFs with variants only (Mutect2, sarek `bcftools call -mv`),
#     where homozygous-reference AIMs are absent (biased; see docs/ANCESTRY.md)
#   - per-sample AIM genotypes from scripts/genotype_aims.sh (0/0 included)
#
# Each VCF is read, standardized and featurized exactly as in `tumoronly run`
# (read_variant_input -> standardize_variant_table -> attach_sample_metadata
# -> add_basic_features) and passed to infer_ancestry(), one file at a time
# to bound memory. Unlike `tumoronly run`, hotspot/driver annotation is not
# computed, so those loci are not excluded (irrelevant for AIM-only input).
#
# Usage:
#   Rscript scripts/run_ancestry_only.R --config <config.yml> --output <dir> <a.vcf.gz> [b.vcf.gz ...]
#
# Uses the config's `ancestry:` section (marker_panel, qc, thresholds) and
# its input/technical settings; `ancestry.enabled` is not required.
# Writes ancestry_summary.tsv, ancestry_qc.tsv, ancestry_projection.tsv,
# ancestry_snps.tsv.gz and plots/ancestry/ under <dir>.

script_root <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  dirname(dirname(normalizePath(file_arg[1])))
}

# Same dev loader as exec/tumoronly: source R/ (utils/io/config/input first).
root <- script_root()
r_dir <- file.path(root, "R")
files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
first <- file.path(r_dir, c("utils.R", "io.R", "config.R", "input.R"))
for (f in c(intersect(first, files), setdiff(files, first))) sys.source(f, envir = globalenv())
Sys.setenv(TUMORONLY_R_DIR = r_dir)

argv <- commandArgs(trailingOnly = TRUE)
get_opt <- function(name) {
  i <- which(argv == name)
  if (length(i) && length(argv) >= i[1] + 1) argv[i[1] + 1] else NULL
}
config <- get_opt("--config"); out_dir <- get_opt("--output")
opt_idx <- which(argv %in% c("--config", "--output"))
vcfs <- argv[-c(opt_idx, opt_idx + 1)]
if (is.null(config) || is.null(out_dir) || length(vcfs) == 0) {
  cat("Usage: Rscript scripts/run_ancestry_only.R --config <config.yml> --output <dir> <a.vcf.gz> [b.vcf.gz ...]\n")
  quit(status = 1)
}
missing <- vcfs[!file.exists(vcfs)]
if (length(missing)) stop("VCF(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)

cfg <- resolve_config(read_config(config))
cfg$input$vcf <- vcfs
build <- resolve_genome_build(cfg, vcfs)
if (is.null(load_aims_panel(cfg, build)))
  stop("ancestry.marker_panel.path is missing or unreadable in ", config, call. = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

parts <- list()
for (f in vcfs) {
  log_step("ancestry-only", "reading", file = f)
  ing <- read_variant_input(f, format = cfg_get(cfg, c("input", "format"), "auto"),
                            sample_metadata = cfg_get(cfg, c("input", "sample_metadata"), NULL),
                            chunk_size = cfg_get(cfg, c("input", "chunk_size"), NULL))
  v <- standardize_variant_table(ing$variants, cfg)
  rm(ing); gc()
  v <- attach_sample_metadata(v, cfg)
  v <- add_basic_features(v, cfg)
  log_step("ancestry-only", "inferring", file = basename(f), n = nrow(v))
  parts[[f]] <- infer_ancestry(v, cfg, build)
  rm(v); gc()
}

bind <- function(k) {
  d <- Filter(function(z) !is.null(z) && nrow(z) > 0, lapply(parts, `[[`, k))
  if (length(d)) data.table::rbindlist(d, fill = TRUE) |> as.data.frame() else data.frame()
}
anc <- list(summary = bind("summary"), snps = bind("snps"), qc = bind("qc"),
            projection = bind("projection"), panel = parts[[1]]$panel)
rownames(anc$summary) <- NULL

write_tsv(anc$summary, file.path(out_dir, "ancestry_summary.tsv"))
write_tsv(anc$qc, file.path(out_dir, "ancestry_qc.tsv"))
write_tsv(anc$projection, file.path(out_dir, "ancestry_projection.tsv"))
if (nrow(anc$snps) > 0) write_tsv(anc$snps, file.path(out_dir, "ancestry_snps.tsv.gz"))
tryCatch(create_ancestry_plots(anc, file.path(out_dir, "plots"), cfg),
         error = function(e) log_step("ancestry-only", "plots failed", error = conditionMessage(e)))
log_step("ancestry-only", "done", output = out_dir, samples = nrow(anc$summary))
