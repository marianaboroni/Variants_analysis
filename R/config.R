# Configuration reading, validation and default resolution.
# read_config() (thin yaml::read_yaml wrapper) and cfg_get() live in io.R.

VALID_BUILDS <- c("GRCh37", "GRCh38")

#' Validate a tumoronly configuration.
#'
#' Checks required fields and reports, for each problem: the missing/invalid
#' field path, the value received, the value expected, and a correction example.
#' Errors are collected and reported together.
#'
#' @param cfg a configuration list (from [read_config()]).
#' @param require_input whether `input.vcf`/`input.variants` must exist on disk.
#' @return `cfg` invisibly if valid; otherwise stops with an actionable message.
#' @keywords internal
validate_config <- function(cfg, require_input = TRUE) {
  errs <- character()
  add_err <- function(field, received, expected, example) {
    errs[[length(errs) + 1]] <<- sprintf(
      "  - %s\n      received: %s\n      expected: %s\n      example:  %s",
      field, received, expected, example)
  }

  # --- input --------------------------------------------------------------
  vcf <- cfg_get(cfg, c("input", "path"),
                 cfg_get(cfg, c("input", "vcf"), cfg_get(cfg, c("input", "variants"), NULL)))
  if (is.null(vcf) || is.na(vcf)) {
    add_err("input.path", "missing", "path to a VCF/MAF/TSV variant file",
            "input:\n    path: data/sample.vep.vcf.gz")
  } else if (require_input && !file.exists(vcf)) {
    add_err("input.path", vcf, "an existing file path", "input:\n    path: data/sample.vep.vcf.gz")
  }

  build <- cfg_get(cfg, c("input", "genome_build"), NULL)
  if (!is.null(build) && !is.na(build) && !(build %in% VALID_BUILDS)) {
    add_err("input.genome_build", build, paste(VALID_BUILDS, collapse = " or "),
            "input:\n    genome_build: GRCh38")
  }

  # --- cosmic (optional, but if present validate builds) ------------------
  cosmic <- cfg_get(cfg, "cosmic", NULL)
  if (!is.null(cosmic)) {
    sb <- cfg_get(cfg, c("cosmic", "source_build"), NULL)
    tb <- cfg_get(cfg, c("cosmic", "target_build"), build)
    if (!is.null(sb) && !is.na(sb) && !(sb %in% VALID_BUILDS)) {
      add_err("cosmic.source_build", sb, paste(VALID_BUILDS, collapse = " or "),
              "cosmic:\n    source_build: GRCh37")
    }
    if (!is.null(tb) && !is.na(tb) && !(tb %in% VALID_BUILDS)) {
      add_err("cosmic.target_build", tb, paste(VALID_BUILDS, collapse = " or "),
              "cosmic:\n    target_build: GRCh38")
    }
  }

  # --- output -------------------------------------------------------------
  outdir <- cfg_get(cfg, c("analysis", "output_dir"), cfg_get(cfg, c("output", "dir"), NULL))
  if (is.null(outdir) || is.na(outdir)) {
    add_err("analysis.output_dir", "missing", "a directory path for run outputs",
            "analysis:\n    output_dir: results")
  }

  if (length(errs) > 0) {
    stop(sprintf("Invalid configuration (%d problem%s):\n%s",
                 length(errs), if (length(errs) == 1) "" else "s",
                 paste(errs, collapse = "\n")), call. = FALSE)
  }
  invisible(cfg)
}

#' Resolve configuration defaults into a fully-specified list.
#' Unifies the legacy `output.dir` / new `analysis.output_dir` and
#' `input.variants` / `input.vcf` spellings.
#' @keywords internal
resolve_config <- function(cfg) {
  cfg$input$vcf <- cfg_get(cfg, c("input", "path"),
                           cfg_get(cfg, c("input", "vcf"), cfg_get(cfg, c("input", "variants"), NULL)))
  cfg$input$format <- cfg_get(cfg, c("input", "format"), "auto")
  if (is.null(cfg$analysis)) cfg$analysis <- list()
  cfg$analysis$output_dir <- cfg_get(cfg, c("analysis", "output_dir"),
                                     cfg_get(cfg, c("output", "dir"), "results"))
  cfg$analysis$run_id <- cfg_get(cfg, c("analysis", "run_id"), NULL)
  cfg$analysis$threads <- cfg_get(cfg, c("analysis", "threads"), 1L)
  cfg
}

#' Write the resolved config as YAML into the run directory.
#' @keywords internal
write_resolved_config <- function(cfg, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(cfg, path)
  invisible(path)
}

#' Deterministic run id from config + input checksum when none is supplied.
#' @keywords internal
resolve_run_id <- function(cfg) {
  rid <- cfg_get(cfg, c("analysis", "run_id"), NULL)
  if (!is.null(rid) && !is.na(rid) && nzchar(rid)) return(as.character(rid))
  vcf <- cfg$input$vcf
  seed <- paste(basename(vcf %||% "na"), file_sha256(vcf), sep = "|")
  paste0("run_", substr(digest::digest(seed, algo = "sha256"), 1, 12))
}
