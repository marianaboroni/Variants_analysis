# run_tumor_only(): the single tumor-only pipeline orchestrator.
#
# VEP-VCF/TSV -> build resolution -> standardize -> features -> technical filters
# -> scoring -> COSMIC (canonical key) -> driver/hotspot -> guideline -> classify
# -> decision trail -> standardized run directory (tables, MAF, plots, manifest).
# OncoKB is NOT part of this flow (post-hoc; see annotate_oncokb()).

#' Run the tumor-only variant analysis.
#'
#' @param config path to a YAML config, or an already-read config list.
#' @param run_id optional explicit run id (otherwise from config, else derived
#'   deterministically from the input checksum).
#' @return (invisibly) a list with the annotated `variants`, the resolved
#'   `run_dir`, the `manifest`, and the filter `audit`.
#' @export
#' @examples
#' \dontrun{
#' run_tumor_only("config/example.yml")
#' }
run_tumor_only <- function(config, run_id = NULL) {
  t0 <- Sys.time()
  cfg <- if (is.character(config)) read_config(config) else config
  validate_config(cfg, require_input = TRUE)
  cfg <- resolve_config(cfg)
  if (!is.null(run_id)) cfg$analysis$run_id <- run_id

  rid <- resolve_run_id(cfg)
  run_dir <- file.path(cfg$analysis$output_dir, rid)
  dirs <- setup_run_dir(run_dir)
  cfg$input$genome_build <- resolve_genome_build(cfg, cfg$input$vcf)
  build <- cfg$input$genome_build
  write_resolved_config(cfg, file.path(run_dir, "config.resolved.yml"))

  log_step("run", "reading variants", file = cfg$input$vcf, run_id = rid)
  variants <- read_variants(cfg$input$vcf, cfg_get(cfg, c("input", "delimiter"), "\t"))
  n_input <- nrow(variants)
  variants <- standardize_variant_table(variants, cfg)
  variants <- attach_sample_metadata(variants, cfg)

  log_step("run", "feature engineering + technical filters", n = nrow(variants))
  variants <- add_basic_features(variants, cfg)
  hf <- add_sample_qc_and_hard_filters(variants, cfg)
  variants <- hf$variants
  sample_qc <- hf$sample_qc
  recurrence <- compute_cohort_recurrence(variants, cfg)
  variants <- merge_recurrence_features(variants, recurrence)

  log_step("run", "scoring + COSMIC evidence + classification")
  variants <- score_variants(variants, cfg)
  variants <- annotate_cosmic(variants, cfg, build)
  variants <- add_driver_annotation_layers(variants, cfg)
  variants <- add_guideline_classifications(variants, cfg)
  variants <- classify_variants(variants, cfg)
  variants <- score_and_classify_driver_layers(variants, cfg)

  variants <- add_decision_trail(variants)

  # ---- outputs ----------------------------------------------------------
  audit <- write_run_tables(variants, dirs)
  maf_path <- file.path(dirs$maf, "filtered.maf.gz")
  retained <- variants[variants$filter_status %in% c("PASS", "REVIEW"), , drop = FALSE]
  maf <- NULL
  if (nrow(retained) > 0) {
    maf <- tryCatch(create_maf(retained, maf_path), error = function(e) {
      log_step("maf", "MAF creation failed", error = conditionMessage(e)); NULL })
    if (!is.null(maf)) make_oncoplot(maf_path, dirs$plots, cfg)
  } else {
    log_step("maf", "no retained variants; MAF and oncoplot skipped")
  }

  manifest <- build_run_manifest(cfg, run_dir, rid, build, variants, n_input, audit, t0)
  jsonlite::write_json(manifest, file.path(run_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")

  render_ok <- tryCatch({ render_tumor_only_report(run_dir); TRUE },
                        error = function(e) { log_step("report", "report skipped",
                          error = conditionMessage(e)); FALSE })

  log_step("run", "done", run_dir = run_dir,
           retained = nrow(retained), excluded = sum(variants$filter_status == "FAIL"),
           seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  invisible(list(variants = variants, run_dir = run_dir, manifest = manifest,
                 audit = audit, report_rendered = render_ok))
}

#' @keywords internal
setup_run_dir <- function(run_dir) {
  dirs <- list(
    root = run_dir,
    tables = file.path(run_dir, "tables"),
    maf = file.path(run_dir, "maf"),
    plots = file.path(run_dir, "plots"),
    report = file.path(run_dir, "report"),
    logs = file.path(run_dir, "logs"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

#' Write the standardized run tables; returns the filter audit.
#' @keywords internal
write_run_tables <- function(variants, dirs) {
  split <- split_by_status(variants)
  write_tsv(variants, file.path(dirs$tables, "variants_all.tsv.gz"))
  write_tsv(split$retained, file.path(dirs$tables, "variants_retained.tsv.gz"))
  write_tsv(split$excluded, file.path(dirs$tables, "variants_excluded.tsv.gz"))
  audit <- build_filter_audit(variants)
  write_tsv(audit, file.path(dirs$tables, "filter_audit.tsv.gz"))
  cosmic_cols <- intersect(c("sample_id", "variant_id", "chrom", "pos", "ref", "alt",
    "gene", "tumor_type", "tumor_subtype", "COSMIC_MATCH", "COSMIC_RELEASE", "COSMIC_MUTATION_IDS",
    "COSMIC_TOTAL_OCCURRENCES", "COSMIC_MATCHING_TUMOR_OCCURRENCES", "COSMIC_OTHER_TUMOR_OCCURRENCES",
    "COSMIC_MATCHING_TUMOR_TYPES", "COSMIC_OTHER_TUMOR_TYPES",
    "COSMIC_SITE_MATCH", "COSMIC_HISTOLOGY_MATCH", "COSMIC_SUBTYPE_MATCH", "COSMIC_MATCH_LEVEL",
    "COSMIC_TUMOR_CONTEXT_STATUS", "COSMIC_CONTEXT_EVALUABLE", "COSMIC_CONTEXT_INTERPRETATION",
    "COSMIC_GLOBAL_PANCANCER_RECURRENT", "COSMIC_PANCANCER_TUMOR_COUNT",
    "COSMIC_PANCANCER_INCLUDES_SAMPLE_TUMOR",
    "OS4_CONTEXTUAL", "OM4_CONTEXTUAL", "COSMIC_GLOBAL_RECURRENT", "COSMIC_GLOBAL_EVIDENCE_ONLY",
    "COSMIC_TOTAL_OCCURRENCES",
    "COSMIC_GLOBAL_RECURRENCE_SCORE", "COSMIC_MATCHED_TUMOR_RECURRENCE_SCORE",
    "COSMIC_TUMOR_SPECIFICITY_SCORE", "COSMIC_CONTEXT_SUPPORT_SCORE",
    "COSMIC_TUMOR_CONTEXT_REASON", "COSMIC_TUMOR_TYPES", "COSMIC_EVIDENCE_SUMMARY"), names(variants))
  cm <- variants[isTRUE_vec(variants$COSMIC_MATCH), cosmic_cols, drop = FALSE]
  write_tsv(cm, file.path(dirs$tables, "cosmic_matches.tsv.gz"))
  audit
}

#' @keywords internal
build_run_manifest <- function(cfg, run_dir, rid, build, variants, n_input, audit, t0) {
  status_counts <- as.list(table(factor(variants$filter_status, c("PASS", "REVIEW", "FAIL"))))
  list(
    run_id = rid,
    run_dir = normalizePath(run_dir, mustWork = FALSE),
    package_version = tumoronly_version(),
    date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    genome_build = build,
    input_file = normalizePath(cfg$input$vcf, mustWork = FALSE),
    input_sha256 = file_sha256(cfg$input$vcf),
    fasta_file = cfg_get(cfg, c("reference", "fasta"), NA_character_),
    cosmic_release = cfg_get(cfg, c("cosmic", "release"), NA_character_),
    n_input_variants = n_input,
    n_standardized = nrow(variants),
    status_counts = status_counts,
    n_cosmic_matches = sum(isTRUE_vec(variants$COSMIC_MATCH)),
    r_version = R.version.string,
    packages = installed_pkg_versions(c("data.table", "maftools", "GenomicRanges")),
    runtime_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
  )
}

#' @keywords internal
installed_pkg_versions <- function(pkgs) {
  out <- lapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                           error = function(e) NA_character_))
  stats::setNames(out, pkgs)
}

# ---- dev loader --------------------------------------------------------------

#' Source every R/ file in dependency-safe order for dev use (when the package
#' is not installed). Used by the CLI as a SINGLE call instead of many source()s.
#' @keywords internal
load_all_r <- function(r_dir = NULL) {
  if (is.null(r_dir)) {
    r_dir <- Sys.getenv("TUMORONLY_R_DIR", unset = "")
    if (r_dir == "") r_dir <- file.path(getwd(), "R")
  }
  files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  # utils/io first (define helpers), workflow last
  first <- file.path(r_dir, c("utils.R", "io.R", "config.R", "input.R"))
  files <- c(intersect(first, files), setdiff(files, first))
  for (f in files) sys.source(f, envir = globalenv())
  invisible(files)
}
