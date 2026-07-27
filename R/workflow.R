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

  log_step("run", "reading variants (multi-format ingestion)",
           file = if (length(cfg$input$vcf) > 1)
             sprintf("cohort of %d files", length(cfg$input$vcf)) else cfg$input$vcf,
           run_id = rid)
  ing <- read_variant_input(
    cfg$input$vcf,
    format = cfg_get(cfg, c("input", "format"), "auto"),
    sample_metadata = cfg_get(cfg, c("input", "sample_metadata"), NULL),
    column_map = cfg_get(cfg, c("input", "column_map"), NULL),
    chunk_size = cfg_get(cfg, c("input", "chunk_size"), NULL))
  variants <- ing$variants
  n_input <- nrow(variants)
  write_tsv(ing$schema_report, file.path(dirs$tables, "input_schema_report.tsv"))
  log_step("run", "input ingested", format = ing$format, n = n_input)

  variants <- standardize_variant_table(variants, cfg)
  variants <- attach_sample_metadata(variants, cfg)

  log_step("run", "feature engineering + technical filters", n = nrow(variants))
  variants <- add_basic_features(variants, cfg)
  hf <- add_sample_qc_and_hard_filters(variants, cfg)
  variants <- hf$variants
  sample_qc <- hf$sample_qc
  recurrence <- compute_cohort_recurrence(variants, cfg)
  variants <- merge_recurrence_features(variants, recurrence)

  adaptive_thresholds <- NULL
  if (isTRUE(cfg_get(cfg, c("adaptive_filtering", "enabled"), TRUE))) {
    log_step("run", "adaptive per-sample filtering (experimental)")
    variants <- add_adaptive_filtering(variants, cfg)
    adaptive_thresholds <- attr(variants, "adaptive_thresholds")
  }

  log_step("run", "scoring + COSMIC evidence + classification")
  variants <- score_variants(variants, cfg)
  variants <- annotate_cosmic(variants, cfg, build)
  variants <- add_driver_annotation_layers(variants, cfg)
  variants <- add_guideline_classifications(variants, cfg)
  variants <- classify_variants(variants, cfg)
  variants <- score_and_classify_driver_layers(variants, cfg)

  log_step("run", "functional predictors + prioritization + authenticity")
  variants <- add_predictor_evidence(variants, cfg)
  predictor_inventory <- attr(variants, "predictor_inventory")
  variants <- add_prioritization(variants, cfg)
  variants <- add_variant_authenticity(variants, cfg)

  log_step("run", "Brazilian-aware population evidence")
  variants <- add_brazilian_population_evidence(variants, cfg, build)

  variants <- add_decision_trail(variants)
  variants <- add_unified_labels(variants)
  variants <- apply_brazilian_population_rules(variants, cfg)   # documented population REVIEW rule
  variants <- apply_adaptive_review(variants, cfg)              # documented adaptive-technical REVIEW rule

  # ---- outputs ----------------------------------------------------------
  audit <- write_run_tables(variants, dirs)
  write_tsv(predictor_inventory %||% data.frame(note = "no predictors"),
            file.path(dirs$tables, "predictor_inventory.tsv"))
  write_summary_tables(variants, dirs)
  export_review_template(variants, file.path(dirs$tables, "review_template.tsv"))
  if (!is.null(adaptive_thresholds))
    write_tsv(adaptive_thresholds, file.path(dirs$tables, "adaptive_thresholds.tsv"))
  maf_path <- file.path(dirs$maf, "filtered.maf.gz")
  retained <- variants[variants$filter_status %in% c("PASS", "REVIEW"), , drop = FALSE]
  plot_manifest <- NULL
  if (nrow(retained) > 0) {
    maf <- tryCatch(create_maf(retained, maf_path), error = function(e) {
      log_step("maf", "MAF creation failed", error = conditionMessage(e)); NULL })
    if (!is.null(maf))
      plot_manifest <- tryCatch(create_maftools_plots(maf_path, dirs$plots, cfg,
        sample_metadata = cfg_get(cfg, c("input", "sample_metadata"), NULL)),
        error = function(e) { log_step("plot", "plots failed", error = conditionMessage(e)); NULL })
  } else {
    log_step("maf", "no retained variants; MAF and plots skipped")
  }
  if (!is.null(plot_manifest)) write_tsv(plot_manifest, file.path(dirs$plots, "plot_manifest.tsv"))

  # ---- optional, experimental genetic-ancestry inference ------------------
  if (isTRUE(cfg_get(cfg, c("ancestry", "enabled"), FALSE))) {
    log_step("run", "genetic-ancestry inference (experimental)")
    anc <- tryCatch(infer_ancestry(variants, cfg, build),
                    error = function(e) { log_step("ancestry", "failed", error = conditionMessage(e)); NULL })
    if (!is.null(anc)) {
      write_tsv(anc$summary, file.path(dirs$tables, "ancestry_summary.tsv"))
      if (nrow(anc$snps) > 0) write_tsv(anc$snps, file.path(dirs$tables, "ancestry_snps.tsv.gz"))
      write_tsv(anc$qc, file.path(dirs$tables, "ancestry_qc.tsv"))
      write_tsv(anc$projection, file.path(dirs$tables, "ancestry_projection.tsv"))
      tryCatch(create_ancestry_plots(anc, dirs$plots, cfg),
               error = function(e) log_step("ancestry", "plots failed", error = conditionMessage(e)))
    }
  }

  manifest <- build_run_manifest(cfg, run_dir, rid, build, variants, n_input, audit, t0)
  manifest$input_format <- ing$format
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

#' Write the review-oriented summary tables (top-prioritized, manual review,
#' predictor conflicts, cosmic context matches, and a compact review summary).
#' @keywords internal
write_summary_tables <- function(variants, dirs) {
  keep <- function(cols) intersect(cols, names(variants))
  review_cols <- keep(c("sample_id", "gene", "variant_id", "consequence", "protein_change", "HGVSP",
    "vaf", "dp", "alt_count", "filter_status", "final_class", "PRIORITY_CATEGORY", "PRIORITY_SCORE_BASE",
    "VARIANT_AUTHENTICITY_CATEGORY", "VARIANT_AUTHENTICITY_SCORE",
    "COMPUTATIONAL_EVIDENCE_CATEGORY", "COMPUTATIONAL_CONFLICT_FLAG", "hotspot_match", "gene_driver_match",
    "COSMIC_TUMOR_CONTEXT_STATUS", "COSMIC_MATCHING_TUMOR_OCCURRENCES",
    "CONFIDENCE_CATEGORY_BASE", "filter_reasons"))
  summ <- variants[, review_cols, drop = FALSE]
  if ("PRIORITY_SCORE_BASE" %in% names(summ)) summ <- summ[order(-summ$PRIORITY_SCORE_BASE), ]
  write_tsv(summ, file.path(dirs$tables, "variant_review_summary.tsv"))

  retained <- variants[variants$filter_status %in% c("PASS", "REVIEW"), , drop = FALSE]
  top <- retained
  if ("PRIORITY_SCORE_BASE" %in% names(top)) top <- top[order(-top$PRIORITY_SCORE_BASE), ]
  write_tsv(utils::head(top[, review_cols, drop = FALSE], 500),
            file.path(dirs$tables, "top_prioritized_variants.tsv"))

  mr <- variants[variants$filter_status == "REVIEW" |
    (("final_class" %in% names(variants)) & variants$final_class == "manual_review_required"), review_cols, drop = FALSE]
  write_tsv(mr, file.path(dirs$tables, "manual_review_variants.tsv"))

  conf <- if ("COMPUTATIONAL_CONFLICT_FLAG" %in% names(variants))
    variants[isTRUE_vec(variants$COMPUTATIONAL_CONFLICT_FLAG),
             keep(c("sample_id", "gene", "variant_id", "consequence", "COMPUTATIONAL_EVIDENCE_DETAILS",
                    "COMPUTATIONAL_PREDICTORS_DAMAGING", "COMPUTATIONAL_PREDICTORS_BENIGN",
                    "COMPUTATIONAL_CONFLICT_REASON")), drop = FALSE]
    else variants[0, , drop = FALSE]
  write_tsv(conf, file.path(dirs$tables, "predictor_conflicts.tsv"))

  ctx <- if ("COSMIC_MATCH" %in% names(variants))
    variants[isTRUE_vec(variants$COSMIC_MATCH),
             keep(c("sample_id", "gene", "variant_id", "tumor_type", "COSMIC_TUMOR_CONTEXT_STATUS",
                    "COSMIC_MATCH_LEVEL", "COSMIC_MATCHING_TUMOR_OCCURRENCES", "COSMIC_OTHER_TUMOR_OCCURRENCES",
                    "COSMIC_CONTEXT_EVALUABLE", "COSMIC_CONTEXT_SUPPORT_SCORE")), drop = FALSE]
    else variants[0, , drop = FALSE]
  write_tsv(ctx, file.path(dirs$tables, "cosmic_context_matches.tsv"))

  auth_col <- "VARIANT_AUTHENTICITY_CATEGORY"
  if (auth_col %in% names(variants)) {
    write_tsv(variants[variants[[auth_col]] == "LIKELY_ARTIFACT", review_cols, drop = FALSE],
              file.path(dirs$tables, "likely_artifacts.tsv"))
    write_tsv(variants[variants[[auth_col]] == "LIKELY_TRUE", review_cols, drop = FALSE],
              file.path(dirs$tables, "likely_true_variants.tsv"))
  }

  # ---- Brazilian population evidence tables -------------------------------
  if ("POPULATION_EVIDENCE_STATUS" %in% names(variants)) {
    bcols <- keep(c("sample_id", "gene", "variant_id", "consequence", "vaf",
      "GLOBAL_POPULATION_STATUS", "BRAZILIAN_POPULATION_STATUS", "POPULATION_EVIDENCE_STATUS",
      "POPULATION_EVIDENCE_CONFIDENCE", "MAX_GLOBAL_AF", "MAX_ANCESTRY_SPECIFIC_AF",
      "BRAZILIAN_AF", "BRAZILIAN_AC", "BRAZILIAN_AN", "BRAZILIAN_HOM_COUNT",
      "LOCAL_CONTROL_AF", "BRAZILIAN_DB_ABSENCE_INTERPRETABLE",
      "POSSIBLE_BRAZILIAN_GERMLINE", "BRAZILIAN_GERMLINE_EVIDENCE",
      "BRAZILIAN_POPULATION_EVIDENCE_SCORE", "STATUS", "filter_status"))
    write_tsv(variants[, bcols, drop = FALSE], file.path(dirs$tables, "brazilian_population_evidence.tsv"))
    write_tsv(variants[isTRUE_vec(variants$POSSIBLE_BRAZILIAN_GERMLINE), bcols, drop = FALSE],
              file.path(dirs$tables, "possible_brazilian_germline.tsv"))
    write_tsv(variants[isTRUE_vec(variants$POSSIBLE_BRAZILIAN_GERMLINE) |
                       variants$POPULATION_EVIDENCE_STATUS == "insufficient_coverage", bcols, drop = FALSE],
              file.path(dirs$tables, "brazilian_population_review.tsv"))
    conflict <- variants$POPULATION_EVIDENCE_STATUS %in%
      c("absent_global_present_brazilian", "absent_brazilian_present_global")
    write_tsv(variants[conflict, bcols, drop = FALSE],
              file.path(dirs$tables, "global_brazilian_frequency_conflicts.tsv"))
    lc <- variants[!is.na(variants$LOCAL_CONTROL_AF) & variants$LOCAL_CONTROL_AF > 0, bcols, drop = FALSE]
    write_tsv(lc, file.path(dirs$tables, "local_normal_recurrence.tsv"))
  }
  invisible(NULL)
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
