# v2 output contract: traceable evidence columns, user-facing table aliases,
# reproducibility files, and publication-oriented QC figures. These helpers only
# materialize already-computed evidence; they do not change scientific rules or
# final classifications.

#' Add user-facing evidence columns for every classified variant.
#'
#' @param variants data.frame after the audited classification workflow.
#' @param cfg resolved tumoronly configuration.
#' @return `variants` with supporting, opposing, missing-evidence, and textual
#'   explanation columns. Existing scientific classes are not modified.
#' @export
add_tumoronly_evidence_columns <- function(variants, cfg = NULL) {
  if (!is.data.frame(variants) || nrow(variants) == 0) {
    variants$evidence_supporting_classification <- character()
    variants$evidence_against_classification <- character()
    variants$missing_evidence <- character()
    variants$classification_explanation <- character()
    return(variants)
  }

  n <- nrow(variants)
  chr <- function(col, default = NA_character_) {
    if (col %in% names(variants)) as.character(variants[[col]]) else rep(default, n)
  }
  num <- function(col, default = NA_real_) {
    if (col %in% names(variants)) suppressWarnings(as.numeric(variants[[col]])) else rep(default, n)
  }
  lg <- function(col, default = FALSE) {
    if (!(col %in% names(variants))) return(rep(default, n))
    variants[[col]] %in% c(TRUE, "TRUE", "true", "1", 1)
  }
  has_value <- function(x) !is.na(x) & nzchar(as.character(x)) & as.character(x) != "NA"
  add_part <- function(parts, flag, label) {
    if (length(flag) == 1L) flag <- rep(flag, length(parts))
    flag[is.na(flag)] <- FALSE
    label <- as.character(label)
    if (length(label) == 1L) label <- rep(label, length(parts))
    parts[flag] <- ifelse(nzchar(parts[flag]),
                           paste(parts[flag], label[flag], sep = ";"),
                           label[flag])
    parts
  }
  fmt_num <- function(x, digits = 3) {
    ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE, scientific = FALSE))
  }

  fc <- chr("final_class", "uncertain_tumor_only")
  status <- chr("filter_status", "REVIEW")
  primary <- chr("primary_reason", "")
  hard <- chr("hard_filter_reason", "")
  pop <- chr("population_category", "")
  art <- chr("artifact_category", "")
  rec <- chr("recurrence_category", "")
  somatic_score <- num("somatic_score")
  technical_score <- num("technical_evidence_score")
  artifact_score <- num("artifact_score")
  germline_score <- num("germline_score")
  max_af <- num("max_pop_af")
  dp <- num("dp")
  alt_count <- num("alt_count")
  vaf <- num("vaf")
  mbq <- num("mbq")
  mmq <- num("mmq")
  tlod <- num("tlod")
  strand <- num("strand_artifact")
  gene <- chr("gene", "")
  cosmic_match <- lg("COSMIC_MATCH") | lg("cosmic_match")
  hotspot <- lg("hotspot_match")
  gene_driver <- lg("gene_driver_match")
  pon <- lg("pon_flag")
  predictor_cat <- chr("COMPUTATIONAL_EVIDENCE_CATEGORY", "")
  is_somatic_candidate <- fc %in% c("high_confidence_somatic", "probable_somatic")
  is_review_or_uncertain <- fc %in% c("manual_review_required", "uncertain_tumor_only")
  is_germline_like <- fc == "likely_germline"
  is_artifact_like <- fc %in% c("likely_artifact", "technical_fail")

  supporting <- rep("", n)
  supporting <- add_part(supporting, has_value(primary), paste0("primary_rule=", primary))
  supporting <- add_part(supporting, status == "PASS", "retained_by_filtering_workflow")
  supporting <- add_part(supporting, status == "REVIEW", "flagged_for_manual_review")
  supporting <- add_part(supporting, is_somatic_candidate & pop == "population_rare_or_absent", "population_frequency_rare_or_absent")
  supporting <- add_part(supporting, is_germline_like & pop %in% c("population_common", "population_low_frequency"), paste0("population_evidence=", pop))
  supporting <- add_part(supporting, is_germline_like & !is.na(germline_score) & germline_score >= 0.60, paste0("germline_score=", fmt_num(germline_score)))
  supporting <- add_part(supporting, is_germline_like & !is.na(vaf) & (abs(vaf - 0.5) <= 0.10 | abs(vaf - 1.0) <= 0.10), "germline_like_VAF")
  supporting <- add_part(supporting, is_artifact_like & hard != "PASS" & has_value(hard), paste0("technical_filter=", hard))
  supporting <- add_part(supporting, art == "artifact_strong", "strong_artifact_evidence")
  supporting <- add_part(supporting, is_artifact_like & !is.na(artifact_score) & artifact_score >= 0.60, paste0("artifact_score=", fmt_num(artifact_score)))
  supporting <- add_part(supporting, (is_somatic_candidate | is_review_or_uncertain) & !is.na(technical_score) & technical_score >= 0.55, paste0("technical_score=", fmt_num(technical_score)))
  supporting <- add_part(supporting, is_somatic_candidate & !is.na(somatic_score), paste0("somatic_score=", fmt_num(somatic_score)))
  supporting <- add_part(supporting, (is_somatic_candidate | is_review_or_uncertain) & !is.na(vaf) & !is.na(dp) & !is.na(alt_count), paste0("VAF=", fmt_num(vaf), "_DP=", fmt_num(dp, 0), "_ALT=", fmt_num(alt_count, 0)))
  supporting <- add_part(supporting, is_somatic_candidate & cosmic_match, "COSMIC_coordinate_match")
  supporting <- add_part(supporting, is_somatic_candidate & hotspot, "hotspot_reference_match")
  supporting <- add_part(supporting, is_somatic_candidate & gene_driver, "driver_gene_reference_match")
  supporting <- add_part(supporting, is_somatic_candidate & rec == "recurrence_tumor_type_supported", "tumor_type_recurrence_support")
  supporting <- add_part(supporting, is_somatic_candidate & predictor_cat %in% c("strong_support", "moderate_support"), paste0("computational_predictors=", predictor_cat))
  supporting[!nzchar(supporting)] <- "no_positive_rule_evidence_recorded"

  against <- rep("", n)
  against <- add_part(against, fc %in% c("high_confidence_somatic", "probable_somatic", "manual_review_required", "uncertain_tumor_only") & pop %in% c("population_common", "population_low_frequency"), paste0("population_evidence=", pop))
  against <- add_part(against, fc %in% c("high_confidence_somatic", "probable_somatic", "manual_review_required", "uncertain_tumor_only") & art %in% c("artifact_strong", "artifact_possible"), paste0("artifact_evidence=", art))
  against <- add_part(against, hard != "PASS" & has_value(hard), paste0("technical_filter=", hard))
  against <- add_part(against, rec == "recurrence_artifact_suspected", "cohort_recurrence_artifact_signal")
  against <- add_part(against, pon, "panel_of_normals_flag")
  against <- add_part(against, !is.na(max_af) & max_af >= cfg_get(cfg, c("population_filters", "rare_af"), 0.001), paste0("max_population_AF=", fmt_num(max_af, 5)))
  against <- add_part(against, !is.na(artifact_score) & artifact_score >= 0.60, paste0("artifact_score=", fmt_num(artifact_score)))
  against <- add_part(against, !is.na(germline_score) & germline_score >= 0.60, paste0("germline_score=", fmt_num(germline_score)))
  against <- add_part(against, is_artifact_like & !is.na(somatic_score) & somatic_score >= 0.55, paste0("somatic_score=", fmt_num(somatic_score)))
  against <- add_part(against, is_artifact_like & !is.na(technical_score) & technical_score >= 0.55, paste0("technical_score=", fmt_num(technical_score)))
  against[!nzchar(against)] <- "no_contrary_rule_evidence_recorded"

  total_samples <- num("total_samples")
  missing <- rep("matched_normal_absent", n)
  missing <- add_part(missing, is.na(tlod), "TLOD_missing")
  missing <- add_part(missing, is.na(strand), "strand_bias_score_missing")
  missing <- add_part(missing, is.na(mmq), "mapping_quality_missing")
  missing <- add_part(missing, is.na(mbq), "base_quality_missing")
  missing <- add_part(missing, is.na(max_af), "population_frequency_missing")
  missing <- add_part(missing, !has_value(gene), "gene_annotation_missing")
  missing <- add_part(missing, is.na(total_samples) | total_samples < cfg_get(cfg, c("cohort", "min_cohort_size_for_recurrence"), 10), "cohort_recurrence_not_evaluable")
  missing <- add_part(missing, is.null(cfg_get(cfg, c("cosmic", "processed_db"), NULL)) || is.na(cfg_get(cfg, c("cosmic", "processed_db"), NA_character_)), "COSMIC_database_not_configured")
  missing <- add_part(missing, is.null(cfg_get(cfg, c("driver_resources", "driver_genes"), NULL)) && is.null(cfg_get(cfg, c("driver_resources", "hotspots"), NULL)), "driver_resources_not_configured")
  missing <- add_part(missing, predictor_cat %in% c("", "insufficient", "not_applicable"), "functional_predictor_support_limited")

  variants$evidence_supporting_classification <- supporting
  variants$evidence_against_classification <- against
  variants$missing_evidence <- missing
  variants$classification_explanation <- vapply(seq_len(n), function(i) {
    lead <- switch(fc[[i]],
      high_confidence_somatic = "Classified as high-confidence somatic candidate",
      probable_somatic = "Classified as probable somatic candidate",
      likely_germline = "Consistent with a likely germline event",
      likely_artifact = "Flagged as a likely technical artifact",
      technical_fail = "Excluded by technical fail-safe filters",
      manual_review_required = "Requires manual review because evidence is conflicting or incomplete",
      uncertain_tumor_only = "Remains uncertain in the tumor-only setting",
      "Classified by the tumor-only rules")
    paste0(
      lead,
      ". Supporting evidence: ", supporting[[i]],
      ". Contrary evidence: ", against[[i]],
      ". Missing or limited evidence: ", missing[[i]],
      ". Tumor-only results are probabilistic and require matched-normal or orthogonal confirmation for clinical interpretation."
    )
  }, character(1))
  variants
}

#' Write v2 output artifacts for a completed run.
#' @keywords internal
write_tumoronly_v2_artifacts <- function(variants, cfg, dirs, manifest, audit) {
  for (d in c(dirs$figures, dirs$figure_data, dirs$logs)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  write_resolved_config(cfg, file.path(dirs$root, "config_used.yaml"))
  write_session_info_file(file.path(dirs$root, "session_info.txt"))

  write_v2_tables(variants, dirs, audit)
  warning_table <- build_v2_warning_table(variants, cfg)
  write_tsv(warning_table, file.path(dirs$logs, "warnings.tsv"))

  figure_manifest <- create_tumoronly_figures(variants, audit, dirs$figures,
                                              dirs$figure_data, cfg = cfg)
  write_tsv(figure_manifest, file.path(dirs$figures, "figure_manifest.tsv"))

  run_manifest <- build_v2_run_manifest(cfg, dirs$root, manifest, warning_table,
                                        figure_manifest)
  jsonlite::write_json(run_manifest, file.path(dirs$root, "run_manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")

  invisible(list(
    warnings = warning_table,
    figures = figure_manifest,
    summary = list(
      warning_count = nrow(warning_table),
      figure_count = sum(figure_manifest$status == "ok"),
      artifact_contract = "v2"
    )
  ))
}

write_v2_tables <- function(variants, dirs, audit) {
  tables <- dirs$tables
  key_cols <- v2_table_columns(variants)
  classified <- variants[, key_cols, drop = FALSE]
  write_tsv(variants, file.path(tables, "all_variants.tsv"))
  write_tsv(classified, file.path(tables, "classified_variants.tsv"))
  write_tsv(variants[variants$final_class == "high_confidence_somatic", key_cols, drop = FALSE],
            file.path(tables, "high_confidence_somatic.tsv"))
  write_tsv(variants[variants$final_class == "likely_germline", key_cols, drop = FALSE],
            file.path(tables, "likely_germline.tsv"))
  write_tsv(variants[variants$final_class == "likely_artifact", key_cols, drop = FALSE],
            file.path(tables, "likely_artifact.tsv"))
  driver_classes <- c("known_driver", "probable_driver", "possible_driver",
                      "uncertain_possible_driver")
  drv <- if ("driver_class" %in% names(variants))
    variants[variants$driver_class %in% driver_classes, key_cols, drop = FALSE]
    else variants[0, key_cols, drop = FALSE]
  write_tsv(drv, file.path(tables, "known_drivers.tsv"))
  write_tsv(build_v2_sample_summary(variants), file.path(tables, "sample_summary.tsv"))
  write_tsv(audit, file.path(tables, "filter_audit.tsv"))
  invisible(NULL)
}

v2_table_columns <- function(variants) {
  intersect(c(
    "sample_id", "tumor_type", "variant_id", "chrom", "pos", "ref", "alt",
    "gene", "consequence", "protein_change", "HGVSP", "vaf", "dp",
    "alt_count", "ref_count", "FILTER", "mutect_filter", "tlod", "mbq", "mmq",
    "strand_artifact", "max_pop_af", "abraom_af", "sabe_af", "gnomad_af",
    "filter_status", "final_class", "somatic_score", "technical_evidence_score",
    "germline_score", "artifact_score", "population_category",
    "artifact_category", "recurrence_category", "oncogenic_category",
    "driver_class", "driver_score", "driver_evidence", "hotspot_match",
    "gene_driver_match", "COSMIC_MATCH", "COSMIC_TOTAL_OCCURRENCES",
    "COSMIC_TUMOR_CONTEXT_STATUS", "PRIORITY", "PRIORITY_SCORE_BASE",
    "AUTHENTICITY", "VARIANT_AUTHENTICITY_SCORE",
    "tmb_countable", "purity_used", "total_cn_used", "multiplicity_used",
    "copy_number_available", "multiplicity_available", "ccf_estimate",
    "ccf_capped", "clonality_class", "clonality_method",
    "clonality_cluster_id", "clonality_cluster_center",
    "clonality_cluster_label", "ml_status", "ml_objective", "ml_model_id",
    "ml_true_positive_probability", "ml_true_variant_probability",
    "ml_prediction_explanation",
    "evidence_supporting_classification", "evidence_against_classification",
    "missing_evidence", "classification_explanation", "filter_reasons",
    "recommended_action"
  ), names(variants))
}

build_v2_sample_summary <- function(variants) {
  if (nrow(variants) == 0) {
    return(data.frame(sample_id = character(), n_variants = integer(),
                      stringsAsFactors = FALSE))
  }
  dt <- data.table::as.data.table(variants)
  cls <- dt[, .N, by = .(sample_id, final_class)]
  class_cols <- as.data.frame(data.table::dcast(cls, sample_id ~ final_class,
                                                value.var = "N", fill = 0))
  base <- as.data.frame(dt[, .(
    n_variants = .N,
    n_pass = sum(filter_status == "PASS", na.rm = TRUE),
    n_review = sum(filter_status == "REVIEW", na.rm = TRUE),
    n_fail = sum(filter_status == "FAIL", na.rm = TRUE),
    median_depth = safe_median(dp),
    median_vaf = safe_median(vaf),
    missing_depth_fraction = round(mean(is.na(dp)), 4),
    missing_vaf_fraction = round(mean(is.na(vaf)), 4),
    missing_population_af_fraction = round(mean(is.na(max_pop_af)), 4)
  ), by = sample_id])
  merge(base, class_cols, by = "sample_id", all.x = TRUE, sort = FALSE)
}

build_v2_warning_table <- function(variants, cfg) {
  n <- nrow(variants)
  sample_n <- if ("sample_id" %in% names(variants))
    length(unique(stats::na.omit(variants$sample_id))) else NA_integer_
  rows <- list()
  add <- function(code, severity, message, n_affected = NA_integer_) {
    rows[[length(rows) + 1]] <<- data.frame(
      code = code, severity = severity, n_affected = n_affected,
      message = message, stringsAsFactors = FALSE)
  }
  min_cohort <- cfg_get(cfg, c("cohort", "min_cohort_size_for_recurrence"), 10)
  if (is.na(sample_n) || sample_n < min_cohort) {
    add("W_COHORT_RECURRENCE_NOT_EVALUABLE", "warning",
        sprintf("Cohort recurrence evidence requires at least %d samples; this run has %s.",
                min_cohort, as.character(sample_n)), n)
  }
  if ("tlod" %in% names(variants) && all(is.na(variants$tlod))) {
    add("W_TLOD_MISSING", "warning",
        "TLOD is absent for all variants; caller-specific quality evidence is incomplete.", n)
  }
  if ("max_pop_af" %in% names(variants) && all(is.na(variants$max_pop_af))) {
    add("W_POPULATION_AF_MISSING", "warning",
        "No configured population-frequency evidence was detected.", n)
  }
  cosmic_db <- cfg_get(cfg, c("cosmic", "processed_db"), NULL)
  if (is.null(cosmic_db) || is.na(cosmic_db) || !file.exists(cosmic_db)) {
    add("W_COSMIC_NOT_CONFIGURED", "warning",
        "No processed COSMIC database was configured; COSMIC evidence is unavailable.", n)
  }
  driver_genes <- cfg_get(cfg, c("driver_resources", "driver_genes"), NULL)
  hotspots <- cfg_get(cfg, c("driver_resources", "hotspots"), NULL)
  if ((is.null(driver_genes) || is.na(driver_genes) || !file.exists(driver_genes)) &&
      (is.null(hotspots) || is.na(hotspots) || !file.exists(hotspots))) {
    add("W_DRIVER_RESOURCES_NOT_CONFIGURED", "warning",
        "No driver-gene or hotspot resource was configured; driver classification is limited.", n)
  }
  if ("gene" %in% names(variants)) {
    miss_gene <- sum(is.na(variants$gene) | variants$gene == "")
    if (miss_gene > 0) add("W_GENE_ANNOTATION_MISSING", "warning",
                           "Some variants lack gene annotation.", miss_gene)
  }
  if (isTRUE(cfg_get(cfg, c("tmb", "enabled"), TRUE))) {
    callable <- suppressWarnings(as.numeric(cfg_get(cfg, c("tmb", "callable_mb"), NA_real_)))
    require_callable <- isTRUE(cfg_get(cfg, c("tmb", "require_callable_mb"), TRUE))
    if (require_callable && (length(callable) != 1L || is.na(callable) || callable <= 0)) {
      add("W_TMB_CALLABLE_MB_MISSING", "warning",
          "TMB module ran as countable variant burden only because tmb.callable_mb is missing or invalid.", n)
    }
  }
  if ("clonality_method" %in% names(variants) &&
      any(variants$clonality_method == "vaf_proxy_no_purity", na.rm = TRUE)) {
    add("W_CLONALITY_VAF_PROXY_ONLY", "warning",
        "Clonality module used VAF-proxy classes for at least one variant because purity was unavailable.",
        sum(variants$clonality_method == "vaf_proxy_no_purity", na.rm = TRUE))
  }
  if ("ml_status" %in% names(variants) &&
      any(variants$ml_status == "not_evaluable_no_active_model", na.rm = TRUE)) {
    add("W_ML_NO_ACTIVE_MODEL", "warning",
        "ML module is enabled but no activated model was available; predictions were not generated.", n)
  }
  if (isTRUE(cfg_get(cfg, c("ancestry", "enabled"), FALSE))) {
    panel <- cfg_get(cfg, c("ancestry", "marker_panel", "path"), NULL)
    if (is.null(panel) || is.na(panel) || !file.exists(panel)) {
      add("W_ANCESTRY_PANEL_MISSING", "warning",
          "Ancestry module is enabled but no AIMs marker panel was available; ancestry is not evaluable.", n)
    }
  }
  if (length(rows) == 0) return(data.frame(
    code = character(), severity = character(), n_affected = integer(),
    message = character(), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

write_session_info_file <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(utils::sessionInfo()), path)
  invisible(path)
}

build_v2_run_manifest <- function(cfg, run_dir, manifest, warnings, figures) {
  list(
    tool = list(name = "tumoronly", version = tumoronly_version(),
                commit = git_commit_id()),
    run = list(
      id = manifest$run_id,
      status = if (nrow(warnings) > 0) "completed_with_warnings" else "completed",
      date = manifest$date,
      runtime_seconds = manifest$runtime_seconds,
      run_dir = normalizePath(run_dir, mustWork = FALSE)
    ),
    input = list(
      files = normalizePath(manifest$input_file, mustWork = FALSE),
      sha256 = manifest$input_sha256,
      format = manifest$input_format,
      genome_build = manifest$genome_build
    ),
    parameters = list(
      technical_filters = cfg_get(cfg, "technical_filters", list()),
      hard_filters = cfg_get(cfg, "hard_filters", list()),
      population_filters = cfg_get(cfg, "population_filters", list()),
      scoring = cfg_get(cfg, "scoring", list()),
      cohort = cfg_get(cfg, "cohort", list())
    ),
    external_resources = list(
      cosmic_release = cfg_get(cfg, c("cosmic", "release"), NA_character_),
      cosmic_processed_db = cfg_get(cfg, c("cosmic", "processed_db"), NA_character_),
      driver_genes = cfg_get(cfg, c("driver_resources", "driver_genes"), NA_character_),
      hotspots = cfg_get(cfg, c("driver_resources", "hotspots"), NA_character_),
      brazilian_population_db = cfg_get(cfg, c("population", "brazilian", "db"), NA_character_)
    ),
    counts = list(
      n_input_variants = manifest$n_input_variants,
      n_standardized = manifest$n_standardized,
      status_counts = manifest$status_counts,
      n_cosmic_matches = manifest$n_cosmic_matches
    ),
    outputs = list(
      report = file.path(normalizePath(run_dir, mustWork = FALSE), "report.html"),
      tables = file.path(normalizePath(run_dir, mustWork = FALSE), "tables"),
      figures = file.path(normalizePath(run_dir, mustWork = FALSE), "figures"),
      figure_data = file.path(normalizePath(run_dir, mustWork = FALSE), "figure_data"),
      module_status = file.path(normalizePath(run_dir, mustWork = FALSE), "tables", "module_status.tsv"),
      tmb_summary = file.path(normalizePath(run_dir, mustWork = FALSE), "tables", "tmb_summary.tsv"),
      clonality_summary = file.path(normalizePath(run_dir, mustWork = FALSE), "tables", "clonality_summary.tsv"),
      ml_status = file.path(normalizePath(run_dir, mustWork = FALSE), "tables", "ml_status.tsv")
    ),
    warnings = warnings,
    figures = figures,
    environment = list(
      r_version = R.version.string,
      platform = R.version$platform,
      system = as.list(Sys.info()),
      packages = installed_pkg_versions(c("data.table", "yaml", "jsonlite",
                                          "digest", "ggplot2", "maftools",
                                          "rmarkdown"))
    )
  )
}

update_tumoronly_run_manifests <- function(run_dir, manifest, report_rendered = NA) {
  jsonlite::write_json(manifest, file.path(run_dir, "manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")
  path <- file.path(run_dir, "run_manifest.json")
  if (file.exists(path)) {
    v2 <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE),
                   error = function(e) NULL)
    if (!is.null(v2)) {
      v2$run$runtime_seconds <- manifest$runtime_seconds
      v2$run$report_rendered <- isTRUE(report_rendered)
      v2$outputs$report <- file.path(normalizePath(run_dir, mustWork = FALSE), "report.html")
      jsonlite::write_json(v2, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    }
  }
  invisible(NULL)
}

git_commit_id <- function() {
  root <- normalizePath(getwd(), mustWork = FALSE)
  if (!nzchar(Sys.which("git")) || !dir.exists(file.path(root, ".git"))) {
    return(NA_character_)
  }
  out <- suppressWarnings(tryCatch(
    system2("git", c("-C", root, "rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    warning = function(w) NA_character_,
    error = function(e) NA_character_))
  if (length(out) == 0 || is.na(out[[1]]) || !nzchar(out[[1]])) NA_character_ else out[[1]]
}

#' Generate v2 QC figures and figure-data tables.
#'
#' @param variants classified variant table.
#' @param audit filter audit table.
#' @param figures_dir destination for PDF/SVG/PNG figures.
#' @param figure_data_dir destination for per-figure TSV data.
#' @param cfg resolved config.
#' @return a manifest with one row per attempted figure.
#' @export
create_tumoronly_figures <- function(variants, audit, figures_dir, figure_data_dir,
                                     cfg = NULL) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  add <- function(id, title, status, files = character(), data = NA_character_,
                  reason = NA_character_) {
    manifest[[length(manifest) + 1]] <<- data.frame(
      figure_id = id, title = title, status = status,
      files = paste(files, collapse = ";"), figure_data = data,
      reason = reason, stringsAsFactors = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    add("all", "v2 figures", "skipped", reason = "ggplot2 is not installed")
    return(do.call(rbind, manifest))
  }
  palette <- v2_class_palette()
  theme_pub <- function() {
    ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(color = "grey30"),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9),
        legend.key.width = grid::unit(0.8, "lines"),
        plot.margin = ggplot2::margin(10, 35, 10, 10)
      )
  }
  save_figure <- function(id, title, plot, data, width = 7.5, height = 5.2) {
    data_path <- file.path(figure_data_dir, paste0(id, ".tsv"))
    write_tsv(data, data_path)
    base <- file.path(figures_dir, id)
    files <- c(paste0(base, ".pdf"), paste0(base, ".svg"), paste0(base, ".png"))
    ok_files <- character()
    for (f in files) {
      dev_opened <- FALSE
      ok <- tryCatch({
        ext <- tolower(tools::file_ext(f))
        if (ext == "pdf") {
          grDevices::pdf(f, width = width, height = height)
        } else if (ext == "svg") {
          suppressWarnings(grDevices::svg(f, width = width, height = height))
        } else if (ext == "png") {
          grDevices::png(f, width = width, height = height, units = "in", res = 300, type = "cairo")
        } else {
          stop("unsupported figure extension")
        }
        dev_opened <- TRUE
        print(plot)
        grDevices::dev.off()
        dev_opened <- FALSE
        TRUE
      }, error = function(e) FALSE, finally = {
        if (isTRUE(dev_opened) && grDevices::dev.cur() > 1) {
          try(grDevices::dev.off(), silent = TRUE)
        }
      })
      if (tolower(tools::file_ext(f)) == "svg" &&
          (!isTRUE(ok) || !file.exists(f) || file.info(f)$size == 0)) {
        ok <- tryCatch({
          write_basic_svg_from_data(f, title, data, width = width, height = height)
          TRUE
        }, error = function(e) FALSE)
      }
      if (ok && file.exists(f) && file.info(f)$size > 0) ok_files <- c(ok_files, f)
    }
    if (length(ok_files) == 0) {
      add(id, title, "failed", data = data_path, reason = "no export format succeeded")
    } else {
      add(id, title, "ok", files = ok_files, data = data_path)
    }
  }

  audit_fig <- audit[audit$filter %in% c("technical", "population_common",
    "technical_artifact", "cohort_recurrence_artifact", "final_PASS",
    "final_REVIEW", "final_FAIL"), , drop = FALSE]
  if (nrow(audit_fig) > 0) {
    audit_fig$filter <- factor(audit_fig$filter, levels = rev(audit_fig$filter))
    p <- ggplot2::ggplot(audit_fig, ggplot2::aes(filter, n_variants)) +
      ggplot2::geom_col(fill = "#0072B2", width = 0.72) +
      ggplot2::geom_text(ggplot2::aes(label = n_variants), hjust = -0.12, size = 3) +
      ggplot2::labs(title = "A. Filtering workflow",
        subtitle = "Counts are audit categories; retained variants are PASS or REVIEW.",
        x = NULL, y = "Variants") +
      ggplot2::coord_flip(ylim = c(0, max(audit_fig$n_variants, na.rm = TRUE) * 1.12),
                          clip = "off") +
      theme_pub()
    save_figure("figure_01_filtering_workflow", "Filtering workflow", p, audit_fig)
  } else {
    add("figure_01_filtering_workflow", "Filtering workflow", "skipped",
        reason = "filter audit is empty")
  }

  qc_cols <- intersect(c("dp", "vaf", "alt_count", "mbq", "mmq", "tlod",
                         "somatic_score"), names(variants))
  qc <- data.frame()
  for (col in qc_cols) {
    val <- suppressWarnings(as.numeric(variants[[col]]))
    val <- val[!is.na(val) & is.finite(val)]
    if (length(val) > 0) qc <- rbind(qc, data.frame(metric = col, value = val))
  }
  if (nrow(qc) > 0) {
    p <- ggplot2::ggplot(qc, ggplot2::aes(value)) +
      ggplot2::geom_histogram(bins = 40, fill = "#009E73", color = "white", size = 0.15) +
      ggplot2::facet_wrap(~ metric, scales = "free", ncol = 2) +
      ggplot2::labs(title = "B. Quality-control and score distributions",
        subtitle = "Only available, numeric metric values are plotted.",
        x = "Metric value", y = "Variants") +
      theme_pub()
    save_figure("figure_02_qc_overview", "QC overview", p, qc, width = 8, height = 6)
  } else {
    add("figure_02_qc_overview", "QC overview", "skipped",
        reason = "no numeric QC metrics available")
  }

  cls <- as.data.frame(table(final_class = variants$final_class), stringsAsFactors = FALSE)
  names(cls) <- c("final_class", "n_variants")
  cls <- cls[order(-cls$n_variants), , drop = FALSE]
  if (nrow(cls) > 0) {
    cls_palette <- palette[intersect(names(palette), cls$final_class)]
    unknown_cls <- setdiff(cls$final_class, names(cls_palette))
    if (length(unknown_cls) > 0) cls_palette <- c(cls_palette, stats::setNames(rep("#666666", length(unknown_cls)), unknown_cls))
    p <- ggplot2::ggplot(cls, ggplot2::aes(stats::reorder(final_class, n_variants), n_variants,
                                           fill = final_class)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = cls_palette, breaks = names(cls_palette), na.value = "grey70") +
      ggplot2::labs(title = "C. Final tumor-only classification",
        subtitle = "Classes are inference categories, not clinical confirmation.",
        x = NULL, y = "Variants") +
      ggplot2::guides(fill = "none") +
      theme_pub()
    save_figure("figure_03_classification", "Final classification", p, cls)
  }

  scatter <- variants[!is.na(variants$vaf) & !is.na(variants$dp) & variants$dp > 0, , drop = FALSE]
  if (nrow(scatter) > 0) {
    max_points <- cfg_get(cfg, c("visualization", "max_plot_points"),
                          cfg_get(cfg, c("report", "max_plot_points"), 200000))
    if (nrow(scatter) > max_points) {
      set.seed(cfg_get(cfg, c("visualization", "seed"), 20260621))
      scatter <- scatter[sample(seq_len(nrow(scatter)), max_points), , drop = FALSE]
    }
    dat <- scatter[, intersect(c("sample_id", "variant_id", "final_class", "vaf",
                                 "dp", "alt_count", "gene"), names(scatter)), drop = FALSE]
    scatter_classes <- unique(as.character(dat$final_class))
    scatter_palette <- palette[intersect(names(palette), scatter_classes)]
    unknown_scatter <- setdiff(scatter_classes, names(scatter_palette))
    if (length(unknown_scatter) > 0) {
      scatter_palette <- c(scatter_palette,
                           stats::setNames(rep("#666666", length(unknown_scatter)),
                                           unknown_scatter))
    }
    p <- ggplot2::ggplot(dat, ggplot2::aes(dp, vaf, color = final_class)) +
      ggplot2::geom_point(alpha = 0.45, size = 0.75) +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::scale_color_manual(values = scatter_palette, breaks = names(scatter_palette),
                                  na.value = "grey50") +
      ggplot2::labs(title = "D. VAF by depth",
        subtitle = "Depth is shown on a log10 scale to keep high-depth outliers visible.",
        x = "Depth (log10)", y = "Variant allele fraction", color = "Class") +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE,
        override.aes = list(alpha = 1, size = 2))) +
      theme_pub()
    save_figure("figure_04_vaf_depth", "VAF by depth", p, dat, width = 7.2, height = 5.4)
  } else {
    add("figure_04_vaf_depth", "VAF by depth", "skipped",
        reason = "depth and VAF are not jointly available")
  }

  sample_cls <- as.data.frame(data.table::as.data.table(variants)[, .N,
    by = .(sample_id, final_class)])
  if (nrow(sample_cls) > 0) {
    sample_classes <- unique(as.character(sample_cls$final_class))
    sample_palette <- palette[intersect(names(palette), sample_classes)]
    unknown_sample <- setdiff(sample_classes, names(sample_palette))
    if (length(unknown_sample) > 0) {
      sample_palette <- c(sample_palette,
                          stats::setNames(rep("#666666", length(unknown_sample)),
                                          unknown_sample))
    }
    p <- ggplot2::ggplot(sample_cls, ggplot2::aes(sample_id, N, fill = final_class)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = sample_palette, breaks = names(sample_palette),
                                 na.value = "grey70") +
      ggplot2::labs(title = "E. Classification burden by sample",
        subtitle = "Counts are variant counts; this is not tumor mutational burden.",
        x = NULL, y = "Variants", fill = "Class") +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
      theme_pub()
    save_figure("figure_05_sample_class_distribution",
                "Classification burden by sample", p, sample_cls, width = 7.5,
                height = max(4.5, min(10, 0.22 * length(unique(sample_cls$sample_id)) + 3)))
  } else {
    add("figure_05_sample_class_distribution", "Classification burden by sample",
        "skipped", reason = "no sample/class counts available")
  }

  miss <- build_missing_evidence_counts(variants)
  if (nrow(miss) > 0) {
    p <- ggplot2::ggplot(miss, ggplot2::aes(stats::reorder(evidence, n_variants), n_variants)) +
      ggplot2::geom_col(fill = "#D55E00", width = 0.72) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "F. Missing or limited evidence",
        subtitle = "Each row may contribute more than one missing-evidence item.",
        x = NULL, y = "Variants") +
      theme_pub()
    save_figure("figure_06_missing_evidence", "Missing evidence", p, miss)
  } else {
    add("figure_06_missing_evidence", "Missing evidence", "skipped",
        reason = "missing_evidence column is empty")
  }

  do.call(rbind, manifest)
}

v2_class_palette <- function() {
  c(
    high_confidence_somatic = "#0072B2",
    probable_somatic = "#56B4E9",
    likely_germline = "#009E73",
    likely_artifact = "#D55E00",
    technical_fail = "#CC79A7",
    manual_review_required = "#E69F00",
    uncertain_tumor_only = "#999999"
  )
}

build_missing_evidence_counts <- function(variants) {
  if (!"missing_evidence" %in% names(variants) || nrow(variants) == 0)
    return(data.frame(evidence = character(), n_variants = integer()))
  pieces <- strsplit(as.character(variants$missing_evidence), ";", fixed = TRUE)
  long <- data.frame(evidence = trimws(unlist(pieces, use.names = FALSE)),
                     stringsAsFactors = FALSE)
  long <- long[nzchar(long$evidence), , drop = FALSE]
  if (nrow(long) == 0) return(data.frame(evidence = character(), n_variants = integer()))
  out <- as.data.frame(table(evidence = long$evidence), stringsAsFactors = FALSE)
  names(out) <- c("evidence", "n_variants")
  out[order(-out$n_variants), , drop = FALSE]
}

write_basic_svg_from_data <- function(path, title, data, width = 7.5, height = 5.2) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  w <- 900
  h <- max(420, as.integer(height / width * w))
  margin <- list(left = 90, right = 40, top = 70, bottom = 80)
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  svg <- c(
    sprintf("<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' viewBox='0 0 %d %d'>", w, h, w, h),
    "<rect width='100%' height='100%' fill='white'/>",
    "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#222}.axis{stroke:#333;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}.label{font-size:12px}.title{font-size:20px;font-weight:bold}.subtitle{font-size:12px;fill:#555}</style>",
    sprintf("<text class='title' x='%d' y='34'>%s</text>", margin$left, esc(title)),
    "<text class='subtitle' x='90' y='54'>Fallback vector SVG generated from saved figure_data.</text>"
  )
  plot_w <- w - margin$left - margin$right
  plot_h <- h - margin$top - margin$bottom
  x0 <- margin$left
  y0 <- margin$top + plot_h
  svg <- c(svg,
           sprintf("<line class='axis' x1='%d' y1='%d' x2='%d' y2='%d'/>", x0, y0, x0 + plot_w, y0),
           sprintf("<line class='axis' x1='%d' y1='%d' x2='%d' y2='%d'/>", x0, margin$top, x0, y0))

  if (is.null(data) || nrow(data) == 0) {
    svg <- c(svg, sprintf("<text class='label' x='%d' y='%d'>No rows available for this figure.</text>",
                          x0 + 20, margin$top + 40))
    writeLines(c(svg, "</svg>"), path)
    return(invisible(path))
  }
  num_cols <- names(data)[vapply(data, is.numeric, logical(1))]
  if (length(num_cols) >= 2) {
    x <- data[[num_cols[[1]]]]
    y <- data[[num_cols[[2]]]]
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) > 2500) {
      set.seed(20260621)
      idx <- sample(seq_along(x), 2500)
      x <- x[idx]; y <- y[idx]
    }
    if (length(x) > 0) {
      xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
      if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
      if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
      sx <- x0 + (x - xr[1]) / diff(xr) * plot_w
      sy <- y0 - (y - yr[1]) / diff(yr) * plot_h
      pts <- sprintf("<circle cx='%.1f' cy='%.1f' r='2.1' fill='#0072B2' opacity='0.45'/>", sx, sy)
      svg <- c(svg, pts,
               sprintf("<text class='label' x='%d' y='%d'>%s</text>", x0, h - 38, esc(num_cols[[1]])),
               sprintf("<text class='label' transform='translate(24,%d) rotate(-90)'>%s</text>", y0, esc(num_cols[[2]])))
    }
  } else if (length(num_cols) >= 1) {
    ycol <- if ("n_variants" %in% num_cols) "n_variants" else num_cols[[1]]
    y <- data[[ycol]]
    char_cols <- names(data)[!vapply(data, is.numeric, logical(1))]
    if (nrow(data) > 50 && ycol == "value") {
      br <- pretty(y[is.finite(y)], n = 25)
      counts <- hist(y[is.finite(y)], breaks = br, plot = FALSE)
      labels <- paste0(head(counts$breaks, -1), "-", tail(counts$breaks, -1))
      plot_data <- data.frame(label = labels, y = counts$counts)
    } else {
      label_col <- if (length(char_cols) > 0) char_cols[[1]] else NULL
      labels <- if (!is.null(label_col)) as.character(data[[label_col]]) else as.character(seq_along(y))
      plot_data <- data.frame(label = labels, y = y)
      plot_data <- plot_data[order(-plot_data$y), , drop = FALSE]
      plot_data <- utils::head(plot_data, 30)
    }
    plot_data <- plot_data[is.finite(plot_data$y), , drop = FALSE]
    if (nrow(plot_data) > 0) {
      ymax <- max(plot_data$y, na.rm = TRUE)
      if (!is.finite(ymax) || ymax <= 0) ymax <- 1
      gap <- 4
      bar_w <- max(3, (plot_w - gap * (nrow(plot_data) - 1)) / max(nrow(plot_data), 1))
      for (i in seq_len(nrow(plot_data))) {
        bh <- plot_data$y[[i]] / ymax * plot_h
        bx <- x0 + (i - 1) * (bar_w + gap)
        by <- y0 - bh
        svg <- c(svg, sprintf("<rect x='%.1f' y='%.1f' width='%.1f' height='%.1f' fill='#0072B2'/>",
                              bx, by, bar_w, bh))
      }
      show_every <- max(1, ceiling(nrow(plot_data) / 12))
      for (i in seq(1, nrow(plot_data), by = show_every)) {
        bx <- x0 + (i - 1) * (bar_w + gap) + bar_w / 2
        lab <- substr(plot_data$label[[i]], 1, 18)
        svg <- c(svg, sprintf("<text class='label' transform='translate(%.1f,%d) rotate(45)' text-anchor='start'>%s</text>",
                              bx, y0 + 14, esc(lab)))
      }
      svg <- c(svg, sprintf("<text class='label' transform='translate(24,%d) rotate(-90)'>%s</text>",
                            y0, esc(ycol)))
    }
  } else {
    svg <- c(svg, sprintf("<text class='label' x='%d' y='%d'>No numeric data available for a chart.</text>",
                          x0 + 20, margin$top + 40))
  }
  writeLines(c(svg, "</svg>"), path)
  invisible(path)
}

copy_report_to_v2_root <- function(run_dir, report_path = NULL) {
  if (is.null(report_path) || is.na(report_path) || !file.exists(report_path)) {
    report_path <- file.path(run_dir, "report", "tumor_only_report.html")
  }
  out <- file.path(run_dir, "report.html")
  if (file.exists(report_path)) {
    file.copy(report_path, out, overwrite = TRUE)
    return(invisible(out))
  }
  invisible(NULL)
}
