# v2 user-facing wrappers: config template, input validation, and a thin API
# around the existing audited run_tumor_only() orchestrator.

#' Return the v2 default configuration as a list.
#' @export
tumoronly_default_config <- function() {
  list(
    analysis = list(run_id = NULL, output_dir = "results", threads = 1L),
    input = list(
      path = NULL,
      format = "auto",
      genome_build = NULL,
      sample_metadata = NULL,
      chunk_size = NULL,
      column_map = NULL
    ),
    validation = list(strict = TRUE),
    reference = list(fasta = NULL, chain_file = NULL),
    cosmic = list(
      processed_db = NULL, raw_file = NULL, release = NULL,
      source_build = NULL, target_build = NULL, cache_dir = "db/cosmic",
      rebuild = FALSE, tumor_type_mapping = "config/cosmic_tumor_type_mapping.tsv"
    ),
    cohort = list(
      cancer_type = "PANCANCER",
      recurrent_variant_fraction_artifact = 0.30,
      recurrent_locus_fraction_artifact = 0.40,
      min_cohort_size_for_recurrence = 10,
      min_recurrent_samples = 3,
      min_recurrent_locus_samples = 3
    ),
    cancer = list(default_tumor_type = "PANCANCER"),
    hard_filters = list(
      enabled = TRUE, assay = "WGS",
      adaptive_depth_fraction = 0.35,
      adaptive_q25_depth_fraction = 0.50,
      missing_required_metrics_fail = FALSE,
      caller_filter_accepted_values = c("PASS")
    ),
    adaptive_filtering = list(enabled = TRUE, escalate_to_review = TRUE),
    technical_filters = list(
      min_depth = 20,
      min_alt_count_snv = 5,
      min_alt_count_indel = 8,
      min_af_snv = 0.03,
      min_af_indel = 0.05,
      min_tlod = 6,
      min_mbq = 25,
      min_mmq = 40
    ),
    population_filters = list(
      common_af = 0.005,
      rare_af = 0.001,
      use_population_columns = c(
        "GMAF", "gnomADe_AF", "gnomADg_AF", "ABraOM_AF",
        "AbraOM_AF", "ABraOM_MAF", "BIPMed_AF")
    ),
    scoring = list(
      high_confidence_somatic = 0.75,
      probable_somatic = 0.55,
      probable_germline = 0.70,
      probable_artifact = 0.70,
      validation_rescue_somatic_min = 0.45
    ),
    guideline_classification = list(
      enabled = TRUE,
      somatic = list(
        benign_maf_very_strong = 0.05,
        benign_maf_strong = 0.01,
        absent_population_af = 0.00001,
        cosmic_moderate_count = 3,
        cosmic_high_count = 10
      ),
      germline = list(ba1_af = 0.05, bs1_af = 0.01, pm2_af = 0.00001)
    ),
    driver_resources = list(
      driver_genes = NULL, driver_genes_delimiter = "\t",
      hotspots = NULL, hotspots_delimiter = "\t"
    ),
    driver_scoring = list(
      known_driver = 0.75, probable_driver = 0.60,
      possible_driver = 0.40, min_tumor_type_samples_for_recurrence_driver = 5
    ),
    predictors = list(
      detect_automatically = TRUE,
      use_for_prioritization = TRUE,
      use_for_filtering = FALSE
    ),
    authenticity = list(likely_true_min = 0.66, likely_artifact_max = 0.40),
    population = list(
      global = list(common_af = 0.01, rare_af = 0.001),
      brazilian = list(
        db = NULL, raw_file = NULL, cache_dir = "db/abraom",
        name = "ABraOM_SABE", genome_build = NULL,
        common_af = 0.01, review_af = 0.001,
        minimum_allele_count = 2, cohort_size = NULL,
        escalate_to_review = TRUE
      ),
      local_controls = list(path = NULL, artifact_or_germline_af = 0.005,
                            minimum_observations = 2)
    ),
    ancestry = list(enabled = FALSE),
    plots = list(
      minimum_samples_for_oncoplot = 2,
      minimum_samples_for_interactions = 20,
      minimum_mutations_for_rainfall = 50,
      minimum_protein_changes_for_lollipop = 3,
      oncoplot = list(top = 20, sort_by_annotation = TRUE,
                      remove_non_mutated = FALSE)
    ),
    report = list(format = c("html"), theme = "publication", top_variants = 30)
  )
}

#' Write a v2 tumoronly YAML configuration template.
#' @export
write_tumoronly_config_template <- function(output, force = FALSE) {
  if (missing(output) || is.null(output) || !nzchar(output))
    stop("write_tumoronly_config_template(): provide `output`.", call. = FALSE)
  if (file.exists(output) && !isTRUE(force))
    stop("Config file already exists: ", output, ". Use force=TRUE to overwrite.", call. = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(tumoronly_default_config(), output)
  invisible(normalizePath(output, mustWork = FALSE))
}

#' Validate an input file against the tumoronly canonical schema.
#' @export
validate_tumoronly_input <- function(input = NULL, config = NULL, strict = TRUE,
                                     output = NULL) {
  cfg <- load_v2_config(config)
  if (!is.null(input)) cfg$input$path <- input
  if (is.null(cfg$analysis)) cfg$analysis <- list(output_dir = "results")
  cfg <- resolve_config(cfg)

  issues <- list()
  add_issue <- function(severity, code, field, message) {
    issues[[length(issues) + 1]] <<- data.frame(
      severity = severity, code = code, field = field, message = message,
      stringsAsFactors = FALSE)
  }

  paths <- cfg$input$vcf
  if (is.null(paths) || length(paths) == 0 || any(is.na(paths))) {
    add_issue("error", "E_INPUT_MISSING", "input.path", "No input file was provided.")
  } else {
    missing <- paths[!file.exists(paths)]
    if (length(missing) > 0)
      add_issue("error", "E_INPUT_NOT_FOUND", "input.path",
                paste("Input file(s) not found:", paste(missing, collapse = ", ")))
  }

  genome_build <- NA_character_
  if (length(issues) == 0) {
    genome_build <- tryCatch(resolve_genome_build(cfg, paths), error = function(e) {
      add_issue("error", "E_GENOME_BUILD", "input.genome_build", conditionMessage(e))
      NA_character_
    })
  }

  ing <- NULL
  variants <- standardized <- featured <- NULL
  if (length(issues) == 0) {
    ing <- tryCatch(read_variant_input(
      paths,
      format = cfg_get(cfg, c("input", "format"), "auto"),
      sample_metadata = cfg_get(cfg, c("input", "sample_metadata"), NULL),
      column_map = cfg_get(cfg, c("input", "column_map"), NULL),
      chunk_size = cfg_get(cfg, c("input", "chunk_size"), NULL)
    ), error = function(e) {
      add_issue("error", "E_INPUT_PARSE", "input.path", conditionMessage(e))
      NULL
    })
  }

  if (!is.null(ing)) {
    variants <- ing$variants
    standardized <- tryCatch(standardize_variant_table(variants, cfg), error = function(e) {
      add_issue("error", "E_STANDARDIZE", "input.schema", conditionMessage(e))
      NULL
    })
    if (!is.null(standardized)) {
      featured <- tryCatch(add_basic_features(standardized, cfg), error = function(e) {
        add_issue("error", "E_FEATURES", "input.metrics", conditionMessage(e))
        NULL
      })
    }
  }

  if (!is.null(featured)) {
    if (all(is.na(featured$sample_id) | featured$sample_id == ""))
      add_issue("error", "E_SAMPLE_ID_MISSING", "SAMPLE_ID",
                "No usable sample identifier was detected.")
    if (isTRUE(strict) && all(is.na(featured$dp)))
      add_issue("error", "E_DEPTH_MISSING", "DP",
                "No usable depth evidence was detected or derivable.")
    if (isTRUE(strict) && all(is.na(featured$alt_count)))
      add_issue("error", "E_ALT_COUNT_MISSING", "AD_ALT",
                "No usable alternate-read evidence was detected.")
    if (isTRUE(strict) && all(is.na(featured$vaf)))
      add_issue("error", "E_VAF_MISSING", "VAF",
                "No usable VAF evidence was detected or derivable.")

    if (length(unique(stats::na.omit(featured$sample_id))) <
        cfg_get(cfg, c("cohort", "min_cohort_size_for_recurrence"), 10)) {
      add_issue("warning", "W_COHORT_RECURRENCE_NOT_EVALUABLE", "sample_id",
                "Cohort recurrence evidence is not evaluable at this sample count.")
    }
    if (all(is.na(featured$tlod)))
      add_issue("warning", "W_TLOD_MISSING", "TLOD",
                "TLOD is absent; caller-specific quality evidence is incomplete.")
    if (all(is.na(featured$max_pop_af)))
      add_issue("warning", "W_POPULATION_AF_MISSING", "population_frequency",
                "No configured population frequency evidence was detected.")
    if (all(is.na(featured$gene) | featured$gene == ""))
      add_issue("warning", "W_GENE_MISSING", "GENE",
                "No gene annotation was detected.")
  }

  schema <- if (!is.null(ing)) ing$schema_report else data.frame()
  summary <- validation_summary_table(ing, standardized, featured, genome_build)
  issue_df <- if (length(issues)) do.call(rbind, issues) else data.frame(
    severity = character(), code = character(), field = character(),
    message = character(), stringsAsFactors = FALSE)
  ok <- !any(issue_df$severity == "error")

  result <- list(
    ok = ok,
    strict = isTRUE(strict),
    input = paths,
    format = if (!is.null(ing)) ing$format else NA_character_,
    genome_build = genome_build,
    summary = summary,
    issues = issue_df,
    schema_report = schema,
    column_mapping = if (!is.null(ing)) ing$column_mapping else NULL
  )
  class(result) <- c("tumoronly_validation", class(result))
  write_validation_output(result, output)
  result
}

#' Run tumoronly with v2-style input/output overrides.
#' @export
run_tumoronly <- function(input = NULL, config, output_dir = NULL, run_id = NULL,
                          strict = TRUE, validate = TRUE) {
  cfg <- load_v2_config(config)
  if (!is.null(input)) cfg$input$path <- input
  if (!is.null(output_dir)) {
    if (is.null(cfg$analysis)) cfg$analysis <- list()
    cfg$analysis$output_dir <- output_dir
  }
  if (!is.null(run_id)) {
    if (is.null(cfg$analysis)) cfg$analysis <- list()
    cfg$analysis$run_id <- run_id
  }
  validation <- NULL
  if (isTRUE(validate)) {
    validation <- validate_tumoronly_input(config = cfg, strict = strict)
    if (isTRUE(strict) && !isTRUE(validation$ok)) {
      stop("Input validation failed:\n",
           paste(sprintf("- %s: %s", validation$issues$code, validation$issues$message),
                 collapse = "\n"), call. = FALSE)
    }
  }
  run <- run_tumor_only(cfg, run_id = run_id)
  run$validation <- validation
  run
}

load_v2_config <- function(config = NULL) {
  if (is.null(config)) return(tumoronly_default_config())
  if (is.character(config)) return(read_config(config))
  config
}

validation_summary_table <- function(ing, standardized, featured, genome_build) {
  variants <- if (!is.null(ing)) ing$variants else NULL
  n <- if (!is.null(variants)) nrow(variants) else 0L
  sample_count <- if (!is.null(standardized) && "sample_id" %in% names(standardized))
    length(unique(stats::na.omit(standardized$sample_id))) else NA_integer_
  variant_count <- if (!is.null(standardized) && "variant_id" %in% names(standardized))
    length(unique(standardized$variant_id)) else NA_integer_
  duplicate_sample_variant <- if (!is.null(standardized) && "sample_variant_key" %in% names(standardized))
    nrow(standardized) - length(unique(stats::na.omit(standardized$sample_variant_key))) else NA_integer_
  data.frame(
    metric = c("rows", "samples", "unique_variants", "duplicate_sample_variants",
               "genome_build", "format", "depth_available", "vaf_available"),
    value = c(
      as.character(n),
      as.character(sample_count),
      as.character(variant_count),
      as.character(duplicate_sample_variant),
      as.character(genome_build),
      if (!is.null(ing)) ing$format else NA_character_,
      as.character(!is.null(featured) && any(!is.na(featured$dp))),
      as.character(!is.null(featured) && any(!is.na(featured$vaf)))
    ),
    stringsAsFactors = FALSE
  )
}

write_validation_output <- function(result, output = NULL) {
  if (is.null(output) || is.na(output) || !nzchar(output)) return(invisible(NULL))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(output))
  if (ext == "json") {
    jsonlite::write_json(result, output, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else if (ext %in% c("tsv", "txt")) {
    write_tsv(result$issues, output)
  } else if (ext %in% c("html", "htm")) {
    write_validation_html(result, output)
  } else {
    stop("Unsupported validation output extension: ", ext, call. = FALSE)
  }
  invisible(output)
}

write_validation_html <- function(result, output) {
  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", as.character(x)))
  table_html <- function(df) {
    if (is.null(df) || nrow(df) == 0) return("<p><em>none</em></p>")
    hdr <- paste0("<tr>", paste(sprintf("<th>%s</th>", esc(names(df))), collapse = ""), "</tr>")
    rows <- apply(df, 1, function(r)
      paste0("<tr>", paste(sprintf("<td>%s</td>", esc(r)), collapse = ""), "</tr>"))
    paste0("<table>", hdr, paste(rows, collapse = "\n"), "</table>")
  }
  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'><title>tumoronly input validation</title>",
    "<style>body{font-family:Arial,sans-serif;margin:2rem;max-width:1100px}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 6px;font-size:12px}th{background:#eee}</style>",
    "</head><body>",
    sprintf("<h1>tumoronly input validation: %s</h1>", if (isTRUE(result$ok)) "PASS" else "FAIL"),
    "<h2>Summary</h2>", table_html(result$summary),
    "<h2>Issues</h2>", table_html(result$issues),
    "<h2>Schema</h2>", table_html(result$schema_report),
    "</body></html>")
  writeLines(html, output)
}
