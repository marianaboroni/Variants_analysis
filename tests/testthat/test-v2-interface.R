minimal_v2_tsv <- function() {
  path <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(
    Sample_Barcode = c("S1", "S1", "S1"),
    CHROM = c("chr1", "chr1", "chr2"),
    START = c(101, 202, 303),
    REF = c("A", "G", "C"),
    ALT = c("T", "A", "CT"),
    FILTER = "PASS",
    ref_count = c(40, 30, 60),
    alt_count = c(20, 12, 15),
    AF = c(0.333, 0.286, 0.200),
    MBQ = c(35, 34, 33),
    MMQ = c(60, 55, 50),
    SYMBOL = c("TP53", "KRAS", "BRAF"),
    Consequence = c("missense_variant", "synonymous_variant", "inframe_insertion"),
    stringsAsFactors = FALSE
  ), path, sep = "\t")
  path
}

minimal_v2_config <- function(input = NULL, output_dir = tempfile("tumoronly_results")) {
  cfg <- tumoronly_default_config()
  cfg$input$path <- input
  cfg$input$genome_build <- "GRCh38"
  cfg$analysis$output_dir <- output_dir
  cfg$analysis$run_id <- "v2_test"
  cfg$cosmic$processed_db <- NULL
  cfg$cosmic$release <- NULL
  cfg$driver_resources$driver_genes <- NULL
  cfg$driver_resources$hotspots <- NULL
  cfg$ancestry$enabled <- FALSE
  cfg$plots$minimum_mutations_for_rainfall <- 1000
  cfg
}

test_that("write_tumoronly_config_template writes and refuses accidental overwrite", {
  cfg_file <- tempfile(fileext = ".yml")
  expect_false(file.exists(cfg_file))
  write_tumoronly_config_template(cfg_file)
  expect_true(file.exists(cfg_file))
  cfg <- yaml::read_yaml(cfg_file)
  expect_true("technical_filters" %in% names(cfg))
  expect_error(write_tumoronly_config_template(cfg_file), "already exists")
  expect_silent(write_tumoronly_config_template(cfg_file, force = TRUE))
})

test_that("validate_tumoronly_input accepts a minimal real TSV contract", {
  input <- minimal_v2_tsv()
  cfg <- minimal_v2_config(input)
  res <- validate_tumoronly_input(config = cfg)
  expect_true(res$ok)
  expect_equal(res$format, "tsv")
  expect_equal(res$genome_build, "GRCh38")
  expect_true("W_TLOD_MISSING" %in% res$issues$code)
  expect_equal(res$summary$value[res$summary$metric == "samples"], "1")
  expect_equal(res$summary$value[res$summary$metric == "depth_available"], "TRUE")
})

test_that("validate_tumoronly_input blocks TSVs without explicit genome build", {
  input <- minimal_v2_tsv()
  cfg <- minimal_v2_config(input)
  cfg$input$genome_build <- NULL
  res <- validate_tumoronly_input(config = cfg)
  expect_false(res$ok)
  expect_true("E_GENOME_BUILD" %in% res$issues$code)
})

test_that("CLI validate shares the package validation implementation", {
  input <- minimal_v2_tsv()
  cfg_file <- tempfile(fileext = ".yml")
  yaml::write_yaml(minimal_v2_config(input), cfg_file)
  out_html <- tempfile(fileext = ".html")
  root <- get(".tumoronly_root", envir = globalenv())
  script <- file.path(root, "exec", "tumoronly")
  cmd <- c(script, "validate", "--input", input, "--config", cfg_file,
           "--output", out_html)
  res <- system2("Rscript", cmd, stdout = TRUE, stderr = TRUE)
  expect_null(attr(res, "status"))
  expect_true(any(grepl("Validation: PASS", res, fixed = TRUE)))
  expect_true(file.exists(out_html))
})

test_that("run_tumoronly executes the workflow through the v2 API wrapper", {
  skip_if_not_installed("maftools")
  input <- minimal_v2_tsv()
  cfg <- minimal_v2_config(input)
  cfg$analysis$run_id <- "api_wrapper"
  capture.output({
    res <- suppressWarnings(suppressMessages(run_tumoronly(config = cfg, strict = TRUE)))
  })
  expect_true(res$validation$ok)
  expect_true(file.exists(file.path(res$run_dir, "tables", "variants_all.tsv.gz")))
  expect_true(file.exists(file.path(res$run_dir, "report", "tumor_only_report.html")))
})

test_that("add_tumoronly_evidence_columns adds traceability without changing classes", {
  x <- data.frame(
    sample_id = "S1", final_class = "probable_somatic", filter_status = "PASS",
    primary_reason = "probable_somatic_rule", hard_filter_reason = "PASS",
    population_category = "population_rare_or_absent",
    artifact_category = "artifact_not_detected",
    recurrence_category = "recurrence_non_informative",
    somatic_score = 0.72, technical_evidence_score = 0.8,
    dp = 100, alt_count = 35, vaf = 0.35, tlod = NA_real_,
    mbq = 35, mmq = 60, max_pop_af = NA_real_, gene = "TP53",
    stringsAsFactors = FALSE
  )
  cfg <- tumoronly_default_config()
  out <- add_tumoronly_evidence_columns(x, cfg)
  expect_equal(out$final_class, x$final_class)
  expect_true(all(c("evidence_supporting_classification",
                    "evidence_against_classification",
                    "missing_evidence",
                    "classification_explanation") %in% names(out)))
  expect_match(out$evidence_supporting_classification, "somatic_score")
  expect_match(out$missing_evidence, "matched_normal_absent")
  expect_match(out$classification_explanation, "probabilistic")
})

test_that("run_tumoronly writes the v2 output contract", {
  skip_if_not_installed("ggplot2")
  input <- minimal_v2_tsv()
  cfg <- minimal_v2_config(input)
  cfg$analysis$run_id <- "v2_output_contract"
  cfg$analysis$output_dir <- tempfile("tumoronly_v2_contract")
  capture.output({
    res <- suppressWarnings(suppressMessages(run_tumoronly(config = cfg, strict = TRUE)))
  })

  expected <- c(
    "report.html", "run_manifest.json", "config_used.yaml", "session_info.txt",
    file.path("tables", "all_variants.tsv"),
    file.path("tables", "classified_variants.tsv"),
    file.path("tables", "high_confidence_somatic.tsv"),
    file.path("tables", "likely_germline.tsv"),
    file.path("tables", "likely_artifact.tsv"),
    file.path("tables", "known_drivers.tsv"),
    file.path("tables", "sample_summary.tsv"),
    file.path("tables", "filter_audit.tsv"),
    file.path("figures", "figure_manifest.tsv"),
    file.path("figure_data", "figure_01_filtering_workflow.tsv"),
    file.path("logs", "warnings.tsv")
  )
  expect_true(all(file.exists(file.path(res$run_dir, expected))))

  classified <- read_variants(file.path(res$run_dir, "tables", "classified_variants.tsv"), "\t")
  expect_true(all(c("evidence_supporting_classification",
                    "evidence_against_classification",
                    "missing_evidence",
                    "classification_explanation") %in% names(classified)))
  expect_equal(nrow(classified), nrow(res$variants))

  manifest <- jsonlite::read_json(file.path(res$run_dir, "run_manifest.json"),
                                  simplifyVector = TRUE)
  expect_equal(manifest$tool$name, "tumoronly")
  expect_true(manifest$run$status %in% c("completed", "completed_with_warnings"))
  expect_equal(manifest$input$genome_build, "GRCh38")

  figman <- read_variants(file.path(res$run_dir, "figures", "figure_manifest.tsv"), "\t")
  ok <- figman[figman$figure_id == "figure_03_classification", , drop = FALSE]
  expect_equal(ok$status, "ok")
  expect_match(ok$files, "[.]pdf")
  expect_match(ok$files, "[.]svg")
  expect_match(ok$files, "[.]png")
})
