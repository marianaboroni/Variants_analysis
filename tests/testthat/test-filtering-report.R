test_that("decision trail assigns PASS/REVIEW/FAIL and never drops variants", {
  x <- data.frame(
    final_class = c("high_confidence_somatic", "manual_review_required",
                    "likely_germline", "technical_fail", "probable_somatic"),
    primary_reason = c("a", "b", "c", "d", "e"),
    hard_filter_reason = c("PASS", "PASS", "PASS", "low_depth", "PASS"),
    population_category = c("population_rare_or_absent", "population_uninformative",
                            "population_common", "population_uninformative", "population_rare_or_absent"),
    artifact_category = rep("artifact_not_detected", 5),
    recurrence_category = rep("recurrence_non_informative", 5),
    somatic_score = c(0.9, 0.5, 0.1, 0.0, 0.6),               # base (no COSMIC context)
    somatic_score_validated = c(0.95, 0.5, 0.1, 0.0, 0.65),   # experimental (with COSMIC)
    stringsAsFactors = FALSE)
  out <- add_decision_trail(x)
  expect_equal(out$filter_status, c("PASS", "REVIEW", "FAIL", "FAIL", "PASS"))
  expect_equal(nrow(out), 5)                        # nothing removed
  # main confidence == BASE (no COSMIC context)
  expect_equal(out$confidence_category[1], "high")
  expect_equal(out$CONFIDENCE_CATEGORY_BASE[1], "high")
  expect_equal(out$confidence_score[1], 0.9)
  expect_equal(out$CONFIDENCE_SCORE_WITH_COSMIC_EXPERIMENTAL[1], 0.95)
  expect_equal(out$confidence_category[3], "excluded")
  sp <- split_by_status(out)
  expect_equal(nrow(sp$retained) + nrow(sp$excluded), nrow(out))
})

test_that("report renders (or falls back) even with no retained variants", {
  run_dir <- tempfile("emptyrun"); dir.create(file.path(run_dir, "tables"), recursive = TRUE)
  empty <- data.frame(sample_id = character(), gene = character(),
                      filter_status = character(), stringsAsFactors = FALSE)
  write_tsv(empty, file.path(run_dir, "tables", "variants_retained.tsv.gz"))
  write_tsv(data.frame(filter = "final_PASS", n_variants = 0),
            file.path(run_dir, "tables", "filter_audit.tsv.gz"))
  write_tsv(data.frame(note = "no variants"), file.path(run_dir, "tables", "variants_all.tsv.gz"))
  jsonlite::write_json(list(run_id = "empty", genome_build = "GRCh38"),
                       file.path(run_dir, "manifest.json"), auto_unbox = TRUE)
  html <- render_fallback_report(run_dir, file.path(run_dir, "report.html"))
  expect_true(file.exists(html))
  expect_true(any(grepl("No retained variants", readLines(html))))
})

test_that("oncoplot degrades gracefully with no variants", {
  skip_if_not_installed("maftools")
  empty_maf <- tempfile(fileext = ".maf")
  write_tsv(data.frame(Hugo_Symbol = character(), Tumor_Sample_Barcode = character()), empty_maf)
  res <- make_oncoplot(empty_maf, tempfile("plots"))
  expect_true(res$status %in% c("no_variants", "no_genes", "read_maf_failed"))
})
