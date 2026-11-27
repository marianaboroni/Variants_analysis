module_variants <- function() {
  data.frame(
    sample_id = c("S1", "S1", "S2"),
    tumor_type = c("UCEC", "UCEC", "BRCA"),
    final_class = c("high_confidence_somatic", "likely_germline", "probable_somatic"),
    consequence = c("missense_variant", "missense_variant", "synonymous_variant"),
    vaf = c(0.40, 0.50, 0.10),
    dp = c(100, 80, 60),
    alt_count = c(40, 40, 6),
    tumor_purity = c(0.80, 0.80, NA),
    driver_class = c("known_driver", "not_evaluable_as_driver", "likely_passenger"),
    stringsAsFactors = FALSE
  )
}

test_that("TMB requires explicit callable territory but always reports burden counts", {
  v <- module_variants()
  cfg <- tumoronly_default_config()
  cfg$tmb$callable_mb <- NULL
  cfg$tmb$require_callable_mb <- TRUE
  res <- calculate_tmb(v, cfg)
  expect_true("tmb_countable" %in% names(res$variants))
  expect_true(any(res$summary$module_status == "not_evaluable_missing_callable_mb"))
  expect_true(all(is.na(res$summary$tmb_mut_per_mb)))
  expect_equal(res$summary$n_tmb_countable[res$summary$sample_id == "S1"], 1)

  cfg$tmb$callable_mb <- 30
  res2 <- calculate_tmb(v, cfg)
  expect_equal(res2$summary$tmb_mut_per_mb[res2$summary$sample_id == "S1"], 1 / 30)
  expect_equal(res2$summary$burden_label[res2$summary$sample_id == "S1"], "tumor_mutational_burden")
})

test_that("clonality uses CCF when purity is present and labels VAF proxy when missing", {
  v <- module_variants()
  cfg <- tumoronly_default_config()
  res <- add_clonality_estimates(v, cfg)
  expect_true(all(c("ccf_estimate", "clonality_class", "clonality_method") %in% names(res$variants)))
  expect_equal(res$variants$clonality_method[1], "purity_adjusted_copy_neutral")
  expect_false(is.na(res$variants$ccf_estimate[1]))
  expect_equal(res$variants$clonality_method[3], "vaf_proxy_no_purity")
  expect_match(res$summary$limitation[res$summary$sample_id == "S2"], "VAF-proxy")
})

write_purple_cnv <- function(path) {
  writeLines(c(
    paste("chromosome", "start", "end", "copyNumber", "bafCount", "observedBAF", "baf",
          "minorAlleleCopyNumber", "majorAlleleCopyNumber", sep = "\t"),
    paste("chr1", 1, 1000, 4.02, 10, 0.5, 0.5, 2.01, 2.01, sep = "\t"),
    paste("chr1", 1001, 5000, 2.95, 10, 0.66, 0.66, 1.0, 1.95, sep = "\t"),
    paste("chr2", 1, 5000, -0.01, 0, 0, 0, 0, 0, sep = "\t")
  ), path)
}

cn_module_variants <- function(cn_file) {
  data.frame(
    sample_id = c("S1", "S1", "S1", "S1", "S1"),
    tumor_type = "HGSOC",
    final_class = c("high_confidence_somatic", "probable_somatic", "likely_germline",
                    "high_confidence_somatic", "high_confidence_somatic"),
    chrom = c("1", "chr1", "1", "1", "2"),
    pos = c(500, 2000, 600, 9000, 100),
    vaf = c(0.40, 0.30, 0.50, 0.35, 0.20),
    tumor_purity = 0.80,
    copy_number_file = cn_file,
    stringsAsFactors = FALSE
  )
}

test_that("clonality takes total copy number from a PURPLE segment file", {
  cn_file <- tempfile(fileext = ".purple.cnv.somatic.tsv")
  write_purple_cnv(cn_file)
  v <- cn_module_variants(cn_file)
  res <- add_clonality_estimates(v, tumoronly_default_config())$variants

  # somatic variants inside segments use the segment CN (chr prefix ignored)
  expect_equal(res$total_cn_used[1:2], c(4.02, 2.95))
  expect_equal(res$minor_cn_used[1:2], c(2.01, 1.0))
  expect_equal(res$major_cn_used[1:2], c(2.01, 1.95))
  expect_equal(res$copy_number_source[1:2], c("segment_file", "segment_file"))
  expect_equal(res$clonality_method[1:2], rep("purity_copy_number_adjusted", 2))
  expect_equal(res$ccf_estimate[1], 0.40 * (0.8 * 4.02 + 0.2 * 2) / 0.8)

  # germline rows are not looked up; uncovered positions and non-positive CN
  # segments fall back to the copy-neutral default
  expect_equal(res$copy_number_source[3:5], rep("default_copy_neutral", 3))
  expect_equal(res$total_cn_used[4:5], c(2, 2))
  expect_equal(res$clonality_method[4], "purity_adjusted_copy_neutral")
})

test_that("per-variant total_cn takes precedence over the segment file", {
  cn_file <- tempfile(fileext = ".tsv")
  write_purple_cnv(cn_file)
  v <- cn_module_variants(cn_file)
  v$total_cn <- c(3, NA, NA, NA, NA)
  res <- add_clonality_estimates(v, tumoronly_default_config())$variants
  expect_equal(res$total_cn_used[1:2], c(3, 2.95))
  expect_equal(res$copy_number_source[1:2], c("variant_column", "segment_file"))
})

test_that("missing copy-number file degrades to copy-neutral without failing", {
  v <- cn_module_variants(file.path(tempdir(), "does_not_exist.cnv.tsv"))
  res <- add_clonality_estimates(v, tumoronly_default_config())
  expect_equal(unique(res$variants$copy_number_source), "default_copy_neutral")
  expect_equal(res$variants$clonality_method[1], "purity_adjusted_copy_neutral")
  expect_match(res$summary$copy_number_sources, "default_copy_neutral")
})

test_that("sample metadata carries copy_number_file, resolving relative paths", {
  dir <- tempfile("meta_dir"); dir.create(dir)
  meta <- file.path(dir, "meta.tsv")
  writeLines(c("sample_id\ttumor_purity\tcopy_number_file",
               "S1\t0.62\tS1.purple.cnv.somatic.tsv",
               "S2\t0.40\t/abs/S2.purple.cnv.somatic.tsv"), meta)
  cfg <- tumoronly_default_config()
  cfg$input$sample_metadata <- meta
  x <- data.frame(sample_id = c("S1", "S2", "S3"), tumor_type = "HGSOC", stringsAsFactors = FALSE)
  out <- attach_sample_metadata(x, cfg)
  out <- out[order(out$sample_id), ]
  expect_equal(out$copy_number_file[1], file.path(dir, "S1.purple.cnv.somatic.tsv"))
  expect_equal(out$copy_number_file[2], "/abs/S2.purple.cnv.somatic.tsv")
  expect_true(is.na(out$copy_number_file[3]))
  expect_equal(out$tumor_purity[1:2], c(0.62, 0.40))
})

test_that("ML module degrades explicitly when no active model is available", {
  v <- module_variants()
  cfg <- tumoronly_default_config()
  cfg$ml$db_dir <- tempfile("missing_ml_db")
  res <- apply_ml_predictions(v, cfg)
  expect_true("ml_status" %in% names(res$variants))
  expect_equal(unique(res$variants$ml_status), "not_evaluable_no_active_model")
  expect_equal(res$status$status, "not_evaluable_no_active_model")
})
