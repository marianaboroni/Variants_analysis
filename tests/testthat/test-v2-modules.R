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

test_that("ML module degrades explicitly when no active model is available", {
  v <- module_variants()
  cfg <- tumoronly_default_config()
  cfg$ml$db_dir <- tempfile("missing_ml_db")
  res <- apply_ml_predictions(v, cfg)
  expect_true("ml_status" %in% names(res$variants))
  expect_equal(unique(res$variants$ml_status), "not_evaluable_no_active_model")
  expect_equal(res$status$status, "not_evaluable_no_active_model")
})
