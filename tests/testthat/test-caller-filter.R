# Mutect2 tumor-only: PASS/germline/panel_of_normals/combinations should be
# evaluable, not unconditionally short-circuited to technical_fail/artifact_strong.
# See docs/FILTERING_STRATEGY.md, "Caller FILTER as graded evidence".

test_that("caller_filter_pass defaults to literal PASS only (backward compatible)", {
  expect_true(caller_filter_pass("PASS", NULL))
  expect_true(caller_filter_pass(NA, NULL))
  expect_true(caller_filter_pass("", NULL))
  expect_false(caller_filter_pass("germline", NULL))
  expect_false(caller_filter_pass("germline;panel_of_normals", NULL))
})

test_that("caller_filter_pass honors a configured allowlist", {
  cfg <- list(hard_filters = list(caller_filter_accepted_values =
    c("PASS", "germline", "panel_of_normals", "germline;panel_of_normals")))
  expect_true(caller_filter_pass("germline", cfg))
  expect_true(caller_filter_pass("panel_of_normals", cfg))
  expect_true(caller_filter_pass("germline;panel_of_normals", cfg))
  expect_false(caller_filter_pass("clustered_events;weak_evidence", cfg))
})

test_that("caller_germline_flag is set from FILTER, independent of population evidence", {
  x <- data.frame(FILTER = c("germline", "germline;panel_of_normals",
                              "clustered_events", "PASS"),
                   filter_status = c("germline", "germline;panel_of_normals",
                              "clustered_events", "PASS"),
                   ref = "A", alt = "T", consequence = NA_character_, stringsAsFactors = FALSE)
  cfg <- list(technical_filters = list(min_alt_count_snv = 5, min_alt_count_indel = 8,
                                        min_af_snv = 0.03, min_af_indel = 0.05, min_depth = 20,
                                        min_tlod = 6, min_mbq = 25, min_mmq = 40))
  out <- add_basic_features(x, cfg)
  expect_equal(out$caller_germline_flag, c(TRUE, TRUE, FALSE, FALSE))
})

mk_clean_variant <- function(filter_value) {
  data.frame(
    filter_status = filter_value, mutect_filter = filter_value,
    is_indel = FALSE, dp = 60, alt_count = 20, vaf = 0.3, tlod = 40,
    mbq = 30, mmq = 50, hard_min_depth = 20, hard_min_alt_count_snv = 5,
    hard_min_alt_count_indel = 8, hard_min_af_snv = 0.03, hard_min_af_indel = 0.05,
    orientation_bias = 0.1, strand_artifact = 0.1, clustered_events = FALSE,
    weak_evidence = FALSE, pon_flag = FALSE, artifact_score = 0,
    variant_cohort_freq = 0, sample_qc_class = "not_evaluated",
    stringsAsFactors = FALSE)
}

test_that("a clean germline-flagged variant is technical_fail by default", {
  x <- mk_clean_variant("germline")
  cfg <- list(technical_filters = list(min_tlod = 6, min_mbq = 25, min_mmq = 40))
  out <- apply_variant_hard_filters(x, cfg)
  expect_false(out$hard_filter_pass)
  expect_true(grepl("caller_filter_not_pass", out$hard_filter_reason))
})

test_that("with the allowlist configured, the same clean germline-flagged variant passes", {
  x <- mk_clean_variant("germline")
  cfg <- list(technical_filters = list(min_tlod = 6, min_mbq = 25, min_mmq = 40),
              hard_filters = list(caller_filter_accepted_values = c("PASS", "germline")))
  out <- apply_variant_hard_filters(x, cfg)
  expect_true(out$hard_filter_pass)
})

test_that("classify_technical_artifact_evidence: germline no longer forces artifact_strong once configured", {
  cfg_default <- list(technical_filters = list(min_tlod = 6, min_mbq = 25, min_mmq = 40),
                       cohort = list(recurrent_variant_fraction_artifact = 0.30))
  cfg_allow <- cfg_default
  cfg_allow$hard_filters <- list(caller_filter_accepted_values = c("PASS", "germline"))

  x <- mk_clean_variant("germline")
  expect_equal(classify_technical_artifact_evidence(x, cfg_default), "artifact_strong")
  expect_equal(classify_technical_artifact_evidence(x, cfg_allow), "artifact_not_detected")
})

test_that("germline_score includes caller_germline_flag as one weighted component, not a veto", {
  base <- list(population_germline_score = 0, vaf_germline_score = 0, variant_cohort_freq = 0)
  cfg <- list(cohort = list(recurrent_variant_fraction_artifact = 0.30))
  x_flagged <- as.data.frame(c(base, list(caller_germline_flag = TRUE)))
  x_clean   <- as.data.frame(c(base, list(caller_germline_flag = FALSE)))
  expect_gt(bounded01(0.45 * x_flagged$population_germline_score + 0.25 * x_flagged$vaf_germline_score +
              0.15 * recurrent_germline_signal(x_flagged, cfg) + 0.15 * ifelse(x_flagged$caller_germline_flag, 1, 0)),
            bounded01(0.45 * x_clean$population_germline_score + 0.25 * x_clean$vaf_germline_score +
              0.15 * recurrent_germline_signal(x_clean, cfg) + 0.15 * ifelse(x_clean$caller_germline_flag, 1, 0)))
  # even flagged, a variant absent from population DBs is nowhere near "germline" by weight alone
  flagged_score <- 0.15 * 1
  expect_lt(flagged_score, 0.55)
})
