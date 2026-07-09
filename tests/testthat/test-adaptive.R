mk_cohort <- function() {
  set.seed(1)
  hi <- data.frame(sample_id = "HIGH", dp = round(rnorm(40, 500, 40)),
    vaf = 0.4, alt_count = 200, orientation_bias = 0.1, strand_artifact = 0.1,
    pon_flag = FALSE, stringsAsFactors = FALSE)
  lo <- data.frame(sample_id = "LOW", dp = round(rnorm(40, 30, 5)),
    vaf = 0.4, alt_count = 12, orientation_bias = 0.1, strand_artifact = 0.1,
    pon_flag = FALSE, stringsAsFactors = FALSE)
  rbind(hi, lo)
}

test_that("adaptive min DP differs by sample and is never below the absolute floor", {
  x <- add_adaptive_filtering(mk_cohort(), list(adaptive_filtering = list(
    depth = list(absolute_floor = 8, method = "percentile_mad", percentile = 0.05, mad_multiplier = 3))))
  thr <- attr(x, "adaptive_thresholds")
  hi <- thr$ADAPTIVE_MIN_DP[thr$sample_id == "HIGH"]
  lo <- thr$ADAPTIVE_MIN_DP[thr$sample_id == "LOW"]
  expect_gt(hi, lo)                    # high-coverage sample gets a higher adaptive floor
  expect_true(all(thr$ADAPTIVE_MIN_DP >= 8))   # GUARDRAIL
})

test_that("Panel-of-Normals and strong bias force FAIL; deterministic", {
  base <- mk_cohort()
  base$pon_flag[1] <- TRUE
  base$strand_artifact[2] <- 0.95
  o1 <- add_adaptive_filtering(base, NULL)
  o2 <- add_adaptive_filtering(base, NULL)
  expect_equal(o1$ADAPTIVE_FILTER_STATUS[1], "FAIL")   # PoN
  expect_equal(o1$ADAPTIVE_FILTER_STATUS[2], "FAIL")   # strong bias
  expect_identical(o1$ADAPTIVE_FILTER_STATUS, o2$ADAPTIVE_FILTER_STATUS)  # deterministic
})

test_that("low alt-read count below absolute floor FAILs (AD_ALT=1 not rescued)", {
  x <- mk_cohort()
  x$alt_count[1] <- 1; x$dp[1] <- 6      # below absolute floors
  o <- add_adaptive_filtering(x, NULL)
  expect_equal(o$ADAPTIVE_FILTER_STATUS[1], "FAIL")
})

test_that("three zones are produced and PASS/REVIEW/FAIL are all reachable", {
  x <- mk_cohort()
  o <- add_adaptive_filtering(x, NULL)
  expect_true(all(o$ADAPTIVE_FILTER_STATUS %in% c("PASS", "REVIEW", "FAIL")))
  expect_true(all(c("SAMPLE_NOISE_FLOOR", "AD_ALT_SIGNAL_TO_NOISE", "ALT_READ_ERROR_PVALUE",
    "EXPECTED_MIN_DETECTABLE_VAF", "ADAPTIVE_MIN_DP") %in% names(o)))
})

test_that("globally poor sample is flagged (not rescued by lowering thresholds)", {
  bad <- data.frame(sample_id = "BAD", dp = round(rnorm(40, 6, 2)), vaf = 0.3,
    alt_count = 2, orientation_bias = 0.1, strand_artifact = 0.1, pon_flag = FALSE,
    stringsAsFactors = FALSE)
  thr <- attr(add_adaptive_filtering(bad, NULL), "adaptive_thresholds")
  expect_true(thr$SAMPLE_QC_STATUS %in% c("FAIL", "REVIEW"))
  expect_gte(thr$ADAPTIVE_MIN_DP, 8)   # still respects floor
})

test_that("adaptive REVIEW escalation is documented and does not alter FAIL/CONFIDENCE_SCORE_BASE", {
  x <- data.frame(final_class = "high_confidence_somatic", primary_reason = "x",
    hard_filter_reason = "PASS", population_category = "population_rare_or_absent",
    artifact_category = "artifact_not_detected", recurrence_category = "recurrence_non_informative",
    somatic_score = 0.9, sample_id = "S", dp = 12, vaf = 0.4, alt_count = 6,
    orientation_bias = 0.1, strand_artifact = 0.1, pon_flag = FALSE, stringsAsFactors = FALSE)
  x <- add_adaptive_filtering(x, NULL)
  base <- add_decision_trail(x)
  before_conf <- base$CONFIDENCE_SCORE_BASE
  lab <- add_unified_labels(base)
  esc <- apply_adaptive_review(lab, NULL)
  expect_identical(esc$CONFIDENCE_SCORE_BASE, before_conf)  # base confidence unchanged
})
