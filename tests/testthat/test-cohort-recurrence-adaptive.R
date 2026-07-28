# Cohort-recurrence-based artifact evidence must adapt to cohort size: a
# fraction alone is meaningless at N=2 (1/2 samples sharing a variant is
# already 50%, trivially over any sane threshold), but is exactly the right
# signal at N=100. See docs/FILTERING_STRATEGY.md, "Cohort-size-adaptive
# recurrence", and cohort_recurrent_flag() (utils.R).

cfg_default <- list(cohort = list(
  recurrent_variant_fraction_artifact = 0.30,
  recurrent_locus_fraction_artifact = 0.40,
  min_cohort_size_for_recurrence = 10,
  min_recurrent_samples = 3,
  min_recurrent_locus_samples = 3
))

test_that("cohort_recurrent_flag: NA (not FALSE) when the cohort is too small to evaluate", {
  x <- data.frame(variant_n_samples = 5, variant_cohort_freq = 0.9, total_samples = 2)
  expect_true(is.na(cohort_recurrent_flag(x, cfg_default, "variant")))
})

test_that("cohort_recurrent_flag: fraction met but below the absolute sample-count floor -> FALSE", {
  x <- data.frame(variant_n_samples = 2, variant_cohort_freq = 0.30, total_samples = 10)
  expect_false(isTRUE(cohort_recurrent_flag(x, cfg_default, "variant")))
})

test_that("cohort_recurrent_flag: fraction met and floor met, cohort large enough -> TRUE", {
  x <- data.frame(variant_n_samples = 3, variant_cohort_freq = 0.30, total_samples = 10)
  expect_true(cohort_recurrent_flag(x, cfg_default, "variant"))
})

test_that("cohort_recurrent_flag: missing recurrence columns are treated as not-evaluable, not zero", {
  x <- data.frame(chrom = "1", pos = 100)  # no variant_n_samples/variant_cohort_freq/total_samples at all
  flag <- cohort_recurrent_flag(x, cfg_default, "variant")
  expect_true(is.na(flag))
})

test_that("cohort_recurrent_flag: locus level uses its own floor/threshold keys", {
  x <- data.frame(locus_n_samples = 5, locus_cohort_freq = 0.50, total_samples = 10)
  expect_true(cohort_recurrent_flag(x, cfg_default, "locus"))
  x2 <- data.frame(locus_n_samples = 5, locus_cohort_freq = 0.35, total_samples = 10)  # below 0.40 locus threshold
  expect_false(isTRUE(cohort_recurrent_flag(x2, cfg_default, "locus")))
})

test_that("cohort_artifact_score ramps up FROM the threshold (rare=low, recurrent=high), not inverted", {
  x <- data.frame(variant_cohort_freq = c(0.05, 0.60), locus_cohort_freq = c(0.05, 0.05),
                   variant_n_samples = c(20, 20), locus_n_samples = c(20, 20), total_samples = c(100, 100),
                   variant_median_vaf = c(0.3, 0.3))
  s <- cohort_artifact_score(x, cfg_default)
  expect_lt(s[1], 0.1)   # rare variant: near-zero artifact evidence
  expect_gt(s[2], 0.9)   # recurrent variant: near-maximum artifact evidence
})

test_that("cohort_artifact_score contributes nothing when the cohort is too small, however high the freq", {
  x <- data.frame(variant_cohort_freq = 1.0, locus_cohort_freq = 1.0,
                   variant_n_samples = 2, locus_n_samples = 2, total_samples = 2,
                   variant_median_vaf = 0.3)
  expect_equal(cohort_artifact_score(x, cfg_default), 0)
})

test_that("cohort_artifact_score does not use low median VAF as recurrence evidence in a tiny cohort", {
  x <- data.frame(variant_cohort_freq = 1.0, locus_cohort_freq = 1.0,
                   variant_n_samples = 1, locus_n_samples = 1, total_samples = 1,
                   variant_median_vaf = 0.03)
  expect_equal(cohort_artifact_score(x, cfg_default), 0)
})

test_that("recurrent_germline_signal: N=2 contributes 0 regardless of freq; adequate N+count triggers it", {
  small <- data.frame(variant_n_samples = 1, variant_cohort_freq = 1.0, total_samples = 2)
  expect_equal(recurrent_germline_signal(small, cfg_default), 0)

  big <- data.frame(variant_n_samples = 5, variant_cohort_freq = 0.5, total_samples = 10)
  expect_equal(recurrent_germline_signal(big, cfg_default), 1)
})

test_that("classify_technical_artifact_evidence: a 2-sample cohort never reaches artifact_possible via recurrence alone", {
  x <- data.frame(
    filter_status = "PASS", mutect_filter = "PASS", is_indel = FALSE,
    dp = 60, alt_count = 20, vaf = 0.5, tlod = 40, mbq = 30, mmq = 50,
    hard_min_depth = 20, hard_min_alt_count_snv = 5, hard_min_alt_count_indel = 8,
    orientation_bias = 0.1, strand_artifact = 0.1, clustered_events = FALSE, weak_evidence = FALSE,
    pon_flag = FALSE, artifact_score = 0,
    variant_n_samples = 1, variant_cohort_freq = 1.0, total_samples = 2,
    stringsAsFactors = FALSE)
  cfg <- list(technical_filters = list(min_tlod = 6, min_mbq = 25, min_mmq = 40),
              cohort = cfg_default$cohort)
  expect_equal(classify_technical_artifact_evidence(x, cfg), "artifact_not_detected")
})

test_that("classify_technical_artifact_evidence: the same recurrence pattern at adequate cohort size IS flagged", {
  x <- data.frame(
    filter_status = "PASS", mutect_filter = "PASS", is_indel = FALSE,
    dp = 60, alt_count = 20, vaf = 0.5, tlod = 40, mbq = 30, mmq = 50,
    hard_min_depth = 20, hard_min_alt_count_snv = 5, hard_min_alt_count_indel = 8,
    orientation_bias = 0.1, strand_artifact = 0.1, clustered_events = FALSE, weak_evidence = FALSE,
    pon_flag = FALSE, artifact_score = 0,
    variant_n_samples = 5, variant_cohort_freq = 0.5, total_samples = 10,
    stringsAsFactors = FALSE)
  cfg <- list(technical_filters = list(min_tlod = 6, min_mbq = 25, min_mmq = 40),
              cohort = cfg_default$cohort)
  expect_equal(classify_technical_artifact_evidence(x, cfg), "artifact_possible")
})

test_that("classify_internal_recurrence: not-evaluable cohort -> recurrence_non_informative, not recurrence_artifact_suspected", {
  x <- data.frame(variant_cohort_freq = 1.0, locus_cohort_freq = 1.0,
                   variant_n_samples = 1, locus_n_samples = 1, total_samples = 2,
                   variant_tumor_type_freq = NA_real_)
  expect_equal(classify_internal_recurrence(x, cfg_default), "recurrence_non_informative")
})

test_that("classify_internal_recurrence: one sample is not tumor-type recurrence support", {
  x <- data.frame(variant_cohort_freq = 1.0, locus_cohort_freq = 1.0,
                   variant_n_samples = 1, locus_n_samples = 1, total_samples = 1,
                   variant_tumor_type_n_samples = 1, tumor_type_total_samples = 1,
                   variant_tumor_type_freq = 1.0)
  expect_equal(classify_internal_recurrence(x, cfg_default), "recurrence_non_informative")
})

test_that("classify_internal_recurrence: tumor-type support needs evaluable sample counts", {
  x <- data.frame(variant_cohort_freq = 0.10, locus_cohort_freq = 0.10,
                   variant_n_samples = 2, locus_n_samples = 2, total_samples = 20,
                   variant_tumor_type_n_samples = 3, tumor_type_total_samples = 20,
                   variant_tumor_type_freq = 0.15)
  expect_equal(classify_internal_recurrence(x, cfg_default), "recurrence_tumor_type_supported")

  x$variant_tumor_type_n_samples <- 2
  expect_equal(classify_internal_recurrence(x, cfg_default), "recurrence_non_informative")
})

test_that("classify_internal_recurrence: adequate cohort size correctly flags recurrence_artifact_suspected", {
  x <- data.frame(variant_cohort_freq = 0.5, locus_cohort_freq = 0.1,
                   variant_n_samples = 5, locus_n_samples = 1, total_samples = 10,
                   variant_tumor_type_freq = NA_real_)
  expect_equal(classify_internal_recurrence(x, cfg_default), "recurrence_artifact_suspected")
})
