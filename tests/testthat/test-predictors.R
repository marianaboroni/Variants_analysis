mk <- function(consequence, ...) data.frame(consequence = consequence, ..., stringsAsFactors = FALSE)

test_that("SIFT 'deleterious(0.01)' is split into label + score", {
  ev <- extract_pred_value(c("deleterious(0.01)", "tolerated(0.4)", "0.998", "benign"))
  expect_equal(ev$label, c("deleterious", "tolerated", NA, "benign"))
  expect_equal(ev$score[1], 0.01); expect_equal(ev$score[3], 0.998)
})

test_that("SIFT + PolyPhen are one family (correlated): concordant damaging = weak support", {
  x <- mk("missense_variant", SIFT = "deleterious(0.01)", PolyPhen = "probably_damaging(0.99)")
  out <- add_predictor_evidence(x)
  # both protein_function -> 1 family -> not counted as 2 independent families
  expect_equal(out$COMPUTATIONAL_EVIDENCE_CATEGORY, "weak_support")
  expect_equal(out$SIFT_PRED, "deleterious"); expect_equal(out$SIFT_SCORE, 0.01)
  expect_false(out$COMPUTATIONAL_CONFLICT_FLAG)
})

test_that("two independent families concordant -> strong/moderate support", {
  x <- mk("missense_variant", SIFT = "deleterious(0.01)", REVEL = 0.9)  # protein_function + ensemble_missense
  out <- add_predictor_evidence(x)
  expect_true(out$COMPUTATIONAL_EVIDENCE_CATEGORY %in% c("strong_support", "moderate_support"))
})

test_that("SIFT damaging + PolyPhen benign -> conflicting", {
  x <- mk("missense_variant", SIFT = "deleterious(0.01)", PolyPhen = "benign(0.02)")
  out <- add_predictor_evidence(x)
  expect_true(out$COMPUTATIONAL_CONFLICT_FLAG)
  expect_equal(out$COMPUTATIONAL_EVIDENCE_CATEGORY, "conflicting")
})

test_that("SpliceAI high on a splice variant is applicable and damaging", {
  x <- mk("splice_donor_variant", SpliceAI_pred_DS_AG = 0.9, SpliceAI_pred_DS_DL = 0.1)
  out <- add_predictor_evidence(x)
  expect_equal(out$PREDICTOR_APPLICABILITY_STATUS, "applicable")
  expect_gt(out$COMPUTATIONAL_EVIDENCE_SCORE, 0)
})

test_that("missense predictors are NOT applicable to a frameshift", {
  x <- mk("frameshift_variant", SIFT = "deleterious(0.01)")
  out <- add_predictor_evidence(x)
  expect_equal(out$PREDICTOR_APPLICABILITY_STATUS, "not_applicable")
})

test_that("all predictors missing -> applicable_but_missing, not benign", {
  x <- mk("missense_variant")
  out <- add_predictor_evidence(x)
  expect_equal(out$PREDICTOR_APPLICABILITY_STATUS, "applicable_but_missing")
  expect_equal(out$COMPUTATIONAL_PREDICTORS_DAMAGING, 0L)
  expect_equal(out$COMPUTATIONAL_EVIDENCE_CATEGORY, "insufficient")
})

test_that("out-of-range predictor scores are dropped (not trusted)", {
  x <- mk("missense_variant", REVEL = 5)  # REVEL range 0..1
  out <- add_predictor_evidence(x)
  expect_true(is.na(out$REVEL_SCORE))
})

test_that("adding predictors does NOT change filter_status or final_class", {
  base <- data.frame(final_class = c("high_confidence_somatic", "likely_germline"),
    primary_reason = "x", hard_filter_reason = "PASS",
    population_category = c("population_rare_or_absent", "population_common"),
    artifact_category = "artifact_not_detected", recurrence_category = "recurrence_non_informative",
    somatic_score = c(0.9, 0.1), consequence = c("missense_variant", "missense_variant"),
    SIFT = c("deleterious(0.01)", "tolerated(0.5)"), stringsAsFactors = FALSE)
  before <- add_decision_trail(base)$filter_status
  withpred <- add_decision_trail(add_predictor_evidence(base))$filter_status
  expect_identical(before, withpred)
})
