# Authenticity = TECHNICAL only; biology must never mask a technical artifact.

tech_only <- function(tech, artifact, ...) data.frame(
  technical_evidence_score = tech, artifact_score = artifact, ..., stringsAsFactors = FALSE)

test_that("technical authenticity ignores biology (same tech, different biology => same score)", {
  a <- tech_only(0.9, 0.1, hotspot_match = TRUE, driver_score = 0.9,
                 COSMIC_GLOBAL_RECURRENCE_SCORE = 1, impact_rank = 3)
  b <- tech_only(0.9, 0.1, hotspot_match = FALSE, driver_score = 0,
                 COSMIC_GLOBAL_RECURRENCE_SCORE = 0, impact_rank = 0)
  oa <- add_variant_authenticity(a); ob <- add_variant_authenticity(b)
  expect_equal(oa$TECHNICAL_AUTHENTICITY_SCORE, ob$TECHNICAL_AUTHENTICITY_SCORE)
  expect_gt(oa$BIOLOGICAL_SUPPORT_SCORE, ob$BIOLOGICAL_SUPPORT_SCORE)  # biology differs
})

test_that("hotspot with weak technical evidence is NOT masked -> not LIKELY_TRUE", {
  x <- tech_only(0.05, 0.85, hotspot_match = TRUE, driver_score = 0.95,
                 impact_rank = 3, COMPUTATIONAL_EVIDENCE_SCORE = 0.9)
  out <- add_variant_authenticity(x)
  expect_true(out$VARIANT_AUTHENTICITY_CATEGORY %in% c("LIKELY_ARTIFACT", "UNCERTAIN"))
  expect_gte(out$BIOLOGICAL_SUPPORT_SCORE, 0.5)          # biologically interesting
  expect_gt(out$VARIANT_REVIEW_SCORE, 0.5)               # surfaced for review
})

test_that("hotspot with strong strand bias stays LIKELY_ARTIFACT", {
  x <- tech_only(0.4, 0.2, hotspot_match = TRUE, strand_artifact = 0.95, driver_score = 0.9)
  out <- add_variant_authenticity(x)
  expect_equal(out$VARIANT_AUTHENTICITY_CATEGORY, "LIKELY_ARTIFACT")
})

test_that("technically excellent non-hotspot is LIKELY_TRUE", {
  x <- tech_only(0.95, 0.02, hotspot_match = FALSE, driver_score = 0)
  out <- add_variant_authenticity(x)
  expect_equal(out$VARIANT_AUTHENTICITY_CATEGORY, "LIKELY_TRUE")
})

test_that("Panel-of-Normals reinforces artifact evidence", {
  x <- tech_only(0.6, 0.1, pon_flag = TRUE)
  out <- add_variant_authenticity(x)
  expect_gte(out$ARTIFACT_EVIDENCE_SCORE, 0.6)
})

test_that("conflicting predictors are down-weighted vs concordant in biological support", {
  conc <- tech_only(0.5, 0.1, impact_rank = 2, COMPUTATIONAL_EVIDENCE_SCORE = 0.8,
                    COMPUTATIONAL_CONFLICT_FLAG = FALSE)
  conf <- tech_only(0.5, 0.1, impact_rank = 2, COMPUTATIONAL_EVIDENCE_SCORE = 0.8,
                    COMPUTATIONAL_CONFLICT_FLAG = TRUE)
  expect_gt(add_variant_authenticity(conc)$BIOLOGICAL_SUPPORT_SCORE,
            add_variant_authenticity(conf)$BIOLOGICAL_SUPPORT_SCORE)
})

test_that("robust to missing columns; deterministic", {
  x <- data.frame(consequence = "missense_variant", stringsAsFactors = FALSE)
  o1 <- add_variant_authenticity(x); o2 <- add_variant_authenticity(x)
  expect_false(is.na(o1$TECHNICAL_AUTHENTICITY_SCORE))
  expect_identical(o1$VARIANT_AUTHENTICITY_SCORE, o2$VARIANT_AUTHENTICITY_SCORE)
})

test_that("authenticity does not change filter_status / filter_reasons / CONFIDENCE_SCORE_BASE", {
  x <- data.frame(final_class = c("high_confidence_somatic", "likely_artifact"),
    primary_reason = "x", hard_filter_reason = c("PASS", "low_depth"),
    population_category = "population_rare_or_absent", artifact_category = "artifact_not_detected",
    recurrence_category = "recurrence_non_informative", somatic_score = c(0.9, 0.1),
    technical_evidence_score = c(0.9, 0.1), artifact_score = c(0.1, 0.8),
    hotspot_match = c(TRUE, TRUE), stringsAsFactors = FALSE)
  before <- add_decision_trail(x)
  after <- add_decision_trail(add_variant_authenticity(x))
  expect_identical(before$filter_status, after$filter_status)
  expect_identical(before$filter_reasons, after$filter_reasons)
  expect_identical(before$CONFIDENCE_SCORE_BASE, after$CONFIDENCE_SCORE_BASE)
})
