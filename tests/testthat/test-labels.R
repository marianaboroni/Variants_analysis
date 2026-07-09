test_that("unified labels map the four DISTINCT concepts correctly", {
  x <- data.frame(
    filter_status = c("PASS", "REVIEW", "FAIL"),
    PRIORITY_CATEGORY = c("high", "moderate", "low"),
    COMPUTATIONAL_EVIDENCE_CATEGORY = c("strong_support", "conflicting", "not_applicable"),
    VARIANT_AUTHENTICITY_CATEGORY = c("LIKELY_TRUE", "UNCERTAIN", "LIKELY_ARTIFACT"),
    filter_reasons = c("a", "b", "c"), filters_failed = c("", "", "low_depth"),
    PRIORITY_COMPONENTS = "x", COMPUTATIONAL_EVIDENCE_DETAILS = "y",
    COMPUTATIONAL_CONFLICT_FLAG = c(FALSE, TRUE, FALSE),
    BIOLOGICAL_SUPPORT_SCORE = c(0.2, 0.7, 0.1),
    TECHNICAL_AUTHENTICITY_SCORE = c(0.9, 0.3, 0.1),
    stringsAsFactors = FALSE)
  out <- add_unified_labels(x)
  expect_equal(out$STATUS, c("PASS", "REVIEW", "FAIL"))
  expect_equal(out$PRIORITY, c("HIGH", "MEDIUM", "NOT_PRIORITIZED"))  # FAIL -> NOT_PRIORITIZED
  expect_equal(out$EVIDENCE, c("STRONG", "CONFLICTING", "NOT_APPLICABLE"))
  expect_equal(out$AUTHENTICITY, c("LIKELY_TRUE", "UNCERTAIN", "LIKELY_ARTIFACT"))
  expect_true(all(c("status_reason", "failed_filters", "priority_reasons",
                    "evidence_summary", "review_reasons") %in% names(out)))
  # REVIEW row: conflicting predictors + biologically interesting but technically weak
  expect_true(grepl("conflicting_predictors", out$review_reasons[2]))
  expect_true(grepl("biologically_interesting_but_technically_weak", out$review_reasons[2]))
})

test_that("STATUS/PRIORITY/EVIDENCE/AUTHENTICITY are independent columns", {
  x <- data.frame(filter_status = "PASS", PRIORITY_CATEGORY = "low",
    COMPUTATIONAL_EVIDENCE_CATEGORY = "insufficient",
    VARIANT_AUTHENTICITY_CATEGORY = "LIKELY_ARTIFACT", stringsAsFactors = FALSE)
  out <- add_unified_labels(x)
  # a PASS variant can still be LOW priority, INSUFFICIENT evidence, LIKELY_ARTIFACT
  expect_equal(out$STATUS, "PASS"); expect_equal(out$PRIORITY, "LOW")
  expect_equal(out$EVIDENCE, "INSUFFICIENT"); expect_equal(out$AUTHENTICITY, "LIKELY_ARTIFACT")
})
