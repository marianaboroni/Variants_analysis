# Decision trail. Turns the (unchanged) scientific classification into the
# unified, auditable contract required by the brief. No variant is removed here:
# every variant keeps its full trail in variants_all; retained/excluded are
# derived views. filter_status depends ONLY on the technical gate and the
# rule-based final_class, never on OncoKB (see docs/REFACTOR_AUDIT.md RISK-1).

# Final classes that are kept as somatic candidates / for review vs. excluded.
STATUS_PASS_CLASSES   <- c("high_confidence_somatic", "probable_somatic")
STATUS_REVIEW_CLASSES <- c("manual_review_required", "uncertain_tumor_only")
STATUS_FAIL_CLASSES   <- c("technical_fail", "likely_artifact", "likely_germline")

#' Add the unified decision-trail columns to a classified variant table.
#'
#' @param x variant data.frame after [classify_variants()] (has final_class,
#'   primary_reason, hard_filter_pass, hard_filter_reason, *_category columns).
#' @return `x` with `filter_status`, `filter_reasons`, `filters_failed`,
#'   `filters_passed`, `confidence_category`, `confidence_score`.
#' @keywords internal
add_decision_trail <- function(x) {
  fc <- x$final_class
  status <- rep("REVIEW", nrow(x))
  status[fc %in% STATUS_PASS_CLASSES] <- "PASS"
  status[fc %in% STATUS_FAIL_CLASSES] <- "FAIL"
  x$filter_status <- status

  # per-filter pass/fail matrix (documented, machine-readable)
  hf <- if ("hard_filter_reason" %in% names(x)) x$hard_filter_reason else rep("PASS", nrow(x))
  technical_failed <- hf != "PASS" & !is.na(hf)

  pop_cat <- if ("population_category" %in% names(x)) x$population_category else rep(NA_character_, nrow(x))
  art_cat <- if ("artifact_category" %in% names(x)) x$artifact_category else rep(NA_character_, nrow(x))
  rec_cat <- if ("recurrence_category" %in% names(x)) x$recurrence_category else rep(NA_character_, nrow(x))

  population_failed <- pop_cat %in% c("population_common")
  artifact_failed   <- art_cat %in% c("artifact_strong")
  recurrence_failed <- rec_cat %in% c("recurrence_artifact_suspected")

  failed_list <- mapply(function(tech, tech_reason, pop, art, rec) {
    f <- character()
    if (isTRUE(tech)) f <- c(f, paste0("technical:", tech_reason))
    if (isTRUE(pop)) f <- c(f, "population_common")
    if (isTRUE(art)) f <- c(f, "technical_artifact")
    if (isTRUE(rec)) f <- c(f, "cohort_recurrence_artifact")
    if (length(f) == 0) "" else paste(f, collapse = "|")
  }, technical_failed, hf, population_failed, artifact_failed, recurrence_failed, USE.NAMES = FALSE)

  passed_list <- mapply(function(tech, pop, art, rec) {
    p <- character()
    if (!isTRUE(tech)) p <- c(p, "technical")
    if (!isTRUE(pop)) p <- c(p, "population")
    if (!isTRUE(art)) p <- c(p, "artifact")
    if (!isTRUE(rec)) p <- c(p, "recurrence")
    paste(p, collapse = "|")
  }, technical_failed, population_failed, artifact_failed, recurrence_failed, USE.NAMES = FALSE)

  x$filters_failed <- failed_list
  x$filters_passed <- passed_list

  # human-readable combined reasons (primary + failed filters)
  primary <- if ("primary_reason" %in% names(x)) x$primary_reason else rep(NA_character_, nrow(x))
  x$filter_reasons <- mapply(function(pr, ff) {
    parts <- c(if (!is.na(pr) && nzchar(pr)) pr, if (nzchar(ff)) ff)
    if (length(parts) == 0) "pass_all_filters" else paste(parts, collapse = "; ")
  }, primary, failed_list, USE.NAMES = FALSE)

  # confidence - BASE (no COSMIC tumor context) is the MAIN output/filter driver.
  # The COSMIC-incorporated score is explicitly EXPERIMENTAL until weights are
  # calibrated on real data (see docs/COSMIC_DECISION_RULES.md P4).
  base_score <- if ("somatic_score" %in% names(x)) x$somatic_score else rep(NA_real_, nrow(x))
  exp_score <- if ("somatic_score_validated" %in% names(x)) x$somatic_score_validated else base_score
  x$CONFIDENCE_SCORE_BASE <- round(base_score, 4)
  x$CONFIDENCE_CATEGORY_BASE <- confidence_category_from(status, fc, base_score)
  x$COSMIC_CONTEXT_SUPPORT_SCORE <- if ("COSMIC_CONTEXT_SUPPORT_SCORE" %in% names(x))
    x$COSMIC_CONTEXT_SUPPORT_SCORE else rep(NA_real_, nrow(x))
  x$CONFIDENCE_SCORE_WITH_COSMIC_EXPERIMENTAL <- round(exp_score, 4)
  x$CONFIDENCE_CATEGORY_WITH_COSMIC_EXPERIMENTAL <- confidence_category_from(status, fc, exp_score)
  # main/primary confidence = BASE (backward-compatible column names)
  x$confidence_score <- x$CONFIDENCE_SCORE_BASE
  x$confidence_category <- x$CONFIDENCE_CATEGORY_BASE
  x
}

#' @keywords internal
confidence_category_from <- function(status, final_class, score) {
  cat <- rep("low", length(status))
  cat[status == "FAIL"] <- "excluded"
  cat[status == "REVIEW"] <- "review"
  is_pass <- status == "PASS"
  cat[is_pass & !is.na(score) & score >= 0.75] <- "high"
  cat[is_pass & !is.na(score) & score >= 0.55 & score < 0.75] <- "moderate"
  cat[is_pass & (is.na(score) | score < 0.55)] <- "low"
  cat
}

#' Add the simplified, unified interface labels. Four DISTINCT concepts, never
#' conflated: STATUS (technical/filtering decision), PRIORITY (biological
#' interest), EVIDENCE (computational predictor consensus), AUTHENTICITY
#' (technical real-vs-artifact). Detail stays in the *_reason columns.
#' @keywords internal
add_unified_labels <- function(x) {
  n <- nrow(x)
  gv <- function(col, d = NA_character_) if (col %in% names(x)) as.character(x[[col]]) else rep(d, n)

  x$STATUS <- gv("filter_status", "REVIEW")   # PASS / REVIEW / FAIL

  pc <- gv("PRIORITY_CATEGORY", "")
  prio <- c(high = "HIGH", moderate = "MEDIUM", low = "LOW", very_low = "NOT_PRIORITIZED")[pc]
  prio[is.na(prio)] <- "NOT_PRIORITIZED"
  prio[x$STATUS == "FAIL"] <- "NOT_PRIORITIZED"   # excluded variants are not prioritized
  x$PRIORITY <- unname(prio)

  ec <- gv("COMPUTATIONAL_EVIDENCE_CATEGORY", "")
  ev <- c(strong_support = "STRONG", moderate_support = "MODERATE", weak_support = "WEAK",
          conflicting = "CONFLICTING", insufficient = "INSUFFICIENT",
          not_applicable = "NOT_APPLICABLE")[ec]
  ev[is.na(ev)] <- "INSUFFICIENT"
  x$EVIDENCE <- unname(ev)

  x$AUTHENTICITY <- gv("VARIANT_AUTHENTICITY_CATEGORY", "UNCERTAIN")  # LIKELY_TRUE/UNCERTAIN/LIKELY_ARTIFACT

  # detail / reason columns (kept out of the four primary labels)
  x$status_reason <- gv("filter_reasons", "")
  x$failed_filters <- gv("filters_failed", "")
  x$priority_reasons <- gv("PRIORITY_COMPONENTS", "")
  x$evidence_summary <- gv("COMPUTATIONAL_EVIDENCE_DETAILS", "")
  # why a variant warrants human review (technical/biological conflict, etc.)
  conflict_pred <- if ("COMPUTATIONAL_CONFLICT_FLAG" %in% names(x)) x$COMPUTATIONAL_CONFLICT_FLAG %in% c(TRUE, "TRUE") else rep(FALSE, n)
  bio_hi_tech_lo <- suppressWarnings(as.numeric(gv("BIOLOGICAL_SUPPORT_SCORE", "0"))) >= 0.5 &
    suppressWarnings(as.numeric(gv("TECHNICAL_AUTHENTICITY_SCORE", "1"))) < 0.5
  x$review_reasons <- mapply(function(st, auth, cp, bt) {
    r <- character()
    if (identical(st, "REVIEW")) r <- c(r, "flagged_for_review")
    if (identical(auth, "LIKELY_ARTIFACT")) r <- c(r, "technical_artifact_evidence")
    if (isTRUE(cp)) r <- c(r, "conflicting_predictors")
    if (isTRUE(bt)) r <- c(r, "biologically_interesting_but_technically_weak")
    if (length(r) == 0) "" else paste(r, collapse = ";")
  }, x$STATUS, x$AUTHENTICITY, conflict_pred, bio_hi_tech_lo, USE.NAMES = FALSE)
  x
}

#' Split the annotated table into retained / excluded views (no data loss:
#' both are subsets of variants_all).
#' @keywords internal
split_by_status <- function(x) {
  list(
    retained = x[x$filter_status %in% c("PASS", "REVIEW"), , drop = FALSE],
    excluded = x[x$filter_status == "FAIL", , drop = FALSE]
  )
}

#' Build a per-filter audit table (definition + counts + loss at each stage).
#' @keywords internal
build_filter_audit <- function(x) {
  total <- nrow(x)
  count_fail <- function(pattern) sum(grepl(pattern, x$filters_failed, fixed = TRUE))
  data.frame(
    filter = c("technical", "population_common", "technical_artifact",
               "cohort_recurrence_artifact", "final_PASS", "final_REVIEW", "final_FAIL"),
    category = c("quality/technical", "population_frequency", "artifact_evidence",
                 "recurrence_evidence", "classification", "classification", "classification"),
    definition = c(
      "adaptive depth/alt/VAF/TLOD/MBQ/MMQ + caller FILTER + sample QC gate",
      "population MAF >= common threshold (likely germline)",
      "strong technical artifact evidence (bias/PoN/low-quality)",
      "recurrent in cohort above artifact fraction",
      "retained somatic candidate (high/probable somatic)",
      "flagged for manual review (uncertain/conflicting)",
      "excluded (technical_fail / likely_artifact / likely_germline)"),
    n_variants = c(
      count_fail("technical:"), count_fail("population_common"),
      count_fail("technical_artifact"), count_fail("cohort_recurrence_artifact"),
      sum(x$filter_status == "PASS"), sum(x$filter_status == "REVIEW"),
      sum(x$filter_status == "FAIL")),
    pct_of_total = NA_real_,
    stringsAsFactors = FALSE
  ) -> audit
  audit$pct_of_total <- round(100 * audit$n_variants / max(total, 1), 2)
  audit
}
