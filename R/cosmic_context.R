# Tumor-type-stratified interpretation of COSMIC evidence.
#
# COSMIC evidence is NOT reduced to global presence and NOT reduced to a single
# score. We compute, separately and traceably:
#   - global recurrence  (from total occurrences across ALL tumor types)
#   - matched-tumor recurrence (from occurrences in the SAMPLE's tumor type only)
#   - a controlled context status + granular match level
#   - a BANDED context-support score whose tiers cannot be crossed by raw counts
# Harmonization uses an explicit, versioned table - never fuzzy text matching.
# See docs/COSMIC_DECISION_RULES.md for the rules, bands and weights.

# Controlled context statuses.
CONTEXT_STATUSES <- c(
  "exact_match", "compatible_match",
  "pan_cancer_with_sample_tumor", "pan_cancer_without_sample_tumor",
  "other_tumor_only", "tumor_type_unknown",
  "cosmic_tumor_type_missing", "no_cosmic_match")

# Granular match levels (best category quality between sample and COSMIC).
MATCH_LEVELS <- c("exact_site_histology_subtype", "exact_site_histology",
                  "compatible_histology", "broad_site", "global_only", "none")

# ordinal rank of harmonization mapping levels (higher = more specific)
MAP_LEVEL_RANK <- c(exact = 3L, compatible = 2L, broad = 1L)

# tier specificity value (pure, recurrence-independent) per status
CONTEXT_TIER_SPECIFICITY <- c(
  exact_match = 1.00, compatible_match = 0.65,
  pan_cancer_with_sample_tumor = 0.45,
  other_tumor_only = 0.15, pan_cancer_without_sample_tumor = 0.15)

# Non-overlapping support bands [floor, ceiling] guaranteeing, by construction,
# exact_match > compatible_match > pan_cancer_with_sample_tumor > other_tumor_only
# for ALL occurrence counts. Recurrence only modulates WITHIN a band.
CONTEXT_SUPPORT_BAND <- list(
  exact_match = c(0.80, 1.00),
  compatible_match = c(0.55, 0.75),
  pan_cancer_with_sample_tumor = c(0.35, 0.50),
  other_tumor_only = c(0.05, 0.25),
  pan_cancer_without_sample_tumor = c(0.05, 0.25))

#' Base-size-aware transform of a raw recurrence count. NOTE: this is a
#' transformed RAW RECURRENCE COUNT, not a frequency/prevalence and NOT a
#' correction for the uneven number of samples represented per tumor type in
#' COSMIC (that would require per-tumor denominators, which raw COSMIC does not
#' provide). See docs/COSMIC_DECISION_RULES.md.
#' @keywords internal
recurrence_transform <- function(count) {
  count[is.na(count)] <- 0
  pmin(1, log10(count + 1) / 3)
}

#' Load the versioned tumor-type harmonization table.
#' @return a normalized data.frame, or NULL if none is configured/found.
#' @keywords internal
load_tumor_type_mapping <- function(cfg) {
  path <- cfg_get(cfg, c("cosmic", "tumor_type_mapping"), NULL)
  if (is.null(path) || is.na(path)) {
    dflt <- file.path("config", "cosmic_tumor_type_mapping.tsv")
    if (file.exists(dflt)) path <- dflt else return(NULL)
  }
  if (!file.exists(path)) return(NULL)
  m <- read_variants(path, "\t")
  m$input_tumor_type <- toupper(trimws(as.character(m$input_tumor_type)))
  st <- if ("input_tumor_subtype" %in% names(m)) as.character(m$input_tumor_subtype) else rep("", nrow(m))
  st[is.na(st)] <- ""
  m$input_tumor_subtype <- toupper(trimws(st))
  m$cosmic_primary_site <- norm_cosmic_cat(m$cosmic_primary_site)
  m$cosmic_histology <- norm_cosmic_cat(m$cosmic_histology)
  m$match_level <- tolower(trimws(as.character(m$match_level)))
  m$match_level[!m$match_level %in% names(MAP_LEVEL_RANK)] <- "broad"
  attr(m, "version") <- paste(unique(stats::na.omit(m$mapping_version)), collapse = ",")
  m
}

#' Allowed COSMIC categories (site/histology/level/via_subtype) for a sample.
#' `via_subtype` marks rows that came from a subtype-specific mapping entry.
#' @return data.frame; 0 rows if the type is not harmonized.
#' @keywords internal
harmonize_sample_tumor <- function(tumor_type, tumor_subtype, mapping) {
  if (is.null(mapping)) return(mapping_empty())
  clean <- function(z) { z <- as.character(z); z[is.na(z)] <- ""; toupper(trimws(z)) }
  tt <- clean(tumor_type); st <- clean(tumor_subtype)
  if (tt %in% c("", "NA", "PANCANCER", "UNKNOWN")) return(mapping_empty())
  keep <- !is.na(mapping$input_tumor_type) & mapping$input_tumor_type == tt &
    (mapping$input_tumor_subtype == "" | mapping$input_tumor_subtype == st)
  keep[is.na(keep)] <- FALSE
  sel <- mapping[keep, , drop = FALSE]
  if (nrow(sel) == 0) return(mapping_empty())
  data.frame(site = sel$cosmic_primary_site, hist = sel$cosmic_histology,
             level = sel$match_level,
             via_subtype = sel$input_tumor_subtype != "" & st != "",
             stringsAsFactors = FALSE)
}
#' @keywords internal
mapping_empty <- function() data.frame(site = character(), hist = character(),
  level = character(), via_subtype = logical(), stringsAsFactors = FALSE)

#' Classify the tumor context of a single COSMIC-matched variant from its
#' long-format category rows.
#'
#' @param long data.frame(site, hist, subtype, count) for the variant (one row
#'   per COSMIC tumor category).
#' @param allowed harmonized sample categories from [harmonize_sample_tumor()].
#' @param sample_known whether the sample tumor type was harmonized.
#' @param sample_subtype_available whether the sample provided a subtype.
#' @param pan_threshold min distinct tumor categories to call pan-cancer.
#' @return one-row data.frame with the full COSMIC context column set.
#' @keywords internal
classify_cosmic_tumor_context <- function(long, allowed, sample_known,
                                          sample_subtype_available = FALSE,
                                          pan_threshold = 5,
                                          os4_min = 10, om4_min = 3) {
  total <- sum(long$count)
  n_distinct <- nrow(long)
  has_info <- any(nzchar(long$site) | nzchar(long$hist))
  global_pan <- n_distinct >= pan_threshold
  global_rec <- recurrence_transform(total)
  lab <- function(df) if (nrow(df) == 0) "" else paste(ifelse(nzchar(df$hist),
    paste(df$site, df$hist, sep = "/"), df$site), collapse = ";")

  # --- not evaluable: no tumor info, or sample tumor type unknown -----------
  if (!has_info) {
    return(ctx_row("cosmic_tumor_type_missing", "global_only", FALSE, FALSE, FALSE,
      matching = 0, other = total, mtypes = "", otypes = lab(long),
      evaluable = FALSE, matched_rec = 0, global_rec = global_rec,
      specificity = NA_real_, support = NA_real_,
      global_pan = global_pan, pan_count = n_distinct, includes_sample = FALSE,
      total = total, os4_min = os4_min, om4_min = om4_min,
      reason = "COSMIC record(s) lack usable tumor-context information"))
  }
  if (!sample_known) {
    return(ctx_row("tumor_type_unknown", "global_only", FALSE, FALSE, FALSE,
      matching = 0, other = total, mtypes = "", otypes = lab(long),
      evaluable = FALSE, matched_rec = 0, global_rec = global_rec,
      specificity = NA_real_, support = NA_real_,
      global_pan = global_pan, pan_count = n_distinct, includes_sample = FALSE,
      total = total, os4_min = os4_min, om4_min = om4_min,
      reason = "sample tumor type not provided or not harmonized; global evidence only"))
  }

  # --- per-category harmonization ------------------------------------------
  lvl <- integer(n_distinct); hist_exact <- logical(n_distinct); via_st <- logical(n_distinct)
  for (i in seq_len(n_distinct)) {
    cand <- allowed[!is.na(allowed$site) & allowed$site == long$site[i], , drop = FALSE]
    if (nrow(cand) == 0) next
    applies <- cand$hist == "" | cand$hist == long$hist[i]
    cand <- cand[applies, , drop = FALSE]
    if (nrow(cand) == 0) next
    ranks <- MAP_LEVEL_RANK[cand$level]
    lvl[i] <- as.integer(max(c(0L, ranks[!is.na(ranks)])))
    exact_rows <- cand$level == "exact" & cand$hist != "" & cand$hist == long$hist[i]
    hist_exact[i] <- any(exact_rows)
    via_st[i] <- any(cand$via_subtype[exact_rows])
  }
  matched <- lvl > 0
  matching_occ <- sum(long$count[matched])
  other_occ <- total - matching_occ
  includes_sample <- matching_occ > 0

  site_match <- includes_sample
  histology_match <- any(hist_exact[matched])
  subtype_match <- sample_subtype_available && any(matched & via_st & nzchar(long$subtype))

  match_level <- if (!includes_sample) "global_only"
    else if (subtype_match) "exact_site_histology_subtype"
    else if (histology_match) "exact_site_histology"
    else if (max(lvl[matched]) == MAP_LEVEL_RANK[["compatible"]]) "compatible_histology"
    else "broad_site"

  if (includes_sample) {
    if (global_pan) { status <- "pan_cancer_with_sample_tumor"
      reason <- sprintf("recurrent across %d tumor categories including the sample type", n_distinct) }
    else if (match_level %in% c("exact_site_histology_subtype", "exact_site_histology")) {
      status <- "exact_match"; reason <- "observed in the same tumor type/histology" }
    else { status <- "compatible_match"; reason <- "observed in a compatible tumor category" }
  } else {
    if (global_pan) { status <- "pan_cancer_without_sample_tumor"
      reason <- sprintf("recurrent across %d tumor categories but NOT the sample type (global evidence only)", n_distinct) }
    else { status <- "other_tumor_only"; reason <- "present in COSMIC only in other tumor type(s)" }
  }

  matched_rec <- recurrence_transform(matching_occ)
  specificity <- unname(CONTEXT_TIER_SPECIFICITY[status])
  band <- CONTEXT_SUPPORT_BAND[[status]]
  # recurrence modulates within the (contextual) band; tumor-matched tiers use
  # matched-tumor recurrence, non-matched tiers use global recurrence.
  rec_for_band <- if (includes_sample) matched_rec else global_rec
  support <- band[1] + (band[2] - band[1]) * rec_for_band

  ctx_row(status, match_level, site_match, histology_match, subtype_match,
    matching = matching_occ, other = other_occ,
    mtypes = lab(long[matched, , drop = FALSE]), otypes = lab(long[!matched, , drop = FALSE]),
    evaluable = TRUE, matched_rec = matched_rec, global_rec = global_rec,
    specificity = specificity, support = round(support, 4),
    global_pan = global_pan, pan_count = n_distinct, includes_sample = includes_sample,
    total = total, os4_min = os4_min, om4_min = om4_min, reason = reason)
}

#' @keywords internal
ctx_row <- function(status, match_level, site_match, histology_match, subtype_match,
                    matching, other, mtypes, otypes, evaluable, matched_rec, global_rec,
                    specificity, support, global_pan, pan_count, includes_sample,
                    total, os4_min, om4_min, reason) {
  # OS4/OM4 are attributable ONLY from tumor-CONTEXT recurrence (occurrences in
  # the sample's tumor type), never from global-only recurrence. See
  # docs/COSMIC_OS4_OM4_AUDIT.md and COSMIC_DECISION_RULES.md.
  os4_contextual <- isTRUE(evaluable) && isTRUE(includes_sample) && matching >= os4_min
  om4_contextual <- isTRUE(evaluable) && isTRUE(includes_sample) &&
    matching >= om4_min && matching < os4_min
  global_evidence_only <- !isTRUE(includes_sample)   # ctx_row is only called on genomic matches
  global_recurrent <- total >= om4_min
  interpretation <- if (isTRUE(evaluable)) "evaluable" else "not_evaluable"
  data.frame(
    COSMIC_TUMOR_CONTEXT_STATUS = status,
    COSMIC_MATCH_LEVEL = match_level,
    COSMIC_SITE_MATCH = site_match,
    COSMIC_HISTOLOGY_MATCH = histology_match,
    COSMIC_SUBTYPE_MATCH = subtype_match,
    COSMIC_MATCHING_TUMOR_OCCURRENCES = matching,
    COSMIC_OTHER_TUMOR_OCCURRENCES = other,
    COSMIC_MATCHING_TUMOR_TYPES = mtypes,
    COSMIC_OTHER_TUMOR_TYPES = otypes,
    COSMIC_CONTEXT_EVALUABLE = evaluable,
    COSMIC_GLOBAL_RECURRENCE_SCORE = round(global_rec, 4),
    COSMIC_MATCHED_TUMOR_RECURRENCE_SCORE = round(matched_rec, 4),
    COSMIC_TUMOR_SPECIFICITY_SCORE = specificity,
    COSMIC_CONTEXT_SUPPORT_SCORE = support,
    COSMIC_GLOBAL_PANCANCER_RECURRENT = global_pan,
    COSMIC_PANCANCER_TUMOR_COUNT = pan_count,
    COSMIC_PANCANCER_INCLUDES_SAMPLE_TUMOR = includes_sample,
    OS4_CONTEXTUAL = os4_contextual,
    OM4_CONTEXTUAL = om4_contextual,
    COSMIC_GLOBAL_RECURRENT = global_recurrent,
    COSMIC_GLOBAL_EVIDENCE_ONLY = global_evidence_only,
    COSMIC_CONTEXT_INTERPRETATION = interpretation,
    COSMIC_TUMOR_CONTEXT_REASON = reason,
    stringsAsFactors = FALSE)
}

#' Exploratory weight-sensitivity analysis for the COSMIC context score.
#'
#' Recomputes an ALTERNATIVE linear score `w_rec * global_recurrence +
#' w_ctx * specificity` (NOT the shipped banded score) under several plausible
#' recurrence/context weightings, and reports ranking stability and category
#' churn. Used to justify keeping the score exploratory. Operates on an
#' already-annotated variant table.
#'
#' @param x annotated variants (must have the COSMIC context columns).
#' @param weights list of c(recurrence, context) pairs summing to 1.
#' @return a list with per-weighting categories, Spearman rank correlations vs a
#'   reference weighting, and the count of variants that change category.
#' @keywords internal
cosmic_weight_sensitivity <- function(x, weights = list(c(0.2, 0.8), c(0.4, 0.6), c(0.5, 0.5))) {
  ev <- isTRUE_vec(x$COSMIC_CONTEXT_EVALUABLE)
  rec <- ifelse(is.na(x$COSMIC_GLOBAL_RECURRENCE_SCORE), 0, x$COSMIC_GLOBAL_RECURRENCE_SCORE)[ev]
  spec <- ifelse(is.na(x$COSMIC_TUMOR_SPECIFICITY_SCORE), 0, x$COSMIC_TUMOR_SPECIFICITY_SCORE)[ev]
  cat_of <- function(s) cut(s, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf),
                            labels = c("low", "moderate", "high", "very_high"))
  scores <- lapply(weights, function(w) w[1] * rec + w[2] * spec)
  names(scores) <- vapply(weights, function(w) sprintf("%g/%g", w[1] * 100, w[2] * 100), character(1))
  cats <- lapply(scores, cat_of)
  ref <- scores[[length(scores)]]
  cors <- vapply(scores, function(s) if (length(s) < 2 || stats::sd(s) == 0) NA_real_
                 else stats::cor(s, ref, method = "spearman"), numeric(1))
  changed <- if (length(cats) >= 2) sum(as.character(cats[[1]]) != as.character(cats[[length(cats)]])) else 0L
  list(n_evaluable = sum(ev), weightings = names(scores),
       rank_correlation_vs_5050 = cors, n_category_changes_20_80_vs_50_50 = changed,
       category_tables = lapply(cats, function(c) as.list(table(c))))
}
