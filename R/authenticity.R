# Authenticity vs biological relevance -- two DISTINCT concepts kept separate:
#
#   TECHNICAL_AUTHENTICITY_SCORE : is the CALL technically real (not an artifact)?
#                                   depends ONLY on technical/quality signals.
#   BIOLOGICAL_SUPPORT_SCORE     : is the variant biologically relevant?
#                                   depends ONLY on consequence/predictors/hotspot/
#                                   driver/COSMIC/clinical evidence.
#   VARIANT_REVIEW_SCORE         : triage score -- surfaces variants that are
#                                   biologically interesting BUT technically weak
#                                   (the ones a human should look at).
#
# Crucially, VARIANT_AUTHENTICITY_CATEGORY derives ONLY from technical evidence, so
# strong biological evidence can NEVER mask a strong technical artifact
# (a hotspot with severe strand bias stays LIKELY_ARTIFACT). None of these scores
# set filter_status. All weights/thresholds are HEURISTIC and configurable
# (see docs/VARIANT_AUTHENTICITY.md).

#' Add the technical-authenticity, biological-support and review scores.
#'
#' @param x variant data.frame after scoring/cosmic/driver/predictor layers.
#' @param cfg resolved config (reads `authenticity:` thresholds).
#' @return `x` with `ARTIFACT_EVIDENCE_SCORE`, `TECHNICAL_AUTHENTICITY_SCORE`,
#'   `BIOLOGICAL_SUPPORT_SCORE`, `VARIANT_REVIEW_SCORE`,
#'   `VARIANT_AUTHENTICITY_SCORE` (alias of technical), `VARIANT_AUTHENTICITY_CATEGORY`
#'   (LIKELY_TRUE/UNCERTAIN/LIKELY_ARTIFACT, technical only),
#'   `VARIANT_AUTHENTICITY_COMPONENTS`.
#' @export
#' @examples
#' add_variant_authenticity(data.frame(technical_evidence_score = 0.9, artifact_score = 0.1))
add_variant_authenticity <- function(x, cfg = NULL) {
  n <- nrow(x)
  getn <- function(col, d = NA_real_) if (col %in% names(x)) suppressWarnings(as.numeric(x[[col]])) else rep(d, n)
  getl <- function(col) if (col %in% names(x)) (x[[col]] %in% c(TRUE, "TRUE")) else rep(FALSE, n)

  # ---- TECHNICAL ONLY: artifact evidence (higher = more likely artifact) --
  artifact <- zero_na(getn("artifact_score", 0))
  bias <- pmax(zero_na(getn("orientation_bias")), zero_na(getn("strand_artifact")))
  reinforce <- pmax(
    ifelse(getl("pon_flag"), 0.6, 0),
    ifelse(getl("clustered_events"), 0.4, 0),
    ifelse(getl("weak_evidence"), 0.5, 0),
    ifelse(!is.na(bias) & bias > 0.8, 0.75, 0))  # strong strand/orientation bias => strong artifact signal
  artifact_evidence <- bounded01(pmax(artifact, reinforce))

  # ---- TECHNICAL ONLY: authenticity (higher = more likely a real call) ----
  tech_quality <- zero_na(getn("technical_evidence_score", 0))
  technical_authenticity <- bounded01(0.6 * tech_quality + 0.4 * (1 - artifact_evidence))

  # ---- BIOLOGICAL ONLY: relevance support (never affects authenticity) ----
  consequence_comp <- pmin(zero_na(getn("impact_rank")) / 3, 1)
  recurrence <- bounded01(pmax(
    zero_na(getn("COSMIC_GLOBAL_RECURRENCE_SCORE")),
    zero_na(getn("COSMIC_CONTEXT_SUPPORT_SCORE")),
    ifelse(getl("hotspot_match"), 1, 0),
    zero_na(getn("driver_score"))))
  computational <- zero_na(getn("COMPUTATIONAL_EVIDENCE_SCORE"))
  computational[getl("COMPUTATIONAL_CONFLICT_FLAG")] <- computational[getl("COMPUTATIONAL_CONFLICT_FLAG")] * 0.5
  clinical <- ifelse(grepl("pathogenic", tolower(as.character(
    if ("clinvar_significance" %in% names(x)) x$clinvar_significance else rep(NA, n)))) &
    !grepl("benign", tolower(as.character(
    if ("clinvar_significance" %in% names(x)) x$clinvar_significance else rep(NA, n)))), 1, 0)
  biological_support <- bounded01(0.30 * consequence_comp + 0.30 * recurrence +
    0.25 * computational + 0.15 * clinical)

  # ---- review triage: biologically interesting BUT technically weak -------
  review_score <- bounded01(0.5 * (1 - technical_authenticity) + 0.5 * biological_support)

  hi <- cfg_get(cfg, c("authenticity", "likely_true_min"), 0.66)
  lo <- cfg_get(cfg, c("authenticity", "likely_artifact_max"), 0.40)
  # authenticity category is TECHNICAL ONLY -> biology cannot mask an artifact
  category <- rep("UNCERTAIN", n)
  category[technical_authenticity >= hi] <- "LIKELY_TRUE"
  category[technical_authenticity < lo] <- "LIKELY_ARTIFACT"

  x$ARTIFACT_EVIDENCE_SCORE <- round(artifact_evidence, 4)
  x$TECHNICAL_AUTHENTICITY_SCORE <- round(technical_authenticity, 4)
  x$BIOLOGICAL_SUPPORT_SCORE <- round(biological_support, 4)
  x$VARIANT_REVIEW_SCORE <- round(review_score, 4)
  x$VARIANT_AUTHENTICITY_SCORE <- round(technical_authenticity, 4)   # authenticity == technical
  x$VARIANT_AUTHENTICITY_CATEGORY <- category
  x$VARIANT_AUTHENTICITY_COMPONENTS <- sprintf(
    "technical_quality=%.2f;artifact_evidence=%.2f | biology[consequence=%.2f;recurrence=%.2f;computational=%.2f;clinical=%.2f]",
    tech_quality, artifact_evidence, consequence_comp, recurrence, computational, clinical)
  x
}

#' @keywords internal
zero_na <- function(v) { v[is.na(v)] <- 0; v }
