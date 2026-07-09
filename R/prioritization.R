# Decomposed, heuristic prioritization of variants for review. Prioritization
# RANKS variants; it does not filter. Each component is exposed separately so the
# score is never an opaque number. Weights are HEURISTIC until calibrated
# (see docs/FILTERING_STRATEGY.md). PRIORITY uses the BASE (non-COSMIC-context)
# evidence; COSMIC global recurrence contributes a small informational term.

# component weights (documented, heuristic)
PRIORITY_WEIGHTS <- c(consequence = 0.25, computational = 0.20, hotspot = 0.25,
                      driver = 0.15, clinical = 0.10, cosmic_global = 0.05)

#' @keywords internal
priority_consequence_component <- function(x) {
  ir <- if ("impact_rank" %in% names(x)) x$impact_rank else rep(NA_real_, nrow(x))
  comp <- pmin(zero_if_na_pri(ir) / 3, 1)
  # truncating / splice-disruptive get full weight regardless of rank granularity
  if ("is_truncating" %in% names(x)) comp[isTRUE_vec(x$is_truncating)] <- 1
  if ("is_splice_disruptive" %in% names(x)) comp[isTRUE_vec(x$is_splice_disruptive)] <- pmax(comp[isTRUE_vec(x$is_splice_disruptive)], 0.9)
  comp
}
#' @keywords internal
priority_computational_component <- function(x) {
  if (!"COMPUTATIONAL_EVIDENCE_SCORE" %in% names(x)) return(rep(0, nrow(x)))
  s <- x$COMPUTATIONAL_EVIDENCE_SCORE; s[is.na(s)] <- 0
  # conflicting evidence is down-weighted
  if ("COMPUTATIONAL_CONFLICT_FLAG" %in% names(x)) s[isTRUE_vec(x$COMPUTATIONAL_CONFLICT_FLAG)] <- s[isTRUE_vec(x$COMPUTATIONAL_CONFLICT_FLAG)] * 0.5
  s
}
#' @keywords internal
priority_hotspot_component <- function(x) {
  h <- rep(0, nrow(x))
  if ("hotspot_match" %in% names(x)) h[isTRUE_vec(x$hotspot_match)] <- 0.85
  if ("hotspot_tumor_specific" %in% names(x)) h[isTRUE_vec(x$hotspot_tumor_specific)] <- 1
  h
}
#' @keywords internal
priority_driver_component <- function(x) {
  d <- rep(0, nrow(x))
  if ("gene_driver_match" %in% names(x)) d[isTRUE_vec(x$gene_driver_match)] <- 0.6
  if ("gene_driver_tumor_specific" %in% names(x)) d[isTRUE_vec(x$gene_driver_tumor_specific)] <- 0.85
  if ("driver_class" %in% names(x)) d[x$driver_class %in% c("known_driver")] <- 1
  d
}
#' @keywords internal
priority_clinical_component <- function(x) {
  cv <- tolower(as.character(if ("clinvar_significance" %in% names(x)) x$clinvar_significance
                             else if ("CLIN_SIG" %in% names(x)) x$CLIN_SIG else rep(NA, nrow(x))))
  c0 <- rep(0, nrow(x))
  c0[grepl("pathogenic", cv) & !grepl("benign|conflicting", cv)] <- 1
  c0[grepl("likely_pathogenic", cv)] <- pmax(c0[grepl("likely_pathogenic", cv)], 0.8)
  c0[grepl("uncertain|vus", cv)] <- pmax(c0[grepl("uncertain|vus", cv)], 0.2)
  c0
}
#' @keywords internal
priority_cosmic_global_component <- function(x) {
  if (!"COSMIC_GLOBAL_RECURRENCE_SCORE" %in% names(x)) return(rep(0, nrow(x)))
  s <- x$COSMIC_GLOBAL_RECURRENCE_SCORE; s[is.na(s)] <- 0; s
}

#' Add decomposed prioritization columns to a variant table.
#' @keywords internal
add_prioritization <- function(x, cfg = NULL) {
  cc <- priority_consequence_component(x)
  comp <- priority_computational_component(x)
  hs <- priority_hotspot_component(x)
  dr <- priority_driver_component(x)
  cl <- priority_clinical_component(x)
  cg <- priority_cosmic_global_component(x)

  x$PRIORITY_CONSEQUENCE_COMPONENT <- round(cc, 4)
  x$PRIORITY_COMPUTATIONAL_COMPONENT <- round(comp, 4)
  x$PRIORITY_HOTSPOT_COMPONENT <- round(hs, 4)
  x$PRIORITY_DRIVER_COMPONENT <- round(dr, 4)
  x$PRIORITY_CLINICAL_COMPONENT <- round(cl, 4)
  x$PRIORITY_COSMIC_GLOBAL_COMPONENT <- round(cg, 4)

  w <- PRIORITY_WEIGHTS
  score <- w["consequence"] * cc + w["computational"] * comp + w["hotspot"] * hs +
    w["driver"] * dr + w["clinical"] * cl + w["cosmic_global"] * cg
  x$PRIORITY_SCORE_BASE <- round(as.numeric(score), 4)
  x$PRIORITY_CATEGORY <- cut(x$PRIORITY_SCORE_BASE, breaks = c(-Inf, 0.15, 0.35, 0.6, Inf),
    labels = c("very_low", "low", "moderate", "high"), right = FALSE) |> as.character()
  x$PRIORITY_COMPONENTS <- sprintf(
    "consequence=%.2f;computational=%.2f;hotspot=%.2f;driver=%.2f;clinical=%.2f;cosmic_global=%.2f",
    cc, comp, hs, dr, cl, cg)
  x
}

#' @keywords internal
zero_if_na_pri <- function(v) { v[is.na(v)] <- 0; v }
