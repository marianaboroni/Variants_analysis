# Functional-predictor layer: registry-driven detection, normalization,
# applicability-by-consequence, and family-aware computational evidence.
# Predictors NEVER set filter_status. They contribute one decomposed component
# to prioritization. See docs/PREDICTOR_RULES.md and inst/config/predictor_registry.yml.

PREDICTOR_FAMILIES <- c("protein_function", "ensemble_missense", "splicing",
                        "conservation", "clinical_database")

#' Load the predictor registry (installed inst/config or dev path).
#' @keywords internal
load_predictor_registry <- function(path = NULL) {
  if (is.null(path)) {
    p <- system.file("config", "predictor_registry.yml", package = "tumoronly")
    if (!nzchar(p) || !file.exists(p)) {
      # dev-mode: search upward from cwd for inst/config/predictor_registry.yml
      dir <- getwd()
      for (i in 1:5) {
        cand <- file.path(dir, "inst", "config", "predictor_registry.yml")
        if (file.exists(cand)) { p <- cand; break }
        dir <- dirname(dir)
      }
    }
    path <- p
  }
  if (is.null(path) || !nzchar(path) || !file.exists(path))
    stop("predictor registry not found (looked in inst/config/predictor_registry.yml)", call. = FALSE)
  yaml::read_yaml(path)
}

#' Union of all predictor alias column names from the registry, so the VCF
#' reader can surface EVERY predictor present in CSQ (not a fixed subset).
#' Falls back to a known set if the registry cannot be read.
#' @keywords internal
predictor_alias_union <- function() {
  reg <- tryCatch(load_predictor_registry(), error = function(e) NULL)
  if (is.null(reg)) {
    return(c("SIFT", "PolyPhen", "REVEL", "REVEL_score", "CADD_PHRED", "CADD_phred",
             "AlphaMissense_score", "AlphaMissense_pred", "MetaLR_score", "MetaLR_pred",
             "MetaRNN_score", "MetaRNN_pred", "MutationTaster_score", "MutationTaster_pred",
             "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL"))
  }
  unique(unlist(lapply(reg, function(s) s$aliases), use.names = FALSE))
}

#' Classify each variant's consequence into a predictor scope.
#' @keywords internal
predictor_variant_scope <- function(consequence) {
  c0 <- tolower(as.character(consequence))
  scope <- rep("other", length(c0))
  scope[grepl("synonymous|stop_retained", c0)] <- "synonymous"
  scope[grepl("missense|protein_altering|inframe", c0)] <- "missense"
  scope[grepl("splice", c0)] <- "splicing"
  scope[grepl("frameshift|stop_gained|start_lost|stop_lost|nonsense|transcript_ablation", c0)] <- "truncating"
  scope
}

#' Extract a numeric score and a categorical label from a raw predictor value,
#' handling the VEP `label(score)` form (e.g. "deleterious(0.01)").
#' @keywords internal
extract_pred_value <- function(v) {
  v <- as.character(v)
  has_paren <- grepl("\\(.*\\)", v)
  score <- rep(NA_real_, length(v))
  label <- rep(NA_character_, length(v))
  # label(score)
  score[has_paren] <- suppressWarnings(as.numeric(sub(".*\\(([-0-9.eE]+)\\).*", "\\1", v[has_paren])))
  label[has_paren] <- tolower(trimws(sub("\\(.*$", "", v[has_paren])))
  # plain numeric
  plain_num <- !has_paren & !is.na(suppressWarnings(as.numeric(v)))
  score[plain_num] <- suppressWarnings(as.numeric(v[plain_num]))
  # plain label
  plain_lab <- !has_paren & !plain_num & !is.na(v) & v != "" & v != "."
  label[plain_lab] <- tolower(trimws(v[plain_lab]))
  list(score = score, label = label)
}

#' Normalize one predictor from the variant table into score/label/present.
#' @keywords internal
normalize_one_predictor <- function(x, spec) {
  cols <- spec$aliases[spec$aliases %in% names(x)]
  if (length(cols) == 0) return(NULL)
  n <- nrow(x)
  score <- rep(NA_real_, n); label <- rep(NA_character_, n)
  agg_max <- identical(spec$aggregate, "max")
  for (col in cols) {
    ev <- extract_pred_value(x[[col]])
    if (agg_max) {
      score <- pmax(score, ev$score, na.rm = TRUE)
    } else {
      score[is.na(score)] <- ev$score[is.na(score)]
    }
    label[is.na(label)] <- ev$label[is.na(label)]
  }
  score[is.infinite(score)] <- NA_real_
  # range validation (out-of-range -> NA, recorded as a QC concern by caller)
  n_oor <- 0L
  if (!is.null(spec$score_min) && !is.null(spec$score_max)) {
    oor <- !is.na(score) & (score < spec$score_min | score > spec$score_max)
    n_oor <- sum(oor); score[oor] <- NA_real_
  }
  present <- !is.na(score) | !is.na(label)
  list(score = score, label = label, present = present, columns = cols, n_out_of_range = n_oor)
}

#' Is a predictor value "damaging"? Uses score+direction when available, else the
#' categorical label. Returns logical (NA when not determinable).
#' @keywords internal
predictor_is_damaging <- function(score, label, spec) {
  dmg <- rep(NA, length(score))
  have_score <- !is.na(score)
  if (identical(spec$direction, "lower_is_more_damaging")) {
    dmg[have_score] <- score[have_score] <= spec$damaging_threshold
  } else {
    dmg[have_score] <- score[have_score] >= spec$damaging_threshold
  }
  # fall back to label where score missing
  labs <- tolower(as.character(spec$pred_damaging_labels %||% character()))
  use_label <- is.na(dmg) & !is.na(label)
  if (length(labs) > 0) dmg[use_label] <- label[use_label] %in% labs
  dmg
}

#' Normalized damaging magnitude in [0,1] (0 = benign end, 1 = damaging end).
#' @keywords internal
predictor_magnitude <- function(score, spec) {
  lo <- spec$score_min; hi <- spec$score_max
  if (is.null(lo) || is.null(hi) || hi == lo) return(rep(NA_real_, length(score)))
  m <- (score - lo) / (hi - lo)
  if (identical(spec$direction, "lower_is_more_damaging")) m <- 1 - m
  pmax(0, pmin(1, m))
}

#' Add the functional-predictor evidence layer to a variant table.
#'
#' Detects predictors via the registry, creates canonical `*_SCORE`/`*_PRED`
#' columns, computes applicability by consequence, and produces the
#' `COMPUTATIONAL_EVIDENCE_*` columns using a family-aware consensus (correlated
#' predictors within a family do not count as independent). Does not touch
#' `filter_status`.
#'
#' @param x variant data.frame (needs `consequence`).
#' @param cfg resolved config (reads `predictors:`).
#' @return `x` with predictor columns; the inventory is attached as
#'   `attr(x, "predictor_inventory")`.
#' @keywords internal
add_predictor_evidence <- function(x, cfg = NULL) {
  registry <- tryCatch(load_predictor_registry(), error = function(e) NULL)
  n <- nrow(x)
  scope <- predictor_variant_scope(if ("consequence" %in% names(x)) x$consequence else rep(NA, n))
  min_applic <- cfg_get(cfg, c("predictors", "minimum_applicable_predictors"), 2)
  min_dmg    <- cfg_get(cfg, c("predictors", "minimum_damaging_predictors"), 2)

  inv <- list()
  # per (variant x family) accumulators
  fam_present <- matrix(FALSE, n, length(PREDICTOR_FAMILIES), dimnames = list(NULL, PREDICTOR_FAMILIES))
  fam_dmg     <- matrix(FALSE, n, length(PREDICTOR_FAMILIES), dimnames = list(NULL, PREDICTOR_FAMILIES))
  fam_ben     <- matrix(FALSE, n, length(PREDICTOR_FAMILIES), dimnames = list(NULL, PREDICTOR_FAMILIES))
  fam_mag     <- matrix(NA_real_, n, length(PREDICTOR_FAMILIES), dimnames = list(NULL, PREDICTOR_FAMILIES))
  avail <- rep(0L, n); dmg_cnt <- rep(0L, n); ben_cnt <- rep(0L, n)
  details <- vector("list", n)
  any_applicable_family <- rep(FALSE, n)

  if (!is.null(registry)) for (pname in names(registry)) {
    spec <- registry[[pname]]
    norm <- normalize_one_predictor(x, spec)
    scopes <- unlist(spec$variant_scope)
    applic <- scope %in% scopes
    any_applicable_family <- any_applicable_family | applic
    # canonical columns
    x[[paste0(toupper(pname), "_SCORE")]] <- if (!is.null(norm)) norm$score else rep(NA_real_, n)
    x[[paste0(toupper(pname), "_PRED")]]  <- if (!is.null(norm)) norm$label else rep(NA_character_, n)
    inv[[pname]] <- build_predictor_inventory_row(pname, spec, norm, applic, scope)
    if (is.null(norm)) next
    fam <- spec$family
    dmg <- predictor_is_damaging(norm$score, norm$label, spec)
    mag <- predictor_magnitude(norm$score, spec)
    use <- applic & norm$present & !is.na(dmg)      # applicable + present + determinable
    avail[use] <- avail[use] + 1L
    dmg_cnt[use & dmg] <- dmg_cnt[use & dmg] + 1L
    ben_cnt[use & !dmg] <- ben_cnt[use & !dmg] + 1L
    fam_present[use, fam] <- TRUE
    fam_dmg[use & dmg, fam] <- TRUE
    fam_ben[use & !dmg, fam] <- TRUE
    fm <- ifelse(use & dmg & !is.na(mag), mag, NA_real_)
    fam_mag[, fam] <- pmax(fam_mag[, fam], fm, na.rm = TRUE)
    for (i in which(use)) details[[i]] <- c(details[[i]],
      sprintf("%s=%s(%s)", pname, ifelse(dmg[i], "D", "N"),
              ifelse(is.na(norm$score[i]), norm$label[i], format(round(norm$score[i], 3)))))
  }
  fam_mag[is.infinite(fam_mag)] <- NA_real_

  n_fam_present <- rowSums(fam_present)
  n_fam_dmg <- rowSums(fam_dmg)
  n_fam_ben <- rowSums(fam_ben)
  # family-aware score: mean over families present of their damaging magnitude (0 if none damaging)
  fm2 <- fam_mag; fm2[is.na(fm2) & fam_present] <- 0
  fm2[!fam_present] <- NA_real_
  comp_score <- rowMeans(fm2, na.rm = TRUE); comp_score[is.nan(comp_score)] <- NA_real_

  conflict <- (dmg_cnt > 0 & ben_cnt > 0)
  applicability <- rep("unsupported_consequence", n)
  applicability[any_applicable_family] <- "applicable_but_missing"
  applicability[any_applicable_family & avail > 0] <- "applicable"
  applicability[!any_applicable_family & scope == "truncating"] <- "not_applicable"  # use consequence directly

  category <- rep("insufficient", n)
  category[!any_applicable_family] <- "not_applicable"
  category[any_applicable_family & avail == 0] <- "insufficient"
  frac <- ifelse(avail > 0, dmg_cnt / avail, NA_real_)
  is_conf <- avail > 0 & conflict & !is.na(frac) & frac >= 0.34 & frac <= 0.66
  category[avail >= 1 & dmg_cnt >= 1 & !is_conf] <- "weak_support"
  category[n_fam_dmg >= min_dmg & !is_conf & n_fam_ben == 0] <- "strong_support"
  category[n_fam_dmg >= min_dmg & !is_conf & n_fam_ben > 0] <- "moderate_support"
  category[is_conf] <- "conflicting"

  x$COMPUTATIONAL_EVIDENCE_SCORE <- round(comp_score, 4)
  x$COMPUTATIONAL_EVIDENCE_CATEGORY <- category
  x$COMPUTATIONAL_EVIDENCE_DETAILS <- vapply(details, function(d) if (is.null(d)) "" else paste(d, collapse = ";"), character(1))
  x$COMPUTATIONAL_PREDICTORS_AVAILABLE <- avail
  x$COMPUTATIONAL_PREDICTORS_DAMAGING <- dmg_cnt
  x$COMPUTATIONAL_PREDICTORS_BENIGN <- ben_cnt
  x$COMPUTATIONAL_PREDICTORS_CONFLICTING <- ifelse(conflict, pmin(dmg_cnt, ben_cnt), 0L)
  x$COMPUTATIONAL_CONFLICT_FLAG <- conflict
  x$COMPUTATIONAL_CONFLICT_REASON <- ifelse(conflict,
    sprintf("%d damaging vs %d benign among applicable predictors", dmg_cnt, ben_cnt), NA_character_)
  x$PREDICTOR_APPLICABILITY_STATUS <- applicability

  attr(x, "predictor_inventory") <- if (length(inv)) do.call(rbind, inv) else
    data.frame(predictor = character(), stringsAsFactors = FALSE)
  x
}

#' @keywords internal
build_predictor_inventory_row <- function(pname, spec, norm, applic, scope) {
  present <- if (is.null(norm)) rep(FALSE, length(applic)) else norm$present
  n_app <- sum(applic)
  score <- if (is.null(norm)) NA_real_ else norm$score
  data.frame(
    predictor = pname,
    detected_column = if (is.null(norm)) NA_character_ else paste(norm$columns, collapse = ";"),
    family = spec$family,
    n_non_missing = if (is.null(norm)) 0L else sum(present),
    fraction_non_missing = if (is.null(norm)) 0 else round(mean(present), 4),
    n_applicable_variants = n_app,
    fraction_applicable_annotated = if (n_app > 0) round(sum(present & applic) / n_app, 4) else 0,
    score_min_observed = if (all(is.na(score))) NA_real_ else min(score, na.rm = TRUE),
    score_max_observed = if (all(is.na(score))) NA_real_ else max(score, na.rm = TRUE),
    n_out_of_range = if (is.null(norm)) 0L else norm$n_out_of_range,
    status = if (is.null(norm)) "not_detected" else if (sum(present) == 0) "detected_all_missing" else "detected",
    stringsAsFactors = FALSE)
}
