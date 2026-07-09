# Adaptive per-sample filtering (EXPERIMENTAL). Complements the universal hard
# floors with per-sample technical profiles, an empirical noise floor, a
# probabilistic alt-read error test, and a three-zone (PASS/REVIEW/FAIL) decision.
#
# Guardrails (enforced): adaptive thresholds NEVER go below the absolute safety
# floors; biological relevance NEVER relaxes technical requirements; borderline
# calls go to REVIEW (not silently PASS); globally poor samples are flagged
# (SAMPLE_QC_STATUS) rather than rescued by lowering thresholds; the module never
# learns TP/FP labels from its own PASS/FAIL output. All thresholds are heuristic
# and assay-dependent (see docs/ADAPTIVE_FILTERING.md).

#' Compute a robust per-sample technical profile.
#' @keywords internal
compute_sample_technical_profile <- function(x, cfg = NULL) {
  by <- split(seq_len(nrow(x)), x$sample_id)
  callable_dp <- cfg_get(cfg, c("adaptive_filtering", "callable_min_dp"),
                         cfg_get(cfg, c("adaptive_filtering", "hard_safety_limits", "minimum_dp"), 8))
  rows <- lapply(names(by), function(s) {
    idx <- by[[s]]
    dp <- suppressWarnings(as.numeric(x$dp[idx])); dp <- dp[!is.na(dp)]
    vaf <- suppressWarnings(as.numeric(x$vaf[idx])); vaf <- vaf[!is.na(vaf)]
    altd <- suppressWarnings(as.numeric(x$alt_count[idx])); altd <- altd[!is.na(altd)]
    q <- function(p) if (length(dp)) as.numeric(stats::quantile(dp, p, names = FALSE)) else NA_real_
    # empirical noise floor: central tendency of the low-VAF tail (likely error/subclonal)
    lowv <- vaf[vaf > 0 & vaf < 0.10]
    bg <- if (length(lowv) >= 5) stats::median(lowv) else NA_real_
    noise_floor <- max(0.001, min(0.02, if (is.na(bg)) 0.005 else bg))
    data.frame(sample_id = s,
      SAMPLE_MEDIAN_DP = safe_median(dp), SAMPLE_MEAN_DP = if (length(dp)) mean(dp) else NA_real_,
      SAMPLE_DP_MAD = if (length(dp)) stats::mad(dp) else NA_real_,
      SAMPLE_DP_IQR = if (length(dp)) stats::IQR(dp) else NA_real_,
      SAMPLE_DP_P05 = q(0.05), SAMPLE_DP_P10 = q(0.10), SAMPLE_DP_P25 = q(0.25),
      SAMPLE_DP_P75 = q(0.75), SAMPLE_DP_P90 = q(0.90), SAMPLE_DP_P95 = q(0.95),
      SAMPLE_MEDIAN_ALT_DEPTH = safe_median(altd),
      SAMPLE_BACKGROUND_VAF = bg, SAMPLE_NOISE_FLOOR = noise_floor,
      SAMPLE_CALLABLE_FRACTION = if (length(dp)) mean(dp >= callable_dp) else NA_real_,
      SAMPLE_LOW_COVERAGE_FRACTION = if (length(dp)) mean(dp < callable_dp) else NA_real_,
      SAMPLE_N_VARIANTS = length(idx), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Derive per-sample adaptive thresholds, never below absolute safety floors.
#' @keywords internal
compute_adaptive_thresholds <- function(profile, cfg = NULL) {
  floor_dp <- cfg_get(cfg, c("adaptive_filtering", "depth", "absolute_floor"),
                      cfg_get(cfg, c("adaptive_filtering", "hard_safety_limits", "minimum_dp"), 8))
  floor_alt <- cfg_get(cfg, c("adaptive_filtering", "alt_depth", "absolute_floor"),
                       cfg_get(cfg, c("adaptive_filtering", "hard_safety_limits", "minimum_alt_depth"), 2))
  method <- cfg_get(cfg, c("adaptive_filtering", "depth", "method"), "percentile_mad")
  pct <- cfg_get(cfg, c("adaptive_filtering", "depth", "percentile"), 0.05)
  madk <- cfg_get(cfg, c("adaptive_filtering", "depth", "mad_multiplier"), 3)
  above_noise <- cfg_get(cfg, c("adaptive_filtering", "vaf", "minimum_above_noise_factor"), 3)

  p <- profile
  adaptive_from_pct <- p[[sprintf("SAMPLE_DP_P%02d", round(pct * 100))]]
  if (is.null(adaptive_from_pct)) adaptive_from_pct <- p$SAMPLE_DP_P05
  adaptive_from_mad <- p$SAMPLE_MEDIAN_DP - madk * p$SAMPLE_DP_MAD
  adaptive_dp <- switch(method,
    percentile = adaptive_from_pct,
    robust_zscore = adaptive_from_mad,
    percentile_mad = pmin(adaptive_from_pct, adaptive_from_mad, na.rm = TRUE),
    adaptive_from_pct)
  adaptive_dp[is.na(adaptive_dp)] <- floor_dp
  p$ABSOLUTE_MIN_DP <- floor_dp
  p$ABSOLUTE_MIN_ALT_DEPTH <- floor_alt
  p$ADAPTIVE_MIN_DP <- pmax(floor_dp, floor(adaptive_dp))          # GUARDRAIL: never below floor
  p$ADAPTIVE_MIN_ALT_DEPTH <- floor_alt                            # per-variant SNR test does the rest
  p$ADAPTIVE_MIN_VAF <- pmax(p$SAMPLE_NOISE_FLOOR * above_noise, 0)
  p$THRESHOLD_METHOD <- method
  p$THRESHOLD_VERSION <- "adaptive_v1"
  # global sample QC (do not rescue a bad sample by lowering thresholds)
  min_callable <- cfg_get(cfg, c("adaptive_filtering", "sample_qc", "min_callable_fraction"), 0.5)
  min_median_dp <- cfg_get(cfg, c("adaptive_filtering", "sample_qc", "min_median_dp"), floor_dp * 2)
  p$SAMPLE_QC_STATUS <- ifelse(!is.na(p$SAMPLE_CALLABLE_FRACTION) & p$SAMPLE_CALLABLE_FRACTION < min_callable * 0.5, "FAIL",
    ifelse((!is.na(p$SAMPLE_CALLABLE_FRACTION) & p$SAMPLE_CALLABLE_FRACTION < min_callable) |
           (!is.na(p$SAMPLE_MEDIAN_DP) & p$SAMPLE_MEDIAN_DP < min_median_dp), "REVIEW", "PASS"))
  p$SAMPLE_QC_WARNINGS <- ifelse(p$SAMPLE_QC_STATUS == "PASS", "",
    paste0("callable_fraction=", round(p$SAMPLE_CALLABLE_FRACTION, 3), ";median_dp=", round(p$SAMPLE_MEDIAN_DP, 1)))
  p
}

#' Add adaptive per-sample filtering columns + a three-zone ADAPTIVE_FILTER_STATUS.
#' Attaches the per-sample threshold manifest as attr(x, "adaptive_thresholds").
#' @keywords internal
add_adaptive_filtering <- function(x, cfg = NULL) {
  profile <- compute_sample_technical_profile(x, cfg)
  thr <- compute_adaptive_thresholds(profile, cfg)
  keepcols <- c("sample_id", "SAMPLE_MEDIAN_DP", "SAMPLE_DP_MAD", "SAMPLE_NOISE_FLOOR",
    "SAMPLE_CALLABLE_FRACTION", "ADAPTIVE_MIN_DP", "ADAPTIVE_MIN_ALT_DEPTH", "ADAPTIVE_MIN_VAF",
    "ABSOLUTE_MIN_DP", "ABSOLUTE_MIN_ALT_DEPTH", "SAMPLE_QC_STATUS")
  x <- merge(x, thr[, keepcols, drop = FALSE], by = "sample_id", all.x = TRUE, sort = FALSE)

  dp <- suppressWarnings(as.numeric(x$dp)); altc <- suppressWarnings(as.numeric(x$alt_count))
  vaf <- suppressWarnings(as.numeric(x$vaf)); nf <- x$SAMPLE_NOISE_FLOOR
  review_margin <- cfg_get(cfg, c("adaptive_filtering", "decision", "review_margin"), 0.15)
  snr_min <- cfg_get(cfg, c("adaptive_filtering", "alt_depth", "minimum_signal_to_noise"), 3)

  x$EXPECTED_ERROR_ALT_READS <- round(zero_na2(dp) * zero_na2(nf), 3)
  x$AD_ALT_SIGNAL_TO_NOISE <- round(zero_na2(altc) / (x$EXPECTED_ERROR_ALT_READS + 1e-6), 3)
  x$ALT_READ_ERROR_PVALUE <- mapply(function(a, d, p) {
    if (is.na(a) || is.na(d) || d <= 0 || is.na(p)) return(NA_real_)
    stats::pbinom(max(0, a - 1), size = round(d), prob = p, lower.tail = FALSE)
  }, altc, dp, nf)
  x$EXPECTED_MIN_DETECTABLE_VAF <- round(pmax(x$ADAPTIVE_MIN_ALT_DEPTH / pmax(dp, 1), zero_na2(nf)), 4)
  x$VAF_ABOVE_SAMPLE_NOISE <- !is.na(vaf) & !is.na(nf) & vaf >= nf * snr_min
  x$VAF_DETECTION_CONFIDENCE <- ifelse(is.na(vaf) | is.na(x$EXPECTED_MIN_DETECTABLE_VAF), "unknown",
    ifelse(vaf >= x$EXPECTED_MIN_DETECTABLE_VAF * 2, "high",
      ifelse(vaf >= x$EXPECTED_MIN_DETECTABLE_VAF, "moderate", "low")))

  bias <- pmax(zero_na2(suppressWarnings(as.numeric(x$orientation_bias))),
               zero_na2(suppressWarnings(as.numeric(x$strand_artifact))))
  pon <- if ("pon_flag" %in% names(x)) x$pon_flag %in% c(TRUE, "TRUE") else rep(FALSE, nrow(x))

  # ---- three-zone decision (technical only) -------------------------------
  fail <- (!is.na(dp) & dp < x$ABSOLUTE_MIN_DP) |
          (!is.na(altc) & altc < x$ABSOLUTE_MIN_ALT_DEPTH) |
          (!is.na(x$ALT_READ_ERROR_PVALUE) & x$ALT_READ_ERROR_PVALUE > 0.05 & !is.na(altc) & altc > 0) |
          (bias > 0.9) | pon
  review <- !fail & (
          (!is.na(dp) & dp < x$ADAPTIVE_MIN_DP) |
          (!is.na(altc) & x$AD_ALT_SIGNAL_TO_NOISE < snr_min) |
          (!is.na(vaf) & !is.na(x$ADAPTIVE_MIN_VAF) & vaf < x$ADAPTIVE_MIN_VAF & x$VAF_ABOVE_SAMPLE_NOISE) |
          (!is.na(dp) & dp < x$ADAPTIVE_MIN_DP * (1 + review_margin)) |
          (bias > 0.8) | (x$SAMPLE_QC_STATUS %in% c("FAIL", "REVIEW")))
  status <- ifelse(fail, "FAIL", ifelse(review, "REVIEW", "PASS"))
  x$ADAPTIVE_FILTER_STATUS <- status
  x$ADAPTIVE_FILTER_REASON <- mapply(function(f, r, d, a, sn, pv, bi, po, sq) {
    z <- character()
    if (isTRUE(f)) {
      if (!is.na(d) && d < 999999) z <- c(z, "below_absolute_or_noise")
      if (isTRUE(po)) z <- c(z, "panel_of_normals")
      if (!is.na(bi) && bi > 0.9) z <- c(z, "strong_bias")
      if (!is.na(pv) && pv > 0.05) z <- c(z, "compatible_with_noise")
    } else if (isTRUE(r)) {
      if (!is.na(sn) && sn < 3) z <- c(z, "low_signal_to_noise")
      if (sq %in% c("FAIL", "REVIEW")) z <- c(z, "sample_qc_" )
      if (length(z) == 0) z <- "near_adaptive_threshold"
    }
    if (length(z) == 0) "" else paste(z, collapse = ";")
  }, fail, review, dp, altc, x$AD_ALT_SIGNAL_TO_NOISE, x$ALT_READ_ERROR_PVALUE, bias, pon, x$SAMPLE_QC_STATUS,
     USE.NAMES = FALSE)

  attr(x, "adaptive_thresholds") <- thr
  x
}

#' Documented adaptive-filtering REVIEW escalation (borderline technical -> REVIEW).
#' Never rescues below absolute floors; configurable.
#' @keywords internal
apply_adaptive_review <- function(x, cfg = NULL) {
  if (!isTRUE(cfg_get(cfg, c("adaptive_filtering", "escalate_to_review"), TRUE))) return(x)
  if (!all(c("ADAPTIVE_FILTER_STATUS", "STATUS") %in% names(x))) return(x)
  esc <- x$ADAPTIVE_FILTER_STATUS == "REVIEW" & x$STATUS == "PASS"
  if (any(esc)) {
    x$STATUS[esc] <- "REVIEW"
    if ("filter_status" %in% names(x)) x$filter_status[esc] <- "REVIEW"
    x$status_reason[esc] <- ifelse(is.na(x$status_reason[esc]) | x$status_reason[esc] == "",
      "adaptive_borderline_technical", paste(x$status_reason[esc], "adaptive_borderline_technical", sep = "; "))
    log_step("adaptive", "adaptive borderline REVIEW escalation applied", escalated = sum(esc))
  }
  x
}

#' @keywords internal
zero_na2 <- function(v) { v[is.na(v)] <- 0; v }
