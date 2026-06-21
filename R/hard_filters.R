add_sample_qc_and_hard_filters <- function(x, cfg) {
  if (is.null(cfg$hard_filters) || isFALSE(cfg$hard_filters$enabled)) {
    x$hard_min_depth <- cfg$technical_filters$min_depth
    x$hard_min_alt_count_snv <- cfg$technical_filters$min_alt_count_snv
    x$hard_min_alt_count_indel <- cfg$technical_filters$min_alt_count_indel
    x$hard_min_af_snv <- cfg$technical_filters$min_af_snv
    x$hard_min_af_indel <- cfg$technical_filters$min_af_indel
    x$sample_qc_class <- "not_evaluated"
    x$hard_filter_reason <- ifelse(x$technical_pass_floor, "PASS", "base_technical_filter")
    x$hard_filter_pass <- x$technical_pass_floor
    return(list(variants = x, sample_qc = compute_sample_qc_from_variants(x, cfg)))
  }

  sample_qc <- compute_sample_qc_from_variants(x, cfg)
  external_qc_path <- cfg_get(cfg, c("input", "sample_qc"), NULL)
  if (!is.null(external_qc_path) && !is.na(external_qc_path) && file.exists(external_qc_path)) {
    sample_qc <- merge_external_sample_qc(sample_qc, external_qc_path)
  }
  sample_qc <- derive_sample_thresholds(sample_qc, cfg)
  x <- merge(x, sample_qc, by = c("sample_id", "tumor_type"), all.x = TRUE)
  x <- apply_variant_hard_filters(x, cfg)
  list(variants = x, sample_qc = sample_qc)
}

compute_sample_qc_from_variants <- function(x, cfg) {
  split_idx <- split(seq_len(nrow(x)), x$sample_id)
  out <- lapply(names(split_idx), function(sample_id) {
    idx <- split_idx[[sample_id]]
    dp <- x$dp[idx]
    filter_status <- x$filter_status[idx]
    pass_filter <- is.na(filter_status) | filter_status %in% c("PASS", ".", "NA") | filter_status == ""
    data.frame(
      sample_id = sample_id,
      tumor_type = unique(x$tumor_type[idx])[1],
      sample_n_variants = length(idx),
      sample_variant_median_dp = safe_median(dp),
      sample_variant_q10_dp = safe_quantile(dp, 0.10),
      sample_variant_q25_dp = safe_quantile(dp, 0.25),
      sample_variant_q75_dp = safe_quantile(dp, 0.75),
      sample_pass_filter_fraction = mean(pass_filter, na.rm = TRUE),
      sample_missing_dp_fraction = mean(is.na(dp)),
      stringsAsFactors = FALSE
    )
  })
  sample_qc <- do.call(rbind, out)
  sample_qc$external_median_depth <- NA_real_
  sample_qc$external_mean_depth <- NA_real_
  sample_qc$external_pct_20x <- NA_real_
  sample_qc$external_pct_30x <- NA_real_
  sample_qc$external_pct_100x <- NA_real_
  sample_qc$external_contamination <- NA_real_
  sample_qc$external_duplicate_rate <- NA_real_
  sample_qc
}

merge_external_sample_qc <- function(sample_qc, path) {
  ext <- read_variants(path, "\t")
  ext$sample_id <- as.character(coalesce_columns(
    ext,
    c("Tumor_Sample_Barcode", "Sample_Barcode", "Tumor_Sample", "Sample", "sample", "sample_id")
  ))
  ext$external_median_depth <- to_numeric_safe(coalesce_columns(
    ext,
    c("median_depth", "Median_Depth", "MEAN_COVERAGE", "mean_coverage", "coverage_median")
  ))
  ext$external_mean_depth <- to_numeric_safe(coalesce_columns(
    ext,
    c("mean_depth", "Mean_Depth", "mean_coverage", "MEAN_COVERAGE")
  ))
  ext$external_pct_20x <- normalize_fraction(to_numeric_safe(coalesce_columns(
    ext,
    c("pct_20x", "PCT_20X", "pct_bases_20x", "coverage_20x")
  )))
  ext$external_pct_30x <- normalize_fraction(to_numeric_safe(coalesce_columns(
    ext,
    c("pct_30x", "PCT_30X", "pct_bases_30x", "coverage_30x")
  )))
  ext$external_pct_100x <- normalize_fraction(to_numeric_safe(coalesce_columns(
    ext,
    c("pct_100x", "PCT_100X", "pct_bases_100x", "coverage_100x")
  )))
  ext$external_contamination <- normalize_fraction(to_numeric_safe(coalesce_columns(
    ext,
    c("contamination", "CONTAMINATION", "contamination_fraction")
  )))
  ext$external_duplicate_rate <- normalize_fraction(to_numeric_safe(coalesce_columns(
    ext,
    c("duplicate_rate", "PERCENT_DUPLICATION", "duplication_rate")
  )))
  keep <- unique(ext[, c(
    "sample_id",
    "external_median_depth",
    "external_mean_depth",
    "external_pct_20x",
    "external_pct_30x",
    "external_pct_100x",
    "external_contamination",
    "external_duplicate_rate"
  )])
  merge(sample_qc, keep, by = "sample_id", all.x = TRUE)
}

derive_sample_thresholds <- function(sample_qc, cfg) {
  profile <- assay_profile(cfg)
  tf <- cfg$technical_filters
  hf <- cfg$hard_filters
  depth_fraction <- cfg_get(cfg, c("hard_filters", "adaptive_depth_fraction"), 0.35)
  q25_fraction <- cfg_get(cfg, c("hard_filters", "adaptive_q25_depth_fraction"), 0.50)

  external_depth <- coalesce_numeric(sample_qc$external_median_depth, sample_qc$external_mean_depth)
  observed_depth <- coalesce_numeric(external_depth, sample_qc$sample_variant_median_dp)
  observed_q25 <- coalesce_numeric(sample_qc$sample_variant_q25_dp, observed_depth)

  adaptive_from_depth <- floor(observed_depth * depth_fraction)
  adaptive_from_q25 <- floor(observed_q25 * q25_fraction)
  adaptive_depth <- pmax(tf$min_depth, adaptive_from_depth, adaptive_from_q25, na.rm = TRUE)
  adaptive_depth[is.infinite(adaptive_depth)] <- tf$min_depth
  adaptive_depth <- pmin(adaptive_depth, profile$max_adaptive_min_depth)

  sample_qc$assay <- profile$assay
  sample_qc$observed_depth_for_filter <- observed_depth
  sample_qc$hard_min_depth <- adaptive_depth
  sample_qc$hard_min_alt_count_snv <- pmax(
    tf$min_alt_count_snv,
    ceiling(sample_qc$hard_min_depth * tf$min_af_snv)
  )
  sample_qc$hard_min_alt_count_indel <- pmax(
    tf$min_alt_count_indel,
    ceiling(sample_qc$hard_min_depth * tf$min_af_indel)
  )
  sample_qc$hard_min_af_snv <- tf$min_af_snv
  sample_qc$hard_min_af_indel <- tf$min_af_indel

  callable_fraction <- choose_callable_fraction(sample_qc, profile)
  median_depth_for_qc <- coalesce_numeric(sample_qc$external_median_depth, sample_qc$sample_variant_median_dp)
  sample_qc$sample_qc_reason <- "PASS"
  sample_qc$sample_qc_class <- "pass"

  low_depth <- !is.na(median_depth_for_qc) & median_depth_for_qc < profile$min_median_depth
  low_callable <- !is.na(callable_fraction) & callable_fraction < profile$min_callable_fraction
  high_contam <- !is.na(sample_qc$external_contamination) &
    sample_qc$external_contamination > profile$max_contamination
  high_dup <- !is.na(sample_qc$external_duplicate_rate) &
    sample_qc$external_duplicate_rate > profile$max_duplicate_rate

  warn_depth <- !is.na(median_depth_for_qc) & median_depth_for_qc < profile$warn_median_depth
  sample_qc$sample_qc_class[warn_depth] <- "warn"
  sample_qc$sample_qc_reason[warn_depth] <- "low_depth_warning"
  fail <- low_depth | low_callable | high_contam | high_dup
  sample_qc$sample_qc_class[fail] <- "fail"
  fail_reasons <- paste_reasons(
    low_depth, "low_median_depth",
    low_callable, "low_callable_fraction",
    high_contam, "high_contamination",
    high_dup, "high_duplicate_rate"
  )
  sample_qc$sample_qc_reason[fail] <- fail_reasons[fail]
  sample_qc
}

apply_variant_hard_filters <- function(x, cfg) {
  tf <- cfg$technical_filters
  missing_required_fail <- isTRUE(cfg_get(cfg, c("hard_filters", "missing_required_metrics_fail"), FALSE))
  min_alt <- ifelse(x$is_indel, x$hard_min_alt_count_indel, x$hard_min_alt_count_snv)
  min_af <- ifelse(x$is_indel, x$hard_min_af_indel, x$hard_min_af_snv)
  pass_filter <- is.na(x$filter_status) | x$filter_status %in% c("PASS", ".", "NA") | x$filter_status == ""

  depth_fail <- ifelse(is.na(x$dp), missing_required_fail, x$dp < x$hard_min_depth)
  alt_fail <- ifelse(is.na(x$alt_count), missing_required_fail, x$alt_count < min_alt)
  af_fail <- ifelse(is.na(x$vaf), missing_required_fail, x$vaf < min_af)
  tlod_fail <- ifelse(is.na(x$tlod), missing_required_fail, x$tlod < tf$min_tlod)
  mbq_fail <- ifelse(is.na(x$mbq), missing_required_fail, x$mbq < tf$min_mbq)
  mmq_fail <- ifelse(is.na(x$mmq), missing_required_fail, x$mmq < tf$min_mmq)
  sample_fail <- x$sample_qc_class == "fail"
  filter_fail <- !pass_filter

  x$hard_min_alt_count <- min_alt
  x$hard_min_af <- min_af
  x$hard_filter_reason <- paste_reasons(
    depth_fail, "low_depth",
    alt_fail, "low_alt_count",
    af_fail, "low_vaf",
    tlod_fail, "low_tlod",
    mbq_fail, "low_mbq",
    mmq_fail, "low_mmq",
    filter_fail, "caller_filter_not_pass",
    sample_fail, "sample_qc_fail"
  )
  x$hard_filter_pass <- x$hard_filter_reason == "PASS"
  x$technical_pass_floor <- x$hard_filter_pass
  x
}

assay_profile <- function(cfg) {
  assay <- toupper(as.character(cfg_get(cfg, c("hard_filters", "assay"), "WGS")))
  profiles <- list(
    WGS = list(
      assay = "WGS",
      min_median_depth = 25,
      warn_median_depth = 30,
      min_callable_fraction = 0.70,
      callable_column = "external_pct_20x",
      max_contamination = 0.05,
      max_duplicate_rate = 0.35,
      max_adaptive_min_depth = 60
    ),
    WES = list(
      assay = "WES",
      min_median_depth = 40,
      warn_median_depth = 50,
      min_callable_fraction = 0.75,
      callable_column = "external_pct_20x",
      max_contamination = 0.05,
      max_duplicate_rate = 0.45,
      max_adaptive_min_depth = 100
    ),
    PANEL = list(
      assay = "PANEL",
      min_median_depth = 250,
      warn_median_depth = 350,
      min_callable_fraction = 0.80,
      callable_column = "external_pct_100x",
      max_contamination = 0.03,
      max_duplicate_rate = 0.60,
      max_adaptive_min_depth = 500
    )
  )
  profile <- profiles[[assay]]
  if (is.null(profile)) profile <- profiles$WGS

  overrides <- cfg_get(cfg, c("hard_filters", "sample_qc"), list())
  for (nm in names(overrides)) {
    profile[[nm]] <- overrides[[nm]]
  }
  profile
}

choose_callable_fraction <- function(sample_qc, profile) {
  col <- profile$callable_column
  if (!is.null(col) && col %in% names(sample_qc)) {
    return(sample_qc[[col]])
  }
  rep(NA_real_, nrow(sample_qc))
}

safe_median <- function(x) {
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

safe_quantile <- function(x, p) {
  if (all(is.na(x))) NA_real_ else as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

normalize_fraction <- function(x) {
  x <- ifelse(x > 1 & x <= 100, x / 100, x)
  x
}

coalesce_numeric <- function(...) {
  vals <- list(...)
  out <- vals[[1]]
  if (length(vals) == 1) return(out)
  for (z in vals[-1]) {
    idx <- is.na(out)
    out[idx] <- z[idx]
  }
  out
}

paste_reasons <- function(...) {
  args <- list(...)
  flags <- args[seq(1, length(args), by = 2)]
  labels <- args[seq(2, length(args), by = 2)]
  n <- length(flags[[1]])
  out <- rep("PASS", n)
  for (i in seq_len(n)) {
    hit <- character()
    for (j in seq_along(flags)) {
      flag <- flags[[j]][i]
      if (!is.na(flag) && isTRUE(flag)) hit <- c(hit, labels[[j]])
    }
    if (length(hit) > 0) out[[i]] <- paste(hit, collapse = ";")
  }
  out
}
