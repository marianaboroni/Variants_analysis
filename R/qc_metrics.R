build_variant_call_qc_summary <- function(x, sample_qc = NULL) {
  dt <- data.table::as.data.table(x)
  out <- dt[, .(
    n_raw_variants = .N,
    n_hard_filter_fail = sum(final_class == "hard_filter_fail", na.rm = TRUE),
    n_probable_germline = sum(final_class == "probable_germline", na.rm = TRUE),
    n_probable_artifact = sum(final_class == "probable_artifact", na.rm = TRUE),
    n_probable_somatic = sum(final_class == "probable_somatic", na.rm = TRUE),
    n_high_confidence_somatic = sum(final_class == "high_confidence_somatic", na.rm = TRUE),
    n_known_driver = sum(driver_class == "known_driver", na.rm = TRUE),
    n_probable_driver = sum(driver_class == "probable_driver", na.rm = TRUE),
    median_dp = safe_median(dp),
    median_alt_count = safe_median(alt_count),
    median_vaf = safe_median(vaf),
    median_mbq = safe_median(mbq),
    median_mmq = safe_median(mmq),
    recurrent_variant_fraction = mean(!is.na(variant_cohort_freq) & variant_cohort_freq >= 0.05, na.rm = TRUE)
  ), by = .(sample_id, tumor_type)]
  out[, hard_filter_fail_fraction := n_hard_filter_fail / n_raw_variants]
  out[, germline_fraction := n_probable_germline / n_raw_variants]
  out[, artifact_fraction := n_probable_artifact / n_raw_variants]
  out[, somatic_fraction := (n_probable_somatic + n_high_confidence_somatic) / n_raw_variants]
  if (!is.null(sample_qc)) {
    out <- merge(out, data.table::as.data.table(sample_qc), by = c("sample_id", "tumor_type"), all.x = TRUE)
  }
  as.data.frame(out)
}

