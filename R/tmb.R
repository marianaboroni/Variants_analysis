calculate_tmb <- function(x, cfg) {
  callable_mb <- cfg_get(cfg, c("tmb", "callable_mb"), NA_real_)
  if (is.null(callable_mb) || is.na(callable_mb)) {
    callable_mb <- assay_default_callable_mb(cfg_get(cfg, c("hard_filters", "assay"), "WGS"))
  }
  use_ml <- isTRUE(cfg_get(cfg, c("tmb", "use_ml_filter"), FALSE))
  ml_cutoff <- cfg_get(cfg, c("tmb", "ml_true_positive_cutoff"), 0.60)

  eligible <- x$final_class %in% c("high_confidence_somatic", "probable_somatic") &
    is_tmb_countable_consequence(x$consequence)
  if (use_ml && "ml_true_positive_probability" %in% names(x)) {
    eligible <- eligible & !is.na(x$ml_true_positive_probability) & x$ml_true_positive_probability >= ml_cutoff
  }

  dt <- data.table::as.data.table(x)
  dt[, tmb_countable := eligible]
  out <- dt[, .(
    n_somatic_variants = sum(final_class %in% c("high_confidence_somatic", "probable_somatic"), na.rm = TRUE),
    n_tmb_countable = sum(tmb_countable, na.rm = TRUE),
    n_known_driver = sum(driver_class == "known_driver", na.rm = TRUE),
    n_probable_driver = sum(driver_class == "probable_driver", na.rm = TRUE),
    median_somatic_vaf = safe_median(vaf[final_class %in% c("high_confidence_somatic", "probable_somatic")])
  ), by = .(sample_id, tumor_type)]
  out[, callable_mb := callable_mb]
  out[, tmb_mut_per_mb := n_tmb_countable / callable_mb]
  out[, tmb_category := tmb_category(tmb_mut_per_mb, cfg)]
  as.data.frame(out)
}

is_tmb_countable_consequence <- function(consequence) {
  z <- tolower(as.character(consequence))
  grepl("missense|frameshift|stop_gained|stop_lost|start_lost|splice_acceptor|splice_donor|inframe|protein_altering", z)
}

assay_default_callable_mb <- function(assay) {
  assay <- toupper(as.character(assay))
  if (assay == "PANEL") return(1)
  if (assay == "WGS") return(30)
  30
}

tmb_category <- function(tmb, cfg) {
  high <- cfg_get(cfg, c("tmb", "high_threshold"), 10)
  intermediate <- cfg_get(cfg, c("tmb", "intermediate_threshold"), 5)
  ifelse(tmb >= high, "TMB_high",
         ifelse(tmb >= intermediate, "TMB_intermediate", "TMB_low"))
}

