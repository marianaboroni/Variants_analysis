# Tumor mutational burden / variant-burden auxiliary module.
#
# This module is intentionally conservative: it only reports TMB in mutations/Mb
# when an explicit callable territory is configured. Otherwise it still reports
# countable candidate counts as a burden summary, but marks TMB as not evaluable.

#' Calculate tumor mutational burden and countable variant burden.
#'
#' @param variants Classified tumoronly variant table.
#' @param cfg Resolved tumoronly configuration.
#' @return A list with updated `variants` and per-sample `summary`.
#' @export
calculate_tmb <- function(variants, cfg = NULL) {
  if (!is.data.frame(variants)) stop("calculate_tmb(): `variants` must be a data.frame.", call. = FALSE)
  enabled <- isTRUE(cfg_get(cfg, c("tmb", "enabled"), TRUE))
  if (!enabled) {
    variants$tmb_countable <- FALSE
    return(list(variants = variants, summary = tmb_disabled_summary(variants)))
  }

  callable_mb <- suppressWarnings(as.numeric(cfg_get(cfg, c("tmb", "callable_mb"), NA_real_)))
  require_callable <- isTRUE(cfg_get(cfg, c("tmb", "require_callable_mb"), TRUE))
  has_callable <- length(callable_mb) == 1L && !is.na(callable_mb) && is.finite(callable_mb) && callable_mb > 0
  if (!has_callable && !require_callable) {
    callable_mb <- assay_default_callable_mb(cfg_get(cfg, c("hard_filters", "assay"), "WGS"))
    has_callable <- TRUE
  }

  somatic_classes <- cfg_get(cfg, c("tmb", "include_final_classes"),
                             c("high_confidence_somatic", "probable_somatic"))
  use_ml <- isTRUE(cfg_get(cfg, c("tmb", "use_ml_filter"), FALSE))
  ml_cutoff <- cfg_get(cfg, c("tmb", "ml_true_positive_cutoff"), 0.60)

  eligible <- variants$final_class %in% somatic_classes &
    is_tmb_countable_consequence(variants$consequence)
  if (use_ml && "ml_true_positive_probability" %in% names(variants)) {
    eligible <- eligible & !is.na(variants$ml_true_positive_probability) &
      variants$ml_true_positive_probability >= ml_cutoff
  }
  eligible[is.na(eligible)] <- FALSE
  variants$tmb_countable <- eligible

  dt <- data.table::as.data.table(variants)
  out <- dt[, .(
    module_status = if (has_callable) "evaluated" else "not_evaluable_missing_callable_mb",
    n_somatic_candidate_variants = sum(final_class %in% somatic_classes, na.rm = TRUE),
    n_tmb_countable = sum(tmb_countable, na.rm = TRUE),
    n_known_driver = sum(driver_class == "known_driver", na.rm = TRUE),
    n_probable_driver = sum(driver_class == "probable_driver", na.rm = TRUE),
    median_somatic_vaf = safe_median(vaf[final_class %in% somatic_classes])
  ), by = .(sample_id, tumor_type)]
  out[, callable_mb := if (has_callable) callable_mb else NA_real_]
  out[, tmb_mut_per_mb := if (has_callable) n_tmb_countable / callable_mb else NA_real_]
  out[, tmb_category := if (has_callable) tmb_category(tmb_mut_per_mb, cfg) else "not_evaluable"]
  out[, burden_label := if (has_callable) "tumor_mutational_burden" else "countable_variant_burden_only"]
  out[, interpretation := if (has_callable)
    "TMB is reported using the configured callable territory; interpret as tumor-only exploratory output."
    else "Callable territory is missing, so only countable candidate burden is reported; do not interpret as TMB."]
  out[, limitation := paste(
    "Tumor-only TMB is exploratory and depends on upstream calling, filtering, consequence annotation,",
    "and a correct callable territory denominator."
  )]
  list(variants = variants, summary = as.data.frame(out))
}

is_tmb_countable_consequence <- function(consequence) {
  z <- tolower(as.character(consequence))
  grepl("missense|frameshift|stop_gained|stop_lost|start_lost|splice_acceptor|splice_donor|inframe|protein_altering", z)
}

assay_default_callable_mb <- function(assay) {
  assay <- toupper(as.character(assay))
  if (assay == "PANEL") return(1)
  if (assay %in% c("WES", "EXOME")) return(30)
  if (assay == "WGS") return(30)
  30
}

tmb_category <- function(tmb, cfg) {
  high <- cfg_get(cfg, c("tmb", "high_threshold"), 10)
  intermediate <- cfg_get(cfg, c("tmb", "intermediate_threshold"), 5)
  ifelse(tmb >= high, "TMB_high",
         ifelse(tmb >= intermediate, "TMB_intermediate", "TMB_low"))
}

tmb_disabled_summary <- function(variants) {
  samples <- module_sample_frame(variants)
  samples$module_status <- "disabled"
  samples$n_somatic_candidate_variants <- NA_integer_
  samples$n_tmb_countable <- NA_integer_
  samples$n_known_driver <- NA_integer_
  samples$n_probable_driver <- NA_integer_
  samples$median_somatic_vaf <- NA_real_
  samples$callable_mb <- NA_real_
  samples$tmb_mut_per_mb <- NA_real_
  samples$tmb_category <- "disabled"
  samples$burden_label <- "disabled"
  samples$interpretation <- "TMB module disabled."
  samples$limitation <- "No TMB or burden result was calculated."
  samples
}

module_sample_frame <- function(variants) {
  if (!is.data.frame(variants) || nrow(variants) == 0 || !"sample_id" %in% names(variants)) {
    return(data.frame(sample_id = character(), tumor_type = character(), stringsAsFactors = FALSE))
  }
  tt <- if ("tumor_type" %in% names(variants)) variants$tumor_type else rep(NA_character_, nrow(variants))
  unique(data.frame(sample_id = as.character(variants$sample_id),
                    tumor_type = as.character(tt), stringsAsFactors = FALSE))
}
