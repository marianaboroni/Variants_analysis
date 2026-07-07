add_confidence_ranking <- function(x, cfg) {
  x$somatic_confidence_score <- somatic_confidence_score(x, cfg)
  x$germline_confidence_score <- germline_confidence_score_rank(x, cfg)
  x$polymorphism_confidence_score <- polymorphism_confidence_score(x, cfg)
  x$artifact_confidence_score <- artifact_confidence_score_rank(x, cfg)
  x$manual_review_priority_score <- manual_review_priority_score(x, cfg)

  x$somatic_confidence_tier <- confidence_tier(x$somatic_confidence_score)
  x$germline_confidence_tier <- confidence_tier(x$germline_confidence_score)
  x$polymorphism_confidence_tier <- confidence_tier(x$polymorphism_confidence_score)
  x$artifact_confidence_tier <- confidence_tier(x$artifact_confidence_score)

  x$confidence_bucket <- assign_confidence_bucket(x)
  x$recommended_variant_action <- recommended_variant_action(x)
  x$confidence_rank <- rank_confidence_variants(x)
  x
}

somatic_confidence_score <- function(x, cfg) {
  ml <- if ("ml_true_positive_probability" %in% names(x)) x$ml_true_positive_probability else NA_real_
  ml[is.na(ml)] <- 0.5
  onc <- somatic_oncogenicity_component(x)
  driver <- driver_component(x)
  bounded01(
    0.35 * zero_if_na_rank(x$somatic_score_validated) +
      0.18 * (1 - zero_if_na_rank(x$germline_score)) +
      0.17 * (1 - zero_if_na_rank(x$artifact_score)) +
      0.12 * onc +
      0.08 * driver +
      0.10 * ml
  )
}

germline_confidence_score_rank <- function(x, cfg) {
  acmg <- germline_acmg_component(x)
  bounded01(
    0.40 * zero_if_na_rank(x$germline_score) +
      0.25 * zero_if_na_rank(x$vaf_germline_score) +
      0.20 * population_component(x) +
      0.15 * acmg
  )
}

polymorphism_confidence_score <- function(x, cfg) {
  common <- cfg_get(cfg, c("population_filters", "common_af_threshold"), cfg_get(cfg, c("population_filters", "common_af"), 0.01))
  low <- cfg_get(cfg, c("population_filters", "low_frequency_af_threshold"), cfg_get(cfg, c("population_filters", "rare_af"), 0.001))
  pop <- saturating_score(x$max_pop_af, low, common)
  pop[is.na(pop)] <- 0
  benign <- ifelse(safe_rank_chr(x, "germline_acmg_class") %in% c("Benign", "Likely Benign"), 1, 0)
  recurrence <- saturating_score(x$variant_cohort_freq, cfg_get(cfg, c("cohort", "recurrent_variant_fraction_artifact"), 0.05), 0.20)
  recurrence[is.na(recurrence)] <- 0
  bounded01(0.55 * pop + 0.25 * benign + 0.20 * recurrence)
}

artifact_confidence_score_rank <- function(x, cfg) {
  bounded01(
    0.60 * zero_if_na_rank(x$artifact_score) +
      0.25 * ifelse(safe_rank_chr(x, "artifact_category") == "artifact_strong", 1, 0) +
      0.15 * ifelse(safe_rank_chr(x, "final_class") %in% c("likely_artifact", "technical_fail"), 1, 0)
  )
}

manual_review_priority_score <- function(x, cfg) {
  conflict <- (
    x$somatic_confidence_score >= 0.60 &
      (x$germline_confidence_score >= 0.55 | x$polymorphism_confidence_score >= 0.55 | x$artifact_confidence_score >= 0.55)
  )
  pathogenic_germline <- safe_rank_chr(x, "germline_disposition") == "possible_pathogenic_germline_review"
  rescue <- safe_rank_logical(x, "validation_rescue_candidate")
  uncertain <- safe_rank_chr(x, "final_class") %in% c("manual_review_required", "uncertain_tumor_only")
  ml_prob <- if ("ml_true_positive_probability" %in% names(x)) x$ml_true_positive_probability else rep(NA_real_, nrow(x))
  ml_uncertain <- !is.na(ml_prob) & abs(ml_prob - cfg_get(cfg, c("ml_filter", "true_positive_cutoff"), 0.60)) <= 0.10
  bounded01(
    0.35 * as.numeric(conflict) +
      0.25 * as.numeric(pathogenic_germline) +
      0.20 * as.numeric(uncertain) +
      0.10 * as.numeric(rescue) +
      0.10 * as.numeric(ml_uncertain)
  )
}

assign_confidence_bucket <- function(x) {
  out <- rep("manual_review", nrow(x))
  out[safe_rank_chr(x, "final_class") %in% c("likely_artifact", "technical_fail")] <- "likely_artifact_remove"
  out[safe_rank_chr(x, "final_class") == "likely_germline"] <- "likely_polymorphism_remove"
  out[x$artifact_confidence_score >= 0.70] <- "likely_artifact_remove"
  out[x$polymorphism_confidence_score >= 0.70 & x$germline_confidence_score >= 0.55] <- "likely_polymorphism_remove"
  out[x$germline_confidence_score >= 0.70 & x$polymorphism_confidence_score < 0.70] <- "likely_germline_check"
  out[x$somatic_confidence_score >= 0.80 & x$artifact_confidence_score < 0.50 & x$polymorphism_confidence_score < 0.50] <- "likely_somatic_high_confidence"
  out[x$somatic_confidence_score >= 0.60 & out == "manual_review"] <- "likely_somatic_moderate_confidence"
  out[x$manual_review_priority_score >= 0.50] <- "manual_review_priority"
  out[safe_rank_chr(x, "germline_disposition") == "possible_pathogenic_germline_review"] <- "possible_pathogenic_germline_review"
  out
}

recommended_variant_action <- function(x) {
  bucket <- x$confidence_bucket
  ifelse(bucket %in% c("likely_somatic_high_confidence", "likely_somatic_moderate_confidence"), "keep_somatic_candidate",
         ifelse(bucket %in% c("likely_polymorphism_remove", "likely_artifact_remove"), "remove_from_somatic_list",
                ifelse(bucket %in% c("likely_germline_check", "possible_pathogenic_germline_review", "manual_review_priority"), "manual_review", "manual_review")))
}

rank_confidence_variants <- function(x) {
  priority <- ifelse(x$recommended_variant_action == "keep_somatic_candidate", x$somatic_confidence_score,
                     ifelse(x$recommended_variant_action == "remove_from_somatic_list",
                            pmax(x$polymorphism_confidence_score, x$artifact_confidence_score, na.rm = TRUE),
                            x$manual_review_priority_score))
  rank(-priority, ties.method = "first")
}

confidence_tier <- function(score) {
  ifelse(score >= 0.80, "very_high",
         ifelse(score >= 0.60, "high",
                ifelse(score >= 0.40, "moderate",
                       ifelse(score >= 0.20, "low", "very_low"))))
}

somatic_oncogenicity_component <- function(x) {
  cls <- safe_rank_chr(x, "somatic_oncogenicity_class")
  ifelse(cls == "Oncogenic", 1,
         ifelse(cls == "Likely Oncogenic", 0.8,
                ifelse(cls == "VUS", 0.4, 0)))
}

germline_acmg_component <- function(x) {
  cls <- safe_rank_chr(x, "germline_acmg_class")
  ifelse(cls == "Pathogenic", 1,
         ifelse(cls == "Likely Pathogenic", 0.85,
                ifelse(cls == "VUS", 0.45,
                       ifelse(cls == "Likely Benign", 0.35,
                              ifelse(cls == "Benign", 0.25, 0)))))
}

population_component <- function(x) {
  cat <- safe_rank_chr(x, "population_category")
  ifelse(cat == "population_common", 1,
         ifelse(cat == "population_low_frequency", 0.65,
                ifelse(cat == "population_rare_or_absent", 0.1, 0)))
}

driver_component <- function(x) {
  cls <- safe_rank_chr(x, "driver_class")
  ifelse(cls == "known_driver", 1,
         ifelse(cls == "probable_driver", 0.8,
                ifelse(cls == "possible_driver", 0.55,
                       ifelse(cls == "uncertain_possible_driver", 0.45, 0))))
}

confidence_ranking_table <- function(x) {
  cols <- c(
    "confidence_rank", "recommended_variant_action", "confidence_bucket",
    "sample_id", "tumor_type", "chrom", "pos", "ref", "alt", "gene", "protein_change",
    "final_class", "primary_reason", "driver_class",
    "somatic_confidence_score", "somatic_confidence_tier",
    "germline_confidence_score", "germline_confidence_tier",
    "polymorphism_confidence_score", "polymorphism_confidence_tier",
    "artifact_confidence_score", "artifact_confidence_tier",
    "manual_review_priority_score",
    "somatic_oncogenicity_class", "somatic_oncogenicity_criteria",
    "germline_acmg_class", "germline_acmg_criteria", "germline_disposition",
    "ml_true_positive_probability", "somatic_score_validated", "germline_score",
    "artifact_score", "max_pop_af", "variant_cohort_freq", "variant_tumor_type_freq",
    "dp", "alt_count", "vaf"
  )
  cols <- cols[cols %in% names(x)]
  out <- x[, cols, drop = FALSE]
  out[order(out$confidence_rank), , drop = FALSE]
}

confidence_summary <- function(x) {
  data.table::as.data.table(x)[, .N, by = .(recommended_variant_action, confidence_bucket)][order(recommended_variant_action, confidence_bucket)] |>
    as.data.frame()
}

zero_if_na_rank <- function(x) {
  x[is.na(x)] <- 0
  x
}

safe_rank_chr <- function(x, col) {
  if (!(col %in% names(x))) return(rep(NA_character_, nrow(x)))
  as.character(x[[col]])
}

safe_rank_logical <- function(x, col) {
  if (!(col %in% names(x))) return(rep(FALSE, nrow(x)))
  z <- x[[col]]
  z[is.na(z)] <- FALSE
  as.logical(z)
}
