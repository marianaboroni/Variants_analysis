score_variants <- function(x, cfg) {
  x$technical_evidence_score <- technical_evidence_score(x, cfg)
  x$population_germline_score <- population_germline_score(x, cfg)
  x$vaf_germline_score <- vaf_germline_score(x)
  x$cohort_artifact_score <- cohort_artifact_score(x, cfg)
  x$quality_artifact_score <- quality_artifact_score(x, cfg)

  x$germline_score <- bounded01(
    0.55 * x$population_germline_score +
      0.30 * x$vaf_germline_score +
      0.15 * recurrent_germline_signal(x, cfg)
  )
  x$artifact_score <- bounded01(
    0.45 * x$cohort_artifact_score +
      0.45 * x$quality_artifact_score +
      0.10 * ifelse(x$pon_flag, 1, 0)
  )
  x$somatic_score <- bounded01(
    0.45 * x$technical_evidence_score +
      0.25 * (1 - x$germline_score) +
      0.20 * (1 - x$artifact_score) +
      0.10 * pmin(x$impact_rank / 3, 1)
  )
  x
}

technical_evidence_score <- function(x, cfg) {
  tf <- cfg$technical_filters
  min_depth <- metric_threshold(x, "hard_min_depth", tf$min_depth)
  min_alt_snv <- metric_threshold(x, "hard_min_alt_count_snv", tf$min_alt_count_snv)
  min_alt_indel <- metric_threshold(x, "hard_min_alt_count_indel", tf$min_alt_count_indel)
  min_af_snv <- metric_threshold(x, "hard_min_af_snv", tf$min_af_snv)
  min_af_indel <- metric_threshold(x, "hard_min_af_indel", tf$min_af_indel)
  min_alt <- ifelse(x$is_indel, min_alt_indel, min_alt_snv)
  min_af <- ifelse(x$is_indel, min_af_indel, min_af_snv)

  depth_s <- saturating_score(x$dp, min_depth, min_depth * 4)
  alt_s <- saturating_score(x$alt_count, min_alt, min_alt * 4)
  vaf_s <- saturating_score(x$vaf, min_af, 0.35)
  tlod_s <- saturating_score(x$tlod, tf$min_tlod, 20)
  mbq_s <- saturating_score(x$mbq, tf$min_mbq, 35)
  mmq_s <- saturating_score(x$mmq, tf$min_mmq, 60)

  bounded01(rowMeans(cbind(depth_s, alt_s, vaf_s, tlod_s, mbq_s, mmq_s), na.rm = TRUE))
}

population_germline_score <- function(x, cfg) {
  common <- cfg_get(cfg, c("population_filters", "common_af_threshold"), cfg_get(cfg, c("population_filters", "common_af"), 0.01))
  rare <- cfg_get(cfg, c("population_filters", "low_frequency_af_threshold"), cfg_get(cfg, c("population_filters", "rare_af"), 0.001))
  af <- x$max_pop_af
  out <- rep(0, length(af))
  out[!is.na(af) & af >= common] <- 1
  mid <- !is.na(af) & af >= rare & af < common
  out[mid] <- 0.5 + 0.5 * ((af[mid] - rare) / (common - rare))
  out
}

vaf_germline_score <- function(x) {
  score_05 <- ifelse(!is.na(x$distance_vaf_05) & x$distance_vaf_05 <= 0.08, 0.85,
                     ifelse(!is.na(x$distance_vaf_05) & x$distance_vaf_05 <= 0.15, 0.45, 0))
  score_10 <- ifelse(!is.na(x$distance_vaf_10) & x$distance_vaf_10 <= 0.08, 0.75, 0)
  high_depth <- ifelse(!is.na(x$dp) & x$dp >= 30, 1, 0.75)
  pmax(score_05, score_10) * high_depth
}

recurrent_germline_signal <- function(x, cfg) {
  freq <- x$variant_cohort_freq
  out <- rep(0, length(freq))
  out[!is.na(freq) & freq >= cfg$cohort$recurrent_variant_fraction_artifact] <- 0.8
  out[!is.na(freq) & freq >= 0.15] <- 1
  out
}

cohort_artifact_score <- function(x, cfg) {
  vf <- x$variant_cohort_freq
  lf <- x$locus_cohort_freq
  a <- saturating_score(vf, cfg$cohort$recurrent_variant_fraction_artifact, 0.20)
  b <- saturating_score(lf, cfg$cohort$recurrent_locus_fraction_artifact, 0.30)
  low_median_vaf <- ifelse(!is.na(x$variant_median_vaf) & x$variant_median_vaf < 0.08, 0.2, 0)
  bounded01(pmax(a, b) + low_median_vaf)
}

quality_artifact_score <- function(x, cfg) {
  tf <- cfg$technical_filters
  min_depth <- metric_threshold(x, "hard_min_depth", tf$min_depth)
  min_alt_snv <- metric_threshold(x, "hard_min_alt_count_snv", tf$min_alt_count_snv)
  min_alt_indel <- metric_threshold(x, "hard_min_alt_count_indel", tf$min_alt_count_indel)
  min_alt <- ifelse(x$is_indel, min_alt_indel, min_alt_snv)

  bad_depth <- inverse_score(x$dp, min_depth, min_depth * 0.5)
  bad_alt <- inverse_score(x$alt_count, min_alt, 1)
  bad_tlod <- inverse_score(x$tlod, tf$min_tlod, 2)
  bad_mbq <- inverse_score(x$mbq, tf$min_mbq, 15)
  bad_mmq <- inverse_score(x$mmq, tf$min_mmq, 20)
  failed_filter <- !(is.na(x$filter_status) | x$filter_status %in% c("PASS", ".", "NA") | x$filter_status == "")

  bounded01(rowMeans(cbind(bad_depth, bad_alt, bad_tlod, bad_mbq, bad_mmq), na.rm = TRUE) +
              ifelse(failed_filter, 0.4, 0))
}

saturating_score <- function(value, low, high) {
  out <- (value - low) / (high - low)
  out[is.na(out)] <- NA_real_
  bounded01(out)
}

inverse_score <- function(value, good, bad) {
  out <- (good - value) / (good - bad)
  out[is.na(out)] <- NA_real_
  bounded01(out)
}

bounded01 <- function(x) {
  x <- ifelse(is.nan(x), NA_real_, x)
  pmax(0, pmin(1, x))
}

metric_threshold <- function(x, col, fallback) {
  if (!(col %in% names(x))) return(rep(fallback, nrow(x)))
  out <- x[[col]]
  out[is.na(out)] <- fallback
  out
}

classify_variants <- function(x, cfg) {
  sc <- cfg$scoring
  validation_support <- ifelse(is.na(x$validation_support_score), 0, x$validation_support_score)
  x$somatic_score_validated <- bounded01(x$somatic_score + 0.10 * validation_support)
  hard_pass <- if ("hard_filter_pass" %in% names(x)) x$hard_filter_pass else rep(TRUE, nrow(x))

  x$population_category <- classify_population_germline_evidence(x, cfg)
  x$artifact_category <- classify_technical_artifact_evidence(x, cfg)
  x$recurrence_category <- classify_internal_recurrence(x, cfg)
  x$oncogenic_category <- classify_oncogenic_evidence(x, cfg)

  x$final_class <- "uncertain_tumor_only"
  x$primary_reason <- NA_character_
  x$secondary_reasons <- NA_character_

  x$final_class[!hard_pass] <- "technical_fail"
  x$primary_reason[!hard_pass] <- "hard_filter_fail"

  strong_artifact <- x$artifact_category == "artifact_strong" | x$pon_flag
  x$final_class[!hard_pass & strong_artifact] <- "likely_artifact"
  x$primary_reason[!hard_pass & strong_artifact] <- "strong_artifact_or_pon"

  population_common <- x$population_category == "population_common"
  x$final_class[population_common & hard_pass] <- "likely_germline"
  x$primary_reason[population_common & hard_pass] <- "population_common"
  x$final_class[population_common & !hard_pass] <- "technical_fail"
  x$primary_reason[population_common & !hard_pass] <- "technical_fail_with_population_common"

  population_low <- x$population_category == "population_low_frequency"
  vaf_near_germline <- !is.na(x$distance_vaf_05) & x$distance_vaf_05 <= 0.10
  x$final_class[population_low & vaf_near_germline & hard_pass] <- "likely_germline"
  x$primary_reason[population_low & vaf_near_germline & hard_pass] <- "population_low_vaf_germline_like"
  x$final_class[population_low & !vaf_near_germline & hard_pass] <- "manual_review_required"
  x$primary_reason[population_low & !vaf_near_germline & hard_pass] <- "population_low_conflict"
  x$final_class[population_low & !hard_pass] <- "technical_fail"
  x$primary_reason[population_low & !hard_pass] <- "technical_fail_with_population_low"

  artifact_possible <- x$artifact_category == "artifact_possible"
  x$final_class[artifact_possible & x$final_class %in% c("uncertain_tumor_only", "manual_review_required")] <- "manual_review_required"
  x$primary_reason[artifact_possible & x$final_class == "manual_review_required"] <- "possible_artifact"

  artifact_not_detected <- x$artifact_category == "artifact_not_detected"
  somatic_support <- x$oncogenic_category == "oncogenic_exact" & x$population_category == "population_rare_or_absent"
  probable_support <- x$technical_evidence_score >= 0.55 & x$population_category == "population_rare_or_absent" & x$artifact_category %in% c("artifact_not_detected", "artifact_uninformative")

  x$final_class[somatic_support & hard_pass] <- "high_confidence_somatic"
  x$primary_reason[somatic_support & hard_pass] <- "oncogenic_exact_rare_good_quality"

  x$final_class[!somatic_support & probable_support & hard_pass] <- "probable_somatic"
  x$primary_reason[!somatic_support & probable_support & hard_pass] <- "probable_somatic_rule"

  high_conflict <- x$final_class %in% c("high_confidence_somatic", "probable_somatic") & (population_common | x$artifact_category == "artifact_strong" | x$pon_flag)
  x$final_class[high_conflict] <- "manual_review_required"
  x$primary_reason[high_conflict] <- "conflicting_evidence_with_somatic"

  uncertain_conflict <- x$final_class == "uncertain_tumor_only" & (artifact_possible | (!population_low & x$oncogenic_category == "oncogenic_exact" & x$population_category != "population_rare_or_absent"))
  x$final_class[uncertain_conflict] <- "manual_review_required"
  x$primary_reason[uncertain_conflict] <- "uncertain_conflict_requires_review"

  x$primary_reason[x$final_class == "manual_review_required" & is.na(x$primary_reason)] <- "needs_manual_review"

  x$rule_based_class <- x$final_class
  x$evidence_for_somatic <- ifelse(x$final_class %in% c("high_confidence_somatic", "probable_somatic"), "rule_support", "")
  x$evidence_against_somatic <- ifelse(x$final_class %in% c("likely_germline", "likely_artifact", "technical_fail"), "rule_against", "")
  x$suggested_label <- map_final_class_to_suggested_label(x$final_class)
  x$recommended_action <- map_final_class_to_action(x$final_class)

  x
}

classify_population_germline_evidence <- function(x, cfg) {
  common <- cfg_get(cfg, c("population_filters", "common_af_threshold"), cfg_get(cfg, c("population_filters", "common_af"), 0.01))
  low <- cfg_get(cfg, c("population_filters", "low_frequency_af_threshold"), cfg_get(cfg, c("population_filters", "rare_af"), 0.001))
  max_af <- pmax(x$max_pop_af, x$abraom_af, x$sabe_af, x$gnomad_af, na.rm = TRUE)
  max_af[is.infinite(max_af)] <- NA_real_

  category <- rep("population_uninformative", nrow(x))
  category[!is.na(max_af) & max_af >= common] <- "population_common"
  category[!is.na(max_af) & max_af >= low & max_af < common] <- "population_low_frequency"
  category[!is.na(max_af) & max_af < low] <- "population_rare_or_absent"

  has_coverage_info <- !is.na(x$dp) & !is.na(x$alt_count)
  category[is.na(max_af) & has_coverage_info] <- "population_rare_or_absent"
  category
}

classify_technical_artifact_evidence <- function(x, cfg) {
  tf <- cfg$technical_filters
  low_tlod <- !is.na(x$tlod) & x$tlod < tf$min_tlod
  low_mbq <- !is.na(x$mbq) & x$mbq < tf$min_mbq
  low_mmq <- !is.na(x$mmq) & x$mmq < tf$min_mmq
  low_depth <- !is.na(x$dp) & x$dp < x$hard_min_depth
  low_alt <- !is.na(x$alt_count) & x$alt_count < ifelse(x$is_indel, x$hard_min_alt_count_indel, x$hard_min_alt_count_snv)
  filter_fail <- !is.na(x$mutect_filter) & !(x$mutect_filter %in% c("PASS", ".", "NA", ""))
  strong_bias <- !is.na(x$orientation_bias) & x$orientation_bias > 0.8
  strand_bias <- !is.na(x$strand_artifact) & x$strand_artifact > 0.8
  clustered <- x$clustered_events
  weak <- x$weak_evidence
  pon <- x$pon_flag

  strong_artifact <- filter_fail | pon | strong_bias | strand_bias | clustered | weak | low_tlod | low_mbq | low_mmq
  possible_artifact <- (!strong_artifact & (low_depth | low_alt | x$artifact_score >= 0.6 | x$variant_cohort_freq >= cfg$cohort$recurrent_variant_fraction_artifact))
  high_confidence <- !strong_artifact & !possible_artifact & !filter_fail &
    (is.na(x$dp) | x$dp >= x$hard_min_depth) &
    (is.na(x$alt_count) | x$alt_count >= ifelse(x$is_indel, x$hard_min_alt_count_indel, x$hard_min_alt_count_snv)) &
    (is.na(x$tlod) | x$tlod >= tf$min_tlod)
  uninformative <- is.na(x$dp) | is.na(x$alt_count) | is.na(x$tlod)

  category <- rep("artifact_uninformative", nrow(x))
  category[strong_artifact] <- "artifact_strong"
  category[possible_artifact & !strong_artifact] <- "artifact_possible"
  category[high_confidence] <- "artifact_not_detected"
  category[uninformative & !strong_artifact & !possible_artifact] <- "artifact_uninformative"
  category
}

classify_internal_recurrence <- function(x, cfg) {
  variant_artifact <- !is.na(x$variant_cohort_freq) & x$variant_cohort_freq >= cfg$cohort$recurrent_variant_fraction_artifact
  locus_artifact <- !is.na(x$locus_cohort_freq) & x$locus_cohort_freq >= cfg$cohort$recurrent_locus_fraction_artifact
  tumor_supported <- !is.na(x$variant_tumor_type_freq) & x$variant_tumor_type_freq >= 0.05 & x$variant_tumor_type_freq >= x$variant_cohort_freq
  gene_only <- !is.na(x$variant_cohort_freq) & x$variant_cohort_freq >= 0.02

  category <- rep("recurrence_non_informative", nrow(x))
  category[variant_artifact | locus_artifact] <- "recurrence_artifact_suspected"
  category[tumor_supported & !variant_artifact & !locus_artifact] <- "recurrence_tumor_type_supported"
  category[gene_only & category == "recurrence_non_informative"] <- "recurrence_non_informative"
  category
}

classify_oncogenic_evidence <- function(x, cfg) {
  oncogenic <- rep("oncogenic_none", nrow(x))
  exact_oncogenic <- x$hotspot_match & x$hotspot_tumor_specific
  supportive_oncogenic <- x$hotspot_match & !x$hotspot_tumor_specific
  gene_only <- x$oncokb_match & !x$hotspot_match

  oncogenic[exact_oncogenic] <- "oncogenic_exact"
  oncogenic[supportive_oncogenic] <- "oncogenic_supportive"
  oncogenic[gene_only] <- "oncogenic_gene_only"
  oncogenic[!is.na(x$cosmic_match) & x$cosmic_match & oncogenic == "oncogenic_none"] <- "oncogenic_supportive"
  oncogenic
}

map_final_class_to_suggested_label <- function(final_class) {
  ifelse(final_class %in% c("high_confidence_somatic", "probable_somatic"), "somatic_like",
         ifelse(final_class == "likely_germline", "germline_like",
                ifelse(final_class %in% c("likely_artifact", "technical_fail"), "artifact_like",
                       ifelse(final_class == "uncertain_tumor_only", "uncertain", "review_priority"))))
}

map_final_class_to_action <- function(final_class) {
  ifelse(final_class %in% c("high_confidence_somatic", "probable_somatic"), "review_as_somatic_like",
         ifelse(final_class == "likely_germline", "review_as_germline_like",
                ifelse(final_class %in% c("likely_artifact", "technical_fail"), "remove_as_artifact",
                       ifelse(final_class == "uncertain_tumor_only", "review_as_uncertain", "manual_review_required"))))
}

