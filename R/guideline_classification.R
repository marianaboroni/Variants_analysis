add_guideline_classifications <- function(x, cfg) {
  if (!isTRUE(cfg_get(cfg, c("guideline_classification", "enabled"), TRUE))) {
    x$somatic_oncogenicity_score <- NA_real_
    x$somatic_oncogenicity_class <- "not_evaluated"
    x$somatic_oncogenicity_criteria <- NA_character_
    x$germline_acmg_points <- NA_real_
    x$germline_acmg_class <- "not_evaluated"
    x$germline_acmg_criteria <- NA_character_
    x$germline_disposition <- "not_evaluated"
    return(x)
  }

  x <- add_somatic_oncogenicity_classification(x, cfg)
  x <- add_germline_acmg_classification(x, cfg)
  x
}

add_somatic_oncogenicity_classification <- function(x, cfg) {
  n <- nrow(x)
  points <- rep(0, n)
  criteria <- vector("list", n)

  add_points <- function(flag, value, code) {
    idx <- which(!is.na(flag) & flag)
    if (length(idx) == 0) return(invisible(NULL))
    points[idx] <<- points[idx] + value
    criteria[idx] <<- mapply(
      function(old, i) c(old, sprintf("%s:%+g", code, value)),
      criteria[idx],
      idx,
      SIMPLIFY = FALSE
    )
    invisible(NULL)
  }

  af <- x$max_pop_af
  add_points(!is.na(af) & af > cfg_get(cfg, c("guideline_classification", "somatic", "benign_maf_very_strong"), 0.05), -8, "SBVS1_population_maf_gt_5pct")
  add_points(!is.na(af) & af > cfg_get(cfg, c("guideline_classification", "somatic", "benign_maf_strong"), 0.01) & af <= 0.05, -4, "SBS1_population_maf_gt_1pct")
  add_points(!is.na(af) & af <= cfg_get(cfg, c("guideline_classification", "somatic", "absent_population_af"), 1e-5), 1, "OP4_absent_or_ultrarare_population")

  add_points(is_silent_no_splice(x), -1, "SBP1_silent_no_splice_prediction")
  add_points(has_benign_clinvar(x), -1, "SBP2_benign_reputable_database")
  add_points(has_pathogenic_clinvar(x), 1, "OP5_pathogenic_reputable_database")

  tsg_lof <- is_loss_of_function_variant(x) & is_tumor_suppressor_like(x)
  add_points(tsg_lof, 8, "OVS1_null_in_tumor_suppressor")

  add_points(safe_logical_guideline(x, "hotspot_match") & safe_logical_guideline(x, "hotspot_tumor_specific"), 4, "OS3_tumor_specific_hotspot")
  add_points(safe_logical_guideline(x, "hotspot_match") & !safe_logical_guideline(x, "hotspot_tumor_specific"), 2, "OM3_pan_cancer_hotspot")
  # OS4/OM4 use TUMOR-CONTEXT recurrence only (occurrences in the sample's tumor
  # type), computed in R/cosmic_context.R. Global-only recurrence (present only in
  # other tumor types, or non-evaluable context) MUST NOT grant OS4/OM4 and thus
  # cannot change final_class. See docs/COSMIC_OS4_OM4_AUDIT.md.
  add_points(safe_logical_guideline(x, "OS4_CONTEXTUAL"), 4, "OS4_high_somatic_recurrence_tumor_context")
  add_points(safe_logical_guideline(x, "OM4_CONTEXTUAL"), 2, "OM4_moderate_somatic_recurrence_tumor_context")
  # NOTE: the OncoKB OS5 oncogenicity point was intentionally removed. OncoKB is a
  # post-hoc annotation only (see R/oncokb.R) and MUST NOT influence filtering or
  # oncogenicity classification. See docs/REFACTOR_AUDIT.md (RISK-1).

  add_points(functional_prediction_supports_oncogenicity(x), 1, "OP1_computational_oncogenic_support")
  add_points(functional_prediction_supports_benignity(x), -1, "SBP3_computational_benign_support")
  add_points(safe_logical_guideline(x, "structural_hotspot"), 1, "OP2_functional_region")

  x$somatic_oncogenicity_score <- points
  x$somatic_oncogenicity_class <- classify_somatic_oncogenicity_points(points)
  x$somatic_oncogenicity_criteria <- vapply(criteria, collapse_criteria, character(1))
  x
}

classify_somatic_oncogenicity_points <- function(points) {
  ifelse(points >= 10, "Oncogenic",
         ifelse(points >= 6, "Likely Oncogenic",
                ifelse(points >= 0, "VUS",
                       ifelse(points >= -6, "Likely Benign", "Benign"))))
}

add_germline_acmg_classification <- function(x, cfg) {
  n <- nrow(x)
  points <- rep(0, n)
  criteria <- vector("list", n)

  add_points <- function(flag, value, code) {
    idx <- which(!is.na(flag) & flag)
    if (length(idx) == 0) return(invisible(NULL))
    points[idx] <<- points[idx] + value
    criteria[idx] <<- mapply(
      function(old, i) c(old, sprintf("%s:%+g", code, value)),
      criteria[idx],
      idx,
      SIMPLIFY = FALSE
    )
    invisible(NULL)
  }

  af <- x$max_pop_af
  ba1 <- !is.na(af) & af >= cfg_get(cfg, c("guideline_classification", "germline", "ba1_af"), 0.05)
  bs1 <- !is.na(af) & af >= cfg_get(cfg, c("guideline_classification", "germline", "bs1_af"), 0.01) & !ba1
  pm2 <- is.na(af) | af <= cfg_get(cfg, c("guideline_classification", "germline", "pm2_af"), 1e-5)

  add_points(ba1, -8, "BA1_common_population_af")
  add_points(bs1, -4, "BS1_af_too_high_for_disorder")
  add_points(pm2, 1, "PM2_supporting_absent_or_ultrarare")
  add_points(has_benign_clinvar(x), -1, "BP6_benign_reputable_database")
  add_points(has_pathogenic_clinvar(x), 1, "PP5_pathogenic_reputable_database")
  add_points(is_silent_no_splice(x), -1, "BP7_silent_no_splice_prediction")
  add_points(functional_prediction_supports_benignity(x), -1, "BP4_computational_benign_support")
  add_points(functional_prediction_supports_pathogenicity(x), 1, "PP3_computational_pathogenic_support")
  add_points(is_loss_of_function_variant(x) & is_tumor_suppressor_like(x), 8, "PVS1_null_in_lof_relevant_gene")
  add_points(safe_logical_guideline(x, "hotspot_match") | safe_logical_guideline(x, "structural_hotspot"), 2, "PM1_hotspot_or_functional_domain")
  add_points(is_missense_variant(x) & safe_logical_guideline(x, "hotspot_match"), 2, "PM5_missense_at_established_residue")

  x$germline_acmg_points <- points
  x$germline_acmg_class <- classify_acmg_points(points)
  x$germline_acmg_criteria <- vapply(criteria, collapse_criteria, character(1))
  x$germline_disposition <- classify_germline_disposition(x)
  x
}

classify_acmg_points <- function(points) {
  ifelse(points >= 10, "Pathogenic",
         ifelse(points >= 6, "Likely Pathogenic",
                ifelse(points <= -7, "Benign",
                       ifelse(points <= -1, "Likely Benign", "VUS"))))
}

classify_germline_disposition <- function(x) {
  germline_like_vaf <- (!is.na(x$distance_vaf_05) & x$distance_vaf_05 <= 0.12) |
    (!is.na(x$distance_vaf_10) & x$distance_vaf_10 <= 0.08)
  out <- rep("germline_uncertain", nrow(x))
  out[x$germline_acmg_class %in% c("Benign", "Likely Benign") & germline_like_vaf] <- "likely_benign_germline_polymorphism"
  out[x$germline_acmg_class %in% c("Benign", "Likely Benign") & !germline_like_vaf] <- "benign_germline_evidence_no_vaf_support"
  out[x$germline_acmg_class %in% c("Pathogenic", "Likely Pathogenic") & germline_like_vaf] <- "possible_pathogenic_germline_review"
  out[x$germline_acmg_class %in% c("Pathogenic", "Likely Pathogenic") & !germline_like_vaf] <- "pathogenic_evidence_without_germline_vaf"
  out
}

is_loss_of_function_variant <- function(x) {
  z <- tolower(as.character(x$consequence))
  grepl("frameshift|stop_gained|splice_acceptor|splice_donor|start_lost|nonsense", z)
}

is_missense_variant <- function(x) {
  grepl("missense", tolower(as.character(x$consequence)))
}

is_silent_no_splice <- function(x) {
  z <- tolower(as.character(x$consequence))
  silent <- grepl("synonymous|silent", z)
  splice <- grepl("splice", z)
  silent & !splice & zero_if_na_guideline(x$spliceai_max_score) < 0.20
}

is_tumor_suppressor_like <- function(x) {
  role <- tolower(as.character(if ("gene_driver_role" %in% names(x)) x$gene_driver_role else NA_character_))
  grepl("tumou?r suppressor|suppressor|loss|lof|inactiv", role)
}

functional_prediction_supports_oncogenicity <- function(x) {
  zero_if_na_guideline(x$predictor_functional_score) >= 0.70 |
    zero_if_na_guideline(x$spliceai_max_score) >= 0.50 |
    zero_if_na_guideline(x$alphamissense_score) >= 0.80 |
    zero_if_na_guideline(x$revel_score) >= 0.75 |
    zero_if_na_guideline(x$cadd_phred) >= 25
}

functional_prediction_supports_pathogenicity <- function(x) {
  functional_prediction_supports_oncogenicity(x)
}

functional_prediction_supports_benignity <- function(x) {
  has_predictor <- !is.na(x$alphamissense_score) |
    !is.na(x$revel_score) |
    !is.na(x$cadd_phred) |
    !is.na(x$spliceai_max_score) |
    !is.na(x$meta_predictor_score)
  has_predictor &
    zero_if_na_guideline(x$predictor_functional_score) <= 0.20 &
    zero_if_na_guideline(x$alphamissense_score) <= 0.34 &
    zero_if_na_guideline(x$revel_score) <= 0.35 &
    zero_if_na_guideline(x$cadd_phred) <= 15 &
    zero_if_na_guideline(x$spliceai_max_score) < 0.20
}

has_pathogenic_clinvar <- function(x) {
  cv <- tolower(as.character(if ("clinvar_significance" %in% names(x)) x$clinvar_significance else NA_character_))
  grepl("pathogenic", cv) & !grepl("benign", cv)
}

has_benign_clinvar <- function(x) {
  cv <- tolower(as.character(if ("clinvar_significance" %in% names(x)) x$clinvar_significance else NA_character_))
  grepl("benign", cv) & !grepl("pathogenic", cv)
}

safe_logical_guideline <- function(x, col) {
  if (!(col %in% names(x))) return(rep(FALSE, nrow(x)))
  z <- x[[col]]
  z[is.na(z)] <- FALSE
  as.logical(z)
}

zero_if_na_guideline <- function(x) {
  if (is.null(x)) return(0)
  x[is.na(x)] <- 0
  x
}

collapse_criteria <- function(z) {
  z <- unique(stats::na.omit(as.character(z)))
  if (length(z) == 0) "" else paste(z, collapse = ";")
}

guideline_summary <- function(x) {
  data.frame(
    metric = c(
      "somatic_oncogenicity_oncogenic",
      "somatic_oncogenicity_likely_oncogenic",
      "somatic_oncogenicity_vus",
      "somatic_oncogenicity_likely_benign_or_benign",
      "germline_acmg_pathogenic_or_likely_pathogenic",
      "germline_acmg_vus",
      "germline_acmg_likely_benign_or_benign",
      "possible_pathogenic_germline_review"
    ),
    value = c(
      sum(x$somatic_oncogenicity_class == "Oncogenic", na.rm = TRUE),
      sum(x$somatic_oncogenicity_class == "Likely Oncogenic", na.rm = TRUE),
      sum(x$somatic_oncogenicity_class == "VUS", na.rm = TRUE),
      sum(x$somatic_oncogenicity_class %in% c("Likely Benign", "Benign"), na.rm = TRUE),
      sum(x$germline_acmg_class %in% c("Pathogenic", "Likely Pathogenic"), na.rm = TRUE),
      sum(x$germline_acmg_class == "VUS", na.rm = TRUE),
      sum(x$germline_acmg_class %in% c("Likely Benign", "Benign"), na.rm = TRUE),
      sum(x$germline_disposition == "possible_pathogenic_germline_review", na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}
