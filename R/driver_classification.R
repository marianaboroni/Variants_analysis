add_driver_layers <- function(x, cfg) {
  x <- add_validation_keys(x)
  x <- attach_driver_gene_reference(x, cfg)
  x <- attach_hotspot_reference(x, cfg)
  x <- add_functional_predictor_features(x)
  x <- score_driver_evidence(x, cfg)
  x <- classify_driver_status(x, cfg)
  x
}

attach_driver_gene_reference <- function(x, cfg) {
  x$gene_driver_match <- FALSE
  x$gene_driver_tumor_specific <- FALSE
  x$gene_driver_source <- NA_character_
  x$gene_driver_role <- NA_character_
  x$gene_driver_confidence <- NA_character_

  path <- cfg_get(cfg, c("driver_resources", "driver_genes"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) return(x)

  ref <- read_variants(path, cfg_get(cfg, c("driver_resources", "driver_genes_delimiter"), "\t"))
  ref <- standardize_driver_gene_reference(ref)
  x$tumor_gene_key <- make_tumor_gene_key(x)
  x$gene_key <- make_gene_key(x)

  tumor_ref <- collapse_driver_gene_reference(ref, "tumor_gene_key")
  pan_ref <- collapse_driver_gene_reference(ref, "gene_key")

  idx_tumor <- match(x$tumor_gene_key, tumor_ref$tumor_gene_key)
  idx_pan <- match(x$gene_key, pan_ref$gene_key)
  hit_tumor <- !is.na(idx_tumor)
  hit_pan <- !hit_tumor & !is.na(idx_pan)

  x$gene_driver_match[hit_tumor | hit_pan] <- TRUE
  x$gene_driver_tumor_specific[hit_tumor] <- TRUE
  x$gene_driver_source[hit_tumor] <- tumor_ref$source[idx_tumor[hit_tumor]]
  x$gene_driver_role[hit_tumor] <- tumor_ref$role[idx_tumor[hit_tumor]]
  x$gene_driver_confidence[hit_tumor] <- tumor_ref$confidence[idx_tumor[hit_tumor]]
  x$gene_driver_source[hit_pan] <- pan_ref$source[idx_pan[hit_pan]]
  x$gene_driver_role[hit_pan] <- pan_ref$role[idx_pan[hit_pan]]
  x$gene_driver_confidence[hit_pan] <- pan_ref$confidence[idx_pan[hit_pan]]
  x
}

attach_hotspot_reference <- function(x, cfg) {
  x$hotspot_match <- FALSE
  x$hotspot_tumor_specific <- FALSE
  x$hotspot_source <- NA_character_
  x$hotspot_evidence <- NA_character_

  path <- cfg_get(cfg, c("driver_resources", "hotspots"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) return(x)

  ref <- read_variants(path, cfg_get(cfg, c("driver_resources", "hotspots_delimiter"), "\t"))
  ref <- standardize_hotspot_reference(ref)

  match <- prioritized_reference_lookup(
    x,
    ref,
    value_cols = c("source", "evidence"),
    prefixes = c("tumor_coord", "tumor_gene_protein", "coord", "gene_protein")
  )
  x$hotspot_match <- !is.na(match$source) | !is.na(match$evidence)
  x$hotspot_tumor_specific <- match$match_scope %in% c("tumor_coord", "tumor_gene_protein")
  x$hotspot_source <- match$source
  x$hotspot_evidence <- match$evidence
  x
}

standardize_driver_gene_reference <- function(ref) {
  out <- data.frame(
    gene = as.character(coalesce_columns(ref, c("gene", "Hugo_Symbol", "SYMBOL", "Gene"))),
    tumor_type = as.character(coalesce_columns(
      ref,
      c("tumor_type", "ONCOTREE_CODE", "Cancer_Type", "Tumor_Type", "Primary_Site"),
      default = "PANCANCER"
    )),
    role = as.character(coalesce_columns(
      ref,
      c("role", "mode_of_action", "Mode_of_Action", "MOA", "cancer_role"),
      default = "unknown"
    )),
    source = as.character(coalesce_columns(ref, c("source", "db", "database"), default = "driver_reference")),
    confidence = as.character(coalesce_columns(ref, c("confidence", "tier", "Tier", "level"), default = "curated")),
    stringsAsFactors = FALSE
  )
  out$gene <- toupper(trimws(out$gene))
  out$tumor_type[is_missing_value(out$tumor_type)] <- "PANCANCER"
  out$tumor_type <- toupper(trimws(out$tumor_type))
  out$role <- tolower(trimws(out$role))
  out$tumor_gene_key <- make_tumor_gene_key(out)
  out$gene_key <- make_gene_key(out)
  out[!is_missing_value(out$gene), ]
}

standardize_hotspot_reference <- function(ref) {
  out <- standardize_reference_core(ref)
  out$source <- as.character(coalesce_columns(ref, c("source", "db", "database"), default = "hotspot_reference"))
  out$evidence <- as.character(coalesce_columns(ref, c("evidence", "level", "confidence"), default = "curated_hotspot"))
  out
}

collapse_driver_gene_reference <- function(ref, key_col) {
  keep <- ref[!is_missing_value(ref[[key_col]]), c(key_col, "role", "source", "confidence"), drop = FALSE]
  aggregate(
    keep[c("role", "source", "confidence")],
    by = list(lookup_key = keep[[key_col]]),
    FUN = function(z) {
      z <- unique(na.omit(as.character(z)))
      if (length(z) == 0) NA_character_ else paste(z[1:min(length(z), 5)], collapse = ";")
    }
  ) |>
    stats::setNames(c(key_col, "role", "source", "confidence"))
}

make_gene_key <- function(x) {
  toupper(trimws(as.character(x$gene)))
}

make_tumor_gene_key <- function(x) {
  paste(toupper(trimws(as.character(x$tumor_type))), make_gene_key(x), sep = "|")
}

add_functional_predictor_features <- function(x) {
  x$spliceai_max_score <- max_numeric_columns(x, c(
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
    "SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL",
    "DS_AG", "DS_AL", "DS_DG", "DS_DL"
  ))
  x$alphamissense_score <- to_numeric_safe(coalesce_columns(x, c("AlphaMissense_score", "alphamissense_score")))
  x$revel_score <- to_numeric_safe(coalesce_columns(x, c("REVEL", "REVEL_score", "revel_score")))
  x$cadd_phred <- to_numeric_safe(coalesce_columns(x, c("CADD_PHRED", "CADD_phred", "CADD", "CADD_score")))
  x$meta_predictor_score <- max_numeric_columns(x, c(
    "MetaRNN_score", "MetaLR_score", "MutationTaster_score", "M-CAP_score", "MPC_score"
  ))

  x$is_truncating <- grepl(
    "frameshift|stop_gained|splice_acceptor|splice_donor|start_lost|stop_lost|nonsense",
    tolower(as.character(x$consequence))
  )
  x$is_splice_disruptive <- (!is.na(x$spliceai_max_score) & x$spliceai_max_score >= 0.5) |
    grepl("splice_acceptor|splice_donor", tolower(as.character(x$consequence)))
  x$is_missense <- grepl("missense", tolower(as.character(x$consequence)))
  x$is_inframe <- grepl("inframe", tolower(as.character(x$consequence)))
  x$structural_annotation <- as.character(coalesce_columns(
    x,
    c("structure_feature", "functional_region", "AlphaFold_feature", "COSMIC3D",
      "cosmic3d", "UniProt_domain", "DOMAINS", "Domain", "protein_domain"),
    default = NA
  ))
  x$domain_annotated <- !is_missing_value(x$structural_annotation)
  x$structural_hotspot <- grepl(
    "active|binding|interface|pocket|catalytic|hotspot|cosmic.?3d|functional",
    tolower(x$structural_annotation)
  )

  x$predictor_functional_score <- bounded01(rowMeans(cbind(
    predictor_threshold_score(x$alphamissense_score, low = 0.34, high = 0.80),
    predictor_threshold_score(x$revel_score, low = 0.50, high = 0.80),
    predictor_threshold_score(x$cadd_phred, low = 20, high = 30),
    predictor_threshold_score(x$spliceai_max_score, low = 0.20, high = 0.80),
    predictor_threshold_score(x$meta_predictor_score, low = 0.50, high = 0.80)
  ), na.rm = TRUE))
  x$predictor_functional_score[is.na(x$predictor_functional_score)] <- 0
  x$functional_impact_score <- bounded01(
    0.35 * zero_if_na(pmin(x$impact_rank / 3, 1)) +
      0.30 * zero_if_na(x$predictor_functional_score) +
      0.20 * ifelse(x$is_truncating | x$is_splice_disruptive, 1, 0) +
      0.10 * ifelse(x$domain_annotated, 1, 0) +
      0.10 * ifelse(x$structural_hotspot, 1, 0) +
      0.05 * ifelse(x$is_inframe, 1, 0)
  )
  x
}

max_numeric_columns <- function(x, cols) {
  cols <- cols[cols %in% names(x)]
  if (length(cols) == 0) return(rep(NA_real_, nrow(x)))
  mat <- do.call(cbind, lapply(cols, function(nm) to_numeric_safe(x[[nm]])))
  suppressWarnings(apply(mat, 1, function(z) {
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  }))
}

predictor_threshold_score <- function(x, low, high) {
  out <- (x - low) / (high - low)
  out[is.na(out)] <- NA_real_
  bounded01(out)
}

score_driver_evidence <- function(x, cfg) {
  x$oncokb_driver_score <- oncokb_driver_score(x)
  x$cosmic_driver_score <- cosmic_driver_score(x)
  x$hotspot_driver_score <- ifelse(x$hotspot_match, ifelse(x$hotspot_tumor_specific, 1, 0.85), 0)
  x$variant_driver_evidence_score <- pmax(
    x$oncokb_driver_score,
    x$cosmic_driver_score,
    x$hotspot_driver_score,
    na.rm = TRUE
  )
  x$variant_driver_evidence_score[is.infinite(x$variant_driver_evidence_score)] <- 0

  x$gene_driver_score <- gene_driver_score(x)
  x$mechanism_compatibility_score <- mechanism_compatibility_score(x)
  x$cohort_driver_signal_score <- cohort_driver_signal_score(x, cfg)
  x$driver_penalty <- driver_penalty_score(x)

  x$driver_score <- bounded01(
    0.40 * zero_if_na(x$variant_driver_evidence_score) +
      0.20 * zero_if_na(x$gene_driver_score) +
      0.20 * zero_if_na(x$functional_impact_score) +
      0.10 * zero_if_na(x$mechanism_compatibility_score) +
      0.10 * zero_if_na(x$cohort_driver_signal_score) -
      zero_if_na(x$driver_penalty)
  )
  x$driver_evidence <- driver_evidence_text(x)
  x
}

zero_if_na <- function(x) {
  x[is.na(x)] <- 0
  x
}

oncokb_driver_score <- function(x) {
  onc <- toupper(as.character(x$oncokb_oncogenic))
  score <- rep(0, nrow(x))
  score[grepl("ONCOGENIC", onc)] <- 1
  score[grepl("LIKELY ONCOGENIC", onc)] <- 0.9
  score[grepl("RESISTANCE", onc)] <- pmax(score[grepl("RESISTANCE", onc)], 0.85)
  score[x$oncokb_match & score == 0] <- 0.45
  score
}

cosmic_driver_score <- function(x) {
  count <- x$cosmic_count
  score <- rep(0, nrow(x))
  score[x$cosmic_match] <- 0.35
  score[!is.na(count) & count >= 3] <- 0.55
  score[!is.na(count) & count >= 10] <- 0.75
  score[!is.na(count) & count >= 50] <- 0.90
  score
}

gene_driver_score <- function(x) {
  score <- rep(0, nrow(x))
  score[x$gene_driver_match] <- 0.50
  score[x$gene_driver_tumor_specific] <- 0.75
  conf <- tolower(as.character(x$gene_driver_confidence))
  score[grepl("tier.?1|canonical|high|curated", conf) & x$gene_driver_match] <-
    pmax(score[grepl("tier.?1|canonical|high|curated", conf) & x$gene_driver_match], 0.85)
  score
}

mechanism_compatibility_score <- function(x) {
  role <- tolower(as.character(x$gene_driver_role))
  oncogene <- grepl("oncogene|activa|gain|gof", role)
  tsg <- grepl("tumou?r suppressor|suppressor|loss|lof|inactiv", role)
  both <- grepl("both|dual", role)
  hotspot_like <- x$hotspot_match | x$oncokb_driver_score >= 0.85 |
    (!is.na(x$cosmic_count) & x$cosmic_count >= 10)

  score <- rep(0, nrow(x))
  score[oncogene & (x$is_missense | x$is_inframe) & hotspot_like] <- 1
  score[oncogene & (x$is_missense | x$is_inframe) & !hotspot_like] <- 0.55
  score[tsg & (x$is_truncating | x$is_splice_disruptive)] <- 1
  score[tsg & x$is_missense & x$functional_impact_score >= 0.65] <- 0.65
  score[both & score == 0] <- pmax(score[both & score == 0], 0.60)
  score[x$gene_driver_match & score == 0 & x$functional_impact_score >= 0.70] <- 0.45
  score
}

cohort_driver_signal_score <- function(x, cfg) {
  min_samples <- cfg_get(cfg, c("driver_scoring", "min_tumor_type_samples_for_recurrence_driver"), 5)
  tumor_n <- if ("tumor_type_total_samples" %in% names(x)) x$tumor_type_total_samples else rep(NA_real_, nrow(x))
  tumor_freq <- if ("variant_tumor_type_freq" %in% names(x)) x$variant_tumor_type_freq else rep(NA_real_, nrow(x))
  cohort_freq <- if ("variant_cohort_freq" %in% names(x)) x$variant_cohort_freq else rep(NA_real_, nrow(x))
  score <- rep(0, nrow(x))
  enough <- !is.na(tumor_n) & tumor_n >= min_samples
  score[enough & !is.na(tumor_freq) & tumor_freq >= 0.05 & x$gene_driver_match] <- 0.45
  score[enough & !is.na(tumor_freq) & tumor_freq >= 0.10 & (x$hotspot_match | x$oncokb_driver_score >= 0.85)] <- 0.75
  score[!is.na(cohort_freq) & cohort_freq >= 0.20 & !x$gene_driver_match & !x$hotspot_match] <- 0
  score
}

driver_penalty_score <- function(x) {
  penalty <- rep(0, nrow(x))
  penalty[x$final_class %in% c("hard_filter_fail", "probable_artifact", "probable_germline")] <- 0.70
  penalty[x$final_class %in% c("uncertain")] <- 0.20
  penalty[!is.na(x$max_pop_af) & x$max_pop_af >= 0.005] <- pmax(
    penalty[!is.na(x$max_pop_af) & x$max_pop_af >= 0.005],
    0.50
  )
  penalty
}

classify_driver_status <- function(x, cfg) {
  known_cutoff <- cfg_get(cfg, c("driver_scoring", "known_driver"), 0.75)
  probable_cutoff <- cfg_get(cfg, c("driver_scoring", "probable_driver"), 0.60)
  possible_cutoff <- cfg_get(cfg, c("driver_scoring", "possible_driver"), 0.40)
  somatic_ok <- x$final_class %in% c("high_confidence_somatic", "probable_somatic")

  x$driver_class <- "not_evaluable_as_driver"
  x$driver_class[somatic_ok & x$driver_score < possible_cutoff] <- "likely_passenger"
  x$driver_class[somatic_ok & x$driver_score >= possible_cutoff] <- "possible_driver"
  x$driver_class[somatic_ok & x$driver_score >= probable_cutoff &
                   ((x$gene_driver_score > 0 & x$functional_impact_score >= 0.35) |
                      x$variant_driver_evidence_score >= 0.85)] <- "probable_driver"
  x$driver_class[somatic_ok &
                   x$driver_score >= probable_cutoff &
                   x$variant_driver_evidence_score >= 0.85 &
                   (x$hotspot_match | x$oncokb_driver_score >= 0.85)] <- "known_driver"
  x$driver_class[somatic_ok &
                   x$driver_score >= known_cutoff &
                   x$variant_driver_evidence_score >= 0.85] <- "known_driver"
  x$driver_class[x$final_class == "uncertain" & x$driver_score >= probable_cutoff] <- "uncertain_possible_driver"
  x
}

driver_evidence_text <- function(x) {
  out <- rep("", nrow(x))
  add <- function(flag, label) {
    idx <- !is.na(flag) & flag
    out[idx] <<- ifelse(out[idx] == "", label, paste(out[idx], label, sep = ";"))
  }
  add(x$oncokb_driver_score >= 0.85, "OncoKB_oncogenic")
  add(x$cosmic_driver_score >= 0.75, "COSMIC_recurrent")
  add(x$hotspot_match, "hotspot_reference")
  add(x$gene_driver_tumor_specific, "tumor_type_driver_gene")
  add(x$gene_driver_match & !x$gene_driver_tumor_specific, "pan_cancer_driver_gene")
  add(x$mechanism_compatibility_score >= 0.75, "mechanism_compatible")
  add(x$functional_impact_score >= 0.70, "strong_functional_prediction")
  add(x$is_splice_disruptive, "splice_disruptive")
  add(x$structural_hotspot, "structural_functional_region")
  add(!is.na(x$max_pop_af) & x$max_pop_af >= 0.005, "population_AF_penalty")
  out[out == ""] <- "no_driver_evidence"
  out
}
