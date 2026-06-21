run_ml_filter <- function(x, cfg) {
  if (!isTRUE(cfg_get(cfg, c("ml_filter", "enabled"), TRUE))) {
    x <- initialize_ml_columns(x)
    return(list(
      variants = x,
      metrics = ml_metric_table("global", "PANCANCER", c(status = "disabled")),
      training_labels = empty_training_label_audit(),
      active_learning_candidates = empty_active_learning_candidates()
    ))
  }

  x <- attach_review_labels(x, cfg)
  x$ml_pseudo_label <- make_pseudo_labels(x, cfg)
  x$ml_manual_training_label <- map_reviewed_label_to_ml(x$reviewed_label)
  x$ml_label <- ifelse(!is.na(x$ml_manual_training_label), x$ml_manual_training_label, x$ml_pseudo_label)
  x$ml_label_source <- ifelse(
    !is.na(x$ml_manual_training_label),
    "manual_review",
    ifelse(!is.na(x$ml_pseudo_label), "pseudo_rules", NA_character_)
  )
  x$ml_label_weight <- ml_label_weights(x, cfg)

  feature_matrix <- build_ml_feature_matrix(x)
  model_result <- train_incremental_ml_models(x, feature_matrix, cfg)
  x <- model_result$variants

  training_labels <- build_training_label_audit(x)
  active_learning_candidates <- build_active_learning_candidates(x, cfg)

  list(
    variants = x,
    metrics = model_result$metrics,
    training_labels = training_labels,
    active_learning_candidates = active_learning_candidates
  )
}

initialize_ml_columns <- function(x) {
  x$reviewed_label <- NA_character_
  x$review_match_scope <- NA_character_
  x$review_evidence <- NA_character_
  x$reviewer <- NA_character_
  x$review_date <- NA_character_
  x$review_notes <- NA_character_
  x$ml_pseudo_label <- NA_character_
  x$ml_manual_training_label <- NA_character_
  x$ml_label <- NA_character_
  x$ml_label_source <- NA_character_
  x$ml_label_weight <- NA_real_
  x$ml_true_positive_probability <- NA_real_
  x$ml_somatic_probability <- NA_real_
  x$ml_germline_probability <- NA_real_
  x$ml_artifact_probability <- NA_real_
  x$ml_uncertainty <- NA_real_
  x$ml_blocked_reason <- NA_character_
  x$ml_pancancer_true_positive_probability <- NA_real_
  x$ml_tumor_type_true_positive_probability <- NA_real_
  x$ml_predicted_class <- NA_character_
  x$ml_model_scope <- NA_character_
  x$ml_model_id <- NA_character_
  x$ml_model_version <- NA_character_
  x
}

attach_review_labels <- function(x, cfg) {
  x$reviewed_label <- NA_character_
  x$review_match_scope <- NA_character_
  x$review_evidence <- NA_character_
  x$reviewer <- NA_character_
  x$review_date <- NA_character_
  x$review_notes <- NA_character_

  review <- load_review_labels(cfg)
  if (nrow(review) == 0) return(x)

  x <- add_validation_keys(x)
  prefixes <- cfg_get(
    cfg,
    c("ml_filter", "review_match_priority"),
    c("sample_coord", "sample_gene_protein", "tumor_coord", "tumor_gene_protein", "coord", "gene_protein")
  )
  match <- prioritized_reference_lookup(
    x,
    review,
    value_cols = c("reviewed_label_src", "review_evidence_src", "reviewer_src", "review_date_src", "review_notes_src"),
    prefixes = prefixes
  )

  x$reviewed_label <- match$reviewed_label_src
  x$review_match_scope <- match$match_scope
  x$review_evidence <- match$review_evidence_src
  x$reviewer <- match$reviewer_src
  x$review_date <- match$review_date_src
  x$review_notes <- match$review_notes_src
  x
}

load_review_labels <- function(cfg) {
  path <- cfg_get(cfg, c("ml_filter", "review_labels"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    return(data.frame())
  }
  delimiter <- cfg_get(cfg, c("ml_filter", "review_labels_delimiter"), "\t")
  review <- read_variants(path, delimiter)
  standardize_review_label_table(review)
}

standardize_review_label_table <- function(review) {
  review <- standardize_validation_table(review)
  review$reviewed_label_src <- as.character(coalesce_columns(
    review,
    c("reviewed_label", "manual_label", "curated_label", "truth_label", "label", "classification")
  ))
  review$review_evidence_src <- as.character(coalesce_columns(
    review,
    c("evidence", "review_evidence", "validation_evidence", "source"),
    default = NA
  ))
  review$reviewer_src <- as.character(coalesce_columns(
    review,
    c("reviewer", "curator", "analyst"),
    default = NA
  ))
  review$review_date_src <- as.character(coalesce_columns(
    review,
    c("review_date", "date", "curation_date"),
    default = NA
  ))
  review$review_notes_src <- as.character(coalesce_columns(
    review,
    c("notes", "comment", "comments", "review_notes"),
    default = NA
  ))
  review[!is_missing_value(review$reviewed_label_src), ]
}

map_reviewed_label_to_ml <- function(label) {
  z <- tolower(trimws(as.character(label)))
  z[is.na(label) | z == "" | z == "na"] <- NA_character_

  positive <- grepl(
    "true_somatic|validated_somatic|somatic_true|true_positive|\\btp\\b|real_somatic|confirmed_somatic",
    z
  )
  negative <- grepl(
    "false_positive|\\bfp\\b|artifact|artefact|technical|germline|polymorphism|polymorphism|snp|not_somatic|contamination",
    z
  )
  uncertain <- grepl("uncertain|unknown|review|ambiguous|exclude|low_confidence|not_evaluable", z)

  out <- rep(NA_character_, length(z))
  out[positive & !negative & !uncertain] <- "true_positive"
  out[negative & !positive & !uncertain] <- "false_positive"
  out
}

ml_label_weights <- function(x, cfg) {
  manual_weight <- cfg_get(cfg, c("ml_filter", "manual_label_weight"), 5)
  pseudo_weight <- cfg_get(cfg, c("ml_filter", "pseudo_label_weight"), 1)
  weights <- rep(NA_real_, nrow(x))
  weights[x$ml_label_source == "manual_review"] <- manual_weight
  weights[x$ml_label_source == "pseudo_rules"] <- pseudo_weight
  weights
}

make_pseudo_labels <- function(x, cfg) {
  label <- rep(NA_character_, nrow(x))
  strong_validated <- zero_if_na(x$validation_support_score) >= 0.8 |
    x$driver_class %in% c("known_driver", "probable_driver")
  high_somatic <- x$final_class == "high_confidence_somatic" &
    x$hard_filter_pass &
    x$germline_score < 0.45 &
    x$artifact_score < 0.45
  probable_fp <- x$final_class %in% c("hard_filter_fail", "probable_artifact", "probable_germline") &
    (x$artifact_score >= 0.70 | x$germline_score >= 0.70 | !is.na(x$max_pop_af) & x$max_pop_af >= 0.005)

  label[high_somatic & (strong_validated | x$somatic_score_validated >= 0.85)] <- "true_positive"
  label[probable_fp] <- "false_positive"
  label
}

train_incremental_ml_models <- function(x, feature_matrix, cfg) {
  x$ml_true_positive_probability <- NA_real_
  x$ml_pancancer_true_positive_probability <- NA_real_
  x$ml_tumor_type_true_positive_probability <- NA_real_
  x$ml_germline_probability <- NA_real_
  x$ml_artifact_probability <- NA_real_
  x$ml_uncertainty <- NA_real_
  x$ml_blocked_reason <- NA_character_
  x$ml_predicted_class <- NA_character_
  x$ml_model_scope <- NA_character_
  x$ml_model_id <- NA_character_
  x$ml_model_version <- NA_character_

  cutoff <- cfg_get(cfg, c("ml_filter", "true_positive_cutoff"), 0.60)
  seed <- cfg_get(cfg, c("ml_filter", "seed"), 20260621)
  set.seed(seed)

  metrics <- empty_ml_metrics()

  pan_idx <- which(!is.na(x$ml_label))
  pan_model <- fit_ml_scope_model(x, feature_matrix, pan_idx, cfg, scope = "pancancer", tumor_type = "PANCANCER")
  metrics <- rbind(metrics, pan_model$metrics)

  if (isTRUE(pan_model$trained)) {
    pan_prob <- predict_ml_model(pan_model$model, feature_matrix)
    x$ml_pancancer_true_positive_probability <- pan_prob
    x$ml_true_positive_probability <- pan_prob
    x$ml_model_scope[!is.na(pan_prob)] <- "pancancer"
    x$ml_model_id[!is.na(pan_prob)] <- "glm_PANCANCER"
    x$ml_model_version[!is.na(pan_prob)] <- model_version_string(pan_model$model)
  }

  use_tumor_models <- isTRUE(cfg_get(cfg, c("ml_filter", "tumor_type_models"), TRUE))
  if (use_tumor_models) {
    tumor_types <- sort(unique(na.omit(as.character(x$tumor_type))))
    tumor_types <- tumor_types[!is_missing_value(tumor_types)]
    for (tt in tumor_types) {
      row_idx <- which(x$tumor_type == tt)
      train_idx <- row_idx[!is.na(x$ml_label[row_idx])]
      tumor_model <- fit_ml_scope_model(
        x,
        feature_matrix,
        train_idx,
        cfg,
        scope = "tumor_type",
        tumor_type = tt,
        per_tumor = TRUE
      )
      metrics <- rbind(metrics, tumor_model$metrics)
      if (!isTRUE(tumor_model$trained)) next

      tumor_prob <- predict_ml_model(tumor_model$model, feature_matrix[row_idx, , drop = FALSE])
      x$ml_tumor_type_true_positive_probability[row_idx] <- tumor_prob
      x$ml_true_positive_probability[row_idx] <- tumor_prob
      x$ml_model_scope[row_idx] <- "tumor_type"
      x$ml_model_id[row_idx] <- paste0("glm_", tt)
      x$ml_model_version[row_idx] <- model_version_string(tumor_model$model)
    }
  }

  prob <- x$ml_true_positive_probability
  x$ml_somatic_probability <- prob
  x$ml_predicted_class <- ifelse(prob >= cutoff, "ml_likely_true_positive", "ml_likely_false_positive")
  x$ml_predicted_class[is.na(prob)] <- NA_character_
  x$ml_germline_probability <- 1 - pmin(pmax(prob, 0), 1)
  x$ml_artifact_probability <- x$ml_germline_probability
  x$ml_uncertainty <- ifelse(is.na(prob), NA_real_, 1 - abs(prob - 0.5) * 2)

  hard_pass <- if ("hard_filter_pass" %in% names(x)) x$hard_filter_pass else rep(TRUE, nrow(x))
  blocked <- !isTRUE(cfg_get(cfg, c("ml_filter", "allow_ml_override_hard_filters"), FALSE)) & !hard_pass
  x$ml_blocked_reason <- NA_character_
  x$ml_blocked_reason[blocked] <- "hard_filter_fail"

  if (!any(metrics$metric == "status" & metrics$value == "trained")) {
    x$ml_model_scope[] <- NA_character_
    x$ml_model_id[] <- NA_character_
    x$ml_model_version[] <- NA_character_
    x$ml_predicted_class[] <- NA_character_
    x$ml_true_positive_probability[] <- NA_real_
    x$ml_somatic_probability[] <- NA_real_
    x$ml_germline_probability[] <- NA_real_
    x$ml_artifact_probability[] <- NA_real_
    x$ml_uncertainty[] <- NA_real_
    x$ml_pancancer_true_positive_probability[] <- NA_real_
    x$ml_tumor_type_true_positive_probability[] <- NA_real_
  }

  list(variants = x, metrics = metrics)
}

fit_ml_scope_model <- function(x, feature_matrix, candidate_idx, cfg, scope, tumor_type, per_tumor = FALSE) {
  min_pos <- cfg_get(
    cfg,
    if (per_tumor) c("ml_filter", "per_tumor_min_positive_labels") else c("ml_filter", "min_positive_labels"),
    cfg_get(cfg, c("ml_filter", "min_positive_labels"), 50)
  )
  min_neg <- cfg_get(
    cfg,
    if (per_tumor) c("ml_filter", "per_tumor_min_negative_labels") else c("ml_filter", "min_negative_labels"),
    cfg_get(cfg, c("ml_filter", "min_negative_labels"), 50)
  )
  max_rows <- cfg_get(cfg, c("ml_filter", "max_training_rows"), 200000)

  label_counts <- ml_label_count_summary(x, candidate_idx)
  status_values <- c(
    status = "insufficient_training_labels",
    model = "glm_binomial_incremental",
    n_positive_labels = label_counts$n_pos,
    n_negative_labels = label_counts$n_neg,
    n_manual_positive_labels = label_counts$n_manual_pos,
    n_manual_negative_labels = label_counts$n_manual_neg,
    n_pseudo_positive_labels = label_counts$n_pseudo_pos,
    n_pseudo_negative_labels = label_counts$n_pseudo_neg,
    n_training_rows = 0,
    training_auc = NA,
    true_positive_cutoff = cfg_get(cfg, c("ml_filter", "true_positive_cutoff"), 0.60)
  )

  if (label_counts$n_pos < min_pos || label_counts$n_neg < min_neg) {
    return(list(
      trained = FALSE,
      model = NULL,
      metrics = ml_metric_table(scope, tumor_type, status_values)
    ))
  }

  train_idx <- stratified_sample_training_idx(x$ml_label, candidate_idx, max_rows)
  y <- ifelse(x$ml_label[train_idx] == "true_positive", 1, 0)
  weights <- x$ml_label_weight[train_idx]
  weights[is.na(weights) | weights <= 0] <- 1
  train <- as.data.frame(feature_matrix[train_idx, , drop = FALSE])
  train$y <- y

  model <- tryCatch(
    suppressWarnings(stats::glm(y ~ ., data = train, family = stats::binomial(), weights = weights)),
    error = function(e) structure(list(message = conditionMessage(e)), class = "ml_model_error")
  )
  if (inherits(model, "ml_model_error")) {
    status_values["status"] <- paste0("training_failed: ", model$message)
    return(list(
      trained = FALSE,
      model = NULL,
      metrics = ml_metric_table(scope, tumor_type, status_values)
    ))
  }

  train_prob <- predict_ml_model(model, feature_matrix[train_idx, , drop = FALSE])
  status_values["status"] <- "trained"
  status_values["n_training_rows"] <- length(train_idx)
  status_values["training_auc"] <- round(binary_auc(y, train_prob), 4)
  list(
    trained = TRUE,
    model = model,
    metrics = ml_metric_table(scope, tumor_type, status_values)
  )
}

model_version_string <- function(model) {
  if (is.null(model) || !inherits(model, "glm")) return(NA_character_)
  version <- tryCatch(
    sprintf("%s_%s", model$family$family, format(Sys.Date(), "%Y%m%d")),
    error = function(e) NA_character_
  )
  version
}

ml_label_count_summary <- function(x, idx) {
  labels <- x$ml_label[idx]
  sources <- x$ml_label_source[idx]
  list(
    n_pos = sum(labels == "true_positive", na.rm = TRUE),
    n_neg = sum(labels == "false_positive", na.rm = TRUE),
    n_manual_pos = sum(labels == "true_positive" & sources == "manual_review", na.rm = TRUE),
    n_manual_neg = sum(labels == "false_positive" & sources == "manual_review", na.rm = TRUE),
    n_pseudo_pos = sum(labels == "true_positive" & sources == "pseudo_rules", na.rm = TRUE),
    n_pseudo_neg = sum(labels == "false_positive" & sources == "pseudo_rules", na.rm = TRUE)
  )
}

empty_ml_metrics <- function() {
  data.frame(
    model_scope = character(),
    tumor_type = character(),
    metric = character(),
    value = character(),
    stringsAsFactors = FALSE
  )
}

ml_metric_table <- function(scope, tumor_type, values) {
  data.frame(
    model_scope = scope,
    tumor_type = tumor_type,
    metric = names(values),
    value = as.character(values),
    stringsAsFactors = FALSE
  )
}

predict_ml_model <- function(model, feature_matrix) {
  out <- tryCatch(
    suppressWarnings(as.numeric(stats::predict(model, newdata = as.data.frame(feature_matrix), type = "response"))),
    error = function(e) rep(NA_real_, nrow(feature_matrix))
  )
  out[!is.finite(out)] <- NA_real_
  out
}

build_ml_feature_matrix <- function(x) {
  feature_cols <- c(
    "dp", "alt_count", "vaf", "tlod", "mbq", "mmq", "strandq",
    "popaf_mutect", "contq", "max_pop_af", "variant_cohort_freq",
    "locus_cohort_freq", "variant_tumor_type_freq", "locus_tumor_type_freq",
    "technical_evidence_score", "population_germline_score", "vaf_germline_score",
    "germline_score", "artifact_score", "somatic_score_validated",
    "functional_impact_score", "driver_score", "spliceai_max_score",
    "alphamissense_score", "revel_score", "cadd_phred", "meta_predictor_score",
    "ccf_estimate", "clonality_cluster_center"
  )
  feature_cols <- feature_cols[feature_cols %in% names(x)]
  if (length(feature_cols) == 0) {
    mat <- matrix(0, nrow = nrow(x), ncol = 1)
    colnames(mat) <- "constant_feature"
    return(mat)
  }
  mat <- do.call(cbind, lapply(feature_cols, function(col) to_numeric_safe(x[[col]])))
  colnames(mat) <- feature_cols
  mat <- impute_numeric_matrix(mat)
  mat
}

impute_numeric_matrix <- function(mat) {
  if (ncol(mat) == 0) return(mat)
  for (j in seq_len(ncol(mat))) {
    z <- mat[, j]
    if (all(is.na(z))) {
      mat[, j] <- 0
    } else {
      med <- stats::median(z, na.rm = TRUE)
      z[is.na(z)] <- med
      mat[, j] <- z
    }
  }
  mat
}

stratified_sample_training_idx <- function(labels, candidate_idx, max_rows) {
  idx_pos <- candidate_idx[labels[candidate_idx] == "true_positive"]
  idx_neg <- candidate_idx[labels[candidate_idx] == "false_positive"]
  per_class <- floor(max_rows / 2)
  idx_pos <- if (length(idx_pos) > per_class) sample(idx_pos, per_class) else idx_pos
  idx_neg <- if (length(idx_neg) > per_class) sample(idx_neg, per_class) else idx_neg
  sample(c(idx_pos, idx_neg))
}

binary_auc <- function(y, prob) {
  pos <- prob[y == 1]
  neg <- prob[y == 0]
  if (length(pos) == 0 || length(neg) == 0 || all(is.na(prob))) return(NA_real_)
  ranks <- rank(c(pos, neg), na.last = "keep")
  n_pos <- length(pos)
  n_neg <- length(neg)
  (sum(ranks[seq_len(n_pos)], na.rm = TRUE) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

build_training_label_audit <- function(x) {
  keep <- !is.na(x$ml_label)
  if (!any(keep)) return(empty_training_label_audit())
  cols <- c(
    "sample_id", "tumor_type", "chrom", "pos", "ref", "alt", "gene", "protein_change",
    "final_class", "driver_class", "reviewed_label", "review_match_scope", "review_evidence",
    "reviewer", "review_date", "ml_pseudo_label", "ml_manual_training_label",
    "ml_label", "ml_label_source", "ml_label_weight", "somatic_score", "germline_score",
    "artifact_score", "max_pop_af", "variant_cohort_freq", "variant_tumor_type_freq"
  )
  cols <- cols[cols %in% names(x)]
  x[keep, cols, drop = FALSE]
}

empty_training_label_audit <- function() {
  data.frame(
    sample_id = character(),
    tumor_type = character(),
    chrom = character(),
    pos = numeric(),
    ref = character(),
    alt = character(),
    gene = character(),
    protein_change = character(),
    final_class = character(),
    driver_class = character(),
    reviewed_label = character(),
    review_match_scope = character(),
    review_evidence = character(),
    reviewer = character(),
    review_date = character(),
    ml_pseudo_label = character(),
    ml_manual_training_label = character(),
    ml_label = character(),
    ml_label_source = character(),
    ml_label_weight = numeric(),
    somatic_score = numeric(),
    germline_score = numeric(),
    artifact_score = numeric(),
    max_pop_af = numeric(),
    variant_cohort_freq = numeric(),
    variant_tumor_type_freq = numeric(),
    stringsAsFactors = FALSE
  )
}

build_active_learning_candidates <- function(x, cfg) {
  if (!isTRUE(cfg_get(cfg, c("active_learning", "enabled"), TRUE))) {
    return(empty_active_learning_candidates())
  }

  max_candidates <- cfg_get(cfg, c("active_learning", "max_candidates"), 500)
  cutoff <- cfg_get(cfg, c("ml_filter", "true_positive_cutoff"), 0.60)
  include_reviewed <- isTRUE(cfg_get(cfg, c("active_learning", "include_reviewed"), FALSE))

  prob <- x$ml_true_positive_probability
  uncertainty <- rep(0, nrow(x))
  has_prob <- !is.na(prob)
  uncertainty[has_prob] <- 1 - pmin(abs(prob[has_prob] - cutoff) / max(cutoff, 1 - cutoff), 1)

  somatic_by_rules <- x$final_class %in% c("high_confidence_somatic", "probable_somatic")
  nonsomatic_by_rules <- x$final_class %in% c("probable_germline", "probable_artifact", "hard_filter_fail")
  ml_somatic <- has_prob & prob >= cutoff
  discordance <- (somatic_by_rules & has_prob & !ml_somatic) | (nonsomatic_by_rules & ml_somatic)

  validation_rescue <- safe_logical_column(x, "validation_rescue_candidate")
  driver_interest <- safe_character_column(x, "driver_class") %in% c("known_driver", "probable_driver", "possible_driver", "uncertain_possible_driver")
  min_recurrence_samples <- cfg_get(
    cfg,
    c("active_learning", "min_samples_for_recurrence"),
    cfg_get(cfg, c("cohort", "min_samples_for_recurrence_model"), 20)
  )
  enough_recurrence_samples <- !is.na(x$total_samples) & x$total_samples >= min_recurrence_samples
  recurrent_unvalidated <- enough_recurrence_samples &
    !is.na(x$variant_cohort_freq) & x$variant_cohort_freq >= cfg_get(cfg, c("active_learning", "recurrent_variant_fraction"), 0.05) &
    !safe_logical_column(x, "oncokb_match") & !safe_logical_column(x, "cosmic_match")
  population_conflict <- somatic_by_rules & !is.na(x$max_pop_af) & x$max_pop_af >= cfg_get(cfg, c("population_filters", "rare_af"), 0.001)
  high_impact_uncertain <- x$final_class == "uncertain" & zero_if_na(x$impact_rank) >= 2

  priority <- 0.35 * uncertainty +
    0.25 * as.numeric(discordance) +
    0.15 * as.numeric(x$final_class == "uncertain") +
    0.10 * as.numeric(validation_rescue) +
    0.10 * as.numeric(driver_interest) +
    0.10 * as.numeric(recurrent_unvalidated) +
    0.05 * as.numeric(population_conflict) +
    0.05 * as.numeric(high_impact_uncertain)

  reviewed <- !is.na(x$ml_manual_training_label)
  candidate <- priority > 0
  if (!include_reviewed) candidate <- candidate & !reviewed

  if (!any(candidate)) return(empty_active_learning_candidates())

  reasons <- active_learning_reasons(
    uncertainty >= cfg_get(cfg, c("active_learning", "uncertainty_score_min"), 0.70),
    discordance,
    x$final_class == "uncertain",
    validation_rescue,
    driver_interest,
    recurrent_unvalidated,
    population_conflict,
    high_impact_uncertain
  )

  out <- data.frame(
    review_rank = NA_integer_,
    active_learning_priority = priority,
    active_learning_reason = reasons,
    sample_id = x$sample_id,
    tumor_type = x$tumor_type,
    chrom = x$chrom,
    pos = x$pos,
    ref = x$ref,
    alt = x$alt,
    gene = x$gene,
    protein_change = x$protein_change,
    final_class = x$final_class,
    driver_class = safe_character_column(x, "driver_class"),
    ml_true_positive_probability = x$ml_true_positive_probability,
    ml_model_scope = x$ml_model_scope,
    ml_model_id = x$ml_model_id,
    somatic_score = x$somatic_score,
    germline_score = x$germline_score,
    artifact_score = x$artifact_score,
    max_pop_af = x$max_pop_af,
    variant_cohort_freq = x$variant_cohort_freq,
    variant_tumor_type_freq = x$variant_tumor_type_freq,
    dp = x$dp,
    alt_count = x$alt_count,
    vaf = x$vaf,
    stringsAsFactors = FALSE
  )
  out <- out[candidate, , drop = FALSE]
  out <- out[order(-out$active_learning_priority, out$sample_id, out$chrom, out$pos), , drop = FALSE]
  out <- head(out, max_candidates)
  out$review_rank <- seq_len(nrow(out))
  out
}

empty_active_learning_candidates <- function() {
  data.frame(
    review_rank = integer(),
    active_learning_priority = numeric(),
    active_learning_reason = character(),
    sample_id = character(),
    tumor_type = character(),
    chrom = character(),
    pos = numeric(),
    ref = character(),
    alt = character(),
    gene = character(),
    protein_change = character(),
    final_class = character(),
    driver_class = character(),
    ml_true_positive_probability = numeric(),
    ml_model_scope = character(),
    ml_model_id = character(),
    somatic_score = numeric(),
    germline_score = numeric(),
    artifact_score = numeric(),
    max_pop_af = numeric(),
    variant_cohort_freq = numeric(),
    variant_tumor_type_freq = numeric(),
    dp = numeric(),
    alt_count = numeric(),
    vaf = numeric(),
    stringsAsFactors = FALSE
  )
}

safe_logical_column <- function(x, col) {
  if (!(col %in% names(x))) return(rep(FALSE, nrow(x)))
  z <- x[[col]]
  z[is.na(z)] <- FALSE
  as.logical(z)
}

safe_character_column <- function(x, col) {
  if (!(col %in% names(x))) return(rep(NA_character_, nrow(x)))
  as.character(x[[col]])
}

active_learning_reasons <- function(
    uncertain_probability,
    discordance,
    uncertain_class,
    validation_rescue,
    driver_interest,
    recurrent_unvalidated,
    population_conflict,
    high_impact_uncertain) {
  n <- length(discordance)
  out <- rep("", n)
  add_reason <- function(flag, label) {
    idx <- !is.na(flag) & flag
    out[idx] <<- ifelse(out[idx] == "", label, paste(out[idx], label, sep = ";"))
  }
  add_reason(uncertain_probability, "ml_uncertain_probability")
  add_reason(discordance, "rule_ml_discordance")
  add_reason(uncertain_class, "uncertain_final_class")
  add_reason(validation_rescue, "validation_rescue_candidate")
  add_reason(driver_interest, "driver_or_hotspot_interest")
  add_reason(recurrent_unvalidated, "recurrent_unvalidated_variant")
  add_reason(population_conflict, "somatic_population_af_conflict")
  add_reason(high_impact_uncertain, "high_impact_uncertain")
  out[out == ""] <- "low_priority"
  out
}
