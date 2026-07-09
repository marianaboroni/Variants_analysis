# Supervised continuous-learning scaffold (EXPERIMENTAL). Two SEPARATE objectives
# are kept distinct and never merged into one opaque score:
#   P_TRUE_VARIANT        -- technical: real call vs artifact (technical features only)
#   P_BIOLOGICAL_RELEVANCE-- biological relevance (biological features only)
# Models NEVER alter filter_status; they only rank / suggest review. Training uses
# ONLY human-reviewed labels, group-aware CV (no row-level leakage), explicit
# approval to activate, and versioned rollback. Interpretable models only.

ML_TECHNICAL_FEATURES <- c("dp", "alt_count", "vaf", "tlod", "mbq", "mmq",
  "strand_artifact", "orientation_bias", "pon_flag", "clustered_events", "weak_evidence")
ML_BIOLOGICAL_FEATURES <- c("impact_rank", "COMPUTATIONAL_EVIDENCE_SCORE", "hotspot_match",
  "driver_score", "COSMIC_GLOBAL_RECURRENCE_SCORE", "COSMIC_CONTEXT_SUPPORT_SCORE")

#' Report the state of the training database (label counts, readiness).
#' @keywords internal
ml_training_status <- function(db_dir = "db/variant_evidence") {
  rp <- file.path(db_dir, "reviewed_variants.tsv.gz")
  if (!file.exists(rp)) { cat("No reviewed variants yet. Import reviews first.\n"); return(invisible(list(n = 0))) }
  rv <- read_variants(rp, "\t")
  tab <- table(rv$review_label)
  n_tp <- sum(rv$review_label == "TRUE_POSITIVE", na.rm = TRUE)
  n_fp <- sum(rv$review_label == "FALSE_POSITIVE", na.rm = TRUE)
  min_each <- 20
  ready <- n_tp >= min_each && n_fp >= min_each
  cat("Training database status\n========================\n")
  print(tab)
  cat(sprintf("\nTRUE_POSITIVE=%d  FALSE_POSITIVE=%d  (min %d each to train)\n", n_tp, n_fp, min_each))
  cat(sprintf("Ready to train P_TRUE_VARIANT: %s\n", if (ready) "YES" else "NO"))
  models <- list.files(file.path(db_dir, "model_registry"), pattern = "[.]rds$")
  active <- file.path(db_dir, "model_registry", "ACTIVE")
  cat(sprintf("Candidate models: %d  | Active: %s\n", length(models),
              if (file.exists(active)) readLines(active)[1] else "none"))
  invisible(list(n = nrow(rv), n_tp = n_tp, n_fp = n_fp, ready = ready))
}

#' Build a training set by joining human labels to run features (no leakage: join
#' by variant key; group column = sample_id).
#' @keywords internal
ml_build_training_set <- function(db_dir = "db/variant_evidence", run_dir = NULL,
                                  objective = "P_TRUE_VARIANT") {
  rv <- read_variants(file.path(db_dir, "reviewed_variants.tsv.gz"), "\t")
  rv <- rv[rv$review_label %in% c("TRUE_POSITIVE", "FALSE_POSITIVE"), , drop = FALSE]
  if (nrow(rv) == 0) stop("No confirmed TRUE_POSITIVE/FALSE_POSITIVE labels to build a training set.", call. = FALSE)
  if (is.null(run_dir)) stop("Provide run_dir with tables/variants_all.tsv.gz to join features.", call. = FALSE)
  va <- read_variants(file.path(run_dir, "tables", "variants_all.tsv.gz"), "\t")
  va$variant_key <- ifelse(!is.na(va$sample_id) & va$sample_id != "",
                           paste(va$sample_id, va$variant_id, sep = "|"), va$variant_id)
  feats <- if (objective == "P_BIOLOGICAL_RELEVANCE") ML_BIOLOGICAL_FEATURES else ML_TECHNICAL_FEATURES
  feats <- intersect(feats, names(va))
  m <- merge(rv[, c("variant_key", "review_label")], va[, c("variant_key", feats), drop = FALSE],
             by = "variant_key")
  m$y <- as.integer(m$review_label == "TRUE_POSITIVE")
  m$group <- sub("[|].*$", "", m$variant_key)   # sample_id for group-aware CV
  out <- file.path(db_dir, sprintf("training_set_%s.tsv.gz", objective))
  write_tsv(m, out)
  log_step("ml", "training set built", objective = objective, n = nrow(m),
           tp = sum(m$y == 1), fp = sum(m$y == 0), features = length(feats))
  invisible(list(path = out, n = nrow(m), features = feats))
}

#' Train an interpretable (regularization-free logistic) candidate model with
#' group-aware (by-sample) cross-validation. Writes a candidate to the registry;
#' does NOT activate it.
#' @keywords internal
ml_train_model <- function(db_dir = "db/variant_evidence", objective = "P_TRUE_VARIANT",
                           model_id = NULL) {
  ts <- file.path(db_dir, sprintf("training_set_%s.tsv.gz", objective))
  if (!file.exists(ts)) stop("No training set. Run build-training-set first.", call. = FALSE)
  d <- read_variants(ts, "\t")
  feats <- setdiff(names(d), c("variant_key", "review_label", "y", "group"))
  for (f in feats) d[[f]] <- suppressWarnings(as.numeric(d[[f]]))
  d[is.na(d)] <- 0
  if (length(unique(d$y)) < 2) stop("Training set needs both classes (TP and FP).", call. = FALSE)

  form <- stats::as.formula(paste("y ~", paste(feats, collapse = " + ")))
  # group-aware CV (leave-one-sample-group-out, capped at 5 folds)
  groups <- unique(d$group); set.seed(1)
  folds <- split(groups, seq_along(groups) %% min(5, length(groups)))
  preds <- rep(NA_real_, nrow(d))
  for (f in folds) {
    te <- d$group %in% f; if (all(te) || !any(te)) next
    fit <- suppressWarnings(stats::glm(form, data = d[!te, ], family = stats::binomial()))
    preds[te] <- suppressWarnings(stats::predict(fit, d[te, ], type = "response"))
  }
  ok <- !is.na(preds)
  metrics <- ml_metrics(d$y[ok], preds[ok])
  final <- suppressWarnings(stats::glm(form, data = d, family = stats::binomial()))

  mid <- model_id %||% paste0(objective, "_", substr(digest::digest(list(feats, nrow(d))), 1, 8))
  reg <- file.path(db_dir, "model_registry")
  saveRDS(list(objective = objective, features = feats, model = final,
               metrics = metrics, n = nrow(d), trained = format(Sys.time()),
               package_version = tumoronly_version(), status = "candidate"),
          file.path(reg, paste0(mid, ".rds")))
  jsonlite::write_json(c(list(model_id = mid), metrics), file.path(reg, paste0(mid, ".metrics.json")),
                       auto_unbox = TRUE, pretty = TRUE)
  log_step("ml", "candidate model trained (NOT active)", model_id = mid,
           auc = round(metrics$roc_auc, 3), n = nrow(d))
  invisible(list(model_id = mid, metrics = metrics))
}

#' @keywords internal
ml_metrics <- function(y, p) {
  # ROC-AUC via Mann-Whitney
  pos <- p[y == 1]; neg <- p[y == 0]
  auc <- if (length(pos) && length(neg)) mean(outer(pos, neg, ">") + 0.5 * outer(pos, neg, "==")) else NA_real_
  pred <- as.integer(p >= 0.5)
  tp <- sum(pred == 1 & y == 1); fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1); tn <- sum(pred == 0 & y == 0)
  list(roc_auc = auc,
       precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
       recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
       specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_,
       brier = mean((p - y)^2),
       confusion = list(tp = tp, fp = fp, fn = fn, tn = tn), n = length(y))
}

#' @keywords internal
ml_evaluate_model <- function(db_dir = "db/variant_evidence", model_id) {
  mp <- file.path(db_dir, "model_registry", paste0(model_id, ".rds"))
  if (!file.exists(mp)) stop("Model not found: ", model_id, call. = FALSE)
  m <- readRDS(mp); print(m$metrics); invisible(m$metrics)
}

#' Activate a candidate model. Requires explicit approval and minimum evidence;
#' preserves the previous active model for rollback.
#' @keywords internal
ml_activate_model <- function(db_dir = "db/variant_evidence", model_id, approve = FALSE,
                              min_labels = 40) {
  reg <- file.path(db_dir, "model_registry")
  mp <- file.path(reg, paste0(model_id, ".rds"))
  if (!file.exists(mp)) stop("Model not found: ", model_id, call. = FALSE)
  m <- readRDS(mp)
  if (m$n < min_labels) stop(sprintf("Refusing to activate: only %d labeled variants (min %d).", m$n, min_labels), call. = FALSE)
  if (!isTRUE(approve)) stop("Activation requires explicit approval (approve = TRUE / --approve).", call. = FALSE)
  active <- file.path(reg, "ACTIVE")
  if (file.exists(active)) writeLines(readLines(active)[1], file.path(reg, "PREVIOUS_ACTIVE"))
  writeLines(model_id, active)
  log_step("ml", "model activated (experimental; does not alter filter_status)", model_id = model_id)
  invisible(model_id)
}

#' @keywords internal
ml_rollback_model <- function(db_dir = "db/variant_evidence") {
  reg <- file.path(db_dir, "model_registry")
  prev <- file.path(reg, "PREVIOUS_ACTIVE")
  if (!file.exists(prev)) stop("No previous active model to roll back to.", call. = FALSE)
  writeLines(readLines(prev)[1], file.path(reg, "ACTIVE"))
  log_step("ml", "rolled back to previous active model", model_id = readLines(prev)[1])
  invisible(readLines(prev)[1])
}

#' Textual explanation for a single prediction (top contributing features).
#' @keywords internal
ml_explain_prediction <- function(model, feature_row) {
  co <- stats::coef(model$model); co <- co[names(co) != "(Intercept)"]
  x <- unlist(feature_row[names(co)]); x[is.na(x)] <- 0
  contrib <- co * x
  ord <- order(-abs(contrib))
  top <- head(names(contrib)[ord], 4)
  sprintf("Driven by: %s", paste(sprintf("%s(%+.2f)", top, contrib[top]), collapse = ", "))
}
