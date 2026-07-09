# Human review system. Reviews are the ONLY source of confirmed labels
# (TRUE_POSITIVE / FALSE_POSITIVE); the pipeline never derives them from
# PASS/FAIL. Storage is append-only: history is never overwritten.

REVIEW_LABELS <- c("TRUE_POSITIVE", "FALSE_POSITIVE", "UNCERTAIN", "NOT_REVIEWED")
REVIEW_TEMPLATE_COLS <- c("sample_id", "variant_id", "gene", "consequence",
  "STATUS", "PRIORITY", "AUTHENTICITY", "VARIANT_REVIEW_SCORE",
  "review_label", "review_confidence", "review_reason", "reviewer",
  "review_date", "validation_method", "validation_result", "notes")

#' Write a human-review template for a run's variants (pre-filled with current
#' labels; review columns left blank for the reviewer to fill).
#' @keywords internal
export_review_template <- function(variants, path) {
  keep <- intersect(c("sample_id", "variant_id", "gene", "consequence", "STATUS",
    "PRIORITY", "AUTHENTICITY", "VARIANT_REVIEW_SCORE"), names(variants))
  tmpl <- variants[, keep, drop = FALSE]
  tmpl$review_label <- "NOT_REVIEWED"
  for (c0 in c("review_confidence", "review_reason", "reviewer", "review_date",
               "validation_method", "validation_result", "notes")) tmpl[[c0]] <- NA_character_
  # prioritize dubious variants at the top of the template
  if ("VARIANT_REVIEW_SCORE" %in% names(tmpl)) tmpl <- tmpl[order(-tmpl$VARIANT_REVIEW_SCORE), ]
  write_tsv(tmpl, path)
  invisible(path)
}

#' Import human variant reviews into the append-only evidence database.
#'
#' Confirmed labels come ONLY from the review file; `PASS`/`FAIL` are never
#' auto-converted to `TRUE_POSITIVE`/`FALSE_POSITIVE`. History is append-only and
#' never overwritten; `reviewed_variants` keeps the latest review per variant.
#'
#' @param file review TSV (from `export-review-template`, filled in).
#' @param db_dir evidence database directory (default `db/variant_evidence`).
#' @param reviewer optional default reviewer if absent in the file.
#' @return (invisibly) a summary list (counts by label, n_appended, paths).
#' @export
#' @examples
#' \dontrun{ import_variant_reviews("reviews.tsv") }
import_variant_reviews <- function(file, db_dir = "db/variant_evidence", reviewer = NULL) {
  if (!file.exists(file)) stop("Review file not found: ", file, call. = FALSE)
  ev <- ensure_variant_evidence_db(db_dir)
  rv <- read_variants(file, "\t")

  # variant key
  if (!"variant_id" %in% names(rv)) {
    need <- c("chrom", "pos", "ref", "alt"); need <- need[!need %in% tolower(names(rv))]
    if (length(need)) stop("Review file needs `variant_id` or chrom/pos/ref/alt columns.", call. = FALSE)
  }
  rv$variant_id <- as.character(coalesce_columns(rv, c("variant_id", "VARIANT_ID")))
  rv$sample_id <- as.character(coalesce_columns(rv, c("sample_id", "Sample", "Tumor_Sample_Barcode"), default = NA_character_))
  rv$review_label <- toupper(trimws(as.character(coalesce_columns(rv, c("review_label", "label")))))
  rv$reviewer <- as.character(coalesce_columns(rv, c("reviewer", "Reviewer"), default = reviewer %||% NA_character_))
  rv$review_date <- as.character(coalesce_columns(rv, c("review_date", "date"), default = NA_character_))
  rv$review_reason <- as.character(coalesce_columns(rv, c("review_reason", "reason"), default = NA_character_))
  rv$validation_method <- as.character(coalesce_columns(rv, c("validation_method"), default = NA_character_))

  # only rows with an actual (non-default) label are ingested
  rv <- rv[!is.na(rv$review_label) & rv$review_label != "" & rv$review_label != "NOT_REVIEWED", , drop = FALSE]
  if (nrow(rv) == 0) { log_step("review", "no reviewed rows to import (all NOT_REVIEWED)"); return(invisible(list(n_appended = 0))) }

  bad <- !rv$review_label %in% REVIEW_LABELS
  if (any(bad)) stop(sprintf("Invalid review_label(s): %s. Allowed: %s",
    paste(unique(rv$review_label[bad]), collapse = ", "), paste(REVIEW_LABELS, collapse = ", ")), call. = FALSE)
  # confirmed labels require human evidence (reviewer)
  confirmed <- rv$review_label %in% c("TRUE_POSITIVE", "FALSE_POSITIVE")
  if (any(confirmed & (is.na(rv$reviewer) | rv$reviewer == "")))
    stop("TRUE_POSITIVE / FALSE_POSITIVE reviews require a `reviewer` (human/independent evidence).", call. = FALSE)

  rv$variant_key <- ifelse(!is.na(rv$sample_id) & rv$sample_id != "",
                           paste(rv$sample_id, rv$variant_id, sep = "|"), rv$variant_id)
  rv$import_timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  keep <- c("variant_key", "sample_id", "variant_id", "review_label", "review_reason",
            "reviewer", "review_date", "validation_method", "import_timestamp")
  new_rows <- rv[, intersect(keep, names(rv)), drop = FALSE]

  # append-only history (store as text; coerce to character so re-read columns
  # with all-NA -> logical never clash with character on append)
  new_rows[] <- lapply(new_rows, as.character)
  hist_path <- file.path(db_dir, "review_history.tsv.gz")
  if (file.exists(hist_path)) {
    old <- read_variants(hist_path, "\t"); old[] <- lapply(old, as.character)
    hist <- data.table::rbindlist(list(old, new_rows), fill = TRUE, use.names = TRUE)
  } else hist <- new_rows
  write_tsv(as.data.frame(hist), hist_path)

  # reviewed_variants = latest review per key (by import order in history)
  h <- as.data.frame(hist)
  latest <- h[!duplicated(h$variant_key, fromLast = TRUE), , drop = FALSE]
  write_tsv(latest, file.path(db_dir, "reviewed_variants.tsv.gz"))

  summ <- list(n_appended = nrow(new_rows), n_total_history = nrow(h),
               by_label = as.list(table(new_rows$review_label)),
               history = hist_path)
  log_step("review", "reviews imported (append-only)", appended = nrow(new_rows),
           total_history = nrow(h))
  invisible(summ)
}

#' Ensure the append-only variant-evidence database layout exists.
#' @keywords internal
ensure_variant_evidence_db <- function(db_dir = "db/variant_evidence") {
  dir.create(file.path(db_dir, "model_registry"), recursive = TRUE, showWarnings = FALSE)
  schema <- file.path(db_dir, "feature_schema.json")
  if (!file.exists(schema))
    jsonlite::write_json(list(
      technical = c("dp", "alt_count", "vaf", "tlod", "mbq", "mmq", "strand_artifact",
                    "orientation_bias", "pon_flag", "clustered_events", "weak_evidence"),
      biological = c("impact_rank", "COMPUTATIONAL_EVIDENCE_SCORE", "hotspot_match",
                     "driver_score", "COSMIC_GLOBAL_RECURRENCE_SCORE", "COSMIC_CONTEXT_SUPPORT_SCORE")),
      schema, auto_unbox = TRUE, pretty = TRUE)
  ld <- file.path(db_dir, "label_dictionary.yml")
  if (!file.exists(ld))
    writeLines(c("# Confirmed labels come ONLY from human/independent review.",
      "TRUE_POSITIVE: confirmed real variant",
      "FALSE_POSITIVE: confirmed artifact / not real",
      "UNCERTAIN: reviewed, inconclusive",
      "NOT_REVIEWED: not yet reviewed"), ld)
  db_dir
}
