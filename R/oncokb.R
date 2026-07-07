# OncoKB annotation as a SEPARATE, post-hoc step. It reads frozen run results,
# adds ONCOKB_* columns, and writes a standalone annotation table. It NEVER
# modifies filter_status / filter_reasons and never rewrites variants_all
# (see docs/REFACTOR_AUDIT.md RISK-1 and the brief section 12).
#
# The API token comes ONLY from the ONCOKB_TOKEN environment variable and is
# never logged, serialized, or included in outputs or errors.

ONCOKB_COLS <- c("ONCOKB_ANNOTATED", "ONCOKB_ONCOGENIC", "ONCOKB_MUTATION_EFFECT",
                 "ONCOKB_HIGHEST_LEVEL", "ONCOKB_HIGHEST_SENSITIVE_LEVEL",
                 "ONCOKB_HIGHEST_RESISTANCE_LEVEL", "ONCOKB_KNOWN_EFFECT",
                 "ONCOKB_QUERY_DATE", "ONCOKB_QUERY_STATUS")

#' Annotate a completed run with OncoKB (optional, non-filtering).
#'
#' Reads `tables/variants_all.tsv.gz` from `run_dir`, annotates the selected
#' scope against OncoKB, and writes `tables/oncokb_annotations.tsv.gz`. The
#' primary run results are never modified; if the API is unavailable the step
#' still completes and marks variants as not annotated.
#'
#' @param run_dir a completed run directory (from [run_tumor_only()]).
#' @param scope "retained" (default) or "all".
#' @param cfg optional resolved config (for cache_dir/endpoint); read from the
#'   run's `config.resolved.yml` when omitted.
#' @return the annotation data.frame (invisibly).
#' @export
#' @examples
#' \dontrun{
#' Sys.setenv(ONCOKB_TOKEN = "...")
#' annotate_oncokb("results/example_run")
#' }
annotate_oncokb <- function(run_dir, scope = c("retained", "all"), cfg = NULL) {
  scope <- match.arg(scope)
  all_path <- file.path(run_dir, "tables", "variants_all.tsv.gz")
  if (!file.exists(all_path))
    stop(sprintf("No frozen results at %s. Run `tumoronly run` first.", all_path), call. = FALSE)
  if (is.null(cfg)) {
    cpath <- file.path(run_dir, "config.resolved.yml")
    cfg <- if (file.exists(cpath)) read_config(cpath) else list()
  }
  variants <- read_variants(all_path, "\t")

  in_scope <- if (scope == "all") rep(TRUE, nrow(variants))
              else variants$filter_status %in% c("PASS", "REVIEW")
  build <- cfg_get(cfg, c("input", "genome_build"), "GRCh38")

  ann <- init_oncokb_columns(variants)
  ann$ONCOKB_QUERY_STATUS[!in_scope] <- "not_in_scope"

  token <- Sys.getenv("ONCOKB_TOKEN", unset = "")
  cache_dir <- cfg_get(cfg, c("oncokb", "cache_dir"), file.path(run_dir, "..", ".oncokb_cache"))
  endpoint <- cfg_get(cfg, c("oncokb", "endpoint"), "https://www.oncokb.org/api/v1")

  idx <- which(in_scope)
  if (length(idx) > 0) {
    res <- oncokb_annotate_variants(variants[idx, , drop = FALSE], build, token, cache_dir, endpoint)
    for (col in ONCOKB_COLS) ann[[col]][idx] <- res[[col]]
  }

  out_cols <- c(intersect(c("sample_id", "variant_id", "chrom", "pos", "ref", "alt",
                            "gene", "filter_status"), names(ann)), ONCOKB_COLS)
  out <- ann[, out_cols, drop = FALSE]
  out_path <- file.path(run_dir, "tables", "oncokb_annotations.tsv.gz")
  write_tsv(out, out_path)

  # SAFETY: prove filter_status is untouched relative to the frozen file.
  frozen_status <- read_variants(all_path, "\t")$filter_status
  stopifnot(identical(frozen_status, variants$filter_status))
  log_step("oncokb", "OncoKB annotation written (filter_status unchanged)",
           scope = scope, annotated = sum(ann$ONCOKB_ANNOTATED, na.rm = TRUE),
           in_scope = sum(in_scope))
  invisible(out)
}

#' @keywords internal
init_oncokb_columns <- function(x) {
  x$ONCOKB_ANNOTATED <- FALSE
  x$ONCOKB_ONCOGENIC <- NA_character_
  x$ONCOKB_MUTATION_EFFECT <- NA_character_
  x$ONCOKB_HIGHEST_LEVEL <- NA_character_
  x$ONCOKB_HIGHEST_SENSITIVE_LEVEL <- NA_character_
  x$ONCOKB_HIGHEST_RESISTANCE_LEVEL <- NA_character_
  x$ONCOKB_KNOWN_EFFECT <- NA_character_
  x$ONCOKB_QUERY_DATE <- NA_character_
  x$ONCOKB_QUERY_STATUS <- "pending"
  x
}

#' Annotate a set of in-scope variants, using a per-variant on-disk cache and
#' (when possible) the OncoKB genomic-change API. Returns a data.frame with the
#' ONCOKB_* columns in row order.
#' @keywords internal
oncokb_annotate_variants <- function(v, build, token, cache_dir, endpoint) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  keys <- canonical_variant_key(build, v$chrom, v$pos, v$ref, v$alt)
  query_date <- format(Sys.Date())
  res <- init_oncokb_columns(v)

  have_api <- nzchar(token) && requireNamespace("httr2", quietly = TRUE)
  to_query <- integer()
  for (i in seq_len(nrow(v))) {
    cached <- oncokb_cache_get(cache_dir, keys[i])
    if (!is.null(cached)) {
      for (col in ONCOKB_COLS) res[[col]][i] <- cached[[col]] %||% res[[col]][i]
    } else if (have_api) {
      to_query <- c(to_query, i)
    } else {
      res$ONCOKB_QUERY_STATUS[i] <- if (!nzchar(token)) "no_token" else "no_api_client"
    }
  }

  if (length(to_query) > 0) {
    fetched <- tryCatch(
      oncokb_api_batch(v[to_query, , drop = FALSE], build, token, endpoint),
      error = function(e) {
        log_step("oncokb", "OncoKB API unavailable; leaving variants not annotated",
                 error = conditionMessage(e))
        NULL
      })
    for (j in seq_along(to_query)) {
      i <- to_query[j]
      if (is.null(fetched)) { res$ONCOKB_QUERY_STATUS[i] <- "api_unavailable"; next }
      rec <- fetched[[j]]
      rec$ONCOKB_QUERY_DATE <- query_date
      for (col in ONCOKB_COLS) res[[col]][i] <- rec[[col]] %||% res[[col]][i]
      oncokb_cache_put(cache_dir, keys[i], rec)
    }
  }
  res[, ONCOKB_COLS, drop = FALSE]
}

#' @keywords internal
oncokb_cache_path <- function(cache_dir, key) {
  file.path(cache_dir, paste0(digest::digest(key, algo = "sha256"), ".json"))
}
#' @keywords internal
oncokb_cache_get <- function(cache_dir, key) {
  p <- oncokb_cache_path(cache_dir, key)
  if (!file.exists(p)) return(NULL)
  tryCatch(jsonlite::read_json(p, simplifyVector = TRUE), error = function(e) NULL)
}
#' @keywords internal
oncokb_cache_put <- function(cache_dir, key, rec) {
  jsonlite::write_json(rec, oncokb_cache_path(cache_dir, key), auto_unbox = TRUE, null = "null")
}

#' Batch call to the OncoKB genomic-change endpoint with retries + backoff.
#' Requires httr2 and a token. Never logs the token. Returns a list of ONCOKB_*
#' record lists, one per input variant.
#' @keywords internal
oncokb_api_batch <- function(v, build, token, endpoint, max_retries = 4, timeout = 30) {
  body <- lapply(seq_len(nrow(v)), function(i) {
    list(genomicLocation = sprintf("%s,%s,%s,%s,%s",
      normalize_chrom(v$chrom[i]), v$pos[i], v$pos[i] + nchar(v$ref[i]) - 1L,
      v$ref[i], v$alt[i]),
      referenceGenome = build)
  })
  url <- paste0(sub("/$", "", endpoint), "/annotate/mutations/byGenomicChange")
  attempt <- 0
  repeat {
    attempt <- attempt + 1
    resp <- tryCatch({
      req <- httr2::request(url)
      req <- httr2::req_headers(req, Authorization = paste("Bearer", token),
                                `Content-Type` = "application/json")
      req <- httr2::req_body_json(req, body)
      req <- httr2::req_timeout(req, timeout)
      httr2::req_perform(req)
    }, error = function(e) e)
    if (!inherits(resp, "error")) {
      code <- httr2::resp_status(resp)
      if (code == 200) {
        parsed <- httr2::resp_body_json(resp, simplifyVector = FALSE)
        return(lapply(parsed, oncokb_parse_record))
      }
      if (code != 429 && code < 500) stop(sprintf("OncoKB API returned HTTP %d", code))
    }
    if (attempt > max_retries) stop("OncoKB API failed after retries.")
    Sys.sleep(min(2^attempt, 30))  # exponential backoff, capped
  }
}

#' Parse one OncoKB API record into the ONCOKB_* schema.
#' @keywords internal
oncokb_parse_record <- function(rec) {
  g <- function(k) { v <- rec[[k]]; if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1] }
  list(
    ONCOKB_ANNOTATED = isTRUE(rec$geneExist) || !is.null(rec$oncogenic),
    ONCOKB_ONCOGENIC = g("oncogenic"),
    ONCOKB_MUTATION_EFFECT = if (!is.null(rec$mutationEffect)) g0(rec$mutationEffect$knownEffect) else NA_character_,
    ONCOKB_HIGHEST_LEVEL = g("highestSensitiveLevel") %||% NA_character_,
    ONCOKB_HIGHEST_SENSITIVE_LEVEL = g("highestSensitiveLevel"),
    ONCOKB_HIGHEST_RESISTANCE_LEVEL = g("highestResistanceLevel"),
    ONCOKB_KNOWN_EFFECT = if (!is.null(rec$mutationEffect)) g0(rec$mutationEffect$knownEffect) else NA_character_,
    ONCOKB_QUERY_DATE = format(Sys.Date()),
    ONCOKB_QUERY_STATUS = "annotated"
  )
}
#' @keywords internal
g0 <- function(v) if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
