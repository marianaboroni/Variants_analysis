# Shared low-level primitives for tumoronly.
# I/O helpers (read_variants, write_tsv, coalesce_columns, to_numeric_safe,
# cfg_get, is_missing_value) live in io.R and are reused here.

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- structured logging ------------------------------------------------------

#' Emit a structured, single-line log message.
#'
#' @param step short step identifier
#' @param msg human-readable message
#' @param ... named fields appended as key=value (e.g. n=1000, sample="S1")
#' @keywords internal
log_step <- function(step, msg, ...) {
  fields <- list(...)
  kv <- ""
  if (length(fields) > 0) {
    kv <- paste0(" | ", paste(sprintf("%s=%s", names(fields),
      vapply(fields, function(v) as.character(v)[1], character(1))), collapse = " "))
  }
  message(sprintf("[%s] %s%s", step, msg, kv))
}

# ---- checksums / hashing -----------------------------------------------------

#' SHA-256 of a file (streamed). Vectorized over `path` (one hash per file, in
#' the same order) so a cohort of input files hashes as cleanly as a single one.
#' @keywords internal
file_sha256 <- function(path) {
  if (is.null(path)) return(NA_character_)
  vapply(path, function(p) {
    if (is.na(p) || !nzchar(p) || !file.exists(p)) NA_character_
    else digest::digest(file = p, algo = "sha256")
  }, character(1), USE.NAMES = FALSE)
}

# ---- chromosome and allele normalization ------------------------------------

#' Normalize a chromosome name to a build-neutral form: strip the `chr` prefix,
#' map `M`/`chrM` to `MT`. Vectorized.
#' @keywords internal
normalize_chrom <- function(chrom) {
  x <- toupper(trimws(as.character(chrom)))
  x <- sub("^CHR", "", x)
  x[x %in% c("M", "MT")] <- "MT"
  x
}

#' Whether a caller FILTER value passes the technical gate. Defaults to
#' literal "PASS" (+ missing/"."/""), matching every caller unless configured
#' otherwise. A Mutect2 tumor-only analysis that has already restricted its
#' input to PASS/germline/panel_of_normals/combinations upstream can set
#' `hard_filters.caller_filter_accepted_values` so those values are not all
#' collapsed into a single "not PASS" technical failure (see
#' docs/FILTERING_STRATEGY.md, "Caller FILTER as graded evidence"). Vectorized;
#' shared by every call site that checks the caller FILTER so the allowlist is
#' consistent everywhere it is consulted.
#' @keywords internal
#' Whether a variant's (or locus') cohort-wide recurrence is trustworthy
#' evidence of a systematic artifact, rather than small-cohort noise. Requires
#' ALL of: the cohort itself large enough to make a fraction meaningful
#' (`cohort.min_cohort_size_for_recurrence`; below it, returns NA - "not
#' assessable at this cohort size", never silently FALSE or TRUE), an
#' absolute floor on how many samples actually carry it
#' (`cohort.min_recurrent_samples`/`min_recurrent_locus_samples` - a fraction
#' alone is meaningless at N=2, where one shared variant is already 50%), and
#' the existing fraction-of-cohort threshold
#' (`cohort.recurrent_variant_fraction_artifact`/`recurrent_locus_fraction_artifact`).
#' Vectorized; missing `total_samples`/`variant_n_samples`/`locus_n_samples`
#' columns (e.g. classification called outside the full pipeline, before
#' `merge_recurrence_features()`) are treated as "not evaluable", not as zero
#' recurrence. See docs/FILTERING_STRATEGY.md, "Cohort-size-adaptive recurrence".
#' @keywords internal
cohort_recurrent_flag <- function(x, cfg, level = c("variant", "locus")) {
  level <- match.arg(level)
  col <- function(name) if (name %in% names(x)) x[[name]] else rep(NA_real_, nrow(x))
  n <- col(if (level == "variant") "variant_n_samples" else "locus_n_samples")
  freq <- col(if (level == "variant") "variant_cohort_freq" else "locus_cohort_freq")
  total <- col("total_samples")

  min_n <- cfg_get(cfg, c("cohort",
    if (level == "variant") "min_recurrent_samples" else "min_recurrent_locus_samples"), 3)
  frac_thresh <- cfg_get(cfg, c("cohort",
    if (level == "variant") "recurrent_variant_fraction_artifact" else "recurrent_locus_fraction_artifact"),
    if (level == "variant") 0.30 else 0.40)
  min_cohort <- cfg_get(cfg, c("cohort", "min_cohort_size_for_recurrence"), 10)

  evaluable <- !is.na(total) & total >= min_cohort
  flag <- evaluable & !is.na(n) & n >= min_n & !is.na(freq) & freq >= frac_thresh
  ifelse(evaluable, flag, NA)
}

caller_filter_pass <- function(filter_status, cfg = NULL) {
  accepted <- cfg_get(cfg, c("hard_filters", "caller_filter_accepted_values"), "PASS")
  is.na(filter_status) | filter_status %in% c(accepted, ".", "NA") | filter_status == ""
}

#' Trim shared leading/trailing bases from a REF/ALT pair and left-normalize the
#' position, matching bcftools/VEP-style minimal representation. Purely local
#' (no reference sequence); handles SNVs, insertions and deletions. Vectorized
#' across ALL rows per trim step (not a per-row loop): each step trims one
#' character from every row that still needs it, so cost is
#' O(longest_allele_in_the_whole_input) vectorized passes, not O(n) R-level
#' loop iterations - the difference between seconds and hours at COSMIC-v104
#' scale (~10^8 rows). Semantically identical to trimming one row fully before
#' moving to the next (suffix trim always completes before prefix trim starts).
#'
#' @return a data.frame with normalized `pos`, `ref`, `alt`.
#' @keywords internal
normalize_alleles <- function(pos, ref, alt) {
  pos <- as.integer(pos)
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  out_pos <- pos; out_ref <- ref; out_alt <- alt
  valid <- !is.na(out_ref) & !is.na(out_alt) & out_ref != "" & out_alt != ""

  max_len <- suppressWarnings(max(nchar(ref), nchar(alt), na.rm = TRUE))
  if (!is.finite(max_len)) max_len <- 0L

  # trim shared suffix (keep at least 1 base on each side)
  for (i in seq_len(max_len)) {
    nr <- nchar(out_ref); nalt <- nchar(out_alt)
    idx <- which(valid & nr > 1L & nalt > 1L)
    if (length(idx) == 0L) break
    same <- substr(out_ref[idx], nr[idx], nr[idx]) == substr(out_alt[idx], nalt[idx], nalt[idx])
    if (!any(same)) break
    idx <- idx[same]
    out_ref[idx] <- substr(out_ref[idx], 1L, nchar(out_ref[idx]) - 1L)
    out_alt[idx] <- substr(out_alt[idx], 1L, nchar(out_alt[idx]) - 1L)
  }
  # trim shared prefix, advancing position
  for (i in seq_len(max_len)) {
    nr <- nchar(out_ref); nalt <- nchar(out_alt)
    idx <- which(valid & nr > 1L & nalt > 1L)
    if (length(idx) == 0L) break
    same <- substr(out_ref[idx], 1L, 1L) == substr(out_alt[idx], 1L, 1L)
    if (!any(same)) break
    idx <- idx[same]
    out_ref[idx] <- substr(out_ref[idx], 2L, nchar(out_ref[idx]))
    out_alt[idx] <- substr(out_alt[idx], 2L, nchar(out_alt[idx]))
    out_pos[idx] <- out_pos[idx] + 1L
  }
  data.frame(pos = out_pos, ref = out_ref, alt = out_alt, stringsAsFactors = FALSE)
}

#' Build the canonical variant key `BUILD|CHROM|POS|REF|ALT` with normalized
#' chrom and alleles. Vectorized.
#' @keywords internal
canonical_variant_key <- function(build, chrom, pos, ref, alt) {
  chrom <- normalize_chrom(chrom)
  na <- normalize_alleles(pos, ref, alt)
  paste(build, chrom, na$pos, na$ref, na$alt, sep = "|")
}

#' Split multiallelic ALT (comma-separated) into one row per ALT allele.
#' Returns a data.frame; requires a column named `alt`.
#' @keywords internal
split_multiallelic <- function(df, alt_col = "alt") {
  alt <- as.character(df[[alt_col]])
  has_multi <- grepl(",", alt, fixed = TRUE)
  if (!any(has_multi)) return(df)
  parts <- strsplit(alt, ",", fixed = TRUE)
  reps <- lengths(parts)
  out <- df[rep(seq_len(nrow(df)), reps), , drop = FALSE]
  out[[alt_col]] <- unlist(parts, use.names = FALSE)
  rownames(out) <- NULL
  out
}

# ---- directory lock ----------------------------------------------------------

#' Acquire a directory-based lock (portable, no filelock dependency).
#' Uses atomic dir.create. Blocks up to `timeout` seconds.
#' @return a lock handle (the lock dir path) to pass to release_lock().
#' @keywords internal
acquire_lock <- function(lock_dir, timeout = 300, poll = 0.5) {
  deadline <- Sys.time() + timeout
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE)) {
      writeLines(as.character(Sys.getpid()), file.path(lock_dir, "pid"))
      return(lock_dir)
    }
    if (Sys.time() > deadline) {
      stop(sprintf("Could not acquire build lock at %s within %d s. If no other build is running, remove the directory manually.", lock_dir, timeout))
    }
    Sys.sleep(poll)
  }
}

#' @keywords internal
release_lock <- function(handle) {
  if (!is.null(handle) && dir.exists(handle)) unlink(handle, recursive = TRUE, force = TRUE)
  invisible(NULL)
}

#' Package version string, robust to dev-loading (no installed DESCRIPTION).
#' @keywords internal
tumoronly_version <- function() {
  v <- tryCatch(as.character(utils::packageVersion("tumoronly")), error = function(e) NA_character_)
  if (is.na(v)) "0.1.0-dev" else v
}
