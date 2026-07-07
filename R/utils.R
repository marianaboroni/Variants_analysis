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

#' SHA-256 of a file (streamed).
#' @keywords internal
file_sha256 <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
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

#' Trim shared leading/trailing bases from a REF/ALT pair and left-normalize the
#' position, matching bcftools/VEP-style minimal representation. Purely local
#' (no reference sequence); handles SNVs, insertions and deletions. Vectorized.
#'
#' @return a data.frame with normalized `pos`, `ref`, `alt`.
#' @keywords internal
normalize_alleles <- function(pos, ref, alt) {
  pos <- as.integer(pos)
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  n <- length(ref)
  out_pos <- pos; out_ref <- ref; out_alt <- alt
  for (i in seq_len(n)) {
    r <- out_ref[i]; a <- out_alt[i]; p <- out_pos[i]
    if (is.na(r) || is.na(a) || r == "" || a == "") next
    # trim shared suffix (keep at least 1 base on each side)
    while (nchar(r) > 1 && nchar(a) > 1 &&
           substr(r, nchar(r), nchar(r)) == substr(a, nchar(a), nchar(a))) {
      r <- substr(r, 1, nchar(r) - 1)
      a <- substr(a, 1, nchar(a) - 1)
    }
    # trim shared prefix, advancing position
    while (nchar(r) > 1 && nchar(a) > 1 &&
           substr(r, 1, 1) == substr(a, 1, 1)) {
      r <- substr(r, 2, nchar(r))
      a <- substr(a, 2, nchar(a))
      p <- p + 1L
    }
    out_ref[i] <- r; out_alt[i] <- a; out_pos[i] <- p
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
