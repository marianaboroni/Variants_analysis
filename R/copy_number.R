# Per-variant copy number from copy-number segment files.
#
# Each sample may point (via the `copy_number_file` column of
# `input.sample_metadata`) to a tab-separated segment table such as PURPLE's
# `<sample>.purple.cnv.somatic.tsv`. Variants are matched to the segment that
# contains their position; the segment's total copy number (and allele-specific
# minor/major copy number when present) is returned for use by the clonality
# module. Unreadable files or unmatched variants yield NA, so callers fall back
# to the configured copy-neutral default.

#' Look up segment copy number for variant rows.
#'
#' @param variants Variant table with `chrom`, `pos` and `copy_number_file`.
#' @param rows Row indices to annotate (others are returned as NA).
#' @return A list of numeric vectors `total_cn`, `minor_cn`, `major_cn`, each
#'   of length `nrow(variants)`.
#' @keywords internal
segment_copy_number_for_variants <- function(variants, rows = seq_len(nrow(variants))) {
  n <- nrow(variants)
  out <- list(total_cn = rep(NA_real_, n), minor_cn = rep(NA_real_, n), major_cn = rep(NA_real_, n))
  if (!all(c("copy_number_file", "chrom", "pos") %in% names(variants)) || length(rows) == 0) {
    return(out)
  }
  files <- as.character(variants$copy_number_file)
  pos <- suppressWarnings(as.numeric(variants$pos))
  rows <- rows[!is_missing_value(files[rows]) & !is.na(pos[rows]) & !is_missing_value(variants$chrom[rows])]
  if (length(rows) == 0) return(out)

  for (idx in split(rows, files[rows])) {
    path <- files[[idx[[1]]]]
    seg <- read_copy_number_segments(path)
    if (is.null(seg) || nrow(seg) == 0) next
    q <- data.table::data.table(
      row_i = idx,
      q_chrom = normalize_cn_chrom(variants$chrom[idx]),
      q_pos = pos[idx]
    )
    hit <- seg[q, on = .(seg_chrom == q_chrom, seg_start <= q_pos, seg_end >= q_pos),
               mult = "first", nomatch = NULL,
               .(row_i = i.row_i, total_cn = x.total_cn, minor_cn = x.minor_cn, major_cn = x.major_cn)]
    if (nrow(hit) == 0) next
    out$total_cn[hit$row_i] <- hit$total_cn
    out$minor_cn[hit$row_i] <- hit$minor_cn
    out$major_cn[hit$row_i] <- hit$major_cn
  }
  out
}

#' Read a copy-number segment table (PURPLE, ASCAT-like or generic).
#'
#' Recognised columns: chromosome (`chromosome`/`chrom`/`chr`), start
#' (`start`/`startpos`/`loc.start`), end (`end`/`endpos`/`loc.end`), total copy
#' number (`copyNumber`/`total_cn`/`tcn`/`CN`), and optionally allele-specific
#' copy number (`minorAlleleCopyNumber`/`minor_cn`/`nMinor`,
#' `majorAlleleCopyNumber`/`major_cn`/`nMajor`).
#'
#' @param path Path to a tab-separated segment file (may be gzipped).
#' @return A data.table with `seg_chrom`, `seg_start`, `seg_end`, `total_cn`,
#'   `minor_cn`, `major_cn`, or NULL when the file is missing or unusable.
#' @keywords internal
read_copy_number_segments <- function(path) {
  if (!file.exists(path)) {
    log_step("clonality", "copy-number segment file not found; using copy-neutral default", file = path)
    return(NULL)
  }
  seg <- data.table::fread(path, sep = "\t", showProgress = FALSE)
  pick <- function(candidates) {
    hit <- candidates[candidates %in% names(seg)]
    if (length(hit)) hit[[1]] else NA_character_
  }
  chrom_col <- pick(c("chromosome", "chrom", "Chromosome", "chr", "CHROM"))
  start_col <- pick(c("start", "Start", "startpos", "loc.start"))
  end_col <- pick(c("end", "End", "endpos", "loc.end"))
  cn_col <- pick(c("copyNumber", "total_cn", "Total_CN", "tcn", "CN", "copy_number"))
  if (anyNA(c(chrom_col, start_col, end_col, cn_col))) {
    log_step("clonality", "copy-number segment file lacks chromosome/start/end/copy-number columns; using copy-neutral default",
             file = path)
    return(NULL)
  }
  minor_col <- pick(c("minorAlleleCopyNumber", "minor_cn", "nMinor"))
  major_col <- pick(c("majorAlleleCopyNumber", "major_cn", "nMajor"))
  col_or_na <- function(col) if (is.na(col)) rep(NA_real_, nrow(seg)) else to_numeric_safe(seg[[col]])

  out <- data.table::data.table(
    seg_chrom = normalize_cn_chrom(seg[[chrom_col]]),
    seg_start = to_numeric_safe(seg[[start_col]]),
    seg_end = to_numeric_safe(seg[[end_col]]),
    total_cn = to_numeric_safe(seg[[cn_col]]),
    minor_cn = col_or_na(minor_col),
    major_cn = col_or_na(major_col)
  )
  out[!is.na(seg_start) & !is.na(seg_end)]
}

normalize_cn_chrom <- function(x) {
  sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)
}
