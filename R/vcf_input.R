# VCF reading with dynamic VEP CSQ parsing and configurable chunked (block)
# reading so whole-genome VCFs are not loaded into memory at once.
# Low-level CSQ/INFO/FORMAT helpers live in io.R and are reused here.

#' Read a VCF body in blocks, applying `fn` to each parsed block.
#'
#' Each block is a data.frame with the fixed VCF columns + per-sample columns,
#' already INFO/CSQ-annotated and sample-expanded (via the io.R helpers). Blocks
#' are read with base connections + `data.table::fread(text=)`, so peak memory is
#' bounded by `chunk_size`, not by the file size.
#'
#' @param path VCF path (.vcf/.vcf.gz).
#' @param chunk_size lines per block (NULL/Inf = single block).
#' @param fn function(block_df, block_index) -> value; values are returned as a list.
#' @return list of `fn` return values (one per block).
#' @keywords internal
read_vcf_blocks <- function(path, chunk_size = NULL, fn = identity) {
  meta <- read_vcf_meta(path)
  csq_fields <- parse_vcf_annotation_format(meta, "CSQ")
  ann_fields <- parse_vcf_annotation_format(meta, "ANN")
  header_line <- meta[startsWith(meta, "#CHROM")]
  if (length(header_line) == 0) stop("VCF has no #CHROM header line: ", path, call. = FALSE)
  col_names <- strsplit(sub("^#", "", header_line[[1]]), "\t", fixed = TRUE)[[1]]

  con <- open_text(path); on.exit(close(con), add = TRUE)
  # skip to first body line
  repeat { z <- readLines(con, n = 1L, warn = FALSE); if (length(z) == 0 || startsWith(z, "#CHROM")) break }
  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size <= 0) chunk_size <- 5e6

  out <- list(); bi <- 0L
  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0) break
    bi <- bi + 1L
    df <- data.table::fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
                            data.table = FALSE, na.strings = c("", ".", "NA"), quote = "")
    names(df) <- col_names[seq_len(ncol(df))]
    if ("#CHROM" %in% names(df)) names(df)[names(df) == "#CHROM"] <- "CHROM"
    if (!"CHROM" %in% names(df) && "CHROM" %in% col_names) names(df)[1] <- "CHROM"
    df <- expand_vcf_samples(df, path)
    df <- add_info_annotations(df, csq_fields, ann_fields)
    out[[bi]] <- fn(df, bi)
    if (length(lines) < chunk_size) break
  }
  out
}

#' Read a full VCF into one data.frame (block-wise under the hood).
#' @keywords internal
read_vcf_canonical <- function(path, chunk_size = NULL) {
  blocks <- read_vcf_blocks(path, chunk_size, fn = function(df, i) df)
  if (length(blocks) == 1) blocks[[1]] else data.table::rbindlist(blocks, fill = TRUE, use.names = TRUE) |>
    as.data.frame()
}
