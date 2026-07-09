# Header/metadata parsing for all input formats. Lines beginning with `#` are
# read and interpreted before the tabular body; nothing is discarded silently.

#' Open a possibly-gzipped text connection.
#' @keywords internal
open_text <- function(path) {
  if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
}

#' Read the leading comment/header block and locate the first tabular line.
#'
#' @return list(header_lines, column_line, column_line_index, delimiter).
#' @keywords internal
read_header_block <- function(path, max_lines = 20000) {
  con <- open_text(path); on.exit(close(con), add = TRUE)
  header <- character(); column_line <- NA_character_; idx <- NA_integer_; vcf_chrom <- FALSE
  i <- 0L
  repeat {
    z <- readLines(con, n = 1L, warn = FALSE)
    if (length(z) == 0) break
    i <- i + 1L
    if (startsWith(z, "##")) { header <- c(header, z); next }
    if (startsWith(z, "#CHROM")) { column_line <- sub("^#", "", z); idx <- i; vcf_chrom <- TRUE; break }  # VCF
    if (startsWith(z, "#")) { header <- c(header, z); next }                            # MAF/TSV comment
    # first non-comment line is the column header (MAF/TSV)
    column_line <- z; idx <- i; break
  }
  delim <- if (grepl("\t", column_line %||% "")) "\t" else if (grepl(",", column_line %||% "")) "," else "\t"
  list(header_lines = header, column_line = column_line, column_line_index = idx,
       delimiter = delim, vcf_chrom = vcf_chrom)
}

#' Parse VCF `##` metadata into a named structure (reference, VEP, contigs, ...).
#' @keywords internal
parse_vcf_header_meta <- function(header_lines) {
  get1 <- function(prefix) {
    hit <- header_lines[startsWith(header_lines, prefix)]
    if (length(hit) == 0) NA_character_ else sub(prefix, "", hit[[1]])
  }
  list(
    fileformat = get1("##fileformat="),
    reference = get1("##reference="),
    vep = { v <- header_lines[grepl("^##VEP", header_lines)]; if (length(v)) v[[1]] else NA_character_ },
    snpeff = { v <- header_lines[grepl("^##SnpEff", header_lines)]; if (length(v)) v[[1]] else NA_character_ },
    n_contigs = sum(grepl("^##contig", header_lines)),
    annotation_tool = if (any(grepl("^##VEP", header_lines))) "VEP"
                      else if (any(grepl("^##SnpEff", header_lines))) "SnpEff" else NA_character_,
    raw = header_lines)
}

#' Parse `#key value` / `#key=value` comment metadata (MAF/TSV).
#' @keywords internal
parse_comment_meta <- function(header_lines) {
  comments <- header_lines[startsWith(header_lines, "#") & !startsWith(header_lines, "##")]
  meta <- list()
  for (c0 in comments) {
    body <- sub("^#+\\s*", "", c0)
    kv <- regmatches(body, regexec("^([A-Za-z0-9_.]+)\\s*[:=]?\\s*(.*)$", body))[[1]]
    if (length(kv) == 3 && nzchar(kv[2])) meta[[kv[2]]] <- kv[3]
  }
  meta$raw <- comments
  meta
}

#' Detect the input format from extension + content (never extension alone).
#' @keywords internal
detect_input_format <- function(path, declared = "auto") {
  if (!is.null(declared) && !is.na(declared) && declared != "auto") return(tolower(declared))
  hb <- read_header_block(path)
  # a true VCF has ##fileformat=VCF or a real #CHROM header line (not merely a
  # TSV whose first column happens to be named CHROM)
  if (any(grepl("^##fileformat=VCF", hb$header_lines)) || isTRUE(hb$vcf_chrom))
    return("vcf")
  cols <- strsplit(hb$column_line %||% "", hb$delimiter, fixed = TRUE)[[1]]
  maf_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "Tumor_Sample_Barcode")
  if (sum(maf_cols %in% cols) >= 3) return("maf")
  # extension fallback for ambiguous
  if (grepl("[.]maf([.]gz)?$", path, ignore.case = TRUE)) return("maf")
  if (grepl("[.]vcf([.]gz)?$", path, ignore.case = TRUE)) return("vcf")
  "tsv"
}
