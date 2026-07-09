# Single ingestion entry point. Detects the format (VCF/MAF/TSV - not by
# extension alone), interprets `#` headers, maps to the canonical internal schema
# while preserving all original columns, and emits an ingestion (schema) report.
# Build detection/validation lives here too.

CANONICAL_FIELDS <- c("CHROM", "POS", "END", "REF", "ALT", "SAMPLE_ID", "FILTER", "QUAL",
  "DP", "AD_REF", "AD_ALT", "VAF", "GENE", "TRANSCRIPT", "CONSEQUENCE", "HGVSC",
  "HGVSP", "EXON", "INTRON", "GENOME_BUILD")
CANONICAL_REQUIRED <- c("CHROM", "POS", "REF", "ALT")

CHR1_LENGTH <- c(GRCh37 = 249250621L, GRCh38 = 248956422L)

#' Read a tumor-only variant input (VCF / MAF / TSV) into the canonical schema.
#'
#' @param path input file (.vcf/.vcf.gz/.maf/.maf.gz/.tsv/.tsv.gz/.txt/.txt.gz).
#' @param format "auto" (default) or one of "vcf"/"maf"/"tsv".
#' @param sample_metadata optional path/data.frame with per-sample metadata.
#' @param column_map optional named list mapping canonical (lowercase) fields to
#'   source columns (for custom TSVs).
#' @param chunk_size optional block size for VCF reading (bounded memory).
#' @return a list with `variants` (canonical data.frame, originals preserved),
#'   `format`, `header_metadata`, `schema_report`, `column_mapping`.
#' @export
#' @examples
#' \dontrun{
#' res <- read_variant_input("data/sample.vep.vcf.gz")
#' head(res$variants); res$schema_report
#' }
read_variant_input <- function(path, format = "auto", sample_metadata = NULL,
                               column_map = NULL, chunk_size = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  fmt <- detect_input_format(path, format)
  hb <- read_header_block(path)

  if (fmt == "vcf") {
    header_meta <- parse_vcf_header_meta(hb$header_lines)
    df <- read_vcf_canonical(path, chunk_size)
    mapping <- vcf_column_mapping(df)
  } else if (fmt == "maf") {
    header_meta <- parse_comment_meta(hb$header_lines)
    df <- read_maf_table(path); maf_validate_input(df)
    mapping <- maf_column_mapping(df)
  } else {
    header_meta <- parse_comment_meta(hb$header_lines)
    df <- read_tsv_table(path)
    mapping <- resolve_tsv_columns(df, column_map)
  }

  built <- build_canonical(df, mapping, fmt, header_meta)
  variants <- built$variants
  report <- build_schema_report(variants, built$mapping, built$method)

  missing_req <- CANONICAL_REQUIRED[!CANONICAL_REQUIRED %in% names(variants) |
    vapply(CANONICAL_REQUIRED, function(f) f %in% names(variants) && all(is.na(variants[[f]])), logical(1))]
  if (length(missing_req)) {
    stop(sprintf(paste0("Input is missing required canonical field(s): %s.\n",
      "Detected format: %s. Provide input.column_map to map them, e.g.:\n",
      "  input:\n    column_map:\n      chrom: <your_column>\n      pos: <your_column>\n",
      "      ref: <your_column>\n      alt: <your_column>"),
      paste(missing_req, collapse = ", "), fmt), call. = FALSE)
  }

  variants$INPUT_FORMAT <- fmt
  attr(variants, "input_format") <- fmt
  attr(variants, "header_metadata") <- header_meta
  attr(variants, "schema_report") <- report
  attr(variants, "column_mapping") <- built$mapping
  list(variants = variants, format = fmt, header_metadata = header_meta,
       schema_report = report, column_mapping = built$mapping)
}

#' @keywords internal
vcf_column_mapping <- function(df) {
  m <- c(CHROM = "CHROM", POS = "POS", REF = "REF", ALT = "ALT", FILTER = "FILTER",
    QUAL = "QUAL", SAMPLE_ID = "Tumor_Sample_Barcode", GENE = "SYMBOL",
    CONSEQUENCE = "Consequence", HGVSP = "HGVSp", HGVSC = "HGVSc", TRANSCRIPT = "Feature",
    EXON = "EXON", INTRON = "INTRON", DP = "DP", VAF = "AF")
  m[m %in% names(df)]
}
#' @keywords internal
maf_column_mapping <- function(df) {
  m <- c(CHROM = "Chromosome", POS = "Start_Position", END = "End_Position",
    REF = "Reference_Allele", ALT = "Tumor_Seq_Allele2", SAMPLE_ID = "Tumor_Sample_Barcode",
    GENE = "Hugo_Symbol", CONSEQUENCE = "Variant_Classification", HGVSP = "HGVSp_Short",
    HGVSC = "HGVSc", TRANSCRIPT = "Transcript_ID", FILTER = "FILTER",
    DP = "t_depth", AD_REF = "t_ref_count", AD_ALT = "t_alt_count")
  m[m %in% names(df)]
}

#' Build canonical UPPERCASE columns from a source df + mapping, preserving all
#' original columns and useful downstream aliases.
#' @keywords internal
build_canonical <- function(df, mapping, fmt, header_meta) {
  method <- attr(mapping, "method")
  for (canon in names(mapping)) {
    src <- mapping[[canon]]
    if (!is.na(src) && src %in% names(df) && !(canon %in% names(df))) df[[canon]] <- df[[src]]
  }
  if ("CHROM" %in% names(df)) df$CHROM <- as.character(df$CHROM)  # keep chrom names as text (X/Y/MT)
  # derived VAF from allele depths if absent
  if (!"VAF" %in% names(df) && all(c("AD_ALT", "AD_REF") %in% names(df))) {
    a <- to_numeric_safe(df$AD_ALT); r <- to_numeric_safe(df$AD_REF)
    df$VAF <- ifelse((a + r) > 0, a / (a + r), NA_real_)
  }
  if (!"END" %in% names(df) && all(c("POS", "REF") %in% names(df)))
    df$END <- to_numeric_safe(df$POS) + pmax(nchar(as.character(df$REF)) - 1L, 0L)
  # downstream-compatibility aliases (standardize_variant_table/features/predictors)
  alias <- function(dst, src) if (!(dst %in% names(df)) && src %in% names(df)) df[[dst]] <<- df[[src]]
  alias("Tumor_Sample_Barcode", "SAMPLE_ID"); alias("SYMBOL", "GENE")
  alias("Consequence", "CONSEQUENCE"); alias("HGVSp", "HGVSP")
  list(variants = df, mapping = mapping, method = method)
}

#' Build the ingestion / schema report (one row per canonical field).
#' @keywords internal
build_schema_report <- function(variants, mapping, method) {
  rows <- lapply(CANONICAL_FIELDS, function(f) {
    src <- if (f %in% names(mapping)) mapping[[f]] else NA_character_
    det <- if (!is.null(method) && f %in% names(method)) method[[f]] else if (!is.na(src)) "alias" else "none"
    present <- f %in% names(variants) && !all(is.na(variants[[f]]))
    vals <- if (f %in% names(variants)) variants[[f]] else rep(NA, nrow(variants))
    nmiss <- sum(is.na(vals) | vals == "")
    ex <- utils::head(unique(stats::na.omit(as.character(vals[vals != "" & !is.na(vals)]))), 3)
    data.frame(
      canonical_field = f,
      detected_source_column = src %||% NA_character_,
      detection_method = det,
      required = f %in% CANONICAL_REQUIRED,
      present = present,
      data_type = if (f %in% names(variants)) class(variants[[f]])[1] else NA_character_,
      n_missing = nmiss,
      fraction_missing = round(nmiss / max(nrow(variants), 1), 4),
      example_values = paste(ex, collapse = " | "),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- genome build (unchanged policy) ----------------------------------------

#' @keywords internal
detect_genome_build <- function(meta) {
  if (length(meta) == 0) return(NA_character_)
  ref_line <- meta[grepl("^##reference", meta, ignore.case = TRUE)]
  if (length(ref_line) > 0) {
    rl <- tolower(paste(ref_line, collapse = " "))
    if (grepl("grch38|hg38|hs38", rl) && !grepl("grch37|hg19", rl)) return("GRCh38")
    if (grepl("grch37|hg19|hs37", rl) && !grepl("grch38|hg38", rl)) return("GRCh37")
  }
  contig1 <- meta[grepl("^##contig", meta) & grepl("ID=(chr)?1[,>]", meta, ignore.case = TRUE)]
  if (length(contig1) > 0) {
    len <- suppressWarnings(as.integer(sub(".*length=([0-9]+).*", "\\1", contig1[[1]])))
    if (!is.na(len)) { hit <- names(CHR1_LENGTH)[match(len, CHR1_LENGTH)]; if (!is.na(hit)) return(hit) }
  }
  NA_character_
}

#' @keywords internal
resolve_genome_build <- function(cfg, vcf_path = NULL) {
  cfg_build <- cfg_get(cfg, c("input", "genome_build"), NULL)
  if (!is.null(cfg_build) && is.na(cfg_build)) cfg_build <- NULL
  header_build <- NA_character_
  if (!is.null(vcf_path) && file.exists(vcf_path) && grepl("[.]vcf([.]gz)?$", vcf_path, ignore.case = TRUE))
    header_build <- detect_genome_build(read_vcf_meta(vcf_path))
  if (!is.null(cfg_build) && !is.na(header_build) && cfg_build != header_build)
    stop(sprintf("Genome-build conflict: config=%s but VCF header=%s. Resolve before continuing.",
                 cfg_build, header_build), call. = FALSE)
  build <- cfg_build %||% (if (!is.na(header_build)) header_build else NULL)
  if (is.null(build)) stop(paste0("Could not determine GRCh37 vs GRCh38.\n",
    "Set input.genome_build in the configuration file."), call. = FALSE)
  if (!build %in% names(CHR1_LENGTH))
    stop(sprintf("Unsupported genome build '%s' (expected GRCh37 or GRCh38).", build), call. = FALSE)
  log_step("build", "genome build resolved", build = build,
           source = if (!is.null(cfg_build)) "config" else "vcf_header")
  build
}
