# TSV reading (VEP / ANNOVAR / Funcotator / custom). The column line is detected
# after any `#` comments; canonical fields are resolved by alias or an explicit
# YAML column_map, with actionable errors when minimum fields are missing.

# Alias tables for auto-detection of canonical fields from arbitrary TSVs.
TSV_ALIASES <- list(
  CHROM = c("CHROM", "Chromosome", "chr", "chrom", "CHROMOSOME", "Chr"),
  POS = c("POS", "Start_Position", "START", "pos", "Start", "Position"),
  END = c("END", "End_Position", "Stop", "end"),
  REF = c("REF", "Reference_Allele", "ref", "REF_ALLELE", "Ref"),
  ALT = c("ALT", "Tumor_Seq_Allele2", "alt", "ALT_ALLELE", "Alt", "Allele"),
  SAMPLE_ID = c("SAMPLE_ID", "Tumor_Sample_Barcode", "Sample_Barcode", "SAMPLE", "sample", "sample_id", "Sample"),
  GENE = c("GENE", "Hugo_Symbol", "SYMBOL", "Gene", "gene", "Gene.refGene"),
  CONSEQUENCE = c("CONSEQUENCE", "Consequence", "Variant_Classification", "ExonicFunc.refGene", "Func.refGene", "Annotation"),
  HGVSC = c("HGVSC", "HGVSc", "cHGVS"),
  HGVSP = c("HGVSP", "HGVSp", "HGVSp_Short", "pHGVS", "Amino_acids"),
  TRANSCRIPT = c("TRANSCRIPT", "Transcript_ID", "Feature", "transcript"),
  FILTER = c("FILTER", "filter", "Filter"),
  QUAL = c("QUAL", "qual"),
  DP = c("DP", "t_depth", "depth", "TOTAL_DEPTH"),
  AD_REF = c("AD_REF", "t_ref_count", "ref_count", "REF_COUNT"),
  AD_ALT = c("AD_ALT", "t_alt_count", "alt_count", "ALT_COUNT"),
  VAF = c("VAF", "AF", "tumor_vaf", "vaf", "t_vaf")
)

#' Read a TSV/TXT (optionally gzipped) into a data.frame, skipping comments.
#' @keywords internal
read_tsv_table <- function(path, delimiter = NULL) {
  hb <- read_header_block(path)
  skip <- max(0L, (hb$column_line_index %||% 1L) - 1L)
  sep <- delimiter %||% hb$delimiter
  df <- data.table::fread(path, sep = sep, skip = skip, header = TRUE,
                          data.table = FALSE, na.strings = c("", ".", "NA"), quote = "")
  attr(df, "header_block") <- hb
  df
}

#' Resolve canonical fields for a TSV via explicit column_map then aliases.
#' @return named character vector canonical_field -> source_column (NA if unresolved).
#' @keywords internal
resolve_tsv_columns <- function(df, column_map = NULL) {
  cols <- names(df)
  resolved <- setNames(rep(NA_character_, length(TSV_ALIASES)), names(TSV_ALIASES))
  method <- setNames(rep(NA_character_, length(TSV_ALIASES)), names(TSV_ALIASES))
  for (canon in names(TSV_ALIASES)) {
    # explicit map wins
    if (!is.null(column_map) && !is.null(column_map[[tolower(canon)]]) &&
        column_map[[tolower(canon)]] %in% cols) {
      resolved[canon] <- column_map[[tolower(canon)]]; method[canon] <- "column_map"; next
    }
    hit <- TSV_ALIASES[[canon]][TSV_ALIASES[[canon]] %in% cols]
    if (length(hit) >= 1) { resolved[canon] <- hit[[1]]; method[canon] <- "alias" }
  }
  attr(resolved, "method") <- method
  resolved
}
