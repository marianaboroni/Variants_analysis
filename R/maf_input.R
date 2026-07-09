# MAF reading. Comment lines (`#...`) are preserved as metadata; the column line is
# auto-detected. Standard MAF columns are recognized; extra columns are kept.

MAF_CORE_COLS <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                   "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                   "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")

#' Read a MAF (optionally gzipped) into a data.frame, skipping comment lines.
#' @keywords internal
read_maf_table <- function(path) {
  hb <- read_header_block(path)
  skip <- max(0L, (hb$column_line_index %||% 1L) - 1L)
  df <- data.table::fread(path, sep = "\t", skip = skip, header = TRUE,
                          data.table = FALSE, na.strings = c("", ".", "NA"), quote = "")
  attr(df, "header_block") <- hb
  df
}

#' Validate a MAF table with maftools when possible (non-fatal).
#' @keywords internal
maf_validate_input <- function(df) {
  if (!requireNamespace("maftools", quietly = TRUE)) return(invisible(NULL))
  if (!all(c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification") %in% names(df)))
    return(invisible(NULL))
  tryCatch(invisible(maftools::read.maf(maf = df, verbose = FALSE)),
           error = function(e) { warning("maftools could not validate the input MAF: ",
                                          conditionMessage(e), call. = FALSE); invisible(NULL) })
}
