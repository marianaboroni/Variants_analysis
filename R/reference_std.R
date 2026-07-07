# Generic standardizer for driver-gene / hotspot reference tables. Produces the
# sample/tumor/coord/gene-protein keys used by prioritized_reference_lookup.
# (COSMIC and OncoKB standardizers were removed: COSMIC is prepared via
# prepare_cosmic_db() and OncoKB is a post-hoc API step.)

#' @keywords internal
standardize_reference_core <- function(x) {
  out <- data.frame(
    sample_id = as.character(coalesce_columns(
      x, c("sample_id", "Sample_Barcode", "Tumor_Sample_Barcode", "Tumor_Sample", "Sample"))),
    tumor_type = as.character(coalesce_columns(
      x, c("tumor_type", "ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "Primary_Site"),
      default = "PANCANCER")),
    chrom = as.character(coalesce_columns(x, c("CHROM", "Chromosome", "chr", "chrom"))),
    pos = to_numeric_safe(coalesce_columns(x, c("START", "Start_Position", "POS", "pos"))),
    ref = as.character(coalesce_columns(x, c("REF", "Reference_Allele", "ref"))),
    alt = as.character(coalesce_columns(x, c("ALT", "Tumor_Seq_Allele2", "alt"))),
    gene = as.character(coalesce_columns(x, c("Hugo_Symbol", "SYMBOL", "Gene", "gene"))),
    protein_change = as.character(coalesce_columns(
      x, c("HGVSp_Short", "HGVSp", "Protein_Change", "Amino_acids", "Mutation AA"))),
    stringsAsFactors = FALSE
  )
  out$tumor_type[is_missing_value(out$tumor_type)] <- "PANCANCER"
  out$tumor_type <- toupper(trimws(out$tumor_type))
  add_validation_keys(out)
}
