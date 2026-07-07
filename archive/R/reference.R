standardize_oncokb_reference <- function(path, output_path = NULL, delimiter = "\t") {
  x <- read_variants(path, delimiter)
  out <- standardize_reference_core(x)
  out$oncokb_oncogenic <- as.character(coalesce_columns(x, c("ONCOGENIC", "Oncogenic")))
  out$oncokb_mutation_effect <- as.character(coalesce_columns(x, c("MUTATION_EFFECT", "Mutation_Effect")))
  out$oncokb_highest_level <- as.character(coalesce_columns(
    x,
    c("HIGHEST_LEVEL", "Highest_Level", "LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4")
  ))
  out$source <- "OncoKB"
  out <- unique(out)
  if (!is.null(output_path)) write_tsv(out, output_path)
  out
}

standardize_cosmic_reference <- function(path, output_path = NULL, delimiter = "\t") {
  x <- read_variants(path, delimiter)
  out <- standardize_reference_core(x)
  out$cosmic_id <- as.character(coalesce_columns(
    x,
    c("COSMIC_ID", "COSMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "Existing_variation")
  ))
  out$cosmic_count <- to_numeric_safe(coalesce_columns(x, c("COSMIC_COUNT", "CNT", "Count", "Mutation somatic status count")))
  out$cosmic_primary_site <- as.character(coalesce_columns(x, c("Primary site", "PRIMARY_SITE", "primary_site")))
  out$source <- "COSMIC"
  out <- unique(out)
  if (!is.null(output_path)) write_tsv(out, output_path)
  out
}

standardize_reference_core <- function(x) {
  out <- data.frame(
    sample_id = as.character(coalesce_columns(
      x,
      c("sample_id", "Sample_Barcode", "Tumor_Sample_Barcode", "Tumor_Sample", "Sample")
    )),
    tumor_type = as.character(coalesce_columns(
      x,
      c("tumor_type", "ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "Primary_Site"),
      default = "PANCANCER"
    )),
    chrom = as.character(coalesce_columns(x, c("CHROM", "Chromosome", "chr", "chrom"))),
    pos = to_numeric_safe(coalesce_columns(x, c("START", "Start_Position", "POS", "pos"))),
    ref = as.character(coalesce_columns(x, c("REF", "Reference_Allele", "ref"))),
    alt = as.character(coalesce_columns(x, c("ALT", "Tumor_Seq_Allele2", "alt"))),
    gene = as.character(coalesce_columns(x, c("Hugo_Symbol", "SYMBOL", "Gene", "gene"))),
    protein_change = as.character(coalesce_columns(
      x,
      c("HGVSp_Short", "HGVSp", "Protein_Change", "Amino_acids", "Mutation AA")
    )),
    stringsAsFactors = FALSE
  )
  out$tumor_type[is_missing_value(out$tumor_type)] <- "PANCANCER"
  out$tumor_type <- toupper(trimws(out$tumor_type))
  out <- add_validation_keys(out)
  out
}
