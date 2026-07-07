# Reference-key machinery used for DRIVER-GENE and CANCER-HOTSPOT matching
# (gene / tumor-type / protein keys). COSMIC and OncoKB matching were removed
# from this file: COSMIC now uses the canonical genomic key (R/cosmic_match.R)
# and OncoKB is a post-hoc step (R/oncokb.R). See docs/REFACTOR_AUDIT.md.

standardize_validation_table <- function(x) {
  x$sample_id <- as.character(coalesce_columns(
    x,
    c("Tumor_Sample_Barcode", "Sample_Barcode", "Tumor_Sample", "Sample", "sample", "sample_id")
  ))
  x$chrom <- as.character(coalesce_columns(x, c("CHROM", "Chromosome", "chr", "chrom")))
  x$pos <- to_numeric_safe(coalesce_columns(x, c("START", "Start_Position", "POS", "pos")))
  x$ref <- as.character(coalesce_columns(x, c("REF", "Reference_Allele", "ref")))
  x$alt <- as.character(coalesce_columns(x, c("ALT", "Tumor_Seq_Allele2", "alt")))
  x$gene <- as.character(coalesce_columns(x, c("Hugo_Symbol", "SYMBOL", "Gene", "gene")))
  x$protein_change <- as.character(coalesce_columns(
    x,
    c("HGVSp_Short", "HGVSp", "Protein_Change", "Amino_acids", "Mutation AA")
  ))
  x$tumor_type <- as.character(coalesce_columns(
    x,
    c("tumor_type", "ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "Primary_Site"),
    default = "PANCANCER"
  ))
  x$tumor_type[is_missing_value(x$tumor_type)] <- "PANCANCER"
  x$tumor_type <- toupper(trimws(x$tumor_type))
  x <- add_validation_keys(x)
  x$validation_key <- make_validation_key(x)
  x
}

make_validation_key <- function(x) {
  has_coord <- !is.na(x$sample_id) & !is.na(x$chrom) & !is.na(x$pos) &
    !is.na(x$ref) & !is.na(x$alt)
  ifelse(has_coord, x$sample_coord_validation_key, x$sample_gene_protein_validation_key)
}

add_validation_keys <- function(x) {
  x$coord_validation_key <- make_coord_validation_key(x)
  x$gene_protein_validation_key <- make_gene_protein_validation_key(x)
  x$sample_coord_validation_key <- make_sample_coord_validation_key(x)
  x$sample_gene_protein_validation_key <- make_sample_gene_protein_validation_key(x)
  x$tumor_coord_validation_key <- make_tumor_coord_validation_key(x)
  x$tumor_gene_protein_validation_key <- make_tumor_gene_protein_validation_key(x)
  x
}

make_coord_validation_key <- function(x) {
  paste(x$chrom, x$pos, x$ref, x$alt, sep = "|")
}

make_gene_protein_validation_key <- function(x) {
  paste(x$gene, x$protein_change, sep = "|")
}

make_sample_coord_validation_key <- function(x) {
  paste(x$sample_id, x$coord_validation_key, sep = "|")
}

make_sample_gene_protein_validation_key <- function(x) {
  paste(x$sample_id, x$gene_protein_validation_key, sep = "|")
}

make_tumor_coord_validation_key <- function(x) {
  paste(x$tumor_type, x$coord_validation_key, sep = "|")
}

make_tumor_gene_protein_validation_key <- function(x) {
  paste(x$tumor_type, x$gene_protein_validation_key, sep = "|")
}

prioritized_reference_lookup <- function(x, ref, value_cols, prefixes) {
  out <- data.frame(match_scope = rep(NA_character_, nrow(x)), stringsAsFactors = FALSE)
  for (col in value_cols) out[[col]] <- rep(NA_character_, nrow(x))

  for (prefix in prefixes) {
    key_col <- paste0(prefix, "_validation_key")
    if (!(key_col %in% names(x)) || !(key_col %in% names(ref))) next
    ref_keep <- collapse_reference_by_key(ref, key_col, value_cols)
    idx <- match(x[[key_col]], ref_keep[[key_col]])
    still_empty <- is.na(out$match_scope) & !is.na(idx)
    if (!any(still_empty)) next
    out$match_scope[still_empty] <- prefix
    for (col in value_cols) {
      out[[col]][still_empty] <- as.character(ref_keep[[col]][idx[still_empty]])
    }
  }
  out
}

collapse_reference_by_key <- function(ref, key_col, value_cols) {
  keep <- ref[!is_missing_value(ref[[key_col]]), c(key_col, value_cols), drop = FALSE]
  if (nrow(keep) == 0) return(keep)
  agg <- stats::aggregate(
    keep[value_cols],
    by = list(lookup_key = keep[[key_col]]),
    FUN = function(z) {
      z <- unique(stats::na.omit(as.character(z)))
      if (length(z) == 0) NA_character_ else paste(z[1:min(length(z), 5)], collapse = ";")
    }
  )
  names(agg)[names(agg) == "lookup_key"] <- key_col
  agg
}
