add_validation_layers <- function(x, cfg) {
  x <- add_validation_keys(x)
  x$validation_key <- make_validation_key(x)
  x$oncokb_match <- FALSE
  x$oncokb_oncogenic <- NA_character_
  x$oncokb_highest_level <- NA_character_
  x$cosmic_match <- FALSE
  x$cosmic_id <- NA_character_
  x$cosmic_count <- NA_real_

  if (!is.null(cfg$validation$oncokb) && !is.na(cfg$validation$oncokb)) {
    okb <- read_variants(cfg$validation$oncokb, "\t")
    x <- add_oncokb(x, okb)
  }
  if (!is.null(cfg$validation$cosmic) && !is.na(cfg$validation$cosmic)) {
    cosmic <- read_variants(cfg$validation$cosmic, "\t")
    x <- add_cosmic(x, cosmic)
  }

  x$validation_support_score <- validation_support_score(x)
  x
}

add_oncokb <- function(x, okb) {
  okb <- standardize_validation_table(okb)
  okb$oncokb_oncogenic_src <- as.character(coalesce_columns(okb, c("ONCOGENIC", "Oncogenic")))
  okb$oncokb_highest_level_src <- as.character(coalesce_columns(
    okb,
    c("HIGHEST_LEVEL", "Highest_Level", "LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4")
  ))

  match <- prioritized_reference_lookup(
    x,
    okb,
    value_cols = c("oncokb_oncogenic_src", "oncokb_highest_level_src"),
    prefixes = c("sample_coord", "sample_gene_protein", "tumor_coord", "tumor_gene_protein", "coord", "gene_protein")
  )

  hit <- !is.na(match$oncokb_oncogenic_src) | !is.na(match$oncokb_highest_level_src)
  x$oncokb_match <- hit
  x$oncokb_oncogenic <- match$oncokb_oncogenic_src
  x$oncokb_highest_level <- match$oncokb_highest_level_src
  x$oncokb_match_scope <- match$match_scope
  x
}

add_cosmic <- function(x, cosmic) {
  cosmic <- standardize_validation_table(cosmic)
  cosmic$cosmic_id_src <- as.character(coalesce_columns(
    cosmic,
    c("COSMIC_ID", "COSMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "Existing_variation")
  ))
  cosmic$cosmic_count_src <- to_numeric_safe(coalesce_columns(cosmic, c("COSMIC_COUNT", "CNT", "Count")))

  match <- prioritized_reference_lookup(
    x,
    cosmic,
    value_cols = c("cosmic_id_src", "cosmic_count_src"),
    prefixes = c("sample_coord", "sample_gene_protein", "tumor_coord", "tumor_gene_protein", "coord", "gene_protein")
  )
  x$cosmic_match <- !is.na(match$cosmic_id_src) | !is.na(match$cosmic_count_src)
  x$cosmic_id <- match$cosmic_id_src
  x$cosmic_count <- to_numeric_safe(match$cosmic_count_src)
  x$cosmic_match_scope <- match$match_scope
  x
}

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
  agg <- aggregate(
    keep[value_cols],
    by = list(lookup_key = keep[[key_col]]),
    FUN = function(z) {
      z <- unique(na.omit(as.character(z)))
      if (length(z) == 0) NA_character_ else paste(z[1:min(length(z), 5)], collapse = ";")
    }
  )
  names(agg)[names(agg) == "lookup_key"] <- key_col
  agg
}

validation_support_score <- function(x) {
  okb <- rep(0, nrow(x))
  onc <- toupper(as.character(x$oncokb_oncogenic))
  okb[grepl("ONCOGENIC|LIKELY", onc)] <- 1
  okb[x$oncokb_match & okb == 0] <- 0.5

  cosmic <- rep(0, nrow(x))
  cosmic[x$cosmic_match] <- 0.4
  cosmic[!is.na(x$cosmic_count) & x$cosmic_count >= 3] <- 0.7
  cosmic[!is.na(x$cosmic_count) & x$cosmic_count >= 10] <- 0.9
  pmax(okb, cosmic)
}
