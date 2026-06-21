read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with install.packages('yaml').")
  }
  yaml::read_yaml(path)
}

read_variants <- function(path, delimiter = "\t") {
  if (!file.exists(path)) {
    stop("Input variant file not found: ", path)
  }
  data.table::fread(path, sep = delimiter, data.table = FALSE, na.strings = c("", ".", "NA"))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

coalesce_columns <- function(x, candidates, default = NA) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(rep(default, nrow(x)))
  out <- x[[hit[[1]]]]
  if (length(hit) > 1) {
    for (nm in hit[-1]) {
      out[is.na(out) | out == ""] <- x[[nm]][is.na(out) | out == ""]
    }
  }
  out
}

to_numeric_safe <- function(x) {
  if (is.null(x)) return(NA_real_)
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = FALSE)))
}

is_missing_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

cfg_get <- function(cfg, path, default = NULL) {
  if (is.null(cfg)) return(default)
  cur <- cfg
  for (key in path) {
    if (is.null(cur[[key]])) return(default)
    cur <- cur[[key]]
  }
  if (length(cur) == 1 && is.na(cur)) default else cur
}

default_tumor_type <- function(cfg = NULL) {
  value <- cfg_get(cfg, c("cancer", "default_tumor_type"), NULL)
  if (is.null(value)) value <- cfg_get(cfg, c("cohort", "cancer_type"), "PANCANCER")
  as.character(value)
}

standardize_variant_table <- function(x, cfg = NULL) {
  x$sample_id <- as.character(coalesce_columns(
    x,
    c("Tumor_Sample_Barcode", "Sample_Barcode", "Tumor_Sample", "Sample", "sample", "sample_id")
  ))
  x$chrom <- as.character(coalesce_columns(x, c("CHROM", "Chromosome", "chr", "chrom")))
  x$pos <- to_numeric_safe(coalesce_columns(x, c("START", "Start_Position", "POS", "pos")))
  x$ref <- as.character(coalesce_columns(x, c("REF", "Reference_Allele", "ref")))
  x$alt <- as.character(coalesce_columns(x, c("ALT", "Tumor_Seq_Allele2", "alt")))
  x$gene <- as.character(coalesce_columns(x, c("Hugo_Symbol", "SYMBOL", "Gene", "gene")))
  x$consequence <- as.character(coalesce_columns(
    x,
    c("Consequence", "Variant_Classification", "EFFECT", "Annotation")
  ))
  x$protein_change <- as.character(coalesce_columns(
    x,
    c("HGVSp_Short", "HGVSp", "Protein_Change", "Amino_acids")
  ))
  x$tumor_type <- as.character(coalesce_columns(
    x,
    c("ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "tumor_type", "Primary_Site"),
    default = default_tumor_type(cfg)
  ))
  x$filter_status <- as.character(coalesce_columns(x, c("FILTER", "filter", "Filter"), default = "NA"))

  if (any(is.na(x$sample_id))) {
    warning("Some variants have missing sample IDs.")
  }
  required <- c("chrom", "pos", "ref", "alt")
  missing_core <- vapply(required, function(nm) all(is.na(x[[nm]]) | x[[nm]] == ""), logical(1))
  if (any(missing_core)) {
    stop("Missing core coordinate columns after standardization: ",
         paste(required[missing_core], collapse = ", "))
  }
  x$variant_id <- paste(x$chrom, x$pos, x$ref, x$alt, sep = ":")
  x$variant_key <- x$variant_id
  x$sample_variant_key <- ifelse(
    !is.na(x$sample_id) & x$sample_id != "",
    paste(x$sample_id, x$variant_id, sep = "|")
    , NA_character_
  )
  x$locus_id <- paste(x$chrom, x$pos, sep = ":")
  x$tumor_type[is_missing_value(x$tumor_type)] <- default_tumor_type(cfg)
  x$tumor_type <- toupper(trimws(x$tumor_type))
  x
}
