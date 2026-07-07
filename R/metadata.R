attach_sample_metadata <- function(x, cfg) {
  metadata_path <- cfg_get(cfg, c("input", "sample_metadata"), NULL)
  if (is.null(metadata_path) || is.na(metadata_path) || !file.exists(metadata_path)) {
    return(x)
  }

  meta <- read_variants(metadata_path, cfg_get(cfg, c("input", "metadata_delimiter"), "\t"))
  meta$sample_id <- as.character(coalesce_columns(
    meta,
    c("sample_id", "Sample_Barcode", "Tumor_Sample_Barcode", "Tumor_Sample", "Sample")
  ))
  meta$tumor_type_meta <- as.character(coalesce_columns(
    meta,
    c("tumor_type", "ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "Primary_Site")
  ))
  meta$assay_meta <- as.character(coalesce_columns(meta, c("assay", "Assay", "sequencing_assay")))
  meta$tumor_purity_meta <- normalize_fraction(to_numeric_safe(coalesce_columns(
    meta,
    c("tumor_purity", "purity", "PURITY")
  )))

  # tumor-context fields for COSMIC stratification (original terms preserved)
  meta$tumor_subtype_meta <- as.character(coalesce_columns(
    meta, c("tumor_subtype", "Tumor_Subtype", "ONCOTREE_SUBTYPE", "subtype")))
  meta$primary_site_meta <- as.character(coalesce_columns(
    meta, c("primary_site", "Primary_Site", "PRIMARY_SITE", "site")))
  meta$histology_meta <- as.character(coalesce_columns(
    meta, c("histology", "Histology", "HISTOLOGY", "histological_subtype")))

  keep <- unique(meta[, c("sample_id", "tumor_type_meta", "assay_meta", "tumor_purity_meta",
                          "tumor_subtype_meta", "primary_site_meta", "histology_meta")])
  x <- merge(x, keep, by = "sample_id", all.x = TRUE)
  use_meta_tumor <- !is_missing_value(x$tumor_type_meta)
  x$tumor_type[use_meta_tumor] <- toupper(trimws(x$tumor_type_meta[use_meta_tumor]))
  x$assay <- as.character(coalesce_columns(x, c("assay_meta"), default = cfg_get(cfg, c("hard_filters", "assay"), "WGS")))
  x$tumor_purity <- x$tumor_purity_meta
  x$tumor_subtype <- x$tumor_subtype_meta
  x$primary_site <- x$primary_site_meta
  x$histology <- x$histology_meta
  x$tumor_type_meta <- NULL
  x$assay_meta <- NULL
  x$tumor_purity_meta <- NULL
  x$tumor_subtype_meta <- NULL
  x$primary_site_meta <- NULL
  x$histology_meta <- NULL
  x
}

