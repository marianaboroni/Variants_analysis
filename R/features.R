add_basic_features <- function(x, cfg) {
  x$dp <- to_numeric_safe(coalesce_columns(x, c("DP", "t_depth", "depth", "Tumor_Depth")))
  x$alt_count <- to_numeric_safe(coalesce_columns(x, c("alt_count", "t_alt_count", "AD_ALT", "AO")))
  x$ref_count <- to_numeric_safe(coalesce_columns(x, c("ref_count", "t_ref_count", "AD_REF", "RO")))

  if (all(is.na(x$alt_count)) && "AD" %in% names(x)) {
    parsed <- parse_ad(x$AD)
    x$ref_count <- parsed$ref
    x$alt_count <- parsed$alt
  }
  if (all(is.na(x$dp)) && !all(is.na(x$alt_count)) && !all(is.na(x$ref_count))) {
    x$dp <- x$alt_count + x$ref_count
  }

  x$vaf <- to_numeric_safe(coalesce_columns(x, c("AF", "VAF", "tumor_f", "Allele_Fraction")))
  missing_vaf <- is.na(x$vaf) & !is.na(x$alt_count) & !is.na(x$dp) & x$dp > 0
  x$vaf[missing_vaf] <- x$alt_count[missing_vaf] / x$dp[missing_vaf]
  x$vaf <- ifelse(x$vaf > 1 & x$vaf <= 100, x$vaf / 100, x$vaf)

  x$tlod <- to_numeric_safe(coalesce_columns(x, c("TLOD", "Tumor_LOD")))
  x$nlod <- to_numeric_safe(coalesce_columns(x, c("NLOD", "Normal_LOD")))
  x$as_filter_status <- as.character(coalesce_columns(x, c("AS_FilterStatus", "AS_Filter_Status", "AS_Filter"), default = NA))
  x$mapping_quality <- to_numeric_safe(coalesce_columns(x, c("MQ", "mapping_quality", "MAPQ")))
  x$base_quality <- to_numeric_safe(coalesce_columns(x, c("MBQ", "BaseQRankSum", "BQ")))
  x$read_position_bias <- to_numeric_safe(coalesce_columns(x, c("ReadPosRankSum", "ReadPos", "read_pos")))
  x$orientation_bias <- to_numeric_safe(coalesce_columns(x, c("F1R2", "F2R1", "orientation_bias", "OB")))
  x$strand_artifact <- to_numeric_safe(coalesce_columns(x, c("SB", "StrandBias", "strand_bias")))
  x$clustered_events <- as_logical_flag(coalesce_columns(x, c("clustered_events", "clustered_event", "clustered")))
  x$weak_evidence <- as_logical_flag(coalesce_columns(x, c("weak_evidence", "WeakEvidence", "weakEvidence")))
  x$germline_risk <- to_numeric_safe(coalesce_columns(x, c("germline_risk", "GermlineRisk", "germline_probability", "germline_prob")))
  x$contamination <- to_numeric_safe(coalesce_columns(x, c("contamination", "CONTAMINATION", "contamination_fraction", "contam")))
  x$mutect_filter <- as.character(coalesce_columns(x, c("FILTER", "filter", "Filter"), default = ""))

  if ("Tumor_LOD" %in% names(x)) {
    x$tlod <- ifelse(is.na(x$tlod) & !is.na(x$Tumor_LOD), to_numeric_safe(x$Tumor_LOD), x$tlod)
  }
  x$mbq <- extract_alt_metric(coalesce_columns(x, c("MBQ", "BaseQRankSum", "BQ")))
  x$mmq <- extract_alt_metric(coalesce_columns(x, c("MMQ", "MQ", "mapping_quality")))
  x$strandq <- to_numeric_safe(coalesce_columns(x, c("STRANDQ", "STRQ", "SBQ", "SB")))
  x$popaf_mutect <- to_numeric_safe(coalesce_columns(x, c("POPAF", "popaf", "Mutect_PopAF")))
  x$pon_flag <- as_logical_flag(coalesce_columns(x, c("PON", "panel_of_normals", "pon")))

  x$abraom_af <- max_numeric_columns(x, c("ABraOM_AF", "AbraOM_AF", "ABraOM_MAF", "abraom_af"))
  x$sabe_af <- max_numeric_columns(x, c("SABE_AF", "Sabe_AF", "sabe_af"))
  x$gnomad_af <- max_numeric_columns(x, c(
    "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMR_AF", "gnomADg_NFE_AF",
    "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_NFE_AF",
    "gnomad_af", "gnomad_genome_af", "gnomad_exome_af"
  ))
  x$dbsnp_present <- !is_missing_value(coalesce_columns(x, c("dbSNP", "Existing_variation", "RSID", "rsid", "dbsnp")))
  x$clinvar_significance <- as.character(coalesce_columns(x, c("CLINVAR_SIG", "ClinVar_Significance", "clinvar_significance", "CLINVAR", "clinvar"), default = NA))

  x$variant_type <- infer_variant_type(x$ref, x$alt)
  x$is_indel <- x$variant_type %in% c("INS", "DEL", "COMPLEX_INDEL")
  x$impact_rank <- consequence_rank(x$consequence)
  x$max_pop_af <- max_population_af(x, cfg$population_filters$use_population_columns)
  x$distance_vaf_05 <- abs(x$vaf - 0.5)
  x$distance_vaf_10 <- abs(x$vaf - 1.0)
  x$technical_pass_floor <- technical_floor_pass(x, cfg)
  x$assay_type <- toupper(as.character(cfg_get(cfg, c("hard_filters", "assay"), "WGS")))
  x
}

parse_ad <- function(ad) {
  vals <- strsplit(as.character(ad), ",")
  ref <- alt <- rep(NA_real_, length(vals))
  for (i in seq_along(vals)) {
    z <- suppressWarnings(as.numeric(vals[[i]]))
    if (length(z) >= 2) {
      ref[[i]] <- z[[1]]
      alt[[i]] <- z[[2]]
    }
  }
  data.frame(ref = ref, alt = alt)
}

extract_alt_metric <- function(x) {
  vals <- strsplit(as.character(x), ",")
  out <- rep(NA_real_, length(vals))
  for (i in seq_along(vals)) {
    z <- suppressWarnings(as.numeric(vals[[i]]))
    if (length(z) > 0) out[[i]] <- tail(z, 1)
  }
  out
}

as_logical_flag <- function(x) {
  z <- toupper(as.character(x))
  z %in% c("TRUE", "T", "YES", "Y", "1", "PON")
}

infer_variant_type <- function(ref, alt) {
  ref_n <- nchar(as.character(ref))
  alt_n <- nchar(as.character(alt))
  ifelse(ref == "-", "INS",
         ifelse(alt == "-", "DEL",
                ifelse(ref_n == 1 & alt_n == 1, "SNV",
                       ifelse(ref_n < alt_n, "INS",
                              ifelse(ref_n > alt_n, "DEL", "COMPLEX_INDEL")))))
}

consequence_rank <- function(consequence) {
  z <- tolower(as.character(consequence))
  high <- grepl("frameshift|stop_gained|splice_acceptor|splice_donor|start_lost|stop_lost", z)
  moderate <- grepl("missense|inframe|protein_altering|splice_region", z)
  low <- grepl("synonymous|stop_retained|start_retained", z)
  ifelse(high, 3, ifelse(moderate, 2, ifelse(low, 1, 0)))
}

max_population_af <- function(x, cols) {
  cols <- cols[cols %in% names(x)]
  if (length(cols) == 0) return(rep(NA_real_, nrow(x)))
  mat <- do.call(cbind, lapply(cols, function(nm) to_numeric_safe(x[[nm]])))
  mat[mat < 0 | mat > 1] <- NA_real_
  suppressWarnings(apply(mat, 1, function(z) {
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  }))
}

technical_floor_pass <- function(x, cfg) {
  tf <- cfg$technical_filters
  min_alt <- ifelse(x$is_indel, tf$min_alt_count_indel, tf$min_alt_count_snv)
  min_af <- ifelse(x$is_indel, tf$min_af_indel, tf$min_af_snv)
  pass_filter <- is.na(x$filter_status) |
    x$filter_status %in% c("PASS", ".", "NA") |
    x$filter_status == ""
  pass_filter &
    (is.na(x$dp) | x$dp >= tf$min_depth) &
    (is.na(x$alt_count) | x$alt_count >= min_alt) &
    (is.na(x$vaf) | x$vaf >= min_af) &
    (is.na(x$tlod) | x$tlod >= tf$min_tlod) &
    (is.na(x$mbq) | x$mbq >= tf$min_mbq) &
    (is.na(x$mmq) | x$mmq >= tf$min_mmq)
}

compute_cohort_recurrence <- function(x, cfg) {
  dt <- data.table::as.data.table(x)
  total_samples <- data.table::uniqueN(dt$sample_id[!is.na(dt$sample_id)])
  tumor_type_samples <- dt[!is.na(sample_id), .(
    tumor_type_total_samples = data.table::uniqueN(sample_id)
  ), by = tumor_type]

  by_variant <- dt[, .(
    variant_n_samples = data.table::uniqueN(sample_id),
    variant_median_vaf = safe_median(vaf),
    variant_median_alt_count = safe_median(alt_count),
    variant_median_dp = safe_median(dp)
  ), by = .(variant_id, chrom, pos, ref, alt)]
  by_variant[, variant_cohort_freq := variant_n_samples / total_samples]

  by_locus <- dt[, .(
    locus_n_samples = data.table::uniqueN(sample_id)
  ), by = locus_id]
  by_locus[, locus_cohort_freq := locus_n_samples / total_samples]

  by_variant_tumor <- dt[, .(
    variant_tumor_type_n_samples = data.table::uniqueN(sample_id)
  ), by = .(variant_id, tumor_type)]
  by_variant_tumor <- merge(by_variant_tumor, tumor_type_samples, by = "tumor_type", all.x = TRUE)
  by_variant_tumor[, variant_tumor_type_freq := variant_tumor_type_n_samples / tumor_type_total_samples]

  by_locus_tumor <- dt[, .(
    locus_tumor_type_n_samples = data.table::uniqueN(sample_id)
  ), by = .(locus_id, tumor_type)]
  by_locus_tumor <- merge(by_locus_tumor, tumor_type_samples, by = "tumor_type", all.x = TRUE)
  by_locus_tumor[, locus_tumor_type_freq := locus_tumor_type_n_samples / tumor_type_total_samples]

  list(
    by_variant = as.data.frame(by_variant),
    by_locus = as.data.frame(by_locus),
    by_variant_tumor = as.data.frame(by_variant_tumor),
    by_locus_tumor = as.data.frame(by_locus_tumor),
    tumor_type_samples = as.data.frame(tumor_type_samples),
    total_samples = total_samples
  )
}

merge_recurrence_features <- function(x, recurrence) {
  dt <- data.table::as.data.table(x)
  dt <- merge(dt, data.table::as.data.table(recurrence$by_variant),
              by = c("variant_id", "chrom", "pos", "ref", "alt"), all.x = TRUE, sort = FALSE)
  dt <- merge(dt, data.table::as.data.table(recurrence$by_locus),
              by = "locus_id", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, data.table::as.data.table(recurrence$by_variant_tumor),
              by = c("variant_id", "tumor_type"), all.x = TRUE, sort = FALSE)
  dt <- merge(dt, data.table::as.data.table(recurrence$by_locus_tumor),
              by = c("locus_id", "tumor_type"), all.x = TRUE, sort = FALSE)
  dt[, total_samples := recurrence$total_samples]
  as.data.frame(dt)
}
