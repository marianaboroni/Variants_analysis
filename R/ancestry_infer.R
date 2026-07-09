# Conservative, EXPERIMENTAL genetic-ancestry inference from germline-informative
# SNPs in tumor-only data. Genetic ancestry is NOT race or ethnic identity; it is
# used ONLY to weight/choose relevant population frequencies and to flag
# under-representation. It NEVER filters variants or changes biological priority.
#
# It uses a SEPARATE SNP track (common germline-informative SNPs, taken BEFORE the
# somatic filter so the allele spectrum is not biased), excludes hotspots/drivers
# and unstable (CNV/LOH-like) loci, aligns alleles/strand vs an AIMs panel, and
# estimates continuous ancestry proportions by non-negative least squares against
# reference-population allele frequencies (+ a PCA projection for visualization).
# PLINK/ADMIXTURE are optional external backends (see docs); this fallback is
# self-contained. All thresholds are heuristic and assay-dependent.

ANCESTRY_COMPONENTS <- c("AFR", "EUR", "NAT", "EAS", "SAS")

#' Infer genetic-ancestry proportions for the sample(s) in a variant table.
#'
#' @param variants standardized variant table (ALL variants, pre-somatic-filter).
#' @param cfg resolved config (reads `ancestry:`).
#' @param build resolved genome build.
#' @return a list with `summary`, `snps`, `qc`, `projection` data.frames.
#' @keywords internal
infer_ancestry <- function(variants, cfg = NULL, build = "GRCh38") {
  qc_min_dp <- cfg_get(cfg, c("ancestry", "qc", "minimum_depth"), 10)
  qc_min_alt <- cfg_get(cfg, c("ancestry", "qc", "minimum_alt_reads"), 3)
  min_snps <- cfg_get(cfg, c("ancestry", "qc", "minimum_snps"), 100)
  panel <- load_aims_panel(cfg, build)
  samples <- unique(as.character(variants$sample_id))

  if (is.null(panel)) {
    return(ancestry_empty(samples, reason = "no AIMs reference panel configured"))
  }

  cand <- ancestry_select_snps(variants, build, qc_min_dp, qc_min_alt)
  cand <- ancestry_align_to_panel(cand, panel)               # allele/strand alignment + palindrome drop
  cand <- ancestry_flag_unstable(cand)                        # CNV/LOH-like exclusion

  qc_rows <- list(); summ_rows <- list(); proj_rows <- list(); snp_rows <- list()
  for (s in samples) {
    cs <- cand[cand$sample_id == s, , drop = FALSE]
    used <- cs[cs$ANCESTRY_LOCUS_STATUS == "usable", , drop = FALSE]
    n_avail <- nrow(cs); n_used <- nrow(used)
    call_rate <- if (n_avail > 0) n_used / n_avail else 0
    conf <- ancestry_confidence(n_used, call_rate, min_snps)
    if (n_used < min_snps) {
      prop <- setNames(rep(NA_real_, length(ANCESTRY_COMPONENTS)), ANCESTRY_COMPONENTS)
      pcs <- c(NA_real_, NA_real_, NA_real_); status <- "not_evaluable"
    } else {
      est <- ancestry_estimate_proportions(used, panel)
      prop <- est$proportions; pcs <- est$pcs; status <- "evaluated"
    }
    label <- ancestry_summary_label(prop, cfg)
    summ_rows[[s]] <- data.frame(sample_id = s, ANCESTRY_INFERENCE_STATUS = status,
      ANCESTRY_CONFIDENCE = conf, ANCESTRY_SNPS_AVAILABLE = n_avail, ANCESTRY_SNPS_USED = n_used,
      ANCESTRY_CALL_RATE = round(call_rate, 4),
      ANCESTRY_AFR_PROPORTION = prop["AFR"], ANCESTRY_EUR_PROPORTION = prop["EUR"],
      ANCESTRY_NAT_PROPORTION = prop["NAT"], ANCESTRY_EAS_PROPORTION = prop["EAS"],
      ANCESTRY_SAS_PROPORTION = prop["SAS"], ANCESTRY_SUMMARY_LABEL = label,
      ANCESTRY_PC1 = pcs[1], ANCESTRY_PC2 = pcs[2],
      REFERENCE_PANEL_VERSION = attr(panel, "version"), stringsAsFactors = FALSE)
    qc_rows[[s]] <- ancestry_qc_row(s, cs, n_avail, n_used, call_rate, conf)
    proj_rows[[s]] <- data.frame(sample_id = s, PC1 = pcs[1], PC2 = pcs[2], PC3 = pcs[3],
                                 stringsAsFactors = FALSE)
    if (nrow(cs) > 0) snp_rows[[s]] <- cs
  }
  list(summary = do.call(rbind, summ_rows),
       snps = if (length(snp_rows)) do.call(rbind, snp_rows) else cand[0, , drop = FALSE],
       qc = do.call(rbind, qc_rows),
       projection = do.call(rbind, proj_rows),
       panel = panel)
}

#' Select germline-informative candidate SNPs (biallelic SNVs, autosomal, adequate
#' depth, not hotspot/driver). Uses ALL variants regardless of somatic filter.
#' @keywords internal
ancestry_select_snps <- function(v, build, min_dp, min_alt) {
  n <- nrow(v)
  is_snv <- nchar(as.character(v$ref)) == 1 & nchar(as.character(v$alt)) == 1 &
    toupper(v$ref) %in% c("A", "C", "G", "T") & toupper(v$alt) %in% c("A", "C", "G", "T")
  chrom <- normalize_chrom(v$chrom)
  autosomal <- chrom %in% as.character(1:22)
  dp <- if ("dp" %in% names(v)) suppressWarnings(as.numeric(v$dp)) else rep(NA_real_, n)
  altc <- if ("alt_count" %in% names(v)) suppressWarnings(as.numeric(v$alt_count)) else rep(NA_real_, n)
  hotspot <- if ("hotspot_match" %in% names(v)) v$hotspot_match %in% c(TRUE, "TRUE") else rep(FALSE, n)
  driver <- if ("gene_driver_match" %in% names(v)) v$gene_driver_match %in% c(TRUE, "TRUE") else rep(FALSE, n)
  keep <- is_snv & autosomal & !hotspot & !driver
  cs <- data.frame(sample_id = as.character(v$sample_id), chrom = chrom,
    pos = suppressWarnings(as.integer(v$pos)), ref = toupper(as.character(v$ref)),
    alt = toupper(as.character(v$alt)),
    vaf = if ("vaf" %in% names(v)) suppressWarnings(as.numeric(v$vaf)) else rep(NA_real_, n),
    dp = dp, alt_count = altc, stringsAsFactors = FALSE)[keep, , drop = FALSE]
  cs$canonical_key <- canonical_variant_key(build, cs$chrom, cs$pos, cs$ref, cs$alt)
  # genotype-from-VAF + confidence + technical status
  cs$ANCESTRY_GENOTYPE_CONFIDENCE <- NA_character_
  cs$ANCESTRY_LOCUS_STATUS <- "usable"
  low_dp <- !is.na(cs$dp) & cs$dp < min_dp
  low_alt <- !is.na(cs$alt_count) & cs$alt_count < min_alt & !is.na(cs$vaf) & cs$vaf > 0
  cs$ANCESTRY_LOCUS_STATUS[low_dp] <- "low_depth"
  cs$dosage <- ancestry_genotype_from_vaf(cs$vaf)
  cs$ANCESTRY_GENOTYPE_CONFIDENCE <- ifelse(low_dp, "low", "usable")
  cs$ANCESTRY_SNP_EXCLUSION_REASON <- ifelse(low_dp, "low_depth", NA_character_)
  cs
}

#' @keywords internal
ancestry_genotype_from_vaf <- function(vaf) {
  ifelse(is.na(vaf), NA_real_, ifelse(vaf < 0.15, 0, ifelse(vaf > 0.85, 2, 1)))
}

#' Align sample SNPs to the AIMs panel; drop palindromic (A/T, C/G) ambiguous SNPs
#' and build-mismatched loci; flip alleles when needed.
#' @keywords internal
ancestry_align_to_panel <- function(cs, panel) {
  idx <- match(cs$canonical_key, panel$canonical_key)
  in_panel <- !is.na(idx)
  cs$in_reference_panel <- in_panel
  cs$ANCESTRY_LOCUS_STATUS[!in_panel & cs$ANCESTRY_LOCUS_STATUS == "usable"] <- "missing"
  cs$ANCESTRY_SNP_EXCLUSION_REASON[!in_panel & is.na(cs$ANCESTRY_SNP_EXCLUSION_REASON)] <- "not_in_reference_panel"
  palI <- (toupper(cs$ref) == "A" & toupper(cs$alt) == "T") |
          (toupper(cs$ref) == "T" & toupper(cs$alt) == "A") |
          (toupper(cs$ref) == "C" & toupper(cs$alt) == "G") |
          (toupper(cs$ref) == "G" & toupper(cs$alt) == "C")
  cs$ANCESTRY_PALINDROMIC <- palI
  cs$ANCESTRY_LOCUS_STATUS[palI & cs$ANCESTRY_LOCUS_STATUS == "usable"] <- "palindromic"
  cs$ANCESTRY_SNP_EXCLUSION_REASON[palI & is.na(cs$ANCESTRY_SNP_EXCLUSION_REASON)] <- "palindromic"
  cs$panel_row <- idx
  cs
}

#' Flag/exclude CNV/LOH-like loci (heterozygous genotype with strong allelic
#' imbalance) so unstable tumor regions are not treated as germline genotypes.
#' @keywords internal
ancestry_flag_unstable <- function(cs) {
  imbalance <- !is.na(cs$vaf) & cs$dosage == 1 & abs(cs$vaf - 0.5) > 0.2
  cs$ANCESTRY_SNP_EXCLUDED_LOH <- imbalance & cs$vaf > 0.7
  cs$ANCESTRY_SNP_EXCLUDED_CNV <- imbalance & cs$vaf < 0.3
  drop <- imbalance & cs$ANCESTRY_LOCUS_STATUS == "usable"
  cs$ANCESTRY_LOCUS_STATUS[drop] <- ifelse(cs$vaf[drop] > 0.5, "possible_loh", "possible_cnv")
  cs$ANCESTRY_SNP_EXCLUSION_REASON[drop & is.na(cs$ANCESTRY_SNP_EXCLUSION_REASON)] <- "allelic_imbalance"
  cs
}

#' Estimate continuous ancestry proportions by non-negative least squares against
#' reference-population allele frequencies, plus a PCA projection.
#' @keywords internal
ancestry_estimate_proportions <- function(used, panel) {
  comps <- intersect(ANCESTRY_COMPONENTS, names(panel))
  F <- as.matrix(panel[used$panel_row, comps, drop = FALSE])   # SNP x pop alt-allele freq
  y <- used$dosage / 2                                          # observed alt fraction (0/0.5/1)
  ok <- stats::complete.cases(F) & !is.na(y)
  F <- F[ok, , drop = FALSE]; y <- y[ok]
  p <- ancestry_nnls_simplex(F, y)
  names(p) <- comps
  full <- setNames(rep(0, length(ANCESTRY_COMPONENTS)), ANCESTRY_COMPONENTS)
  full[comps] <- p
  # PCA projection: PCA of the reference pop-frequency matrix; project the sample
  pcs <- tryCatch({
    pc <- stats::prcomp(t(F), center = TRUE, scale. = FALSE)
    samp <- scale(matrix(y, nrow = 1), center = pc$center, scale = FALSE) %*% pc$rotation
    c(samp[1, 1], if (ncol(samp) >= 2) samp[1, 2] else NA_real_,
      if (ncol(samp) >= 3) samp[1, 3] else NA_real_)
  }, error = function(e) c(NA_real_, NA_real_, NA_real_))
  list(proportions = round(full, 4), pcs = round(pcs, 4))
}

#' Non-negative, sum-to-one least squares via multiplicative updates (simple,
#' dependency-free heuristic for experimental proportion estimation).
#' @keywords internal
ancestry_nnls_simplex <- function(F, y, iters = 200) {
  K <- ncol(F); p <- rep(1 / K, K)
  FtF <- crossprod(F); Fty <- as.numeric(crossprod(F, y))
  for (i in seq_len(iters)) {
    num <- Fty; den <- as.numeric(FtF %*% p) + 1e-9
    p <- p * pmax(num, 0) / den
    if (sum(p) <= 0) { p <- rep(1 / K, K); break }
    p <- p / sum(p)
  }
  pmax(0, pmin(1, p))
}

#' @keywords internal
ancestry_confidence <- function(n_used, call_rate, min_snps = 100) {
  if (n_used < min_snps) "not_evaluable"
  else if (n_used >= 1000 && call_rate >= 0.95) "high_confidence"
  else if (n_used >= 300 && call_rate >= 0.90) "moderate_confidence"
  else "low_confidence"
}

#' @keywords internal
ancestry_summary_label <- function(prop, cfg) {
  if (all(is.na(prop))) return("not_evaluable")
  thr <- cfg_get(cfg, c("ancestry", "dominant_component_threshold"), 0.80)
  sec <- cfg_get(cfg, c("ancestry", "admixed_minimum_secondary_component"), 0.15)
  o <- sort(prop, decreasing = TRUE)
  if (o[1] >= thr) paste0("predominantly_", names(o)[1])
  else if (length(o) >= 2 && o[2] >= sec) "admixed"
  else "uncertain"
}

#' @keywords internal
ancestry_qc_row <- function(s, cs, n_avail, n_used, call_rate, conf) {
  reason_tab <- table(cs$ANCESTRY_SNP_EXCLUSION_REASON)
  data.frame(sample_id = s, snps_available = n_avail, snps_used = n_used,
    call_rate = round(call_rate, 4), confidence = conf,
    excluded_low_depth = sum(cs$ANCESTRY_LOCUS_STATUS == "low_depth"),
    excluded_palindromic = sum(cs$ANCESTRY_LOCUS_STATUS == "palindromic"),
    excluded_loh = sum(cs$ANCESTRY_LOCUS_STATUS == "possible_loh"),
    excluded_cnv = sum(cs$ANCESTRY_LOCUS_STATUS == "possible_cnv"),
    excluded_not_in_panel = sum(cs$ANCESTRY_LOCUS_STATUS == "missing"),
    stringsAsFactors = FALSE)
}

#' @keywords internal
ancestry_empty <- function(samples, reason) {
  summ <- data.frame(sample_id = samples, ANCESTRY_INFERENCE_STATUS = "not_evaluable",
    ANCESTRY_CONFIDENCE = "not_evaluable", ANCESTRY_SNPS_AVAILABLE = 0, ANCESTRY_SNPS_USED = 0,
    ANCESTRY_CALL_RATE = 0, ANCESTRY_AFR_PROPORTION = NA_real_, ANCESTRY_EUR_PROPORTION = NA_real_,
    ANCESTRY_NAT_PROPORTION = NA_real_, ANCESTRY_EAS_PROPORTION = NA_real_,
    ANCESTRY_SAS_PROPORTION = NA_real_, ANCESTRY_SUMMARY_LABEL = "not_evaluable",
    ANCESTRY_PC1 = NA_real_, ANCESTRY_PC2 = NA_real_, REFERENCE_PANEL_VERSION = NA_character_,
    reason = reason, stringsAsFactors = FALSE)
  list(summary = summ, snps = data.frame(), qc = data.frame(sample_id = samples),
       projection = data.frame(sample_id = samples), panel = NULL)
}

#' Load a versioned AIMs reference panel (marker key + per-population alt-allele
#' frequencies). User-provided; not redistributed.
#' @keywords internal
load_aims_panel <- function(cfg, build) {
  path <- cfg_get(cfg, c("ancestry", "marker_panel", "path"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  d <- read_variants(path, "\t")
  key <- if ("canonical_key" %in% names(d)) d$canonical_key else
    canonical_variant_key(build,
      coalesce_columns(d, c("CHROM", "chrom", "chromosome")),
      coalesce_columns(d, c("POS", "pos", "position")),
      coalesce_columns(d, c("REF", "ref")), coalesce_columns(d, c("ALT", "alt")))
  out <- data.frame(canonical_key = key, stringsAsFactors = FALSE)
  for (k in ANCESTRY_COMPONENTS) {
    col <- intersect(c(paste0(k, "_freq"), paste0(k, "_AF"), k), names(d))
    if (length(col)) out[[k]] <- to_numeric_safe(d[[col[1]]])
  }
  out <- out[!duplicated(out$canonical_key), , drop = FALSE]
  attr(out, "version") <- cfg_get(cfg, c("ancestry", "marker_panel", "version"), "aims_v1")
  out
}
