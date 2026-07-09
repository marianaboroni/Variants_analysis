# Brazilian-population-aware population evidence for tumor-only curation.
#
# Rationale: the Brazilian population is highly admixed and under-represented in
# global databases, so ABSENCE from global databases (e.g. gnomAD) is NOT evidence
# of somaticity. Global and Brazilian evidence are kept SEPARATE, coverage/absence
# interpretability is explicit, and a variant present in a Brazilian database (e.g.
# ABraOM SABE) is treated as possible germline -> REVIEW (documented, auditable
# population rule; configurable). Matching is by the canonical key
# BUILD|CHROM|POS|REF|ALT. The Brazilian DB is user-provided (not redistributed).
# All thresholds are heuristic and cohort-size dependent.

POP_EVIDENCE_STATUSES <- c("common_global", "common_brazilian", "rare_global_and_brazilian",
  "absent_global_present_brazilian", "absent_brazilian_present_global",
  "absent_in_checked_databases", "insufficient_coverage", "not_evaluable")

#' Add Brazilian-aware population evidence columns to a variant table.
#'
#' @param x variant data.frame (after features; uses max_pop_af, vaf, gene, dp).
#' @param cfg resolved config (reads `population:`).
#' @param build resolved genome build.
#' @return `x` with GLOBAL/BRAZILIAN/POPULATION_EVIDENCE_* statuses, Brazilian
#'   AF/AC/AN/HOM, coverage/absence-interpretability, POSSIBLE_BRAZILIAN_GERMLINE,
#'   BRAZILIAN_POPULATION_EVIDENCE_SCORE, and `population_review_flag`.
#' @keywords internal
add_brazilian_population_evidence <- function(x, cfg = NULL, build = "GRCh38") {
  n <- nrow(x)
  g_common <- cfg_get(cfg, c("population", "global", "common_af"),
                      cfg_get(cfg, c("population_filters", "common_af_threshold"), 0.01))
  g_rare <- cfg_get(cfg, c("population", "global", "rare_af"),
                    cfg_get(cfg, c("population_filters", "low_frequency_af_threshold"), 0.001))
  b_common <- cfg_get(cfg, c("population", "brazilian", "common_af"), 0.01)
  b_review <- cfg_get(cfg, c("population", "brazilian", "review_af"), 0.001)
  b_min_ac <- cfg_get(cfg, c("population", "brazilian", "minimum_allele_count"), 2)
  lc_af <- cfg_get(cfg, c("population", "local_controls", "artifact_or_germline_af"), 0.005)

  global_af <- if ("max_pop_af" %in% names(x)) suppressWarnings(as.numeric(x$max_pop_af)) else rep(NA_real_, n)
  x$MAX_GLOBAL_AF <- global_af
  x$MAX_ANCESTRY_SPECIFIC_AF <- brazil_ancestry_max_af(x)

  # ---- Brazilian database (ABraOM etc.) matching by canonical key ---------
  bdb <- load_brazilian_db(cfg, build)
  x$POPULATION_DATABASES_CHECKED <- paste(c("gnomAD_global",
    if (!is.null(bdb)) attr(bdb, "db_name") else NULL), collapse = ";")
  x$BRAZILIAN_AF <- NA_real_; x$BRAZILIAN_AC <- NA_real_; x$BRAZILIAN_AN <- NA_real_
  x$BRAZILIAN_HOM_COUNT <- NA_real_; x$BRAZILIAN_DATABASE_SAMPLE_SIZE <- NA_real_
  x$BRAZILIAN_DB_VARIANT_OBSERVED <- FALSE
  x$BRAZILIAN_DB_REGION_COVERED <- NA
  x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE <- FALSE
  if (!is.null(bdb)) {
    keys <- canonical_variant_key(build, x$chrom, x$pos, x$ref, x$alt)
    idx <- match(keys, bdb$canonical_key)
    hit <- !is.na(idx)
    x$BRAZILIAN_DB_VARIANT_OBSERVED <- hit
    x$BRAZILIAN_AF[hit] <- bdb$af[idx[hit]]
    x$BRAZILIAN_AC[hit] <- bdb$ac[idx[hit]]
    x$BRAZILIAN_AN[hit] <- bdb$an[idx[hit]]
    x$BRAZILIAN_HOM_COUNT[hit] <- bdb$hom[idx[hit]]
    x$BRAZILIAN_DATABASE_SAMPLE_SIZE <- attr(bdb, "sample_size")
    # coverage: if the DB provides a callable/covered flag per key, use it; else
    # treat observed variants as covered and rely on region coverage table if present
    x$BRAZILIAN_DB_REGION_COVERED <- ifelse(hit, TRUE, NA)  # unknown coverage where not observed
    # absence is interpretable only where the region is known-covered
    x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE <- !hit & isTRUE_vec(x$BRAZILIAN_DB_REGION_COVERED %in% TRUE)
  }
  x$LOCAL_CONTROL_AF <- brazil_local_control_af(x, cfg, build)

  # ---- statuses -----------------------------------------------------------
  x$GLOBAL_POPULATION_STATUS <- global_status(global_af, g_common, g_rare, x)
  x$BRAZILIAN_POPULATION_STATUS <- brazilian_status(x, bdb, b_common, b_review, b_min_ac)
  x$POPULATION_COVERAGE_STATUS <- ifelse(is.null(bdb), "brazilian_db_not_configured",
    ifelse(x$BRAZILIAN_DB_VARIANT_OBSERVED, "observed",
      ifelse(x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE, "covered_not_observed", "coverage_unknown")))
  ev <- population_evidence_status(x, has_bdb = !is.null(bdb))
  x$POPULATION_EVIDENCE_STATUS <- ev$status
  x$POPULATION_EVIDENCE_CONFIDENCE <- ev$confidence

  # ---- possible Brazilian germline (documented review rule) ---------------
  br <- brazilian_germline_flags(x, cfg, b_review)
  x$POSSIBLE_BRAZILIAN_GERMLINE <- br$flag
  x$BRAZILIAN_GERMLINE_EVIDENCE <- br$evidence
  x$BRAZILIAN_GERMLINE_REVIEW_REASON <- br$reason
  x$population_review_flag <- br$flag

  x$BRAZILIAN_POPULATION_EVIDENCE_SCORE <- brazilian_population_evidence_score(x, b_common)
  x
}

#' @keywords internal
brazil_ancestry_max_af <- function(x) {
  cols <- intersect(c("AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
    "gnomADg_AFR_AF", "gnomADg_AMR_AF", "gnomADg_EAS_AF", "gnomADg_NFE_AF", "gnomADg_SAS_AF",
    "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_EAS_AF", "gnomADe_NFE_AF", "gnomADe_SAS_AF"), names(x))
  if (length(cols) == 0) return(rep(NA_real_, nrow(x)))
  m <- suppressWarnings(do.call(cbind, lapply(cols, function(c0) as.numeric(x[[c0]]))))
  suppressWarnings(apply(m, 1, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)))
}

#' @keywords internal
brazil_local_control_af <- function(x, cfg, build) {
  path <- cfg_get(cfg, c("population", "local_controls", "path"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) return(rep(NA_real_, nrow(x)))
  db <- read_variants(path, "\t")
  db$canonical_key <- canonical_variant_key(build,
    coalesce_columns(db, c("CHROM", "chrom", "Chromosome")),
    coalesce_columns(db, c("POS", "pos", "Start_Position")),
    coalesce_columns(db, c("REF", "ref", "Reference_Allele")),
    coalesce_columns(db, c("ALT", "alt", "Tumor_Seq_Allele2")))
  db$af <- to_numeric_safe(coalesce_columns(db, c("AF", "af", "control_af")))
  keys <- canonical_variant_key(build, x$chrom, x$pos, x$ref, x$alt)
  db$af[match(keys, db$canonical_key)]
}

#' Load a Brazilian population database (ABraOM/BIPMed/... user-provided TSV).
#' @keywords internal
load_brazilian_db <- function(cfg, build) {
  path <- cfg_get(cfg, c("population", "brazilian", "db"),
                  cfg_get(cfg, c("population", "brazilian", "path"), NULL))
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  db_build <- cfg_get(cfg, c("population", "brazilian", "genome_build"), build)
  if (!identical(db_build, build))
    warning(sprintf("Brazilian DB build (%s) != analysis build (%s); matching may be incomplete. Prepare a lifted copy.",
                    db_build, build), call. = FALSE)
  d <- read_variants(path, "\t")
  out <- data.frame(
    canonical_key = canonical_variant_key(build,
      coalesce_columns(d, c("CHROM", "chrom", "Chromosome", "CHROMOSOME")),
      coalesce_columns(d, c("POS", "pos", "Start_Position", "START", "GENOME_START")),
      coalesce_columns(d, c("REF", "ref", "Reference_Allele", "GENOMIC_WT_ALLELE")),
      coalesce_columns(d, c("ALT", "alt", "Tumor_Seq_Allele2", "GENOMIC_MUT_ALLELE"))),
    af = to_numeric_safe(coalesce_columns(d, c("AF", "af", "ABraOM_AF", "brazilian_af", "AF_sabe"))),
    ac = to_numeric_safe(coalesce_columns(d, c("AC", "ac", "allele_count"))),
    an = to_numeric_safe(coalesce_columns(d, c("AN", "an", "allele_number"))),
    hom = to_numeric_safe(coalesce_columns(d, c("HOM", "hom", "nhomalt", "n_hom"))),
    stringsAsFactors = FALSE)
  out <- out[!duplicated(out$canonical_key), , drop = FALSE]
  attr(out, "db_name") <- cfg_get(cfg, c("population", "brazilian", "name"), "brazilian_db")
  attr(out, "sample_size") <- cfg_get(cfg, c("population", "brazilian", "cohort_size"), NA_real_)
  out
}

#' @keywords internal
global_status <- function(af, common, rare, x) {
  has_cov <- if ("dp" %in% names(x)) !is.na(x$dp) else rep(TRUE, length(af))
  st <- rep("not_evaluable", length(af))
  st[!is.na(af) & af >= common] <- "common_global"
  st[!is.na(af) & af >= rare & af < common] <- "rare_global"
  st[!is.na(af) & af < rare] <- "absent_or_ultrarare_global"
  st[is.na(af) & has_cov] <- "absent_global"
  st
}

#' @keywords internal
brazilian_status <- function(x, bdb, common, review, min_ac) {
  n <- nrow(x)
  if (is.null(bdb)) return(rep("not_evaluable", n))
  af <- x$BRAZILIAN_AF; ac <- x$BRAZILIAN_AC; obs <- isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED)
  st <- rep("absent_in_checked_databases", n)
  st[!obs & !isTRUE_vec(x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE)] <- "insufficient_coverage"
  st[obs & !is.na(af) & af >= common] <- "common_brazilian"
  st[obs & !is.na(af) & af >= review & af < common] <- "rare_brazilian"
  st[obs & (is.na(af) | af < review) & !is.na(ac) & ac >= min_ac] <- "rare_brazilian"
  st[obs & is.na(af) & (is.na(ac) | ac < min_ac)] <- "present_low_count_brazilian"
  st
}

#' @keywords internal
population_evidence_status <- function(x, has_bdb) {
  n <- nrow(x)
  g <- x$GLOBAL_POPULATION_STATUS; b <- x$BRAZILIAN_POPULATION_STATUS
  st <- rep("not_evaluable", n); conf <- rep("low", n)
  st[g == "common_global"] <- "common_global"
  st[b %in% c("common_brazilian")] <- "common_brazilian"
  st[g %in% c("rare_global") & b %in% c("rare_brazilian", "present_low_count_brazilian")] <- "rare_global_and_brazilian"
  st[g %in% c("absent_global", "absent_or_ultrarare_global") & isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED)] <- "absent_global_present_brazilian"
  st[g %in% c("common_global", "rare_global") & has_bdb & !isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED) & isTRUE_vec(x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE)] <- "absent_brazilian_present_global"
  st[g %in% c("absent_global", "absent_or_ultrarare_global") & has_bdb & isTRUE_vec(x$BRAZILIAN_DB_ABSENCE_INTERPRETABLE) & !isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED)] <- "absent_in_checked_databases"
  st[has_bdb & x$POPULATION_COVERAGE_STATUS == "coverage_unknown" & !isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED) & g %in% c("absent_global","absent_or_ultrarare_global")] <- "insufficient_coverage"
  conf[st %in% c("common_global", "common_brazilian", "rare_global_and_brazilian")] <- "high"
  conf[st %in% c("absent_global_present_brazilian", "absent_brazilian_present_global")] <- "moderate"
  conf[st %in% c("insufficient_coverage", "not_evaluable")] <- "low"
  list(status = st, confidence = conf)
}

#' @keywords internal
brazilian_germline_flags <- function(x, cfg, review) {
  n <- nrow(x)
  vaf <- if ("vaf" %in% names(x)) suppressWarnings(as.numeric(x$vaf)) else rep(NA_real_, n)
  het_like <- !is.na(vaf) & vaf >= 0.35 & vaf <= 0.65
  present_br <- isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED)
  absent_global <- x$GLOBAL_POPULATION_STATUS %in% c("absent_global", "absent_or_ultrarare_global")
  local_ctrl <- !is.na(x$LOCAL_CONTROL_AF) & x$LOCAL_CONTROL_AF >= cfg_get(cfg, c("population", "local_controls", "artifact_or_germline_af"), 0.005)
  predisp <- brazil_is_predisposition_gene(x, cfg)
  ev <- vector("list", n); flag <- logical(n)
  for (i in seq_len(n)) {
    e <- character()
    if (absent_global[i] && present_br[i]) e <- c(e, "absent_global_present_brazilian")
    if (present_br[i] && het_like[i]) e <- c(e, "brazilian_db_with_heterozygous_vaf")
    if (isTRUE(local_ctrl[i])) e <- c(e, "present_in_local_controls")
    if (isTRUE(predisp[i]) && (present_br[i] || het_like[i])) e <- c(e, "predisposition_gene_review")
    flag[i] <- length(e) > 0
    ev[[i]] <- e
  }
  list(flag = flag,
       evidence = vapply(ev, function(e) paste(e, collapse = ";"), character(1)),
       reason = vapply(ev, function(e) if (length(e)) "possible_brazilian_germline" else NA_character_, character(1)))
}

#' @keywords internal
brazil_is_predisposition_gene <- function(x, cfg) {
  genes <- cfg_get(cfg, c("population", "predisposition_genes"),
    c("BRCA1", "BRCA2", "TP53", "PALB2", "ATM", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2",
      "APC", "MUTYH", "PTEN", "STK11", "CDH1", "RB1", "VHL", "RET", "NF1", "NF2"))
  g <- toupper(as.character(if ("gene" %in% names(x)) x$gene else rep(NA, nrow(x))))
  g %in% toupper(genes)
}

#' Separate Brazilian population evidence score (0..1). Higher = more germline-like
#' / lower confidence in somaticity. NOT a biological-relevance score.
#' @keywords internal
brazilian_population_evidence_score <- function(x, b_common) {
  s <- rep(0, nrow(x))
  s[isTRUE_vec(x$BRAZILIAN_DB_VARIANT_OBSERVED)] <- 0.5
  s[!is.na(x$BRAZILIAN_AF) & x$BRAZILIAN_AF >= 0.001] <- 0.7
  s[!is.na(x$BRAZILIAN_AF) & x$BRAZILIAN_AF >= b_common] <- 0.95
  s[!is.na(x$LOCAL_CONTROL_AF) & x$LOCAL_CONTROL_AF >= 0.005] <- pmax(s[!is.na(x$LOCAL_CONTROL_AF) & x$LOCAL_CONTROL_AF >= 0.005], 0.6)
  round(s, 4)
}

#' Apply the documented Brazilian population REVIEW rule (escalate PASS -> REVIEW
#' for possible Brazilian germline). Explicit and auditable; configurable.
#' @keywords internal
apply_brazilian_population_rules <- function(x, cfg = NULL) {
  if (!isTRUE(cfg_get(cfg, c("population", "brazilian", "escalate_to_review"), TRUE))) return(x)
  if (!"POSSIBLE_BRAZILIAN_GERMLINE" %in% names(x) || !"STATUS" %in% names(x)) return(x)
  esc <- isTRUE_vec(x$POSSIBLE_BRAZILIAN_GERMLINE) & x$STATUS == "PASS"
  if (any(esc)) {
    x$STATUS[esc] <- "REVIEW"
    if ("filter_status" %in% names(x)) x$filter_status[esc] <- "REVIEW"
    rr <- "present_in_brazilian_population_database;possible_germline"
    x$status_reason[esc] <- ifelse(is.na(x$status_reason[esc]) | x$status_reason[esc] == "",
      rr, paste(x$status_reason[esc], rr, sep = "; "))
    if ("review_reasons" %in% names(x))
      x$review_reasons[esc] <- ifelse(x$review_reasons[esc] == "", "possible_brazilian_germline",
        paste(x$review_reasons[esc], "possible_brazilian_germline", sep = ";"))
    log_step("population", "Brazilian germline REVIEW rule applied", escalated = sum(esc))
  }
  x
}
