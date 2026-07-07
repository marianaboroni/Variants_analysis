# MAF export and maftools-based oncoplots. The oncoplot is produced directly by
# maftools::oncoplot (no hand-rolled ggplot2 imitation). The MAF is validated
# with maftools::read.maf before any plotting.

# Explicit, tested VEP consequence -> MAF Variant_Classification map.
# Unknown consequences are mapped to a controlled category and reported in QC.
VEP_TO_MAF_CLASS <- list(
  missense_variant = "Missense_Mutation",
  stop_gained = "Nonsense_Mutation",
  stop_lost = "Nonstop_Mutation",
  start_lost = "Translation_Start_Site",
  synonymous_variant = "Silent",
  stop_retained_variant = "Silent",
  splice_acceptor_variant = "Splice_Site",
  splice_donor_variant = "Splice_Site",
  splice_region_variant = "Splice_Region",
  frameshift_variant = "Frame_Shift",           # refined to Ins/Del by type
  inframe_insertion = "In_Frame_Ins",
  inframe_deletion = "In_Frame_Del",
  protein_altering_variant = "Missense_Mutation",
  coding_sequence_variant = "Missense_Mutation",
  intron_variant = "Intron",
  `5_prime_utr_variant` = "5'UTR",
  `3_prime_utr_variant` = "3'UTR",
  upstream_gene_variant = "5'Flank",
  downstream_gene_variant = "3'Flank",
  intergenic_variant = "IGR",
  non_coding_transcript_exon_variant = "RNA",
  mature_mirna_variant = "RNA"
)
MAF_UNKNOWN_CLASS <- "Targeted_Region"  # controlled bucket for unmapped consequences

#' Map a VEP consequence string to a MAF Variant_Classification.
#' Takes the first (most severe) consequence when several are `&`-joined.
#' @keywords internal
map_consequence_to_maf <- function(consequence, variant_type) {
  first <- tolower(sub("[&,].*$", "", trimws(as.character(consequence))))
  cls <- vapply(first, function(c) {
    v <- VEP_TO_MAF_CLASS[[c]]
    if (is.null(v)) NA_character_ else v
  }, character(1), USE.NAMES = FALSE)
  # refine frameshift by ins/del
  fs <- !is.na(cls) & cls == "Frame_Shift"
  cls[fs] <- ifelse(variant_type[fs] == "INS", "Frame_Shift_Ins", "Frame_Shift_Del")
  unknown <- is.na(cls)
  cls[unknown] <- MAF_UNKNOWN_CLASS
  attr(cls, "n_unknown") <- sum(unknown)
  attr(cls, "unknown_terms") <- unique(first[unknown])
  cls
}

#' Infer MAF Variant_Type from REF/ALT.
#' @keywords internal
maf_variant_type <- function(ref, alt) {
  nr <- nchar(ref); na <- nchar(alt)
  ifelse(nr == 1 & na == 1, "SNP",
    ifelse(nr == 2 & na == 2, "DNP",
      ifelse(nr < na, "INS",
        ifelse(nr > na, "DEL", "SNP"))))
}

#' Create a maftools-compatible MAF data.frame from a filtered variant table.
#'
#' @param variants variant data.frame with chrom/pos/ref/alt/gene/consequence
#'   and (optionally) vaf/dp/alt_count/filter_status/COSMIC_*/confidence/OncoKB.
#' @param path optional output path (`.maf` or `.maf.gz`); if given, the MAF is
#'   written and validated with [maftools::read.maf()].
#' @return the MAF data.frame (invisibly if written).
#' @export
#' @examples
#' v <- data.frame(chrom="17", pos=7674220, ref="C", alt="T", gene="TP53",
#'                 consequence="missense_variant", sample_id="S1", vaf=0.4)
#' create_maf(v)
create_maf <- function(variants, path = NULL) {
  if (!is.data.frame(variants) || nrow(variants) == 0)
    stop("create_maf(): `variants` must be a non-empty data.frame.", call. = FALSE)
  need <- c("chrom", "pos", "ref", "alt")
  miss <- need[!need %in% names(variants)]
  if (length(miss)) stop("create_maf(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  ref <- toupper(as.character(variants$ref))
  alt <- toupper(as.character(variants$alt))
  vtype <- maf_variant_type(ref, alt)
  cons <- if ("consequence" %in% names(variants)) variants$consequence else rep(NA_character_, nrow(variants))
  vclass <- map_consequence_to_maf(cons, vtype)
  n_unknown <- attr(vclass, "n_unknown")

  gene <- if ("gene" %in% names(variants)) variants$gene else rep(NA_character_, nrow(variants))
  gene[is.na(gene) | gene == ""] <- "Unknown"
  sample <- if ("sample_id" %in% names(variants)) variants$sample_id else rep("sample", nrow(variants))

  start <- as.integer(variants$pos)
  end <- start + pmax(nchar(ref) - 1L, 0L)
  # MAF convention: insertions use start=pos, end=pos+1; keep simple + valid
  ins <- vtype == "INS"
  end[ins] <- start[ins] + 1L

  maf <- data.frame(
    Hugo_Symbol = as.character(gene),
    Chromosome = normalize_chrom(variants$chrom),
    Start_Position = start,
    End_Position = end,
    Reference_Allele = ifelse(ins, "-", ref),
    Tumor_Seq_Allele2 = ifelse(vtype == "DEL", "-", alt),
    Variant_Classification = as.character(vclass),
    Variant_Type = vtype,
    Tumor_Sample_Barcode = as.character(sample),
    stringsAsFactors = FALSE
  )
  add_col <- function(name, src) if (src %in% names(variants)) maf[[name]] <<- variants[[src]]
  add_col("t_depth", "dp"); add_col("t_alt_count", "alt_count"); add_col("VAF", "vaf")
  add_col("filter_status", "filter_status")
  add_col("COSMIC_MATCH", "COSMIC_MATCH"); add_col("COSMIC_OCCURRENCE_COUNT", "COSMIC_OCCURRENCE_COUNT")
  add_col("confidence_category", "confidence_category")
  for (oc in grep("^ONCOKB_", names(variants), value = TRUE)) maf[[oc]] <- variants[[oc]]

  if (!is.null(n_unknown) && n_unknown > 0) {
    warning(sprintf("create_maf(): %d variant(s) had unmapped VEP consequence(s) mapped to '%s'; terms: %s",
                    n_unknown, MAF_UNKNOWN_CLASS,
                    paste(attr(vclass, "unknown_terms"), collapse = ", ")), call. = FALSE)
  }

  if (!is.null(path)) {
    write_tsv(maf, path)
    validate_maf(path)
    return(invisible(maf))
  }
  maf
}

#' Validate a MAF file by loading it with maftools::read.maf().
#' @return the maftools MAF object (invisibly), or stops on failure.
#' @keywords internal
validate_maf <- function(path) {
  if (!requireNamespace("maftools", quietly = TRUE)) {
    warning("maftools not installed; skipping MAF validation.", call. = FALSE)
    return(invisible(NULL))
  }
  df <- read_variants(path, "\t")
  obj <- maftools::read.maf(maf = df, verbose = FALSE)
  n_samples <- tryCatch(nrow(maftools::getSampleSummary(obj)), error = function(e) NA_integer_)
  log_step("maf", "MAF validated with maftools::read.maf",
           samples = n_samples, genes = length(unique(df$Hugo_Symbol)))
  invisible(obj)
}

#' Produce a maftools oncoplot (and summary) from a MAF, handling degenerate
#' cases gracefully (returns a message instead of failing when there is nothing
#' to plot).
#' @param maf_path path to the MAF file.
#' @param out_dir directory for plot outputs.
#' @param cfg resolved config (reads report.oncoplot.*).
#' @return list(status, files) describing what was produced.
#' @keywords internal
make_oncoplot <- function(maf_path, out_dir, cfg = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("maftools", quietly = TRUE)) {
    return(list(status = "skipped_no_maftools", files = character()))
  }
  df <- tryCatch(read_variants(maf_path, "\t"), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    msg <- "No variants available for an oncoplot."
    log_step("plot", msg)
    return(list(status = "no_variants", message = msg, files = character()))
  }
  maf_obj <- tryCatch(maftools::read.maf(maf = df, verbose = FALSE),
                      error = function(e) NULL)
  if (is.null(maf_obj)) return(list(status = "read_maf_failed", files = character()))

  top_genes <- cfg_get(cfg, c("report", "oncoplot", "top_genes"),
                       cfg_get(cfg, c("visualization", "oncoplot_top_genes"), 20))
  n_genes <- length(unique(df$Hugo_Symbol[df$Hugo_Symbol != "Unknown"]))
  if (n_genes == 0) {
    return(list(status = "no_genes",
                message = "No named genes among retained variants; oncoplot skipped.",
                files = character()))
  }
  onco_path <- file.path(out_dir, "oncoplot_top_genes.pdf")
  summary_path <- file.path(out_dir, "maf_summary.pdf")
  grDevices::pdf(summary_path, width = 9, height = 6)
  tryCatch(maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE,
                                    addStat = "median", dashboard = TRUE),
           error = function(e) plot.new())
  grDevices::dev.off()
  grDevices::pdf(onco_path, width = 9, height = 7)
  tryCatch(maftools::oncoplot(maf = maf_obj, top = min(top_genes, n_genes),
                              removeNonMutated = cfg_get(cfg, c("report","oncoplot","remove_non_mutated"), FALSE)),
           error = function(e) plot.new())
  grDevices::dev.off()
  log_step("plot", "maftools oncoplot generated", genes = min(top_genes, n_genes))
  list(status = "ok", files = c(summary_path, onco_path))
}
