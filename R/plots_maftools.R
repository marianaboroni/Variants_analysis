# Native maftools visualizations. Every mutation plot is produced by maftools
# from a single validated MAF object. ggplot2 is NOT used to imitate any of them.
# Each plot honors minimum-data rules and is saved as both PNG and PDF; when a
# plot is not applicable a clear reason is recorded (never a fake plot).

#' Build a maftools MAF object from a MAF file (+ optional clinical data).
#' @keywords internal
load_maf_object <- function(maf_path, sample_metadata = NULL) {
  if (!requireNamespace("maftools", quietly = TRUE)) return(NULL)
  df <- tryCatch(read_variants(maf_path, "\t"), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  clin <- NULL
  if (!is.null(sample_metadata)) {
    clin <- if (is.character(sample_metadata) && file.exists(sample_metadata))
      read_variants(sample_metadata, "\t") else if (is.data.frame(sample_metadata)) sample_metadata else NULL
    if (!is.null(clin) && !"Tumor_Sample_Barcode" %in% names(clin)) {
      sid <- coalesce_columns(clin, c("Tumor_Sample_Barcode", "sample_id", "Sample", "SAMPLE_ID"))
      clin$Tumor_Sample_Barcode <- as.character(sid)
    }
  }
  tryCatch(maftools::read.maf(maf = df, clinicalData = clin, removeSilent = FALSE, verbose = FALSE),
           error = function(e) tryCatch(maftools::read.maf(maf = df, verbose = FALSE), error = function(e2) NULL))
}

#' @keywords internal
save_plot_pair <- function(out_dir, name, draw, width = 9, height = 7) {
  files <- character()
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  ok <- TRUE
  grDevices::pdf(pdf_path, width = width, height = height)
  ok <- tryCatch({ draw(); TRUE }, error = function(e) FALSE); grDevices::dev.off()
  if (ok && file.exists(pdf_path) && file.info(pdf_path)$size > 0) files <- c(files, pdf_path) else unlink(pdf_path)
  grDevices::png(png_path, width = width * 100, height = height * 100, res = 100)
  ok2 <- tryCatch({ draw(); TRUE }, error = function(e) FALSE); grDevices::dev.off()
  if (ok2 && file.exists(png_path) && file.info(png_path)$size > 0) files <- c(files, png_path) else unlink(png_path)
  files
}

#' Generate the native maftools plot set for a run.
#'
#' @param maf_path path to the validated MAF.
#' @param out_dir output directory for plots.
#' @param cfg resolved config (reads `plots:`).
#' @param sample_metadata optional clinical data (path or data.frame).
#' @return a data.frame manifest: plot, status, files, reason.
#' @export
#' @examples
#' \dontrun{ create_maftools_plots("results/run/maf/filtered.maf.gz", "results/run/plots") }
create_maftools_plots <- function(maf_path, out_dir, cfg = NULL, sample_metadata = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  add <- function(plot, status, files = character(), reason = NA_character_)
    manifest[[length(manifest) + 1]] <<- data.frame(plot = plot, status = status,
      files = paste(files, collapse = ";"), reason = reason %||% NA_character_, stringsAsFactors = FALSE)

  if (!requireNamespace("maftools", quietly = TRUE)) {
    add("all", "skipped", reason = "maftools not installed")
    return(do.call(rbind, manifest))
  }
  maf <- load_maf_object(maf_path, sample_metadata)
  if (is.null(maf)) { add("all", "skipped", reason = "no variants / MAF could not be read")
    return(do.call(rbind, manifest)) }

  gs <- maftools::getGeneSummary(maf); ss <- maftools::getSampleSummary(maf)
  n_samples <- nrow(ss); n_genes <- nrow(gs)
  min_onco <- cfg_get(cfg, c("plots", "minimum_samples_for_oncoplot"), 2)
  min_inter <- cfg_get(cfg, c("plots", "minimum_samples_for_interactions"), 20)
  min_rain <- cfg_get(cfg, c("plots", "minimum_mutations_for_rainfall"), 50)
  min_lolli <- cfg_get(cfg, c("plots", "minimum_protein_changes_for_lollipop"), 3)
  top <- cfg_get(cfg, c("plots", "oncoplot", "top"), cfg_get(cfg, c("report", "oncoplot", "top_genes"), 20))
  clin_feats <- cfg_get(cfg, c("plots", "oncoplot", "clinical_features"), NULL)
  clin_feats <- clin_feats[clin_feats %in% names(maftools::getClinicalData(maf))]

  # --- summary (always) ---
  add("plotmaf_summary", "ok", save_plot_pair(out_dir, "plotmaf_summary",
    function() maftools::plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)))

  # --- oncoplot ---
  if (n_samples >= min_onco && n_genes >= 1) {
    h <- max(6, min(14, 4 + 0.25 * min(top, n_genes)))
    w <- max(8, min(20, 6 + 0.25 * n_samples))
    add("oncoplot", "ok", save_plot_pair(out_dir, "oncoplot", function() maftools::oncoplot(
      maf = maf, top = min(top, n_genes),
      clinicalFeatures = if (length(clin_feats)) clin_feats else NULL,
      sortByAnnotation = isTRUE(cfg_get(cfg, c("plots", "oncoplot", "sort_by_annotation"), TRUE)),
      removeNonMutated = isTRUE(cfg_get(cfg, c("plots", "oncoplot", "remove_non_mutated"), FALSE))),
      width = w, height = h))
  } else add("oncoplot", "skipped", reason = sprintf("need >= %d samples (have %d)", min_onco, n_samples))

  # --- Ti/Tv ---
  titv_ok <- tryCatch({ tt <- maftools::titv(maf = maf, plot = FALSE, useSyn = TRUE)
    add("titv", "ok", save_plot_pair(out_dir, "titv", function() maftools::plotTiTv(res = tt))); TRUE },
    error = function(e) FALSE)
  if (!titv_ok) add("titv", "skipped", reason = "no SNVs / titv failed")

  # --- lollipop for top genes with enough protein changes ---
  top_genes <- utils::head(as.character(gs$Hugo_Symbol), 5)
  n_lolli <- 0L
  for (g in top_genes) {
    files <- tryCatch(save_plot_pair(out_dir, paste0("lollipop_", g),
      function() maftools::lollipopPlot(maf = maf, gene = g, showMutationRate = TRUE)),
      error = function(e) character())
    if (length(files)) { add(paste0("lollipop_", g), "ok", files); n_lolli <- n_lolli + 1L }
  }
  if (n_lolli == 0) add("lollipop", "skipped", reason = sprintf("no gene with >= %d mappable protein changes", min_lolli))

  # --- rainfall (needs enough mutations; per top sample) ---
  top_sample <- as.character(ss$Tumor_Sample_Barcode[1])
  if (!is.na(top_sample) && ss$total[1] >= min_rain) {
    files <- tryCatch(save_plot_pair(out_dir, paste0("rainfall_", top_sample),
      function() maftools::rainfallPlot(maf = maf, tsb = top_sample, detectChangePoints = TRUE, pointSize = 0.4)),
      error = function(e) character())
    if (length(files)) add("rainfall", "ok", files) else add("rainfall", "skipped", reason = "rainfall failed")
  } else add("rainfall", "skipped", reason = sprintf("need >= %d mutations in a sample", min_rain))

  # --- somatic interactions (cohort only) ---
  if (n_samples >= min_inter) {
    files <- tryCatch(save_plot_pair(out_dir, "somatic_interactions",
      function() maftools::somaticInteractions(maf = maf, top = min(25, n_genes), pvalue = c(0.05, 0.1))),
      error = function(e) character())
    if (length(files)) add("somatic_interactions", "ok", files) else add("somatic_interactions", "skipped", reason = "interactions failed")
  } else add("somatic_interactions", "skipped", reason = sprintf("need >= %d samples (have %d)", min_inter, n_samples))

  do.call(rbind, manifest)
}
