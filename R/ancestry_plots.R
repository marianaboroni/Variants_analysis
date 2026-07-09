# Informative ancestry plots (ggplot2 is appropriate here: QC/ancestry panels, not
# maftools mutation plots). Each plot is generated only with sufficient data and is
# saved as PNG + PDF; inapplicable plots are skipped with a recorded reason and
# never halt the pipeline. A manifest lists every attempted plot.

#' Generate the ancestry plot set for a run.
#'
#' @param anc result of [infer_ancestry()].
#' @param out_dir output directory (a `ancestry/` subdir is created).
#' @param cfg resolved config (reads `ancestry:plots`).
#' @return a data.frame manifest (plot_id, status, reason, files).
#' @keywords internal
create_ancestry_plots <- function(anc, out_dir, cfg = NULL) {
  dir <- file.path(out_dir, "ancestry"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  man <- list()
  add <- function(id, status, files = character(), reason = NA_character_)
    man[[length(man) + 1]] <<- data.frame(plot_id = id, status = status,
      files = paste(files, collapse = ";"), reason = reason %||% NA_character_,
      generated_at = format(Sys.time()), stringsAsFactors = FALSE)
  has_gg <- requireNamespace("ggplot2", quietly = TRUE)
  if (!has_gg) { add("all", "failed", reason = "ggplot2 not installed"); return(do.call(rbind, man)) }

  save_pair <- function(id, p, w = 8, h = 5) {
    f <- character()
    for (ext in c("png", "pdf")) {
      path <- file.path(dir, paste0(id, ".", ext))
      okp <- tryCatch({ ggplot2::ggsave(path, p, width = w, height = h, dpi = 100); TRUE },
                      error = function(e) FALSE)
      if (okp && file.exists(path) && file.info(path)$size > 0) f <- c(f, path)
    }
    f
  }
  summ <- anc$summary; qc <- anc$qc; snps <- anc$snps
  evaluable <- !is.null(summ) && any(summ$ANCESTRY_INFERENCE_STATUS == "evaluated")

  # 1. proportions (stacked)
  if (evaluable) {
    long <- do.call(rbind, lapply(seq_len(nrow(summ)), function(i) data.frame(
      sample_id = summ$sample_id[i],
      component = ANCESTRY_COMPONENTS,
      proportion = as.numeric(summ[i, paste0("ANCESTRY_", ANCESTRY_COMPONENTS, "_PROPORTION")]),
      stringsAsFactors = FALSE)))
    long <- long[!is.na(long$proportion), ]
    if (nrow(long) > 0) {
      p <- ggplot2::ggplot(long, ggplot2::aes(sample_id, proportion, fill = component)) +
        ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::labs(title = "Inferred genetic-ancestry proportions (experimental)",
          subtitle = "Genetic ancestry is not race or ethnic identity",
          x = NULL, y = "proportion") + ggplot2::theme_minimal()
      add("ancestry_proportions", "generated", save_pair("ancestry_proportions", p))
    } else add("ancestry_proportions", "insufficient_data", reason = "no proportions")
  } else add("ancestry_proportions", "not_applicable", reason = "no evaluable sample")

  # 2. SNP QC funnel + 3. exclusion reasons + 5. by chromosome + 7. VAF dist
  if (!is.null(snps) && nrow(snps) > 0) {
    funnel <- data.frame(stage = c("candidate", "in_panel", "usable"),
      n = c(nrow(snps), sum(snps$in_reference_panel %in% c(TRUE, "TRUE"), na.rm = TRUE),
            sum(snps$ANCESTRY_LOCUS_STATUS == "usable")))
    funnel$stage <- factor(funnel$stage, levels = funnel$stage)
    p2 <- ggplot2::ggplot(funnel, ggplot2::aes(stage, n)) + ggplot2::geom_col(fill = "#2c7fb8") +
      ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.3) +
      ggplot2::labs(title = "Ancestry SNP QC funnel", x = NULL, y = "SNPs") + ggplot2::theme_minimal()
    add("ancestry_snp_qc_funnel", "generated", save_pair("ancestry_snp_qc_funnel", p2))

    ex <- as.data.frame(table(snps$ANCESTRY_SNP_EXCLUSION_REASON))
    if (nrow(ex) > 0) {
      names(ex) <- c("reason", "n"); ex <- ex[order(-ex$n), ]
      p3 <- ggplot2::ggplot(ex, ggplot2::aes(stats::reorder(reason, n), n)) +
        ggplot2::geom_col(fill = "#d95f0e") + ggplot2::coord_flip() +
        ggplot2::labs(title = "Ancestry SNP exclusion reasons", x = NULL, y = "SNPs") + ggplot2::theme_minimal()
      add("ancestry_snp_exclusion_reasons", "generated", save_pair("ancestry_snp_exclusion_reasons", p3))
    } else add("ancestry_snp_exclusion_reasons", "insufficient_data")

    used <- snps[snps$ANCESTRY_LOCUS_STATUS == "usable", , drop = FALSE]
    if (nrow(used) > 0) {
      bychr <- as.data.frame(table(factor(used$chrom, levels = as.character(1:22))))
      names(bychr) <- c("chrom", "n")
      p5 <- ggplot2::ggplot(bychr, ggplot2::aes(chrom, n)) + ggplot2::geom_col(fill = "#31a354") +
        ggplot2::labs(title = "Usable ancestry SNPs by chromosome", x = "chromosome", y = "SNPs") +
        ggplot2::theme_minimal()
      add("ancestry_snps_by_chromosome", "generated", save_pair("ancestry_snps_by_chromosome", p5))

      p7 <- ggplot2::ggplot(used[!is.na(used$vaf), ], ggplot2::aes(vaf)) +
        ggplot2::geom_histogram(bins = 30, fill = "#756bb1") +
        ggplot2::geom_vline(xintercept = c(0, 0.5, 1), linetype = "dashed", color = "grey40") +
        ggplot2::labs(title = "VAF distribution of usable ancestry SNPs",
          subtitle = "Peaks near 0 / 0.5 / 1 expected for germline genotypes", x = "VAF", y = "count") +
        ggplot2::theme_minimal()
      add("ancestry_snp_vaf_distribution", "generated", save_pair("ancestry_snp_vaf_distribution", p7))
    }
  } else {
    for (id in c("ancestry_snp_qc_funnel", "ancestry_snp_exclusion_reasons",
                 "ancestry_snps_by_chromosome", "ancestry_snp_vaf_distribution"))
      add(id, "insufficient_data", reason = "no candidate SNPs")
  }

  # 4. call rate
  if (!is.null(summ) && nrow(summ) > 0) {
    cr <- summ[, c("sample_id", "ANCESTRY_CALL_RATE", "ANCESTRY_CONFIDENCE")]
    p4 <- ggplot2::ggplot(cr, ggplot2::aes(sample_id, ANCESTRY_CALL_RATE, fill = ANCESTRY_CONFIDENCE)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::labs(title = "Ancestry call rate by sample", x = NULL, y = "call rate") +
      ggplot2::theme_minimal()
    add("ancestry_call_rate", "generated", save_pair("ancestry_call_rate", p4))
  }

  # 6. PCA projection
  proj <- anc$projection
  if (!is.null(proj) && any(!is.na(proj$PC1)) && any(!is.na(proj$PC2))) {
    p6 <- ggplot2::ggplot(proj[!is.na(proj$PC1), ], ggplot2::aes(PC1, PC2)) +
      ggplot2::geom_point(color = "#e6550d", size = 3) +
      ggplot2::labs(title = "PCA projection onto reference populations (experimental)",
        subtitle = "Represents genetic similarity to reference panel; not race/ethnicity") +
      ggplot2::theme_minimal()
    add("ancestry_pca_projection", "generated", save_pair("ancestry_pca_projection", p6))
  } else add("ancestry_pca_projection", "insufficient_data", reason = "no PCA coordinates")

  mani <- do.call(rbind, man)
  write_tsv(mani, file.path(dir, "plot_manifest.tsv"))
  mani
}
