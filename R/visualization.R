make_visualization_outputs <- function(
    variants,
    cfg,
    sample_qc = NULL,
    tmb_summary = NULL,
    variant_qc_summary = NULL,
    clonality_summary = NULL) {
  enabled <- isTRUE(cfg_get(cfg, c("visualization", "enabled"), TRUE))
  outdir <- cfg_get(cfg, c("visualization", "dir"), file.path(cfg$output$dir, "figures"))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  manifest <- empty_figure_manifest()
  if (!enabled) {
    manifest <- add_figure_manifest(manifest, NA_character_, "visualization", "skipped", "visualization.enabled is false")
    return(manifest)
  }

  somatic_variants <- visualization_somatic_variants(variants, cfg)
  maf_table <- build_maf_export(somatic_variants)
  maf_path <- file.path(outdir, "somatic_maftools_input.maf")
  write_tsv(maf_table, maf_path)
  manifest <- add_figure_manifest(
    manifest,
    maf_path,
    "maftools_input",
    "maftools_input",
    "written",
    paste(nrow(maf_table), "somatic variants exported")
  )

  manifest <- run_maftools_plots(maf_path, maf_table, cfg, outdir, manifest)
  manifest <- run_fallback_oncoplot(maf_table, cfg, outdir, manifest)
  manifest <- plot_variant_class_stack(variants, cfg, outdir, manifest)
  manifest <- plot_tmb_bar(tmb_summary, cfg, outdir, manifest)
  manifest <- plot_driver_summary(variants, cfg, outdir, manifest)
  manifest <- plot_qc_depth_vaf(variants, cfg, outdir, manifest)
  manifest <- plot_recurrence_population(variants, cfg, outdir, manifest)
  manifest <- plot_clonality_outputs(variants, clonality_summary, cfg, outdir, manifest)

  manifest
}

empty_figure_manifest <- function() {
  data.frame(
    file = character(),
    type = character(),
    plot_category = character(),
    status = character(),
    note = character(),
    stringsAsFactors = FALSE
  )
}

add_figure_manifest <- function(manifest, file, type, plot_category = NA_character_, status, note = "") {
  rbind(
    manifest,
    data.frame(
      file = as.character(file),
      type = as.character(type),
      plot_category = as.character(plot_category),
      status = as.character(status),
      note = as.character(note),
      stringsAsFactors = FALSE
    )
  )
}

visualization_somatic_variants <- function(x, cfg) {
  include_classes <- cfg_get(
    cfg,
    c("visualization", "include_final_classes"),
    c("high_confidence_somatic", "probable_somatic")
  )
  keep <- x$final_class %in% include_classes
  use_ml <- isTRUE(cfg_get(cfg, c("visualization", "use_ml_filter"), FALSE))
  if (use_ml && "ml_true_positive_probability" %in% names(x)) {
    cutoff <- cfg_get(cfg, c("visualization", "ml_true_positive_cutoff"), 0.60)
    keep <- keep & !is.na(x$ml_true_positive_probability) & x$ml_true_positive_probability >= cutoff
  }
  x[keep, , drop = FALSE]
}

build_maf_export <- function(x) {
  if (nrow(x) == 0) {
    return(data.frame(
      Hugo_Symbol = character(),
      Chromosome = character(),
      Start_Position = numeric(),
      End_Position = numeric(),
      Reference_Allele = character(),
      Tumor_Seq_Allele2 = character(),
      Variant_Classification = character(),
      Variant_Type = character(),
      Tumor_Sample_Barcode = character(),
      t_depth = numeric(),
      t_ref_count = numeric(),
      t_alt_count = numeric(),
      tumor_f = numeric(),
      HGVSp_Short = character(),
      tumor_type = character(),
      final_class = character(),
      driver_class = character(),
      stringsAsFactors = FALSE
    ))
  }

  end_pos <- x$pos + pmax(nchar(as.character(x$ref)), 1) - 1
  end_pos[is.na(end_pos)] <- x$pos[is.na(end_pos)]

  data.frame(
    Hugo_Symbol = clean_maf_text(x$gene, "Unknown"),
    Chromosome = gsub("^chr", "", as.character(x$chrom), ignore.case = TRUE),
    Start_Position = as.integer(x$pos),
    End_Position = as.integer(end_pos),
    Reference_Allele = clean_maf_text(x$ref, "-"),
    Tumor_Seq_Allele2 = clean_maf_text(x$alt, "-"),
    Variant_Classification = map_consequence_to_maf_class(x$consequence, x$variant_type),
    Variant_Type = map_variant_type_to_maf(x$variant_type),
    Tumor_Sample_Barcode = clean_maf_text(x$sample_id, "Unknown"),
    t_depth = as.integer(round(x$dp)),
    t_ref_count = as.integer(round(x$ref_count)),
    t_alt_count = as.integer(round(x$alt_count)),
    tumor_f = x$vaf,
    HGVSp_Short = clean_maf_text(x$protein_change, ""),
    tumor_type = clean_maf_text(x$tumor_type, ""),
    final_class = clean_maf_text(x$final_class, ""),
    driver_class = clean_maf_text(if ("driver_class" %in% names(x)) x$driver_class else NA, ""),
    stringsAsFactors = FALSE
  )
}

clean_maf_text <- function(x, default = "") {
  out <- as.character(x)
  out[is.na(out) | out == "NA" | trimws(out) == ""] <- default
  out
}

map_variant_type_to_maf <- function(variant_type) {
  z <- toupper(as.character(variant_type))
  ifelse(z == "SNV", "SNP",
         ifelse(z == "INS", "INS",
                ifelse(z == "DEL", "DEL", "ONP")))
}

map_consequence_to_maf_class <- function(consequence, variant_type = NULL) {
  z <- tolower(as.character(consequence))
  vt <- toupper(as.character(variant_type))
  out <- rep("Targeted_Region", length(z))
  out[grepl("missense", z)] <- "Missense_Mutation"
  out[grepl("synonymous|stop_retained|start_retained", z)] <- "Silent"
  out[grepl("stop_gained|nonsense", z)] <- "Nonsense_Mutation"
  out[grepl("stop_lost", z)] <- "Nonstop_Mutation"
  out[grepl("start_lost", z)] <- "Translation_Start_Site"
  out[grepl("splice_acceptor|splice_donor|splice_site", z)] <- "Splice_Site"
  out[grepl("splice_region", z)] <- "Splice_Region"
  out[grepl("inframe", z) & vt == "DEL"] <- "In_Frame_Del"
  out[grepl("inframe", z) & vt == "INS"] <- "In_Frame_Ins"
  out[grepl("inframe", z) & !(vt %in% c("DEL", "INS"))] <- "In_Frame_Ins"
  out[grepl("frameshift", z) & vt == "DEL"] <- "Frame_Shift_Del"
  out[grepl("frameshift", z) & vt == "INS"] <- "Frame_Shift_Ins"
  out[grepl("frameshift", z) & !(vt %in% c("DEL", "INS"))] <- "Frame_Shift_Del"
  out
}

run_maftools_plots <- function(maf_path, maf_table, cfg, outdir, manifest) {
  top_n <- cfg_get(cfg, c("visualization", "oncoplot_top_genes"), 20)
  if (nrow(maf_table) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "maftools", "maftools", "skipped", "no somatic variants for maftools"))
  }
  if (!requireNamespace("maftools", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "maftools", "maftools", "skipped", "R package maftools is not installed"))
  }

  result <- tryCatch({
    maf_obj <- maftools::read.maf(maf = maf_path, verbose = FALSE)

    summary_pdf <- file.path(outdir, "maftools_summary.pdf")
    grDevices::pdf(summary_pdf, width = 10, height = 8)
    maftools::plotmafSummary(maf = maf_obj, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
    grDevices::dev.off()

    oncoplot_pdf <- file.path(outdir, "maftools_oncoplot_top_genes.pdf")
    grDevices::pdf(oncoplot_pdf, width = 12, height = 8)
    maftools::oncoplot(maf = maf_obj, top = top_n, removeNonMutated = TRUE, draw_titv = FALSE)
    grDevices::dev.off()

    list(summary_pdf = summary_pdf, oncoplot_pdf = oncoplot_pdf)
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    structure(list(message = conditionMessage(e)), class = "maftools_error")
  })

  if (inherits(result, "maftools_error")) {
    return(add_figure_manifest(manifest, NA_character_, "maftools", "maftools", "failed", result$message))
  }
  manifest <- add_figure_manifest(manifest, result$summary_pdf, "maftools_summary", "maftools", "written", "")
  add_figure_manifest(manifest, result$oncoplot_pdf, "maftools_oncoplot", "maftools", "written", paste("top", top_n, "genes"))
}

run_fallback_oncoplot <- function(maf_table, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "oncoplot_fallback", "oncoplot", "skipped", "R package ggplot2 is not installed"))
  }
  if (nrow(maf_table) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "oncoplot_fallback", "oncoplot", "skipped", "no somatic variants"))
  }

  top_n <- cfg_get(cfg, c("visualization", "oncoplot_top_genes"), 20)
  dt <- data.table::as.data.table(maf_table)
  gene_counts <- dt[Hugo_Symbol != "Unknown", .N, by = Hugo_Symbol][order(-N)]
  if (nrow(gene_counts) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "oncoplot_fallback", "oncoplot", "skipped", "no gene symbols available"))
  }
  top_genes <- gene_counts$Hugo_Symbol[seq_len(min(top_n, nrow(gene_counts)))]
  plot_dt <- dt[Hugo_Symbol %in% top_genes]
  plot_dt[, priority := maf_class_priority(Variant_Classification)]
  plot_dt <- plot_dt[order(priority), .(
    Variant_Classification = ifelse(.N > 1, "Multi_Hit", Variant_Classification[1]),
    n_hits = .N
  ), by = .(Tumor_Sample_Barcode, Hugo_Symbol)]

  sample_order <- dt[, .N, by = Tumor_Sample_Barcode][order(-N)]$Tumor_Sample_Barcode
  gene_order <- rev(top_genes)
  plot_dt$Tumor_Sample_Barcode <- factor(plot_dt$Tumor_Sample_Barcode, levels = sample_order)
  plot_dt$Hugo_Symbol <- factor(plot_dt$Hugo_Symbol, levels = gene_order)

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Variant_Classification), color = "white", size = 0.25) +
    ggplot2::scale_fill_brewer(palette = "Set2", na.value = "grey70") +
    ggplot2::labs(x = "Sample", y = "Gene", fill = "Variant class") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )

  path_png <- file.path(outdir, "oncoplot_fallback_top_genes.png")
  path_pdf <- file.path(outdir, "oncoplot_fallback_top_genes.pdf")
  save_ggplot_pair(p, path_png, path_pdf, width = 12, height = 7)
  manifest <- add_figure_manifest(manifest, path_png, "oncoplot_fallback_png", "oncoplot", "written", paste("top", length(top_genes), "genes"))
  add_figure_manifest(manifest, path_pdf, "oncoplot_fallback_pdf", "oncoplot", "written", paste("top", length(top_genes), "genes"))
}

maf_class_priority <- function(x) {
  priority <- c(
    Nonsense_Mutation = 1,
    Frame_Shift_Del = 2,
    Frame_Shift_Ins = 3,
    Splice_Site = 4,
    Translation_Start_Site = 5,
    Missense_Mutation = 6,
    In_Frame_Del = 7,
    In_Frame_Ins = 8,
    Nonstop_Mutation = 9,
    Splice_Region = 10,
    Silent = 11,
    Targeted_Region = 12
  )
  out <- priority[as.character(x)]
  out[is.na(out)] <- 99
  as.numeric(out)
}

plot_variant_class_stack <- function(variants, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "classification_stack", "classification", "skipped", "R package ggplot2 is not installed"))
  }
  dt <- data.table::as.data.table(variants)
  if (nrow(dt) == 0 || !"final_class" %in% names(dt)) {
    return(add_figure_manifest(manifest, NA_character_, "classification_stack", "classification", "skipped", "no final_class data"))
  }
  plot_dt <- dt[, .N, by = .(sample_id, final_class)]
  sample_order <- dt[, .N, by = sample_id][order(-N)]$sample_id
  plot_dt$sample_id <- factor(plot_dt$sample_id, levels = sample_order)

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = sample_id, y = N, fill = final_class)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Sample", y = "Variant count (log10)", fill = "Final class") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  path_png <- file.path(outdir, "variant_classification_stack.png")
  path_pdf <- file.path(outdir, "variant_classification_stack.pdf")
  save_ggplot_pair(p, path_png, path_pdf, width = 12, height = 6)
  manifest <- add_figure_manifest(manifest, path_png, "classification_stack_png", "classification", "written", "")
  add_figure_manifest(manifest, path_pdf, "classification_stack_pdf", "classification", "written", "")
}

plot_tmb_bar <- function(tmb_summary, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "tmb_bar", "tmb", "skipped", "R package ggplot2 is not installed"))
  }
  if (is.null(tmb_summary) || nrow(tmb_summary) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "tmb_bar", "tmb", "skipped", "no TMB summary"))
  }
  dt <- data.table::as.data.table(tmb_summary)
  dt <- dt[order(-tmb_mut_per_mb)]
  dt$sample_id <- factor(dt$sample_id, levels = dt$sample_id)
  high <- cfg_get(cfg, c("tmb", "high_threshold"), 10)
  intermediate <- cfg_get(cfg, c("tmb", "intermediate_threshold"), 5)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = sample_id, y = tmb_mut_per_mb, fill = tmb_category)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_hline(yintercept = intermediate, linetype = "dashed", color = "grey45") +
    ggplot2::geom_hline(yintercept = high, linetype = "dotted", color = "grey25") +
    ggplot2::labs(x = "Sample", y = "TMB (mut/Mb)", fill = "TMB category") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  path_png <- file.path(outdir, "tmb_by_sample.png")
  path_pdf <- file.path(outdir, "tmb_by_sample.pdf")
  save_ggplot_pair(p, path_png, path_pdf, width = 12, height = 6)
  manifest <- add_figure_manifest(manifest, path_png, "tmb_bar_png", "tmb", "written", "")
  add_figure_manifest(manifest, path_pdf, "tmb_bar_pdf", "tmb", "written", "")
}

plot_driver_summary <- function(variants, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "driver_summary", "driver", "skipped", "R package ggplot2 is not installed"))
  }
  if (!"driver_class" %in% names(variants)) {
    return(add_figure_manifest(manifest, NA_character_, "driver_summary", "driver", "skipped", "no driver_class data"))
  }
  dt <- data.table::as.data.table(variants)
  plot_dt <- dt[, .N, by = .(sample_id, driver_class)]
  sample_order <- dt[, .N, by = sample_id][order(-N)]$sample_id
  plot_dt$sample_id <- factor(plot_dt$sample_id, levels = sample_order)

  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = sample_id, y = N, fill = driver_class)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Sample", y = "Variant count (log10)", fill = "Driver class") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  path_png <- file.path(outdir, "driver_class_by_sample.png")
  path_pdf <- file.path(outdir, "driver_class_by_sample.pdf")
  save_ggplot_pair(p, path_png, path_pdf, width = 12, height = 6)
  manifest <- add_figure_manifest(manifest, path_png, "driver_summary_png", "driver", "written", "")
  add_figure_manifest(manifest, path_pdf, "driver_summary_pdf", "driver", "written", "")
}

plot_qc_depth_vaf <- function(variants, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "qc_depth_vaf", "qc", "skipped", "R package ggplot2 is not installed"))
  }
  dt <- data.table::as.data.table(sample_variants_for_plot(variants, cfg))
  dt <- dt[!is.na(dp) & !is.na(vaf)]
  if (nrow(dt) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "qc_depth_vaf", "qc", "skipped", "no DP/VAF data"))
  }

  p_vaf <- ggplot2::ggplot(dt, ggplot2::aes(x = vaf, fill = final_class)) +
    ggplot2::geom_histogram(bins = 50, alpha = 0.70, position = "identity") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "VAF", y = "Variant count", fill = "Final class") +
    ggplot2::theme_minimal(base_size = 10)

  p_depth <- ggplot2::ggplot(dt, ggplot2::aes(x = dp, fill = final_class)) +
    ggplot2::geom_histogram(bins = 80, alpha = 0.75, position = "identity") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Depth (log10)", y = "Variant count", fill = "Final class") +
    ggplot2::theme_minimal(base_size = 10)

  vaf_png <- file.path(outdir, "qc_vaf_histogram_by_class.png")
  vaf_pdf <- file.path(outdir, "qc_vaf_histogram_by_class.pdf")
  depth_png <- file.path(outdir, "qc_depth_histogram_by_class.png")
  depth_pdf <- file.path(outdir, "qc_depth_histogram_by_class.pdf")
  save_ggplot_pair(p_vaf, vaf_png, vaf_pdf, width = 9, height = 6)
  save_ggplot_pair(p_depth, depth_png, depth_pdf, width = 9, height = 6)
  manifest <- add_figure_manifest(manifest, vaf_png, "qc_vaf_histogram_png", "qc", "written", "")
  manifest <- add_figure_manifest(manifest, vaf_pdf, "qc_vaf_histogram_pdf", "qc", "written", "")
  manifest <- add_figure_manifest(manifest, depth_png, "qc_depth_histogram_png", "qc", "written", "")
  add_figure_manifest(manifest, depth_pdf, "qc_depth_histogram_pdf", "qc", "written", "")
}

plot_recurrence_population <- function(variants, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "recurrence_population", "recurrence", "skipped", "R package ggplot2 is not installed"))
  }
  dt <- data.table::as.data.table(sample_variants_for_plot(variants, cfg))
  if (!all(c("variant_cohort_freq", "max_pop_af") %in% names(dt))) {
    return(add_figure_manifest(manifest, NA_character_, "recurrence_population", "recurrence", "skipped", "missing recurrence or population AF"))
  }
  dt <- dt[!is.na(variant_cohort_freq) | !is.na(max_pop_af)]
  if (nrow(dt) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "recurrence_population", "recurrence", "skipped", "no recurrence/population AF data"))
  }
  dt$max_pop_af_plot <- pmax(dt$max_pop_af, 1e-6, na.rm = TRUE)
  dt$max_pop_af_plot[is.infinite(dt$max_pop_af_plot)] <- 1e-6

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = max_pop_af_plot, y = variant_cohort_freq, color = final_class)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Max population AF (log10)", y = "Internal cohort recurrence", color = "Final class") +
    ggplot2::theme_minimal(base_size = 10)

  path_png <- file.path(outdir, "recurrence_vs_population_af.png")
  path_pdf <- file.path(outdir, "recurrence_vs_population_af.pdf")
  save_ggplot_pair(p, path_png, path_pdf, width = 9, height = 6)
  manifest <- add_figure_manifest(manifest, path_png, "recurrence_population_png", "recurrence", "written", "")
  add_figure_manifest(manifest, path_pdf, "recurrence_population_pdf", "recurrence", "written", "")
}

plot_clonality_outputs <- function(variants, clonality_summary, cfg, outdir, manifest) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "clonality", "clonality", "skipped", "R package ggplot2 is not installed"))
  }
  if (!"clonality_class" %in% names(variants)) {
    return(add_figure_manifest(manifest, NA_character_, "clonality", "clonality", "skipped", "no clonality_class data"))
  }
  dt <- data.table::as.data.table(variants)
  somatic_classes <- c("high_confidence_somatic", "probable_somatic")
  dt <- dt[final_class %in% somatic_classes]
  if (nrow(dt) == 0) {
    return(add_figure_manifest(manifest, NA_character_, "clonality", "clonality", "skipped", "no somatic variants"))
  }

  count_dt <- dt[, .N, by = .(sample_id, clonality_class)]
  sample_order <- dt[, .N, by = sample_id][order(-N)]$sample_id
  count_dt$sample_id <- factor(count_dt$sample_id, levels = sample_order)

  p_counts <- ggplot2::ggplot(count_dt, ggplot2::aes(x = sample_id, y = N, fill = clonality_class)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::labs(x = "Sample", y = "Somatic variant count", fill = "Clonality") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

  counts_png <- file.path(outdir, "clonality_class_by_sample.png")
  counts_pdf <- file.path(outdir, "clonality_class_by_sample.pdf")
  save_ggplot_pair(p_counts, counts_png, counts_pdf, width = 12, height = 6)
  manifest <- add_figure_manifest(manifest, counts_png, "clonality_counts_png", "clonality", "written", "")
  manifest <- add_figure_manifest(manifest, counts_pdf, "clonality_counts_pdf", "clonality", "written", "")

  plot_dt <- data.table::as.data.table(sample_variants_for_plot(as.data.frame(dt), cfg))
  plot_dt$plot_value <- ifelse(!is.na(plot_dt$ccf_capped), plot_dt$ccf_capped, plot_dt$vaf)
  plot_dt$sample_id <- factor(plot_dt$sample_id, levels = sample_order)
  plot_dt <- plot_dt[!is.na(plot_value)]
  if (nrow(plot_dt) > 0) {
    p_dist <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = sample_id, y = plot_value, color = clonality_class)) +
      ggplot2::geom_jitter(width = 0.18, height = 0, alpha = 0.45, size = 0.9) +
      ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
      ggplot2::labs(x = "Sample", y = "CCF capped at 1 or VAF proxy", color = "Clonality") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

    dist_png <- file.path(outdir, "clonality_ccf_or_vaf_distribution.png")
    dist_pdf <- file.path(outdir, "clonality_ccf_or_vaf_distribution.pdf")
    save_ggplot_pair(p_dist, dist_png, dist_pdf, width = 12, height = 6)
    manifest <- add_figure_manifest(manifest, dist_png, "clonality_distribution_png", "clonality", "written", "")
    manifest <- add_figure_manifest(manifest, dist_pdf, "clonality_distribution_pdf", "clonality", "written", "")
  }

  if (!is.null(clonality_summary) && nrow(clonality_summary) > 0) {
    path <- file.path(outdir, "clonality_summary_for_figures.tsv")
    write_tsv(clonality_summary, path)
    manifest <- add_figure_manifest(manifest, path, "clonality_summary_copy", "written", "")
  }
  manifest
}

sample_variants_for_plot <- function(x, cfg) {
  max_points <- cfg_get(cfg, c("visualization", "max_plot_points"), 200000)
  if (nrow(x) <= max_points) return(x)
  set.seed(cfg_get(cfg, c("visualization", "seed"), 20260621))
  x[sample.int(nrow(x), max_points), , drop = FALSE]
}

save_ggplot_pair <- function(plot, png_path, pdf_path, width, height) {
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 160)
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in")
}
