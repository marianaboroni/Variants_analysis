run_variant_analysis <- function(config_path) {
  cfg <- read_config(config_path)
  outdir <- cfg$output$dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Reading variants: ", cfg$input$variants)
  variants <- read_variants(cfg$input$variants, cfg$input$delimiter)
  variants <- standardize_variant_table(variants, cfg)
  variants <- attach_sample_metadata(variants, cfg)

  message("Computing technical and population features")
  variants <- add_basic_features(variants, cfg)
  hard_filter_result <- add_sample_qc_and_hard_filters(variants, cfg)
  variants <- hard_filter_result$variants
  sample_qc <- hard_filter_result$sample_qc
  recurrence <- compute_cohort_recurrence(variants, cfg)
  variants <- merge_recurrence_features(variants, recurrence)

  message("Scoring variants")
  variants <- score_variants(variants, cfg)

  message("Adding OncoKB/COSMIC validation layers")
  variants <- add_validation_layers(variants, cfg)

  message("Classifying driver/passenger evidence")
  variants <- add_driver_layers(variants, cfg)

  message("Assigning conservative final classes")
  variants <- classify_variants(variants, cfg)

  message("Training semi-supervised ML filter")
  ml_result <- run_ml_filter(variants, cfg)
  variants <- ml_result$variants
  ml_metrics <- ml_result$metrics
  ml_training_labels <- ml_result$training_labels
  active_learning_candidates <- ml_result$active_learning_candidates

  message("Estimating clonality")
  clonality_result <- add_clonality_estimates(variants, cfg)
  variants <- clonality_result$variants
  clonality_summary <- clonality_result$summary

  message("Calculating TMB and QC summaries")
  tmb_summary <- calculate_tmb(variants, cfg)
  variant_qc_summary <- build_variant_call_qc_summary(variants, sample_qc)

  message("Writing outputs")
  write_outputs(
    variants,
    recurrence,
    cfg,
    sample_qc,
    tmb_summary,
    ml_metrics,
    variant_qc_summary,
    clonality_summary,
    ml_training_labels,
    active_learning_candidates
  )

  message("Generating figures")
  figure_manifest <- make_visualization_outputs(
    variants,
    cfg,
    sample_qc = sample_qc,
    tmb_summary = tmb_summary,
    variant_qc_summary = variant_qc_summary,
    clonality_summary = clonality_summary
  )
  write_tsv(figure_manifest, file.path(outdir, "figure_manifest.tsv"))

  message("Done: ", normalizePath(outdir, mustWork = FALSE))
  invisible(list(
    variants = variants,
    recurrence = recurrence,
    sample_qc = sample_qc,
    tmb_summary = tmb_summary,
    ml_metrics = ml_metrics,
    ml_training_labels = ml_training_labels,
    active_learning_candidates = active_learning_candidates,
    variant_qc_summary = variant_qc_summary,
    clonality_summary = clonality_summary,
    figure_manifest = figure_manifest
  ))
}
