write_outputs <- function(
    variants,
    recurrence,
    cfg,
    sample_qc = NULL,
    tmb_summary = NULL,
    ml_metrics = NULL,
    variant_qc_summary = NULL,
    clonality_summary = NULL,
    ml_training_labels = NULL,
    active_learning_candidates = NULL) {
  outdir <- cfg$output$dir
  variants$validation_key <- make_validation_key(variants)

  write_tsv(variants, file.path(outdir, "variants_final.tsv"))

  removed <- variants[variants$final_class %in% c("likely_germline", "likely_artifact", "technical_fail", "manual_review_required"), ]
  removed$primary_reason <- as.character(removed$primary_reason)
  removed$final_class <- as.character(removed$final_class)
  write_tsv(removed, file.path(outdir, "removed_variants_summary.tsv"))

  if (!is.null(active_learning_candidates)) {
    write_tsv(active_learning_candidates, file.path(outdir, "variants_for_review.tsv"))
  }

  write_tsv(class_summary(variants), file.path(outdir, "classification_summary.tsv"))
  write_tsv(validation_summary(variants), file.path(outdir, "oncokb_cosmic_validation_summary.tsv"))
  if (!is.null(sample_qc)) {
    write_tsv(sample_qc, file.path(outdir, "sample_qc_summary.tsv"))
  }
  if (!is.null(tmb_summary)) {
    write_tsv(tmb_summary, file.path(outdir, "tmb_summary.tsv"))
  }
  if (!is.null(ml_metrics)) {
    write_tsv(ml_metrics, file.path(outdir, "ml_model_report.tsv"))
  }
  if (!is.null(variant_qc_summary)) {
    write_tsv(variant_qc_summary, file.path(outdir, "sample_variant_qc_summary.tsv"))
  }
  if (!is.null(clonality_summary)) {
    write_tsv(clonality_summary, file.path(outdir, "clonality_summary.tsv"))
  }
  if (!is.null(ml_training_labels)) {
    write_tsv(ml_training_labels, file.path(outdir, "ml_training_table.tsv"))
  }

  if (isTRUE(cfg_get(cfg, c("output", "generate_filter_report"), FALSE))) {
    write_filter_report(variants, sample_qc, tmb_summary, outdir)
  }
}

write_filter_report <- function(variants, sample_qc, tmb_summary, outdir) {
  report_path <- file.path(outdir, "filter_report.md")
  counts <- as.list(table(factor(variants$final_class, levels = c(
    "high_confidence_somatic", "probable_somatic", "uncertain_tumor_only",
    "likely_germline", "likely_artifact", "technical_fail", "manual_review_required"
  ))))
  names(counts) <- c(
    "high_confidence_somatic", "probable_somatic", "uncertain_tumor_only",
    "likely_germline", "likely_artifact", "technical_fail", "manual_review_required"
  )

  populational_removals <- sum(variants$final_class %in% c("likely_germline"), na.rm = TRUE)
  artifact_removals <- sum(variants$final_class %in% c("likely_artifact", "technical_fail"), na.rm = TRUE)
  review_required <- sum(variants$final_class == "manual_review_required", na.rm = TRUE)
  high_confidence <- counts$high_confidence_somatic
  probable_somatic <- counts$probable_somatic

  population_sources <- sort(table(variants$population_category), decreasing = TRUE)
  artifact_sources <- sort(table(variants$artifact_category), decreasing = TRUE)
  filter_reasons <- sort(table(variants$primary_reason), decreasing = TRUE)

  lines <- c(
    "# Filter Report",
    "",
    "## Resumo geral",
    sprintf("- Número inicial de variantes: %d", nrow(variants)),
    sprintf("- Número de high_confidence_somatic: %d", high_confidence),
    sprintf("- Número de probable_somatic: %d", probable_somatic),
    sprintf("- Número de likely_germline: %d", counts$likely_germline),
    sprintf("- Número de likely_artifact: %d", counts$likely_artifact),
    sprintf("- Número de technical_fail: %d", counts$technical_fail),
    sprintf("- Número de manual_review_required: %d", review_required),
    sprintf("- Número de uncertain_tumor_only: %d", counts$uncertain_tumor_only),
    "",
    "## Impacto dos filtros",
    sprintf("- Número removido por evidência populacional: %d", populational_removals),
    sprintf("- Número removido por artefato técnico: %d", artifact_removals),
    sprintf("- Número mantido como high_confidence_somatic: %d", high_confidence),
    sprintf("- Número mantido como probable_somatic: %d", probable_somatic),
    sprintf("- Número enviado para revisão manual: %d", review_required),
    "",
    "## Principais categorias populacionais",
    paste0("- ", names(population_sources), ": ", population_sources),
    "",
    "## Principais categorias de artefato",
    paste0("- ", names(artifact_sources), ": ", artifact_sources),
    "",
    "## Principais razões de filtro",
    paste0("- ", names(filter_reasons), ": ", filter_reasons),
    "",
    "## Impacto no TMB",
    if (!is.null(tmb_summary)) {
      sapply(seq_len(nrow(tmb_summary)), function(i) {
        sprintf("- Amostra %s (%s): %0.2f mut/Mb (%s)",
                tmb_summary$sample_id[i], tmb_summary$tumor_type[i], tmb_summary$tmb_mut_per_mb[i], tmb_summary$tmb_category[i])
      })
    } else {
      "- TMB não foi calculado."
    },
    "",
    "## Limitações da análise tumor-only",
    "- Como a análise é tumor-only, as variantes classificadas como high_confidence_somatic ou probable_somatic são inferências baseadas em evidência técnica, populacional e oncológica, não confirmação definitiva de somaticidade.",
    "- Variantes com sinal forte de filtro técnico ou AF populacional alta não são resgatadas automaticamente por COSMIC/OncoKB.",
    "- A confiabilidade depende da qualidade dos metadados, cobertura e bancos populacionais disponíveis.",
    ""
  )
  writeLines(unlist(lines), report_path)
}

class_summary <- function(x) {
  tbl <- as.data.frame(table(x$final_class), stringsAsFactors = FALSE)
  names(tbl) <- c("final_class", "n_variants")
  tbl$fraction <- tbl$n_variants / sum(tbl$n_variants)
  tbl
}

validation_summary <- function(x) {
  data.frame(
    metric = c(
      "n_variants",
      "oncokb_matches",
      "cosmic_matches",
      "oncokb_or_cosmic_matches",
      "validation_rescue_candidates",
      "validated_high_confidence_somatic",
      "validated_probable_germline",
      "validated_probable_artifact"
    ),
    value = c(
      nrow(x),
      sum(x$oncokb_match, na.rm = TRUE),
      sum(x$cosmic_match, na.rm = TRUE),
      sum(x$oncokb_match | x$cosmic_match, na.rm = TRUE),
      sum(x$validation_rescue_candidate, na.rm = TRUE),
      sum((x$oncokb_match | x$cosmic_match) & x$final_class == "high_confidence_somatic", na.rm = TRUE),
      sum((x$oncokb_match | x$cosmic_match) & x$final_class == "probable_germline", na.rm = TRUE),
      sum((x$oncokb_match | x$cosmic_match) & x$final_class == "probable_artifact", na.rm = TRUE)
    )
  )
}

driver_summary <- function(x) {
  tbl <- as.data.frame(table(x$driver_class), stringsAsFactors = FALSE)
  names(tbl) <- c("driver_class", "n_variants")
  tbl$fraction <- tbl$n_variants / sum(tbl$n_variants)
  tbl
}

dndscv_input_table <- function(x) {
  keep <- x[x$final_class %in% c("high_confidence_somatic", "probable_somatic"), , drop = FALSE]
  data.frame(
    sampleID = keep$sample_id,
    chr = gsub("^chr", "", keep$chrom, ignore.case = TRUE),
    pos = keep$pos,
    ref = keep$ref,
    mut = keep$alt,
    gene = keep$gene,
    tumor_type = keep$tumor_type,
    stringsAsFactors = FALSE
  )
}

chasmplus_input_table <- function(x) {
  keep <- x[x$final_class %in% c("high_confidence_somatic", "probable_somatic") &
            grepl("missense", tolower(as.character(x$consequence))), , drop = FALSE]
  data.frame(
    sample_id = keep$sample_id,
    tumor_type = keep$tumor_type,
    gene = keep$gene,
    protein_change = keep$protein_change,
    chrom = keep$chrom,
    pos = keep$pos,
    ref = keep$ref,
    alt = keep$alt,
    driver_class = keep$driver_class,
    driver_score = keep$driver_score,
    stringsAsFactors = FALSE
  )
}
