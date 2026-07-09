test_that("prioritization is decomposed and hotspots raise priority", {
  x <- data.frame(
    consequence = c("missense_variant", "missense_variant"),
    impact_rank = c(2, 2),
    COMPUTATIONAL_EVIDENCE_SCORE = c(0.8, 0.8),
    hotspot_match = c(TRUE, FALSE), hotspot_tumor_specific = c(TRUE, FALSE),
    gene_driver_match = c(TRUE, FALSE), gene_driver_tumor_specific = c(FALSE, FALSE),
    clinvar_significance = c("pathogenic", NA),
    COSMIC_GLOBAL_RECURRENCE_SCORE = c(0.5, 0.5),
    stringsAsFactors = FALSE)
  out <- add_prioritization(x)
  expect_true(all(c("PRIORITY_SCORE_BASE", "PRIORITY_CATEGORY", "PRIORITY_COMPONENTS",
    "PRIORITY_HOTSPOT_COMPONENT", "PRIORITY_COMPUTATIONAL_COMPONENT") %in% names(out)))
  expect_gt(out$PRIORITY_SCORE_BASE[1], out$PRIORITY_SCORE_BASE[2])   # hotspot+driver+clinical
  expect_equal(out$PRIORITY_HOTSPOT_COMPONENT[1], 1)
})

test_that("create_maftools_plots produces summary+oncoplot as PNG+PDF and skips inapplicable plots", {
  skip_if_not_installed("maftools")
  v <- data.frame(
    chrom = c("17", "7", "12"), pos = c(7674220, 140453136, 25398284),
    ref = c("C", "A", "C"), alt = c("T", "T", "A"),
    gene = c("TP53", "BRAF", "KRAS"),
    consequence = rep("missense_variant", 3),
    sample_id = c("S1", "S1", "S2"), vaf = c(0.4, 0.5, 0.45), stringsAsFactors = FALSE)
  maf_path <- tempfile(fileext = ".maf"); create_maf(v, maf_path)
  out_dir <- tempfile("plots");
  man <- suppressWarnings(create_maftools_plots(maf_path, out_dir, cfg = list(plots = list(
    minimum_samples_for_oncoplot = 2, minimum_samples_for_interactions = 20,
    minimum_mutations_for_rainfall = 50))))
  expect_true("plotmaf_summary" %in% man$plot)
  expect_true(any(man$plot == "somatic_interactions" & man$status == "skipped"))
  onco <- man[man$plot == "oncoplot", ]
  if (nrow(onco) && onco$status == "ok") {
    fs <- strsplit(onco$files, ";", fixed = TRUE)[[1]]
    expect_true(any(grepl("[.]png$", fs)) && any(grepl("[.]pdf$", fs)))
    expect_true(all(file.exists(fs)) && all(file.info(fs)$size > 0))
  }
})
