# Hermetic mini-run: a small inline VEP-like TSV + minimal config, so the test
# is self-contained and runs under R CMD check (no repo data/ or config/).
demo_run <- function() {
  variants <- data.frame(
    Tumor_Sample_Barcode = c("S1", "S2", "S3", "S4", "S5"),
    CHROM = c("17", "17", "7", "7", "12"),
    START = c(7673803, 7673803, 140453136, 140453136, 25398284),
    REF = c("C", "C", "A", "A", "C"),
    ALT = c("T", "T", "T", "T", "A"),
    Hugo_Symbol = c("TP53", "TP53", "BRAF", "BRAF", "KRAS"),
    HGVSp_Short = c("p.R175H", "p.R175H", "p.V600E", "p.V600E", "p.G12V"),
    Consequence = "missense_variant",
    DP = c(80, 75, 90, 8, 70), alt_count = c(38, 35, 45, 2, 35),
    AF = c(0.48, 0.47, 0.5, 0.25, 0.5),
    TLOD = c(80, 78, 85, 3, 70), MBQ = 35, MMQ = 60,
    FILTER = c("PASS", "PASS", "PASS", "weak_evidence", "PASS"),
    gnomADg_AF = c(0.01, 0.01, 0.0, 0.0, 0.0),
    stringsAsFactors = FALSE)
  vfile <- tempfile(fileext = ".tsv")
  write_tsv(variants, vfile)
  cfg <- list(
    input = list(vcf = vfile, delimiter = "\t", genome_build = "GRCh38"),
    analysis = list(output_dir = tempfile("run_out"), run_id = "t"),
    cancer = list(default_tumor_type = "PANCANCER"),
    hard_filters = list(enabled = TRUE, assay = "WGS"),
    technical_filters = list(min_depth = 20, min_alt_count_snv = 5, min_alt_count_indel = 8,
                             min_af_snv = 0.03, min_af_indel = 0.05, min_tlod = 6,
                             min_mbq = 25, min_mmq = 40),
    population_filters = list(common_af_threshold = 0.01, low_frequency_af_threshold = 0.001),
    cohort = list(recurrent_variant_fraction_artifact = 0.3, recurrent_locus_fraction_artifact = 0.4),
    guideline_classification = list(enabled = TRUE))
  suppressWarnings(suppressMessages(run_tumor_only(cfg)))
}

test_that("OncoKB annotation does NOT change filter_status (regression guarantee)", {
  res <- demo_run()
  all_path <- file.path(res$run_dir, "tables", "variants_all.tsv.gz")
  before <- read_variants(all_path, "\t")$filter_status
  Sys.unsetenv("ONCOKB_TOKEN")               # no token -> offline path, still completes
  suppressMessages(annotate_oncokb(res$run_dir, scope = "all"))
  after <- read_variants(all_path, "\t")$filter_status
  expect_identical(before, after)            # the explicit invariant from the brief
})

test_that("annotate_oncokb completes without a token and marks not-annotated", {
  res <- demo_run()
  Sys.unsetenv("ONCOKB_TOKEN")
  ann <- suppressMessages(annotate_oncokb(res$run_dir, scope = "retained"))
  expect_true(file.exists(file.path(res$run_dir, "tables", "oncokb_annotations.tsv.gz")))
  expect_true(all(ann$ONCOKB_ANNOTATED %in% c(FALSE, NA)))
  expect_true(all(ann$ONCOKB_QUERY_STATUS[ann$filter_status %in% c("PASS","REVIEW")] %in%
                    c("no_token", "not_in_scope")))
})

test_that("OncoKB cache round-trips a record", {
  cd <- tempfile("okbcache"); dir.create(cd)
  rec <- list(ONCOKB_ANNOTATED = TRUE, ONCOKB_ONCOGENIC = "Oncogenic",
              ONCOKB_QUERY_STATUS = "annotated")
  oncokb_cache_put(cd, "GRCh38|7|140453136|A|T", rec)
  got <- oncokb_cache_get(cd, "GRCh38|7|140453136|A|T")
  expect_equal(got$ONCOKB_ONCOGENIC, "Oncogenic")
})
