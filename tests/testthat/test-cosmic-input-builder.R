cosmic_input_cfg <- function(cache_dir) list(cosmic = list(
  release = "TEST", cache_dir = cache_dir,
  input_sources = list(
    genome_screens_tsv = fixture("cosmic_genome_screens_small.tsv"),
    noncoding_tsv = fixture("cosmic_noncoding_small.tsv"),
    genome_screens_normal_vcf = fixture("cosmic_genome_screens_normal_small.vcf"),
    noncoding_normal_vcf = fixture("cosmic_noncoding_normal_small.vcf"),
    classification_tsv = fixture("cosmic_classification_small.tsv"))))

test_that("build_cosmic_raw_input joins coordinates (from Normal VCF) with tumor context (Classification)", {
  cache <- tempfile("cosmicinput")
  res <- build_cosmic_raw_input(cosmic_input_cfg(cache))
  expect_false(res$reused)
  expect_true(file.exists(res$raw_file_path))
  m <- res$manifest
  expect_equal(m$n_input, 4)
  expect_equal(m$n_unmatched_to_normal_vcf, 1)
  expect_equal(m$n_output, 3)

  out <- read.delim(res$raw_file_path, stringsAsFactors = FALSE)
  # coordinates come from the Normal VCF, not the mutation TSV's own columns
  row <- out[out$GENOMIC_MUTATION_ID == "COSV70830383", ][1, ]
  expect_equal(row$GENOME_START, 10151)
  expect_equal(row$GENOMIC_WT_ALLELE, "T")
  expect_equal(row$GENOMIC_MUT_ALLELE, "A")
  # tumor context joined via COSMIC_PHENOTYPE_ID
  expect_true(all(out$PRIMARY_SITE == "ovary"))
  expect_true(all(out$PRIMARY_HISTOLOGY == "carcinoma"))
})

test_that("mutation rows with no coordinate match are excluded and reported, not dropped", {
  cache <- tempfile("cosmicinput")
  res <- build_cosmic_raw_input(cosmic_input_cfg(cache))
  unmatched <- read.delim(file.path(dirname(res$raw_file_path), "unmatched_to_normal_vcf.tsv.gz"),
                           stringsAsFactors = FALSE)
  expect_equal(nrow(unmatched), 1)
  expect_equal(unmatched$GENOMIC_MUTATION_ID, "COSVFAKE00000001")
})

test_that("the output feeds prepare_cosmic_db() directly (standardize_cosmic_raw compatible)", {
  cache <- tempfile("cosmicinput")
  built <- build_cosmic_raw_input(cosmic_input_cfg(cache))
  cosmic_cache <- tempfile("cosmicdb")
  res <- prepare_cosmic_db(list(
    cosmic = list(raw_file = built$raw_file_path, release = "vTEST",
                  source_build = "GRCh38", target_build = "GRCh38", cache_dir = cosmic_cache),
    reference = list(fasta = NULL)))
  expect_true(file.exists(res$db_path))
  # 3 sample-occurrence rows -> 2 unique canonical keys (COSV70830383 has two
  # occurrences, DDX11L1/WASH7P transcripts of the same genomic mutation) -
  # this is the point of feeding one-row-per-occurrence into prepare_cosmic_db:
  # recurrence is counted, not silently deduplicated to 1.
  expect_equal(res$manifest$n_converted, 2)
  # occurrence counting (no explicit COSMIC_COUNT column -> each row = 1
  # occurrence, summed) lives in cosmic_context_long, the primary structure
  # for tumor-type-stratified matching (docs/COSMIC_DECISION_RULES.md).
  long <- as.data.frame(readRDS(file.path(dirname(res$db_path), "cosmic_context_long.rds")))
  recurrent <- long[long$canonical_key == "GRCh38|1|10151|T|A", ]
  expect_equal(recurrent$occurrence_count, 2)
  expect_equal(recurrent$cosmic_primary_site, "ovary")
})

test_that("preparation is idempotent: second call reuses (no rebuild)", {
  cache <- tempfile("cosmicinput")
  cfg <- cosmic_input_cfg(cache)
  first <- build_cosmic_raw_input(cfg)
  mtime1 <- file.info(first$raw_file_path)$mtime
  second <- build_cosmic_raw_input(cfg)
  expect_true(second$reused)
  expect_equal(file.info(second$raw_file_path)$mtime, mtime1)
})

test_that("missing input_sources file errors clearly (never silently proceeds)", {
  cache <- tempfile("cosmicinput")
  cfg <- cosmic_input_cfg(cache)
  cfg$cosmic$input_sources$classification_tsv <- "/no/such/file.tsv.gz"
  expect_error(build_cosmic_raw_input(cfg), "classification_tsv", ignore.case = TRUE)
})
