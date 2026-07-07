make_cosmic_cfg <- function(cache_dir, raw = fixture("cosmic_raw.tsv")) {
  list(
    input = list(genome_build = "GRCh38"),
    cosmic = list(raw_file = raw, release = "vTEST",
                  source_build = "GRCh38", target_build = "GRCh38",
                  cache_dir = cache_dir),
    reference = list(fasta = NULL))
}

test_that("prepare_cosmic_db builds a same-build DB with manifest and provenance", {
  cache <- file.path(tempfile("cosmic"))
  res <- prepare_cosmic_db(make_cosmic_cfg(cache))
  expect_false(res$reused)
  expect_true(file.exists(res$db_path))
  expect_true(file.exists(res$manifest_path))
  m <- res$manifest
  # 6 input rows: 1 dropped (missing pos), 1 multiallelic -> 2, 1 duplicate key collapsed
  expect_equal(m$n_input, 6)
  expect_equal(m$n_dropped_missing, 1)      # COSV5 has no position
  expect_equal(m$n_duplicates, 1)           # COSV2b duplicates COSV2 key
  expect_true(all(file.exists(file.path(dirname(res$db_path),
    c("unmapped.tsv.gz", "ref_mismatch.tsv.gz", "duplicates.tsv.gz", "dropped_missing.tsv.gz")))))
})

test_that("duplicate canonical keys are collapsed and counts aggregated", {
  cache <- file.path(tempfile("cosmic"))
  res <- prepare_cosmic_db(make_cosmic_cfg(cache))
  db <- as.data.frame(readRDS(res$db_path))
  kras <- db[db$canonical_key == "GRCh38|12|25398284|C|A", ]
  expect_equal(nrow(kras), 1)
  expect_equal(kras$COSMIC_OCCURRENCE_COUNT, 3100)          # 3000 + 100 summed
  expect_true(grepl("COSV2", kras$COSMIC_MUTATION_IDS))     # original IDs preserved
})

test_that("preparation is idempotent: second call reuses (no rebuild)", {
  cache <- file.path(tempfile("cosmic"))
  cfg <- make_cosmic_cfg(cache)
  first <- prepare_cosmic_db(cfg)
  mtime1 <- file.info(first$db_path)$mtime
  second <- prepare_cosmic_db(cfg)
  expect_true(second$reused)
  expect_equal(file.info(second$db_path)$mtime, mtime1)     # file untouched
})

test_that("changed input invalidates the cache and requires --force", {
  cache <- file.path(tempfile("cosmic"))
  cfg <- make_cosmic_cfg(cache)
  prepare_cosmic_db(cfg)
  # simulate a changed raw file by pointing at a modified copy under same release/build
  alt <- tempfile(fileext = ".tsv")
  x <- read_variants(fixture("cosmic_raw.tsv"), "\t"); x <- x[1:2, ]
  write_tsv(x, alt)
  cfg2 <- make_cosmic_cfg(cache, raw = alt)
  expect_error(prepare_cosmic_db(cfg2, force = FALSE), "force", ignore.case = TRUE)
  expect_silent_build <- prepare_cosmic_db(cfg2, force = TRUE)
  expect_false(expect_silent_build$reused)
})

test_that("cross-build liftover with no backend errors clearly (never silent)", {
  cache <- file.path(tempfile("cosmic"))
  cfg <- make_cosmic_cfg(cache)
  cfg$cosmic$source_build <- "GRCh37"
  cfg$cosmic$target_build <- "GRCh38"
  cfg$cosmic$chain_file <- NULL
  skip_if(requireNamespace("rtracklayer", quietly = TRUE),
          "rtracklayer present; cross-build path would run instead of erroring")
  expect_error(prepare_cosmic_db(cfg), "liftover backend", ignore.case = TRUE)
})

test_that("COSMIC matching uses the canonical key (chrom/pos/ref/alt), not gene", {
  cache <- file.path(tempfile("cosmic"))
  prepare_cosmic_db(make_cosmic_cfg(cache))
  x <- data.frame(
    chrom = c("chr7", "12", "17"),
    pos = c(140453136, 25398284, 999999),   # 3rd has no COSMIC record at this pos
    ref = c("A", "C", "G"), alt = c("T", "A", "C"),
    stringsAsFactors = FALSE)
  cfg <- list(cosmic = list(release = "vTEST", source_build = "GRCh38", cache_dir = cache))
  out <- annotate_cosmic(x, cfg, "GRCh38")
  expect_equal(out$COSMIC_MATCH, c(TRUE, TRUE, FALSE))
  expect_equal(out$cosmic_count[2], 3100)
})
