abraom_prep_cfg <- function(cache_dir) list(
  population = list(brazilian = list(
    raw_file = fixture("abraom_raw.tsv"), name = "TEST", genome_build = "GRCh38",
    cache_dir = cache_dir)),
  reference = list(fasta = fixture("chr1_subset.fa")))

test_that("prepare_brazilian_db reconstructs VCF-standard indels and validates against FASTA", {
  skip_if_not_installed("Biostrings")
  cache <- tempfile("abraom")
  res <- prepare_brazilian_db(abraom_prep_cfg(cache))
  expect_false(res$reused)
  expect_true(file.exists(res$db_path))
  m <- res$manifest
  expect_equal(m$n_input, 5)
  expect_equal(m$n_ref_mismatch, 2)   # the two deliberately-wrong rows
  expect_equal(m$n_output, 3)

  db <- read.delim(res$db_path, stringsAsFactors = FALSE)
  snv <- db[db$Variant_Type == "SNV", ]
  expect_equal(snv$POS, 10380); expect_equal(snv$REF, "C"); expect_equal(snv$ALT, "T")

  del <- db[db$Variant_Type == "deletion", ]
  expect_equal(del$POS, 10378); expect_equal(del$REF, "CCCTAA"); expect_equal(del$ALT, "C")

  ins <- db[db$Variant_Type == "insertion", ]
  expect_equal(ins$POS, 10433); expect_equal(ins$REF, "A"); expect_equal(ins$ALT, "AC")
})

test_that("ref_mismatch rows are excluded from the database but never dropped", {
  skip_if_not_installed("Biostrings")
  cache <- tempfile("abraom")
  res <- prepare_brazilian_db(abraom_prep_cfg(cache))
  mismatch <- read.delim(file.path(dirname(res$db_path), "ref_mismatch.tsv.gz"), stringsAsFactors = FALSE)
  expect_equal(nrow(mismatch), 2)
  expect_true(all(c(40000, 10380) %in% mismatch$orig_pos))
})

test_that("preparation is idempotent: second call reuses (no rebuild)", {
  skip_if_not_installed("Biostrings")
  cache <- tempfile("abraom")
  cfg <- abraom_prep_cfg(cache)
  first <- prepare_brazilian_db(cfg)
  mtime1 <- file.info(first$db_path)$mtime
  second <- prepare_brazilian_db(cfg)
  expect_true(second$reused)
  expect_equal(file.info(second$db_path)$mtime, mtime1)
})

test_that("missing reference.fasta errors clearly (never silently proceeds)", {
  cache <- tempfile("abraom")
  cfg <- abraom_prep_cfg(cache)
  cfg$reference$fasta <- NULL
  expect_error(prepare_brazilian_db(cfg), "reference.fasta", ignore.case = TRUE)
})
