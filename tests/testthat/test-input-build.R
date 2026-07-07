test_that("genome build is detected from the VCF header (##reference + contig)", {
  meta <- read_vcf_meta(fixture("sample.vep.vcf"))
  expect_equal(detect_genome_build(meta), "GRCh38")
})

test_that("VEP CSQ format is parsed dynamically from the header", {
  meta <- read_vcf_meta(fixture("sample.vep.vcf"))
  fields <- parse_vcf_annotation_format(meta, "CSQ")
  expect_true(all(c("Consequence", "SYMBOL", "HGVSp", "gnomADe_AF") %in% fields))
  x <- read_annotated_vcf(fixture("sample.vep.vcf"))
  expect_true(all(c("SYMBOL", "Consequence") %in% names(x)))
  expect_true("BRAF" %in% x$SYMBOL)
})

test_that("config build prevails but conflict with header errors", {
  cfg <- list(input = list(genome_build = "GRCh37"))
  expect_error(resolve_genome_build(cfg, fixture("sample.vep.vcf")), "conflict", ignore.case = TRUE)
})

test_that("undeterminable build stops with an actionable message", {
  cfg <- list(input = list())
  # a TSV path (no header) and no config build -> cannot determine
  tmp <- tempfile(fileext = ".tsv"); writeLines("a\tb", tmp)
  expect_error(resolve_genome_build(cfg, tmp), "genome_build", ignore.case = TRUE)
})

test_that("agreeing config and header resolve without error", {
  cfg <- list(input = list(genome_build = "GRCh38"))
  expect_equal(resolve_genome_build(cfg, fixture("sample.vep.vcf")), "GRCh38")
})
