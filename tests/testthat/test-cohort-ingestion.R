# A cohort's caller-FILTER/artifact evidence relies on cohort-wide recurrence
# (compute_cohort_recurrence() / cohort.recurrent_variant_fraction_artifact) to
# be meaningful. That requires more than one distinct sample_id in the table
# handed to run_tumor_only() in a single run. This file covers ingesting a
# cohort straight from a list of single-sample files (input.path as a YAML
# list), as an alternative to pre-merging with `bcftools merge`.
# See docs/FILTERING_STRATEGY.md, "Cohort-wide ingestion".

second_sample_vcf <- function(sample_name = "TUMOR2") {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    "##contig=<ID=1,length=248956422>",
    "##contig=<ID=7,length=159345973>",
    "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSp|gnomADe_AF|CLIN_SIG\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele fraction\">",
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", sample_name),
    # same BRAF variant as sample.vep.vcf's TUMOR1 -> recurrent across the cohort
    "7\t140453136\t.\tA\tT\t500\tPASS\tCSQ=T|missense_variant|MODERATE|BRAF|ENSG00000157764|ENST00000288602|p.Val600Glu|0.0|.\tGT:AD:DP:AF\t0/1:30,28:58:0.48",
    # private to this sample
    "1\t100\t.\tG\tC\t400\tPASS\tCSQ=C|missense_variant|MODERATE|GENEX|ENSG0000000000|ENST00000000000|p.X1Y|0.0|.\tGT:AD:DP:AF\t0/1:20,20:40:0.50"
  ), vcf)
  vcf
}

test_that("as_input_paths unwraps a YAML list into a character vector, passes a scalar through", {
  expect_equal(as_input_paths("a.vcf"), "a.vcf")
  expect_equal(as_input_paths(list("a.vcf", "b.vcf")), c("a.vcf", "b.vcf"))
  expect_null(as_input_paths(NULL))
})

test_that("validate_config accepts input.path as a list and checks every file", {
  present <- tempfile(fileext = ".vcf"); writeLines("x", present)
  cfg_ok <- list(input = list(path = list(present, present)), analysis = list(output_dir = "o"))
  expect_invisible(validate_config(cfg_ok))

  cfg_missing <- list(input = list(path = list(present, "/no/such/file.vcf.gz")),
                       analysis = list(output_dir = "o"))
  expect_error(validate_config(cfg_missing), "no/such/file", fixed = TRUE)
})

test_that("resolve_config normalizes a YAML list into a plain character vector", {
  cfg <- resolve_config(list(input = list(path = list("a.vcf", "b.vcf")), output = list(dir = "res")))
  expect_true(is.character(cfg$input$vcf))
  expect_equal(cfg$input$vcf, c("a.vcf", "b.vcf"))
})

test_that("resolve_run_id is stable across file order and differs from a single-file run", {
  a <- tempfile(fileext = ".vcf"); writeLines("a", a)
  b <- tempfile(fileext = ".vcf"); writeLines("b", b)
  cfg1 <- list(input = list(vcf = c(a, b)), analysis = list())
  cfg2 <- list(input = list(vcf = c(b, a)), analysis = list())
  cfg_single <- list(input = list(vcf = a), analysis = list())
  expect_equal(resolve_run_id(cfg1), resolve_run_id(cfg2))
  expect_false(resolve_run_id(cfg1) == resolve_run_id(cfg_single))
})

test_that("file_sha256 is vectorized: one hash per file, NA for missing", {
  a <- tempfile(fileext = ".vcf"); writeLines("a", a)
  hashes <- file_sha256(c(a, "/no/such/file"))
  expect_length(hashes, 2)
  expect_false(is.na(hashes[1]))
  expect_true(is.na(hashes[2]))
})

test_that("resolve_genome_build detects a shared build across a cohort and errors on conflict", {
  vcf1 <- fixture("sample.vep.vcf")
  vcf2 <- second_sample_vcf()
  expect_equal(resolve_genome_build(list(), c(vcf1, vcf2)), "GRCh38")

  vcf37 <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##reference=GRCh37",
               "##contig=<ID=1,length=249250621>",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "1\t100\t.\tA\tT\t.\tPASS\t."), vcf37)
  expect_error(resolve_genome_build(list(), c(vcf1, vcf37)), "conflict", ignore.case = TRUE)
})

test_that("read_variant_input ingests a list of VCFs as one cohort table (no bcftools merge needed)", {
  vcf1 <- fixture("sample.vep.vcf")
  vcf2 <- second_sample_vcf()
  res <- read_variant_input(c(vcf1, vcf2))

  expect_equal(nrow(res$variants), 4)  # 2 variants/file x 2 files
  expect_setequal(unique(res$variants$SAMPLE_ID), c("TUMOR1", "TUMOR2"))
  expect_true("SOURCE_FILE" %in% names(res$variants))
  expect_setequal(unique(res$variants$SOURCE_FILE), c(vcf1, vcf2))
  expect_equal(res$format, "vcf")
})

test_that("a cohort ingested this way makes cohort-wide recurrence meaningful", {
  vcf1 <- fixture("sample.vep.vcf")
  vcf2 <- second_sample_vcf()
  res <- read_variant_input(c(vcf1, vcf2))
  x <- standardize_variant_table(res$variants, cfg = NULL)
  cfg <- list(technical_filters = list(min_alt_count_snv = 5, min_alt_count_indel = 8,
                                        min_af_snv = 0.03, min_af_indel = 0.05, min_depth = 20,
                                        min_tlod = 6, min_mbq = 25, min_mmq = 40))
  x <- add_basic_features(x, cfg)

  rec <- compute_cohort_recurrence(x, cfg)
  expect_equal(rec$total_samples, 2)

  braf <- rec$by_variant[rec$by_variant$chrom == "7" & rec$by_variant$pos == 140453136, ]
  expect_equal(braf$variant_n_samples, 2)
  expect_equal(braf$variant_cohort_freq, 1)

  private_var <- rec$by_variant[rec$by_variant$chrom == "1" & rec$by_variant$pos == 100, ]
  expect_equal(private_var$variant_n_samples, 1)
  expect_equal(private_var$variant_cohort_freq, 0.5)
})

test_that("a single-file input.path is completely unaffected (no SOURCE_FILE, no behavior change)", {
  res <- read_variant_input(fixture("sample.vep.vcf"))
  expect_false("SOURCE_FILE" %in% names(res$variants))
  expect_equal(nrow(res$variants), 2)
})
