test_that("format is detected by content, not extension alone (VCF)", {
  expect_equal(detect_input_format(fixture("sample.vep.vcf")), "vcf")
  res <- read_variant_input(fixture("sample.vep.vcf"))
  expect_equal(res$format, "vcf")
  expect_true(all(c("CHROM", "POS", "REF", "ALT", "GENE", "CONSEQUENCE", "SAMPLE_ID") %in% names(res$variants)))
  expect_true(nrow(res$schema_report) == length(CANONICAL_FIELDS))
})

test_that("VCF header metadata and build are interpreted", {
  res <- read_variant_input(fixture("sample.vep.vcf"))
  expect_equal(res$header_metadata$fileformat, "VCFv4.2")
  expect_equal(detect_genome_build(res$header_metadata$raw), "GRCh38")
})

test_that("MAF is detected, comments parsed, VAF derived from allele depths", {
  maf <- data.frame(Hugo_Symbol = "TP53", Chromosome = "17", Start_Position = 7674220,
    End_Position = 7674220, Reference_Allele = "C", Tumor_Seq_Allele2 = "T",
    Variant_Classification = "Missense_Mutation", Variant_Type = "SNP",
    Tumor_Sample_Barcode = "S1", t_ref_count = 60, t_alt_count = 40, stringsAsFactors = FALSE)
  mf <- tempfile(fileext = ".maf")
  writeLines(c("#version 2.4", "#genome_build GRCh38"), mf)
  suppressWarnings(data.table::fwrite(maf, mf, sep = "\t", append = TRUE, col.names = TRUE))
  res <- read_variant_input(mf)
  expect_equal(res$format, "maf")
  expect_equal(res$variants$CHROM[1], "17")
  expect_equal(round(res$variants$VAF[1], 2), 0.4)
  expect_equal(res$header_metadata$genome_build, "GRCh38")
})

test_that("custom TSV uses column_map; unresolved minimum fields raise an actionable error", {
  tsvf <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(CHROMOSOME = "7", START = 140453136, REF_ALLELE = "A",
    ALT_ALLELE = "T", SAMPLE = "X", SYMBOL = "BRAF"), tsvf, sep = "\t")
  res <- read_variant_input(tsvf, column_map = list(chrom = "CHROMOSOME", pos = "START",
    ref = "REF_ALLELE", alt = "ALT_ALLELE", sample_id = "SAMPLE", gene = "SYMBOL"))
  expect_equal(res$format, "tsv")
  expect_equal(res$variants$GENE[1], "BRAF")
  expect_true("column_map" %in% res$schema_report$detection_method)

  bad <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(foo = 1, bar = 2), bad, sep = "\t")
  expect_error(read_variant_input(bad), "required canonical", ignore.case = TRUE)
})

test_that("original columns are preserved (no silent loss)", {
  tsvf <- tempfile(fileext = ".tsv")
  data.table::fwrite(data.frame(CHROM = "1", POS = 100, REF = "A", ALT = "G",
    SAMPLE_ID = "s", my_custom_col = "keep_me"), tsvf, sep = "\t")
  res <- read_variant_input(tsvf)
  expect_true("my_custom_col" %in% names(res$variants))
})
