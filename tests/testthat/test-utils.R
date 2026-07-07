test_that("chromosome names are normalized (chr prefix, M/MT)", {
  expect_equal(normalize_chrom(c("chr1", "1", "CHRX", "chrM", "MT")),
               c("1", "1", "X", "MT", "MT"))
})

test_that("allele normalization trims common bases and left-aligns", {
  # SNV unchanged
  s <- normalize_alleles(100, "C", "T")
  expect_equal(c(s$pos, s$ref, s$alt), c("100", "C", "T"))
  # insertion: CAT>CGAT should trim shared prefix C -> A>GA at pos+1 style
  ins <- normalize_alleles(100, "C", "CA")
  expect_true(nchar(ins$alt) >= nchar(ins$ref))
  # deletion with shared prefix: AAT>A  -> minimal representation retains anchor
  del <- normalize_alleles(100, "AAT", "A")
  expect_equal(del$ref, "AAT")   # no interior trim possible while keeping >=1 base
  # shared suffix trimmed: TAG>TG (share G) -> TA>T
  suf <- normalize_alleles(100, "TAG", "TG")
  expect_equal(c(suf$ref, suf$alt), c("TA", "T"))
})

test_that("canonical key uses build|chrom|pos|ref|alt with normalization", {
  k <- canonical_variant_key("GRCh38", "chr17", 7674220, "C", "T")
  expect_equal(k, "GRCh38|17|7674220|C|T")
})

test_that("multiallelic ALT splits into one row per allele", {
  df <- data.frame(chrom = "1", pos = 100, ref = "G", alt = "A,T", stringsAsFactors = FALSE)
  out <- split_multiallelic(df, "alt")
  expect_equal(nrow(out), 2)
  expect_setequal(out$alt, c("A", "T"))
})
