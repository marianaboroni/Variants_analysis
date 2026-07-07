test_that("VEP consequences map to MAF Variant_Classification explicitly", {
  vt <- c("SNP", "SNP", "INS", "DEL")
  cls <- map_consequence_to_maf(
    c("missense_variant", "stop_gained", "frameshift_variant", "frameshift_variant"), vt)
  expect_equal(as.character(cls),
               c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del"))
})

test_that("unknown consequences fall into a controlled bucket and warn", {
  expect_warning(
    v <- create_maf(data.frame(chrom = "1", pos = 100, ref = "A", alt = "G",
                               gene = "X", consequence = "totally_made_up_variant",
                               sample_id = "S1")),
    "unmapped", ignore.case = TRUE)
  expect_equal(v$Variant_Classification, MAF_UNKNOWN_CLASS)
})

test_that("create_maf builds required MAF fields and validates with maftools", {
  skip_if_not_installed("maftools")
  v <- data.frame(
    chrom = c("17", "7"), pos = c(7674220, 140453136),
    ref = c("C", "A"), alt = c("T", "T"),
    gene = c("TP53", "BRAF"),
    consequence = c("missense_variant", "missense_variant"),
    sample_id = c("S1", "S2"), vaf = c(0.4, 0.5), stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".maf")
  maf <- create_maf(v, path)
  req <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
           "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
           "Variant_Type", "Tumor_Sample_Barcode")
  expect_true(all(req %in% names(maf)))
  obj <- maftools::read.maf(maf = read_variants(path, "\t"), verbose = FALSE)
  expect_s4_class(obj, "MAF")
})
