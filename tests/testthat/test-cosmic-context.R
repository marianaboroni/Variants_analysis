melanoma_allowed <- function() {
  data.frame(site = c("skin", "skin"), hist = c("malignant_melanoma", ""),
             level = c("exact", "compatible"), via_subtype = c(FALSE, FALSE),
             stringsAsFactors = FALSE)
}
lr <- function(...) {  # long rows: each c(site, hist, subtype, count)
  rows <- list(...)
  data.frame(
    site = vapply(rows, function(r) r[[1]], character(1)),
    hist = vapply(rows, function(r) r[[2]], character(1)),
    subtype = vapply(rows, function(r) r[[3]], character(1)),
    count = vapply(rows, function(r) as.numeric(r[[4]]), numeric(1)),
    stringsAsFactors = FALSE)
}
cls <- function(long, allowed = melanoma_allowed(), known = TRUE, sub = FALSE, pan = 5)
  classify_cosmic_tumor_context(long, allowed, known, sample_subtype_available = sub, pan_threshold = pan)

# --- the ten required scenarios ---------------------------------------------

test_that("1. same type + histology -> exact_match", {
  r <- cls(lr(c("skin","malignant_melanoma","",40)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "exact_match")
  expect_equal(r$COSMIC_MATCH_LEVEL, "exact_site_histology")
  expect_equal(r$COSMIC_MATCHING_TUMOR_OCCURRENCES, 40)
  expect_equal(r$COSMIC_OTHER_TUMOR_OCCURRENCES, 0)
})
test_that("2. same organ, different histology -> compatible_match", {
  r <- cls(lr(c("skin","basal_cell_carcinoma","",12)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "compatible_match")
  expect_equal(r$COSMIC_MATCH_LEVEL, "compatible_histology")
})
test_that("3. same type once -> exact_match, count 1", {
  r <- cls(lr(c("skin","malignant_melanoma","",1)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "exact_match")
  expect_equal(r$COSMIC_MATCHING_TUMOR_OCCURRENCES, 1)
})
test_that("4. recurrent only in other tumor -> other_tumor_only", {
  r <- cls(lr(c("breast","ductal_carcinoma","",30)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "other_tumor_only")
  expect_equal(r$COSMIC_MATCHING_TUMOR_OCCURRENCES, 0)
  expect_equal(r$COSMIC_OTHER_TUMOR_OCCURRENCES, 30)
})
test_that("5a. pan-cancer including sample -> pan_cancer_with_sample_tumor", {
  r <- cls(lr(c("skin","malignant_melanoma","",20), c("lung","adenocarcinoma","",40),
              c("large_intestine","adenocarcinoma","",60), c("pancreas","ductal_carcinoma","",30),
              c("breast","ductal_carcinoma","",15), c("ovary","serous_carcinoma","",10)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "pan_cancer_with_sample_tumor")
  expect_true(r$COSMIC_PANCANCER_INCLUDES_SAMPLE_TUMOR)
  expect_equal(r$COSMIC_MATCHING_TUMOR_OCCURRENCES, 20)
})
test_that("5b. pan-cancer NOT including sample -> pan_cancer_without_sample_tumor", {
  r <- cls(lr(c("lung","adenocarcinoma","",40), c("large_intestine","adenocarcinoma","",60),
              c("pancreas","ductal_carcinoma","",30), c("breast","ductal_carcinoma","",15),
              c("ovary","serous_carcinoma","",10), c("kidney","clear_cell","",8)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "pan_cancer_without_sample_tumor")
  expect_false(r$COSMIC_PANCANCER_INCLUDES_SAMPLE_TUMOR)
})
test_that("6. COSMIC record without tumor info -> cosmic_tumor_type_missing (not evaluable)", {
  r <- cls(lr(c("","","",5)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "cosmic_tumor_type_missing")
  expect_false(r$COSMIC_CONTEXT_EVALUABLE)
  expect_true(is.na(r$COSMIC_TUMOR_SPECIFICITY_SCORE))
  expect_true(is.na(r$COSMIC_CONTEXT_SUPPORT_SCORE))
})
test_that("7. sample without tumor type -> tumor_type_unknown (not evaluable)", {
  r <- cls(lr(c("skin","malignant_melanoma","",40)), allowed = mapping_empty(), known = FALSE)
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "tumor_type_unknown")
  expect_false(r$COSMIC_CONTEXT_EVALUABLE)
  expect_true(is.na(r$COSMIC_CONTEXT_SUPPORT_SCORE))
})
test_that("8. unrecognized tumor type is not harmonized", {
  m <- load_tumor_type_mapping(list(cosmic = list(tumor_type_mapping = fixture("cosmic_tumor_type_mapping.tsv"))))
  expect_equal(nrow(harmonize_sample_tumor("NOT_A_REAL_TYPE", NA, m)), 0)
})
test_that("9. site matches, histology diverges -> compatible_match", {
  r <- cls(lr(c("skin","adenocarcinoma","",8)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "compatible_match")
})
test_that("10. multiple COSMIC contexts -> exact wins, occurrences split", {
  r <- cls(lr(c("skin","malignant_melanoma","",30), c("breast","ductal_carcinoma","",10)))
  expect_equal(r$COSMIC_TUMOR_CONTEXT_STATUS, "exact_match")
  expect_equal(r$COSMIC_MATCHING_TUMOR_OCCURRENCES, 30)
  expect_equal(r$COSMIC_OTHER_TUMOR_OCCURRENCES, 10)
})

# --- monotonicity & bounds (P3/P11) -----------------------------------------

test_that("support is strictly ordered exact > compatible > pan_with_sample > other", {
  s_exact <- cls(lr(c("skin","malignant_melanoma","",40)))$COSMIC_CONTEXT_SUPPORT_SCORE
  s_comp  <- cls(lr(c("skin","basal_cell_carcinoma","",40)))$COSMIC_CONTEXT_SUPPORT_SCORE
  s_pan   <- cls(lr(c("skin","malignant_melanoma","",20), c("lung","adenocarcinoma","",40),
                    c("large_intestine","adenocarcinoma","",60), c("pancreas","ductal_carcinoma","",30),
                    c("breast","ductal_carcinoma","",15), c("ovary","serous_carcinoma","",10)))$COSMIC_CONTEXT_SUPPORT_SCORE
  s_other <- cls(lr(c("breast","ductal_carcinoma","",40)))$COSMIC_CONTEXT_SUPPORT_SCORE
  expect_gt(s_exact, s_comp)
  expect_gt(s_comp, s_pan)
  expect_gt(s_pan, s_other)
})

test_that("high-recurrence other_tumor_only never outranks a low-recurrence exact_match", {
  s_other_10000 <- cls(lr(c("breast","ductal_carcinoma","",10000)))$COSMIC_CONTEXT_SUPPORT_SCORE
  s_exact_1     <- cls(lr(c("skin","malignant_melanoma","",1)))$COSMIC_CONTEXT_SUPPORT_SCORE
  expect_lt(s_other_10000, s_exact_1)
})

test_that("global and matched recurrence are computed separately", {
  r <- cls(lr(c("skin","malignant_melanoma","",10), c("breast","ductal_carcinoma","",990)))
  expect_gt(r$COSMIC_GLOBAL_RECURRENCE_SCORE, r$COSMIC_MATCHED_TUMOR_RECURRENCE_SCORE)
})

# --- normalized long-table storage: delimiters / missing fields (P6) --------

test_that("build_cosmic_long stores categories as columns (delimiter-safe) and sums counts", {
  rows <- data.frame(canonical_key = c("K1","K1","K1","K2"),
    tumor_site = c("skin; melanoma", "skin; melanoma", "lung", ""),
    tumor_histology = c("a=b", "a=b", "", ""),
    tumor_subtype = c("", "", "", ""),
    cosmic_count = c(3, 2, 5, NA), stringsAsFactors = FALSE)
  long <- build_cosmic_long(rows, "vX")
  k1 <- long[long$canonical_key == "K1" & long$cosmic_primary_site == "skin;_melanoma", ]
  expect_equal(k1$occurrence_count, 5)          # 3 + 2 summed, delimiters preserved as data
  expect_equal(k1$source_row_count, 2L)
  expect_equal(long$occurrence_count[long$canonical_key == "K2"], 1)  # NA count -> 1
})

# --- integration through processed DB + annotate_cosmic ---------------------

context_cfg <- function() {
  cache <- tempfile("cosmicctx")
  cfg <- list(input = list(genome_build = "GRCh38"),
              cosmic = list(raw_file = fixture("cosmic_raw_context.tsv"), release = "vCTX",
                            source_build = "GRCh38", target_build = "GRCh38", cache_dir = cache,
                            tumor_type_mapping = fixture("cosmic_tumor_type_mapping.tsv"),
                            pan_cancer_min_tumor_types = 5),
              reference = list(fasta = NULL))
  suppressWarnings(suppressMessages(prepare_cosmic_db(cfg)))
  cfg
}

test_that("annotate_cosmic stratifies end-to-end using the long table", {
  cfg <- context_cfg()
  x <- data.frame(sample_id = "S1",
    chrom = c("7","12","17"), pos = c(140453136, 25398284, 7674220),
    ref = c("A","C","C"), alt = c("T","A","T"),
    tumor_type = "MELANOMA", tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(annotate_cosmic(x, cfg, "GRCh38")))
  expect_equal(out$COSMIC_TUMOR_CONTEXT_STATUS,
               c("exact_match", "pan_cancer_with_sample_tumor", "other_tumor_only"))
  expect_equal(out$COSMIC_MATCHING_TUMOR_OCCURRENCES[1], 50)
  expect_equal(out$COSMIC_MATCHING_TUMOR_OCCURRENCES[3], 0)
  expect_gt(out$COSMIC_CONTEXT_SUPPORT_SCORE[1], out$COSMIC_CONTEXT_SUPPORT_SCORE[3])
})

run_with_tumor_type <- function(tt, cache_dir, release) {
  v <- data.frame(
    Tumor_Sample_Barcode = c("S1","S2","S3"),
    CHROM = c("7","12","17"), START = c(140453136, 25398284, 7674220),
    REF = c("A","C","C"), ALT = c("T","A","T"),
    Hugo_Symbol = c("BRAF","KRAS","TP53"), HGVSp_Short = c("p.V600E","p.G12V","p.R248Q"),
    Consequence = "missense_variant",
    DP = c(80,75,90), alt_count = c(38,35,45), AF = c(0.48,0.47,0.5),
    TLOD = c(80,78,85), MBQ = 35, MMQ = 60, FILTER = "PASS", gnomADg_AF = 0.0,
    stringsAsFactors = FALSE)
  vf <- tempfile(fileext = ".tsv"); write_tsv(v, vf)
  cfg <- list(
    input = list(vcf = vf, delimiter = "\t", genome_build = "GRCh38"),
    analysis = list(output_dir = tempfile("rt"), run_id = "t"),
    cancer = list(default_tumor_type = tt),
    hard_filters = list(enabled = TRUE, assay = "WGS"),
    technical_filters = list(min_depth = 20, min_alt_count_snv = 5, min_alt_count_indel = 8,
      min_af_snv = 0.03, min_af_indel = 0.05, min_tlod = 6, min_mbq = 25, min_mmq = 40),
    population_filters = list(common_af_threshold = 0.01, low_frequency_af_threshold = 0.001),
    cohort = list(recurrent_variant_fraction_artifact = 0.3, recurrent_locus_fraction_artifact = 0.4),
    guideline_classification = list(enabled = TRUE),
    cosmic = list(release = release, source_build = "GRCh38", cache_dir = cache_dir,
      tumor_type_mapping = fixture("cosmic_tumor_type_mapping.tsv"), pan_cancer_min_tumor_types = 5))
  res <- suppressWarnings(suppressMessages(run_tumor_only(cfg)))
  read_variants(file.path(res$run_dir, "tables", "variants_all.tsv.gz"), "\t")
}

test_that("tumor context NEVER changes filter_status, filter_reasons or CONFIDENCE_SCORE_BASE", {
  cfg <- context_cfg()
  cache <- cfg$cosmic$cache_dir; rel <- cfg$cosmic$release
  mel <- run_with_tumor_type("MELANOMA", cache, rel)
  brc <- run_with_tumor_type("BRCA", cache, rel)
  ord <- function(d) d[order(d$variant_id), ]
  mel <- ord(mel); brc <- ord(brc)
  expect_identical(mel$filter_status, brc$filter_status)
  expect_identical(mel$filter_reasons, brc$filter_reasons)
  expect_identical(mel$CONFIDENCE_SCORE_BASE, brc$CONFIDENCE_SCORE_BASE)
  # but the experimental COSMIC-context columns DO differ (BRAF exact in MELANOMA
  # vs other_tumor_only in BRCA)
  expect_false(identical(mel$COSMIC_TUMOR_CONTEXT_STATUS, brc$COSMIC_TUMOR_CONTEXT_STATUS))
})

test_that("weight-sensitivity analysis runs and reports churn", {
  cfg <- context_cfg()
  x <- data.frame(sample_id = "S1",
    chrom = c("7","12","17"), pos = c(140453136, 25398284, 7674220),
    ref = c("A","C","C"), alt = c("T","A","T"),
    tumor_type = "MELANOMA", tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(annotate_cosmic(x, cfg, "GRCh38")))
  s <- cosmic_weight_sensitivity(out)
  expect_true(all(c("20/80","40/60","50/50") %in% s$weightings))
  expect_true(is.numeric(s$n_category_changes_20_80_vs_50_50))
})
