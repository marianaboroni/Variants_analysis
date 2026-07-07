# Representative validation (P10). Uses an ILLUSTRATIVE, non-restricted fixture
# (not real COSMIC content) covering canonical driver scenarios. For validation
# against a real COSMIC release, point cosmic.raw_file at the downloaded TSV and
# run `tumoronly prepare-cosmic`; only aggregated results should be reported.

representative_cfg <- function() {
  cache <- tempfile("cosmicrep")
  cfg <- list(input = list(genome_build = "GRCh38"),
    cosmic = list(raw_file = fixture("cosmic_raw_representative.tsv"), release = "vREP",
                  source_build = "GRCh38", target_build = "GRCh38", cache_dir = cache,
                  tumor_type_mapping = fixture("cosmic_tumor_type_mapping.tsv"),
                  pan_cancer_min_tumor_types = 5),
    reference = list(fasta = NULL))
  suppressWarnings(suppressMessages(prepare_cosmic_db(cfg)))
  cfg
}

test_that("representative driver scenarios receive the expected context status", {
  cfg <- representative_cfg()
  x <- data.frame(
    sample_id = c("HGSOC1","UCEC1","COREAD1","MEL1","MEL2","MEL3"),
    gene = c("TP53","PIK3CA","KRAS","BRAF","ABCX","NOCTX"),
    chrom = c("17","3","12","7","1","2"),
    pos = c(7676154, 179234297, 25245350, 140753336, 100000, 200000),
    ref = c("C","A","C","A","G","T"), alt = c("T","G","T","T","A","C"),
    tumor_type = c("HGSOC","UCEC","COREAD","MELANOMA","MELANOMA","MELANOMA"),
    tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(annotate_cosmic(x, cfg, "GRCh38")))
  status <- setNames(out$COSMIC_TUMOR_CONTEXT_STATUS, out$gene)
  expect_equal(unname(status["TP53"]), "pan_cancer_with_sample_tumor")   # 6 sites incl ovary
  expect_equal(unname(status["PIK3CA"]), "exact_match")                  # 3 sites incl endometrium (<5)
  expect_equal(unname(status["KRAS"]), "pan_cancer_with_sample_tumor")   # 5 sites incl large_intestine
  expect_equal(unname(status["BRAF"]), "exact_match")                    # skin-dominant, 2 sites
  expect_equal(unname(status["ABCX"]), "other_tumor_only")               # only liver, sample melanoma
  expect_equal(unname(status["NOCTX"]), "cosmic_tumor_type_missing")     # no tumor context
  # aggregated-only reporting (no restricted content)
  agg <- as.data.frame(table(out$COSMIC_TUMOR_CONTEXT_STATUS))
  expect_true(all(agg$Freq >= 1))
})
