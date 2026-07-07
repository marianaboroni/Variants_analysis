# OS4/OM4 must derive ONLY from tumor-context recurrence, never from global-only
# recurrence. Same variant + same global count in different samples must give
# equal global evidence but different contextual evidence and OS4/OM4.

os4_cfg <- function() {
  cache <- tempfile("os4")
  cfg <- list(input = list(genome_build = "GRCh38"),
    cosmic = list(raw_file = fixture("cosmic_raw_representative.tsv"), release = "vOS4",
      source_build = "GRCh38", target_build = "GRCh38", cache_dir = cache,
      tumor_type_mapping = fixture("cosmic_tumor_type_mapping.tsv"),
      pan_cancer_min_tumor_types = 5, os4_min_matched_occurrences = 10,
      om4_min_matched_occurrences = 3),
    reference = list(fasta = NULL))
  suppressWarnings(suppressMessages(prepare_cosmic_db(cfg)))
  cfg
}

annotate_pair <- function(cfg, gene, chrom, pos, ref, alt, tt_match, tt_unrelated) {
  x <- data.frame(sample_id = c("m","u"), chrom = chrom, pos = pos, ref = ref, alt = alt,
                  tumor_type = c(tt_match, tt_unrelated), tumor_subtype = NA_character_,
                  stringsAsFactors = FALSE)
  suppressWarnings(suppressMessages(annotate_cosmic(x, cfg, "GRCh38")))
}

test_that("same variant, different samples: equal global, different context, OS4 only if matched (TP53)", {
  cfg <- os4_cfg()
  out <- annotate_pair(cfg, "TP53", "17", 7676154, "C", "T", "HGSOC", "MELANOMA")
  hgsoc <- out[1, ]; unrelated <- out[2, ]
  expect_equal(hgsoc$COSMIC_TOTAL_OCCURRENCES, unrelated$COSMIC_TOTAL_OCCURRENCES)   # global equal
  expect_gt(hgsoc$COSMIC_MATCHING_TUMOR_OCCURRENCES, unrelated$COSMIC_MATCHING_TUMOR_OCCURRENCES)
  expect_true(hgsoc$OS4_CONTEXTUAL)                 # ovary-matched, matching >= 10
  expect_false(unrelated$OS4_CONTEXTUAL)            # melanoma: not in sample tumor
  expect_false(unrelated$OM4_CONTEXTUAL)
  expect_true(unrelated$COSMIC_GLOBAL_EVIDENCE_ONLY)
})

test_that("BRAF: exact in melanoma grants OS4; absent in HGSOC does not", {
  cfg <- os4_cfg()
  out <- annotate_pair(cfg, "BRAF", "7", 140753336, "A", "T", "MELANOMA", "HGSOC")
  mel <- out[1, ]; hgsoc <- out[2, ]
  expect_equal(mel$COSMIC_TOTAL_OCCURRENCES, hgsoc$COSMIC_TOTAL_OCCURRENCES)
  expect_true(mel$OS4_CONTEXTUAL)
  expect_false(hgsoc$OS4_CONTEXTUAL)
  expect_false(hgsoc$OM4_CONTEXTUAL)
})

test_that("PIK3CA: exact in endometrium grants OS4; unrelated tumor does not", {
  cfg <- os4_cfg()
  out <- annotate_pair(cfg, "PIK3CA", "3", 179234297, "A", "G", "UCEC", "MELANOMA")
  ucec <- out[1, ]; mel <- out[2, ]
  expect_true(ucec$OS4_CONTEXTUAL)
  expect_false(mel$OS4_CONTEXTUAL)
  expect_false(mel$OM4_CONTEXTUAL)
})

test_that("unknown / missing context is not evaluable and grants no OS4/OM4", {
  cfg <- os4_cfg()
  # NOCTX: COSMIC record without tumor info
  noctx <- data.frame(sample_id = "s", chrom = "2", pos = 200000, ref = "T", alt = "C",
                      tumor_type = "MELANOMA", tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  o1 <- suppressWarnings(suppressMessages(annotate_cosmic(noctx, cfg, "GRCh38")))
  expect_false(o1$COSMIC_CONTEXT_EVALUABLE)
  expect_equal(o1$COSMIC_CONTEXT_INTERPRETATION, "not_evaluable")
  expect_false(o1$OS4_CONTEXTUAL); expect_false(o1$OM4_CONTEXTUAL)
  expect_true(o1$COSMIC_GLOBAL_EVIDENCE_ONLY)
  # unknown sample tumor type
  unk <- data.frame(sample_id = "s", chrom = "7", pos = 140753336, ref = "A", alt = "T",
                    tumor_type = "PANCANCER", tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  o2 <- suppressWarnings(suppressMessages(annotate_cosmic(unk, cfg, "GRCh38")))
  expect_equal(o2$COSMIC_TUMOR_CONTEXT_STATUS, "tumor_type_unknown")
  expect_false(o2$COSMIC_CONTEXT_EVALUABLE)
  expect_false(o2$OS4_CONTEXTUAL); expect_false(o2$OM4_CONTEXTUAL)
})

test_that("global-only recurrence does not change filter_status vs a matched context (regression)", {
  # legacy 'global' mode grants OS4 to discordant variants; contextual mode does not.
  cfg <- os4_cfg()
  x <- data.frame(sample_id = "s", chrom = "7", pos = 140753336, ref = "A", alt = "T",
                  tumor_type = "HGSOC", tumor_subtype = NA_character_, stringsAsFactors = FALSE)
  ctx <- suppressWarnings(suppressMessages(annotate_cosmic(x, cfg, "GRCh38")))
  cfg_g <- cfg; cfg_g$cosmic$os4om4_source <- "global"
  glb <- suppressWarnings(suppressMessages(annotate_cosmic(x, cfg_g, "GRCh38")))
  expect_false(ctx$OS4_CONTEXTUAL)   # contextual: BRAF not in ovary -> no OS4
  expect_true(glb$OS4_CONTEXTUAL)    # legacy global: would have granted OS4
})
