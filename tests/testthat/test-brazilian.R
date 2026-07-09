br_cfg <- function(escalate = TRUE) list(population = list(
  global = list(common_af = 0.01, rare_af = 0.001),
  brazilian = list(db = fixture("brazilian_db.tsv"), name = "ABraOM_test",
    genome_build = "GRCh38", common_af = 0.01, review_af = 0.001,
    minimum_allele_count = 2, cohort_size = 1171, escalate_to_review = escalate)))

mkv <- function(chrom, pos, ref, alt, gene, gnomad, vaf = 0.48) data.frame(
  chrom = chrom, pos = pos, ref = ref, alt = alt, gene = gene,
  max_pop_af = gnomad, vaf = vaf, dp = 80, stringsAsFactors = FALSE)

test_that("absent in gnomAD + present in Brazilian DB -> possible germline (NOT auto-somatic)", {
  x <- mkv("7", 140453136, "A", "T", "BRAF", NA)   # in brazilian_db (AF 0.012)
  out <- add_brazilian_population_evidence(x, br_cfg(), "GRCh38")
  expect_equal(out$POPULATION_EVIDENCE_STATUS, "absent_global_present_brazilian")
  expect_true(out$POSSIBLE_BRAZILIAN_GERMLINE)
  expect_equal(out$BRAZILIAN_AF, 0.012)
})

test_that("present globally, absent in Brazilian DB is distinguished", {
  x <- mkv("1", 999999, "G", "C", "XYZ", 0.05)     # common global, not in brazilian_db
  out <- add_brazilian_population_evidence(x, br_cfg(), "GRCh38")
  expect_equal(out$GLOBAL_POPULATION_STATUS, "common_global")
  expect_false(out$BRAZILIAN_DB_VARIANT_OBSERVED)
})

test_that("low Brazilian allele count is flagged as present_low_count / rare, not common", {
  x <- mkv("12", 25398284, "C", "A", "KRAS", NA)   # brazilian_db AF 0.0008 AC 2
  out <- add_brazilian_population_evidence(x, br_cfg(), "GRCh38")
  expect_true(out$BRAZILIAN_POPULATION_STATUS %in% c("rare_brazilian", "present_low_count_brazilian"))
  expect_equal(out$BRAZILIAN_AC, 2)
})

test_that("predisposition gene with heterozygous VAF is flagged for review", {
  x <- mkv("17", 7673803, "C", "T", "TP53", 0.0, vaf = 0.5)   # TP53 predisposition, het
  out <- add_brazilian_population_evidence(x, br_cfg(), "GRCh38")
  expect_true(out$POSSIBLE_BRAZILIAN_GERMLINE)
  expect_true(grepl("predisposition_gene_review|absent_global_present_brazilian", out$BRAZILIAN_GERMLINE_EVIDENCE))
})

test_that("no Brazilian DB configured -> brazilian status not_evaluable, absence not interpretable", {
  x <- mkv("5", 100, "A", "G", "GENE", NA)
  out <- add_brazilian_population_evidence(x, list(population = list()), "GRCh38")
  expect_equal(out$BRAZILIAN_POPULATION_STATUS, "not_evaluable")
  expect_false(out$BRAZILIAN_DB_ABSENCE_INTERPRETABLE)
})

test_that("documented REVIEW rule escalates PASS->REVIEW only when enabled", {
  x <- mkv("7", 140453136, "A", "T", "BRAF", NA)
  ev <- add_brazilian_population_evidence(x, br_cfg(TRUE), "GRCh38")
  ev$STATUS <- "PASS"; ev$filter_status <- "PASS"; ev$status_reason <- ""; ev$review_reasons <- ""
  on_rule <- apply_brazilian_population_rules(ev, br_cfg(TRUE))
  expect_equal(on_rule$STATUS, "REVIEW")
  expect_true(grepl("brazilian_population_database", on_rule$status_reason))
  # disabled -> invariant
  ev2 <- ev; ev2$STATUS <- "PASS"
  off_rule <- apply_brazilian_population_rules(ev2, br_cfg(FALSE))
  expect_equal(off_rule$STATUS, "PASS")
})

test_that("absence in checked databases is NOT labeled somatic", {
  # absent global (covered) + absent brazilian (interpretable) handled as absent_in_checked_databases
  x <- mkv("7", 140453200, "A", "T", "GENE", NA)   # not in brazilian_db
  out <- add_brazilian_population_evidence(x, br_cfg(), "GRCh38")
  expect_false(out$POPULATION_EVIDENCE_STATUS %in% c("common_global", "common_brazilian"))
  expect_true(out$POPULATION_EVIDENCE_STATUS %in%
    c("absent_in_checked_databases", "insufficient_coverage", "absent_global_present_brazilian"))
})
