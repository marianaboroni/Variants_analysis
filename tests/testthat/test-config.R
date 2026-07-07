test_that("config validation reports missing input with an actionable message", {
  expect_error(validate_config(list(analysis = list(output_dir = "x")), require_input = FALSE),
               "input.vcf", ignore.case = TRUE)
})

test_that("config validation rejects an invalid genome build", {
  cfg <- list(input = list(vcf = "x", genome_build = "hg19"),
              analysis = list(output_dir = "o"))
  expect_error(validate_config(cfg, require_input = FALSE), "GRCh3", ignore.case = TRUE)
})

test_that("resolve_config unifies legacy and new spellings", {
  cfg <- resolve_config(list(input = list(variants = "a.tsv"), output = list(dir = "res")))
  expect_equal(cfg$input$vcf, "a.tsv")
  expect_equal(cfg$analysis$output_dir, "res")
})

test_that("run id is deterministic from the input checksum", {
  tmp <- tempfile(fileext = ".tsv"); writeLines("x", tmp)
  cfg <- list(input = list(vcf = tmp), analysis = list())
  expect_equal(resolve_run_id(cfg), resolve_run_id(cfg))
})

test_that("a valid config passes", {
  tmp <- tempfile(fileext = ".tsv"); writeLines("x", tmp)
  cfg <- list(input = list(vcf = tmp, genome_build = "GRCh38"),
              analysis = list(output_dir = "o"))
  expect_invisible(validate_config(cfg))
})
