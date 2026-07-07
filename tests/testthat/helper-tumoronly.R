# Make the package functions available whether tests run under R CMD check
# (package loaded) or standalone via testthat::test_dir (dev mode: source R/).
if (!exists("run_tumor_only", mode = "function")) {
  # locate the package root by walking up from the test dir
  root <- normalizePath(file.path(getwd(), "..", ".."))
  r_dir <- file.path(root, "R")
  if (!dir.exists(r_dir)) r_dir <- file.path(getwd(), "R")  # when run from root
  files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  first <- file.path(r_dir, c("utils.R", "io.R", "config.R", "input.R"))
  files <- c(intersect(first, files), setdiff(files, first))
  for (f in files) sys.source(f, envir = globalenv())
  assign(".tumoronly_root", root, envir = globalenv())
}

fixture <- function(...) {
  root <- if (exists(".tumoronly_root", envir = globalenv())) get(".tumoronly_root", envir = globalenv())
          else normalizePath(file.path(getwd(), "..", ".."))
  candidates <- c(
    file.path(root, "tests", "fixtures", ...),
    file.path(getwd(), "..", "fixtures", ...),
    file.path(getwd(), "tests", "fixtures", ...))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) stop("fixture not found: ", paste(..., sep = "/"))
  normalizePath(hit[1])
}
