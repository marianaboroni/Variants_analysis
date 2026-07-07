# Environment diagnostics: which R packages and external tools are available.

REQUIRED_PKGS <- c("data.table", "yaml", "jsonlite", "digest")
OPTIONAL_PKGS <- c("maftools", "rmarkdown", "ggplot2", "GenomicRanges",
                   "rtracklayer", "Biostrings", "httr2", "testthat")
EXTERNAL_TOOLS <- c("bcftools", "tabix", "CrossMap.py", "liftOver", "samtools")

#' Check the tumoronly installation and report capabilities.
#'
#' Reports which required/optional R packages and external tools are present and
#' what that enables (cross-build COSMIC liftover, REF validation, OncoKB API,
#' oncoplots, reporting). Never fails; returns a structured status.
#'
#' @param quiet suppress console output.
#' @return (invisibly) a list of `packages`, `tools`, `capabilities`, and `ok`.
#' @export
#' @examples
#' check_installation()
check_installation <- function(quiet = FALSE) {
  pkg_status <- function(p) {
    ok <- requireNamespace(p, quietly = TRUE)
    v <- if (ok) as.character(utils::packageVersion(p)) else NA_character_
    list(installed = ok, version = v)
  }
  packages <- stats::setNames(lapply(c(REQUIRED_PKGS, OPTIONAL_PKGS), pkg_status),
                              c(REQUIRED_PKGS, OPTIONAL_PKGS))
  tools <- stats::setNames(lapply(EXTERNAL_TOOLS, function(t) {
    unname(nzchar(Sys.which(t)))
  }), EXTERNAL_TOOLS)

  have <- function(p) isTRUE(packages[[p]]$installed)
  have_tool <- function(t) isTRUE(tools[[t]])
  capabilities <- list(
    run_core = all(vapply(REQUIRED_PKGS, have, logical(1))),
    same_build_cosmic = TRUE,
    cross_build_liftover = (have("rtracklayer") && have("GenomicRanges")) ||
      have_tool("CrossMap.py") || have_tool("liftOver"),
    ref_validation = have("Biostrings"),
    maf_and_oncoplot = have("maftools"),
    oncokb_api = have("httr2"),
    html_report = have("rmarkdown"))

  ok <- capabilities$run_core
  if (!quiet) {
    cat("tumoronly doctor\n================\n\nR packages:\n")
    for (p in names(packages)) {
      req <- if (p %in% REQUIRED_PKGS) "[required]" else "[optional]"
      mark <- if (packages[[p]]$installed) paste0("OK  ", packages[[p]]$version) else "MISSING"
      cat(sprintf("  %-14s %-10s %s\n", p, req, mark))
    }
    cat("\nExternal tools:\n")
    for (t in names(tools)) cat(sprintf("  %-14s %s\n", t, if (tools[[t]]) "found" else "not found"))
    cat("\nCapabilities:\n")
    for (c in names(capabilities)) cat(sprintf("  %-22s %s\n", c, if (isTRUE(capabilities[[c]])) "yes" else "no"))
    if (!capabilities$cross_build_liftover)
      cat("\nNote: cross-build COSMIC liftover is unavailable. Same-build preparation works.\n",
          "     Install rtracklayer (+ chain file), CrossMap, or UCSC liftOver to enable it.\n", sep = "")
    if (!ok) cat("\nERROR: required packages are missing; core pipeline cannot run.\n")
  }
  invisible(list(packages = packages, tools = tools, capabilities = capabilities, ok = ok))
}
