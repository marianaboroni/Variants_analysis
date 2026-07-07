# VCF/TSV input, dynamic VEP CSQ parsing (in io.R), and genome-build resolution.
# The efficient dynamic-CSQ reader lives in io.R (read_annotated_vcf,
# parse_vcf_annotation_format, add_info_annotations). This file adds robust,
# never-guess genome-build detection.

# Length of chr1/1 in each build; the most reliable single discriminator.
CHR1_LENGTH <- c(GRCh37 = 249250621L, GRCh38 = 248956422L)

#' Detect genome build from VCF header metadata lines.
#'
#' Inspects, in order: `##reference` line, then `##contig` length for chr1/1.
#' Returns `NA_character_` if it cannot be determined confidently (never guesses).
#'
#' @param meta character vector of VCF header lines (from [read_vcf_meta()]).
#' @return one of "GRCh37", "GRCh38", or NA.
#' @keywords internal
detect_genome_build <- function(meta) {
  if (length(meta) == 0) return(NA_character_)
  ref_line <- meta[grepl("^##reference", meta, ignore.case = TRUE)]
  if (length(ref_line) > 0) {
    rl <- tolower(paste(ref_line, collapse = " "))
    if (grepl("grch38|hg38|hs38|38", rl) && !grepl("grch37|hg19|37", rl)) return("GRCh38")
    if (grepl("grch37|hg19|hs37|37", rl) && !grepl("grch38|hg38|38", rl)) return("GRCh37")
  }
  contig1 <- meta[grepl("^##contig", meta) &
                    grepl("ID=(chr)?1[,>]", meta, ignore.case = TRUE)]
  if (length(contig1) > 0) {
    len <- suppressWarnings(as.integer(sub(".*length=([0-9]+).*", "\\1", contig1[[1]])))
    if (!is.na(len)) {
      hit <- names(CHR1_LENGTH)[match(len, CHR1_LENGTH)]
      if (!is.na(hit)) return(hit)
    }
  }
  NA_character_
}

#' Resolve the genome build for a run, enforcing the never-guess policy.
#'
#' Configuration (`input.genome_build`) prevails, but any conflict with the VCF
#' header raises a critical error. If neither config nor header determines the
#' build, execution stops with an actionable message.
#'
#' @param cfg resolved config list.
#' @param vcf_path path to the input VCF (used to read the header); may be a TSV,
#'   in which case only the config is used.
#' @return the resolved build string.
#' @keywords internal
resolve_genome_build <- function(cfg, vcf_path = NULL) {
  cfg_build <- cfg_get(cfg, c("input", "genome_build"), NULL)
  if (!is.null(cfg_build) && is.na(cfg_build)) cfg_build <- NULL

  header_build <- NA_character_
  if (!is.null(vcf_path) && file.exists(vcf_path) &&
      grepl("[.]vcf([.]gz)?$", vcf_path, ignore.case = TRUE)) {
    header_build <- detect_genome_build(read_vcf_meta(vcf_path))
  }

  if (!is.null(cfg_build) && !is.na(header_build) && cfg_build != header_build) {
    stop(sprintf(paste0(
      "Genome-build conflict: configuration says input.genome_build=%s but the ",
      "VCF header indicates %s.\nResolve the mismatch (fix the config or re-annotate ",
      "against the correct reference) before continuing."),
      cfg_build, header_build), call. = FALSE)
  }

  build <- cfg_build %||% (if (!is.na(header_build)) header_build else NULL)
  if (is.null(build)) {
    stop(paste0(
      "Could not determine whether the VCF is in GRCh37 or GRCh38.\n",
      "Set input.genome_build in the configuration file, e.g.:\n",
      "  input:\n    genome_build: GRCh38"), call. = FALSE)
  }
  if (!build %in% VALID_BUILDS) {
    stop(sprintf("Unsupported genome build '%s'. Expected one of: %s.",
                 build, paste(VALID_BUILDS, collapse = ", ")), call. = FALSE)
  }
  if (!is.na(header_build) && header_build == build) {
    log_step("build", "genome build confirmed from config and VCF header", build = build)
  } else {
    log_step("build", "genome build resolved", build = build,
             source = if (!is.null(cfg_build)) "config" else "vcf_header")
  }
  build
}
