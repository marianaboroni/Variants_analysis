# Build the flat, denormalized COSMIC input TSV that `prepare_cosmic_db()`
# consumes as `cosmic.raw_file` (see R/cosmic.R::standardize_cosmic_raw()).
#
# COSMIC v104 ships mutation-occurrence data (Cosmic_GenomeScreensMutant /
# Cosmic_NonCodingVariants TSVs), genomic coordinates (their *_Normal.vcf,
# already 5'-shifted per the VCF standard - see docs/COSMIC_DECISION_RULES.md
# and the package README for why the Normal VCF, not the raw TSV's own
# GENOME_START/GENOMIC_WT_ALLELE/GENOMIC_MUT_ALLELE, is the coordinate source),
# and tumor classification (Cosmic_Classification, COSMIC_PHENOTYPE_ID ->
# PRIMARY_SITE/PRIMARY_HISTOLOGY) as SEPARATE files, joined only by ID columns
# (GENOMIC_MUTATION_ID, COSMIC_PHENOTYPE_ID). This function performs that join
# ONCE, producing one row per (mutation x sample occurrence) with VCF-standard
# coordinates and full tumor-context columns - exactly what `standardize_cosmic_raw()`
# expects, and what `build_cosmic_long()` needs to count occurrences per tumor
# category correctly (one row per real sample occurrence, not one per mutation).
#
# Cosmic_Sample_v104_GRCh38.tsv.gz (per-individual attributes: age, MSI, stage,
# ploidy, ...) is NOT integrated here - no part of `tumoronly`'s COSMIC engine
# currently consumes those fields. Left as a documented, deliberate limitation
# (see docs/AUXILIARY_DATABASES.md) rather than bolted on without a consumer.
#
# Column selection is used throughout (`fread(..., select = ...)`) because the
# real v104 files are tens of millions of rows / tens of GB; reading only the
# ~7 columns actually needed cuts I/O and memory by roughly 3-4x versus reading
# every column.

#' Build (or reuse) the flat COSMIC input TSV used as `cosmic.raw_file`.
#'
#' @param config path to a YAML config, or an already-read config list. Reads
#'   `cosmic.input_sources` (genome_screens_tsv, noncoding_tsv,
#'   genome_screens_normal_vcf, noncoding_normal_vcf, classification_tsv) and
#'   `cosmic.release`/`cosmic.cache_dir`.
#' @param force rebuild even if a valid cached input already exists.
#' @return a list with `raw_file_path` (set this as `cosmic.raw_file`),
#'   `manifest_path`, `manifest`, and `reused`.
#' @export
#' @examples
#' \dontrun{
#' build_cosmic_raw_input("config/example.yml")
#' }
build_cosmic_raw_input <- function(config, force = FALSE) {
  cfg <- if (is.character(config)) read_config(config) else config
  cfg <- resolve_config(cfg)
  src <- cfg_get(cfg, c("cosmic", "input_sources"), NULL)
  if (is.null(src)) stop("No `cosmic: input_sources:` section in configuration.", call. = FALSE)

  genome_screens_tsv        <- cfg_get(cfg, c("cosmic", "input_sources", "genome_screens_tsv"), NULL)
  noncoding_tsv              <- cfg_get(cfg, c("cosmic", "input_sources", "noncoding_tsv"), NULL)
  genome_screens_normal_vcf  <- cfg_get(cfg, c("cosmic", "input_sources", "genome_screens_normal_vcf"), NULL)
  noncoding_normal_vcf       <- cfg_get(cfg, c("cosmic", "input_sources", "noncoding_normal_vcf"), NULL)
  classification_tsv         <- cfg_get(cfg, c("cosmic", "input_sources", "classification_tsv"), NULL)
  release   <- as.character(cfg_get(cfg, c("cosmic", "release"), "vUNKNOWN"))
  cache_dir <- cfg_get(cfg, c("cosmic", "cache_dir"), "db/cosmic")

  required <- c(genome_screens_tsv = genome_screens_tsv, noncoding_tsv = noncoding_tsv,
                genome_screens_normal_vcf = genome_screens_normal_vcf,
                noncoding_normal_vcf = noncoding_normal_vcf,
                classification_tsv = classification_tsv)
  missing <- names(required)[vapply(required, function(p) is.null(p) || is.na(p) || !file.exists(p), logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("cosmic.input_sources missing or file not found for: %s", paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  out_dir       <- file.path(cache_dir, release, "raw_input")
  out_path      <- file.path(out_dir, "cosmic_flat_input.tsv.gz")
  manifest_path <- file.path(out_dir, "manifest.json")

  input_sha <- vapply(required, file_sha256, character(1))

  # ---- idempotent reuse -------------------------------------------------
  if (!force && file.exists(out_path) && file.exists(manifest_path)) {
    m <- tryCatch(jsonlite::read_json(manifest_path, simplifyVector = TRUE), error = function(e) NULL)
    if (!is.null(m) && identical(m$input_sha256, as.list(input_sha)) && identical(m$release, release)) {
      log_step("cosmic-input", "reusing cached COSMIC flat input (checksums match)",
               path = out_path, records = m$n_output)
      return(list(raw_file_path = out_path, manifest_path = manifest_path, manifest = m, reused = TRUE))
    }
    if (!force) {
      stop(sprintf(paste0(
        "A COSMIC flat input exists at %s but its manifest does not match the current inputs.\n",
        "Re-run with force = TRUE (CLI: --force) to rebuild it."), out_path), call. = FALSE)
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(file.path(out_dir, ".build.lock"))
  on.exit(release_lock(lock), add = TRUE)

  log_step("cosmic-input", "reading COSMIC Classification (COSMIC_PHENOTYPE_ID -> site/histology)",
           file = classification_tsv)
  classification <- data.table::fread(classification_tsv,
    select = c("COSMIC_PHENOTYPE_ID", "PRIMARY_SITE", "PRIMARY_HISTOLOGY", "HISTOLOGY_SUBTYPE_1"))
  classification <- unique(classification, by = "COSMIC_PHENOTYPE_ID")

  log_step("cosmic-input", "reading normalized genomic coordinates from *_Normal.vcf",
           coding = genome_screens_normal_vcf, noncoding = noncoding_normal_vcf)
  coords <- data.table::rbindlist(list(
    read_cosmic_normal_vcf(genome_screens_normal_vcf),
    read_cosmic_normal_vcf(noncoding_normal_vcf)
  ))
  coords <- unique(coords, by = "GENOMIC_MUTATION_ID")

  log_step("cosmic-input", "reading mutation-occurrence TSVs (coding + non-coding)",
           coding = genome_screens_tsv, noncoding = noncoding_tsv)
  mut_cols <- c("GENE_SYMBOL", "COSMIC_SAMPLE_ID", "COSMIC_PHENOTYPE_ID", "GENOMIC_MUTATION_ID",
                "LEGACY_MUTATION_ID", "COSMIC_STUDY_ID", "PUBMED_PMID", "MUTATION_SOMATIC_STATUS")
  mutations <- data.table::rbindlist(list(
    data.table::fread(genome_screens_tsv, select = mut_cols),
    data.table::fread(noncoding_tsv, select = mut_cols)
  ), use.names = TRUE, fill = TRUE)
  n_input <- nrow(mutations)

  # ---- join 1: mutation occurrence -> VCF-standard coordinates -----------
  mutations <- merge(mutations, coords, by = "GENOMIC_MUTATION_ID", all.x = TRUE)
  unmatched_coords <- mutations[is.na(CHROM)]
  mutations <- mutations[!is.na(CHROM)]

  # ---- join 2: -> tumor context (site/histology) --------------------------
  mutations <- merge(mutations, classification, by = "COSMIC_PHENOTYPE_ID", all.x = TRUE)

  # ---- rename to standardize_cosmic_raw()'s expected schema ---------------
  final <- mutations[, list(
    CHROMOSOME = CHROM, GENOME_START = POS,
    GENOMIC_WT_ALLELE = REF, GENOMIC_MUT_ALLELE = ALT,
    GENOMIC_MUTATION_ID = GENOMIC_MUTATION_ID, LEGACY_MUTATION_ID = LEGACY_MUTATION_ID,
    GENE_SYMBOL = GENE_SYMBOL, COSMIC_SAMPLE_ID = COSMIC_SAMPLE_ID,
    COSMIC_PHENOTYPE_ID = COSMIC_PHENOTYPE_ID, COSMIC_STUDY_ID = COSMIC_STUDY_ID,
    PUBMED_PMID = PUBMED_PMID, MUTATION_SOMATIC_STATUS = MUTATION_SOMATIC_STATUS,
    PRIMARY_SITE = PRIMARY_SITE, PRIMARY_HISTOLOGY = PRIMARY_HISTOLOGY,
    HISTOLOGY_SUBTYPE_1 = HISTOLOGY_SUBTYPE_1
  )]

  data.table::fwrite(final, out_path, sep = "\t", na = "NA", quote = FALSE, compress = "gzip")
  write_provenance(out_dir, "unmatched_to_normal_vcf", as.data.frame(unmatched_coords))

  manifest <- list(
    release = release,
    input_sources = as.list(required),
    input_sha256 = as.list(input_sha),
    build_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = tumoronly_version(),
    coordinate_source = "*_Normal.vcf (already 5'-shifted per VCF standard); mutation TSVs contribute only relational columns",
    n_input = n_input,
    n_unmatched_to_normal_vcf = nrow(unmatched_coords),
    n_output = nrow(final),
    n_classification_rows = nrow(classification),
    n_coords_rows = nrow(coords),
    cosmic_sample_tsv_integrated = FALSE
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  log_step("cosmic-input", "COSMIC flat input written", path = out_path,
           input = n_input, output = nrow(final), unmatched = nrow(unmatched_coords))

  list(raw_file_path = out_path, manifest_path = manifest_path, manifest = manifest, reused = FALSE)
}

#' Read CHROM/POS/GENOMIC_MUTATION_ID(ID)/REF/ALT from a COSMIC `*_Normal.vcf`
#' (no genotype columns; ID is the COSV genomic mutation identifier).
#' @keywords internal
read_cosmic_normal_vcf <- function(path) {
  dt <- data.table::fread(path, skip = "#CHROM", select = c(1, 2, 3, 4, 5),
                           col.names = c("CHROM", "POS", "GENOMIC_MUTATION_ID", "REF", "ALT"),
                           na.strings = c(".", ""))
  dt
}
