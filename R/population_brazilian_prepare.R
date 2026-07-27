# Preparation of a local, VCF-standard-normalized Brazilian population database
# (e.g. ABraOM/SABE) from a raw ANNOVAR-format TSV.
#
# ANNOVAR-exported population databases represent indels without a VCF anchor
# base. From real ABraOM rows cross-checked against Ensembl GRCh38 coordinates
# (rs56289060: chr1 10434-10436 CCC>CCCC; rs749280879: chr1 51865-51877 poly-A):
#   deletion  (Alt == "-"): Start = first deleted base.
#              VCF: pos = Start - 1, ref = anchor + Ref, alt = anchor.
#   insertion (Ref == "-"): Start = the anchor position itself (no shift).
#              VCF: pos = Start,     ref = anchor,       alt = anchor + Alt.
# This asymmetry is real (confirmed against two independent rsID lookups, not
# assumed) - deletions need one more backward step to reach the anchor,
# insertions do not. It is never trusted blindly: every reconstructed REF is
# re-fetched from the FASTA and compared; mismatches go to a provenance table
# instead of the database (see reconstruct_vcf_representation()).
#
# Preparation happens ONCE here; `population.brazilian.db` in later configs
# points at the output, matched by canonical key BUILD|CHROM|POS|REF|ALT
# (see R/population_brazilian.R::load_brazilian_db()).

#' Prepare (or reuse) a local, VCF-standard-normalized Brazilian population
#' database (ABraOM/BIPMed/... ANNOVAR-format TSV) for use as
#' `population.brazilian.db`.
#'
#' @param config path to a YAML config, or an already-read config list.
#' @param force rebuild even if a valid cached database exists.
#' @return a list with `db_path`, `manifest_path`, `manifest`, and `reused`.
#' @export
#' @examples
#' \dontrun{
#' prepare_brazilian_db("config/example.yml")
#' }
prepare_brazilian_db <- function(config, force = FALSE) {
  cfg <- if (is.character(config)) read_config(config) else config
  cfg <- resolve_config(cfg)
  bz <- cfg_get(cfg, c("population", "brazilian"), NULL)
  if (is.null(bz)) stop("No `population: brazilian:` section in configuration.", call. = FALSE)

  raw_file  <- cfg_get(cfg, c("population", "brazilian", "raw_file"), NULL)
  build     <- as.character(cfg_get(cfg, c("population", "brazilian", "genome_build"), "GRCh38"))
  name      <- as.character(cfg_get(cfg, c("population", "brazilian", "name"), "brazilian_db"))
  cache_dir <- cfg_get(cfg, c("population", "brazilian", "cache_dir"), "db/abraom")
  fasta     <- cfg_get(cfg, c("reference", "fasta"), NULL)

  if (is.null(raw_file) || is.na(raw_file))
    stop("population.brazilian.raw_file is required to prepare the Brazilian population database.", call. = FALSE)
  if (!file.exists(raw_file))
    stop(sprintf("population.brazilian.raw_file not found: %s", raw_file), call. = FALSE)
  if (is.null(fasta) || is.na(fasta) || !file.exists(fasta))
    stop(paste0("reference.fasta is required to prepare the Brazilian population database ",
                "(indel positions are reconstructed against the reference)."), call. = FALSE)
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required. Install with BiocManager::install('Biostrings').", call. = FALSE)

  processed_dir <- file.path(cache_dir, paste0(name, "_", build))
  db_path       <- file.path(processed_dir, "abraom_db.tsv.gz")
  manifest_path <- file.path(processed_dir, "manifest.json")

  raw_sha   <- file_sha256(raw_file)
  fasta_sha <- file_sha256(fasta)

  # ---- idempotent reuse -------------------------------------------------
  if (!force && file.exists(db_path) && file.exists(manifest_path)) {
    m <- tryCatch(jsonlite::read_json(manifest_path, simplifyVector = TRUE), error = function(e) NULL)
    norm <- function(z) if (is.null(z) || length(z) == 0 || is.na(z)) NA_character_ else as.character(z)
    if (!is.null(m) && identical(norm(m$raw_sha256), norm(raw_sha)) &&
        identical(norm(m$fasta_sha256), norm(fasta_sha)) &&
        identical(norm(m$genome_build), norm(build))) {
      log_step("abraom", "reusing cached Brazilian population database (checksums match)",
               db = db_path, records = m$n_output)
      return(list(db_path = db_path, manifest_path = manifest_path, manifest = m, reused = TRUE))
    }
    if (!force) {
      stop(sprintf(paste0(
        "A processed Brazilian population database exists at %s but its manifest does not match ",
        "the current inputs (raw file, FASTA, or build changed).\n",
        "Re-run with force = TRUE (CLI: --force) to rebuild it."), processed_dir), call. = FALSE)
    }
  }

  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(file.path(processed_dir, ".build.lock"))
  on.exit(release_lock(lock), add = TRUE)

  log_step("abraom", "reading raw Brazilian population TSV (ANNOVAR format)", file = raw_file)
  raw <- data.table::fread(raw_file, sep = "\t", na.strings = c("", ".", "NA"))
  n_input <- nrow(raw)

  std <- standardize_brazilian_raw(raw)
  essential_ok <- !is.na(std$chrom) & !is.na(std$pos) & !is.na(std$ref) & !is.na(std$alt) &
    std$ref != "" & std$alt != ""
  dropped_missing <- std[!essential_ok]
  std <- std[essential_ok]
  std[, chrom := normalize_chrom(chrom)]

  log_step("abraom", "loading reference FASTA", file = fasta)
  fa <- Biostrings::readDNAStringSet(fasta)
  names(fa) <- normalize_chrom(sub("\\s.*$", "", names(fa)))

  log_step("abraom", "reconstructing VCF-standard representation from ANNOVAR indel notation", n = nrow(std))
  rec <- reconstruct_vcf_representation(std, fa)
  converted    <- rec$converted
  ref_mismatch <- rec$ref_mismatch
  unresolved   <- rec$unresolved

  converted[, canonical_key := canonical_variant_key(build, chrom, pos, ref, alt)]
  dup_flag <- duplicated(converted$canonical_key)
  duplicates <- converted[dup_flag]
  converted <- converted[!dup_flag]

  out <- converted[, list(
    CHROM = chrom, POS = pos, REF = ref, ALT = alt,
    AF = af, AC = ac, AN = an, HOM = hom,
    Variant_Type = variant_type, rsID = rsid, Source_FILTER = filter_flag,
    Original_CHROM = orig_chrom, Original_POS = orig_pos,
    Original_REF = orig_ref, Original_ALT = orig_alt
  )]

  write_tsv(as.data.frame(out), db_path)
  write_provenance(processed_dir, "ref_mismatch", as.data.frame(ref_mismatch))
  write_provenance(processed_dir, "unresolved", as.data.frame(unresolved))
  write_provenance(processed_dir, "dropped_missing", as.data.frame(dropped_missing))
  write_provenance(processed_dir, "duplicates", as.data.frame(duplicates))

  manifest <- list(
    name = name, genome_build = build,
    raw_file = normalizePath(raw_file, mustWork = FALSE), raw_sha256 = raw_sha,
    fasta_file = normalizePath(fasta, mustWork = FALSE), fasta_sha256 = fasta_sha,
    build_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = tumoronly_version(),
    indel_anchor_convention = paste(
      "deletion: pos=Start-1, ref=anchor+Ref, alt=anchor;",
      "insertion: pos=Start, ref=anchor, alt=anchor+Alt",
      "(empirically confirmed against Ensembl GRCh38 coordinates for rs56289060 and rs749280879,",
      "see R/population_brazilian_prepare.R header)"),
    n_input = n_input,
    n_dropped_missing = nrow(dropped_missing),
    n_snv = sum(out$Variant_Type == "SNV"),
    n_mnp = sum(out$Variant_Type == "MNP"),
    n_insertion = sum(out$Variant_Type == "insertion"),
    n_deletion = sum(out$Variant_Type == "deletion"),
    n_ref_mismatch = nrow(ref_mismatch),
    n_unresolved = nrow(unresolved),
    n_duplicates = nrow(duplicates),
    n_output = nrow(out),
    ref_validated = TRUE
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  log_step("abraom", "processed Brazilian population database written", db = db_path,
           input = n_input, output = nrow(out), ref_mismatch = nrow(ref_mismatch),
           unresolved = nrow(unresolved))

  list(db_path = db_path, manifest_path = manifest_path, manifest = manifest, reused = FALSE)
}

#' Standardize a raw ANNOVAR-format population TSV (ABraOM/BIPMed/...) into
#' chrom/pos/ref/alt/af/ac/an/hom/rsid/filter_flag. Flexible to column-name
#' variants across ANNOVAR-derived exports.
#' @keywords internal
standardize_brazilian_raw <- function(x) {
  x <- data.table::as.data.table(x)
  data.table::data.table(
    chrom  = as.character(coalesce_columns(x, c("Chr", "CHROM", "Chromosome", "chr", "chrom"))),
    pos    = to_numeric_safe(coalesce_columns(x, c("Start", "POS", "pos"))),
    ref    = as.character(coalesce_columns(x, c("Ref", "REF", "ref"))),
    alt    = as.character(coalesce_columns(x, c("Alt", "ALT", "alt"))),
    an     = to_numeric_safe(coalesce_columns(x, c("Allele_number", "AN", "an"))),
    ac     = to_numeric_safe(coalesce_columns(x, c("Allele_ALT_count", "AC", "ac"))),
    af     = to_numeric_safe(coalesce_columns(x, c("Frequencies", "AF", "af"))),
    hom    = to_numeric_safe(coalesce_columns(x, c("HomozygousALT_count", "HOM", "hom"))),
    rsid   = as.character(coalesce_columns(x, c("avsnp150", "rsID", "rsid"))),
    filter_flag = as.character(coalesce_columns(x, c("FILTER", "filter")))
  )
}

#' Reconstruct VCF-standard CHROM/POS/REF/ALT from ANNOVAR-format rows,
#' anchoring indels against the reference FASTA, then re-validating every
#' resulting REF against the same FASTA (never trusted silently). Rows whose
#' chromosome is absent from the FASTA, or whose reconstructed REF does not
#' match the reference, are excluded from `converted` and reported separately.
#'
#' @param std data.table with chrom/pos/ref/alt (+ af/ac/an/hom/rsid/filter_flag).
#' @param fa a named `DNAStringSet` (names = normalized chrom) from the target FASTA.
#' @return list(converted, ref_mismatch, unresolved).
#' @keywords internal
reconstruct_vcf_representation <- function(std, fa) {
  std <- data.table::copy(std)
  std[, `:=`(orig_chrom = chrom, orig_pos = pos, orig_ref = ref, orig_alt = alt)]

  is_del <- std$alt == "-"
  is_ins <- std$ref == "-"
  std[, variant_type := ifelse(is_del, "deletion",
                        ifelse(is_ins, "insertion",
                        ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "SNV", "MNP")))]

  # anchor position: deletion needs the base BEFORE the first deleted base
  # (Start - 1); insertion's Start already IS the anchor (no shift) - see
  # file header for the empirical derivation.
  std[, anchor_pos := data.table::fifelse(is_del, pos - 1L,
                      data.table::fifelse(is_ins, pos, NA_integer_))]

  unresolved <- std[(is_del | is_ins) & !(chrom %in% names(fa))]
  std <- std[!((is_del | is_ins) & !(chrom %in% names(fa)))]
  is_del <- std$alt == "-"; is_ins <- std$ref == "-"

  anchor_base <- rep(NA_character_, nrow(std))
  idx_indel <- which(is_del | is_ins)
  for (chr in unique(std$chrom[idx_indel])) {
    sel <- idx_indel[std$chrom[idx_indel] == chr]
    anchor_base[sel] <- toupper(as.character(
      Biostrings::Views(fa[[chr]], start = std$anchor_pos[sel], width = 1L)))
  }
  std[, anchor_base := anchor_base]

  std[is_del, `:=`(pos = anchor_pos, ref = paste0(anchor_base, ref), alt = anchor_base)]
  std[is_ins, `:=`(pos = anchor_pos, ref = anchor_base, alt = paste0(anchor_base, alt))]
  std[, `:=`(anchor_pos = NULL, anchor_base = NULL)]

  # ---- self-validation: re-fetch REF from the FASTA at (chrom,pos,width)
  # and compare. Meaningful for SNV/MNP/deletion (multi-base match); for pure
  # insertions it only catches out-of-bounds/chrom errors, not a systematic
  # off-by-one - that was ruled out empirically (see file header), not by this
  # check alone. Never trusted implicitly: mismatches are excluded and kept.
  ok <- rep(FALSE, nrow(std))
  for (chr in unique(std$chrom)) {
    sel <- which(std$chrom == chr)
    if (!chr %in% names(fa)) next
    seqs <- toupper(as.character(
      Biostrings::Views(fa[[chr]], start = std$pos[sel], width = nchar(std$ref[sel]))))
    ok[sel] <- seqs == toupper(std$ref[sel])
  }
  list(converted = std[ok], ref_mismatch = std[!ok], unresolved = unresolved)
}
