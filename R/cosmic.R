# Idempotent preparation of a local, coordinate-based COSMIC database.
#
# Liftover and normalization happen ONCE here, never per sample. The processed
# database is keyed on the canonical variant key BUILD|CHROM|POS|REF|ALT and is
# reused across analyses via a checksummed manifest. No variant is ever dropped
# silently: every excluded record is written to a provenance table.

#' Prepare (or reuse) the local processed COSMIC database.
#'
#' Reads the raw COSMIC TSV declared in the config, validates its build, lifts it
#' over to the target build if needed (via a pluggable backend), normalizes
#' chromosome names and alleles, validates REF against the target FASTA when
#' available, builds the canonical key, deduplicates, and writes the processed
#' database plus a manifest and provenance tables. Idempotent: if a valid,
#' checksum-matching database already exists it is reused without re-running the
#' liftover, unless `force = TRUE`.
#'
#' @param config path to a YAML config, or an already-read config list.
#' @param force rebuild even if a valid cached database exists.
#' @return a list with `db_path`, `manifest_path`, `manifest`, and `reused`.
#' @export
#' @examples
#' \dontrun{
#' prepare_cosmic_db("config/example.yml")
#' }
prepare_cosmic_db <- function(config, force = FALSE) {
  cfg <- if (is.character(config)) read_config(config) else config
  cfg <- resolve_config(cfg)
  cs <- cfg_get(cfg, "cosmic", NULL)
  if (is.null(cs)) stop("No `cosmic:` section in configuration.", call. = FALSE)

  raw_file    <- cfg_get(cfg, c("cosmic", "raw_file"), NULL)
  release     <- as.character(cfg_get(cfg, c("cosmic", "release"), "vUNKNOWN"))
  source_build<- as.character(cfg_get(cfg, c("cosmic", "source_build"), "GRCh38"))
  target_build<- as.character(cfg_get(cfg, c("cosmic", "target_build"),
                                      cfg_get(cfg, c("input", "genome_build"), "GRCh38")))
  cache_dir   <- cfg_get(cfg, c("cosmic", "cache_dir"), "db/cosmic")
  chain_file  <- cfg_get(cfg, c("cosmic", "chain_file"), cfg_get(cfg, c("reference", "chain_file"), NULL))
  fasta       <- cfg_get(cfg, c("reference", "fasta"), NULL)

  if (is.null(raw_file) || is.na(raw_file))
    stop("cosmic.raw_file is required to prepare the COSMIC database.", call. = FALSE)
  if (!file.exists(raw_file))
    stop(sprintf("cosmic.raw_file not found: %s", raw_file), call. = FALSE)

  processed_dir <- file.path(cache_dir, release, sprintf("%s_to_%s", source_build, target_build))
  db_path       <- file.path(processed_dir, "cosmic_db.rds")
  manifest_path <- file.path(processed_dir, "manifest.json")

  raw_sha   <- file_sha256(raw_file)
  chain_sha <- file_sha256(chain_file)
  fasta_sha <- file_sha256(fasta)

  # ---- idempotent reuse -------------------------------------------------
  if (!force && file.exists(db_path) && file.exists(manifest_path)) {
    m <- tryCatch(jsonlite::read_json(manifest_path, simplifyVector = TRUE),
                  error = function(e) NULL)
    norm <- function(z) if (is.null(z) || length(z) == 0 || is.na(z)) NA_character_ else as.character(z)
    if (!is.null(m) && identical(norm(m$raw_sha256), norm(raw_sha)) &&
        identical(norm(m$chain_sha256), norm(chain_sha)) &&
        identical(norm(m$source_build), norm(source_build)) &&
        identical(norm(m$target_build), norm(target_build)) &&
        identical(norm(m$release), norm(release))) {
      log_step("cosmic", "reusing cached COSMIC database (checksums match)",
               db = db_path, records = m$n_converted)
      return(list(db_path = db_path, manifest_path = manifest_path,
                  manifest = m, reused = TRUE))
    }
    if (!force) {
      stop(sprintf(paste0(
        "A processed COSMIC database exists at %s but its manifest does not match ",
        "the current inputs (raw file, chain, or build changed).\n",
        "Re-run with force = TRUE (CLI: --force) to rebuild it."), processed_dir),
        call. = FALSE)
    }
  }

  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  lock <- acquire_lock(file.path(processed_dir, ".build.lock"))
  on.exit(release_lock(lock), add = TRUE)

  log_step("cosmic", "reading raw COSMIC file", file = raw_file, release = release)
  raw <- read_variants(raw_file, "\t")
  n_input <- nrow(raw)

  std <- standardize_cosmic_raw(raw)
  # provenance: records dropped for missing essential fields
  essential_ok <- !is.na(std$chrom) & !is.na(std$pos) & !is.na(std$ref) & !is.na(std$alt) &
    std$ref != "" & std$alt != ""
  dropped_missing <- std[!essential_ok, , drop = FALSE]
  std <- std[essential_ok, , drop = FALSE]
  std$chrom <- normalize_chrom(std$chrom)

  # split multiallelic before liftover/normalization
  std <- split_multiallelic(std, "alt")

  # ---- liftover (only when builds differ) -------------------------------
  lifted <- run_cosmic_liftover(std, source_build, target_build, chain_file)
  converted <- lifted$converted
  unmapped  <- lifted$unmapped
  multi     <- lifted$multi_mapped

  # ---- normalize alleles (trim + left-align) ----------------------------
  na <- normalize_alleles(converted$pos, converted$ref, converted$alt)
  converted$pos <- na$pos; converted$ref <- na$ref; converted$alt <- na$alt

  # ---- validate REF against target FASTA (never silently drop) ----------
  ref_mismatch <- converted[0, , drop = FALSE]
  vres <- validate_ref_against_fasta(converted, fasta, target_build)
  if (!is.null(vres)) {
    ref_mismatch <- converted[!vres, , drop = FALSE]
    converted <- converted[vres, , drop = FALSE]
  }

  # ---- canonical key + dedup --------------------------------------------
  converted$canonical_key <- canonical_variant_key(
    target_build, converted$chrom, converted$pos, converted$ref, converted$alt)
  dup_flag <- duplicated(converted$canonical_key)
  duplicates <- converted[dup_flag, , drop = FALSE]
  db <- collapse_cosmic_by_key(converted)
  db$COSMIC_RELEASE <- release
  db$COSMIC_TARGET_BUILD <- target_build

  # normalized long context table (P6): the PRIMARY analytical structure.
  long <- build_cosmic_long(converted, release)

  # ---- write processed db + long table + provenance ---------------------
  dt <- data.table::as.data.table(db)
  data.table::setkey(dt, canonical_key)
  saveRDS(dt, db_path, compress = "gzip")
  write_tsv(as.data.frame(dt), file.path(processed_dir, "cosmic_db.tsv.gz"))
  lt <- data.table::as.data.table(long)
  data.table::setkey(lt, canonical_key)
  saveRDS(lt, file.path(processed_dir, "cosmic_context_long.rds"), compress = "gzip")
  write_tsv(as.data.frame(lt), file.path(processed_dir, "cosmic_context_long.tsv.gz"))

  write_provenance(processed_dir, "converted", db)
  write_provenance(processed_dir, "unmapped", unmapped)
  write_provenance(processed_dir, "multi_mapped", multi)
  write_provenance(processed_dir, "ref_mismatch", ref_mismatch)
  write_provenance(processed_dir, "dropped_missing", dropped_missing)
  write_provenance(processed_dir, "duplicates", duplicates)

  manifest <- list(
    release = release,
    raw_file = normalizePath(raw_file, mustWork = FALSE),
    raw_sha256 = raw_sha,
    source_build = source_build,
    target_build = target_build,
    chain_file = if (is.null(chain_file)) NA_character_ else normalizePath(chain_file, mustWork = FALSE),
    chain_sha256 = chain_sha %||% NA_character_,
    fasta_file = if (is.null(fasta)) NA_character_ else normalizePath(fasta, mustWork = FALSE),
    fasta_sha256 = fasta_sha %||% NA_character_,
    build_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = tumoronly_version(),
    liftover_backend = lifted$backend,
    tool_versions = lifted$tool_versions,
    n_input = n_input,
    n_dropped_missing = nrow(dropped_missing),
    n_converted = nrow(db),
    n_unmapped = nrow(unmapped),
    n_multi_mapped = nrow(multi),
    n_ref_mismatch = nrow(ref_mismatch),
    n_duplicates = nrow(duplicates),
    ref_validated = !is.null(vres),
    schema = names(db),
    context_long_schema = names(long),
    n_context_long_rows = nrow(long)
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  log_step("cosmic", "processed COSMIC database written", db = db_path,
           input = n_input, converted = nrow(db), unmapped = nrow(unmapped),
           ref_mismatch = nrow(ref_mismatch), duplicates = nrow(duplicates))

  list(db_path = db_path, manifest_path = manifest_path, manifest = manifest, reused = FALSE)
}

#' Standardize raw COSMIC columns into chrom/pos/ref/alt/id/count/site/gene.
#' Flexible to COSMIC's varying release schemas.
#' @keywords internal
standardize_cosmic_raw <- function(x) {
  data.frame(
    chrom = as.character(coalesce_columns(x, c("CHROMOSOME", "Chromosome", "CHROM", "chr", "chrom"))),
    pos   = to_numeric_safe(coalesce_columns(x, c("GENOME_START", "GENOMIC_START", "Start_Position", "START", "POS", "pos"))),
    ref   = as.character(coalesce_columns(x, c("GENOMIC_WT_ALLELE", "WT_ALLELE", "REF", "Reference_Allele", "ref"))),
    alt   = as.character(coalesce_columns(x, c("GENOMIC_MUT_ALLELE", "MUT_ALLELE", "ALT", "Tumor_Seq_Allele2", "alt"))),
    cosmic_id = as.character(coalesce_columns(x, c("GENOMIC_MUTATION_ID", "COSMIC_MUTATION_ID", "LEGACY_MUTATION_ID", "COSMIC_ID", "Existing_variation"))),
    cosmic_count = to_numeric_safe(coalesce_columns(x, c("COSMIC_COUNT", "CNT", "Count", "MUTATION_COUNT"))),
    gene  = as.character(coalesce_columns(x, c("GENE_SYMBOL", "GENE_NAME", "Hugo_Symbol", "SYMBOL", "Gene", "gene"))),
    tumor_site = as.character(coalesce_columns(x, c("PRIMARY_SITE", "Primary site", "primary_site", "TUMOUR_SITE"))),
    tumor_histology = as.character(coalesce_columns(x, c("PRIMARY_HISTOLOGY", "Primary histology", "HISTOLOGY", "histology"))),
    tumor_subtype = as.character(coalesce_columns(x, c("HISTOLOGY_SUBTYPE_1", "HISTOLOGY_SUBTYPE", "PRIMARY_HISTOLOGY_SUBTYPE", "subtype"))),
    stringsAsFactors = FALSE
  )
}

#' Build the normalized long-format COSMIC context table: one row per
#' (canonical_key, site, histology, subtype) with summed occurrences and the
#' number of source rows collapsed. This is the PRIMARY analytical structure
#' for tumor-type-stratified matching (see R/cosmic_context.R). Delimiters are
#' irrelevant here because categories are stored as columns, not concatenated.
#' @keywords internal
build_cosmic_long <- function(rows, release) {
  site <- norm_cosmic_cat(rows$tumor_site)
  hist <- norm_cosmic_cat(rows$tumor_histology)
  sub  <- norm_cosmic_cat(rows$tumor_subtype)
  cnt  <- rows$cosmic_count; cnt[is.na(cnt)] <- 1
  grp  <- paste(rows$canonical_key, site, hist, sub, sep = "\r")
  occ  <- tapply(cnt, grp, sum)
  srcn <- tapply(cnt, grp, length)
  first <- !duplicated(grp)
  data.frame(
    canonical_key = rows$canonical_key[first],
    cosmic_primary_site = site[first],
    cosmic_histology = hist[first],
    cosmic_subtype = sub[first],
    occurrence_count = as.numeric(occ[grp[first]]),
    cosmic_release = release,
    source_row_count = as.integer(srcn[grp[first]]),
    stringsAsFactors = FALSE)
}

#' Normalize a COSMIC tumor category (site/histology) to a lowercase key.
#' @keywords internal
norm_cosmic_cat <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[is.na(x) | x %in% c("", "ns", "na", "none", "unknown")] <- ""
  gsub("[[:space:]]+", "_", x)
}

#' Collapse rows sharing a canonical key, aggregating IDs/counts and building a
#' per-tumor-category occurrence breakdown so each key appears once. Preserves
#' original COSMIC identifiers (deduplicated). The breakdown string encodes
#' "site::histology=count" pairs joined by ";", enabling tumor-type-stratified
#' evidence downstream (see R/cosmic_context.R). Counts with no tumor category
#' are aggregated under the empty category.
#' @keywords internal
collapse_cosmic_by_key <- function(rows) {
  key <- rows$canonical_key
  rows$cat <- paste(norm_cosmic_cat(rows$tumor_site), norm_cosmic_cat(rows$tumor_histology), sep = "::")
  agg_ids  <- tapply(rows$cosmic_id, key, function(z) paste(unique(stats::na.omit(z)), collapse = ";"))
  agg_cnt  <- tapply(rows$cosmic_count, key, function(z) if (all(is.na(z))) NA_real_ else sum(z, na.rm = TRUE))
  agg_site <- tapply(rows$tumor_site, key, function(z) paste(unique(stats::na.omit(z[z != ""])), collapse = ";"))
  # per-category breakdown: sum counts within each (site::histology) per key
  safe <- function(s) gsub("[;=]", "_", s)   # keep the export string parseable
  agg_bd <- tapply(seq_len(nrow(rows)), key, function(idx) {
    cats <- safe(rows$cat[idx]); cnt <- rows$cosmic_count[idx]
    cnt[is.na(cnt)] <- 1  # a record with no count still counts as one occurrence
    by_cat <- tapply(cnt, cats, sum)
    paste(sprintf("%s=%s", names(by_cat), format(as.numeric(by_cat), trim = TRUE, scientific = FALSE)),
          collapse = ";")
  })
  first <- rows[!duplicated(key), , drop = FALSE]
  k <- first$canonical_key
  data.frame(
    canonical_key = k,
    chrom = first$chrom, pos = first$pos, ref = first$ref, alt = first$alt,
    gene = first$gene,
    COSMIC_MUTATION_IDS = as.character(agg_ids[k]),
    COSMIC_OCCURRENCE_COUNT = as.numeric(agg_cnt[k]),
    COSMIC_TUMOR_TYPES = as.character(agg_site[k]),
    COSMIC_TUMOR_BREAKDOWN = as.character(agg_bd[k]),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
write_provenance <- function(dir, name, df) {
  path <- file.path(dir, paste0(name, ".tsv.gz"))
  if (is.null(df) || nrow(df) == 0) {
    write_tsv(data.frame(note = sprintf("no %s records", name)), path)
  } else {
    write_tsv(as.data.frame(df), path)
  }
  invisible(path)
}

#' Run (or skip) COSMIC liftover between builds through a pluggable backend.
#'
#' If source and target builds are identical, this is a passthrough (no liftover).
#' For a cross-build conversion it requires one of: rtracklayer (+ chain file),
#' CrossMap, or UCSC liftOver. If none is available it stops with an actionable
#' message - it never silently drops or leaves variants unconverted.
#'
#' @return list(converted, unmapped, multi_mapped, backend, tool_versions).
#' @keywords internal
run_cosmic_liftover <- function(std, source_build, target_build, chain_file) {
  empty <- std[0, , drop = FALSE]
  if (identical(source_build, target_build)) {
    log_step("cosmic", "no liftover needed (source build == target build)",
             build = target_build)
    return(list(converted = std, unmapped = empty, multi_mapped = empty,
                backend = "none (same build)", tool_versions = list()))
  }

  if (requireNamespace("rtracklayer", quietly = TRUE) &&
      requireNamespace("GenomicRanges", quietly = TRUE) &&
      !is.null(chain_file) && !is.na(chain_file) && file.exists(chain_file)) {
    return(liftover_rtracklayer(std, chain_file))
  }

  stop(sprintf(paste0(
    "Cross-build liftover %s -> %s requested but no liftover backend is available.\n",
    "Install one of the following and re-run `tumoronly prepare-cosmic`:\n",
    "  - R package 'rtracklayer' + a UCSC chain file (cosmic.chain_file)\n",
    "  - CrossMap (command-line)\n",
    "  - UCSC liftOver (command-line)\n",
    "Run `tumoronly doctor` to see what is detected."),
    source_build, target_build), call. = FALSE)
}

#' Liftover via rtracklayer + a UCSC chain file, tracking unmapped and
#' multi-mapped records as provenance.
#' @keywords internal
liftover_rtracklayer <- function(std, chain_file) {
  chain <- rtracklayer::import.chain(chain_file)
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", std$chrom),
    ranges = IRanges::IRanges(start = std$pos, width = nchar(std$ref)))
  gr$row <- seq_len(nrow(std))
  lifted <- rtracklayer::liftOver(gr, chain)
  n_hits <- lengths(lifted)
  unmapped <- std[n_hits == 0, , drop = FALSE]
  multi    <- std[n_hits > 1, , drop = FALSE]
  ok_idx   <- which(n_hits == 1)
  flat <- unlist(lifted[ok_idx])
  converted <- std[ok_idx, , drop = FALSE]
  converted$chrom <- normalize_chrom(as.character(GenomicRanges::seqnames(flat)))
  converted$pos   <- GenomicRanges::start(flat)
  list(converted = converted, unmapped = unmapped, multi_mapped = multi,
       backend = "rtracklayer",
       tool_versions = list(rtracklayer = as.character(utils::packageVersion("rtracklayer"))))
}

#' Validate REF alleles against the target FASTA (requires Biostrings + indexed
#' FASTA). Returns a logical vector of matches, or NULL if validation could not
#' be performed (recorded in the manifest as ref_validated = FALSE). Never drops
#' silently: the caller writes mismatches to the ref_mismatch provenance table.
#' @keywords internal
validate_ref_against_fasta <- function(df, fasta, target_build) {
  if (is.null(fasta) || is.na(fasta) || !file.exists(fasta)) {
    warning("No FASTA available for REF validation; skipping (manifest ref_validated=FALSE). ",
            "Provide reference.fasta to enable REF-allele checking.", call. = FALSE)
    return(NULL)
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    warning("Biostrings not installed; skipping REF-allele validation. ",
            "Install Biostrings to enable it.", call. = FALSE)
    return(NULL)
  }
  fa <- Biostrings::readDNAStringSet(fasta)
  names(fa) <- normalize_chrom(sub("\\s.*$", "", names(fa)))
  ok <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(df))) {
    chr <- df$chrom[i]
    if (!chr %in% names(fa)) next
    seqi <- as.character(Biostrings::subseq(fa[[chr]], start = df$pos[i], width = nchar(df$ref[i])))
    ok[i] <- toupper(seqi) == toupper(df$ref[i])
  }
  ok
}
