# Annotate variants against the processed COSMIC database using ONLY the
# canonical key BUILD|CHROM|POS|REF|ALT, then stratify the evidence by the
# sample's tumor type using the normalized long context table
# (see R/cosmic_context.R). The database + long table are loaded once.

#' Annotate a variant table with tumor-type-stratified COSMIC evidence.
#'
#' @param x variant data.frame (must contain chrom/pos/ref/alt; uses tumor_type
#'   and, if present, tumor_subtype).
#' @param cfg resolved config.
#' @param build the resolved genome build of the variants.
#' @return `x` with the full COSMIC_* evidence + context columns and the internal
#'   `cosmic_match`/`cosmic_count` columns. `validation_support_score` is set to
#'   the tumor-context COSMIC_CONTEXT_SUPPORT_SCORE (NA -> 0); it feeds only the
#'   EXPERIMENTAL confidence layer, never filter_status (see R/filtering.R).
#' @keywords internal
annotate_cosmic <- function(x, cfg, build) {
  x <- init_cosmic_columns(x)
  loaded <- load_cosmic_db(cfg, build)
  if (is.null(loaded)) {
    log_step("cosmic", "no processed COSMIC database configured; skipping COSMIC annotation")
    x$validation_support_score <- 0
    return(x)
  }
  db <- loaded$db; long <- loaded$long

  release <- unique(db$COSMIC_RELEASE)[1]
  db_build <- unique(db$COSMIC_TARGET_BUILD)[1]
  if (!is.na(db_build) && !identical(db_build, build)) {
    stop(sprintf(paste0(
      "COSMIC database build (%s) does not match the analysis build (%s). ",
      "Re-run prepare-cosmic targeting %s, or set input.genome_build accordingly."),
      db_build, build, build), call. = FALSE)
  }

  mapping <- load_tumor_type_mapping(cfg)
  pan_threshold <- cfg_get(cfg, c("cosmic", "pan_cancer_min_tumor_types"), 5)
  os4_min <- cfg_get(cfg, c("cosmic", "os4_min_matched_occurrences"), 10)  # heuristic
  om4_min <- cfg_get(cfg, c("cosmic", "om4_min_matched_occurrences"), 3)   # heuristic

  keys <- canonical_variant_key(build, x$chrom, x$pos, x$ref, x$alt)
  idx <- match(keys, db$canonical_key)
  hit <- !is.na(idx)

  x$cosmic_match <- hit
  x$COSMIC_MATCH <- hit
  x$COSMIC_RELEASE <- ifelse(hit, release, NA_character_)
  x$COSMIC_MUTATION_IDS <- db$COSMIC_MUTATION_IDS[idx]
  x$cosmic_id <- db$COSMIC_MUTATION_IDS[idx]
  x$COSMIC_TOTAL_OCCURRENCES <- db$COSMIC_OCCURRENCE_COUNT[idx]
  x$COSMIC_OCCURRENCE_COUNT <- db$COSMIC_OCCURRENCE_COUNT[idx]
  x$cosmic_count <- db$COSMIC_OCCURRENCE_COUNT[idx]        # global count (guideline uses this)
  x$COSMIC_TUMOR_TYPES <- db$COSMIC_TUMOR_TYPES[idx]
  x$COSMIC_TUMOR_BREAKDOWN <- db$COSMIC_TUMOR_BREAKDOWN[idx]

  subtype <- as.character(coalesce_columns(
    x, c("tumor_subtype", "Tumor_Subtype", "ONCOTREE_SUBTYPE"), default = NA_character_))

  matched <- which(hit)
  if (length(matched) > 0) {
    tt <- as.character(x$tumor_type); st <- subtype
    # long rows split by key for fast lookup
    long_by_key <- split(seq_len(nrow(long)), long$canonical_key)
    cache <- new.env(parent = emptyenv())
    ctx_list <- lapply(matched, function(i) {
      ck <- paste(db$canonical_key[idx[i]], tt[i], st[i], sep = "\r")
      if (!is.null(cache[[ck]])) return(cache[[ck]])
      rows <- long[long_by_key[[db$canonical_key[idx[i]]]], , drop = FALSE]
      lr <- data.frame(site = rows$cosmic_primary_site, hist = rows$cosmic_histology,
                       subtype = rows$cosmic_subtype, count = rows$occurrence_count,
                       stringsAsFactors = FALSE)
      allowed <- harmonize_sample_tumor(tt[i], st[i], mapping)
      ctx <- classify_cosmic_tumor_context(lr, allowed, sample_known = nrow(allowed) > 0,
        sample_subtype_available = !is.na(st[i]) && nzchar(st[i]), pan_threshold = pan_threshold,
        os4_min = os4_min, om4_min = om4_min)
      cache[[ck]] <- ctx
      ctx
    })
    ctx_df <- do.call(rbind, ctx_list)
    for (col in names(ctx_df)) x[[col]][matched] <- ctx_df[[col]]
  }

  # Regression aid: reproduce the LEGACY behaviour where OS4/OM4 came from the
  # GLOBAL total (default is the contextual behaviour). Documented, off by default.
  if (identical(cfg_get(cfg, c("cosmic", "os4om4_source"), "context"), "global")) {
    tot <- x$COSMIC_TOTAL_OCCURRENCES
    x$OS4_CONTEXTUAL <- !is.na(tot) & tot >= os4_min
    x$OM4_CONTEXTUAL <- !is.na(tot) & tot >= om4_min & tot < os4_min
  }

  x$COSMIC_EVIDENCE_SUMMARY <- ifelse(hit,
    sprintf("COSMIC %s: %s total (%s in matched tumor type); status=%s, level=%s",
            release, ifelse(is.na(x$COSMIC_TOTAL_OCCURRENCES), "?", format(x$COSMIC_TOTAL_OCCURRENCES, trim = TRUE)),
            ifelse(is.na(x$COSMIC_MATCHING_TUMOR_OCCURRENCES), "?", format(x$COSMIC_MATCHING_TUMOR_OCCURRENCES, trim = TRUE)),
            x$COSMIC_TUMOR_CONTEXT_STATUS, x$COSMIC_MATCH_LEVEL),
    NA_character_)

  # tumor-context support feeds the EXPERIMENTAL confidence only (never filters)
  x$validation_support_score <- ifelse(is.na(x$COSMIC_CONTEXT_SUPPORT_SCORE), 0, x$COSMIC_CONTEXT_SUPPORT_SCORE)
  log_step("cosmic", "tumor-type-stratified COSMIC annotation complete",
           matched = sum(hit),
           exact = sum(x$COSMIC_TUMOR_CONTEXT_STATUS == "exact_match", na.rm = TRUE),
           other_only = sum(x$COSMIC_TUMOR_CONTEXT_STATUS == "other_tumor_only", na.rm = TRUE),
           pan_with = sum(x$COSMIC_TUMOR_CONTEXT_STATUS == "pan_cancer_with_sample_tumor", na.rm = TRUE),
           not_evaluable = sum(!isTRUE_vec(x$COSMIC_CONTEXT_EVALUABLE) & hit), total = nrow(x))
  x
}

#' @keywords internal
init_cosmic_columns <- function(x) {
  x$cosmic_match <- FALSE
  x$cosmic_id <- NA_character_
  x$cosmic_count <- NA_real_
  x$COSMIC_MATCH <- FALSE
  x$COSMIC_RELEASE <- NA_character_
  x$COSMIC_MUTATION_IDS <- NA_character_
  x$COSMIC_TOTAL_OCCURRENCES <- NA_real_
  x$COSMIC_OCCURRENCE_COUNT <- NA_real_
  x$COSMIC_TUMOR_TYPES <- NA_character_
  x$COSMIC_TUMOR_BREAKDOWN <- NA_character_
  x$COSMIC_EVIDENCE_SUMMARY <- NA_character_
  x$COSMIC_TUMOR_CONTEXT_STATUS <- "no_cosmic_match"
  x$COSMIC_MATCH_LEVEL <- "none"
  x$COSMIC_SITE_MATCH <- FALSE
  x$COSMIC_HISTOLOGY_MATCH <- FALSE
  x$COSMIC_SUBTYPE_MATCH <- FALSE
  x$COSMIC_MATCHING_TUMOR_OCCURRENCES <- 0
  x$COSMIC_OTHER_TUMOR_OCCURRENCES <- 0
  x$COSMIC_MATCHING_TUMOR_TYPES <- NA_character_
  x$COSMIC_OTHER_TUMOR_TYPES <- NA_character_
  x$COSMIC_CONTEXT_EVALUABLE <- FALSE
  x$COSMIC_GLOBAL_RECURRENCE_SCORE <- 0
  x$COSMIC_MATCHED_TUMOR_RECURRENCE_SCORE <- 0
  x$COSMIC_TUMOR_SPECIFICITY_SCORE <- NA_real_
  x$COSMIC_CONTEXT_SUPPORT_SCORE <- NA_real_
  x$COSMIC_GLOBAL_PANCANCER_RECURRENT <- FALSE
  x$COSMIC_PANCANCER_TUMOR_COUNT <- 0
  x$COSMIC_PANCANCER_INCLUDES_SAMPLE_TUMOR <- FALSE
  x$OS4_CONTEXTUAL <- FALSE
  x$OM4_CONTEXTUAL <- FALSE
  x$COSMIC_GLOBAL_RECURRENT <- FALSE
  x$COSMIC_GLOBAL_EVIDENCE_ONLY <- FALSE
  x$COSMIC_CONTEXT_INTERPRETATION <- "no_cosmic_match"
  x$COSMIC_TUMOR_CONTEXT_REASON <- "no valid COSMIC genomic match"
  x
}

#' Locate the processed COSMIC directory and load the db + long table (once).
#' @return list(db, long) or NULL if not found.
#' @keywords internal
load_cosmic_db <- function(cfg, build) {
  pd <- cfg_get(cfg, c("cosmic", "processed_db"), NULL)
  candidates <- character()
  if (!is.null(pd) && !is.na(pd)) candidates <- c(pd, file.path(pd, "cosmic_db.rds"))
  release <- cfg_get(cfg, c("cosmic", "release"), NULL)
  source_build <- cfg_get(cfg, c("cosmic", "source_build"), build)
  cache_dir <- cfg_get(cfg, c("cosmic", "cache_dir"), "db/cosmic")
  if (!is.null(release) && !is.na(release)) {
    candidates <- c(candidates, file.path(cache_dir, release,
      sprintf("%s_to_%s", source_build, build), "cosmic_db.rds"))
  }
  hit <- candidates[grepl("[.]rds$", candidates) & file.exists(candidates)]
  if (length(hit) == 0) return(NULL)
  db_path <- hit[[1]]
  long_path <- file.path(dirname(db_path), "cosmic_context_long.rds")
  db <- as.data.frame(readRDS(db_path))
  long <- if (file.exists(long_path)) as.data.frame(readRDS(long_path))
          else data.frame(canonical_key = character(), cosmic_primary_site = character(),
                          cosmic_histology = character(), cosmic_subtype = character(),
                          occurrence_count = numeric(), stringsAsFactors = FALSE)
  list(db = db, long = long)
}

#' @keywords internal
isTRUE_vec <- function(v) !is.na(v) & v
