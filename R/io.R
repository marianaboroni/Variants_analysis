read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install with install.packages('yaml').")
  }
  yaml::read_yaml(path)
}

read_variants <- function(path, delimiter = "\t", select = NULL) {
  if (!file.exists(path)) {
    stop("Input variant file not found: ", path)
  }
  if (grepl("[.]vcf([.]gz)?$", path, ignore.case = TRUE)) {
    return(read_annotated_vcf(path))
  }
  if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    return(tryCatch(
      data.table::fread(path, sep = delimiter, select = select, data.table = FALSE, na.strings = c("", ".", "NA")),
      error = function(e) {
        if (!grepl("memory map", conditionMessage(e), ignore.case = TRUE)) stop(e)
        # fread's native gz decompression mmaps the decompressed buffer, which can
        # fail with "could not memory map it" on clusters with restrictive virtual
        # address space limits or fragmented address space late in a long run.
        # Piping through zcat avoids that mmap path entirely.
        data.table::fread(cmd = paste("zcat", shQuote(path)), sep = delimiter, select = select,
                           data.table = FALSE, na.strings = c("", ".", "NA"))
      }))
  }
  data.table::fread(path, sep = delimiter, select = select, data.table = FALSE, na.strings = c("", ".", "NA"))
}

read_annotated_vcf <- function(path) {
  meta <- read_vcf_meta(path)
  csq_fields <- parse_vcf_annotation_format(meta, "CSQ")
  ann_fields <- parse_vcf_annotation_format(meta, "ANN")
  x <- data.table::fread(path, skip = "#CHROM", data.table = FALSE, na.strings = c("", ".", "NA"))
  names(x)[names(x) == "#CHROM"] <- "CHROM"
  if (!"INFO" %in% names(x)) {
    stop("VCF input does not contain an INFO column: ", path)
  }
  x <- expand_vcf_samples(x, path)
  x <- add_info_annotations(x, csq_fields, ann_fields)
  x
}

read_vcf_meta <- function(path, n = 10000) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  out <- character()
  repeat {
    z <- readLines(con, n = 1, warn = FALSE)
    if (length(z) == 0) break
    out <- c(out, z)
    if (startsWith(z, "#CHROM") || length(out) >= n) break
  }
  out
}

parse_vcf_annotation_format <- function(meta, id) {
  line <- meta[grepl(paste0("ID=", id, ","), meta, fixed = TRUE)]
  if (length(line) == 0) return(character())
  line <- line[[1]]
  fmt <- sub(".*Format: ", "", line)
  fmt <- sub("[\">].*", "", fmt)
  fields <- strsplit(fmt, "[|]", perl = TRUE)[[1]]
  trimws(fields)
}

expand_vcf_samples <- function(x, path) {
  if (!"FORMAT" %in% names(x)) {
    x$Tumor_Sample_Barcode <- tools::file_path_sans_ext(basename(path))
    return(x)
  }
  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  sample_cols <- setdiff(names(x), fixed_cols)
  if (length(sample_cols) == 0) {
    x$Tumor_Sample_Barcode <- tools::file_path_sans_ext(basename(path))
    return(x)
  }

  pieces <- lapply(sample_cols, function(sample_col) {
    y <- x[, fixed_cols[fixed_cols %in% names(x)], drop = FALSE]
    y$Tumor_Sample_Barcode <- sample_col
    fmt_str <- as.character(y$FORMAT)
    val_str <- as.character(x[[sample_col]])
    y$GT <- extract_format_field(fmt_str, val_str, "GT")
    y$AD <- extract_format_field(fmt_str, val_str, "AD")
    y$DP <- to_numeric_safe(extract_format_field(fmt_str, val_str, "DP"))
    y$AF <- to_numeric_safe(extract_format_field(fmt_str, val_str, "AF"))
    y$TLOD <- to_numeric_safe(extract_format_field(fmt_str, val_str, "TLOD"))
    y
  })
  do.call(rbind, pieces)
}

#' Extract one FORMAT field (e.g. "DP") from colon-delimited FORMAT/genotype
#' strings. Groups rows by their (usually 1-2 distinct, for a whole VCF)
#' FORMAT string and extracts each group with one vectorized
#' data.table::tstrsplit() call, instead of a per-row match()+index - needed
#' at WGS scale (millions of rows), where a per-row R loop turns minutes of
#' work into hours.
#' @param fmt_str character vector, the FORMAT column (e.g. "GT:AD:DP:AF").
#' @param val_str character vector, the matching genotype column.
#' @keywords internal
extract_format_field <- function(fmt_str, val_str, key) {
  fmt_str <- as.character(fmt_str); val_str <- as.character(val_str)
  out <- rep(NA_character_, length(fmt_str))
  for (pattern in unique(fmt_str)) {
    if (is.na(pattern)) next
    idx <- match(key, strsplit(pattern, ":", fixed = TRUE)[[1]])
    if (is.na(idx)) next
    sel <- which(fmt_str == pattern)
    parts <- data.table::tstrsplit(val_str[sel], ":", fixed = TRUE)
    if (idx <= length(parts)) out[sel] <- parts[[idx]]
  }
  out
}

add_info_annotations <- function(x, csq_fields, ann_fields) {
  info <- as.character(x$INFO)
  x$CSQ <- extract_info_value(info, "CSQ")
  x$ANN <- extract_info_value(info, "ANN")

  if (any(!is.na(x$CSQ)) && length(csq_fields) > 0) {
    csq <- parse_pipe_annotation(x$CSQ, csq_fields)
    x$SYMBOL <- first_nonmissing_annotation(csq, c("SYMBOL", "Gene", "Feature"))
    x$Consequence <- first_nonmissing_annotation(csq, c("Consequence"))
    x$HGVSp <- first_nonmissing_annotation(csq, c("HGVSp", "Protein_position", "Amino_acids"))
    x$Existing_variation <- first_nonmissing_annotation(csq, c("Existing_variation"))
    x$CLINVAR_SIG <- first_nonmissing_annotation(csq, c("CLIN_SIG", "ClinVar_CLNSIG"))
    x$CLIN_SIG <- first_nonmissing_annotation(csq, c("CLIN_SIG"))
    x$gnomADe_AF <- to_numeric_safe(first_nonmissing_annotation(csq, c("gnomADe_AF", "gnomADg_AF", "gnomAD_AF", "MAX_AF")))
    # Surface ALL predictor CSQ fields present (dynamic, from the registry) plus
    # transcript/population fields, from the first/most-severe VEP annotation, so
    # EVERY predictor that exists is available to the predictor layer + MAF export.
    pred_aliases <- tryCatch(predictor_alias_union(), error = function(e) character())
    surface <- unique(c(pred_aliases, "IMPACT", "HGVSc", "EXON", "INTRON", "Feature",
                 "CANONICAL", "MANE_SELECT", "BIOTYPE", "Amino_acids", "gnomADg_AF", "MAX_AF"))
    for (f in intersect(surface, csq_fields)) if (!(f %in% names(x))) x[[f]] <- csq[[f]]
  }

  if (any(!is.na(x$ANN)) && length(ann_fields) > 0) {
    ann <- parse_pipe_annotation(x$ANN, ann_fields)
    x$SYMBOL <- coalesce_value(x$SYMBOL, first_nonmissing_annotation(ann, c("Gene_Name", "Gene_ID")))
    x$Consequence <- coalesce_value(x$Consequence, first_nonmissing_annotation(ann, c("Annotation")))
    x$HGVSp <- coalesce_value(x$HGVSp, first_nonmissing_annotation(ann, c("HGVS.p", "HGVS_p")))
  }

  x$Gene.refGene <- extract_info_value(info, "Gene.refGene")
  x$ExonicFunc.refGene <- extract_info_value(info, "ExonicFunc.refGene")
  x$Func.refGene <- extract_info_value(info, "Func.refGene")
  x$AAChange.refGene <- extract_info_value(info, "AAChange.refGene")
  x$SYMBOL <- coalesce_value(x$SYMBOL, x$Gene.refGene)
  x$Consequence <- coalesce_value(x$Consequence, x$ExonicFunc.refGene)
  x$Consequence <- coalesce_value(x$Consequence, x$Func.refGene)
  x$HGVSp <- coalesce_value(x$HGVSp, x$AAChange.refGene)

  pop_keys <- c("AF", "AF_popmax", "gnomAD_AF", "gnomADg_AF", "gnomADe_AF", "GMAF")
  for (key in pop_keys) {
    if (!(key %in% names(x))) x[[key]] <- to_numeric_safe(extract_info_value(info, key))
  }
  x
}

extract_info_value <- function(info, key) {
  prefix <- paste0(key, "=")
  vapply(strsplit(info, ";", fixed = TRUE), function(fields) {
    hit <- fields[startsWith(fields, prefix)]
    if (length(hit) == 0) return(NA_character_)
    sub(prefix, "", hit[[1]], fixed = TRUE)
  }, character(1))
}

#' Parse the FIRST pipe-delimited transcript annotation (VEP CSQ / snpEff ANN)
#' out of every row at once with data.table::tstrsplit(), instead of a
#' per-row x per-field double loop - the double loop is O(n_rows x n_fields)
#' and is the single largest ingestion bottleneck at WGS scale (millions of
#' rows x 20+ CSQ fields).
#' @keywords internal
parse_pipe_annotation <- function(values, fields) {
  first <- sub(",.*$", "", as.character(values))
  empty <- is.na(first) | first == ""
  parts <- data.table::tstrsplit(first, "|", fixed = TRUE, fill = NA_character_)
  n <- min(length(parts), length(fields))
  out <- setNames(vector("list", length(fields)), fields)
  for (j in seq_along(fields)) {
    out[[j]] <- if (j <= n) parts[[j]] else rep(NA_character_, length(values))
    if (any(empty)) out[[j]][empty] <- NA_character_
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

first_nonmissing_annotation <- function(x, candidates) {
  candidates <- candidates[candidates %in% names(x)]
  if (length(candidates) == 0) return(rep(NA_character_, nrow(x)))
  out <- x[[candidates[[1]]]]
  for (nm in candidates[-1]) out <- coalesce_value(out, x[[nm]])
  out
}

coalesce_value <- function(a, b) {
  if (is.null(a)) return(b)
  out <- a
  replace <- is.na(out) | out == "" | out == "."
  out[replace] <- b[replace]
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

coalesce_columns <- function(x, candidates, default = NA) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) return(rep(default, nrow(x)))
  out <- x[[hit[[1]]]]
  if (length(hit) > 1) {
    for (nm in hit[-1]) {
      out[is.na(out) | out == ""] <- x[[nm]][is.na(out) | out == ""]
    }
  }
  out
}

to_numeric_safe <- function(x) {
  if (is.null(x)) return(NA_real_)
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = FALSE)))
}

is_missing_value <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

cfg_get <- function(cfg, path, default = NULL) {
  if (is.null(cfg)) return(default)
  cur <- cfg
  for (key in path) {
    if (is.null(cur[[key]])) return(default)
    cur <- cur[[key]]
  }
  if (length(cur) == 1 && is.na(cur)) default else cur
}

default_tumor_type <- function(cfg = NULL) {
  value <- cfg_get(cfg, c("cancer", "default_tumor_type"), NULL)
  if (is.null(value)) value <- cfg_get(cfg, c("cohort", "cancer_type"), "PANCANCER")
  as.character(value)
}

standardize_variant_table <- function(x, cfg = NULL) {
  x$sample_id <- as.character(coalesce_columns(
    x,
    c("Tumor_Sample_Barcode", "Sample_Barcode", "Tumor_Sample", "Sample", "sample", "sample_id")
  ))
  x$chrom <- as.character(coalesce_columns(x, c("CHROM", "Chromosome", "chr", "chrom")))
  x$pos <- to_numeric_safe(coalesce_columns(x, c("START", "Start_Position", "POS", "pos")))
  x$ref <- as.character(coalesce_columns(x, c("REF", "Reference_Allele", "ref")))
  x$alt <- as.character(coalesce_columns(x, c("ALT", "Tumor_Seq_Allele2", "alt")))
  x$gene <- as.character(coalesce_columns(x, c("Hugo_Symbol", "SYMBOL", "Gene", "gene")))
  x$consequence <- as.character(coalesce_columns(
    x,
    c("Consequence", "Variant_Classification", "EFFECT", "Annotation")
  ))
  x$protein_change <- as.character(coalesce_columns(
    x,
    c("HGVSp_Short", "HGVSp", "Protein_Change", "Amino_acids")
  ))
  x$tumor_type <- as.character(coalesce_columns(
    x,
    c("ONCOTREE_CODE", "Oncotree_Code", "Cancer_Type", "Tumor_Type", "tumor_type", "Primary_Site"),
    default = default_tumor_type(cfg)
  ))
  x$filter_status <- as.character(coalesce_columns(x, c("FILTER", "filter", "Filter"), default = "NA"))

  if (any(is.na(x$sample_id))) {
    warning("Some variants have missing sample IDs.")
  }
  required <- c("chrom", "pos", "ref", "alt")
  missing_core <- vapply(required, function(nm) all(is.na(x[[nm]]) | x[[nm]] == ""), logical(1))
  if (any(missing_core)) {
    stop("Missing core coordinate columns after standardization: ",
         paste(required[missing_core], collapse = ", "))
  }
  x$variant_id <- paste(x$chrom, x$pos, x$ref, x$alt, sep = ":")
  x$variant_key <- x$variant_id
  x$sample_variant_key <- ifelse(
    !is.na(x$sample_id) & x$sample_id != "",
    paste(x$sample_id, x$variant_id, sep = "|")
    , NA_character_
  )
  x$locus_id <- paste(x$chrom, x$pos, sep = ":")
  x$tumor_type[is_missing_value(x$tumor_type)] <- default_tumor_type(cfg)
  x$tumor_type <- toupper(trimws(x$tumor_type))
  x
}
