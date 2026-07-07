run_ancestry_analysis <- function(x, cfg) {
  enabled <- isTRUE(cfg_get(cfg, c("ancestry", "enabled"), FALSE))
  if (!enabled) {
    return(empty_ancestry_result("disabled"))
  }

  outdir <- cfg$output$dir
  figdir <- file.path(outdir, "figures")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  snps <- select_ancestry_snps(x, cfg)
  min_snps <- cfg_get(cfg, c("ancestry", "min_snps"), 200)
  if (nrow(snps) < min_snps) {
    result <- empty_ancestry_result("insufficient_snps")
    result$snp_set <- snps
    result$sample_qc <- ancestry_sample_qc_from_snps(snps)
    result$snp_qc <- ancestry_snp_qc_from_snps(snps)
    return(result)
  }

  matrix_result <- build_ancestry_genotype_matrix(snps, cfg)
  sample_qc <- matrix_result$sample_qc
  snp_qc <- matrix_result$snp_qc
  genotype_matrix <- matrix_result$genotype_matrix
  snp_info <- matrix_result$snp_info

  if (nrow(genotype_matrix) == 0 || ncol(genotype_matrix) < min_snps) {
    result <- empty_ancestry_result("insufficient_snps_after_qc")
    result$snp_set <- snps
    result$sample_qc <- sample_qc
    result$snp_qc <- snp_qc
    return(result)
  }

  genotype_matrix <- ancestry_ld_prune(genotype_matrix, snp_info, cfg)
  snp_info <- snp_info[snp_info$snp_id %in% colnames(genotype_matrix), , drop = FALSE]

  reference <- load_ancestry_reference(cfg)
  if (nrow(reference) > 0) {
    pca_result <- ancestry_reference_projection(genotype_matrix, snp_info, reference, cfg)
  } else {
    pca_result <- ancestry_cohort_pca(genotype_matrix, snp_info, cfg)
  }

  figure_manifest <- plot_ancestry_outputs(pca_result, sample_qc, cfg, figdir)
  status <- ancestry_status_table(
    status = pca_result$status,
    method = pca_result$method,
    n_samples = nrow(genotype_matrix),
    n_snps = ncol(genotype_matrix),
    n_reference_samples = nrow(unique_reference_samples(reference)),
    n_reference_snps = length(unique(reference$snp_id))
  )

  list(
    status = status,
    snp_set = snps,
    sample_qc = sample_qc,
    snp_qc = snp_qc,
    pca_scores = pca_result$pca_scores,
    reference_scores = pca_result$reference_scores,
    reference_centroids = pca_result$reference_centroids,
    assignments = pca_result$assignments,
    figure_manifest = figure_manifest
  )
}

empty_ancestry_result <- function(status) {
  list(
    status = ancestry_status_table(status = status),
    snp_set = data.frame(),
    sample_qc = data.frame(),
    snp_qc = data.frame(),
    pca_scores = empty_ancestry_pca_scores(),
    reference_scores = empty_ancestry_reference_scores(),
    reference_centroids = empty_ancestry_reference_centroids(),
    assignments = empty_ancestry_assignments(),
    figure_manifest = empty_figure_manifest()
  )
}

ancestry_status_table <- function(
    status,
    method = NA_character_,
    n_samples = NA_integer_,
    n_snps = NA_integer_,
    n_reference_samples = NA_integer_,
    n_reference_snps = NA_integer_) {
  data.frame(
    metric = c("status", "method", "n_samples", "n_snps", "n_reference_samples", "n_reference_snps"),
    value = as.character(c(status, method, n_samples, n_snps, n_reference_samples, n_reference_snps)),
    stringsAsFactors = FALSE
  )
}

select_ancestry_snps <- function(x, cfg) {
  autosomal <- is_autosomal_chrom(x$chrom)
  snv <- is_biallelic_acgt_snv(x$ref, x$alt)
  min_depth <- cfg_get(cfg, c("ancestry", "min_depth"), 10)
  min_alt_count <- cfg_get(cfg, c("ancestry", "min_alt_count"), 3)
  min_pop_af <- cfg_get(cfg, c("ancestry", "min_population_af"), 0.01)
  max_pop_af <- cfg_get(cfg, c("ancestry", "max_population_af"), 0.99)
  het_min <- cfg_get(cfg, c("ancestry", "het_vaf_min"), 0.25)
  hom_alt_min <- cfg_get(cfg, c("ancestry", "hom_alt_vaf_min"), 0.80)
  require_pop_af <- isTRUE(cfg_get(cfg, c("ancestry", "require_population_af"), TRUE))
  include_classes <- cfg_get(
    cfg,
    c("ancestry", "include_final_classes"),
    c(
      "likely_germline", "manual_review_required", "uncertain_tumor_only",
      "probable_germline", "uncertain"
    )
  )

  technical_ok <- (is.na(x$dp) | x$dp >= min_depth) &
    (is.na(x$alt_count) | x$alt_count >= min_alt_count) &
    !is.na(x$vaf) & x$vaf >= het_min
  pop_ok <- !is.na(x$max_pop_af) & x$max_pop_af >= min_pop_af & x$max_pop_af <= max_pop_af
  if (!require_pop_af) {
    pop_ok <- pop_ok | is.na(x$max_pop_af)
  }
  class_ok <- x$final_class %in% include_classes
  not_driver_hotspot <- !(safe_character_column(x, "driver_class") %in% c("known_driver", "probable_driver"))

  keep <- autosomal & snv & technical_ok & pop_ok & class_ok & not_driver_hotspot
  if (!any(keep)) {
    return(empty_ancestry_snp_set())
  }

  out <- x[keep, , drop = FALSE]
  out$snp_id <- make_snp_id(out)
  out$ancestry_dosage <- ancestry_vaf_to_dosage(out$vaf, hom_alt_min)
  out$ancestry_genotype_class <- ifelse(out$ancestry_dosage == 2, "hom_alt", "het")

  cols <- c(
    "sample_id", "tumor_type", "chrom", "pos", "ref", "alt", "snp_id",
    "gene", "final_class", "dp", "alt_count", "ref_count", "vaf",
    "max_pop_af", "ancestry_dosage", "ancestry_genotype_class",
    "variant_cohort_freq", "variant_tumor_type_freq"
  )
  cols <- cols[cols %in% names(out)]
  out <- out[, cols, drop = FALSE]
  out <- cap_ancestry_candidate_snps(out, cfg)
  out
}

empty_ancestry_snp_set <- function() {
  data.frame(
    sample_id = character(),
    tumor_type = character(),
    chrom = character(),
    pos = numeric(),
    ref = character(),
    alt = character(),
    snp_id = character(),
    gene = character(),
    final_class = character(),
    dp = numeric(),
    alt_count = numeric(),
    ref_count = numeric(),
    vaf = numeric(),
    max_pop_af = numeric(),
    ancestry_dosage = numeric(),
    ancestry_genotype_class = character(),
    stringsAsFactors = FALSE
  )
}

is_autosomal_chrom <- function(chrom) {
  z <- gsub("^chr", "", as.character(chrom), ignore.case = TRUE)
  suppressWarnings(!is.na(as.integer(z)) & as.integer(z) >= 1 & as.integer(z) <= 22)
}

is_biallelic_acgt_snv <- function(ref, alt) {
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  nchar(ref) == 1 & nchar(alt) == 1 &
    ref %in% c("A", "C", "G", "T") &
    alt %in% c("A", "C", "G", "T") &
    !grepl(",", alt, fixed = TRUE)
}

make_snp_id <- function(x) {
  paste(gsub("^chr", "", as.character(x$chrom), ignore.case = TRUE), x$pos, x$ref, x$alt, sep = ":")
}

ancestry_vaf_to_dosage <- function(vaf, hom_alt_min = 0.80) {
  ifelse(vaf >= hom_alt_min, 2, 1)
}

cap_ancestry_candidate_snps <- function(snps, cfg) {
  max_candidate_snps <- cfg_get(cfg, c("ancestry", "max_candidate_snps"), 200000)
  n_unique <- length(unique(snps$snp_id))
  if (n_unique <= max_candidate_snps) return(snps)

  dt <- data.table::as.data.table(snps)
  snp_rank <- dt[, .(
    n_samples = data.table::uniqueN(sample_id),
    median_depth = safe_median(dp),
    maf_proxy = safe_median(pmin(ancestry_dosage, 2 - ancestry_dosage) / 2),
    chrom_num = suppressWarnings(as.integer(gsub("^chr", "", as.character(chrom), ignore.case = TRUE))),
    pos_num = min(pos, na.rm = TRUE)
  ), by = snp_id]
  snp_rank[, maf_balance := abs(maf_proxy - 0.25)]
  snp_rank <- snp_rank[order(-n_samples, maf_balance, -median_depth, chrom_num, pos_num)]
  keep_snps <- snp_rank$snp_id[seq_len(max_candidate_snps)]
  as.data.frame(dt[snp_id %in% keep_snps])
}

build_ancestry_genotype_matrix <- function(snps, cfg) {
  if (nrow(snps) == 0) {
    return(list(
      genotype_matrix = matrix(nrow = 0, ncol = 0),
      snp_info = data.frame(),
      sample_qc = data.frame(),
      snp_qc = data.frame()
    ))
  }

  dt <- data.table::as.data.table(snps)
  dt[, depth_for_choice := ifelse(is.na(dp), 0, dp)]
  data.table::setorder(dt, sample_id, snp_id, -depth_for_choice)
  dt <- dt[, .SD[1], by = .(sample_id, snp_id)]

  wide <- data.table::dcast(dt, sample_id ~ snp_id, value.var = "ancestry_dosage")
  samples <- wide$sample_id
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- samples
  storage.mode(mat) <- "numeric"

  snp_info <- unique(dt[, .(
    snp_id,
    chrom = gsub("^chr", "", as.character(chrom), ignore.case = TRUE),
    pos = as.numeric(pos),
    ref = as.character(ref),
    alt = as.character(alt)
  )])
  snp_info <- as.data.frame(snp_info[match(colnames(mat), snp_id)])

  sample_qc <- ancestry_sample_qc_from_matrix(mat, snps, cfg)
  snp_qc <- ancestry_snp_qc_from_matrix(mat, snp_info)

  min_snp_call_rate <- cfg_get(cfg, c("ancestry", "min_snp_call_rate"), 0.20)
  min_sample_snps <- cfg_get(cfg, c("ancestry", "min_snps_per_sample"), 100)
  min_maf <- cfg_get(cfg, c("ancestry", "min_maf"), 0.01)
  max_missing <- cfg_get(cfg, c("ancestry", "max_sample_missing_rate"), 0.95)

  keep_snp <- snp_qc$call_rate >= min_snp_call_rate & snp_qc$maf >= min_maf & snp_qc$maf <= (1 - min_maf)
  keep_snp[is.na(keep_snp)] <- FALSE
  mat <- mat[, keep_snp, drop = FALSE]
  snp_info <- snp_info[keep_snp, , drop = FALSE]
  snp_qc <- snp_qc[keep_snp, , drop = FALSE]

  sample_called <- rowSums(!is.na(mat))
  sample_missing <- rowMeans(is.na(mat))
  keep_sample <- sample_called >= min_sample_snps & sample_missing <= max_missing
  keep_sample[is.na(keep_sample)] <- FALSE
  mat <- mat[keep_sample, , drop = FALSE]
  sample_qc <- sample_qc[match(rownames(mat), sample_qc$sample_id), , drop = FALSE]

  imputation <- cfg_get(cfg, c("ancestry", "missing_genotype_imputation"), "mean")
  mat <- impute_ancestry_genotypes(mat, method = imputation)

  list(
    genotype_matrix = mat,
    snp_info = snp_info,
    sample_qc = sample_qc,
    snp_qc = snp_qc
  )
}

impute_ancestry_genotypes <- function(mat, method = "mean") {
  if (length(mat) == 0) return(mat)
  if (method == "zero") {
    mat[is.na(mat)] <- 0
    return(mat)
  }
  for (j in seq_len(ncol(mat))) {
    z <- mat[, j]
    fill <- if (all(is.na(z))) 0 else mean(z, na.rm = TRUE)
    z[is.na(z)] <- fill
    mat[, j] <- z
  }
  mat
}

ancestry_sample_qc_from_snps <- function(snps) {
  if (nrow(snps) == 0) return(data.frame())
  dt <- data.table::as.data.table(snps)
  as.data.frame(dt[, .(
    n_candidate_snps = data.table::uniqueN(snp_id),
    median_snp_depth = safe_median(dp),
    median_snp_vaf = safe_median(vaf),
    n_het_like = sum(ancestry_dosage == 1, na.rm = TRUE),
    n_hom_alt_like = sum(ancestry_dosage == 2, na.rm = TRUE)
  ), by = .(sample_id, tumor_type)])
}

ancestry_snp_qc_from_snps <- function(snps) {
  if (nrow(snps) == 0) return(data.frame())
  dt <- data.table::as.data.table(snps)
  as.data.frame(dt[, .(
    n_samples_observed = data.table::uniqueN(sample_id),
    median_depth = safe_median(dp),
    median_vaf = safe_median(vaf),
    max_pop_af = safe_median(max_pop_af)
  ), by = .(snp_id, chrom, pos, ref, alt)])
}

ancestry_sample_qc_from_matrix <- function(mat, snps, cfg) {
  sample_meta <- unique(snps[, c("sample_id", "tumor_type"), drop = FALSE])
  out <- data.frame(
    sample_id = rownames(mat),
    n_ancestry_snps_called = rowSums(!is.na(mat)),
    ancestry_missing_rate = rowMeans(is.na(mat)),
    ancestry_heterozygosity_proxy = rowMeans(mat == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  merge(out, sample_meta, by = "sample_id", all.x = TRUE)
}

ancestry_snp_qc_from_matrix <- function(mat, snp_info) {
  if (ncol(mat) == 0) return(data.frame())
  call_rate <- colMeans(!is.na(mat))
  alt_freq <- colMeans(mat, na.rm = TRUE) / 2
  maf <- pmin(alt_freq, 1 - alt_freq)
  data.frame(
    snp_info,
    call_rate = call_rate,
    alt_freq = alt_freq,
    maf = maf,
    stringsAsFactors = FALSE
  )
}

ancestry_ld_prune <- function(mat, snp_info, cfg) {
  enabled <- isTRUE(cfg_get(cfg, c("ancestry", "ld_prune"), TRUE))
  if (!enabled || ncol(mat) < 2) return(mat)

  max_snps <- cfg_get(cfg, c("ancestry", "max_snps_for_pca"), 50000)
  r2_threshold <- cfg_get(cfg, c("ancestry", "ld_r2_threshold"), 0.20)
  window_bp <- cfg_get(cfg, c("ancestry", "ld_window_bp"), 500000)

  info <- snp_info
  info$chrom_num <- suppressWarnings(as.integer(info$chrom))
  info$pos <- as.numeric(info$pos)
  keep <- rep(FALSE, nrow(info))
  kept_indices <- integer()
  order_idx <- order(info$chrom_num, info$pos)

  for (idx in order_idx) {
    same_window <- kept_indices[
      info$chrom_num[kept_indices] == info$chrom_num[idx] &
        abs(info$pos[kept_indices] - info$pos[idx]) <= window_bp
    ]
    if (length(same_window) == 0) {
      keep[idx] <- TRUE
      kept_indices <- c(kept_indices, idx)
    } else {
      r2 <- vapply(same_window, function(j) {
        r <- suppressWarnings(stats::cor(mat[, idx], mat[, j], use = "pairwise.complete.obs"))
        if (is.na(r)) 0 else r^2
      }, numeric(1))
      if (all(r2 < r2_threshold)) {
        keep[idx] <- TRUE
        kept_indices <- c(kept_indices, idx)
      }
    }
    if (sum(keep) >= max_snps) break
  }

  mat[, keep, drop = FALSE]
}

load_ancestry_reference <- function(cfg) {
  path <- cfg_get(cfg, c("ancestry", "reference_genotypes"), NULL)
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    return(data.frame())
  }
  delimiter <- cfg_get(cfg, c("ancestry", "reference_delimiter"), "\t")
  ref <- read_variants(path, delimiter)
  standardize_ancestry_reference(ref)
}

standardize_ancestry_reference <- function(ref) {
  out <- data.frame(
    sample_id = as.character(coalesce_columns(ref, c("sample_id", "Sample", "IID", "id"))),
    population = as.character(coalesce_columns(ref, c("population", "Population", "pop", "POP"), default = NA)),
    superpopulation = as.character(coalesce_columns(ref, c("superpopulation", "Superpopulation", "continental_group", "superpop", "SUPERPOP"), default = NA)),
    chrom = as.character(coalesce_columns(ref, c("CHROM", "Chromosome", "chr", "chrom"))),
    pos = to_numeric_safe(coalesce_columns(ref, c("POS", "START", "Start_Position", "pos"))),
    ref = as.character(coalesce_columns(ref, c("REF", "Reference_Allele", "ref"))),
    alt = as.character(coalesce_columns(ref, c("ALT", "Tumor_Seq_Allele2", "alt"))),
    dosage = to_numeric_safe(coalesce_columns(ref, c("dosage", "genotype_dosage", "GT_dosage", "alt_dosage"))),
    stringsAsFactors = FALSE
  )
  out$chrom <- gsub("^chr", "", out$chrom, ignore.case = TRUE)
  out$snp_id <- make_snp_id(out)
  out <- out[!is_missing_value(out$sample_id) & !is.na(out$dosage) & out$dosage >= 0 & out$dosage <= 2, ]
  out
}

unique_reference_samples <- function(reference) {
  if (nrow(reference) == 0) return(data.frame())
  unique(reference[, c("sample_id", "population", "superpopulation"), drop = FALSE])
}

ancestry_reference_projection <- function(cohort_mat, snp_info, reference, cfg) {
  ref_samples <- unique_reference_samples(reference)
  common_snps <- intersect(colnames(cohort_mat), reference$snp_id)
  min_overlap <- cfg_get(cfg, c("ancestry", "min_reference_overlap_snps"), 200)
  if (length(common_snps) < min_overlap) {
    fallback <- ancestry_cohort_pca(cohort_mat, snp_info, cfg)
    fallback$status <- "insufficient_reference_overlap_used_cohort_pca"
    return(fallback)
  }

  ref_dt <- data.table::as.data.table(reference[reference$snp_id %in% common_snps, ])
  ref_wide <- data.table::dcast(ref_dt, sample_id ~ snp_id, value.var = "dosage", fun.aggregate = mean)
  ref_mat <- as.matrix(ref_wide[, common_snps, with = FALSE])
  rownames(ref_mat) <- ref_wide$sample_id
  storage.mode(ref_mat) <- "numeric"
  ref_mat <- impute_ancestry_genotypes(ref_mat, "mean")

  cohort_mat <- cohort_mat[, common_snps, drop = FALSE]
  cohort_mat <- impute_ancestry_genotypes(cohort_mat, "mean")

  variable <- matrix_variable_columns(ref_mat)
  ref_mat <- ref_mat[, variable, drop = FALSE]
  cohort_mat <- cohort_mat[, variable, drop = FALSE]
  if (ncol(ref_mat) < min_overlap) {
    fallback <- ancestry_cohort_pca(cohort_mat, snp_info, cfg)
    fallback$status <- "insufficient_variable_reference_snps_used_cohort_pca"
    return(fallback)
  }

  n_pcs <- cfg_get(cfg, c("ancestry", "n_pcs"), 10)
  n_pcs <- min(n_pcs, nrow(ref_mat) - 1, ncol(ref_mat))
  pca <- stats::prcomp(ref_mat, center = TRUE, scale. = TRUE, rank. = n_pcs)

  cohort_scaled <- scale(cohort_mat, center = pca$center, scale = pca$scale)
  cohort_scores <- cohort_scaled %*% pca$rotation[, seq_len(n_pcs), drop = FALSE]
  colnames(cohort_scores) <- paste0("PC", seq_len(ncol(cohort_scores)))

  ref_scores <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  ref_scores$sample_id <- rownames(pca$x)
  ref_scores <- merge(ref_scores, ref_samples, by = "sample_id", all.x = TRUE)
  ref_scores$sample_type <- "reference"

  cohort_scores_df <- as.data.frame(cohort_scores)
  cohort_scores_df$sample_id <- rownames(cohort_mat)
  cohort_scores_df$sample_type <- "cohort"

  centroids <- reference_centroids(ref_scores, cfg)
  assignments <- assign_ancestry_by_centroid(cohort_scores_df, centroids, cfg)

  list(
    status = "trained_reference_projection",
    method = "reference_pca_projection",
    pca_scores = cohort_scores_df,
    reference_scores = ref_scores,
    reference_centroids = centroids,
    assignments = assignments
  )
}

ancestry_cohort_pca <- function(cohort_mat, snp_info, cfg) {
  variable <- matrix_variable_columns(cohort_mat)
  cohort_mat <- cohort_mat[, variable, drop = FALSE]
  n_pcs <- cfg_get(cfg, c("ancestry", "n_pcs"), 10)
  n_pcs <- min(n_pcs, nrow(cohort_mat) - 1, ncol(cohort_mat))
  if (n_pcs < 1) {
    return(list(
      status = "insufficient_samples_for_pca",
      method = "cohort_pca_no_reference",
      pca_scores = empty_ancestry_pca_scores(),
      reference_scores = empty_ancestry_reference_scores(),
      reference_centroids = empty_ancestry_reference_centroids(),
      assignments = empty_ancestry_assignments()
    ))
  }

  pca <- stats::prcomp(cohort_mat, center = TRUE, scale. = TRUE, rank. = n_pcs)
  scores <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  scores$sample_id <- rownames(cohort_mat)
  scores$sample_type <- "cohort"
  scores$ancestry_note <- "cohort_pca_without_reference_not_a_population_assignment"

  list(
    status = "cohort_pca_no_reference",
    method = "cohort_pca_no_reference",
    pca_scores = scores,
    reference_scores = data.frame(),
    reference_centroids = data.frame(),
    assignments = data.frame(
      sample_id = rownames(cohort_mat),
      assignment_status = "no_reference_panel",
      nearest_population = NA_character_,
      nearest_superpopulation = NA_character_,
      ancestry_confidence = NA_real_,
      ancestry_note = "PCA sem painel de referencia mostra estrutura interna, mas nao atribui ancestralidade.",
      stringsAsFactors = FALSE
    )
  )
}

empty_ancestry_pca_scores <- function() {
  data.frame(
    sample_id = character(),
    sample_type = character(),
    PC1 = numeric(),
    PC2 = numeric(),
    ancestry_note = character(),
    stringsAsFactors = FALSE
  )
}

empty_ancestry_reference_scores <- function() {
  data.frame(
    sample_id = character(),
    population = character(),
    superpopulation = character(),
    sample_type = character(),
    PC1 = numeric(),
    PC2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_ancestry_reference_centroids <- function() {
  data.frame(
    population = character(),
    superpopulation = character(),
    PC1 = numeric(),
    PC2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_ancestry_assignments <- function() {
  data.frame(
    sample_id = character(),
    assignment_status = character(),
    nearest_population = character(),
    nearest_superpopulation = character(),
    ancestry_confidence = numeric(),
    nearest_distance = numeric(),
    second_nearest_population = character(),
    second_nearest_superpopulation = character(),
    second_nearest_distance = numeric(),
    ancestry_note = character(),
    stringsAsFactors = FALSE
  )
}

matrix_variable_columns <- function(mat) {
  if (ncol(mat) == 0) return(logical(0))
  vars <- apply(mat, 2, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  vars > 0
}

reference_centroids <- function(ref_scores, cfg) {
  pc_cols <- grep("^PC[0-9]+$", names(ref_scores), value = TRUE)
  n_assign_pcs <- cfg_get(cfg, c("ancestry", "assignment_pcs"), 5)
  pc_cols <- pc_cols[seq_len(min(length(pc_cols), n_assign_pcs))]
  dt <- data.table::as.data.table(ref_scores)
  cent <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(population, superpopulation), .SDcols = pc_cols]
  as.data.frame(cent)
}

assign_ancestry_by_centroid <- function(cohort_scores, centroids, cfg) {
  pc_cols <- intersect(grep("^PC[0-9]+$", names(cohort_scores), value = TRUE), grep("^PC[0-9]+$", names(centroids), value = TRUE))
  if (length(pc_cols) == 0 || nrow(centroids) == 0) return(data.frame())
  out <- lapply(seq_len(nrow(cohort_scores)), function(i) {
    z <- as.numeric(cohort_scores[i, pc_cols, drop = TRUE])
    d <- apply(centroids[, pc_cols, drop = FALSE], 1, function(cn) sqrt(sum((z - as.numeric(cn))^2, na.rm = TRUE)))
    ord <- order(d)
    nearest <- ord[[1]]
    second <- if (length(ord) >= 2) ord[[2]] else ord[[1]]
    confidence <- if (length(ord) >= 2 && d[[second]] > 0) 1 - (d[[nearest]] / d[[second]]) else NA_real_
    data.frame(
      sample_id = cohort_scores$sample_id[[i]],
      assignment_status = "assigned_by_reference_pca_centroid",
      nearest_population = centroids$population[[nearest]],
      nearest_superpopulation = centroids$superpopulation[[nearest]],
      ancestry_confidence = confidence,
      nearest_distance = d[[nearest]],
      second_nearest_population = centroids$population[[second]],
      second_nearest_superpopulation = centroids$superpopulation[[second]],
      second_nearest_distance = d[[second]],
      ancestry_note = "Atribuicao por centroide em PCs de referencia; use como ancestralidade global aproximada.",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

plot_ancestry_outputs <- function(pca_result, sample_qc, cfg, figdir) {
  manifest <- empty_figure_manifest()
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(add_figure_manifest(manifest, NA_character_, "ancestry", "ancestry", "skipped", "R package ggplot2 is not installed"))
  }

  if (nrow(pca_result$pca_scores) > 0) {
    pca_png <- file.path(figdir, "ancestry_pca.png")
    pca_pdf <- file.path(figdir, "ancestry_pca.pdf")
    p <- ancestry_pca_plot(pca_result)
    save_ggplot_pair(p, pca_png, pca_pdf, width = 8, height = 6)
    manifest <- add_figure_manifest(manifest, pca_png, "ancestry_pca_png", "ancestry", "written", pca_result$method)
    manifest <- add_figure_manifest(manifest, pca_pdf, "ancestry_pca_pdf", "ancestry", "written", pca_result$method)
  }

  if (nrow(pca_result$assignments) > 0 && any(!is.na(pca_result$assignments$nearest_superpopulation))) {
    assign_png <- file.path(figdir, "ancestry_assignment.png")
    assign_pdf <- file.path(figdir, "ancestry_assignment.pdf")
    p <- ancestry_assignment_plot(pca_result$assignments)
    save_ggplot_pair(p, assign_png, assign_pdf, width = 10, height = 5)
    manifest <- add_figure_manifest(manifest, assign_png, "ancestry_assignment_png", "ancestry", "written", "")
    manifest <- add_figure_manifest(manifest, assign_pdf, "ancestry_assignment_pdf", "ancestry", "written", "")
  }

  if (nrow(sample_qc) > 0 && "n_ancestry_snps_called" %in% names(sample_qc)) {
    qc_png <- file.path(figdir, "ancestry_snp_qc_by_sample.png")
    qc_pdf <- file.path(figdir, "ancestry_snp_qc_by_sample.pdf")
    p <- ggplot2::ggplot(sample_qc, ggplot2::aes(x = sample_id, y = n_ancestry_snps_called, fill = tumor_type)) +
      ggplot2::geom_col(width = 0.85) +
      ggplot2::labs(x = "Sample", y = "Ancestry SNPs called", fill = "Tumor type") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    save_ggplot_pair(p, qc_png, qc_pdf, width = 10, height = 5)
    manifest <- add_figure_manifest(manifest, qc_png, "ancestry_snp_qc_png", "ancestry", "written", "")
    manifest <- add_figure_manifest(manifest, qc_pdf, "ancestry_snp_qc_pdf", "ancestry", "written", "")
  }

  manifest
}

ancestry_pca_plot <- function(pca_result) {
  cohort <- pca_result$pca_scores
  ref <- pca_result$reference_scores
  cohort <- ensure_pc2_for_plot(cohort)
  ref <- ensure_pc2_for_plot(ref)
  if (nrow(ref) > 0) {
    ref$plot_group <- ifelse(is_missing_value(ref$superpopulation), ref$population, ref$superpopulation)
    cohort$plot_group <- "COHORT"
    ggplot2::ggplot() +
      ggplot2::geom_point(data = ref, ggplot2::aes(x = PC1, y = PC2, color = plot_group), alpha = 0.35, size = 1.4) +
      ggplot2::geom_point(data = cohort, ggplot2::aes(x = PC1, y = PC2), color = "black", fill = "#e15759", shape = 21, size = 3) +
      ggplot2::labs(x = "PC1", y = "PC2", color = "Reference group") +
      ggplot2::theme_minimal(base_size = 11)
  } else {
    ggplot2::ggplot(cohort, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::geom_point(color = "#2f6f9f", size = 3, alpha = 0.85) +
      ggplot2::labs(x = "PC1", y = "PC2") +
      ggplot2::theme_minimal(base_size = 11)
  }
}

ancestry_assignment_plot <- function(assignments) {
  assignments <- assignments[order(assignments$nearest_superpopulation, -assignments$ancestry_confidence), ]
  assignments$sample_id <- factor(assignments$sample_id, levels = assignments$sample_id)
  ggplot2::ggplot(assignments, ggplot2::aes(x = sample_id, y = ancestry_confidence, fill = nearest_superpopulation)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Sample", y = "Assignment confidence", fill = "Nearest reference group") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
}

ensure_pc2_for_plot <- function(x) {
  if (nrow(x) == 0) return(x)
  if (!"PC1" %in% names(x)) x$PC1 <- 0
  if (!"PC2" %in% names(x)) x$PC2 <- 0
  x
}
