# Clonality auxiliary module.
#
# CCF is estimated only when VAF and purity are available. Copy-number and
# multiplicity improve the estimate when supplied; otherwise the module records
# the fallback method. Without purity, the output is explicitly a VAF proxy.

#' Add clonality estimates to classified variants.
#'
#' @param variants Classified tumoronly variant table.
#' @param cfg Resolved tumoronly configuration.
#' @return A list with updated `variants` and per-sample `summary`.
#' @export
add_clonality_estimates <- function(variants, cfg = NULL) {
  if (!is.data.frame(variants)) stop("add_clonality_estimates(): `variants` must be a data.frame.", call. = FALSE)
  enabled <- isTRUE(cfg_get(cfg, c("clonality", "enabled"), TRUE))
  if (!enabled) {
    variants <- add_empty_clonality_columns(variants, method = "disabled")
    return(list(variants = variants, summary = build_clonality_summary(variants)))
  }

  purity_raw <- to_numeric_safe(coalesce_columns(
    variants,
    c("tumor_purity", "purity", "PURITY", "Tumor_Purity", "tumour_purity"),
    default = NA
  ))
  purity <- normalize_fraction(purity_raw)

  total_cn_raw <- to_numeric_safe(coalesce_columns(
    variants,
    c("total_cn", "Total_CN", "TOTAL_CN", "CN", "copy_number", "local_cn", "tcn"),
    default = NA
  ))
  multiplicity_raw <- to_numeric_safe(coalesce_columns(
    variants,
    c("multiplicity", "mutation_multiplicity", "mut_cn", "mutation_cn"),
    default = NA
  ))

  default_total_cn <- cfg_get(cfg, c("clonality", "default_total_cn"), 2)
  default_multiplicity <- cfg_get(cfg, c("clonality", "default_multiplicity"), 1)
  normal_cn <- cfg_get(cfg, c("clonality", "normal_cn"), 2)

  total_cn <- total_cn_raw
  total_cn[is.na(total_cn) | total_cn <= 0] <- default_total_cn
  multiplicity <- multiplicity_raw
  multiplicity[is.na(multiplicity) | multiplicity <= 0] <- default_multiplicity

  somatic_ok <- variants$final_class %in% c("high_confidence_somatic", "probable_somatic")
  has_vaf <- !is.na(variants$vaf) & variants$vaf >= 0 & variants$vaf <= 1
  has_purity <- !is.na(purity) & purity > 0 & purity <= 1
  has_cn <- !is.na(total_cn_raw) & total_cn_raw > 0
  has_multiplicity <- !is.na(multiplicity_raw) & multiplicity_raw > 0

  denom <- purity * total_cn + (1 - purity) * normal_cn
  ccf <- (variants$vaf * denom) / (purity * multiplicity)
  ccf[!somatic_ok | !has_vaf | !has_purity] <- NA_real_
  ccf[!is.finite(ccf)] <- NA_real_

  variants$purity_used <- purity
  variants$total_cn_used <- total_cn
  variants$multiplicity_used <- multiplicity
  variants$copy_number_available <- has_cn
  variants$multiplicity_available <- has_multiplicity
  variants$ccf_estimate <- ccf
  variants$ccf_capped <- pmin(ccf, 1)
  variants$clonality_method <- "not_evaluable"
  variants$clonality_method[somatic_ok & has_vaf & has_purity & has_cn] <- "purity_copy_number_adjusted"
  variants$clonality_method[somatic_ok & has_vaf & has_purity & !has_cn] <- "purity_adjusted_copy_neutral"
  variants$clonality_method[somatic_ok & has_vaf & !has_purity] <- "vaf_proxy_no_purity"

  variants$clonality_class <- classify_clonality(variants, cfg)
  variants <- add_clonality_clusters(variants, cfg)

  list(variants = variants, summary = build_clonality_summary(variants))
}

add_empty_clonality_columns <- function(variants, method = "not_evaluated") {
  variants$purity_used <- NA_real_
  variants$total_cn_used <- NA_real_
  variants$multiplicity_used <- NA_real_
  variants$copy_number_available <- FALSE
  variants$multiplicity_available <- FALSE
  variants$ccf_estimate <- NA_real_
  variants$ccf_capped <- NA_real_
  variants$clonality_class <- "not_evaluable"
  variants$clonality_method <- method
  variants$clonality_cluster_id <- NA_integer_
  variants$clonality_cluster_center <- NA_real_
  variants$clonality_cluster_label <- NA_character_
  variants
}

classify_clonality <- function(variants, cfg) {
  clonal_ccf <- cfg_get(cfg, c("clonality", "clonal_ccf_cutoff"), 0.85)
  subclonal_ccf <- cfg_get(cfg, c("clonality", "subclonal_ccf_cutoff"), 0.55)
  clonal_vaf <- cfg_get(cfg, c("clonality", "clonal_vaf_proxy_cutoff"), 0.30)
  subclonal_vaf <- cfg_get(cfg, c("clonality", "subclonal_vaf_proxy_cutoff"), 0.12)

  out <- rep("not_evaluable", nrow(variants))
  somatic_ok <- variants$final_class %in% c("high_confidence_somatic", "probable_somatic")
  has_ccf <- somatic_ok & !is.na(variants$ccf_estimate)
  has_vaf_proxy <- somatic_ok & is.na(variants$ccf_estimate) & !is.na(variants$vaf)

  out[has_ccf & variants$ccf_estimate >= clonal_ccf] <- "clonal"
  out[has_ccf & variants$ccf_estimate < clonal_ccf & variants$ccf_estimate > subclonal_ccf] <- "intermediate"
  out[has_ccf & variants$ccf_estimate <= subclonal_ccf] <- "subclonal"

  out[has_vaf_proxy & variants$vaf >= clonal_vaf] <- "clonal_like_high_vaf"
  out[has_vaf_proxy & variants$vaf < clonal_vaf & variants$vaf > subclonal_vaf] <- "intermediate_vaf"
  out[has_vaf_proxy & variants$vaf <= subclonal_vaf] <- "subclonal_like_low_vaf"
  out
}

add_clonality_clusters <- function(variants, cfg) {
  variants$clonality_cluster_id <- NA_integer_
  variants$clonality_cluster_center <- NA_real_
  variants$clonality_cluster_label <- NA_character_

  max_k <- cfg_get(cfg, c("clonality", "max_clusters_per_sample"), 3)
  min_variants <- cfg_get(cfg, c("clonality", "min_variants_for_clustering"), 10)
  somatic_ok <- variants$final_class %in% c("high_confidence_somatic", "probable_somatic")
  cluster_value <- ifelse(!is.na(variants$ccf_capped), variants$ccf_capped, variants$vaf)
  eligible <- somatic_ok & !is.na(cluster_value)

  for (sid in unique(variants$sample_id[eligible])) {
    idx <- which(eligible & variants$sample_id == sid)
    if (length(idx) < min_variants) next

    values <- pmin(pmax(cluster_value[idx], 0), 1)
    k <- choose_clonality_k(values, max_k)
    set.seed(cfg_get(cfg, c("clonality", "seed"), 20260621))
    fit <- stats::kmeans(values, centers = k, nstart = 25)
    centers <- as.numeric(fit$centers)
    rank_high_to_low <- rank(-centers, ties.method = "first")

    variants$clonality_cluster_id[idx] <- rank_high_to_low[fit$cluster]
    variants$clonality_cluster_center[idx] <- centers[fit$cluster]
    variants$clonality_cluster_label[idx] <- cluster_label_from_center(centers[fit$cluster])
  }
  variants
}

choose_clonality_k <- function(values, max_k) {
  values <- values[!is.na(values)]
  n <- length(values)
  max_k <- max(1, min(max_k, n, length(unique(round(values, 4)))))
  if (max_k <= 1) return(1)

  scores <- rep(Inf, max_k)
  for (k in seq_len(max_k)) {
    if (k == 1) {
      wss <- sum((values - mean(values))^2)
    } else {
      fit <- stats::kmeans(values, centers = k, nstart = 10)
      wss <- fit$tot.withinss
    }
    wss <- max(wss, .Machine$double.eps)
    scores[[k]] <- n * log(wss / n) + k * log(n)
  }
  which.min(scores)
}

cluster_label_from_center <- function(center) {
  ifelse(center >= 0.85, "high_ccf_cluster",
         ifelse(center <= 0.55, "low_ccf_cluster", "intermediate_cluster"))
}

build_clonality_summary <- function(variants) {
  if (!"clonality_class" %in% names(variants)) {
    return(data.frame())
  }
  dt <- data.table::as.data.table(variants)
  somatic_classes <- c("high_confidence_somatic", "probable_somatic")
  out <- dt[, .(
    n_somatic_variants = sum(final_class %in% somatic_classes, na.rm = TRUE),
    n_clonality_evaluable = sum(final_class %in% somatic_classes & clonality_class != "not_evaluable", na.rm = TRUE),
    n_clonal = sum(clonality_class %in% c("clonal", "clonal_like_high_vaf"), na.rm = TRUE),
    n_intermediate = sum(clonality_class %in% c("intermediate", "intermediate_vaf"), na.rm = TRUE),
    n_subclonal = sum(clonality_class %in% c("subclonal", "subclonal_like_low_vaf"), na.rm = TRUE),
    median_somatic_vaf = safe_median(vaf[final_class %in% somatic_classes]),
    median_ccf = safe_median(ccf_estimate[final_class %in% somatic_classes]),
    median_purity_used = safe_median(purity_used),
    clonality_methods = collapse_unique(clonality_method[final_class %in% somatic_classes])
  ), by = .(sample_id, tumor_type)]
  out[, clonal_fraction := ifelse(n_clonality_evaluable > 0, n_clonal / n_clonality_evaluable, NA_real_)]
  out[, subclonal_fraction := ifelse(n_clonality_evaluable > 0, n_subclonal / n_clonality_evaluable, NA_real_)]
  out[, limitation := ifelse(grepl("vaf_proxy_no_purity", clonality_methods %||% ""),
                             "No purity was available; clonality classes are VAF-proxy labels, not CCF estimates.",
                             "CCF interpretation depends on purity, copy number, and mutation multiplicity accuracy.")]
  as.data.frame(out)
}

collapse_unique <- function(x, max_items = 8) {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(x[seq_len(min(length(x), max_items))], collapse = ";")
}
