add_clonality_estimates <- function(x, cfg) {
  enabled <- isTRUE(cfg_get(cfg, c("clonality", "enabled"), TRUE))
  if (!enabled) {
    x$ccf_estimate <- NA_real_
    x$ccf_capped <- NA_real_
    x$clonality_class <- "not_evaluated"
    x$clonality_method <- "disabled"
    x$clonality_cluster_id <- NA_integer_
    x$clonality_cluster_center <- NA_real_
    x$clonality_cluster_label <- NA_character_
    return(list(variants = x, summary = build_clonality_summary(x)))
  }

  purity_raw <- to_numeric_safe(coalesce_columns(
    x,
    c("tumor_purity", "purity", "PURITY", "Tumor_Purity", "tumour_purity"),
    default = NA
  ))
  purity <- normalize_fraction(purity_raw)

  total_cn_raw <- to_numeric_safe(coalesce_columns(
    x,
    c("total_cn", "Total_CN", "TOTAL_CN", "CN", "copy_number", "local_cn", "tcn"),
    default = NA
  ))
  multiplicity_raw <- to_numeric_safe(coalesce_columns(
    x,
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

  somatic_ok <- x$final_class %in% c("high_confidence_somatic", "probable_somatic")
  has_vaf <- !is.na(x$vaf) & x$vaf >= 0 & x$vaf <= 1
  has_purity <- !is.na(purity) & purity > 0 & purity <= 1
  has_cn <- !is.na(total_cn_raw) & total_cn_raw > 0
  has_multiplicity <- !is.na(multiplicity_raw) & multiplicity_raw > 0

  denom <- purity * total_cn + (1 - purity) * normal_cn
  ccf <- (x$vaf * denom) / (purity * multiplicity)
  ccf[!somatic_ok | !has_vaf | !has_purity] <- NA_real_
  ccf[!is.finite(ccf)] <- NA_real_

  x$purity_used <- purity
  x$total_cn_used <- total_cn
  x$multiplicity_used <- multiplicity
  x$copy_number_available <- has_cn
  x$multiplicity_available <- has_multiplicity
  x$ccf_estimate <- ccf
  x$ccf_capped <- pmin(ccf, 1)
  x$clonality_method <- "not_evaluable"
  x$clonality_method[somatic_ok & has_vaf & has_purity & has_cn] <- "purity_copy_number_adjusted"
  x$clonality_method[somatic_ok & has_vaf & has_purity & !has_cn] <- "purity_adjusted_copy_neutral"
  x$clonality_method[somatic_ok & has_vaf & !has_purity] <- "vaf_proxy_no_purity"

  x$clonality_class <- classify_clonality(x, cfg)
  x <- add_clonality_clusters(x, cfg)

  list(
    variants = x,
    summary = build_clonality_summary(x)
  )
}

classify_clonality <- function(x, cfg) {
  clonal_ccf <- cfg_get(cfg, c("clonality", "clonal_ccf_cutoff"), 0.85)
  subclonal_ccf <- cfg_get(cfg, c("clonality", "subclonal_ccf_cutoff"), 0.55)
  clonal_vaf <- cfg_get(cfg, c("clonality", "clonal_vaf_proxy_cutoff"), 0.30)
  subclonal_vaf <- cfg_get(cfg, c("clonality", "subclonal_vaf_proxy_cutoff"), 0.12)

  out <- rep("not_evaluable", nrow(x))
  somatic_ok <- x$final_class %in% c("high_confidence_somatic", "probable_somatic")
  has_ccf <- somatic_ok & !is.na(x$ccf_estimate)
  has_vaf_proxy <- somatic_ok & is.na(x$ccf_estimate) & !is.na(x$vaf)

  out[has_ccf & x$ccf_estimate >= clonal_ccf] <- "clonal"
  out[has_ccf & x$ccf_estimate < clonal_ccf & x$ccf_estimate > subclonal_ccf] <- "intermediate"
  out[has_ccf & x$ccf_estimate <= subclonal_ccf] <- "subclonal"

  out[has_vaf_proxy & x$vaf >= clonal_vaf] <- "clonal_like_high_vaf"
  out[has_vaf_proxy & x$vaf < clonal_vaf & x$vaf > subclonal_vaf] <- "intermediate_vaf"
  out[has_vaf_proxy & x$vaf <= subclonal_vaf] <- "subclonal_like_low_vaf"
  out
}

add_clonality_clusters <- function(x, cfg) {
  x$clonality_cluster_id <- NA_integer_
  x$clonality_cluster_center <- NA_real_
  x$clonality_cluster_label <- NA_character_

  max_k <- cfg_get(cfg, c("clonality", "max_clusters_per_sample"), 3)
  min_variants <- cfg_get(cfg, c("clonality", "min_variants_for_clustering"), 10)
  somatic_ok <- x$final_class %in% c("high_confidence_somatic", "probable_somatic")
  cluster_value <- ifelse(!is.na(x$ccf_capped), x$ccf_capped, x$vaf)
  eligible <- somatic_ok & !is.na(cluster_value)

  for (sample_id in unique(x$sample_id[eligible])) {
    idx <- which(eligible & x$sample_id == sample_id)
    if (length(idx) < min_variants) next

    values <- cluster_value[idx]
    values <- pmin(pmax(values, 0), 1)
    k <- choose_clonality_k(values, max_k)
    set.seed(cfg_get(cfg, c("clonality", "seed"), 20260621))
    fit <- stats::kmeans(values, centers = k, nstart = 25)
    centers <- as.numeric(fit$centers)
    rank_high_to_low <- rank(-centers, ties.method = "first")

    x$clonality_cluster_id[idx] <- rank_high_to_low[fit$cluster]
    x$clonality_cluster_center[idx] <- centers[fit$cluster]
    x$clonality_cluster_label[idx] <- cluster_label_from_center(centers[fit$cluster])
  }
  x
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

build_clonality_summary <- function(x) {
  if (!"clonality_class" %in% names(x)) {
    return(data.frame())
  }
  dt <- data.table::as.data.table(x)
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
  as.data.frame(out)
}

collapse_unique <- function(x, max_items = 8) {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(x[seq_len(min(length(x), max_items))], collapse = ";")
}
