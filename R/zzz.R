# Mark the package as data.table-aware so that non-standard evaluation
# (.(), .N, by=) inside `[.data.table` resolves correctly when the package is
# loaded via library() rather than dev-sourced into the global environment.
.datatable.aware <- TRUE

# Silence R CMD check "no visible binding for global variable" NOTEs arising
# from data.table NSE column references used across the pipeline.
utils::globalVariables(c(
  ".", ".N", ":=", "sample_id", "tumor_type", "vaf", "alt_count", "dp",
  "variant_id", "chrom", "pos", "ref", "alt", "locus_id", "canonical_key",
  "Hugo_Symbol", "Tumor_Sample_Barcode",
  "variant_cohort_freq", "variant_n_samples", "locus_cohort_freq",
  "locus_n_samples", "variant_tumor_type_freq", "variant_tumor_type_n_samples",
  "tumor_type_total_samples", "locus_tumor_type_freq",
  "variant_median_vaf", "variant_median_dp", "variant_median_alt_count",
  "locus_tumor_type_n_samples", "total_samples"
))
