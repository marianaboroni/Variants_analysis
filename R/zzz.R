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
  "locus_tumor_type_n_samples", "total_samples",
  # ggplot2 aes() symbols used in ancestry/QC plots
  "stage", "n", "reason", "PC1", "PC2", "proportion", "component", "vaf", "chrom",
  "ANCESTRY_CALL_RATE", "ANCESTRY_CONFIDENCE", "sample_id",
  # data.table/ggplot2 NSE symbols used in preparation and v2 output layers
  ".cat", ".cnt", ".cnt1", ".hist", ".site", ".sub", "ac", "af", "ALT", "an",
  "anchor_pos", "CHROM", "cosmic_count", "cosmic_id", "COSMIC_MUTATION_IDS",
  "COSMIC_OCCURRENCE_COUNT", "COSMIC_PHENOTYPE_ID", "COSMIC_SAMPLE_ID",
  "COSMIC_STUDY_ID", "COSMIC_TUMOR_BREAKDOWN", "COSMIC_TUMOR_TYPES",
  "evidence", "filter", "filter_flag", "filter_status", "final_class",
  "gene", "GENE_SYMBOL", "GENOMIC_MUTATION_ID", "HISTOLOGY_SUBTYPE_1",
  "hom", "LEGACY_MUTATION_ID", "max_pop_af", "MUTATION_SOMATIC_STATUS",
  "N", "n_variants", "orig_alt", "orig_chrom", "orig_pos", "orig_ref",
  "POS", "PRIMARY_HISTOLOGY", "PRIMARY_SITE", "PUBMED_PMID", "REF", "rsid",
  "tumor_histology", "tumor_site", "tumor_subtype", "value", "variant_type",
  "burden_label", "ccf_estimate", "clonal_fraction", "clonality_class",
  "clonality_method", "clonality_methods", "driver_class", "interpretation",
  "limitation", "n_clonal", "n_clonality_evaluable", "n_subclonal",
  "n_tmb_countable", "purity_used", "subclonal_fraction", "tmb_countable",
  "tmb_mut_per_mb"
))
