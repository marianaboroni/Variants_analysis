# How it works

This page describes the current v2 scientific workflow. The user-facing path is
`init -> validate -> run -> report`. Optional external resources can add
evidence, but they do not turn tumor-only candidates into clinical
confirmation.

## 1. Input validation

- Objective: determine whether the file can be analyzed safely.
- Input: VCF, MAF, TSV, TXT, or gzipped equivalent.
- Operation: detect format, map columns to the canonical schema, resolve genome
  build, derive basic metrics when possible.
- Justification: tumor-only classification is sensitive to missing depth, VAF,
  population AF, and quality evidence.
- Output: validation summary, schema report, and actionable issues.
- Limitation: validation cannot prove that upstream variant calling was correct.

## 2. Normalization

- Objective: make variant keys comparable across inputs and resources.
- Input: chromosome, position, REF, ALT, genome build.
- Operation: normalize chromosome names, split multiallelic ALT where applicable,
  trim REF/ALT to a canonical local representation.
- Output: `variant_id`, `locus_id`, canonical keys.
- Limitation: local allele normalization does not replace full reference-backed
  normalization for all complex indels.

## 3. Feature engineering

- Objective: expose evidence used by filters and scoring.
- Input: caller fields, read counts, population columns, annotations.
- Operation: derive depth, alternate reads, VAF, variant type, population AF,
  quality metrics, caller flags, and consequence severity.
- Output: standardized evidence columns.

## 4. Hard filters

- Objective: enforce fail-safe technical gates.
- Input: depth, alt reads, VAF, TLOD, MBQ, MMQ, caller `FILTER`, sample QC.
- Operation: evaluate absolute and adaptive thresholds.
- Justification: low-support variants should not be promoted by biological
  annotations.
- Output: `hard_filter_pass`, `hard_filter_reason`.

## 5. Cohort recurrence

- Objective: identify recurrence patterns that may support artifacts or, in
  adequate cohorts, tumor-type recurrence.
- Input: sample IDs and canonical variant/locus keys.
- Operation: count samples carrying each variant/locus.
- Assumption: recurrence requires enough samples to make a fraction meaningful.
- Limitation: recurrence is not evaluable in single-sample runs.

## 6. Scoring

- Objective: combine evidence into interpretable scores.
- Input: technical evidence, population evidence, VAF pattern, artifact signals,
  recurrence, and consequence.
- Operation: calculate `technical_evidence_score`, `germline_score`,
  `artifact_score`, and `somatic_score`.
- Output: bounded scores in `[0,1]`.
- Limitation: weights are documented heuristics unless calibrated on a validated
  dataset.

## 7. Classification

- Objective: assign a conservative tumor-only inference class.
- Input: scores and rule categories.
- Operation: assign `final_class` once in the audited classifier.
- Output classes: `high_confidence_somatic`, `probable_somatic`,
  `likely_germline`, `likely_artifact`, `technical_fail`,
  `manual_review_required`, `uncertain_tumor_only`.
- Interpretation: these are research classifications, not clinical truth.

## 8. Auxiliary modules

TMB/countable burden:

- Objective: summarize per-sample countable somatic-candidate burden and TMB
  when a callable territory is configured.
- Input: final class, consequence, sample ID, tumor type, `tmb.callable_mb`.
- Operation: count eligible somatic-candidate coding/splice consequences.
- Output: `tables/tmb_summary.tsv` and `tmb_countable`.
- Limitation: without `callable_mb`, TMB is not evaluable and only burden counts
  are reported.

Clonality:

- Objective: classify somatic candidates as clonal/intermediate/subclonal when
  evidence supports it.
- Input: VAF, purity, copy number, mutation multiplicity.
- Operation: estimate CCF when possible; otherwise use explicit VAF-proxy
  labels.
- Output: `tables/clonality_summary.tsv`, `ccf_estimate`,
  `clonality_class`, `clonality_method`.
- Limitation: VAF-proxy classes are not CCF estimates.

Genetic ancestry:

- Objective: estimate broad genetic-ancestry proportions for population
  frequency interpretation when a user supplies an AIMs panel.
- Input: common SNP-like loci, VAF, depth, AIMs panel.
- Operation: select germline-informative SNPs, exclude unstable loci, and
  estimate continuous proportions.
- Output: `ancestry_summary.tsv`, `ancestry_qc.tsv`, `ancestry_snps.tsv.gz`,
  and ancestry plots when enabled.
- Limitation: experimental, non-diagnostic, not race or ethnicity, and never a
  filter.

ML-assisted review:

- Objective: add review-prioritization probabilities from human-reviewed labels.
- Input: an activated model in `ml.db_dir`.
- Operation: apply the active model to run features.
- Output: `ml_status.tsv`, `ml_true_positive_probability`,
  `ml_prediction_explanation`.
- Limitation: predictions never alter filters, final classes, or scientific
  scores.

## 9. Evidence trail

- Objective: make every classification auditable.
- Input: computed evidence columns and final class.
- Operation: write supporting, opposing, and missing evidence strings plus a
  textual explanation.
- Output: `evidence_supporting_classification`,
  `evidence_against_classification`, `missing_evidence`,
  `classification_explanation`.

## 10. Outputs and report

- Objective: make the run reproducible and interpretable.
- Output: tables, figures, figure data, manifests, config, session info, and
  HTML report.
- Limitation: the current v2 figure set covers QC, classification, and module
  status through tables/report; additional publication figures are still
  planned.
