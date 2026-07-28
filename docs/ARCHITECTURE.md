# tumoronly v2 Architecture

The v2 architecture keeps one scientific workflow and exposes it through thin
user-facing wrappers.

## Entry Points

CLI:

```text
Rscript exec/tumoronly init
Rscript exec/tumoronly validate
Rscript exec/tumoronly run
Rscript exec/tumoronly report
```

R API:

```r
tumoronly_default_config()
write_tumoronly_config_template()
validate_tumoronly_input()
run_tumoronly()
```

The CLI calls these R functions and does not duplicate scientific logic.

## Core Flow

```text
read_variant_input()
  -> standardize_variant_table()
  -> add_basic_features()
  -> add_sample_qc_and_hard_filters()
  -> compute_cohort_recurrence()
  -> add_adaptive_filtering()
  -> score_variants()
  -> annotate_cosmic()
  -> add_driver_annotation_layers()
  -> add_guideline_classifications()
  -> classify_variants()
  -> add_prioritization()
  -> add_variant_authenticity()
  -> add_tumoronly_evidence_columns()
  -> write tables, figures, manifests, report
```

`add_tumoronly_evidence_columns()` is a v2 traceability layer. It does not
change `final_class`, scores, thresholds, or filter decisions.

## Input Layer

`read_variant_input()` accepts VCF, MAF, TSV, TXT, and gzipped variants of those
formats. It emits:

- the preserved source columns;
- canonical fields such as `CHROM`, `POS`, `REF`, `ALT`, `SAMPLE_ID`;
- a schema report;
- column mapping metadata.

`validate_tumoronly_input()` uses the same ingestion code as `run_tumoronly()`.
Validation is therefore a true preflight of the run path, not a separate parser.

## Output Layer

The original workflow still writes compatibility tables such as
`variants_all.tsv.gz`. The v2 layer additionally writes the user-facing contract:

```text
report.html
tables/all_variants.tsv
tables/classified_variants.tsv
tables/high_confidence_somatic.tsv
tables/likely_germline.tsv
tables/likely_artifact.tsv
tables/known_drivers.tsv
tables/sample_summary.tsv
tables/filter_audit.tsv
figures/*.pdf|*.svg|*.png
figure_data/*.tsv
logs/warnings.tsv
config_used.yaml
session_info.txt
run_manifest.json
```

## Scientific Boundaries

- Tumor-only classes are inference candidates, not clinical truth.
- Matched-normal absence is always recorded as missing evidence.
- Cohort recurrence is disabled below the configured minimum sample count.
- Missing COSMIC, driver, TLOD, population-frequency, mapping-quality, or
  strand-bias evidence remains visible in logs and per-variant explanation
  columns.
- Archived TMB, clonality, ancestry, dashboard, and ML workflows are not part of
  the v2 core user workflow.

## Real-Data Gate

Before v2 development continued, the workflow was run end to end on
`data/test/WGS_all_patients_first100k.tsv`. The validation report is
`docs/REAL_DATA_VALIDATION.md`.
