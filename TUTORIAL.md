# tumoronly v2 Tutorial

This tutorial runs the current v2 workflow end to end. It uses only the
documented API and CLI:

```text
init -> validate -> run -> report
```

The example is intentionally small. It is useful for learning the workflow and
checking installation, not for validating scientific performance.

## 1. Concept

Tumor-only analysis starts from tumor variant calls without a matched normal.
That makes interpretation difficult because a technically credible variant can
be somatic, germline-like, recurrent technical noise, or simply under-supported.

`tumoronly` therefore reports conservative candidates:

- `high_confidence_somatic` or `probable_somatic` means compatible with a
  tumor-only somatic-candidate rule.
- `likely_germline` means consistent with germline-like evidence, usually
  population frequency or VAF pattern.
- `likely_artifact` or `technical_fail` means technical evidence argues against
  using the call as a somatic candidate.
- `manual_review_required` and `uncertain_tumor_only` mean the evidence is
  incomplete or conflicting.

Do not rewrite these classes as confirmed clinical calls.

## 2. Source-Checkout Tutorial

From the repository root:

```bash
Rscript exec/tumoronly validate \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2_validation.html
```

Validation should return `PASS`. Warnings are possible in the demo because the
dataset is intentionally tiny and has no matched normal or COSMIC database.

Run the workflow:

```bash
Rscript exec/tumoronly run \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2 \
  --run-id tutorial_demo
```

Re-render or copy the report:

```bash
Rscript exec/tumoronly report \
  --run-dir results/demo_v2/tutorial_demo \
  --output results/demo_v2/tutorial_demo/report.html
```

Open the report:

```bash
open results/demo_v2/tutorial_demo/report.html
```

On Linux, use `xdg-open`; on Windows, open the file from Explorer.

## 3. Installed-Package Tutorial

After installing the package, the example data are available through
`system.file()`. This R script writes a runnable config, validates the input,
runs the CLI, and prints the run directory.

```r
library(tumoronly)
library(yaml)

variants <- system.file("extdata", "demo_variants.tsv", package = "tumoronly")
metadata <- system.file("extdata", "demo_sample_metadata.tsv", package = "tumoronly")

cfg <- tumoronly_default_config()
cfg$input$path <- variants
cfg$input$genome_build <- "GRCh38"
cfg$input$sample_metadata <- metadata
cfg$analysis$output_dir <- "results/installed_demo"
cfg$cohort$min_cohort_size_for_recurrence <- 10
cfg$cancer$default_tumor_type <- "UCEC"
yaml::write_yaml(cfg, "installed_demo_config.yml")

validation <- validate_tumoronly_input(
  config = "installed_demo_config.yml",
  output = "installed_demo_validation.html"
)
stopifnot(isTRUE(validation$ok))

result <- run_tumoronly(
  config = "installed_demo_config.yml",
  run_id = "installed_demo"
)

result$run_dir
file.path(result$run_dir, "report.html")
```

The installed CLI can run the same config:

```bash
TO_CLI="$(Rscript -e 'cat(system.file("exec", "tumoronly", package = "tumoronly"))')"
Rscript "$TO_CLI" validate --config installed_demo_config.yml --output installed_cli_validation.html
Rscript "$TO_CLI" run --config installed_demo_config.yml --run-id installed_cli_demo
```

## 4. What Validation Checks

Validation reads the file without changing scientific values. It reports:

- detected format;
- genome build;
- mapped canonical columns;
- number of rows, samples, variants, and duplicate sample-variant keys;
- availability of depth and VAF;
- missing indispensable fields;
- warnings for missing evidence such as TLOD, population frequency, gene
  annotation, or cohort recurrence.

Strict validation blocks execution when indispensable evidence is absent.
Permissive validation records issues but does not fill in missing data:

```bash
Rscript exec/tumoronly validate \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --permissive
```

## 5. What `run` Does

The run performs these steps:

1. Ingest VCF, MAF, TSV, TXT, or gzipped inputs.
2. Resolve and check genome build.
3. Standardize chromosomes, positions, REF, ALT, sample ID, gene, consequence,
   and metrics.
4. Derive depth, alternate reads, VAF, population AF, variant type, and quality
   evidence where possible.
5. Apply fail-safe technical filters.
6. Evaluate cohort recurrence only when enough samples exist.
7. Calculate technical, germline, artifact, and somatic evidence scores.
8. Apply the audited classifier once.
9. Add driver and hotspot evidence when resources are configured.
10. Write evidence strings explaining every classification.
11. Generate tables, figures, logs, config snapshots, manifests, and report.

## 6. Output Tour

Start with:

```text
results/<run_id>/report.html
results/<run_id>/run_manifest.json
results/<run_id>/logs/warnings.tsv
results/<run_id>/tables/classified_variants.tsv
results/<run_id>/tables/sample_summary.tsv
results/<run_id>/tables/filter_audit.tsv
```

Then inspect per-class tables:

```text
tables/high_confidence_somatic.tsv
tables/likely_germline.tsv
tables/likely_artifact.tsv
tables/known_drivers.tsv
```

Each row in `classified_variants.tsv` carries:

- `evidence_supporting_classification`;
- `evidence_against_classification`;
- `missing_evidence`;
- `classification_explanation`.

These columns are the audit trail. If a row has missing COSMIC, TLOD, mapping
quality, or cohort evidence, that limitation remains visible.

## 7. Figures

The current v2 figure set writes PDF, SVG, PNG, and source data:

```text
figures/figure_01_filtering_workflow.*
figures/figure_02_qc_overview.*
figures/figure_03_classification.*
figures/figure_04_vaf_depth.*
figures/figure_05_sample_class_distribution.*
figures/figure_06_missing_evidence.*
figure_data/figure_01_filtering_workflow.tsv
...
```

The TSV files in `figure_data/` are the reproducible source for custom plotting.

## 8. Real-Data Gate

The demo is not the validation gate. The v2 implementation was first tested on
`data/test/WGS_all_patients_first100k.tsv`. See
`docs/REAL_DATA_VALIDATION.md` for counts, warnings, runtime, memory, and
traceable example variants.

The real-data gate passed for ingestion and full workflow execution. It remains
limited because the 100k-row slice contains one active sample, so real
multi-sample recurrence behavior still requires an additional real cohort.

## 9. Common Problems

- `E_INPUT_MISSING`: pass `--input` or set `input.path`.
- `E_GENOME_BUILD`: set `input.genome_build` to `GRCh37` or `GRCh38`.
- `E_STANDARDIZE`: add `input.column_map`.
- `E_DEPTH_MISSING`, `E_ALT_COUNT_MISSING`, or `E_VAF_MISSING`: provide
  indispensable technical evidence or deliberately rerun permissively.
- `W_COHORT_RECURRENCE_NOT_EVALUABLE`: expected for single-sample or tiny demo
  runs.
- `W_COSMIC_NOT_CONFIGURED`: COSMIC evidence was not available; classification
  still runs, but this evidence is marked missing.

