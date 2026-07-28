# V2 final report - current status

Status: not ready for version 2.0 release.

This file records the current implementation status on branch
`feature/cohort-ingestion-graded-filters`. It is intentionally conservative:
the real-data gate has passed, and the first v2 user interface/output layers are
implemented, but final clean-environment installation and the complete
publication figure/documentation acceptance suite are not finished.

## Commits in this work segment

| Commit | Purpose |
|---|---|
| `bc9f5c9` | validated the real WGS TSV and fixed Phase 1 blockers |
| `173cdc3` | documented the v2 workflow and input contract |
| `b05cd2c` | added v2 CLI/API entry points |
| `f17774b` | added traceable v2 outputs, manifest, evidence columns, and figures |
| `d6f3acc` | added installation docs, Docker/Conda/CI files, man pages, and check fixes |
| `300833d` | rewrote README/tutorial/quickstart/docs around the current v2 workflow only |
| auxiliary module stage | restored TMB/burden, clonality, ancestry status, and activated-model ML as functional v2 auxiliary modules |

## Real-data validation

Input:

```text
data/test/WGS_all_patients_first100k.tsv
```

Latest auxiliary-module validation run:

```bash
/usr/bin/time -l Rscript exec/tumoronly run \
  --input data/test/WGS_all_patients_first100k.tsv \
  --config config/wgs_all_patients_subset_config.yml \
  --output results/wgs_all_patients_subset_filter \
  --run-id real_data_modules_v2_callable_gate
```

Results:

| Metric | Value |
|---|---:|
| Input rows | 100,000 |
| Standardized variants | 100,000 |
| PASS | 51,531 |
| REVIEW | 4,996 |
| FAIL | 43,473 |
| Pipeline runtime in `run_manifest.json` | 744.62 s |
| Wall time | 762.28 s |
| Maximum resident set size | 853,487,616 bytes |
| Peak memory footprint | 3,082,560,448 bytes |

Warnings:

| Code | Affected variants |
|---|---:|
| `W_COHORT_RECURRENCE_NOT_EVALUABLE` | 100,000 |
| `W_TLOD_MISSING` | 100,000 |
| `W_COSMIC_NOT_CONFIGURED` | 100,000 |
| `W_DRIVER_RESOURCES_NOT_CONFIGURED` | 100,000 |
| `W_GENE_ANNOTATION_MISSING` | 33,275 |
| `W_TMB_CALLABLE_MB_MISSING` | 100,000 |
| `W_CLONALITY_VAF_PROXY_ONLY` | 51,879 |
| `W_ML_NO_ACTIVE_MODEL` | 100,000 |
| `W_ANCESTRY_PANEL_MISSING` | 100,000 |

Auxiliary module statuses on the real file:

| Module | Status | Interpretation |
|---|---|---|
| TMB/burden | `not_evaluable_missing_callable_mb` | 106 countable coding candidates; no mutations/Mb value was reported because no validated callable territory denominator was configured |
| Clonality | `evaluated` | 51,879 somatic-candidate rows received VAF-proxy labels because purity was absent |
| ML | `not_evaluable_no_active_model` | no activated reviewed-label model was available |
| Ancestry | `not_evaluable` | no AIMs marker panel was configured |

The auxiliary modules did not change `final_class`, `filter_status`,
`somatic_score`, or the fail-safe scientific rules.

## Implemented interface changes

- `tumoronly init`
- `tumoronly validate`
- `tumoronly run`
- `tumoronly report`
- `run_tumoronly()`
- `validate_tumoronly_input()`
- `write_tumoronly_config_template()`
- `tumoronly_default_config()`

The CLI calls the R API and does not duplicate scientific classification logic.

## Implemented output contract

Each completed run now writes:

- `report.html`
- `run_manifest.json`
- `config_used.yaml`
- `session_info.txt`
- `logs/warnings.tsv`
- `tables/all_variants.tsv`
- `tables/classified_variants.tsv`
- `tables/high_confidence_somatic.tsv`
- `tables/likely_germline.tsv`
- `tables/likely_artifact.tsv`
- `tables/known_drivers.tsv`
- `tables/sample_summary.tsv`
- `tables/filter_audit.tsv`
- `tables/module_status.tsv`
- `tables/tmb_summary.tsv`
- `tables/clonality_summary.tsv`
- `tables/ml_status.tsv`
- `tables/ancestry_summary.tsv` when ancestry is enabled
- `figures/*.pdf`
- `figures/*.svg`
- `figures/*.png`
- `figure_data/*.tsv`

Every row in `tables/classified_variants.tsv` includes:

- `evidence_supporting_classification`
- `evidence_against_classification`
- `missing_evidence`
- `classification_explanation`

## Figures implemented

The first v2 figure set is generated as PDF, SVG, and 300 dpi PNG, with matching
TSV data:

1. filtering workflow;
2. QC and score distributions;
3. final classification distribution;
4. VAF by depth;
5. classification burden by sample;
6. missing or limited evidence.

This is not yet the complete publication figure set requested for v2.0.

## Tests run

Targeted v2 tests:

```text
Rscript -e 'testthat::test_file("tests/testthat/test-v2-interface.R", reporter = "summary")'
Result: 0 failures

Rscript -e 'testthat::test_file("tests/testthat/test-v2-modules.R", reporter = "summary")'
Result: 0 failures

Rscript -e 'testthat::test_file("tests/testthat/test-review-ml.R", reporter = "summary")'
Result: 0 failures

Rscript -e 'testthat::test_file("tests/testthat/test-ancestry.R", reporter = "summary")'
Result: 0 failures
```

Full development test suite:

```text
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
Result: 0 failures, 1 skipped, 7 known COSMIC/FASTA warnings
```

Package build/check:

```text
R CMD build .
R CMD check --no-manual --no-build-vignettes tumoronly_0.1.0.tar.gz
Result: Status OK
```

Clean local installation from tarball:

```text
R CMD INSTALL --library=<temp_lib> tumoronly_0.1.0.tar.gz
library(tumoronly)
system.file("exec", "tumoronly", package = "tumoronly")
Rscript <installed_cli> validate --input data/test/WGS_all_patients_first100k.tsv --config config/wgs_all_patients_subset_config.yml
Result: package loaded, CLI found, validation PASS with documented warnings
```

Real-data CLI validation:

```text
Rscript exec/tumoronly validate \
  --input data/test/WGS_all_patients_first100k.tsv \
  --config config/wgs_all_patients_subset_config.yml
Result: PASS with warnings
```

## Installation and distribution files added

- `INSTALLATION.md`
- `QUICKSTART.md`
- `TUTORIAL.md`
- `INPUT_FORMAT.md`
- `CONFIGURATION.md`
- `HOW_IT_WORKS.md`
- `INTERPRETING_RESULTS.md`
- `TROUBLESHOOTING.md`
- `CHANGELOG.md`
- `CONTRIBUTING.md`
- `Dockerfile`
- `environment.yml`
- `.github/workflows/R-CMD-check.yaml`
- `config/demo_v2.yml`
- `inst/extdata/demo_variants.tsv`
- `inst/extdata/demo_sample_metadata.tsv`

The README, quick start, tutorial, installation guide, configuration guide,
architecture guide, report guide, and vignette entry point have been rewritten
around `init`, `validate`, `run`, and `report`. TMB/burden, clonality, ancestry,
and ML are documented as functional auxiliary modules that write explicit
status/output tables and do not override the audited core classifier.

These files have been added, but a clean installation from each route still
needs final execution for Conda, Docker, and remote GitHub Actions. Local
tarball installation and installed-CLI validation have passed.

Documentation validation after this stage:

```text
Rscript exec/tumoronly init --output /tmp/tumoronly_doc_config.yml --force
Result: OK

Rscript exec/tumoronly validate --input data/demo/variants.tsv --config config/demo_v2.yml --output results/demo_v2_validation.html
Result: PASS with expected W_COHORT_RECURRENCE_NOT_EVALUABLE warning

Rscript exec/tumoronly run --input data/demo/variants.tsv --config config/demo_v2.yml --output results/demo_v2 --run-id quickstart_demo
Result: completed; report.html, tables, figures, figure_data, logs, and manifests written

Clean tarball install into a temporary library, followed by API and installed-CLI demo using inst/extdata
Result: OK

Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
Result: 0 failures, 1 skipped, 7 known COSMIC/FASTA warnings

R CMD check --no-manual --no-build-vignettes tumoronly_0.1.0.tar.gz
Result: Status OK
```

## Limitations still open

- The real-data file has one active sample, so it does not validate true
  multi-sample cohort recurrence behavior.
- COSMIC and driver resources were not configured in the real-data validation.
- TLOD and parsed strand-bias evidence are missing from the real TSV.
- TMB in mutations/Mb was not evaluable on the real file because no validated
  callable territory denominator was configured; only countable burden was
  reported.
- Clonality on the real file is VAF-proxy only because tumor purity was absent.
- ML was not evaluable on the real file because no reviewed-label model was
  activated.
- Genetic ancestry was not evaluable on the real file because no AIMs marker
  panel was configured.
- The report is improved by the new output contract but has not yet been fully
  rewritten into the complete 16-section v2 report specification.
- The full publication figure set is not complete.
- Conda, Docker, and remote GitHub Actions runs remain to be executed.

## Readiness recommendation

Not ready for v2.0 release or publication yet.

Ready to continue v2 development because the real-data gate passed and the
core CLI/API/output contract now runs on the real file without changing the
audited scientific classification counts.
