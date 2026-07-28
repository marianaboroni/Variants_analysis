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
| packaging/docs stage | adds installation docs, Docker/Conda/CI files, man pages, and check fixes |

## Real-data validation

Input:

```text
data/test/WGS_all_patients_first100k.tsv
```

Final post-commit validation run:

```bash
/usr/bin/time -l Rscript exec/tumoronly run \
  --input data/test/WGS_all_patients_first100k.tsv \
  --config config/wgs_all_patients_subset_config.yml \
  --output results/wgs_all_patients_subset_filter \
  --run-id real_data_phase4_v2_outputs_final
```

Results:

| Metric | Value |
|---|---:|
| Input rows | 100,000 |
| Standardized variants | 100,000 |
| PASS | 51,531 |
| REVIEW | 4,996 |
| FAIL | 43,473 |
| Pipeline runtime in `run_manifest.json` | 130.59 s |
| Wall time | 139.53 s |
| Maximum resident set size | 1,439,297,536 bytes |
| Peak memory footprint | 3,066,564,288 bytes |

Warnings:

| Code | Affected variants |
|---|---:|
| `W_COHORT_RECURRENCE_NOT_EVALUABLE` | 100,000 |
| `W_TLOD_MISSING` | 100,000 |
| `W_COSMIC_NOT_CONFIGURED` | 100,000 |
| `W_DRIVER_RESOURCES_NOT_CONFIGURED` | 100,000 |
| `W_GENE_ANNOTATION_MISSING` | 33,275 |

The final `run_manifest.json` records commit `f17774b`.

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

These files have been added, but a clean installation from each route still
needs final execution for Conda, Docker, and remote GitHub Actions. Local
tarball installation and installed-CLI validation have passed.

## Limitations still open

- The real-data file has one active sample, so it does not validate true
  multi-sample cohort recurrence behavior.
- COSMIC and driver resources were not configured in the real-data validation.
- TLOD and parsed strand-bias evidence are missing from the real TSV.
- The report is improved by the new output contract but has not yet been fully
  rewritten into the complete 16-section v2 report specification.
- The full publication figure set is not complete.
- Conda, Docker, and remote GitHub Actions runs remain to be executed.

## Readiness recommendation

Not ready for v2.0 release or publication yet.

Ready to continue v2 development because the real-data gate passed and the
core CLI/API/output contract now runs on the real file without changing the
audited scientific classification counts.
