# Changelog

## Unreleased - v2 development branch

- Added real-data validation gate for `data/test/WGS_all_patients_first100k.tsv`.
- Fixed small-cohort recurrence so one-sample data are not treated as cohort
  recurrence evidence.
- Added v2 CLI entry points: `init`, `validate`, `run`, and `report`.
- Added R API wrapper `run_tumoronly()`.
- Added input validation reporting and config template generation.
- Added v2 traceability columns for supporting, opposing, and missing evidence.
- Added v2 output aliases, `run_manifest.json`, `config_used.yaml`,
  `session_info.txt`, warnings log, figure data, and QC/classification figures.
- Added regression tests for v2 interface and output contract.
- Rewrote README, quick start, installation, configuration, architecture, report
  guide, and tutorial entry points around the v2 workflow only.
- Added `config/demo_v2.yml` and installed example data under `inst/extdata` for
  reproducible source-checkout and installed-package tutorials.

This branch is not yet declared version 2.0 complete. Installation, CI,
documentation, clean-environment validation, and the full publication figure set
still require final verification.
