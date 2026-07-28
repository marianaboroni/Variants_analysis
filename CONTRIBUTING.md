# Contributing

## Scientific safety rules

- Do not silently change classification thresholds.
- Do not reintroduce archived TMB modules or out-of-scope features.
- Do not present tumor-only output as clinical confirmation.
- Preserve per-variant supporting, opposing, and missing evidence.
- Add regression tests for every scientific fix.

## Development setup

```bash
R CMD INSTALL .
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

Run focused tests while editing:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-v2-interface.R")'
```

## Before committing

- Run relevant targeted tests.
- Run the full test suite for shared behavior.
- Update documentation when user-facing behavior changes.
- Keep commits small and descriptive.

## Pull request checklist

- No changes to scientific rules without documented rationale.
- Input validation still blocks missing required fields.
- CLI and R API remain equivalent for shared parameters.
- `run_manifest.json` records parameters, input checksums, warnings, and status.
- Real-data validation remains documented in `docs/REAL_DATA_VALIDATION.md`.

