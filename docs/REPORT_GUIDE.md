# v2 Report Guide

Every v2 run writes a root-level report:

```text
results/<run_id>/report.html
```

The legacy template path `report/tumor_only_report.html` may still be present
inside a run directory for compatibility, but `report.html` is the user-facing
v2 entry point.

## Expected Sections

The current report summarizes:

- input data and cohort;
- validation issues and warnings;
- quality-control metrics;
- filtering workflow;
- classification counts;
- auxiliary module status for TMB/burden, clonality, ancestry, and ML;
- prioritized variants;
- per-variant evidence columns;
- limitations and reproducibility metadata.

## Evidence Language

Report text must use candidate language:

- `classified as high-confidence somatic candidate`;
- `consistent with a likely germline event`;
- `flagged as a likely technical artifact`;
- `requires matched-normal or orthogonal confirmation for clinical use`.

It must not claim that a tumor-only call is clinically confirmed somatic or
germline.

## Source Files

The report should be interpreted with:

- `tables/classified_variants.tsv`
- `tables/filter_audit.tsv`
- `tables/sample_summary.tsv`
- `tables/module_status.tsv`
- `tables/tmb_summary.tsv`
- `tables/clonality_summary.tsv`
- `tables/ml_status.tsv`
- `tables/ancestry_summary.tsv`, when ancestry is enabled
- `logs/warnings.tsv`
- `run_manifest.json`
- `config_used.yaml`
- `session_info.txt`

Figure source data are stored in `figure_data/*.tsv`.

## Current Limitation

The current v2 report is functional and generated automatically, but the full
publication-style 16-section narrative report remains a release-readiness item
tracked in `docs/V2_FINAL_REPORT.md`.
