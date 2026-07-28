# Interpreting results

Start with `report.html`, then inspect the tables.

## Main classes

| Class | Meaning |
|---|---|
| `high_confidence_somatic` | high-confidence somatic candidate under tumor-only rules |
| `probable_somatic` | probable somatic candidate with supportive technical and population evidence |
| `likely_germline` | consistent with a likely germline-like event |
| `likely_artifact` | flagged as likely technical artifact |
| `technical_fail` | excluded by fail-safe technical filters |
| `manual_review_required` | conflicting or incomplete evidence |
| `uncertain_tumor_only` | insufficient evidence for a more specific class |

Do not write "this variant is somatic" or "this variant is germline" from a
tumor-only run. Use candidate/inference language.

## Per-variant evidence

Read these columns together:

- `evidence_supporting_classification`
- `evidence_against_classification`
- `missing_evidence`
- `classification_explanation`

A useful result is not one without missing evidence. A useful result is one
where the missing evidence is explicit.

## Warnings

Run-level warnings are in `logs/warnings.tsv` and `run_manifest.json`.

Examples:

| Warning | Interpretation |
|---|---|
| `W_COHORT_RECURRENCE_NOT_EVALUABLE` | too few samples for cohort recurrence |
| `W_TLOD_MISSING` | caller-specific quality evidence incomplete |
| `W_COSMIC_NOT_CONFIGURED` | COSMIC evidence unavailable |
| `W_DRIVER_RESOURCES_NOT_CONFIGURED` | driver interpretation limited |

## Figures

Each v2 figure has:

- PDF
- SVG
- 300 dpi PNG
- a matching `figure_data/*.tsv`

The figure data are authoritative for reproduction and custom plotting.

## Clinical limitations

Tumor-only analysis cannot confirm somatic or germline status. Confirmation
requires matched normal, orthogonal validation, or a clinically validated
workflow appropriate to the question.

