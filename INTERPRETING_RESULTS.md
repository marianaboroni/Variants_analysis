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
- `tmb_countable`
- `clonality_class`
- `clonality_method`
- `ml_status`

A useful result is not one without missing evidence. A useful result is one
where the missing evidence is explicit.

## Auxiliary modules

Start with `tables/module_status.tsv`.

| Module | Main table | How to interpret |
|---|---|---|
| TMB/burden | `tmb_summary.tsv` | TMB in mutations/Mb only when `callable_mb` is configured; otherwise burden counts only |
| Clonality | `clonality_summary.tsv` | CCF-based when purity/CN evidence exists; VAF-proxy labels are explicitly marked |
| Ancestry | `ancestry_summary.tsv` | experimental, non-diagnostic genetic ancestry from user-supplied AIMs panel |
| ML | `ml_status.tsv` | activated-model predictions for review prioritization only |

None of these modules should be used to override `final_class` without an
explicit scientific change and new validation.

## Warnings

Run-level warnings are in `logs/warnings.tsv` and `run_manifest.json`.

Examples:

| Warning | Interpretation |
|---|---|
| `W_COHORT_RECURRENCE_NOT_EVALUABLE` | too few samples for cohort recurrence |
| `W_TLOD_MISSING` | caller-specific quality evidence incomplete |
| `W_COSMIC_NOT_CONFIGURED` | COSMIC evidence unavailable |
| `W_DRIVER_RESOURCES_NOT_CONFIGURED` | driver interpretation limited |
| `W_TMB_CALLABLE_MB_MISSING` | burden was counted, but TMB is not evaluable |
| `W_CLONALITY_VAF_PROXY_ONLY` | clonality used VAF-proxy labels because purity was missing |
| `W_ML_NO_ACTIVE_MODEL` | ML is enabled but no activated model was available |
| `W_ANCESTRY_PANEL_MISSING` | ancestry is enabled but the AIMs panel is unavailable |

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
