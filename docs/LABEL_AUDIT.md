# Label audit — four distinct concepts, never conflated

The interface uses four independent labels; detail lives in `*_reason` columns.

| Label | Values | Meaning | Source column |
|---|---|---|---|
| **STATUS** | PASS / REVIEW / FAIL | filtering decision (technical + rule-based class) | `filter_status` |
| **PRIORITY** | HIGH / MEDIUM / LOW / NOT_PRIORITIZED | biological interest (ranking) | `PRIORITY_CATEGORY` |
| **EVIDENCE** | STRONG / MODERATE / WEAK / CONFLICTING / INSUFFICIENT / NOT_APPLICABLE | computational predictor consensus | `COMPUTATIONAL_EVIDENCE_CATEGORY` |
| **AUTHENTICITY** | LIKELY_TRUE / UNCERTAIN / LIKELY_ARTIFACT | technical real-vs-artifact | `VARIANT_AUTHENTICITY_CATEGORY` |

A variant can be e.g. `STATUS=PASS`, `PRIORITY=LOW`, `EVIDENCE=INSUFFICIENT`,
`AUTHENTICITY=LIKELY_ARTIFACT` — the four are orthogonal by design (tested).

## Reason / detail columns

`status_reason` (= filter_reasons), `failed_filters`, `review_reasons`
(why REVIEW: flagged_for_review / technical_artifact_evidence /
conflicting_predictors / biologically_interesting_but_technically_weak /
possible_brazilian_germline / adaptive_borderline_technical), `priority_reasons`
(= PRIORITY_COMPONENTS), `evidence_summary` (= COMPUTATIONAL_EVIDENCE_DETAILS).

## Prior (internal) labels retained for traceability

`final_class`, `COSMIC_TUMOR_CONTEXT_STATUS`, `POPULATION_EVIDENCE_STATUS`,
`ADAPTIVE_FILTER_STATUS`, `VARIANT_AUTHENTICITY_SCORE`, etc. remain in
`variants_all.tsv.gz` for auditability; the four labels above are the primary
interface. Confirmed review labels (TRUE_POSITIVE/FALSE_POSITIVE/UNCERTAIN/
NOT_REVIEWED) come only from human review (see docs/REPORT_GUIDE.md and R/review.R).
