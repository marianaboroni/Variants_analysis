# Technical authenticity vs biological relevance

`R/authenticity.R` keeps two concepts strictly separate so **biological evidence
can never mask a technical artifact**:

| score | question | inputs (ONLY) |
|---|---|---|
| `TECHNICAL_AUTHENTICITY_SCORE` | is the CALL technically real? | technical_evidence_score, artifact_score, PoN, strand/orientation bias, clustered, weak_evidence |
| `BIOLOGICAL_SUPPORT_SCORE` | is the variant biologically relevant? | consequence/impact, computational predictors, hotspot, driver, COSMIC, ClinVar |
| `VARIANT_REVIEW_SCORE` | should a human look? | high biological support **and** low technical authenticity |

- `VARIANT_AUTHENTICITY_CATEGORY` ∈ {LIKELY_TRUE, UNCERTAIN, LIKELY_ARTIFACT} is
  derived **only** from `TECHNICAL_AUTHENTICITY_SCORE`. A hotspot with severe
  strand bias stays LIKELY_ARTIFACT (verified by test).
- `ARTIFACT_EVIDENCE_SCORE` is technical only.
- None of these change `filter_status`, `filter_reasons`, or `CONFIDENCE_SCORE_BASE`
  (regression-tested).

Formulas (heuristic, configurable, not clinically validated):

```
technical_authenticity = 0.6 * technical_quality + 0.4 * (1 - artifact_evidence)
biological_support     = 0.30*consequence + 0.30*recurrence + 0.25*computational + 0.15*clinical
variant_review_score   = 0.5*(1 - technical_authenticity) + 0.5*biological_support
```

Category thresholds: `authenticity.likely_true_min` (0.66), `authenticity.likely_artifact_max` (0.40).

The **adaptive per-sample** layer (`R/adaptive_filtering.R`) provides the technical
inputs (noise floor, signal-to-noise, three-zone status). See
[ADAPTIVE_FILTERING.md](ADAPTIVE_FILTERING.md).
