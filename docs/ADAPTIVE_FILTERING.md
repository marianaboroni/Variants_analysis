# Adaptive per-sample filtering (EXPERIMENTAL)

`R/adaptive_filtering.R`. Complements universal hard floors with per-sample
technical profiles, an empirical noise floor, a probabilistic alt-read error test,
and a three-zone decision. Technical only — biological relevance never relaxes it.

## Per-sample technical profile

`SAMPLE_MEDIAN_DP`, `SAMPLE_MEAN_DP`, `SAMPLE_DP_MAD/IQR`, `SAMPLE_DP_P05..P95`,
`SAMPLE_MEDIAN_ALT_DEPTH`, `SAMPLE_BACKGROUND_VAF`, `SAMPLE_NOISE_FLOOR`,
`SAMPLE_CALLABLE_FRACTION`, `SAMPLE_LOW_COVERAGE_FRACTION` (robust stats;
noise floor from the low-VAF tail, bounded to [0.001, 0.02]).

## Adaptive thresholds (never below absolute floors)

```
ADAPTIVE_MIN_DP = max(absolute_floor, min(DP_P05, median - k*MAD))   # method-dependent
ADAPTIVE_MIN_ALT_DEPTH = absolute_floor                              # per-variant SNR test does the rest
ADAPTIVE_MIN_VAF = SAMPLE_NOISE_FLOOR * minimum_above_noise_factor
```

Per variant: `EXPECTED_ERROR_ALT_READS = dp*noise_floor`,
`AD_ALT_SIGNAL_TO_NOISE`, `ALT_READ_ERROR_PVALUE` (upper-tail binomial vs noise),
`EXPECTED_MIN_DETECTABLE_VAF`, `VAF_ABOVE_SAMPLE_NOISE`, `VAF_DETECTION_CONFIDENCE`.

## Three zones (`ADAPTIVE_FILTER_STATUS`)

- **FAIL**: below absolute floor, compatible with noise (binomial p > 0.05),
  strong bias (>0.9), or PoN.
- **REVIEW**: below adaptive DP, low signal-to-noise, low VAF but above noise
  (possible subclonal), within `review_margin` of a threshold, moderate bias, or a
  REVIEW/FAIL `SAMPLE_QC_STATUS`.
- **PASS**: above adaptive thresholds with adequate support.

## Guardrails (enforced / tested)

Adaptive thresholds never go below `hard_safety_limits`; `AD_ALT=1`/`DP=6` still
FAIL; biology never relaxes technical requirements; globally poor samples get
`SAMPLE_QC_STATUS ∈ {FAIL,REVIEW}` rather than lowered thresholds; the module never
learns TP/FP from its own PASS/FAIL. The documented escalation
(`apply_adaptive_review`) only moves borderline PASS → REVIEW and never changes
`CONFIDENCE_SCORE_BASE` (tested).

## Outputs

`tables/adaptive_thresholds.tsv` (per-sample manifest: median/MAD/percentiles,
noise floor, adaptive & absolute thresholds, method, version, QC, warnings).

## Config

See `config/example.yml` `adaptive_filtering:` (method `percentile_mad` by
default). Start with this simple, auditable baseline; mixture-model / isolation-
forest / supervised alternatives are future comparisons. All thresholds heuristic;
validate against ground truth / matched normal / replicates before production.
