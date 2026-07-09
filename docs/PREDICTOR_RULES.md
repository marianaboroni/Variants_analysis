# Functional predictor rules

Predictors are handled by `R/predictors.R` driven by
`inst/config/predictor_registry.yml`. They contribute **computational evidence**
and one prioritization component; they **never** set `filter_status`
(`use_for_filtering: false` for all).

## Pipeline

1. **Detection** — for each registry predictor, aliases are matched against input
   columns → `tables/predictor_inventory.tsv` (coverage, observed range, status).
2. **Normalization** — canonical `*_SCORE`/`*_PRED` columns; the VEP `label(score)`
   form (e.g. `deleterious(0.01)`) is split; SpliceAI DS_* aggregated as max;
   out-of-range values dropped (recorded, never trusted).
3. **Applicability** — by variant consequence scope:
   - **missense**: SIFT, PolyPhen, REVEL, CADD, AlphaMissense, MetaLR/RNN/SVM,
     M-CAP, PrimateAI, ClinPred, BayesDel, MutationAssessor, FATHMM, PROVEAN.
   - **splicing**: SpliceAI, dbscSNV, MaxEntScan.
   - **nonsense/frameshift/canonical splice** (`truncating`): high-impact
     consequence used directly → `PREDICTOR_APPLICABILITY_STATUS = not_applicable`
     for missense predictors (SIFT/PolyPhen not required).
   - **synonymous**: SpliceAI + conservation + ClinVar.
   - Missing applicable predictors → `applicable_but_missing` (never benign).
4. **Family-aware consensus** — predictors are grouped into families
   (`protein_function`, `ensemble_missense`, `splicing`, `conservation`,
   `clinical_database`). Correlated predictors within a family count as **one**
   family, so SIFT+PolyPhen alone are not two independent lines of evidence.

## Computational evidence outputs

`COMPUTATIONAL_EVIDENCE_SCORE` = mean over present families of their damaging
magnitude (0 if none damaging). `COMPUTATIONAL_EVIDENCE_CATEGORY`:

| category | rule |
|---|---|
| `strong_support` | ≥ `minimum_damaging_predictors` families damaging, no benign family, concordant |
| `moderate_support` | ≥ min damaging families but some benign |
| `weak_support` | ≥1 damaging predictor but < min families |
| `conflicting` | damaging and benign among applicable predictors (≈balanced) |
| `insufficient` | applicable but no usable predictor value |
| `not_applicable` | consequence uses impact directly (e.g. frameshift) |

Also: `COMPUTATIONAL_PREDICTORS_AVAILABLE/DAMAGING/BENIGN/CONFLICTING`,
`COMPUTATIONAL_CONFLICT_FLAG/REASON`, `COMPUTATIONAL_EVIDENCE_DETAILS`.

## Config

```yaml
predictors:
  detect_automatically: true
  use_for_prioritization: true
  use_for_filtering: false        # enforced regardless
  mode: consensus
  minimum_applicable_predictors: 2
  minimum_damaging_predictors: 2  # families, not raw predictors
  missing_policy: ignore
  conflict_policy: report
```

**All thresholds are heuristic** and not clinically validated; predictor output is
never labeled "pathogenic" on its own.
