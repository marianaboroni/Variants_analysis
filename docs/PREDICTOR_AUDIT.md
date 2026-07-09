# Functional predictor audit

Grounded in the real inputs shipped/tested: the Mutect2+VEP VCF
(`data/test/tumor_ATPBR-…VEP.ann.vcf.gz`, 98 CSQ fields) and the VEP TSV
(`data/test/WGS_all_patients_first100k.tsv`). "Used" means the value influences a
downstream column, not merely that a column exists.

## Predictors detected in the real data

| Predictor / field | Origin | Variant scope | Present? | Parsed (before refactor)? | Used (before)? | Affects filter? | Affects prioritization (after refactor)? |
|---|---|---|---|---|---|---|---|
| Consequence | VEP CSQ | all | yes | yes | yes | **yes** (consequence hierarchy → impact_rank, is_truncating) | yes (consequence component) |
| IMPACT | VEP CSQ | all | yes | yes (impact_rank) | yes | indirectly (via classification) | yes |
| SIFT (`deleterious(0.01)`) | VEP CSQ | missense | yes | **no** (not split) | no | no | **yes** (protein_function family) |
| PolyPhen (`probably_damaging(0.99)`) | VEP CSQ | missense | yes | **no** | no | no | **yes** |
| AlphaMissense_pred/score | dbNSFP/VEP plugin | missense | yes | score only | yes (functional_impact_score → driver) | no | yes (ensemble_missense) |
| MetaLR_pred/score | dbNSFP | missense | yes | score only (meta max) | yes (driver) | no | yes |
| MetaRNN_pred/score | dbNSFP | missense | yes | score only (meta max) | yes (driver) | no | yes |
| MutationTaster_pred/score | dbNSFP | missense/splice | yes | score only (meta max) | yes (driver) | no | yes |
| SpliceAI_pred_DS_AG/AL/DG/DL | SpliceAI plugin | splicing | yes | max DS parsed | yes (is_splice_disruptive, driver) | no (`is_splice_disruptive` feeds driver, not filter) | yes (splicing family) |
| CLIN_SIG | VEP CSQ (ClinVar) | all | yes | yes (clinvar_significance) | yes (guideline SBP2/OP5, germline) | **can** (via ACMG germline path) | yes (clinical_database) |
| REVEL, CADD, FATHMM, PROVEAN, PrimateAI, BayesDel, ClinPred, MutationAssessor, M-CAP, MetaSVM, dbscSNV, MaxEntScan, phyloP, phastCons, GERP | dbNSFP/plugins | various | **absent in this dataset** | aliases exist for REVEL/CADD | no | no | supported via registry when present |

## Findings (gaps the refactor closes)

1. **SIFT/PolyPhen were never parsed** — VEP emits `label(score)`; the old code
   never split them, so SIFT/PolyPhen contributed nothing. Fixed: normalization
   splits `deleterious(0.01)` → `SIFT_PRED=deleterious`, `SIFT_SCORE=0.01`.
2. **Predictor use was implicit** — `add_functional_predictor_features`
   (driver_classification.R) hard-codes a subset into `functional_impact_score`
   which feeds `driver_score` only. There was no inventory, no applicability, no
   consensus, no computational-evidence columns.
3. **No applicability gating** — scores were read regardless of consequence
   (e.g. SIFT on a frameshift is meaningless). Fixed with
   `PREDICTOR_APPLICABILITY_STATUS`.
4. **Absence was implicitly benign** in the mean-based score. Fixed: missing
   predictors are `applicable_but_missing`, never counted as benign.
5. **No filter leak** — confirmed predictors do **not** set `filter_status`
   directly (they fed `driver_score`/`functional_impact_score` only). The refactor
   keeps this invariant and adds a regression test
   (`filter_status` identical before/after adding predictors).

## After the refactor

Predictors are handled by a **registry-driven** layer (`R/predictors.R` +
`inst/config/predictor_registry.yml`):
- dynamic detection → `tables/predictor_inventory.tsv`;
- normalization → canonical `*_SCORE`/`*_PRED` columns;
- applicability by consequence → `PREDICTOR_APPLICABILITY_STATUS`;
- family-aware consensus → `COMPUTATIONAL_EVIDENCE_*` (never a filter);
- prioritization uses `COMPUTATIONAL_EVIDENCE_SCORE` as one decomposed component.

Rules and thresholds are in [PREDICTOR_RULES.md](PREDICTOR_RULES.md); all
thresholds are **heuristic** until calibrated.
