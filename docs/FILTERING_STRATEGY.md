# Filtering & prioritization strategy

How `tumoronly` decides `filter_status` (retain/exclude) and `PRIORITY_SCORE_BASE`
(ranking). Filtering removes unsupported/likely-non-somatic calls; prioritization
ranks what remains. Predictors, COSMIC context and OncoKB never remove a variant.

| Stage | Objetivo | Campos/evidências | Estratégia | Efeito | Status |
|---|---|---|---|---|---|
| Technical QC | Remove unsupported calls | FILTER, QUAL, DP, AD, VAF, TLOD, MBQ, MMQ | adaptive + configurable thresholds | **can exclude** (`technical_fail`) | pipeline rule |
| Artifact evidence | Flag suspicious events | caller FILTER, PoN, strand/orientation bias, clustered, weak_evidence | technical rules | **can exclude** (`likely_artifact`) | pipeline rule |
| Population frequency | Reduce likely germline | gnomAD/ExAC/1000G/ABraOM AFs | AF thresholds (common/low) | **can exclude/flag** (`likely_germline`) | pipeline rule |
| Germline disposition | Identify germline-like pattern | AF + VAF + ClinVar (ACMG) | conservative ACMG | exclude or review | guideline-adapted |
| Consequence | Rank functional impact | Consequence, IMPACT, LoF | severity hierarchy | **prioritizes** | pipeline rule |
| Functional predictors | Add computational evidence | SIFT/PolyPhen/REVEL/CADD/AlphaMissense/SpliceAI/… | family-aware consensus | **prioritizes only** | heuristic |
| Driver/hotspot | Prioritize relevant genes/positions | curated drivers, cancer hotspots | curated DBs | **prioritizes** | curated |
| COSMIC global | Record recurrence | COSMIC_TOTAL_OCCURRENCES | global evidence | complementary | evidence |
| COSMIC contextual | Same-tumor recurrence | site/histology/subtype | contextual matching | **experimental** | experimental |
| OncoKB | Post-hoc interpretation | OncoKB API | confirmatory | **never filters** | post-hoc |

## What sets `filter_status`

`filter_status` ∈ {PASS, REVIEW, FAIL} is derived **only** from the rule-based
`final_class` (one decision point, `scoring.R::classify_variants`), which uses:
technical QC gate, population category, artifact category, guideline oncogenicity
(incl. tumor-context OS4/OM4 — never global-only), and germline disposition.

**Never sets `filter_status`:** functional predictors, `COMPUTATIONAL_EVIDENCE_*`,
COSMIC contextual score, `CONFIDENCE_SCORE_WITH_COSMIC_EXPERIMENTAL`, OncoKB,
prioritization.

## Examples

- A variant with low depth, few alt reads, or strand bias is **excluded** (technical).
- A variant with high population AF is **flagged/excluded** as likely germline.
- SIFT and PolyPhen damaging do **not** remove or keep a variant — they add
  computational evidence and raise priority.
- Hotspots and driver genes raise priority but do not bypass technical QC.

## Prioritization (heuristic, decomposed)

`PRIORITY_SCORE_BASE = 0.25·consequence + 0.20·computational + 0.25·hotspot +
0.15·driver + 0.10·clinical + 0.05·cosmic_global`, each component in [0,1] and
exported separately (`PRIORITY_*_COMPONENT`, `PRIORITY_COMPONENTS`). Weights are
**heuristic** until calibrated. Priority ranks retained variants; it never filters.
