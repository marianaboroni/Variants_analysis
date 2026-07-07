# COSMIC decision rules (tumor-type-stratified)

COSMIC evidence is interpreted **in the context of the sample's tumor type**, and
is decomposed into separate, traceable components — never a single opaque score.
This document defines the harmonization, the controlled statuses, the granular
match levels, the decomposed scores (bands + weights), the base-vs-experimental
confidence separation, and what each situation may affect.

## Harmonization

The sample's `tumor_type`/`tumor_subtype` (and, when available, `primary_site`,
`histology`) are mapped to COSMIC categories through an **explicit, versioned**
table — never fuzzy text matching: `config/cosmic_tumor_type_mapping.tsv`
(`input_tumor_type`, `input_tumor_subtype`, `cosmic_primary_site`,
`cosmic_histology`, `match_level` ∈ {`exact`,`compatible`,`broad`},
`mapping_version`, `mapping_source`). Original terms are preserved; normalized
COSMIC terms are used for matching. Different histologies of the same organ are
**not** merged unless an explicit `compatible`/`broad` row says so.

## Storage (P6)

The processed DB stores a **normalized long table** `cosmic_context_long`
(`canonical_key`, `cosmic_primary_site`, `cosmic_histology`, `cosmic_subtype`,
`occurrence_count`, `cosmic_release`, `source_row_count`) — one row per tumor
category. This is the **primary analytical structure** for matching and
aggregation; delimiters are irrelevant because categories are columns, not a
concatenated string. `COSMIC_TUMOR_BREAKDOWN` remains only as a delimiter-safe
summary string for export.

## Global vs matched-tumor recurrence (P1)

Two separate scores, so a high global count can never stand in for absence in the
sample's tumor type:

```
COSMIC_GLOBAL_RECURRENCE_SCORE          = recurrence_transform(COSMIC_TOTAL_OCCURRENCES)
COSMIC_MATCHED_TUMOR_RECURRENCE_SCORE   = recurrence_transform(COSMIC_MATCHING_TUMOR_OCCURRENCES)
recurrence_transform(count)             = min(1, log10(count + 1) / 3)
```

`recurrence_transform` is a **transformed raw recurrence count**. It reduces the
influence of extreme values; it is **not** a frequency, prevalence, or a
correction for the uneven number of samples represented per tumor type in COSMIC
(that would require per-tumor denominators, which raw COSMIC does not provide). If
a COSMIC release supplies denominators, a proper rate should be implemented
separately with its source documented.

## Granular match (P8)

```
COSMIC_SITE_MATCH        COSMIC_HISTOLOGY_MATCH        COSMIC_SUBTYPE_MATCH
COSMIC_MATCH_LEVEL ∈ { exact_site_histology_subtype, exact_site_histology,
                       compatible_histology, broad_site, global_only, none }
```

A subtype-level exact match is declared **only** when the sample subtype and a
COSMIC subtype are both present and concordant via a subtype-specific mapping row.
When subtype is unavailable, no subtype-level exact match is claimed.

## Controlled statuses (P7)

| status | definition |
|---|---|
| `exact_match` | matched in same site+histology (subtype when available), not broadly pan-cancer |
| `compatible_match` | matched in a compatible category (same site, related/broad histology) |
| `pan_cancer_with_sample_tumor` | recurrent across ≥ threshold categories **including** the sample type |
| `pan_cancer_without_sample_tumor` | recurrent across ≥ threshold categories but **not** the sample type (global evidence only) |
| `other_tumor_only` | present only in other tumor type(s) |
| `tumor_type_unknown` | sample tumor type absent/not harmonized — **not evaluable** |
| `cosmic_tumor_type_missing` | COSMIC record lacks tumor context — **not evaluable** |
| `no_cosmic_match` | no valid genomic match |

Pan-cancer is split so a variant recurrent in many tumors **but absent in the
sample's type** never gets the same contextual status as one that **includes** the
sample type; the former remains global evidence only.

## Non-evaluable context (P2)

`tumor_type_unknown` and `cosmic_tumor_type_missing` do **not** receive positive
contextual evidence:

```
COSMIC_TUMOR_SPECIFICITY_SCORE = NA
COSMIC_CONTEXT_SUPPORT_SCORE   = NA
COSMIC_CONTEXT_EVALUABLE       = FALSE
```

Global evidence (`COSMIC_GLOBAL_RECURRENCE_SCORE`, occurrence counts) is still
reported separately, but absence of tumor information is never turned into
positive contextual support.

## Banded support score (P3) — ordering guaranteed for ANY counts

`COSMIC_CONTEXT_SUPPORT_SCORE` uses **non-overlapping bands**; recurrence only
modulates *within* a band, so tiers can never be crossed by raw counts:

| status | band [floor, ceiling] | recurrence used within band |
|---|---|---|
| `exact_match` | [0.80, 1.00] | matched-tumor |
| `compatible_match` | [0.55, 0.75] | matched-tumor |
| `pan_cancer_with_sample_tumor` | [0.35, 0.50] | matched-tumor |
| `other_tumor_only` / `pan_cancer_without_sample_tumor` | [0.05, 0.25] | global |
| not evaluable | NA | — |

`support = floor + (ceiling - floor) * recurrence`. Therefore, always:

```
exact_match > compatible_match > pan_cancer_with_sample_tumor > other_tumor_only
```

and an `other_tumor_only` variant with thousands of occurrences (≤ 0.25) can never
outrank an `exact_match` with a single occurrence (≥ 0.80).

`COSMIC_TUMOR_SPECIFICITY_SCORE` is the pure, recurrence-independent tier value
(`exact_match`=1.0, `compatible_match`=0.65, `pan_cancer_with_sample_tumor`=0.45,
`other_tumor_only`/`pan_cancer_without_sample_tumor`=0.15, non-evaluable=NA).

## Base vs experimental confidence (P4)

The COSMIC context score does **not** silently replace the main confidence:

```
CONFIDENCE_SCORE_BASE                        # main output & filters — NO COSMIC context
CONFIDENCE_CATEGORY_BASE
COSMIC_CONTEXT_SUPPORT_SCORE                 # banded, tumor-context
CONFIDENCE_SCORE_WITH_COSMIC_EXPERIMENTAL    # base + 0.10 * context support (EXPERIMENTAL)
CONFIDENCE_CATEGORY_WITH_COSMIC_EXPERIMENTAL
```

`CONFIDENCE_SCORE_BASE`/`CONFIDENCE_CATEGORY_BASE` and `filter_status` are the
authoritative outputs and are **invariant to tumor context** (regression-tested).
The `_WITH_COSMIC_EXPERIMENTAL` columns are complementary/exploratory.

**Note on what still uses global COSMIC:** the guideline oncogenicity criteria
(OS4/OM4) use the **global** `COSMIC_TOTAL_OCCURRENCES` (genomic recurrence),
which can influence `final_class` (documented, pre-existing). The **tumor-context**
score (matching occurrences, specificity) does **not** feed `final_class` or
`CONFIDENCE_SCORE_BASE`.

## What each situation may affect

| Situação | Evidência genômica | Contexto tumoral | Efeito permitido | Pode alterar `filter_status`? | Pode alterar confiança? |
| --- | --- | --- | --- | --- | --- |
| Mesmo tipo e mesma histologia | Match completo | Exato | Suporte forte | Não (só via OS4/OM4 global, documentado) | Sim (experimental) |
| Mesmo órgão, outra histologia | Match completo | Parcial | Suporte limitado | Não automaticamente | Sim, peso menor (experimental) |
| Apenas outro tipo tumoral | Match completo | Discordante | Informação contextual | Não | No máximo suporte fraco (banda ≤0.25) |
| Pan-câncer incluindo o tumor | Match completo | Multitumoral c/ amostra | Suporte de recorrência | Não automaticamente | Sim (experimental) |
| Pan-câncer sem o tumor | Match completo | Multitumoral s/ amostra | Evidência global apenas | Não | Fraco (banda ≤0.25) |
| Tipo tumoral desconhecido | Match completo | Não avaliável | Evidência global apenas | Não | Não (specificity=NA) |
| Sem match COSMIC | Sem match | Não aplicável | Nenhum suporte COSMIC | Não excluir automaticamente | Não |

## OS4 / OM4 are tumor-context gated (audit)

The only path by which global COSMIC recurrence could reach `final_class` was the
guideline OS4/OM4 criteria, which previously used the **global**
`COSMIC_TOTAL_OCCURRENCES`. They are now granted from **tumor-context** occurrences
only:

```
OS4_CONTEXTUAL = COSMIC_CONTEXT_EVALUABLE & includes_sample &
                 COSMIC_MATCHING_TUMOR_OCCURRENCES >= OS4_MIN_MATCHED_TUMOR_OCCURRENCES (default 10)
OM4_CONTEXTUAL = COSMIC_CONTEXT_EVALUABLE & includes_sample &
                 OM4_MIN_MATCHED <= COSMIC_MATCHING_TUMOR_OCCURRENCES < OS4_MIN (defaults 3..10)
```

`other_tumor_only`, `pan_cancer_without_sample_tumor`, `tumor_type_unknown` and
`cosmic_tumor_type_missing` set `COSMIC_GLOBAL_EVIDENCE_ONLY = TRUE` and never grant
OS4/OM4. See [COSMIC_OS4_OM4_AUDIT.md](COSMIC_OS4_OM4_AUDIT.md) and the regression
in [regression/COSMIC_CONTEXT_CLASSIFICATION_REGRESSION.tsv](regression/COSMIC_CONTEXT_CLASSIFICATION_REGRESSION.tsv).

## Provenance / status of each rule

| Rule | Status |
|---|---|
| ACMG/AMP germline criteria (BA1/BS1/PM2 …) | guideline-based (published) |
| Somatic oncogenicity OVS1/OS3/OM3/OP4/SBVS1 … | guideline-based (ClinGen/CGC-style), adapted |
| **OS4/OM4 tumor-context gating + `OS4/OM4_MIN_MATCHED_TUMOR_OCCURRENCES`** | **heuristic pipeline thresholds** (not clinically validated) |
| Tumor-type harmonization table | curated mapping, versioned |
| `COSMIC_CONTEXT_SUPPORT_SCORE` bands | tumor-only-adapted heuristic |
| `CONFIDENCE_SCORE_WITH_COSMIC_EXPERIMENTAL` | **experimental score** |
| `filter_status`, `CONFIDENCE_SCORE_BASE` | pipeline rule-based (unchanged core thresholds) |

Heuristic thresholds are **not** presented as validated clinical criteria.

## Sensitivity & calibration status (P9/P12)

`cosmic_weight_sensitivity()` recomputes an alternative **linear** score
(`w_rec·global_recurrence + w_ctx·specificity`) under 20/80, 40/60 and 50/50
weightings and reports rank correlation and category churn. Because plausible
weight changes move variants across categories, the COSMIC context score is
**exploratory/complementary**, not a validated component; the main pipeline uses
`CONFIDENCE_SCORE_BASE`.
