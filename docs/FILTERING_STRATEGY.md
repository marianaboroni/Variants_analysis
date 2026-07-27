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

## Caller FILTER as graded evidence (Mutect2 tumor-only)

By default the technical gate treats the caller `FILTER` column as binary:
only the literal string `PASS` (or missing/`.`/empty) passes; anything else
counts as a caller-filter failure everywhere the caller FILTER is consulted
(`utils.R::caller_filter_pass()` — one shared helper used by
`hard_filters.R::apply_variant_hard_filters()` and
`compute_sample_qc_from_variants()`, `scoring.R::quality_artifact_score()` and
`classify_technical_artifact_evidence()`, and
`features.R::technical_floor_pass()`, so the allowlist below is consistent
everywhere, not just at the top-level gate). For Mutect2 tumor-only, this
means a variant flagged `germline` or `panel_of_normals` by `FilterMutectCalls`
can **never** reach `high_confidence_somatic`/`probable_somatic` regardless of
how strong the rest of the evidence is (`classify_variants()` forces
`technical_fail` on `!hard_pass` before any other evidence is consulted).

That default is too blunt for tumor-only Mutect2 calls where `germline`/
`panel_of_normals` are graded evidence, not proof: Mutect2's `germline` tag
comes from GERMQ, a comparison against `--germline-resource` (typically
gnomAD) that can be wrong under **loss of heterozygosity (LOH)** — a somatic
event at a locus that also carries a common germline SNP can push the VAF into
a genotype-like pattern indistinguishable, without a matched normal, from true
germline. This is not an edge case for ovarian HGSOC, where LOH at driver loci
(`TP53`, `BRCA1/2`, `PTEN`) is near-universal.

Two config-driven, backward-compatible extensions (default behavior unchanged
unless configured):

```yaml
hard_filters:
  caller_filter_accepted_values: ["PASS", "germline", "panel_of_normals", "germline;panel_of_normals"]
  caller_filter_germline_values: ["germline", "germline;panel_of_normals"]
```

- `caller_filter_accepted_values` (`hard_filters.R`): which exact FILTER
  strings pass the technical gate. Variants outside this list should already
  have been excluded **upstream** (e.g. a Mutect2-tumor-only harmonization
  step keeping only PASS/germline/panel_of_normals/combinations before this
  point) — this key only controls what the technical gate does with what it
  is handed, it does not itself filter a raw multi-artifact VCF.
- `caller_filter_germline_values` (`features.R::add_basic_features()`): which
  exact FILTER strings set `caller_germline_flag`. This flag is **one weighted
  component of `germline_score`** (`scoring.R::score_variants()`), alongside
  `population_germline_score` (ABraOM/gnomAD AF — the largest weight, since it
  is sample-independent), `vaf_germline_score`, and `recurrent_germline_signal`
  — never a veto on its own. A variant Mutect2 flagged `germline` but that is
  absent from population databases and sits in a driver-gene hotspot can still
  reach `high_confidence_somatic`; a truly common population SNP still lands
  in `likely_germline` regardless of this flag, because population evidence
  dominates the weighted sum.
- `panel_of_normals` already had a graded path before this change:
  `pon_flag` (`features.R`, from an `PON`/`panel_of_normals` column) feeds
  `artifact_score` and several classification rules. It only needed
  `caller_filter_accepted_values` to stop being short-circuited to
  `technical_fail` first.

## Cohort-wide ingestion

`compute_cohort_recurrence()` (`features.R`) is what makes `cohort.
recurrent_variant_fraction_artifact` (used by `recurrent_germline_signal()`,
`quality_artifact_score()`, `classify_technical_artifact_evidence()` and
`classify_variants()`) a real cross-sample signal rather than an inert
config key: it counts, for every variant/locus, how many *distinct*
`sample_id` values in the table handed to that `run` carry it, out of the
total distinct `sample_id` count in that same table (`total_samples`). It runs
**once per `run_tumor_only()` call**, over whatever `read_variant_input()`
returned for that call — so it is only as cohort-wide as the input actually is.

A single tumor-only Mutect2 VCF has exactly one sample column, so running
`tumoronly run` once per sample (the common case — see docs/GETTING_STARTED.md,
Paso 7) makes `total_samples == 1` in every one of those runs:
`variant_cohort_freq` degenerates to 1.0 for everything present and the
recurrence signal is inert, even though the config key is set.

To get real cohort-wide recurrence across several tumor-only VCFs in one run,
`input.path` accepts either a single path (unchanged) or a YAML list of paths:

```yaml
input:
  path:
    - mutect2/sample1/sample1.mutect2.filtered_VEP.ann.vcf.gz
    - mutect2/sample2/sample2.mutect2.filtered_VEP.ann.vcf.gz
    - mutect2/sample3/sample3.mutect2.filtered_VEP.ann.vcf.gz
```

Each file is ingested independently through `read_variant_input()` (its own
format detection, its own column mapping — files are never assumed to share a
layout), then row-bound into one table (`read_variant_input_cohort()`,
`input.R`). Sample identity is whatever ingestion already assigns per file
(the tumor genotype column name for a VCF with a `FORMAT` block — i.e. the
sample name Mutect2 wrote into the VCF header, not the file name — or the
file's basename otherwise); a `SOURCE_FILE` column is added so the originating
file stays traceable in `variants_all.tsv.gz`. `total_samples` is then the
true cohort size, and `variant_cohort_freq`/`locus_cohort_freq` reflect actual
recurrence across samples.

This is an alternative to pre-merging the VCFs yourself with `bcftools merge`
into one multi-sample VCF (which `tumoronly` also reads correctly — a `FORMAT`
block with N sample columns expands to N rows per variant the same way), not
a replacement for it: `bcftools merge` is still the right tool if you need a
merged VCF for other purposes (e.g. feeding a different tool). For `tumoronly`
alone, listing the files directly avoids the extra merge/index step and the
allele-padding edge cases `bcftools merge` can introduce.

All input files in a cohort must share one genome build; `resolve_genome_build()`
detects each file's header build independently and errors on disagreement
instead of silently picking one.

## Cohort-size-adaptive recurrence

`variant_cohort_freq`/`locus_cohort_freq` are fractions (`n_samples /
total_samples`), so their reliability as artifact evidence depends entirely on
`total_samples`. At `total_samples = 2`, one sample sharing a variant with the
other is already a 50% fraction — indistinguishable, by the fraction alone,
from a true systematic artifact recurring across a 100-sample cohort. Applying
`cohort.recurrent_variant_fraction_artifact` (default `0.30`) directly to a
2- or 3-sample cohort pushes most shared variants above threshold and into
`REVIEW`/`likely_artifact`, regardless of whether they are real artifacts.

`cohort_recurrent_flag(x, cfg, level)` (`utils.R`) is the single place that
decides whether a variant's or locus' recurrence counts as evidence, and is
now the only way `variant_cohort_freq`/`locus_cohort_freq` reach
`recurrent_germline_signal()`, `cohort_artifact_score()`,
`classify_technical_artifact_evidence()` and `classify_internal_recurrence()`
(`scoring.R`). It requires, in order:

1. **The cohort is large enough to evaluate at all** —
   `cohort.min_cohort_size_for_recurrence` (default `10`). Below it, the
   function returns `NA` (not `FALSE`): recurrence is "not assessable at this
   cohort size", not "no recurrence detected". `NA` is treated as no evidence
   everywhere it is consulted (`isTRUE_vec()` / explicit `!is.na(flag) & flag`
   checks) — it never silently defaults to a REVIEW-forcing `TRUE`.
2. **An absolute floor on how many samples actually carry it** —
   `cohort.min_recurrent_samples` / `min_recurrent_locus_samples` (default
   `3` each). A fraction alone cannot tell "shared by 1 of 2 samples" apart
   from "shared by 30 of 100"; this floor makes the count itself part of the
   decision, not just the ratio.
3. **The existing fraction-of-cohort threshold** —
   `cohort.recurrent_variant_fraction_artifact` /
   `recurrent_locus_fraction_artifact`, unchanged.

```yaml
cohort:
  min_cohort_size_for_recurrence: 10
  min_recurrent_samples: 3
  min_recurrent_locus_samples: 3
```

Practical effect: a 2-sample run never flags anything via cohort recurrence
(cohort too small — `NA`, whatever the fraction). A 100-sample run behaves
as before once a variant clears both the count floor and the fraction. Runs
in between scale smoothly as `total_samples` grows — see docs/GETTING_STARTED.md
for how to grow a cohort incrementally.

While consolidating these call sites, `cohort_artifact_score()`'s use of
`saturating_score(vf, threshold, 0.20)` was found to be inverted: with
`threshold = 0.30 > 0.20`, the ramp ran backwards — a **rare** variant
(`vf = 0.05`) scored `1` (maximum artifact evidence) and a genuinely
**recurrent** one (`vf = 0.60`) scored `0`. Every other `saturating_score()`
call in `scoring.R` ramps *up from* its threshold (e.g.
`saturating_score(dp, min_depth, min_depth * 4)`); this one now does too
(`saturating_score(vf, threshold, threshold + 0.20)`), independently of the
cohort-size fix above.

## Prioritization (heuristic, decomposed)

`PRIORITY_SCORE_BASE = 0.25·consequence + 0.20·computational + 0.25·hotspot +
0.15·driver + 0.10·clinical + 0.05·cosmic_global`, each component in [0,1] and
exported separately (`PRIORITY_*_COMPONENT`, `PRIORITY_COMPONENTS`). Weights are
**heuristic** until calibrated. Priority ranks retained variants; it never filters.
