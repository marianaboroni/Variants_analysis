# Refactor audit — tumor-only variant analysis

Date: 2026-07-07. Auditor: senior review of the full repository at commit `08b3ed5`.

This document is the Phase 1 deliverable: an inventory of what exists, what is
duplicated, what is dead or unvalidated, what the *real* execution flow is, the
architectural problems, and the scientific risks. It is descriptive; the target
design is in [ARCHITECTURE.md](ARCHITECTURE.md).

## 1. Current structure

The project is **not a formal R package** (no `DESCRIPTION`/`NAMESPACE`). It is a
folder of scripts plus a `R/` directory of loose files that are loaded with
`source()`.

```
R/                     19 files, ~4400 LOC   (the actual logic)
scripts/               3 runner scripts
config/                5 YAML configs
data/demo/             7 small TSV fixtures
data/test/             1 large TSV (100k rows, WGS cohort)
r_libs/                vendored R libraries (should NOT be in git)
references/            empty (oncokb_cosmic subdir, no files)
results/               committed run outputs (should NOT be in git)
figures/               pipeline diagram
vignettes/             2 Rmd + rendered HTML
```

### Input reality vs. the stated goal
The stated goal is a **VEP-annotated VCF** pipeline. The code is in fact
**TSV-first**: `read_variants()` reads a wide TSV where VEP fields are already
flattened into columns (see `data/test/WGS_all_patients_first100k.tsv`: `SYMBOL`,
`HGVSp`, `Consequence`, `gnomADe_AF`, per-sample genotype columns, …). A VCF
reader path exists (`read_annotated_vcf`, io.R:18) and already parses the `CSQ`
header dynamically (io.R:45) and extracts a *fixed subset* of CSQ fields — this
part is actually good and worth keeping. But VCF is a secondary, lightly-tested
path today.

## 2. The three scripts (the requested focus)

| script | lines | what it does |
|---|---|---|
| `run_automated_variant_analysis.R` | 32 | 20× `source("R/*.R")`, then `run_variant_analysis(config)` |
| `run_tumor_only_filter.R` | 32 | **byte-for-byte identical logic** — same 20 `source()`s, same call |
| `build_oncokb_cosmic_reference.R` | 56 | `source()`s 4 files, then standardizes OncoKB + COSMIC + driver genes + hotspots in one script |

**Finding 1 — the first two runners are exact duplicates.** They differ only in
filename and usage string. Both `source()` the same 20 files and call the same
`run_variant_analysis()`. Pure duplication → collapse to one CLI.

**Finding 2 — 20 `source()` calls per executable.** The executables hand-load the
package file-by-file. Order-sensitive, fragile, and the reason a formal package
is needed.

**Finding 3 — `build_oncokb_cosmic_reference.R` mixes four unrelated concerns**
(COSMIC prep, OncoKB prep, driver-gene prep, hotspot prep) in one script. COSMIC
and OncoKB must be *decoupled*; the reference-prep responsibilities should be
split by domain.

## 3. Inventory of `R/` modules and how they connect

`run_variant_analysis()` (workflow.R) is the single orchestrator. Traced call graph:

| module | role | feeds the filter decision? | verdict |
|---|---|---|---|
| `io.R` | config, VCF/TSV reader (dynamic CSQ), standardization, canonical `variant_id` | infra | **keep / split** into config.R + input.R |
| `metadata.R` | merge sample metadata (tumor_type, assay, purity) | yes (tumor_type) | **keep** (fold in) |
| `features.R` | derive dp/vaf/alt_count/is_indel/max_pop_af/pon_flag/… + cohort recurrence | **yes — produces every column the classifier reads** | **keep** (core) |
| `hard_filters.R` | technical + adaptive-depth + sample-QC filters → `hard_filter_pass`, `hard_filter_reason` | **yes — the technical gate** | **keep** (core) |
| `scoring.R` | evidence scores + `classify_variants()` → `final_class`, `primary_reason` | **yes — THE decision** | **keep** (core) |
| `validation.R` | OncoKB + COSMIC matching by 6 keys (coord AND gene/protein AND sample/tumor) | **yes (both)** | **rewrite** (see risks) |
| `reference.R` | `standardize_{oncokb,cosmic}_reference` | build-time | **keep / split** |
| `driver_classification.R` | hotspot/driver-gene annotation, `driver_class` | reads final_class; sets `hotspot_match` used by scoring | **keep** (core, trim) |
| `guideline_classification.R` | somatic oncogenicity + germline ACMG points → `somatic_oncogenicity_class`, `germline_disposition` | **yes — consumed by classify_variants** | **keep** (core, fix OncoKB) |
| `confidence_ranking.R` | confidence scores/tiers/bucket/rank | additive; **consumes ML prob** | **replace** with ML-free confidence |
| `ml_filtering.R` (663 LOC) | semi-supervised GLM, pseudo-labels, active learning | additive — never writes `final_class` | **remove/archive** (unvalidated ML) |
| `ancestry.R` (676 LOC) | SNP-based PCA ancestry assignment | additive analysis | **remove/archive** (out of scope) |
| `clonality.R` | CCF / clonal-subclonal / kmeans clusters | additive analysis | **remove/archive** (out of scope) |
| `tmb.R` | per-sample TMB | additive summary | **keep** (internal, cheap, report-useful) |
| `qc_metrics.R` | per-sample QC counts | additive summary | **keep** (internal) |
| `reporting.R` | ~25 TSV writers + markdown filter report | output | **rewrite** to standardized run dir |
| `visualization.R` (591 LOC) | MAF export + maftools plots + **ggplot2 fallback oncoplot** + many ggplot figures | output | **keep MAF + maftools; drop hand-rolled oncoplot** |
| `dashboard.R` | hand-built HTML dashboard | output | **remove/archive** (replaced by report) |

### Confirmed facts (verified by reading the code)
- **`final_class` is written in exactly one place:** `classify_variants()` (scoring.R:125). Every other module only reads it. Good — one decision point.
- **`ml_filtering.R` and `confidence_ranking.R` never mutate `final_class`.** ML is purely additive; safe to remove.
- **`visualization.R` already produces `somatic_maftools_input.maf` and already calls `maftools::read.maf()` + `maftools::oncoplot()`** (visualization.R:177,186), guarded by `requireNamespace`, with a ggplot2 fallback. The maftools path is the one to keep.

## 4. Dependencies

Available in this environment: `data.table 1.14.2`, `yaml`, `jsonlite`, `digest`,
`GenomicRanges 1.46.1`, `maftools 2.10.5`, `rmarkdown 2.31`, `ggplot2 3.3.6`,
`testthat 3.1.2`.

**Missing** (blocks full end-to-end verification here): `VariantAnnotation`,
`rtracklayer`, `Biostrings`, `BSgenome`, `arrow`, `duckdb`, `RSQLite`, `httr2`,
`roxygen2`, `quarto`, `filelock`. External tools `bcftools`, `tabix`, `CrossMap`,
`liftOver`, `samtools` are all **absent**.

Consequence: **real cross-build COSMIC liftover and the live OncoKB API cannot be
executed in this environment.** The architecture treats them as pluggable
backends behind capability checks (`tumoronly doctor`), and the same-build path,
canonical-key matching, MAF/oncoplot and the OncoKB-independence guarantee are all
fully runnable and tested here.

## 5. Architectural problems

1. **Two identical runners + 20 `source()` calls** — no package boundary.
2. **No public API surface.** Everything is a loose function; nothing is
   exported/documented; no argument validation on entry points.
3. **COSMIC and OncoKB are entangled** in one build script and in the
   classifier.
4. **COSMIC/OncoKB matching uses gene+protein and sample-scoped keys**, not a
   normalized coordinate key. This is fragile and, per the brief, wrong: matching
   must be on normalized `BUILD|CHROM|POS|REF|ALT`.
5. **No genome-build detection or validation.** `chrom` is coalesced but build is
   never checked; a GRCh37 VCF against a GRCh38 COSMIC would silently mismatch.
6. **No idempotent, cached, checksummed COSMIC build.** Reference standardization
   re-runs every time; no manifest, no lock, no provenance.
7. **Outputs sprawl** — ~25 TSVs + PDFs + HTML written flat into one dir, plus
   committed `results/` in git.
8. **`r_libs/` and `results/` are committed** (though `.gitignore` now lists them).
9. **Unvalidated ML in the critical path folder** — 663 LOC of GLM/active-learning
   that the brief explicitly says not to keep unless validated and connected.

## 6. Scientific risks (highest priority)

- **RISK-1 (blocker): OncoKB changes filtering decisions.** OncoKB feeds
  `final_class` through two paths:
  - `validation.R` sets `oncokb_match` → `scoring.R::classify_oncogenic_evidence`
    (`gene_only <- oncokb_match & !hotspot_match`, scoring.R:290) → `oncogenic_category`.
  - `guideline_classification.R:52` grants oncogenicity point **OS5** when
    `oncokb_match & oncogenic`, raising `somatic_oncogenicity_class`, which
    `classify_variants` turns into `high_confidence_somatic` (scoring.R:175–186).

  The brief requires OncoKB to be a **post-hoc annotation that never alters
  filter status**. This must be removed from the classifier and moved to a
  separate `annotate-oncokb` step.
- **RISK-2: matching by gene/protein and by sample.** Sample-scoped and
  gene+protein keys can create spurious or sample-specific "recurrence" evidence.
  COSMIC evidence must be coordinate-based and build-aware.
- **RISK-3: silent build mismatch.** No guard prevents matching a GRCh37 VCF to a
  GRCh38 COSMIC (or vice-versa).
- **RISK-4: COSMIC as binary presence.** Occurrence count / tumor context are
  partially captured but the classifier reduces some paths to `cosmic_match`
  boolean; the brief wants richer, traceable COSMIC evidence fields.
- **RISK-5: no per-variant decision trail.** `final_class` + `primary_reason`
  exist, but there is no unified `filter_status` / `filters_failed` /
  `filters_passed` / `confidence_*` contract, and excluded variants are split into
  separate files rather than retained with a full audit trail.

COSMIC participating in classification is **allowed** by the brief (documented,
traceable) and is retained; OncoKB participation is **removed**.

## 7. What to keep, rewrite, remove

- **Keep (core science, unchanged thresholds):** feature engineering
  (features.R), technical/QC filters (hard_filters.R), scoring/classification
  (scoring.R), guideline classification (guideline_classification.R — minus
  OncoKB), driver/hotspot annotation (driver_classification.R), TMB & QC
  summaries.
- **Rewrite:** reference/COSMIC handling → idempotent build with liftover,
  canonical key, manifest, provenance (`cosmic.R`); OncoKB → separate step
  (`oncokb.R`); outputs → standardized run dir (`maf_reporting.R` + report);
  confidence → ML-free.
- **Remove/archive:** `ml_filtering.R`, `ancestry.R`, `clonality.R`,
  `dashboard.R`, one duplicate runner, the ggplot2 fallback oncoplot.

## 8. Regression safety

`final_class` has a single writer, so a regression test comparing `final_class`
on the demo fixture before/after the refactor is well-defined. The key behavioral
change is intentional: **removing OncoKB from classification.** That change is
isolated and testable (`filter_status` identical with and without the OncoKB
annotation step) and will be documented as an intentional scientific correction,
not a regression.
