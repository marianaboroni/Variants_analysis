# Target architecture — `tumoronly`

Phase 2 deliverable. Defines the public API, the single CLI, the directory
layout, the COSMIC and OncoKB flows, the table schemas, and the migration plan.

## 1. Package identity

- Package name: **`tumoronly`**. One executable: **`tumoronly`**.
- Formal R package: `DESCRIPTION` + `NAMESPACE` (roxygen2 source comments; a
  hand-maintained `NAMESPACE` is committed so the package loads without roxygen2
  present).

## 2. Public API (small and stable)

Everything else is internal (not exported).

```r
check_installation()                       # tumoronly doctor
prepare_cosmic_db(config, force = FALSE)   # tumoronly prepare-cosmic
run_tumor_only(config, run_id = NULL)      # tumoronly run
annotate_oncokb(run_dir, scope = "retained")  # tumoronly annotate-oncokb
render_tumor_only_report(run_dir)          # tumoronly report
create_maf(variants, path = NULL)          # used by run + reusable
```

Each public function: roxygen docs + example, argument validation, informative
errors, structured return value, deterministic behavior (fixed seeds, sorted
outputs).

## 3. CLI

One dispatcher, `inst/exec/tumoronly`, loads the package once (installed
`library(tumoronly)`, or a dev loader that `source()`s `R/` through a **single**
`tumoronly:::load_all_r()` call — never 20 `source()` lines in the executable).

```
tumoronly doctor
tumoronly prepare-cosmic  --config config.yml [--force]
tumoronly run             --config config.yml [--run-id ID]
tumoronly annotate-oncokb --run-dir results/ID [--scope retained|all]
tumoronly report          --run-dir results/ID
```

The two legacy runners become 3-line deprecation wrappers that warn and delegate
to `tumoronly run` (no logic duplicated).

## 4. Directory layout

```
DESCRIPTION  NAMESPACE  LICENSE  README.md
R/
  config.R          # read + validate config, resolve defaults
  input.R           # VCF/TSV read, dynamic CSQ, build detection, standardization
  features.R        # per-variant feature engineering (kept)
  filtering.R       # hard filters + scoring + classification + DECISION TRAIL
  cosmic.R          # idempotent COSMIC build (liftover, key, manifest, lock)
  cosmic_match.R    # canonical-key annotation of variants against processed COSMIC
  oncokb.R          # separate post-hoc OncoKB annotation (never edits filters)
  maf_reporting.R   # create_maf + maftools validation + oncoplot + tables
  report.R          # render_tumor_only_report
  workflow.R        # run_tumor_only orchestrator
  doctor.R          # check_installation
  utils.R           # shared helpers (io, checksums, logging)
inst/
  exec/tumoronly            # single CLI
  report/tumor_only_report.Rmd
config/example.yml
tests/testthat/ + tests/fixtures/
archive/            # de-scoped modules (ml, ancestry, clonality, dashboard), NOT built
docs/
```

Kept core-science files (`hard_filters.R`, `scoring.R`,
`guideline_classification.R`, `driver_classification.R`, `tmb.R`,
`qc_metrics.R`, `metadata.R`) are folded conceptually into `filtering.R`'s
pipeline but stay as separate `R/` files for reviewability; OncoKB references are
surgically removed from them.

## 5. COSMIC flow (build once, reuse forever)

Liftover and normalization happen **only** in `prepare_cosmic_db()`, never per
sample.

```
raw COSMIC tsv.gz
  → validate format + declared source build
  → [cross-build?] liftover backend (rtracklayer | CrossMap | UCSC liftOver)
  → normalize: chr naming (chr1↔1), split multiallelic, trim common bases,
               left-align indels, strand handling
  → validate REF against target FASTA
  → canonical key  BUILD|CHROM|POS|REF|ALT
  → dedup, preserve original COSMIC IDs + provenance
  → write processed db + manifest.json + provenance tables
```

**Liftover backend** is an abstraction (`cosmic_liftover_backend()`): if
`source_build == target_build`, no liftover (fully runnable here); otherwise it
requires one of rtracklayer/CrossMap/liftOver and **errors with an actionable
message if none is available** — it never silently degrades or drops variants.

**Idempotency & cache.** The manifest records: COSMIC release, raw filename +
SHA-256, source/target build, chain SHA-256, FASTA SHA-256/id, build date,
package version, tool versions, N input / converted / unmapped / REF-mismatch /
duplicate, and the output column schema. On reuse: if the processed db exists, the
manifest is intact, checksums match, and the build is compatible → reuse without
re-running liftover. If anything relevant changed → refuse and require `--force`
(never silently use an incompatible db). A directory-based **lock**
(`.build.lock`) prevents two concurrent builds.

**Storage format.** Chosen: **`data.table` keyed table serialized as compressed
`.rds` + a sidecar `.tsv.gz`**. Rationale (benchmark in the final report):
`arrow`/`duckdb`/`RSQLite` are not available in the target environment, and a
keyed `data.table` gives O(log n) canonical-key joins with zero extra
infrastructure — the simplest option that performs well at COSMIC scale. The
match layer loads the db **once** per run.

**Provenance tables** (always written, never silent loss):
`converted.tsv.gz`, `unmapped.tsv.gz`, `multi_mapped.tsv.gz`,
`ref_mismatch.tsv.gz`, `dropped_missing.tsv.gz`, `duplicates.tsv.gz`.

**Processed db location:** `cache_dir/<release>/<source>_to_<target>/`.

## 6. Matching (VEP ↔ COSMIC)

Canonical key only: `BUILD|CHROM|POS|REF|ALT` with normalized chrom (no `chr`
prefix, `MT`→`MT`) and left-normalized REF/ALT. Output columns:
`COSMIC_MATCH`, `COSMIC_RELEASE`, `COSMIC_MUTATION_IDS`,
`COSMIC_OCCURRENCE_COUNT`, `COSMIC_TUMOR_TYPES`, `COSMIC_EVIDENCE_SUMMARY`.
COSMIC may contribute to evidence/priority (documented, traceable in the decision
trail); presence is never sufficient alone, absence never excludes.

## 7. OncoKB flow (separate, post-hoc, non-filtering)

```
run_tumor_only  → results frozen (filter_status written)
                → [optional] annotate_oncokb(run_dir)
                → report enriched with ONCOKB_* columns
```

Guarantees enforced by code + tests:
- `annotate_oncokb` reads frozen `variants_all.tsv.gz`, adds `ONCOKB_*` columns,
  writes `tables/oncokb_annotations.tsv.gz`; it **never** modifies
  `filter_status`, `filter_reasons`, or promotes FAIL→PASS.
- Token only from `ONCOKB_TOKEN` env var; never logged/serialized/committed.
- Batch requests, retries w/ exponential backoff, timeouts, rate-limit handling,
  per-variant cache (`cache_dir`), query date/endpoint/status recorded.
- API unavailable → report still renders; annotation marked `not_annotated`.
- Default scope `retained`; `--scope all` for audit.

## 8. Table schemas (run outputs)

```
results/<run_id>/
  manifest.json                     # provenance, checksums, tool/pkg versions, counts
  config.resolved.yml               # config after defaults
  logs/run.log
  tables/
    variants_all.tsv.gz             # every input variant + full decision trail
    variants_retained.tsv.gz        # filter_status == PASS
    variants_excluded.tsv.gz        # filter_status == FAIL (kept, not dropped)
    filter_audit.tsv.gz             # per-filter definition + counts + loss
    cosmic_matches.tsv.gz
    oncokb_annotations.tsv.gz       # only after annotate-oncokb
  maf/filtered.maf.gz
  plots/
  report/tumor_only_report.html
```

**Decision-trail contract** on every variant (in `variants_all`):
`filter_status` (PASS/FAIL) · `filter_reasons` (delimited) · `filters_failed` ·
`filters_passed` · `confidence_category` · `confidence_score` · plus the retained
scientific columns (`final_class`, `primary_reason`, `COSMIC_*`, driver/hotspot,
guideline). No variant is dropped by failing a filter — it stays in
`variants_all` with its trail.

`filter_status` is derived deterministically from the existing engine:
`PASS` ⇔ `hard_filter_pass` AND `final_class ∈ {high_confidence_somatic,
probable_somatic, manual_review_required, uncertain_tumor_only}` is *not* the
rule — rather `filter_status = PASS` when the variant is retained as a somatic
candidate (`final_class ∈ {high_confidence_somatic, probable_somatic}`) OR flagged
for review; `FAIL` for `technical_fail`, `likely_artifact`, `likely_germline`.
The precise mapping is documented in `filter_audit.tsv.gz` and configurable.

## 9. Genome build

Resolved in `input.R`: config `input.genome_build` **prevails**, but the VCF
header (`##reference`, contig names/lengths) is parsed and any conflict raises a
critical error. If neither config nor header determines the build → stop with an
actionable message (`Set input.genome_build to GRCh37 or GRCh38`).

## 10. Migration plan

1. Package skeleton + single loader + CLI; legacy runners → deprecation wrappers.
   *(Repo stays runnable.)*
2. Surgically remove OncoKB from `scoring.R` + `guideline_classification.R` +
   `validation.R`; add regression test (`filter_status` identical ± OncoKB).
3. `input.R`: build detection/validation; keep the good dynamic-CSQ reader.
4. `cosmic.R` + `cosmic_match.R`: idempotent build, canonical key.
5. `filtering.R`: wrap existing engine, emit the decision trail; standardized run
   dir via `report.R`/`maf_reporting.R`.
6. `oncokb.R`: separate command.
7. `report.R`: Rmd (Quarto not available → R Markdown, which is present and
   mature here).
8. Tests + `R CMD check` (best-effort given missing Suggests) + fixture E2E.
9. Archive de-scoped modules; rewrite README; final report.
```
