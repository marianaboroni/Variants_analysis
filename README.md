# tumoronly — traceable tumor-only somatic variant filtering

`tumoronly` is a reproducible R package for tumor-only somatic variant analysis.
It reads a VEP-annotated VCF, validates the genome build, applies a documented and
auditable set of tumor-only filters, annotates variants against a locally-prepared
COSMIC database (lifted over **once**, at build time, by canonical genomic key),
exports a `maftools`-compatible MAF, draws oncoplots with `maftools`, and renders an
explanatory HTML report. **OncoKB annotation is a separate, optional step that never
alters filtering decisions.**

## 1. What it does

```
VEP VCF ─▶ build check ─▶ tumor-only filters ─▶ COSMIC evidence ─▶ classify
        ─▶ decision trail ─▶ tables + MAF + oncoplot + HTML report
                                            └▶ (optional) OncoKB annotation
```

Every variant keeps a full **decision trail** (`filter_status`, `filter_reasons`,
`filters_failed`, `filters_passed`, `confidence_category`, `confidence_score`) and is
retained in `variants_all` — nothing is dropped silently.

## Tutorial: how the pipeline works (end to end)

This is a walkthrough of one analysis. Commands are the single CLI `tumoronly`
(`exec/tumoronly` in the source tree, or `library(tumoronly)` if installed).

### Step 0 — check your environment

```bash
tumoronly doctor
```

Reports which R packages / external tools are present and what that enables
(cross-build COSMIC liftover, REF validation, OncoKB API, oncoplots, HTML report).
The core pipeline only needs `data.table`, `yaml`, `jsonlite`, `digest`.

### Step 1 — (once) prepare the local COSMIC database

COSMIC is **not** redistributed — you download the TSV and prepare it **once**.
Liftover (if your VCF and COSMIC differ in build), normalization, dedup, the
canonical key `BUILD|CHROM|POS|REF|ALT`, a checksummed manifest and provenance
tables are all produced here, never per sample:

```bash
tumoronly prepare-cosmic --config config/example.yml
```

This is idempotent: re-running reuses the processed DB unless an input changed
(then `--force` is required). It is optional — the pipeline runs without COSMIC,
just without COSMIC evidence.

### Step 2 — run the analysis

```bash
tumoronly run --config config/example.yml
```

Internally, for each run:

1. **Read + build check** — the VEP VCF header is parsed; `input.genome_build`
   prevails but any conflict with the header is a hard error (build is never
   guessed). The `CSQ` field order is read dynamically from the header; only the
   needed fields are extracted.
2. **Feature engineering** — depth, VAF, alt count, variant type, population AF,
   PoN flag, etc. are derived per variant.
3. **Technical filters** — adaptive depth / alt / VAF / TLOD / MBQ / MMQ + caller
   `FILTER` + sample-QC gate → `hard_filter_pass`.
4. **Scoring + classification** — evidence scores → a single conservative
   `final_class` (the one and only decision point).
5. **COSMIC evidence (tumor-type-stratified)** — variants are matched by the
   canonical key; occurrences are split into the sample's tumor type vs. others
   (see *Tumor-type-stratified COSMIC evidence* below). This affects
   **confidence only** (experimental), never `filter_status`.
6. **Decision trail** — every variant gets `filter_status` (PASS / REVIEW / FAIL),
   `filter_reasons`, `filters_failed`, `filters_passed`, `CONFIDENCE_SCORE_BASE`,
   `CONFIDENCE_CATEGORY_BASE`. No variant is dropped — excluded ones stay in
   `variants_all` with their trail.
7. **Outputs** — standardized run directory (tables, validated MAF, maftools
   oncoplot, `manifest.json`, resolved config, HTML report).

### Step 3 — (optional) OncoKB annotation, post-hoc

```bash
export ONCOKB_TOKEN="..."      # token ONLY from the environment; never stored/logged
tumoronly annotate-oncokb --run-dir results/example_run   # scope: retained | all
```

This reads the **frozen** results, adds `ONCOKB_*` columns to a separate table,
and **never** changes `filter_status`. If the API is unavailable the step still
completes (variants marked not annotated).

### Step 4 — (re-)render the report

```bash
tumoronly report --run-dir results/example_run
```

### What is expected as input

- **Primary input (`input.vcf`)**: a **VEP-annotated VCF** (`.vcf`/`.vcf.gz`) — e.g.
  Mutect2 output annotated with Ensembl VEP (a `CSQ` INFO field). A wide,
  VEP-flattened TSV is also accepted. Per-sample genotype columns (`GT:AD:DP:AF`,
  `TLOD`, …) are parsed. Useful CSQ fields: `SYMBOL`, `Consequence`, `HGVSp`,
  `gnomAD*_AF`, `CLIN_SIG`.
- **`input.genome_build`**: `GRCh37` or `GRCh38` (prevails over the header).
- **Sample metadata (optional, `input.sample_metadata`)**: TSV with `sample_id`
  and, for COSMIC tumor-context, `tumor_type`, `tumor_subtype`, `primary_site`,
  `histology`.
- **COSMIC (optional)**: a raw COSMIC TSV (`cosmic.raw_file`) with genomic
  coordinates + `PRIMARY_SITE`/`PRIMARY_HISTOLOGY`/counts, plus the versioned
  harmonization table `config/cosmic_tumor_type_mapping.tsv`.
- **Reference FASTA (optional, `reference.fasta`)**: enables REF-allele validation
  during COSMIC preparation.

### What is produced as output

```
results/<run_id>/
  manifest.json            provenance, checksums, tool/package versions, counts
  config.resolved.yml      the exact configuration used (defaults resolved)
  logs/
  tables/
    variants_all.tsv.gz         every input variant + full decision trail
    variants_retained.tsv.gz    PASS / REVIEW (kept as candidates / for review)
    variants_excluded.tsv.gz    FAIL (kept, not deleted)
    filter_audit.tsv.gz         per-filter definitions + counts + % loss
    cosmic_matches.tsv.gz       global + tumor-context COSMIC evidence per variant
    oncokb_annotations.tsv.gz   only after annotate-oncokb
  maf/filtered.maf.gz            validated with maftools::read.maf
  plots/                         maftools oncoplot + summary
  report/tumor_only_report.html  explanatory report (funnel, provenance, tables, limitations)
```

Key output columns to read first: `filter_status`, `CONFIDENCE_CATEGORY_BASE`,
`final_class`, `filter_reasons`, and the COSMIC context columns
(`COSMIC_TUMOR_CONTEXT_STATUS`, `COSMIC_MATCHING_TUMOR_OCCURRENCES`, …).

## 2. Limitations of tumor-only analysis

- Tumor-only calling **cannot definitively distinguish somatic from germline** variants.
- **Presence in COSMIC is supporting evidence** of recurrence, not proof of somaticity.
- **Absence from COSMIC does not exclude** relevance.
- **OncoKB is used only as a confirmatory annotation** and never changes filtering.
- Results require review by qualified professionals. This is **not a validated clinical test**.

## 3. Installation

```r
# from the package root
install.packages(c("data.table", "yaml", "jsonlite", "digest"))   # required
# optional but recommended:
install.packages(c("maftools", "rmarkdown", "ggplot2"))
# BiocManager::install(c("GenomicRanges", "rtracklayer", "Biostrings"))  # cross-build COSMIC
# install.packages("httr2")                                              # OncoKB API

R CMD INSTALL .
```

You can also run without installing (dev mode): the CLI at `exec/tumoronly`
auto-loads the sources.

## 4. External dependencies (optional)

| capability | needs |
|---|---|
| same-build COSMIC preparation | nothing extra |
| **cross-build** COSMIC liftover | `rtracklayer` + UCSC chain file, **or** CrossMap, **or** UCSC `liftOver` |
| COSMIC REF validation | `Biostrings` + target FASTA |
| OncoKB annotation | `httr2` + `ONCOKB_TOKEN` env var |
| oncoplots / MAF validation | `maftools` |
| HTML report | `rmarkdown` (+ pandoc) — falls back to a self-contained HTML otherwise |

Run `tumoronly doctor` to see exactly what is available.

## 5. Prepare the local COSMIC database

The package **does not redistribute COSMIC**. Download the COSMIC TSV yourself and
point the config at it. Preparation (liftover, normalization, dedup, canonical key,
manifest) runs **once**; later analyses reuse the processed DB.

```yaml
cosmic:
  raw_file: data/COSMIC_updated.tsv.gz
  release: v99
  source_build: GRCh37
  target_build: GRCh38
  cache_dir: db/cosmic
  chain_file: reference/hg19ToHg38.over.chain.gz   # only for cross-build
```

```bash
tumoronly prepare-cosmic --config config/example.yml
```

This writes `db/cosmic/<release>/<source>_to_<target>/` containing `cosmic_db.rds`,
`manifest.json` (release, checksums, builds, tool versions, input/converted/unmapped/
ref-mismatch/duplicate counts, schema) and provenance tables (`unmapped.tsv.gz`,
`ref_mismatch.tsv.gz`, `duplicates.tsv.gz`, …). Re-running reuses the DB unless inputs
changed (then `--force` is required). Matching later uses only the canonical key
`BUILD|CHROM|POS|REF|ALT`.

### Tumor-type-stratified COSMIC evidence

COSMIC evidence is interpreted **in the context of the sample's tumor type**, not as
global presence. The sample's `tumor_type`/`tumor_subtype`/`primary_site`/`histology`
are harmonized to COSMIC categories through the **explicit, versioned** table
`config/cosmic_tumor_type_mapping.tsv` (never fuzzy text matching). Each variant gets a
controlled `COSMIC_TUMOR_CONTEXT_STATUS` (`exact_match`, `compatible_match`,
`pan_cancer_recurrent`, `other_tumor_only`, `tumor_type_unknown`,
`cosmic_tumor_type_missing`, `no_cosmic_match`), matching/other occurrence counts, and
decomposed scores (`cosmic_genomic_match_score`, `cosmic_recurrence_score`,
`cosmic_tumor_specificity_score`, `cosmic_total_score`). Tumor context affects
**confidence only** — never `filter_status` automatically. Rules, weights and thresholds
are documented in [docs/COSMIC_DECISION_RULES.md](docs/COSMIC_DECISION_RULES.md).

## 6. Run one sample

```bash
tumoronly run --config config/example.yml
```

## 7. Run a cohort

Point `input.vcf` at a multi-sample VEP VCF (or a wide annotated TSV with per-sample
genotype columns). Cohort recurrence features are computed automatically.

## 8. Generate the report

The report is produced during `run`. To (re-)render:

```bash
tumoronly report --run-dir results/example_run
```

## 9. Optional OncoKB annotation

```bash
export ONCOKB_TOKEN="..."          # token ONLY via env var; never stored or logged
tumoronly annotate-oncokb --run-dir results/example_run          # scope: retained (default) | all
```

This writes `tables/oncokb_annotations.tsv.gz` and **never** modifies `filter_status`.

## 10. Outputs

```
results/<run_id>/
  manifest.json            provenance, checksums, tool/package versions, counts
  config.resolved.yml
  logs/
  tables/
    variants_all.tsv.gz         every variant + full decision trail
    variants_retained.tsv.gz    PASS / REVIEW
    variants_excluded.tsv.gz    FAIL (kept, not dropped)
    filter_audit.tsv.gz         per-filter definitions + counts
    cosmic_matches.tsv.gz
    oncokb_annotations.tsv.gz   (after annotate-oncokb)
  maf/filtered.maf.gz           validated by maftools::read.maf
  plots/                        maftools oncoplot + summary
  report/tumor_only_report.html
```

## 11. Troubleshooting

- **"Could not determine genome build"** → set `input.genome_build` (GRCh37/GRCh38).
- **"Genome-build conflict"** → the config and the VCF header disagree; fix one.
- **"no liftover backend"** → install rtracklayer (+chain) / CrossMap / liftOver, or
  prepare COSMIC on the same build as your VCF.
- **manifest mismatch on prepare-cosmic** → inputs changed; re-run with `--force`.
- **OncoKB `no_token`** → export `ONCOKB_TOKEN`; the run still completes without it.

## 12. Minimal complete example

```bash
tumoronly doctor
tumoronly prepare-cosmic  --config config/example.yml
tumoronly run             --config config/example.yml
export ONCOKB_TOKEN="..."
tumoronly annotate-oncokb --run-dir results/example_run
tumoronly report          --run-dir results/example_run
```

## Public API

`run_tumor_only()`, `prepare_cosmic_db()`, `annotate_oncokb()`, `create_maf()`,
`render_tumor_only_report()`, `check_installation()`. See `docs/ARCHITECTURE.md` and
`docs/REFACTOR_AUDIT.md` for design and rationale.

## Notes on the refactor

This package replaces an earlier TSV-first pipeline. De-scoped modules (semi-supervised
ML filtering, ancestry PCA, clonality, the bespoke HTML dashboard) were moved to
`archive/` and are **not** part of the package; see `docs/REFACTOR_AUDIT.md`. The most
important scientific change: **OncoKB was removed from the classifier** and is now a
post-hoc annotation only.
