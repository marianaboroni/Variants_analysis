# Real data validation gate for v2

Date: 2026-07-28  
Branch requested: `feature/cohort-ingestion-graded-filters`  
Branch used: `feature/cohort-ingestion-graded-filters`  
Requested remediation commit: `ed2cd31`  
Available branch HEAD at start: `9aea9cb6e1e784e04c07331f6c7a4afc4b2e30fd`

`ed2cd31` was not present in the local clone or in `origin` at validation time.
The validation was therefore performed from the available requested branch HEAD
and the missing commit is recorded as a provenance limitation.

## Decision

**ready for v2 development**

The current branch can ingest and run the required real TSV end-to-end after the
blocking reproducibility issues below were fixed and regression-tested. This is
not a declaration that v2 is complete. The real file is a one-sample WGS slice,
not an active multi-sample cohort, so v2 still needs explicit multi-sample
cohort tests and the new user-facing contract/reporting work.

## Input data

Input file:

```text
data/test/WGS_all_patients_first100k.tsv
```

Read-only inspection:

| Item | Result |
|---|---:|
| Raw rows | 100,000 |
| Columns | 87 |
| Active `Sample_Barcode` values | 1 |
| Active `Patient_ID` values | 1 |
| Sample-like genotype columns | 21 |
| Non-empty sample-like genotype columns | 1 |
| Unique variant keys (`CHROM:START:REF:ALT`) | 100,000 |
| Unique loci (`CHROM:START`) | 100,000 |
| Full-row duplicates | 0 |
| Duplicate sample-variant keys | 0 |
| Chromosomes represented | `chr1`-`chr16` |
| Position range | 10,567-248,945,365 |

Important limitation: although the file contains 21 sample-like columns, only
`ATPBR.054.001.0183.1510.01ABD_tumor_ATPBR.054.001.0183.1510.01ABD` has data.
The actual active sample in this validation is
`tumor_ATPBR-054-001-0183-1510-01ABD`.

## Columns and caller evidence

Core columns available:

```text
Sample_Barcode, CHROM, START, END, REF, ALT, FILTER, ref_count, alt_count, AF,
SYMBOL, Gene, Feature, Consequence, IMPACT, VARIANT_CLASS, MBQ, MMQ, GERMQ,
AS_SB_TABLE, GMAF, gnomADe_AF
```

No explicit caller-name column is present. The field set is consistent with a
Mutect2/FilterMutectCalls-like VCF that has been flattened after VEP annotation
(`FILTER`, `GERMQ`, `AS_FilterStatus`, `AS_SB_TABLE`, `MBQ`, `MMQ`), but the
caller cannot be proven from the TSV alone. `FILTER` is `PASS` in all 100,000
rows, and `TLOD` is absent.

## Variant representation

| Representation check | Result |
|---|---:|
| `chr` prefix present | 100,000 / 100,000 |
| REF/ALT restricted to A/C/G/T/N strings | 100,000 / 100,000 |
| SNV | 58,702 |
| insertion | 16,623 |
| deletion | 12,495 |
| equal-length substitution/DNP-like | 12,180 |
| AF outside [0, 1] | 0 |
| Negative read-count components | 0 |

## Missingness and metric availability

| Field/group | Availability |
|---|---|
| Derived depth (`ref_count + alt_count`) | 100,000 / 100,000 |
| `AF` / VAF | 100,000 / 100,000 |
| `alt_count` / `ref_count` | 100,000 / 100,000 |
| `MBQ`, `MMQ`, `GERMQ` | 100,000 / 100,000 |
| `QUAL` | 0 / 100,000 |
| `TLOD` | column absent |
| `AS_FilterStatus` | 0 / 100,000 |
| `AS_UNIQ_ALT_READ_COUNT`, `MFRL`, `MPOS`, `STR` | 0 / 100,000 |
| Strand-bias evidence | `AS_SB_TABLE` present, but no parsed strand-bias score |
| Mapping-quality evidence | `MMQ` present |
| Population frequency | `GMAF` in 14,547 rows; `gnomADe_AF` in 2,602 rows |
| Functional annotation | `Consequence` and `IMPACT` complete; `SYMBOL` missing in 54,115 rows |

Numeric distributions from read-only inspection:

| Metric | n | Min | Q10 | Q25 | Median | Mean | Q75 | Q90 | Q99 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Derived DP | 100,000 | 1 | 12 | 42 | 92 | 106.15 | 131 | 159 | 453 | 3,816 |
| `alt_count` | 100,000 | 1 | 5 | 13 | 42 | 47.63 | 66 | 95 | 176 | 3,371 |
| `AF` | 100,000 | 0.0047 | 0.123 | 0.390 | 0.515 | 0.561 | 0.857 | 0.981 | 0.993 | 1.000 |
| `GERMQ` | 100,000 | 3 | 9 | 15 | 18 | 36.49 | 66 | 93 | 93 | 93 |
| `MBQ` | 100,000 | 20 | 25 | 31 | 32 | 31.70 | 35 | 35 | 37 | 39 |
| `MMQ` | 100,000 | 21 | 38 | 42 | 60 | 53.51 | 60 | 60 | 60 | 60 |
| `GMAF` | 14,547 | 0 | 0.0038 | 0.0611 | 0.3117 | 0.3508 | 0.5731 | 0.7930 | 0.9872 | 1 |
| `gnomADe_AF` | 2,602 | 0 | 0 | 0.0000009 | 0.00393 | 0.1366 | 0.1875 | 0.4969 | 1 | 1 |

## Commands executed

Read-only inspection:

```bash
Rscript - <<'RS'
suppressPackageStartupMessages(library(data.table))
dt <- fread("data/test/WGS_all_patients_first100k.tsv", sep = "\t")
# column, missingness, duplicate, format, and distribution checks
RS
```

Initial workflow attempt:

```bash
/usr/bin/time -l Rscript exec/tumoronly run \
  --config config/wgs_all_patients_subset_config.yml \
  --run-id real_data_phase1_baseline
```

This failed before ingestion because the TSV has no build metadata and the config
did not define `input.genome_build`.

Final workflow command:

```bash
/usr/bin/time -l Rscript exec/tumoronly run \
  --config config/wgs_all_patients_subset_config.yml \
  --run-id real_data_phase1_final
```

Report rerender check:

```bash
Rscript exec/tumoronly report \
  --run-dir results/wgs_all_patients_subset_filter/real_data_phase1_final
```

Regression tests:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-cohort-recurrence-adaptive.R")'
Rscript -e 'testthat::test_file("tests/testthat/test-plots-prioritization.R")'
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

## Problems found and corrections made

### 1. Required genome build was absent from WGS configs

The real TSV does not carry VCF header metadata, so build inference is
impossible. The initial end-to-end run stopped with:

```text
Could not determine GRCh37 vs GRCh38.
Set input.genome_build in the configuration file.
```

Correction:

- added `input.genome_build: GRCh38` to
  `config/wgs_all_patients_subset_config.yml`;
- added the same explicit build to `config/wgs_all_patients_config.yml`.

This is a reproducibility fix only; it does not alter classification rules.

### 2. CLI used a stale installed package instead of branch source

The first successful run was discarded because `exec/tumoronly` loaded
`r_libs/tumoronly` before the source tree. That stale installation lacked the
cohort-size-adaptive recurrence fix and incorrectly marked all 100,000 variants
as `cohort_recurrence_artifact`.

Correction:

- `exec/tumoronly` now prefers local `R/*.R` files when executed from the source
  tree;
- installed `library(tumoronly)` is used only when no source files are available.

### 3. Single-sample recurrence was still labeled as tumor-type support

After the loader fix, artifact recurrence was no longer falsely applied, but
`classify_internal_recurrence()` labeled every row as
`recurrence_tumor_type_supported` in a one-sample file. This did not change final
classes, but it was misleading evidence.

Correction:

- tumor-type recurrence support now requires evaluable tumor-type sample counts;
- low median VAF contributes to `cohort_artifact_score()` only when the
  variant/locus is recurrent in an evaluable cohort;
- regression tests were added for both cases.

### 4. Report used stale installed template and relative plot paths

The report initially emitted `normalizePath()` and PNG fetch warnings because
`system.file()` found the installed template before the source-tree template,
and the plot manifest stored relative paths.

Correction:

- `report_template_path()` now prefers the source tree via `TUMORONLY_R_DIR`;
- `render_tumor_only_report()` normalizes `run_dir` once at entry;
- the R Markdown template resolves plot paths relative to the run directory.

### 5. maftools rainfall sidecar was written to repository root

`maftools::rainfallPlot()` wrote
`tumor_ATPBR-054-001-0183-1510-01ABD_Kataegis.tsv` to the current working
directory.

Correction:

- `save_plot_pair()` now temporarily sets the working directory to `out_dir`
  while drawing;
- a regression test verifies sidecar files stay inside the plot directory;
- the final run writes the sidecar to:

```text
results/wgs_all_patients_subset_filter/real_data_phase1_final/plots/tumor_ATPBR-054-001-0183-1510-01ABD_Kataegis.tsv
```

## Final run results

Final run directory:

```text
results/wgs_all_patients_subset_filter/real_data_phase1_final
```

Manifest:

| Field | Value |
|---|---|
| Run id | `real_data_phase1_final` |
| Input format | `tsv` |
| Genome build | `GRCh38` |
| Input SHA-256 | `a2e6f6115ee88172caefde47b121de703ac39899eee1e5901839a02b7a23c480` |
| Package version reported | `0.1.0` |
| Internal runtime | 112.64 seconds |
| Wall time from `/usr/bin/time -l` | 141.63 seconds |
| Maximum resident set size | 792,748,032 bytes |
| macOS peak memory footprint | 2,643,244,608 bytes |

Workflow counts:

| Stage | n |
|---|---:|
| Raw rows | 100,000 |
| Ingested variants | 100,000 |
| Standardized variants | 100,000 |
| Duplicate variants removed | 0 |
| Technical hard-filter PASS | 68,700 |
| Technical hard-filter FAIL | 31,300 |
| Final PASS | 51,531 |
| Final REVIEW | 4,996 |
| Final FAIL | 43,473 |

Final class distribution:

| Final class | n |
|---|---:|
| `probable_somatic` | 51,879 |
| `likely_artifact` | 18,375 |
| `technical_fail` | 12,925 |
| `likely_germline` | 12,173 |
| `uncertain_tumor_only` | 4,100 |
| `manual_review_required` | 548 |

Other invariants:

| Check | Result |
|---|---:|
| Unknown final classes | 0 |
| `somatic_score` NA | 0 |
| `somatic_score` outside [0,1] | 0 |
| Recurrence artifact failures in N=1 | 0 |
| `recurrence_category` | 100,000 `recurrence_non_informative` |
| Strong artifacts that hard-passed | 0 |
| COSMIC matches | 0 |
| Driver resources matched | 0 hotspots, 0 driver genes |

Filter audit:

| Filter | n | Percent |
|---|---:|---:|
| technical | 31,300 | 31.30 |
| population_common | 13,834 | 13.83 |
| technical_artifact | 19,425 | 19.42 |
| cohort_recurrence_artifact | 0 | 0.00 |
| final_PASS | 51,531 | 51.53 |
| final_REVIEW | 4,996 | 5.00 |
| final_FAIL | 43,473 | 43.47 |

## Outputs generated

Key outputs:

```text
results/wgs_all_patients_subset_filter/real_data_phase1_final/
├── config.resolved.yml
├── manifest.json
├── maf/filtered.maf.gz
├── plots/
│   ├── plotmaf_summary.pdf
│   ├── plotmaf_summary.png
│   ├── rainfall_tumor_ATPBR-054-001-0183-1510-01ABD.pdf
│   ├── rainfall_tumor_ATPBR-054-001-0183-1510-01ABD.png
│   ├── titv.pdf
│   ├── titv.png
│   └── tumor_ATPBR-054-001-0183-1510-01ABD_Kataegis.tsv
├── report/tumor_only_report.html
└── tables/
    ├── variants_all.tsv.gz
    ├── variants_retained.tsv.gz
    ├── variants_excluded.tsv.gz
    ├── filter_audit.tsv.gz
    ├── input_schema_report.tsv
    ├── variant_review_summary.tsv
    ├── top_prioritized_variants.tsv
    └── ...
```

The HTML report rendered without path/resource warnings after fixes. It is
self-contained, about 2.6 MB, and contains three embedded maftools images.

Figure inspection:

| Figure | Status | Notes |
|---|---|---|
| `plotmaf_summary.png` | non-empty, legible | top-gene percentages are uninformative with one sample |
| `titv.png` | non-empty | includes an empty lower panel from maftools layout; not v2 publication quality |
| `rainfall_*.png` | non-empty, legible | reflects a one-sample slice over chr1-chr16 |
| `ancestry_call_rate.png` | non-empty but weak | ancestry was not evaluable; no AIMs panel configured |

## Warnings in final run

Expected/recorded warnings:

1. `create_maf()` mapped 246 variants with unmapped VEP consequences to
   `Targeted_Region`:
   `regulatory_region_variant`, `splice_polypyrimidine_tract_variant`,
   `splice_donor_region_variant`, `splice_donor_5th_base_variant`.
2. `maftools::titv()` reported non-standard Ti/Tv class `202TRUE` three times.
   The plot was produced, but this needs a v2 plotting/QC redesign.
3. COSMIC annotation was skipped because no processed COSMIC database is
   configured.
4. Ancestry inference was `not_evaluable` because no AIMs reference panel is
   configured.

No final warning remained for `normalizePath()` or missing report resources.

## Traceable examples

These examples come from `tables/variants_all.tsv.gz` in the final run.
They are examples of algorithmic tumor-only classification, not clinical truth.

| Class | Variant | Gene | Evidence supporting | Evidence against | Missing/limited evidence |
|---|---|---|---|---|---|
| `probable_somatic` | `chr5:115212428:AG:A` | `PGGT1B` | `probable_somatic_rule`; VAF 0.603; DP 157; alt 95; `FILTER=PASS`; `somatic_score=0.9719` | no rule-level against flag | no TLOD/QUAL; no COSMIC/driver DB; no matched normal |
| `manual_review_required` | `chr12:113128127:C:T` | `RASAL1` | technical metrics pass; `somatic_score=0.8665` | `population_low_conflict`; low-frequency population evidence | no TLOD/QUAL; no COSMIC/driver DB; no matched normal |
| `uncertain_tumor_only` | `chr5:34191144:C:G` | `ENSG00000215156` | technical metrics pass; VAF 0.073; DP 593 | no decisive somatic/germline/artifact rule | no TLOD/QUAL; no COSMIC/driver DB; no matched normal |
| `likely_germline` | `chr15:96102181:T:C` | `NR2F2-AS1` | high VAF 0.833 and common population frequency (`GMAF=0.2342`) | not clinically confirmed germline | no matched normal; no clinical germline assay |
| `likely_artifact` | `chr13:113403869:G:GACCCTG...` | `ADPRHL1` | `strong_artifact_or_pon`; `low_mbq`; technical artifact flag | none strong enough to retain | no TLOD/QUAL; no orthogonal validation |
| `technical_fail` | `chr15:97487289:G:A` | `LINC00923` | low VAF hard-filter failure (`AF=0.013`) | cannot evaluate biology after technical fail | no TLOD/QUAL; no matched normal |

The current tables contain `evidence_for_somatic` and
`evidence_against_somatic`, but they do not yet contain a complete structured
`missing_evidence` column for every row. That is a required v2 product feature.

## Test evidence

Targeted tests:

```text
test-cohort-recurrence-adaptive.R: 19 pass, 0 fail
test-plots-prioritization.R: 11 pass, 0 fail
```

Full dev test suite:

```text
0 failures
1 skipped test
7 warnings, all from COSMIC fixture runs without reference FASTA validation
```

The warnings are pre-existing fixture/resource warnings and did not come from
the real-data run.

## Limitations

- This file validates a real WGS-like TSV input, but not active multi-sample
  cohort behavior: only one sample has non-empty data.
- No explicit caller-name column is present; caller identity is inferred only
  from Mutect2-like fields and upstream context.
- `TLOD`, `QUAL`, `MPOS`, `STR`, and `AS_UNIQ_ALT_READ_COUNT` are absent, so
  caller-quality evidence is incomplete.
- `AS_SB_TABLE` is present but not parsed into a strand-bias score by the
  current pipeline.
- COSMIC and driver resources are not configured; therefore driver detection and
  COSMIC contextual evidence were not validated in this run.
- Ancestry inference is not evaluable without an AIMs reference panel.
- Current maftools figures are useful as smoke-test outputs, not as v2
  publication-quality figures.
- Tumor-only classifications remain probabilistic candidates. They must not be
  presented as clinically confirmed somatic or germline calls without matched
  normal or orthogonal validation.

## Gate conclusion

The real TSV passes the Phase 1 gate after the corrections above:

- ingestion works on the real file;
- schema/reporting exposes missing core evidence;
- chromosome/position/REF/ALT normalization succeeds;
- no duplicates are introduced or hidden;
- sample grouping correctly identifies one active sample;
- fail-safe and classification rules run without unknown classes;
- `somatic_score` is bounded and complete;
- tables, MAF, plots, manifest, and HTML report are generated;
- stale installed code/template issues are fixed;
- recurrence evidence no longer misrepresents a one-sample file as a cohort.

Proceed to v2 development, with the explicit constraint that v2 must add:

- structured missing-evidence columns for every classification;
- a real multi-sample cohort regression fixture;
- a formal input validation command;
- publication-quality figures and figure data;
- a manifest/schema/report design that does not depend on maftools side effects.
