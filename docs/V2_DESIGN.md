# tumoronly v2 design

Status: design for implementation after the real-data gate in
`docs/REAL_DATA_VALIDATION.md`.

## Gate dependency

v2 development is allowed only because the required real TSV
`data/test/WGS_all_patients_first100k.tsv` was run end-to-end and the Phase 1
gate concluded:

```text
ready for v2 development
```

This design does not declare v2 complete. It defines the product contract that
implementation and final verification must satisfy.

## Product goal

`tumoronly` v2 is a reproducible tumor-only SNV/indel filtering, classification,
and reporting tool for exploratory cancer genomics. It must help users run an
analysis without hiding the assumptions that make tumor-only interpretation
uncertain.

The tool must never state that a variant is clinically confirmed somatic or
germline without matched normal or orthogonal validation. User-facing language
must use terms such as:

- `high-confidence somatic candidate`;
- `probable somatic candidate`;
- `consistent with a likely germline event`;
- `flagged as a likely technical artifact`;
- `requires matched-normal or orthogonal confirmation`.

## Scientific invariants that v2 must preserve

These rules are locked unless a new real-data validation explicitly justifies a
change:

- final classes remain in the audited vocabulary:
  `high_confidence_somatic`, `probable_somatic`,
  `manual_review_required`, `uncertain_tumor_only`, `likely_germline`,
  `likely_artifact`, `technical_fail`;
- tumor-only classes are probabilistic labels, not clinical calls;
- OncoKB annotation remains post-hoc and must not alter filtering decisions;
- COSMIC global recurrence is supporting context, not proof of somaticity;
- tumor-context COSMIC evidence remains separated from the base confidence
  score unless calibrated;
- caller `FILTER` evidence is graded/configurable and must not become a
  hidden binary veto for all non-`PASS` values;
- cohort recurrence is not evaluable below configured cohort/sample floors;
- no TMB module is reintroduced;
- archived TMB/out-of-scope modules stay archived.

## User workflow

The public workflow is intentionally four commands:

```bash
tumoronly init --output config.yaml

tumoronly validate \
  --input variants.tsv \
  --config config.yaml

tumoronly run \
  --input variants.tsv \
  --config config.yaml \
  --output results/

tumoronly report \
  --run-dir results/<run_id>
```

Equivalent R API:

```r
library(tumoronly)

result <- run_tumoronly(
  input = "variants.tsv",
  config = "config.yaml",
  output_dir = "results"
)
```

The CLI must be a thin wrapper. Scientific logic must live in package functions
and be shared by CLI and R API.

## CLI command contracts

### `tumoronly init`

Purpose: create a documented YAML template with safe defaults.

Required behavior:

- writes a config file and refuses to overwrite unless `--force`;
- includes every threshold used by the scientific pipeline;
- includes comments or companion documentation explaining why each group exists;
- defaults to `strict: true` for validation;
- never guesses genome build.

Minimal options:

```text
--output config.yaml
--force
```

### `tumoronly validate`

Purpose: validate input and config before running the scientific workflow.

Required behavior:

- detects format (`vcf`, `maf`, `tsv`) by content where possible;
- maps columns by explicit YAML first, aliases second;
- reports required, optional, detected, missing, and derived fields;
- computes sample counts, variant counts, duplicate counts, coordinate format,
  REF/ALT validity, missingness, and metric distributions;
- distinguishes blockers from warnings;
- writes a validation report when `--output` is supplied;
- exits non-zero on blockers in strict mode;
- in permissive mode, warns but must not silently change values.

Minimal options:

```text
--input variants.tsv
--config config.yaml
--output input_validation.html|tsv|json
--strict
--permissive
```

### `tumoronly run`

Purpose: execute the complete workflow after validation.

Required behavior:

- accepts `--input`, `--config`, and `--output` overrides without editing the
  original config;
- runs validation first unless `--skip-validation` is explicitly provided;
- blocks on validation errors in strict mode;
- writes the resolved config into the run directory;
- writes a manifest with tool version, git commit, input checksums, parameters,
  external database versions, environment, runtime, warnings, and status;
- returns a non-zero exit code if required outputs are missing or empty.

Minimal options:

```text
--input variants.tsv
--config config.yaml
--output results/
--run-id ID
--strict
--permissive
--skip-validation
```

### `tumoronly report`

Purpose: render or rerender the HTML report from an existing run directory.

Required behavior:

- uses only files inside the run directory;
- fails clearly if required tables are missing;
- embeds figures without relying on the caller's current working directory;
- records missing or skipped figures as report warnings.

Minimal options:

```text
--run-dir results/<run_id>
--output report.html
```

## R API contracts

### `run_tumoronly()`

Public convenience wrapper around the existing scientific orchestrator.

Signature:

```r
run_tumoronly <- function(input, config, output_dir = NULL, run_id = NULL,
                          strict = TRUE, validate = TRUE)
```

Required behavior:

- reads the config;
- applies input/output/run-id overrides in memory;
- optionally runs `validate_tumoronly_input()`;
- calls `run_tumor_only()` for the scientific workflow;
- returns a structured list with `run_dir`, `manifest`, `validation`, and
  `report`.

### `validate_tumoronly_input()`

Signature:

```r
validate_tumoronly_input <- function(input, config = NULL, strict = TRUE,
                                     output = NULL)
```

Required behavior:

- wraps the existing ingestion (`read_variant_input()`);
- does not duplicate classification/scoring logic;
- returns a structured validation object;
- can write TSV/JSON/HTML summaries;
- identifies blockers with machine-readable codes.

### `write_tumoronly_config_template()`

Signature:

```r
write_tumoronly_config_template <- function(output, force = FALSE)
```

Required behavior:

- writes the canonical v2 YAML template;
- uses centralized defaults from package data or a single R helper;
- fails if `output` exists and `force = FALSE`.

## Input contract

v2 accepts VCF, MAF, and delimited TSV, but every input must map into the same
canonical fields.

Required canonical fields:

| Field | Meaning |
|---|---|
| `CHROM` | chromosome |
| `POS` | 1-based position |
| `REF` | reference allele |
| `ALT` | alternate allele |
| `genome_build` | `GRCh37` or `GRCh38`; may come from VCF header or config |

Recommended fields:

| Field | Meaning |
|---|---|
| `SAMPLE_ID` | sample identifier |
| `FILTER` | caller filter status |
| `DP` or `AD_REF` + `AD_ALT` | depth evidence |
| `VAF` or derivable VAF | allele fraction |
| `GENE` | gene symbol |
| `CONSEQUENCE` | functional consequence |
| population AF columns | gnomAD, ABraOM/SABE, local population frequency |
| caller quality metrics | e.g. `TLOD`, `MBQ`, `MMQ`, strand/mapping bias |

Supported chromosome conventions:

- `1`, `2`, ..., `22`, `X`, `Y`, `MT`;
- `chr1`, `chr2`, ..., `chr22`, `chrX`, `chrY`, `chrM`;
- mixed conventions should be reported as warnings and normalized only in the
  internal canonical key, not silently rewritten in the raw preserved columns.

Supported variants:

- SNVs;
- insertions;
- deletions;
- equal-length substitutions/DNP-like records;
- multiallelic ALT values must be split explicitly with traceability.

Required validation blockers:

- missing `CHROM`, `POS`, `REF`, or `ALT`;
- missing/unsupported genome build;
- invalid coordinates;
- invalid REF/ALT symbols after allowed normalization;
- no usable sample identifier when cohort/sample outputs are requested;
- no usable depth/VAF/alt-read evidence in strict mode;
- mixed genome builds across cohort files;
- empty input after parsing.

Warnings:

- absent caller-specific metrics such as `TLOD`;
- absent population frequency;
- absent strand-bias or mapping-quality metrics;
- absent functional annotation;
- duplicate variants;
- single-sample input when cohort recurrence is enabled;
- sample-like wide columns that are empty or inconsistent with `SAMPLE_ID`.

## Evidence model

Every classified row must carry:

| Column | Meaning |
|---|---|
| `evidence_supporting_classification` | signals supporting the assigned class |
| `evidence_against_classification` | conflicting or contrary evidence |
| `missing_evidence` | absent evidence that limits interpretation |
| `classification_explanation` | short human-readable explanation |

These columns are additive evidence reporting. They must not alter audited
classification rules unless a separate scientific change is made and validated.

Examples of `missing_evidence` values:

- `missing_matched_normal`;
- `missing_tlod`;
- `missing_population_af`;
- `missing_cosmic_database`;
- `missing_driver_resources`;
- `missing_strand_bias_metric`;
- `cohort_recurrence_not_evaluable_single_sample`.

## Output contract

The v2 run directory must be predictable:

```text
results/
├── report.html
├── tables/
│   ├── all_variants.tsv
│   ├── classified_variants.tsv
│   ├── high_confidence_somatic.tsv
│   ├── likely_germline.tsv
│   ├── likely_artifact.tsv
│   ├── known_drivers.tsv
│   ├── sample_summary.tsv
│   └── filter_audit.tsv
├── figures/
├── figure_data/
├── logs/
├── config_used.yaml
├── session_info.txt
└── run_manifest.json
```

Compatibility note: existing v1 paths (`tables/variants_all.tsv.gz`,
`manifest.json`, `report/tumor_only_report.html`) may be retained during
transition, but v2 names must be produced as stable user-facing outputs.

## Report contract

The HTML report must contain these sections:

1. Executive summary
2. Input data and cohort
3. Input validation
4. Quality control
5. Filtering workflow
6. Classification methodology
7. Cohort overview
8. Sample-level results
9. High-confidence somatic candidates
10. Likely germline variants
11. Likely artifacts
12. Known or putative drivers
13. Evidence supporting each classification
14. Warnings and missing evidence
15. Limitations
16. Reproducibility information

The report must explicitly distinguish:

- robust outputs;
- uncertain exploratory candidates;
- missing evidence;
- clinical non-readiness.

## Figure contract

v2 figures must be package functions, not one-off report code. Each function
must write:

- PDF;
- SVG;
- PNG at publication resolution;
- source data in `figure_data/*.tsv`;
- a manifest row with status and skip reason.

No empty plots are acceptable. If a plot is not applicable, v2 records
`status=skipped` and a reason.

## Implementation sequence

1. Add config defaults and template writer.
2. Add input validation object and report writer.
3. Add CLI `init`, `validate`, `run --input --output`, and `report --output`.
4. Add `run_tumoronly()` API wrapper.
5. Add v2 output aliases and manifest expansion.
6. Add structured evidence/missing-evidence columns.
7. Replace maftools smoke-test plots with v2 scientific figures.
8. Rewrite documentation and tutorials.
9. Add installation, Conda, Docker, CI.
10. Run final clean install, CLI/API equivalence, WGS real data, tests, and
    `R CMD check`.

## Acceptance boundary

Do not tag or describe v2 as complete until:

- the CLI and R API produce equivalent classifications;
- `data/test/WGS_all_patients_first100k.tsv` runs from a clean install;
- every output in the v2 contract exists and has validated content;
- every figure is visually inspected;
- every classification has supporting, contrary, and missing evidence;
- the test suite, `R CMD check`, documentation checks, and clean-install smoke
  tests pass.
