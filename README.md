# tumoronly

Traceable tumor-only variant filtering, classification, and reporting in R.

`tumoronly` ingests VCF, MAF, TSV, TXT, and gzipped variant files; validates the
input schema; normalizes variant coordinates; applies documented fail-safe
technical filters; scores tumor-only evidence; assigns conservative
classification candidates; and writes tables, figures, logs, manifests, and an
HTML report.

Tumor-only results are probabilistic. The package never treats a variant as
clinically confirmed somatic or germline without matched-normal or orthogonal
validation.

## Current v2 Scope

The v2 user workflow is:

```text
tumoronly init -> tumoronly validate -> tumoronly run -> tumoronly report
```

From a source checkout, use the CLI script directly:

```bash
Rscript exec/tumoronly init --output config.yaml

Rscript exec/tumoronly validate \
  --input variants.tsv \
  --config config.yaml \
  --output validation.html

Rscript exec/tumoronly run \
  --input variants.tsv \
  --config config.yaml \
  --output results \
  --run-id my_run

Rscript exec/tumoronly report \
  --run-dir results/my_run \
  --output results/my_run/report.html
```

After installation, locate the same CLI through R:

```bash
TO_CLI="$(Rscript -e 'cat(system.file("exec", "tumoronly", package = "tumoronly"))')"
Rscript "$TO_CLI" validate --config config.yaml
Rscript "$TO_CLI" run --config config.yaml --run-id my_run
```

The same analysis can be launched from R:

```r
library(tumoronly)

result <- run_tumoronly(
  input = "variants.tsv",
  config = "config.yaml",
  output_dir = "results",
  run_id = "my_run"
)
```

The CLI and R API call the same internal workflow. Scientific rules are not
duplicated in the command-line layer.

## Install

From GitHub:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("marianaboroni/Variants_analysis")
library(tumoronly)
```

From a local checkout:

```bash
git clone https://github.com/marianaboroni/Variants_analysis.git
cd Variants_analysis
R CMD INSTALL .
Rscript -e 'library(tumoronly); check_installation()'
```

Core dependencies are `data.table`, `yaml`, `jsonlite`, and `digest`.
`ggplot2`, `rmarkdown`, and `knitr` are recommended for richer figures and
report rendering. The core workflow does not require credentials.

Conda and Docker recipes are included:

```bash
conda env create -f environment.yml
conda activate tumoronly
R CMD INSTALL .
```

```bash
docker build -t tumoronly:dev .
docker run --rm -v "$PWD":/work -w /work tumoronly:dev R CMD INSTALL .
```

## Input Contract

At minimum, every run needs genomic coordinates and genome build:

| Required | Meaning |
|---|---|
| `CHROM` | chromosome, such as `chr17`, `17`, `X`, or `MT` |
| `POS` | 1-based genomic position |
| `REF` | reference allele |
| `ALT` | alternate allele |
| `input.genome_build` | `GRCh37` or `GRCh38` unless safely resolved from VCF header |

Strongly recommended evidence includes sample ID, depth, alternate-read count,
VAF, caller `FILTER`, population frequency, gene, consequence, TLOD, mapping
quality, base quality, and strand-bias metrics.

Common column names are autodetected. Custom layouts should use explicit YAML
mapping:

```yaml
input:
  path: variants.tsv
  genome_build: GRCh38
  column_map:
    chrom: chromosome
    pos: start
    ref: reference
    alt: alternate
    sample_id: tumor_sample
```

Validation is strict by default. Permissive mode records missing evidence but
does not invent or silently alter scientific values.

## Minimal Reproducible Run

The repository includes a small v2 demo config:

```bash
Rscript exec/tumoronly validate \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2_validation.html

Rscript exec/tumoronly run \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2 \
  --run-id demo_v2
```

The installed package also includes a tiny example dataset under `inst/extdata`.
See [TUTORIAL.md](TUTORIAL.md) for a fully reproducible installed-package
example that writes a config, validates, runs, and opens the same output
structure.

## Output Contract

Every v2 run writes a predictable directory:

```text
results/<run_id>/
  report.html
  tables/
    all_variants.tsv
    classified_variants.tsv
    high_confidence_somatic.tsv
    likely_germline.tsv
    likely_artifact.tsv
    known_drivers.tsv
    sample_summary.tsv
    filter_audit.tsv
  figures/
    figure_01_filtering_workflow.pdf/.svg/.png
    figure_02_qc_overview.pdf/.svg/.png
    figure_03_classification.pdf/.svg/.png
    figure_04_vaf_depth.pdf/.svg/.png
    figure_05_sample_class_distribution.pdf/.svg/.png
    figure_06_missing_evidence.pdf/.svg/.png
  figure_data/
  logs/
    warnings.tsv
  config_used.yaml
  config.resolved.yml
  session_info.txt
  run_manifest.json
```

For each classified variant, start with:

| Column | Use |
|---|---|
| `final_class` | conservative tumor-only classification candidate |
| `filter_status` | `PASS`, `REVIEW`, or `FAIL` |
| `somatic_score` | bounded support score for somatic-candidate behavior |
| `evidence_supporting_classification` | evidence favoring the assigned class |
| `evidence_against_classification` | evidence arguing against the assigned class |
| `missing_evidence` | evidence unavailable for interpretation |
| `classification_explanation` | plain-language trace of the decision |

Main classes are `high_confidence_somatic`, `probable_somatic`,
`likely_germline`, `likely_artifact`, `technical_fail`,
`manual_review_required`, and `uncertain_tumor_only`.

## Real-Data Validation Gate

The v2 work was gated on the real file
`data/test/WGS_all_patients_first100k.tsv`, not only on demo data. The gate
validated ingestion, normalization, deduplication, filtering, scoring,
classification, tables, figures, report generation, runtime, memory, and
warnings.

The validation decision is recorded in
[docs/REAL_DATA_VALIDATION.md](docs/REAL_DATA_VALIDATION.md). The file contains
the commands run, row counts, warnings, memory/runtime, class distribution,
traceable examples, and limitations.

Important limitation: the real 100k-row slice contains one active sample, so it
validates real-file ingestion and the full workflow, but not real multi-sample
cohort recurrence behavior.

## What v2 Does Not Claim

- It does not confirm clinical somatic or germline status.
- It does not reintroduce archived TMB, clonality, ancestry, or ML workflows as
  v2 core features.
- It does not use absence from external databases as proof that a variant is
  somatic.
- It does not treat generated files alone as validation; commands and contents
  must be inspected.

## Documentation

- [QUICKSTART.md](QUICKSTART.md)
- [TUTORIAL.md](TUTORIAL.md)
- [INSTALLATION.md](INSTALLATION.md)
- [INPUT_FORMAT.md](INPUT_FORMAT.md)
- [CONFIGURATION.md](CONFIGURATION.md)
- [HOW_IT_WORKS.md](HOW_IT_WORKS.md)
- [INTERPRETING_RESULTS.md](INTERPRETING_RESULTS.md)
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md)
- [docs/V2_DESIGN.md](docs/V2_DESIGN.md)
- [docs/REAL_DATA_VALIDATION.md](docs/REAL_DATA_VALIDATION.md)
- [docs/V2_FINAL_REPORT.md](docs/V2_FINAL_REPORT.md)
