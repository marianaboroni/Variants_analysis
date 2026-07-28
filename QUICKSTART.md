# Quick Start

This guide runs the v2 workflow on the included demo data.

Tumor-only classifications are probabilistic candidates. They are not clinical
confirmation of somatic or germline status without matched-normal or orthogonal
validation.

## 1. Install From a Checkout

```bash
Rscript -e 'install.packages(c("data.table","yaml","jsonlite","digest"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Recommended for report and figures:

```bash
Rscript -e 'install.packages(c("ggplot2","rmarkdown","knitr"), repos="https://cloud.r-project.org")'
```

## 2. Validate

```bash
Rscript exec/tumoronly validate \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2_validation.html
```

Validation must return `PASS` before a strict run. Warnings are expected for
small demo data when evidence such as cohort recurrence or COSMIC is not
available.

## 3. Run

```bash
Rscript exec/tumoronly run \
  --input data/demo/variants.tsv \
  --config config/demo_v2.yml \
  --output results/demo_v2 \
  --run-id quickstart_demo
```

The same run from R:

```r
library(tumoronly)

result <- run_tumoronly(
  input = "data/demo/variants.tsv",
  config = "config/demo_v2.yml",
  output_dir = "results/demo_v2",
  run_id = "quickstart_api_demo"
)
```

## 4. Read the Outputs

Open:

```text
results/demo_v2/quickstart_demo/report.html
```

Start with:

- `report.html`
- `run_manifest.json`
- `logs/warnings.tsv`
- `tables/classified_variants.tsv`
- `tables/sample_summary.tsv`
- `tables/filter_audit.tsv`
- `tables/module_status.tsv`
- `tables/tmb_summary.tsv`
- `tables/clonality_summary.tsv`
- `tables/ml_status.tsv`

For each variant, inspect:

- `final_class`
- `somatic_score`
- `tmb_countable`
- `clonality_class`
- `ml_status`
- `filter_status`
- `evidence_supporting_classification`
- `evidence_against_classification`
- `missing_evidence`
- `classification_explanation`

Continue with [TUTORIAL.md](TUTORIAL.md) for a full walkthrough, including an
installed-package example.
