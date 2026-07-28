# tumoronly quick start

This guide runs the v2 workflow on a TSV, VCF, or MAF input.

Tumor-only classifications are probabilistic. A result can be classified as a
high-confidence somatic candidate, likely germline-like, or likely artifact-like,
but it is not clinical confirmation without matched normal or orthogonal
validation.

## 1. Install

From the repository root:

```bash
Rscript -e 'install.packages(c("data.table","yaml","jsonlite","digest"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Optional but recommended for reports and figures:

```bash
Rscript -e 'install.packages(c("ggplot2","rmarkdown","knitr"), repos="https://cloud.r-project.org")'
```

If working from a checkout, the source-tree CLI also works without installation:

```bash
Rscript exec/tumoronly doctor
```

## 2. Create a config

```bash
Rscript exec/tumoronly init --output config.yaml
```

Edit at least:

```yaml
input:
  path: variants.tsv
  genome_build: GRCh38

analysis:
  output_dir: results
```

## 3. Validate the input

```bash
Rscript exec/tumoronly validate \
  --input variants.tsv \
  --config config.yaml \
  --output validation.html
```

Do not continue if validation reports errors. Warnings mean the run can proceed,
but the report and `logs/warnings.tsv` must be interpreted carefully.

## 4. Run

```bash
Rscript exec/tumoronly run \
  --input variants.tsv \
  --config config.yaml \
  --output results
```

The same workflow is available from R:

```r
library(tumoronly)

result <- run_tumoronly(
  input = "variants.tsv",
  config = "config.yaml",
  output_dir = "results"
)
```

## 5. Open the report

The run writes:

```text
results/<run_id>/report.html
```

If needed, re-render:

```bash
Rscript exec/tumoronly report --run-dir results/<run_id>
```

## 6. Read the key outputs

Start with:

- `report.html`
- `tables/classified_variants.tsv`
- `tables/sample_summary.tsv`
- `tables/filter_audit.tsv`
- `logs/warnings.tsv`
- `run_manifest.json`

For each variant, inspect:

- `final_class`
- `somatic_score`
- `filter_status`
- `evidence_supporting_classification`
- `evidence_against_classification`
- `missing_evidence`
- `classification_explanation`

