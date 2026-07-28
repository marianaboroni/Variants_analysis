# Getting Started With tumoronly v2

Use this page when you already have a checkout of the repository. For the
shortest path, see `QUICKSTART.md`; for the full walkthrough, see
`TUTORIAL.md`.

## 1. Install

```bash
R CMD INSTALL .
Rscript -e 'library(tumoronly); check_installation()'
```

## 2. Create a Config

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

Or start from the runnable demo:

```bash
cp config/demo_v2.yml config/my_analysis.yml
```

## 3. Validate

```bash
Rscript exec/tumoronly validate \
  --input variants.tsv \
  --config config.yaml \
  --output validation.html
```

Strict validation blocks missing indispensable evidence. Use permissive mode
only when you want missing evidence documented rather than blocking execution:

```bash
Rscript exec/tumoronly validate \
  --input variants.tsv \
  --config config.yaml \
  --permissive
```

## 4. Run

```bash
Rscript exec/tumoronly run \
  --input variants.tsv \
  --config config.yaml \
  --output results \
  --run-id my_run
```

## 5. Read Outputs

Open `results/my_run/report.html`, then inspect:

- `run_manifest.json`
- `logs/warnings.tsv`
- `tables/classified_variants.tsv`
- `tables/filter_audit.tsv`
- `tables/sample_summary.tsv`
- `figure_data/*.tsv`

The per-variant evidence trail is in:

- `evidence_supporting_classification`
- `evidence_against_classification`
- `missing_evidence`
- `classification_explanation`

## 6. Real Data Gate

The v2 implementation was gated on
`data/test/WGS_all_patients_first100k.tsv`. Details are in
`docs/REAL_DATA_VALIDATION.md`.
