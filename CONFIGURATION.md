# Configuration

Create a template:

```bash
Rscript exec/tumoronly init --output config.yaml
```

All v2 defaults are centralized in `tumoronly_default_config()`.

## Minimal config

```yaml
input:
  path: variants.tsv
  genome_build: GRCh38

analysis:
  output_dir: results

technical_filters:
  min_depth: 20
  min_alt_count_snv: 5
  min_alt_count_indel: 8
  min_af_snv: 0.03
  min_af_indel: 0.05
  min_tlod: 6
  min_mbq: 25
  min_mmq: 40
```

## Important sections

| Section | Purpose |
|---|---|
| `input` | file paths, genome build, column mapping |
| `analysis` | output directory, run id, threads |
| `hard_filters` | fail-safe technical gate |
| `technical_filters` | depth, VAF, quality thresholds |
| `population_filters` | global AF thresholds and columns |
| `cohort` | recurrence thresholds and minimum cohort size |
| `scoring` | classification score thresholds |
| `driver_resources` | optional driver-gene and hotspot references |
| `cosmic` | optional COSMIC database configuration |
| `report` | report and top-variant display options |

## Reproducibility

Each run writes:

- `config_used.yaml`
- `config.resolved.yml`
- `run_manifest.json`
- `session_info.txt`

These files are the first place to look when reproducing or auditing a run.

