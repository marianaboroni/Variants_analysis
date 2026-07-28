# Configuration

Create a v2 template:

```bash
Rscript exec/tumoronly init --output config.yaml
```

All defaults are centralized in `tumoronly_default_config()`. A runnable demo is
available at `config/demo_v2.yml`.

## Minimal Config

```yaml
input:
  path: variants.tsv
  format: auto
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

## Important Sections

| Section | Purpose |
|---|---|
| `input` | file path, input format, genome build, sample metadata, column mapping |
| `analysis` | output directory, run id, threads |
| `validation` | strict/permissive behavior |
| `hard_filters` | fail-safe technical gate and accepted caller filters |
| `adaptive_filtering` | per-sample threshold review behavior |
| `technical_filters` | depth, alt-read, VAF, TLOD, base-quality, mapping-quality thresholds |
| `population_filters` | global population AF thresholds and source columns |
| `cohort` | recurrence thresholds and minimum sample count |
| `scoring` | score cutoffs for conservative tumor-only classes |
| `driver_resources` | optional driver-gene and hotspot references |
| `driver_scoring` | driver evidence thresholds |
| `cosmic` | optional processed COSMIC resource; missing COSMIC is reported as missing evidence |
| `tmb` | callable territory and burden/TMB thresholds |
| `clonality` | purity/CN/multiplicity defaults and clonality cutoffs |
| `ancestry` | optional AIMs panel and QC thresholds |
| `ml` | optional reviewed-label model registry and objective |
| `report` | report format, theme, and number of prioritized variants |

## Explicit Column Mapping

Autodetection handles common aliases. Use `input.column_map` when column names
are custom:

```yaml
input:
  path: variants.tsv
  genome_build: GRCh38
  column_map:
    chrom: chromosome_name
    pos: start_position
    ref: ref_allele
    alt: alt_allele
    sample_id: sample
```

Mapping identifies columns only. It does not transform allele values, change
positions, or fill missing evidence.

## Reproducibility Files

Each run writes:

- `config_used.yaml`
- `config.resolved.yml`
- `run_manifest.json`
- `session_info.txt`

Use these files to reproduce a run, audit thresholds, and record the package
version, commit, input checksum, warnings, runtime, and environment.

## Auxiliary Module Defaults

TMB requires an explicit callable territory to report mutations/Mb:

```yaml
tmb:
  enabled: true
  callable_mb: 30
  require_callable_mb: true
```

If `callable_mb` is missing, the run writes countable burden but marks TMB as
not evaluable.

Clonality uses purity/CN/multiplicity when available:

```yaml
clonality:
  enabled: true
  clonal_ccf_cutoff: 0.85
  subclonal_ccf_cutoff: 0.55
```

Without purity, the module writes VAF-proxy labels and records that limitation.

Ancestry requires a user-provided AIMs panel:

```yaml
ancestry:
  enabled: true
  marker_panel:
    path: data/aims_panel.tsv
    version: aims_v1
```

ML requires an activated model in the reviewed-evidence registry:

```yaml
ml:
  enabled: true
  db_dir: db/variant_evidence
  objective: P_TRUE_VARIANT
  use_for_filtering: false
```

ML predictions are auxiliary review evidence only.
