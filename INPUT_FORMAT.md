# Input format

`tumoronly` accepts VCF, MAF, TSV, TXT, and gzipped variants of those formats.
The input may represent one sample or a cohort.

## Required biological fields

| Canonical field | Meaning | Examples |
|---|---|---|
| `CHROM` | chromosome | `chr17`, `17`, `X`, `MT` |
| `POS` | 1-based position | `7674220` |
| `REF` | reference allele | `C`, `AG` |
| `ALT` | alternate allele | `T`, `A` |

`input.genome_build` is also required unless it can be safely resolved from a
VCF header. Use `GRCh37` or `GRCh38`.

## Strongly recommended fields

| Field | Why it matters |
|---|---|
| sample identifier | needed for sample grouping and cohort recurrence |
| depth (`DP` or equivalent) | technical quality and fail-safe filters |
| alternate reads | technical quality and VAF derivation |
| VAF | scoring and interpretation |
| caller `FILTER` | caller-level technical evidence |
| gene and consequence | prioritization, report, and driver interpretation |
| population AF | germline-like evidence |
| mapping/base quality | technical artifact evidence |
| TLOD or caller LOD | caller-specific quality evidence |
| strand-bias metric | artifact evidence |

## Common TSV aliases

The validator detects common aliases. Examples:

| Canonical | Accepted aliases |
|---|---|
| `CHROM` | `Chromosome`, `chr`, `chrom` |
| `POS` | `Start_Position`, `START`, `Position` |
| `REF` | `Reference_Allele`, `ref` |
| `ALT` | `Tumor_Seq_Allele2`, `alt` |
| `SAMPLE_ID` | `Tumor_Sample_Barcode`, `Sample_Barcode`, `sample_id` |
| `GENE` | `Hugo_Symbol`, `SYMBOL`, `Gene.refGene` |
| `CONSEQUENCE` | `Consequence`, `Variant_Classification`, `Func.refGene` |

## Explicit YAML mapping

Use `input.column_map` when aliases are not enough:

```yaml
input:
  path: variants.tsv
  genome_build: GRCh38
  column_map:
    chrom: chromosome_name
    pos: start
    ref: reference
    alt: alternate
    sample_id: patient_sample
```

The mapper never silently changes scientific values. It only identifies which
source column should be read.

## Strict and permissive validation

Strict validation is the default:

```bash
Rscript exec/tumoronly validate --input variants.tsv --config config.yaml
```

Permissive mode records issues but allows more incomplete files to pass:

```bash
Rscript exec/tumoronly validate --input variants.tsv --config config.yaml --permissive
```

Permissive mode does not invent missing data. It only changes whether selected
missing-evidence conditions block execution.

## Cohorts

A cohort can be one multi-sample file or a YAML list of single-sample files:

```yaml
input:
  path:
    - sample1.vcf.gz
    - sample2.vcf.gz
  genome_build: GRCh38
```

Recurrence evidence is disabled below `cohort.min_cohort_size_for_recurrence`
because fractions are not meaningful in tiny cohorts.

