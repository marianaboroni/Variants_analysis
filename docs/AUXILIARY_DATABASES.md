# Preparing auxiliary databases (ABraOM, COSMIC) before running `tumoronly`

Two one-time, standalone preparation steps turn raw, differently-represented
population/somatic databases into the exact inputs `tumoronly`'s existing
matching engines (`R/population_brazilian.R`, `R/cosmic.R`) already know how to
consume. Run them **once per database release**; every later `tumoronly run`
across any number of samples/cohorts reuses the result. Neither step touches
your sample VCFs — they only prepare reference databases.

```
raw ANNOVAR TSV (ABraOM/BIPMed)  ──▶ tumoronly prepare-abraom      ──▶ population.brazilian.db
raw COSMIC v104 release files    ──▶ tumoronly build-cosmic-input  ──▶ cosmic.raw_file
                                                                          │
                                                                          ▼
                                                          tumoronly prepare-cosmic (existing)
```

Both are idempotent (checksummed manifest; re-run is a no-op unless an input
changed or `--force` is passed) and never drop rows silently — anything that
could not be converted or matched is written to a `*.tsv.gz` provenance file
next to the output, never merged into the database and never just discarded.

## 1. `tumoronly prepare-abraom` — Brazilian population database

### What it solves

ABraOM (and ANNOVAR-derived population exports in general) represent indels
**without a VCF anchor base**: a deletion gives the deleted sequence as `Ref`
and `-` as `Alt` at the position of the first deleted base; an insertion gives
`-` as `Ref` and the inserted sequence as `Alt`. Handing this directly to
`population.brazilian.db` (which matches by the canonical key
`BUILD|CHROM|POS|REF|ALT`) would silently miss every indel — the columns
(`Chr`/`Start`/`Ref`/`Alt`) don't even have the names `load_brazilian_db()`
looks for, and the representation itself isn't VCF-standard.

`prepare-abraom` reconstructs a VCF-standard representation using the
reference FASTA, **then re-validates every reconstructed REF against the same
FASTA** before trusting it. The exact reconstruction rule (empirically
confirmed against real Ensembl/GRCh38 coordinates for two ABraOM indels with
dbSNP IDs — rs56289060 and rs749280879 — not assumed from generic ANNOVAR
documentation; see the header of `R/population_brazilian_prepare.R`):

| variant type | ANNOVAR `Start` means | VCF reconstruction |
|---|---|---|
| SNV / MNP | the position itself | unchanged |
| deletion (`Alt == "-"`) | first deleted base | `pos = Start - 1`, `ref = anchor + Ref`, `alt = anchor` |
| insertion (`Ref == "-"`) | **the anchor position itself** | `pos = Start`, `ref = anchor`, `alt = anchor + Alt` |

The asymmetry (deletion shifts back one base, insertion does not) is real —
verify it yourself with `Rscript -e 'source("R/utils.R")'`-style exploration
before trusting it blindly on a database you haven't checked, the same way
this derivation was checked here rather than assumed.

### Inputs you must provide

| Config key | What | Example |
|---|---|---|
| `population.brazilian.raw_file` | raw ANNOVAR-format TSV | `databases/abraom/SABE1171.Abraom.clean.tsv` |
| `reference.fasta` | reference FASTA matching `population.brazilian.genome_build` | `referencias/Homo_sapiens_assembly38.fasta` |
| `population.brazilian.genome_build` | build of both the TSV and the FASTA | `GRCh38` |
| `population.brazilian.name` | short label, becomes part of the output path | `ABraOM_SABE` |
| `population.brazilian.cache_dir` | where the processed DB is written | `db/abraom` (default) |

Minimum accepted raw-TSV columns (case-sensitive; ANNOVAR names are the
primary match, generic aliases also work — see `standardize_brazilian_raw()`):
`Chr`, `Start`, `Ref`, `Alt`, and, when present, `Allele_number` (→ `AN`),
`Allele_ALT_count` (→ `AC`), `Frequencies` (→ `AF`), `HomozygousALT_count`
(→ `HOM`), `avsnp150` (→ `rsID`).

### Running it

```bash
tumoronly prepare-abraom --config config/example.yml
```

Requires `Biostrings` (`BiocManager::install("Biostrings")`) — the reference
FASTA is loaded once into memory (`Biostrings::readDNAStringSet`), so expect
memory use on the order of the FASTA size (a full GRCh38 primary assembly is
~3 GB uncompressed).

### Output

```
db/abraom/<name>_<build>/
  abraom_db.tsv.gz        CHROM POS REF ALT AF AC AN HOM Variant_Type rsID
                          Source_FILTER Original_CHROM Original_POS
                          Original_REF Original_ALT
  manifest.json           checksums, counts by variant type, ref_validated
  ref_mismatch.tsv.gz      rows where the reconstructed REF did not match the
                           FASTA — excluded from abraom_db.tsv.gz, never
                           silently dropped
  unresolved.tsv.gz        indels on a chromosome absent from the FASTA
  dropped_missing.tsv.gz   rows missing chrom/pos/ref/alt
  duplicates.tsv.gz        rows collapsing to the same canonical key
```

Point `population.brazilian.db` at `abraom_db.tsv.gz`. **Check
`n_ref_mismatch` in the manifest** — a non-trivial fraction (more than the
occasional edge case) is the signal that something about the FASTA
build/contig-naming doesn't match this database, not something to ignore.

## 2. `tumoronly build-cosmic-input` — COSMIC v104 flat input

### What it solves

`prepare_cosmic_db()` (existing) expects **one flat TSV** with, per row, both
genomic coordinates and tumor-context columns
(`standardize_cosmic_raw()`: `CHROMOSOME`, `GENOME_START`,
`GENOMIC_WT_ALLELE`, `GENOMIC_MUT_ALLELE`, `GENOMIC_MUTATION_ID`,
`GENE_SYMBOL`, `PRIMARY_SITE`, `PRIMARY_HISTOLOGY`, `HISTOLOGY_SUBTYPE_1`).
COSMIC v104 ships these as **separate files**, joined only by ID columns:

- mutation occurrences (`Cosmic_GenomeScreensMutant_v104_GRCh38.tsv`,
  `Cosmic_NonCodingVariants_v104_GRCh38.tsv`) — one row per (mutation ×
  sample), but their own `GENOME_START`/`GENOMIC_WT_ALLELE`/`GENOMIC_MUT_ALLELE`
  are the **pre-normalization** representation (same convention as
  `OLD_VARIANT` in the Normal VCF) — not safe to use directly for matching
  against a `bcftools norm`-normalized sample VCF.
- normalized coordinates (`Cosmic_*_Normal_v104_GRCh38.vcf`) — already
  5'-shifted per the VCF standard (COSMIC's own README cites the same
  normalization reference this project uses), one row per unique
  `GENOMIC_MUTATION_ID` (COSV id, repeated once per overlapping
  transcript/gene — deduplicated by this script).
- tumor context (`Cosmic_Classification_v104_GRCh38.tsv.gz`) —
  `COSMIC_PHENOTYPE_ID → PRIMARY_SITE/PRIMARY_HISTOLOGY/HISTOLOGY_SUBTYPE_*`,
  needed for `docs/COSMIC_DECISION_RULES.md`'s tumor-type-stratified matching;
  without it every variant would be `cosmic_tumor_type_missing`.

`build-cosmic-input` joins these **once**: coordinates come from the Normal
VCF (matches the normalization your sample VCFs should already have), the
mutation TSVs contribute only relational columns (gene, sample, study, PubMed,
somatic status), and Classification adds tumor context — producing exactly
the schema `standardize_cosmic_raw()` expects, with correct row cardinality
(one row per real sample occurrence, which is what
`build_cosmic_long()`/`docs/COSMIC_DECISION_RULES.md`'s recurrence counting
needs — not one row per mutation).

**`Cosmic_Sample_v104_GRCh38.tsv.gz` is deliberately NOT integrated.** No part
of `tumoronly`'s current COSMIC engine consumes its per-individual fields
(age, MSI, stage, ploidy, ...) — bolting it on here without a downstream
consumer would be dead data. If a future feature needs it, join it in then.

### Inputs you must provide

```yaml
cosmic:
  release: v104
  cache_dir: db/cosmic
  input_sources:
    genome_screens_tsv: databases/cosmic/Cosmic_GenomeScreensMutant_v104_GRCh38.tsv
    noncoding_tsv: databases/cosmic/Cosmic_NonCodingVariants_v104_GRCh38.tsv
    genome_screens_normal_vcf: databases/cosmic/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf
    noncoding_normal_vcf: databases/cosmic/Cosmic_NonCodingVariants_Normal_v104_GRCh38.vcf
    classification_tsv: databases/cosmic/Cosmic_Classification_v104_GRCh38.tsv.gz
```

### Running it

```bash
tumoronly build-cosmic-input --config config/example.yml
# then point cosmic.raw_file at the printed path and run:
tumoronly prepare-cosmic --config config/example.yml
```

### Performance — read this before running on the full release

COSMIC v104's genome-screens files are **~60M+ rows / ~45 GB combined** across
the 4 input files. This is a one-time batch job, not an interactive one:

- Only the columns actually needed are read (`data.table::fread(..., select=)`),
  not every column — this cuts I/O/memory by roughly 3-4x versus a naive read,
  but full-scale runs should still be expected to take from tens of minutes to
  a few hours and to need double-digit GB of RAM, depending on your machine.
  Run it with `nohup`/in the background, on a machine with adequate RAM — not
  as a "wait for it in the terminal" step.
- If you only care about specific genes/regions (e.g. a targeted panel or a
  known driver gene list), pre-filter the 4 input files by chromosome/region
  with standard tools (`awk`/`grep`/`bcftools view -r`) **before** pointing
  `input_sources` at them — this script does not do positional filtering
  itself, by design, so that a COSMIC preparation done once stays valid for
  future, larger, or differently-scoped cohorts rather than being scoped to
  today's samples.

### Output

```
db/cosmic/<release>/raw_input/
  cosmic_flat_input.tsv.gz         set this as cosmic.raw_file
  manifest.json                    checksums, counts, coordinate_source
  unmatched_to_normal_vcf.tsv.gz   mutation-TSV rows whose GENOMIC_MUTATION_ID
                                   had no match in either Normal VCF — excluded
                                   from the flat input, never silently dropped
```

## 3. Known limitation: `cosmic_tumor_type_mapping.tsv` histology values

While building and testing `build-cosmic-input` against real
`Cosmic_Classification_v104_GRCh38.tsv.gz` rows, the shipped
`config/cosmic_tumor_type_mapping.tsv` was found to have **incorrect
`cosmic_histology` values for every row except `MELANOMA`**. The matching
engine (`R/cosmic_context.R::classify_cosmic_tumor_context()`) compares
`cosmic_histology` against COSMIC's `PRIMARY_HISTOLOGY` field — and for
carcinomas, COSMIC's `PRIMARY_HISTOLOGY` is the generic `carcinoma`; the
specific descriptor (`serous_carcinoma`, `ductal_carcinoma`,
`endometrioid_carcinoma`, ...) lives in `HISTOLOGY_SUBTYPE_1`, a field the
current matching engine does not compare at all. The template's original
values (verified against real Classification rows for ovary/breast/lung/
endometrium/large_intestine) would never have produced `exact_match` — this
has been corrected (all `cosmic_histology` values set to `carcinoma` for the
carcinoma-type rows; `HGSOC → ovary/carcinoma` is now correct for this
project's ovarian cohort). Subtype-level discrimination (e.g. serous vs. clear
cell ovarian carcinoma) is **not currently implemented** in the matching
engine — `input_tumor_subtype` only marks which mapping row was used, it does
not check `HISTOLOGY_SUBTYPE_1`. Treat any future addition of new tumor types
to this table with the same verification, not by pattern-matching the
existing (now-corrected) rows blindly.
