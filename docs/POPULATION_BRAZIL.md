# Brazilian-population-aware population interpretation

`R/population_brazilian.R`. The Brazilian population is highly admixed and
under-represented in global databases, so **absence in gnomAD is not evidence of
somaticity**. Global and Brazilian evidence are computed and reported separately.

## Outputs (per variant)

`GLOBAL_POPULATION_STATUS`, `BRAZILIAN_POPULATION_STATUS`,
`POPULATION_EVIDENCE_STATUS` (common_global / common_brazilian /
rare_global_and_brazilian / absent_global_present_brazilian /
absent_brazilian_present_global / absent_in_checked_databases /
insufficient_coverage / not_evaluable), `POPULATION_EVIDENCE_CONFIDENCE`,
`POPULATION_DATABASES_CHECKED`, `POPULATION_COVERAGE_STATUS`; `MAX_GLOBAL_AF`,
`MAX_ANCESTRY_SPECIFIC_AF`, `BRAZILIAN_AF/AC/AN/HOM_COUNT`,
`BRAZILIAN_DATABASE_SAMPLE_SIZE`, `LOCAL_CONTROL_AF`;
`BRAZILIAN_DB_VARIANT_OBSERVED/REGION_COVERED/ABSENCE_INTERPRETABLE`;
`POSSIBLE_BRAZILIAN_GERMLINE` + `BRAZILIAN_GERMLINE_EVIDENCE` +
`BRAZILIAN_GERMLINE_REVIEW_REASON`; `BRAZILIAN_POPULATION_EVIDENCE_SCORE`.

## Key rules

- Matching by canonical key `BUILD|CHROM|POS|REF|ALT` (never rsID/gene/position alone).
- **Absence is used as evidence only when interpretable** (`BRAZILIAN_DB_ABSENCE_INTERPRETABLE`);
  otherwise `insufficient_coverage`.
- **Documented REVIEW rule** (`apply_brazilian_population_rules`, configurable via
  `population.brazilian.escalate_to_review`): a variant absent globally but present
  in a Brazilian DB, or a predisposition-gene variant with heterozygous VAF, or one
  present in local controls, is flagged `POSSIBLE_BRAZILIAN_GERMLINE` and its STATUS
  is escalated PASS → REVIEW (not auto-somatic, not auto-FAIL). This is the only
  population path that changes STATUS, and it is explicit and auditable.
- Allele count / cohort size are preserved; a single observation in a small cohort
  is not treated as a precise frequency (`present_low_count_brazilian`).
- `BRAZILIAN_POPULATION_EVIDENCE_SCORE` is separate from biological relevance and
  is never used to declare ancestry.

## Config

```yaml
population:
  global: { common_af: 0.01, rare_af: 0.001 }
  brazilian:
    db: db/abraom_sabe_grch38.tsv      # user-provided (ABraOM/BIPMed); not redistributed
    genome_build: GRCh38               # ABraOM SABE-WGS=GRCh38, SABE-WES=GRCh37
    common_af: 0.01
    review_af: 0.001
    minimum_allele_count: 2
    escalate_to_review: true
  local_controls: { path: null, artifact_or_germline_af: 0.005 }
```

## Outputs / tables

`brazilian_population_evidence.tsv`, `brazilian_population_review.tsv`,
`possible_brazilian_germline.tsv`, `global_brazilian_frequency_conflicts.tsv`,
`local_normal_recurrence.tsv`.

## Status: EXPERIMENTAL

Thresholds are heuristic and cohort-size dependent; validate against a Brazilian
cohort with matched normal or orthogonal validation before production use. Data
handling follows minimization/LGPD principles: no direct identifiers, broad
regions only, local execution (no external calls except optional OncoKB with
minimal data).
