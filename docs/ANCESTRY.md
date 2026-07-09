# Genetic-ancestry inference in tumor-only data (EXPERIMENTAL)

`R/ancestry_infer.R` + `R/ancestry_plots.R`. Genetic ancestry is **not race or
ethnic identity**; it is used only to weight/choose relevant population
frequencies and to flag under-representation. It never filters variants or changes
biological priority.

## Method

1. **Separate SNP track** taken from ALL variants (before the somatic filter), so
   the allele spectrum is not biased by somatic prioritization.
2. **Candidate selection**: biallelic autosomal SNVs, adequate depth, excluding
   hotspots/drivers.
3. **Genotype from VAF** (0/1/2) with confidence; **CNV/LOH-like** heterozygous
   imbalance loci are flagged and excluded (`possible_loh`/`possible_cnv`).
4. **Allele/strand alignment** to a versioned AIMs panel; palindromic (A/T, C/G)
   and non-panel SNPs excluded.
5. **Proportions** by non-negative, sum-to-one least squares against
   reference-population alt-allele frequencies (AFR/EUR/NAT/EAS/SAS) — continuous,
   never a single categorical identity; plus a **PCA projection** for QC/plots.
6. **QC gating**: `ANCESTRY_INFERENCE_STATUS = not_evaluable` when
   `ANCESTRY_SNPS_USED < ancestry.qc.minimum_snps`; confidence tier
   (high/moderate/low) from SNP count + call rate.

PLINK / ADMIXTURE / SNPRelate are **optional external backends** (checked by
`tumoronly doctor`); the shipped fallback is self-contained (base R).

## Outputs

`tables/ancestry_summary.tsv` (proportions + confidence + SNP counts + PCs),
`tables/ancestry_snps.tsv.gz`, `tables/ancestry_qc.tsv`, `tables/ancestry_projection.tsv`.
Plots (`plots/ancestry/`, PNG+PDF, manifest): proportions, SNP QC funnel, exclusion
reasons, SNPs by chromosome, VAF distribution, call rate, PCA projection — each
generated only with sufficient data, else skipped with a reason.

## Config

```yaml
ancestry:
  enabled: true
  marker_panel: { path: db/ancestry/aims_grch38.tsv, version: aims_v1, genome_build: GRCh38 }
  qc: { minimum_depth: 10, minimum_alt_reads: 3, minimum_snps: 100 }
  dominant_component_threshold: 0.80
  admixed_minimum_secondary_component: 0.15
```

## Limitations / status: EXPERIMENTAL

Tumor purity, CNV and LOH distort VAF-based genotypes; targeted panels have too few
AIMs (prefer `not_evaluable` or broad labels); `admixed` is expected and common for
Brazilian samples and is not a lower-quality category. Validate against matched
normal / array / germline WGS before research or production use. Not diagnostic;
not race/ethnicity.
