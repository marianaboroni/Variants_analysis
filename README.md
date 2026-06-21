# Tumor-only Somatic Variant Filtering Pipeline

R package for prioritizing tumor-only somatic variants in WGS/WES/panel sequencing, particularly for cohorts without paired normal samples or panel of normals. Designed with Brazilian cohort context in mind.

## Overview

The pipeline computes:

- Technical features from Mutect2/DRAGEN/VEP
- Internal cohort recurrence (pseudo-panel of normals)
- Recurrence by tumor type (when OncoTree metadata is available)
- Adaptive hard filters by sample and assay type (WGS, WES, PANEL)
- Per-sample variant calling QC
- Semi-supervised ML model for true positive prioritization
- Tumor mutational burden (TMB) per sample
- Clonality/CCF estimates when purity and local copy number are available
- Germline and artifact probability scores
- Somatic probability scores
- Auditable final classification
- External validation with OncoKB and COSMIC
- Driver vs. passenger classification (separate from somatic classification)
- Export for external models (dNdScv, CHASMplus)
- QC, TMB, driver/passenger, recurrence, clonality, and oncoplot figures

When `maftools` is installed, native MAF visualization plots are generated; otherwise, a fallback oncoplot is created using `ggplot2`.

**Important:** OncoKB and COSMIC are used for validation/prioritization, not as standalone proof of somaticity.



## Input Format

A TSV/CSV table with one row per variant per sample. Columns can follow MAF, VEP, or custom formats. The pipeline expects to find or infer:

- **Sample ID:** `Tumor_Sample_Barcode`, `sample`, `Sample`, or similar
- **Coordinates:** `CHROM`, `START`, `REF`, `ALT`
- **Gene:** `Hugo_Symbol`, `SYMBOL`, `Gene`
- **Consequence:** `Consequence`, `Variant_Classification`
- **Depth:** `DP`
- **Alt support:** `alt_count`, `t_alt_count`, `AD`
- **Tumor AF:** `AF`, `VAF`, `tumor_f`
- **Quality:** `TLOD`, `MBQ`, `MMQ`, `STRANDQ`, `STRQ`, `FILTER`
- **Population frequencies (VEP/gnomAD):** `gnomADe_AF`, `gnomADg_AF`, `AFR_AF`, `AMR_AF`, etc.
- **Tumor type:** `tumor_type`, `ONCOTREE_CODE`, `Cancer_Type`, `Tumor_Type`, or via `input.sample_metadata`

Missing fields are treated as `NA` with minimal or neutral penalty. Preserving all technical fields from the VCF improves classification.

For Brazilian cohorts, configure local population frequency columns when available (e.g., `ABraOM_AF`, `BIPMed_AF`). The pipeline uses the maximum AF observed across configured population columns to reduce false positives from variants under-represented in global databases.



## Quick Start

A Bioconductor-style vignette is available in [vignettes/tumor_only_variant_filtering.Rmd](vignettes/tumor_only_variant_filtering.Rmd). It describes pipeline usage, output structure, and how to interpret QC, ML results, TMB, driver classification, clonality, and figures.

A workflow schematic is available in [figures/pipeline_fluxograma.svg](figures/pipeline_fluxograma.svg) (editable Mermaid version at [figures/pipeline_fluxograma.mmd](figures/pipeline_fluxograma.mmd)).

Edit [config/example_config.yml](config/example_config.yml) and run:

```bash
Rscript scripts/run_automated_variant_analysis.R --config config/example_config.yml
```

Or use the explicit alias:

```bash
Rscript scripts/run_tumor_only_filter.R --config config/example_config.yml
```

## Key Output Files

**Classification & QC:**
- `sample_qc_summary.tsv` — Per-sample variant calling QC
- `variant_call_qc_summary.tsv` — Variant-level QC metrics
- `variants_with_features.tsv` — Full feature matrix
- `variants_scored.tsv` — Variants with all scores

**Filtered Sets:**
- `variants_high_confidence_somatic.tsv` — Pass all filters
- `variants_probable_somatic.tsv` — Likely somatic, lower confidence
- `variants_likely_passenger.tsv` — Somatic but unlikely driver
- `variants_driver_candidates.tsv` — Candidate drivers
- `variants_uncertain_for_review.tsv` — Manual review candidates
- `removed_probable_germline.tsv` — Filtered as likely germline
- `removed_probable_artifact.tsv` — Filtered as likely artifact

**ML & Recurrence:**
- `ml_filter_metrics.tsv` — ML model performance
- `ml_training_labels.tsv` — Training labels used
- `active_learning_candidates.tsv` — Candidates for manual curation
- `cohort_recurrence_table.tsv` — Internal variant recurrence
- `tumor_type_recurrence_table.tsv` — Recurrence by tumor type

**External Model Exports:**
- `dndscv_input.tsv` — dNdScv positive selection model
- `chasmplus_missense_input.tsv` — CHASMplus driver prediction
- `oncokb_cosmic_validation_summary.tsv` — OncoKB/COSMIC matches

**TMB, Drivers & Clonality:**
- `tmb_summary.tsv` — Tumor mutational burden per sample
- `driver_classification_summary.tsv` — Driver classification summary
- `clonality_summary.tsv` — CCF/clonality estimates

**Manifest & Figures:**
- `figure_manifest.tsv` — Metadata for all generated figures (file, type, plot category, status, notes)
- `figures/` subdirectory with PNG and PDF versions of QC, TMB, driver, clonality plots, and optional maftools visualizations



## Validation with OncoKB

Preferably use output already annotated by OncoKB Annotator/MafAnnotator or a local reference standardized by [scripts/build_oncokb_cosmic_reference.R](scripts/build_oncokb_cosmic_reference.R). The file should support matching by coordinate or by `sample + gene + protein_change`. If `tumor_type`/OncoTree is available, the pipeline prioritizes tumor-type-specific matches before pan-cancer matches.

Recognized columns:
- `Tumor_Sample_Barcode`, `Sample`
- `Hugo_Symbol`, `SYMBOL`
- `HGVSp_Short`, `Protein_Change`, `HGVSp`
- `CHROM`, `START`, `REF`, `ALT`
- `ONCOGENIC`, `MUTATION_EFFECT`, `HIGHEST_LEVEL`, `LEVEL_1`, `LEVEL_2`, `LEVEL_3A`, `LEVEL_3B`, `LEVEL_4`

## Validation with COSMIC

Use a licensed/exported COSMIC table or VEP annotation with COSMIC identifiers in `Existing_variation`. COSMIC data is not auto-downloaded due to licensing requirements. Recognized columns:
- `COSMIC_ID`, `COSMIC_MUTATION_ID`, `LEGACY_MUTATION_ID`
- `Hugo_Symbol`, `Gene`
- `HGVSp_Short`, `Mutation AA`, `Protein_Change`
- `CHROM`, `START`, `REF`, `ALT`
- `COSMIC_COUNT`, `FATHMM_PREDICTION`, `CNT`

To build standardized local references:

```bash
Rscript scripts/build_oncokb_cosmic_reference.R --config config/reference_build_example.yml
```

This generates `oncokb_reference.tsv` and/or `cosmic_reference.tsv` in the configured directory, which can be used in `validation.oncokb` and `validation.cosmic` config settings.



## Driver Classification

The pipeline separates two distinct questions:

1. **`final_class`**: Does the variant appear real/somatic, or likely germline/artifact?
2. **`driver_class`**: Among high-confidence somatic variants, is there evidence of driver status?

**Classification levels:**
- `known_driver` — Strong evidence in OncoKB/COSMIC/curated hotspots + high-confidence somatic
- `probable_driver` — Driver gene compatible with tumor type + coherent functional consequence
- `possible_driver` — Partial evidence, requires manual review
- `likely_passenger` — High-confidence somatic but no driver support
- `not_evaluable_as_driver` — Non-somatic, likely germline, likely artifact, or hard filter fail

**Score sources (when available):**
- OncoKB: `ONCOGENIC`, `MUTATION_EFFECT`, `HIGHEST_LEVEL`
- COSMIC: `COSMIC_ID`, `COSMIC_COUNT`
- Driver gene references: CGC, IntOGen, NCG, CGI, or custom table
- Hotspot references: OncoKB, COSMIC, Cancer Hotspots, or custom table
- Predictors: `SpliceAI`, `CADD`, `AlphaMissense`, `REVEL`, `MetaRNN`, `MetaLR`, `MutationTaster`, `SIFT`, `PolyPhen`
- Structure/domain: `DOMAINS`, `COSMIC3D`, `AlphaFold_feature`, `functional_region`, `structure_feature`

**Example driver genes file** (`driver_resources.driver_genes`):
```
gene	tumor_type	role	source	confidence
TP53	UCEC	tumor_suppressor	CGC;IntOGen	canonical
PIK3CA	UCEC	oncogene	CGC;IntOGen	canonical
KRAS	PANCANCER	oncogene	CGC;OncoKB	canonical
```

**Example hotspots file** (`driver_resources.hotspots`):
```
tumor_type	Hugo_Symbol	HGVSp_Short	source	evidence
PANCANCER	BRAF	p.V600E	OncoKB;COSMIC	canonical_hotspot
UCEC	KRAS	p.G12V	OncoKB;COSMIC	recurrent_hotspot
```

The `dndscv_input.tsv` and `chasmplus_missense_input.tsv` files are generated from high-confidence and probable somatic variants for external re-evaluation via positive selection and driver prediction models.



## Incremental ML for False Positive Reduction

The pipeline includes an incremental ML layer (`ml_filter`) combining:

- Manual review labels in `ml_filter.review_labels`
- Auditable pseudo-labels from pipeline rules
- Higher weights for manual curation (`manual_label_weight`) vs. pseudo-labels (`pseudo_label_weight`)
- Global pan-cancer model
- Tumor-type-specific models when sufficient labels are available

Benefits:
- Leverage curated manual reviews for learning
- Transparent rules enable interpretation of decisions
- Active learning identifies uncertain candidates for manual review
- Incremental curation improves future predictions

## Technical Architecture

**Workflow stages** ([R/workflow.R](R/workflow.R)):
1. Input standardization + sample metadata
2. Feature extraction (technical + population-based)
3. Hard filtering (sample/assay-adaptive)
4. Scoring (conservative classification + manual review rules)
5. ML layer (auxiliary to hard filters)
6. Clonality/CCF estimation
7. TMB calculation
8. Driver classification
9. Reporting (filter reports, metrics, exports)
10. Visualization (figures + manifest)

**Configuration** ([config/](config/)):
- YAML-driven settings for all stages
- Sample-specific overrides
- Population frequency column mapping for Brazilian cohorts
- OncoKB/COSMIC reference paths

**Figure manifest** (`figure_manifest.tsv`):
Each generated figure is logged with:
- `file` — Output path
- `type` — Figure type identifier
- `plot_category` — High-level category (maftools, oncoplot, classification, tmb, driver, qc, recurrence, clonality)
- `status` — written, skipped, failed
- `note` — Generation notes or failure reason



Rotulos manuais reconhecidos:

- positivos: `true_somatic`, `validated_somatic`, `true_positive`, `TP`;
- negativos: `false_positive`, `artifact`, `germline`, `polymorphism`,
  `technical`, `not_somatic`;
- incertos: `uncertain`, `ambiguous`, `unknown`, `not_evaluable`, que nao entram
  no treino.

Exemplo de `review_labels.tsv`:

```text
sample_id	tumor_type	chrom	pos	ref	alt	gene	protein_change	reviewed_label	evidence	reviewer	review_date	notes
S3	UCEC	7	140453136	A	T	BRAF	p.V600E	true_somatic	IGV_hotspot	MB	2026-06-21	adequate depth and VAF
S4	UCEC	7	140453136	A	T	BRAF	p.V600E	artifact	low_depth_low_quality	MB	2026-06-21	technical failure
```

O pareamento dos rotulos e feito por prioridade:

```text
sample+coord -> sample+gene/protein -> tumor_type+coord ->
tumor_type+gene/protein -> coord -> gene/protein
```

Quando ha rotulos positivos e negativos suficientes, o pipeline treina:

- `glm_PANCANCER`: modelo global usando todos os tipos tumorais;
- `glm_<TUMOR_TYPE>`: modelo especifico, por exemplo `glm_UCEC`, quando aquele
  tumor tem exemplos suficientes.

Na predicao final, o modelo especifico do tipo tumoral e usado primeiro. Se ele
nao existir, o pipeline usa o modelo pan-cancer como fallback.

Saidas principais:

- `ml_label`: rotulo usado no treino;
- `ml_label_source`: `manual_review` ou `pseudo_rules`;
- `ml_true_positive_probability`: probabilidade final usada pelo pipeline;
- `ml_pancancer_true_positive_probability`: probabilidade do modelo global;
- `ml_tumor_type_true_positive_probability`: probabilidade do modelo tumoral;
- `ml_model_scope`: `tumor_type` ou `pancancer`;
- `ml_model_id`: identificador do modelo usado;
- `ml_training_labels.tsv`: auditoria dos exemplos usados no treino;
- `active_learning_candidates.tsv`: variantes sugeridas para a proxima rodada
  de revisao manual;
- `ml_filter_metrics.tsv`: status, contagem de rotulos e AUC de treino por
  escopo.

Se nao houver positivos e negativos suficientes, o pipeline nao forca um modelo:
ele grava `status = insufficient_training_labels` em `ml_filter_metrics.tsv` e
ainda gera `active_learning_candidates.tsv` para orientar a curadoria. Isso evita
treinar um classificador aparentemente sofisticado, mas biologicamente vazio.

O ciclo recomendado e:

```text
rodar pipeline -> revisar active_learning_candidates.tsv ->
atualizar review_labels.tsv -> rerodar pipeline
```

Com o tempo, os modelos especificos por tumor ficam mais fortes; quando varios
tipos tumorais acumulam rotulos confiaveis, o modelo pan-cancer passa a aprender
padroes compartilhados entre tumores.

## TMB

O TMB e calculado a partir de variantes `high_confidence_somatic` e
`probable_somatic` com consequencias proteicas/splice contaveis. Configure o
denominador em `tmb.callable_mb`.

Exemplos:

- WES/coding TMB: geralmente usar o tamanho chamavel do exoma/captura em Mb.
- Painel: usar o tamanho real do painel em Mb.
- WGS: defina explicitamente se quer TMB coding-like ou genome-wide; o default
  do exemplo usa `30 Mb` para uma leitura coding-like.

## Clonalidade

O pipeline estima `ccf_estimate` usando VAF, pureza tumoral e copy number local
quando esses campos existem. Colunas reconhecidas incluem:

- pureza: `tumor_purity`, `purity`, `PURITY`, `Tumor_Purity`
- copy number total/local: `total_cn`, `Total_CN`, `CN`, `copy_number`,
  `local_cn`, `tcn`
- multiplicidade da mutacao: `multiplicity`, `mutation_multiplicity`,
  `mut_cn`, `mutation_cn`

Formula usada:

```text
CCF = VAF * (purity * total_cn + (1 - purity) * normal_cn) /
      (purity * mutation_multiplicity)
```

Se houver pureza, mas nao houver copy number, o pipeline assume regiao
copy-neutral (`total_cn = 2`) e marca o metodo como
`purity_adjusted_copy_neutral`. Se nao houver pureza, ele nao calcula CCF formal:
usa apenas uma classificacao proxy por VAF (`vaf_proxy_no_purity`). Isso e
intencional para preservar a interpretabilidade em tumor-only.

Classes principais:

- `clonal`: CCF acima de `clonality.clonal_ccf_cutoff`.
- `subclonal`: CCF abaixo de `clonality.subclonal_ccf_cutoff`.
- `intermediate`: zona intermediaria.
- `clonal_like_high_vaf` e `subclonal_like_low_vaf`: classes proxy quando nao
  ha pureza.

Tambem sao gerados clusters simples por amostra em
`clonality_cluster_id`, `clonality_cluster_center` e
`clonality_cluster_label`. Para inferencia clonal definitiva, prefira integrar
pureza/CNV de ferramentas dedicadas como ABSOLUTE, FACETS, Sequenza, PURPLE ou
ichorCNA.

## Figuras e maftools

As figuras ficam no subdiretorio `figures/` e sao listadas em
`figure_manifest.tsv`, que tambem informa quando uma figura foi pulada. O
arquivo `figures/somatic_maftools_input.maf` e sempre criado com variantes
`high_confidence_somatic` e `probable_somatic`.

Para gerar os plots nativos do `maftools`, instale o pacote no ambiente R usado
pelo pipeline. Sem `maftools`, o workflow continua funcionando e gera um
oncoplot fallback em `ggplot2`.

## Tipos tumorais

Para aplicar o pipeline a diferentes tumores, forneca `tumor_type` na tabela
principal ou em `input.sample_metadata`. Use preferencialmente codigos OncoTree
ou abreviacoes consistentes, como `UCEC`, `BRCA`, `LUAD`, `COADREAD`.

O pipeline sempre calcula dois sinais de recorrencia:

- pan-coorte: bom para detectar artefatos sistematicos e polimorfismos;
- dentro do tipo tumoral: bom para auditar hotspots recorrentes em um contexto
  biologico especifico.

Mesmo quando OncoKB/COSMIC apoiam uma variante, ela continua passando pelos
hardfilters tecnicos; validacao externa aumenta prioridade, mas nao substitui
profundidade, VAF, qualidade de base/mapeamento e revisao de artefatos.

## Recomendacao de validacao

Depois da primeira rodada, revisar no IGV um conjunto balanceado:

- 50 `high_confidence_somatic`
- 50 `probable_germline`
- 50 `probable_artifact`
- 50 `uncertain`
- todos os casos `validation_rescue_candidate`
