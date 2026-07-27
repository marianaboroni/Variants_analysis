# Tutorial: primera ejecución de `tumoronly` (cohorte de ovario, 7 muestras)

Esta guía asume la disposición real de este proyecto:

```
variant_analysis/
  Variants_analysis/        <- este paquete (raíz de todos los comandos abajo)
  mutect2/                  <- 7 VCF Mutect2+VEP, uno por carpeta de muestra
  databases/
    abraom/SABE1171.Abraom.clean.tsv
    cosmic/Cosmic_*_v104_GRCh38.*
  referencias/
    Homo_sapiens_assembly38.fasta (+ .fai)
```

Todos los comandos se ejecutan desde `Variants_analysis/` salvo que se indique lo contrario.

## Paso 0 — verificar el entorno

`bcftools` (para el filtrado/normalización, Paso 3) vive en un entorno conda separado del R base:

```bash
export PATH=/home/mmparedese/miniconda3/envs/vcfman/bin:$PATH   # da bcftools/tabix/samtools
bcftools --version   # 1.22 esperado
```

El pipeline R en sí (`tumoronly`) no necesita `bcftools` — corre en el R base del sistema. Verificá qué tiene disponible:

```bash
Rscript exec/tumoronly doctor
```

Esperado: los 4 paquetes `required` en OK (`data.table`, `yaml`, `jsonlite`, `digest`) y, para este tutorial, `Biostrings` en OK (lo necesita `prepare-abraom`). Si falta algo: `install.packages(...)` o `BiocManager::install("Biostrings")`.

## Paso 1 — preparar ABraOM (una sola vez)

Editá `config/example.yml` (o copiá a `config/ovarian_cohort.yml` y editá esa copia — recomendado, para no perder el ejemplo original):

```yaml
reference:
  fasta: /home/mmparedese/Desktop/variant_analysis/referencias/Homo_sapiens_assembly38.fasta

population:
  brazilian:
    raw_file: /home/mmparedese/Desktop/variant_analysis/databases/abraom/SABE1171.Abraom.clean.tsv
    name: ABraOM_SABE
    genome_build: GRCh38
    cache_dir: db/abraom
```

```bash
Rscript exec/tumoronly prepare-abraom --config config/ovarian_cohort.yml
```

**Esto es pesado**: el TSV de ABraOM tiene ~76 millones de líneas. Esperá varios minutos y uso de memoria en el orden del tamaño del FASTA (~3 GB, porque se carga completo para reconstruir los indels). Correlo con `nohup ... &` si vas a seguir trabajando en paralelo:

```bash
nohup Rscript exec/tumoronly prepare-abraom --config config/ovarian_cohort.yml > logs/prepare_abraom.log 2>&1 &
```

Al terminar, revisá `db/abraom/ABraOM_SABE_GRCh38/manifest.json` — el campo `n_ref_mismatch` debería ser una fracción pequeña del total, no una fracción grande (si es grande, algo no coincide entre el FASTA y la base, ver `docs/AUXILIARY_DATABASES.md`).

Completá en el config:
```yaml
population:
  brazilian:
    db: db/abraom/ABraOM_SABE_GRCh38/abraom_db.tsv.gz
```

## Paso 2 — preparar COSMIC (una sola vez)

```yaml
cosmic:
  release: v104
  cache_dir: db/cosmic
  source_build: GRCh38
  target_build: GRCh38
  input_sources:
    genome_screens_tsv: /home/mmparedese/Desktop/variant_analysis/databases/cosmic/Cosmic_GenomeScreensMutant_v104_GRCh38.tsv
    noncoding_tsv: /home/mmparedese/Desktop/variant_analysis/databases/cosmic/Cosmic_NonCodingVariants_v104_GRCh38.tsv
    genome_screens_normal_vcf: /home/mmparedese/Desktop/variant_analysis/databases/cosmic/Cosmic_GenomeScreensMutant_Normal_v104_GRCh38.vcf
    noncoding_normal_vcf: /home/mmparedese/Desktop/variant_analysis/databases/cosmic/Cosmic_NonCodingVariants_Normal_v104_GRCh38.vcf
    classification_tsv: /home/mmparedese/Desktop/variant_analysis/databases/cosmic/Cosmic_Classification_v104_GRCh38.tsv.gz
  tumor_type_mapping: config/cosmic_tumor_type_mapping.tsv
```

Dos comandos, en este orden — el primero prepara el TSV plano de entrada, el segundo construye la base final:

```bash
nohup Rscript exec/tumoronly build-cosmic-input --config config/ovarian_cohort.yml > logs/build_cosmic_input.log 2>&1 &
# cuando termine (revisar logs/build_cosmic_input.log): completá cosmic.raw_file con la ruta que imprime, luego:
nohup Rscript exec/tumoronly prepare-cosmic --config config/ovarian_cohort.yml > logs/prepare_cosmic.log 2>&1 &
```

**Esto es aún más pesado que ABraOM**: ~117 millones de filas combinadas entre los 4 archivos de entrada, decenas de GB. Esperá que tome de minutos a algunas horas y uso de RAM de doble dígito en GB, dependiendo de la máquina — no es un paso interactivo (ver `docs/AUXILIARY_DATABASES.md`, sección de rendimiento).

Al terminar, completá:
```yaml
cosmic:
  raw_file: db/cosmic/v104/raw_input/cosmic_flat_input.tsv.gz
  processed_db: db/cosmic/v104/GRCh38_to_GRCh38
```

## Paso 3 — confirmar el subtipo tumoral y los metadatos de muestra

`config/cosmic_tumor_type_mapping.tsv` ya trae una fila `HGSOC → ovary/carcinoma` (corregida y verificada contra COSMIC real). Si tu cohorte es carcinoma seroso de alto grado, no hace falta tocar este archivo. Si hay subtipos no serosos (células claras, endometrioide, mucinoso), usá el mismo `cosmic_primary_site=ovary`/`cosmic_histology=carcinoma` con un `input_tumor_type` propio por subtipo.

Creá `data/ovarian_sample_metadata.tsv` con las 7 muestras (el `sample_id` debe coincidir con la columna de muestra del VCF — para Mutect2 de una sola muestra tumoral, es el nombre de la columna de genotipo, no el nombre del archivo):

```tsv
sample_id	tumor_type	assay
ATPBR-056-001-0367-6535-01ABD_tumor_ATPBR-056-001-0367-6535-01ABD	HGSOC	WGS
ATPBR-056-001-0813-4407-01ABD_tumor_ATPBR-056-001-0813-4407-01ABD	HGSOC	WGS
...  (una fila por cada una de las 7 muestras)
```

Para obtener el nombre exacto de columna de cada VCF:
```bash
bcftools view -h mutect2/tumor_ATPBR-056-001-0367-6535-01ABD/*.vcf.gz | grep "^#CHROM" | cut -f10
```

## Paso 4 — filtrar y normalizar cada VCF de Mutect2

Este paso usa el script que agregamos (`scripts/filter_normalize_mutect2.sh`): selecciona solo `PASS`/`germline`/`panel_of_normals`/`germline;panel_of_normals` y normaliza contra la referencia (`bcftools norm`). Ver `docs/FILTERING_STRATEGY.md` para el porqué de esta selección — no es un simple "quedarse con PASS": `tumoronly` ya sabe evaluar `germline`/`panel_of_normals` como evidencia graduada, no como corte binario, siempre que reciban el config del Paso 5.

```bash
export PATH=/home/mmparedese/miniconda3/envs/vcfman/bin:$PATH
mkdir -p filtered

SAMPLE=tumor_ATPBR-056-001-0367-6535-01ABD
./scripts/filter_normalize_mutect2.sh \
  ../mutect2/$SAMPLE/$SAMPLE.mutect2.filtered_VEP.ann.vcf.gz \
  ../referencias/Homo_sapiens_assembly38.fasta \
  filtered/$SAMPLE
```

Produce `filtered/$SAMPLE.PASS.vcf.gz` (input real para `tumoronly`) y `filtered/$SAMPLE.nonPASS.vcf.gz` (conservado, no se usa, no se borra).

## Paso 5 — el `config.yml` final

Sección `hard_filters` — esto es lo que hace que `germline`/`panel_of_normals` compitan de verdad en vez de ser descartados automáticamente (ver `docs/FILTERING_STRATEGY.md`):

```yaml
hard_filters:
  enabled: true
  assay: "WGS"
  caller_filter_accepted_values: ["PASS", "germline", "panel_of_normals", "germline;panel_of_normals"]
  caller_filter_germline_values: ["germline", "germline;panel_of_normals"]

input:
  path: filtered/tumor_ATPBR-056-001-0367-6535-01ABD.PASS.vcf.gz
  format: vcf
  genome_build: GRCh38
  sample_metadata: data/ovarian_sample_metadata.tsv

cancer:
  default_tumor_type: "HGSOC"
```

El resto de las secciones (`technical_filters`, `scoring`, `guideline_classification`, `predictors`, `adaptive_filtering`, ...) quedan como en `config/example.yml` — son los umbrales científicos, no algo a improvisar sin justificación (ver conversación previa sobre Fase 1).

## Paso 6 — correr con UNA muestra

```bash
Rscript exec/tumoronly run --config config/ovarian_cohort.yml --run-id ATPBR_0367
```

Salida en `results/ATPBR_0367/`:
```
tables/variants_all.tsv.gz          # TODAS las variantes + rastro de decisión completo
tables/variants_retained.tsv.gz     # PASS/REVIEW
tables/variants_excluded.tsv.gz     # FAIL (conservadas, no borradas)
tables/filter_audit.tsv.gz
maf/filtered.maf.gz
plots/
report/tumor_only_report.html
manifest.json
config.resolved.yml
```

**Cómo verificar que la adaptación de `germline`/`panel_of_normals` está funcionando** (columnas clave en `variants_all.tsv.gz`): ninguna variante con `FILTER` original `germline`/`panel_of_normals` debería tener `caller_filter_not_pass` en `hard_filter_reason` — si la ves ahí, revisá que `hard_filters.caller_filter_accepted_values` esté bien escrito en el config resuelto (`config.resolved.yml` del run).

```r
x <- read.delim("results/ATPBR_0367/tables/variants_all.tsv.gz")
table(x$hard_filter_pass, grepl("caller_filter_not_pass", x$hard_filter_reason))
```

## Paso 7 — correr las 7 muestras

Primero filtrá y normalizá las 7 (paso 4, repetido por muestra):

```bash
for dir in ../mutect2/tumor_*/; do
  SAMPLE=$(basename "$dir")
  ./scripts/filter_normalize_mutect2.sh \
    "$dir/$SAMPLE.mutect2.filtered_VEP.ann.vcf.gz" \
    ../referencias/Homo_sapiens_assembly38.fasta \
    filtered/$SAMPLE
done
```

Ahora elegí cómo correr `tumoronly run`:

**Opción A — cohorte conjunta (recomendado si querés recurrencia entre muestras):**
`cohort.recurrent_variant_fraction_artifact` — el filtro que trata una variante
recurrente en muchas muestras como probable artefacto — solo es significativo
si las 7 muestras están en la misma tabla dentro de una única corrida (ver
`docs/FILTERING_STRATEGY.md`, "Cohort-wide ingestion"). Para eso, `input.path`
acepta una lista en vez de un solo archivo — no hace falta `bcftools merge`:

```yaml
input:
  path:
    - filtered/ATPBR-056-001-0367-6535-01ABD.PASS.vcf.gz
    - filtered/ATPBR-056-001-0813-4407-01ABD.PASS.vcf.gz
    # ... las 7 rutas
```

```bash
Rscript exec/tumoronly run --config config/ovarian_cohort.yml --run-id ovarian_cohort
```

Cada archivo se ingiere por separado (su propio muestreo de columnas) y se
combina en una sola tabla; `sample_id` sigue siendo el nombre de la columna de
genotipo de cada VCF (columna `SOURCE_FILE` en `variants_all.tsv.gz` para
rastrear de qué archivo vino cada fila).

**Opción B — una corrida por muestra:** si preferís 7 `run_id`/reportes
separados (y no necesitás la recurrencia entre muestras en esta primera
ejecución):

```bash
for dir in ../mutect2/tumor_*/; do
  SAMPLE=$(basename "$dir")
  Rscript exec/tumoronly run --config config/ovarian_cohort.yml \
    --run-id "$SAMPLE" \
    2>&1 | tee logs/run_$SAMPLE.log
done
```

(Necesitás actualizar `input.path` en el config antes de cada corrida, o pasar
un config distinto por muestra.) En este modo `total_samples == 1` en cada
corrida, así que `variant_cohort_freq` no aporta nada — es la contrapartida de
tener resultados separados por muestra.

## Paso 8 (opcional) — OncoKB y reporte

```bash
export ONCOKB_TOKEN="..."
Rscript exec/tumoronly annotate-oncokb --run-dir results/ATPBR_0367
Rscript exec/tumoronly report --run-dir results/ATPBR_0367
```

## Validado en este entorno

Antes de entregar este tutorial corrí el flujo completo con datos reales: extraje una región pequeña (TP53, una de las 7 muestras reales), la filtré y normalicé con `scripts/filter_normalize_mutect2.sh`, y corrí `tumoronly run` de punta a punta. Confirmé que las 3 variantes `germline` de esa región no se descartan por `caller_filter_not_pass` — solo por razones técnicas genuinas (profundidad/QC de muestra, esperable en una muestra de prueba de 4 variantes). El mecanismo funciona; los tiempos de ABraOM/COSMIC completos no se ejecutaron (son de horas, ver Pasos 1-2).

También corrí `tumoronly run` con `input.path` como lista de dos VCFs de muestra distintos (Opción A del Paso 7): las 4 variantes resultantes conservaron el `sample_id` correcto por archivo (columna `SOURCE_FILE` presente), y la variante compartida entre ambos archivos quedó con `variant_cohort_freq = 1.0` mientras que las privadas de cada muestra quedaron en `0.5` — confirma que la ingesta de lista de VCFs alimenta correctamente `compute_cohort_recurrence()` sin necesitar `bcftools merge`.

## Troubleshooting

- **"Could not determine genome build"** → fijá `input.genome_build: GRCh38` explícitamente.
- **`caller_filter_not_pass` aparece en variantes germline/PON** → revisá `config.resolved.yml` del run: `hard_filters.caller_filter_accepted_values` debe incluir esos valores exactos (comparación de string exacto, no por tag individual).
- **`prepare-abraom`/`build-cosmic-input` fallan por checksum tras cambiar un input** → agregá `--force`.
- **manifest con `n_ref_mismatch` alto en ABraOM** → el FASTA no coincide con el build/versión esperado; ver `docs/AUXILIARY_DATABASES.md`.
- **bcftools "not found"** → activá el entorno: `export PATH=/home/mmparedese/miniconda3/envs/vcfman/bin:$PATH`.
