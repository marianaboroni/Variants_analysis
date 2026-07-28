# Troubleshooting

## `E_INPUT_MISSING`

No input file was provided. Set `input.path` in the YAML or pass `--input`.

## `E_INPUT_NOT_FOUND`

The configured file does not exist from the current working directory. Use an
absolute path or run from the repository/project root.

## `E_GENOME_BUILD`

The genome build could not be resolved or is not supported. Set:

```yaml
input:
  genome_build: GRCh38
```

Accepted values are `GRCh37` and `GRCh38`.

## `E_STANDARDIZE`

The file was read, but canonical fields could not be standardized. Check
`validation.html` and add `input.column_map`.

## `E_DEPTH_MISSING`, `E_ALT_COUNT_MISSING`, `E_VAF_MISSING`

Strict validation could not find indispensable technical evidence. Add the
columns, fix the mapping, or use permissive validation only when you understand
the resulting limitations.

## `W_COHORT_RECURRENCE_NOT_EVALUABLE`

The run has fewer samples than `cohort.min_cohort_size_for_recurrence`. This is
expected for single-sample runs.

## `W_TLOD_MISSING`

The pipeline can still run, but caller-specific quality evidence is incomplete.
The missing evidence is propagated to per-variant explanations.

## `W_COSMIC_NOT_CONFIGURED`

COSMIC evidence is unavailable. This does not invalidate core filtering, but
driver/context interpretation is limited.

## `W_TMB_CALLABLE_MB_MISSING`

The TMB/burden module ran, but `tmb.callable_mb` was missing or invalid. The
output reports countable candidate burden only. Set:

```yaml
tmb:
  callable_mb: 30
```

Use the true callable territory for your assay.

## `W_CLONALITY_VAF_PROXY_ONLY`

Purity was missing for at least one somatic candidate, so clonality labels were
based on VAF proxy rather than CCF. Add tumor purity and, when available, local
copy number and mutation multiplicity.

## `W_ML_NO_ACTIVE_MODEL`

ML is enabled, but no activated model exists in `ml.db_dir`. Import reviewed
labels, build a training set, train a candidate model, and activate it
explicitly before expecting predictions.

## `W_ANCESTRY_PANEL_MISSING`

Ancestry is enabled, but no AIMs marker panel was found at
`ancestry.marker_panel.path`. Provide a versioned panel or disable ancestry for
the run.

## Report did not render with R Markdown

The pipeline falls back to a self-contained HTML report if `rmarkdown` or pandoc
is unavailable. Install `rmarkdown` and pandoc for the template renderer.

## SVG warnings on macOS

Some R installations lack a working Cairo/X11 SVG device. The v2 figure writer
uses a data-driven vector SVG fallback, so SVG files are still produced.
