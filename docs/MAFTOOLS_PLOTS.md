# maftools visualizations

All mutation plots are produced by **native maftools** from a single validated MAF
object (`maftools::read.maf(..., removeSilent = FALSE)` with optional
`clinicalData`). ggplot2 is not used to imitate any maftools plot; the previous
hand-rolled ggplot2 oncoplot was removed.

`R/plots_maftools.R::create_maftools_plots()` produces (each as **PNG + PDF**):

| Plot | maftools function | Minimum data rule (config) |
|---|---|---|
| Summary | `plotmafSummary` | always |
| Oncoplot | `oncoplot` | `plots.minimum_samples_for_oncoplot` (2) |
| Ti/Tv | `titv` + `plotTiTv` | any SNVs |
| Lollipop (top genes) | `lollipopPlot` | `plots.minimum_protein_changes_for_lollipop` (3) |
| Rainfall (top sample) | `rainfallPlot` | `plots.minimum_mutations_for_rainfall` (50) |
| Somatic interactions | `somaticInteractions` | `plots.minimum_samples_for_interactions` (20) |

When a rule is not met, no plot is drawn and the reason is recorded in
`plots/plot_manifest.tsv` (e.g. *"need >= 20 samples (have 1)"*) and surfaced as a
report alert — never a fake/empty plot. Figure dimensions scale with sample and
gene counts. Empty/failed device output is discarded.

Config:

```yaml
plots:
  oncoplot: { top: 20, sort_by_annotation: true, remove_non_mutated: false,
              clinical_features: [tumor_type, tumor_subtype, histology] }
  minimum_samples_for_oncoplot: 2
  minimum_samples_for_interactions: 20
  minimum_mutations_for_rainfall: 50
  minimum_protein_changes_for_lollipop: 3
```

ggplot2 is used **only** for non-maftools panels (filter funnel, predictor
coverage, exclusion reasons, score distributions, QC) — never for oncoplot/Ti-Tv/
rainfall/lollipop.
