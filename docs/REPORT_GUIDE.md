# Report guide

`inst/report/tumor_only_report.Rmd` → `report/tumor_only_report.html`
(self-contained). Designed so a reviewer sees the key findings without opening
multiple tables.

## Structure

1. **Executive summary** (first screen): colored cards — input variants, retained,
   excluded, manual review, high priority, in hotspots, computational support,
   COSMIC same-tumor — plus run/build/format/OncoKB status.
2. **Principais achados**: auto-generated bullet findings.
3. **Alerts**: build unconfirmed, missing fields, low predictor coverage, no tumor
   metadata, COSMIC without context, VCF without CSQ, MAF/TSV mapping issues,
   over-exclusion, no retained variants, plots skipped, OncoKB unavailable, …
4. **Top prioritized variants**: interactive `DT` table (search/sort/filter,
   tooltips, CSV export) — the single consolidated review table.
5. **Details** (tabset): Data quality · Filtering funnel · Functional predictors ·
   COSMIC evidence (global vs contextual) · maftools visualizations · Excluded
   variants · Provenance & parameters.
6. **Limitations**.

Each analytical section carries a short *"What this shows / How to interpret"*
note. Full tables remain downloadable under `tables/` (gzipped TSVs). Friendly
labels are used in prose while technical column names remain in the tables and
tooltips.

If `rmarkdown`/pandoc or `DT` is unavailable, `render_tumor_only_report()` falls
back to a self-contained HTML built directly from the run tables (no interactivity
but the same executive summary, funnel, retained/COSMIC sections and limitations).
