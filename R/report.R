# Render the explanatory HTML report for a completed run.

#' Render the tumor-only HTML report for a run directory.
#'
#' Prefers the Quarto/R Markdown template in `inst/report/`. If rmarkdown or
#' pandoc is unavailable it falls back to a self-contained HTML summary built
#' directly from the run tables, so a report is always produced.
#'
#' @param run_dir a completed run directory (from [run_tumor_only()]).
#' @return path to the rendered HTML (invisibly).
#' @export
#' @examples
#' \dontrun{ render_tumor_only_report("results/example_run") }
render_tumor_only_report <- function(run_dir) {
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  out_html <- file.path(run_dir, "report", "tumor_only_report.html")
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  rmd <- report_template_path()

  can_rmd <- requireNamespace("rmarkdown", quietly = TRUE) &&
    !is.null(rmd) && file.exists(rmd) && rmarkdown::pandoc_available()
  if (can_rmd) {
    ok <- tryCatch({
      rmarkdown::render(rmd, output_file = basename(out_html),
                        output_dir = dirname(out_html),
                        params = list(run_dir = run_dir),
                        quiet = TRUE, envir = new.env())
      TRUE
    }, error = function(e) { log_step("report", "rmarkdown render failed; using fallback",
                                      error = conditionMessage(e)); FALSE })
    if (ok) { log_step("report", "HTML report rendered", file = out_html); return(invisible(out_html)) }
  }
  render_fallback_report(run_dir, out_html)
}

#' @keywords internal
report_template_path <- function() {
  dev_roots <- c()
  r_dir <- Sys.getenv("TUMORONLY_R_DIR", unset = "")
  if (nzchar(r_dir)) dev_roots <- c(dev_roots, dirname(r_dir))
  dev_roots <- unique(c(dev_roots, getwd()))
  for (root in dev_roots) {
    dev <- file.path(root, "inst", "report", "tumor_only_report.Rmd")
    if (file.exists(dev)) return(dev)
  }
  p <- system.file("report", "tumor_only_report.Rmd", package = "tumoronly")
  if (nzchar(p) && file.exists(p)) return(p)
  NULL
}

#' Self-contained HTML summary built from the run tables (no pandoc needed).
#' @keywords internal
render_fallback_report <- function(run_dir, out_html) {
  read_tbl <- function(name) {
    p <- file.path(run_dir, "tables", name)
    if (file.exists(p)) read_variants(p, "\t") else NULL
  }
  manifest <- tryCatch(jsonlite::read_json(file.path(run_dir, "manifest.json"), simplifyVector = TRUE),
                       error = function(e) list())
  audit <- read_tbl("filter_audit.tsv.gz")
  retained <- read_tbl("variants_retained.tsv.gz")
  esc <- function(s) gsub("<", "&lt;", gsub("&", "&amp;", as.character(s)))
  tbl_html <- function(df, max = 200) {
    if (is.null(df) || nrow(df) == 0) return("<p><em>none</em></p>")
    df <- head(df, max)
    hdr <- paste0("<tr>", paste0("<th>", esc(names(df)), "</th>", collapse = ""), "</tr>")
    rows <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", esc(r), "</td>", collapse = ""), "</tr>"))
    paste0("<table border=1 cellspacing=0 cellpadding=3>", hdr, paste(rows, collapse = ""), "</table>")
  }
  show_cols <- intersect(c("sample_id", "gene", "protein_change", "consequence", "vaf", "dp",
    "max_pop_af", "COSMIC_EVIDENCE_SUMMARY", "confidence_category", "final_class",
    "driver_class", "filter_reasons"), names(retained %||% data.frame()))
  html <- c(
    "<!DOCTYPE html><html><head><meta charset='utf-8'><title>Tumor-only report</title>",
    "<style>body{font-family:system-ui,Arial,sans-serif;margin:2rem;max-width:1100px}",
    "table{border-collapse:collapse;font-size:12px}th{background:#f0f0f0}h2{border-bottom:2px solid #ccc}</style></head><body>",
    "<h1>Tumor-only variant analysis report</h1>",
    sprintf("<p><b>Run:</b> %s &nbsp; <b>Build:</b> %s &nbsp; <b>Package:</b> %s &nbsp; <b>Date:</b> %s</p>",
            esc(manifest$run_id %||% "?"), esc(manifest$genome_build %||% "?"),
            esc(manifest$package_version %||% "?"), esc(manifest$date %||% "?")),
    sprintf("<p><b>Input variants:</b> %s &nbsp; <b>COSMIC release:</b> %s &nbsp; <b>COSMIC matches:</b> %s</p>",
            esc(manifest$n_input_variants %||% "?"), esc(manifest$cosmic_release %||% "NA"),
            esc(manifest$n_cosmic_matches %||% "?")),
    "<h2>Filtering funnel / audit</h2>", tbl_html(audit),
    "<h2>Retained variants</h2>",
    if (is.null(retained) || nrow(retained) == 0) "<p>No retained variants.</p>"
    else tbl_html(if (length(show_cols)) retained[, show_cols, drop = FALSE] else retained),
    "<h2>Evid\u00eancia COSMIC no contexto do tipo tumoral</h2>",
    local({
      cm <- read_tbl("cosmic_matches.tsv.gz")
      if (is.null(cm) || nrow(cm) == 0 || !"COSMIC_TUMOR_CONTEXT_STATUS" %in% names(cm))
        return("<p><em>No COSMIC matches.</em></p>")
      st <- as.data.frame(table(cm$COSMIC_TUMOR_CONTEXT_STATUS))
      names(st) <- c("tumor_context_status", "n_variants")
      cc <- intersect(c("sample_id","gene","tumor_type","COSMIC_TUMOR_CONTEXT_STATUS",
        "COSMIC_MATCH_LEVEL","COSMIC_CONTEXT_EVALUABLE","COSMIC_MATCHING_TUMOR_OCCURRENCES",
        "COSMIC_OTHER_TUMOR_OCCURRENCES","COSMIC_CONTEXT_SUPPORT_SCORE",
        "COSMIC_TUMOR_CONTEXT_REASON"), names(cm))
      paste0(tbl_html(st), tbl_html(cm[, cc, drop = FALSE]),
        "<p><em>Bias note: COSMIC representation is uneven across tumor types; absence in a tumor type does not exclude relevance; cross-tumor counts are not directly comparable.</em></p>")
    }),
    "<h2>Limitations</h2><ul>",
    "<li>Tumor-only analysis cannot definitively separate somatic from germline variants.</li>",
    "<li>Presence in COSMIC is supporting evidence of recurrence, not proof of somaticity; absence does not exclude relevance.</li>",
    "<li>OncoKB (if run) is a confirmatory annotation only and never affects filtering.</li>",
    "<li>Results require review by qualified professionals; this is not a validated clinical test.</li></ul>",
    "</body></html>")
  writeLines(html, out_html)
  log_step("report", "HTML report written (fallback builder)", file = out_html)
  invisible(out_html)
}
