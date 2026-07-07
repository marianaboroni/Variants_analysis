build_analysis_dashboard <- function(
    variants,
    cfg,
    figure_manifest = NULL,
    sample_qc = NULL,
    tmb_summary = NULL,
    variant_qc_summary = NULL,
    clonality_summary = NULL,
    ancestry_result = NULL) {
  if (!isTRUE(cfg_get(cfg, c("dashboard", "enabled"), TRUE))) {
    return(invisible(NULL))
  }
  outdir <- cfg$output$dir
  dashboard_path <- file.path(outdir, cfg_get(cfg, c("dashboard", "file"), "variant_filtering_dashboard.html"))
  dir.create(dirname(dashboard_path), recursive = TRUE, showWarnings = FALSE)

  ranking <- if ("confidence_rank" %in% names(variants)) confidence_ranking_table(variants) else data.frame()
  class_counts <- class_summary(variants)
  confidence_counts <- if ("confidence_bucket" %in% names(variants)) confidence_summary(variants) else data.frame()
  guideline_counts <- if ("somatic_oncogenicity_class" %in% names(variants) && "germline_acmg_class" %in% names(variants)) guideline_summary(variants) else data.frame()

  keep_n <- sum(variants$recommended_variant_action == "keep_somatic_candidate", na.rm = TRUE)
  remove_n <- sum(variants$recommended_variant_action == "remove_from_somatic_list", na.rm = TRUE)
  review_n <- sum(variants$recommended_variant_action == "manual_review", na.rm = TRUE)

  figure_html <- dashboard_figures(figure_manifest, dashboard_path)
  top_somatic <- dashboard_top_table(ranking, "keep_somatic_candidate", 25)
  top_review <- dashboard_top_table(ranking, "manual_review", 25)
  top_remove <- dashboard_top_table(ranking, "remove_from_somatic_list", 25)

  lines <- c(
    "<!doctype html>",
    "<html lang=\"pt-BR\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Variant Filtering Dashboard</title>",
    dashboard_css(),
    "</head>",
    "<body>",
    "<main>",
    "<header>",
    "<h1>Variant Filtering Dashboard</h1>",
    sprintf("<p>Output: <code>%s</code></p>", html_escape(normalizePath(outdir, mustWork = FALSE))),
    "</header>",
    "<section class=\"cards\">",
    dashboard_card("Total variantes", nrow(variants), "Todas as chamadas após padronização."),
    dashboard_card("Candidatas somáticas", keep_n, "Manter na lista somática."),
    dashboard_card("Remover", remove_n, "Polimorfismo provável ou artefato provável."),
    dashboard_card("Revisão manual", review_n, "Conflito, incerteza ou possível achado germinativo."),
    "</section>",
    "<section>",
    "<h2>Como interpretar</h2>",
    "<p>O ranking ordena variantes por confiança de uso prático: manter como candidata somática, remover da lista somática ou revisar manualmente. Em tumor-only, nenhuma classe substitui normal pareado, PoN compatível ou validação ortogonal.</p>",
    "</section>",
    "<section class=\"grid two\">",
    dashboard_table_section("Classes finais", class_counts),
    dashboard_table_section("Ações recomendadas", confidence_counts),
    "</section>",
    "<section class=\"grid two\">",
    dashboard_table_section("Diretrizes", guideline_counts),
    dashboard_table_section("TMB", head_table(tmb_summary, 20)),
    "</section>",
    "<section>",
    "<h2>Top candidatas somáticas</h2>",
    dashboard_table(top_somatic),
    "</section>",
    "<section>",
    "<h2>Top revisão manual</h2>",
    dashboard_table(top_review),
    "</section>",
    "<section>",
    "<h2>Top remoção da lista somática</h2>",
    dashboard_table(top_remove),
    "</section>",
    "<section>",
    "<h2>Figuras e oncoplots</h2>",
    figure_html,
    "</section>",
    "<section>",
    "<h2>Arquivos principais</h2>",
    "<ul>",
    "<li><code>variants_final.tsv</code>: tabela completa.</li>",
    "<li><code>mutation_confidence_ranking.tsv</code>: ranking somática/germinativa/polimorfismo/revisão.</li>",
    "<li><code>variants_high_confidence_somatic.tsv</code>: candidatas somáticas estritas.</li>",
    "<li><code>variants_uncertain_for_review.tsv</code>: fila de revisão manual.</li>",
    "<li><code>removed_likely_germline.tsv</code> e <code>removed_likely_artifact.tsv</code>: removidas da lista somática.</li>",
    "</ul>",
    "</section>",
    "</main>",
    "</body>",
    "</html>"
  )
  writeLines(lines, dashboard_path)
  dashboard_path
}

dashboard_top_table <- function(ranking, action, n = 25) {
  if (nrow(ranking) == 0 || !"recommended_variant_action" %in% names(ranking)) return(data.frame())
  cols <- c(
    "confidence_rank", "sample_id", "gene", "protein_change", "final_class",
    "confidence_bucket", "somatic_confidence_score", "germline_confidence_score",
    "polymorphism_confidence_score", "artifact_confidence_score",
    "manual_review_priority_score", "max_pop_af", "vaf", "dp", "primary_reason"
  )
  cols <- cols[cols %in% names(ranking)]
  head(ranking[ranking$recommended_variant_action == action, cols, drop = FALSE], n)
}

dashboard_figures <- function(figure_manifest, dashboard_path) {
  if (is.null(figure_manifest) || nrow(figure_manifest) == 0) return("<p>Nenhuma figura registrada.</p>")
  figs <- figure_manifest[figure_manifest$status == "written" & grepl("[.]png$", figure_manifest$file), , drop = FALSE]
  if (nrow(figs) == 0) return("<p>Nenhuma figura PNG registrada.</p>")
  cards <- vapply(seq_len(nrow(figs)), function(i) {
    rel <- rel_path(figs$file[[i]], dirname(dashboard_path))
    sprintf(
      "<figure><img src=\"%s\" alt=\"%s\"><figcaption>%s</figcaption></figure>",
      html_escape(rel),
      html_escape(figs$type[[i]]),
      html_escape(paste(figs$plot_category[[i]], figs$type[[i]], sep = " - "))
    )
  }, character(1))
  paste0("<div class=\"figures\">", paste(cards, collapse = "\n"), "</div>")
}

dashboard_table_section <- function(title, x) {
  paste0("<div><h2>", html_escape(title), "</h2>", dashboard_table(x), "</div>")
}

dashboard_table <- function(x) {
  if (is.null(x) || nrow(x) == 0) return("<p>Sem dados.</p>")
  x <- as.data.frame(x)
  x <- x[, seq_len(min(ncol(x), 12)), drop = FALSE]
  header <- paste(sprintf("<th>%s</th>", html_escape(names(x))), collapse = "")
  rows <- apply(x, 1, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html_escape(format_dashboard_value(row))), collapse = ""), "</tr>")
  })
  paste0("<div class=\"table-wrap\"><table><thead><tr>", header, "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>")
}

head_table <- function(x, n = 20) {
  if (is.null(x)) return(data.frame())
  head(as.data.frame(x), n)
}

format_dashboard_value <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- ""
  out
}

dashboard_card <- function(title, value, note) {
  sprintf(
    "<article class=\"card\"><span>%s</span><strong>%s</strong><p>%s</p></article>",
    html_escape(title),
    html_escape(format(value, big.mark = ",")),
    html_escape(note)
  )
}

dashboard_css <- function() {
  "<style>
  body { margin: 0; font-family: Arial, Helvetica, sans-serif; background: #f6f7f9; color: #172033; }
  main { max-width: 1280px; margin: 0 auto; padding: 32px 24px 56px; }
  header { margin-bottom: 24px; }
  h1 { margin: 0 0 8px; font-size: 34px; }
  h2 { margin: 24px 0 12px; font-size: 20px; }
  p, li { color: #42526b; line-height: 1.45; }
  code { background: #eef1f5; padding: 2px 5px; border-radius: 4px; }
  .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }
  .card { background: white; border: 1px solid #d8e1ec; border-radius: 8px; padding: 16px; }
  .card span { display: block; color: #53657d; font-size: 13px; }
  .card strong { display: block; font-size: 30px; margin: 8px 0; }
  .card p { margin: 0; font-size: 13px; }
  .grid.two { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 20px; }
  .table-wrap { overflow-x: auto; background: white; border: 1px solid #d8e1ec; border-radius: 8px; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  th, td { padding: 8px 10px; border-bottom: 1px solid #e6ecf3; text-align: left; vertical-align: top; }
  th { background: #eef3f8; color: #2e3b4e; position: sticky; top: 0; }
  .figures { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px; }
  figure { margin: 0; background: white; border: 1px solid #d8e1ec; border-radius: 8px; padding: 12px; }
  img { width: 100%; height: auto; display: block; }
  figcaption { color: #53657d; font-size: 12px; margin-top: 8px; }
  @media (max-width: 900px) { .cards, .grid.two, .figures { grid-template-columns: 1fr; } }
  </style>"
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

rel_path <- function(path, start) {
  path <- normalizePath(path, mustWork = FALSE)
  start <- normalizePath(start, mustWork = FALSE)
  if (requireNamespace("xfun", quietly = TRUE)) {
    return(xfun::relative_path(path, start))
  }
  path
}
