\name{render_tumor_only_report}
\alias{render_tumor_only_report}
\title{Render the tumor-only HTML report for a run directory}
\usage{
render_tumor_only_report(run_dir)
}
\arguments{
\item{run_dir}{a completed run directory (from \code{\link{run_tumor_only}}).}
}
\value{
path to the rendered HTML (invisibly).
}
\description{
Prefers the R Markdown template in \code{inst/report/}; if rmarkdown or pandoc is
unavailable it falls back to a self-contained HTML summary built directly from the
run tables, so a report is always produced.
}
\examples{
\dontrun{
render_tumor_only_report("results/example_run")
}
}
