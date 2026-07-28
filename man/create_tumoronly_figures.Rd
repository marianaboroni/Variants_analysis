\name{create_tumoronly_figures}
\alias{create_tumoronly_figures}
\title{Create v2 QC and classification figures}
\usage{
create_tumoronly_figures(variants, audit, figures_dir, figure_data_dir,
                         cfg = NULL)
}
\arguments{
\item{variants}{classified variant data.frame with v2 evidence columns.}
\item{audit}{filter audit data.frame.}
\item{figures_dir}{directory where PDF, SVG, and PNG figures are written.}
\item{figure_data_dir}{directory where per-figure TSV data files are written.}
\item{cfg}{resolved configuration list.}
}
\value{
A manifest data.frame with one row per attempted figure.
}
\description{
Generates the first v2 figure set: filtering workflow, QC and score
distributions, classification counts, VAF by depth, class burden by sample, and
missing-evidence counts. Each successful figure is exported as PDF, SVG, and
300 dpi PNG, and its source data are saved as TSV. SVG has a data-driven vector
fallback for systems without a working Cairo/SVG device.
}
\examples{
\dontrun{
figures <- create_tumoronly_figures(variants, audit, "figures", "figure_data")
}
}
