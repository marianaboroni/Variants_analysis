\name{create_maftools_plots}
\alias{create_maftools_plots}
\title{Generate the native maftools plot set for a run}
\usage{
create_maftools_plots(maf_path, out_dir, cfg = NULL, sample_metadata = NULL)
}
\arguments{
\item{maf_path}{path to the validated MAF.}
\item{out_dir}{output directory for plots.}
\item{cfg}{resolved config (reads \code{plots:} minimum-data rules and oncoplot options).}
\item{sample_metadata}{optional clinical data (path or data.frame).}
}
\value{
a data.frame manifest with columns \code{plot}, \code{status}, \code{files}, \code{reason}.
}
\description{
Builds one validated maftools MAF object and produces native plots
(plotmafSummary, oncoplot, Ti/Tv, lollipop, rainfall, somaticInteractions), each
as PNG and PDF, honoring minimum-data rules. Plots that are not applicable are
skipped with a recorded reason (never a fake plot). ggplot2 is not used to imitate
maftools plots.
}
\examples{
\dontrun{ create_maftools_plots("results/run/maf/filtered.maf.gz", "results/run/plots") }
}
