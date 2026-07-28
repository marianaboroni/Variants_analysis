\name{add_clonality_estimates}
\alias{add_clonality_estimates}
\title{Add clonality estimates}
\usage{
add_clonality_estimates(variants, cfg = NULL)
}
\arguments{
\item{variants}{Classified tumoronly variant table.}
\item{cfg}{Resolved tumoronly configuration.}
}
\value{
A list with updated \code{variants} and per-sample \code{summary}.
}
\description{
Adds CCF-based clonality estimates when VAF and purity are available. If purity
is missing, variants may receive explicitly labelled VAF-proxy classes. The
module does not alter filtering or final tumor-only classification.
}

