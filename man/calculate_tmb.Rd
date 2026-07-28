\name{calculate_tmb}
\alias{calculate_tmb}
\title{Calculate tumor mutational burden and variant burden}
\usage{
calculate_tmb(variants, cfg = NULL)
}
\arguments{
\item{variants}{Classified tumoronly variant table.}
\item{cfg}{Resolved tumoronly configuration.}
}
\value{
A list with updated \code{variants} and per-sample \code{summary}.
}
\description{
Calculates countable somatic-candidate burden and, when \code{tmb.callable_mb}
is explicitly configured, exploratory tumor mutational burden in mutations/Mb.
Without a valid callable territory, TMB is marked not evaluable.
}

