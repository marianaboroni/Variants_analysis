\name{apply_ml_predictions}
\alias{apply_ml_predictions}
\title{Apply an activated ML model}
\usage{
apply_ml_predictions(variants, cfg = NULL)
}
\arguments{
\item{variants}{Classified tumoronly variant table.}
\item{cfg}{Resolved tumoronly configuration.}
}
\value{
A list with updated \code{variants} and one-row \code{status} table.
}
\description{
Applies an activated model from the configured variant-evidence registry as an
auxiliary prioritization layer. Predictions never alter \code{filter_status},
\code{final_class}, or scientific scores.
}

