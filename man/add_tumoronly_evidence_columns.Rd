\name{add_tumoronly_evidence_columns}
\alias{add_tumoronly_evidence_columns}
\title{Add v2 traceability columns to classified variants}
\usage{
add_tumoronly_evidence_columns(variants, cfg = NULL)
}
\arguments{
\item{variants}{classified variant data.frame after the audited tumoronly
workflow has calculated \code{final_class}.}
\item{cfg}{resolved configuration list.}
}
\value{
The input data.frame with \code{evidence_supporting_classification},
\code{evidence_against_classification}, \code{missing_evidence}, and
\code{classification_explanation}.
}
\description{
Materializes the evidence trail required by the v2 output contract. This is an
output layer only: it does not change \code{final_class}, \code{filter_status},
scores, or thresholds.
}
\examples{
\dontrun{
variants <- add_tumoronly_evidence_columns(variants, cfg)
}
}
