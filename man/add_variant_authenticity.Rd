\name{add_variant_authenticity}
\alias{add_variant_authenticity}
\title{Add technical-authenticity, biological-support and review scores}
\usage{
add_variant_authenticity(x, cfg = NULL)
}
\arguments{
\item{x}{variant data.frame after scoring/cosmic/driver/predictor layers.}
\item{cfg}{resolved config (reads \code{authenticity:} thresholds).}
}
\value{
\code{x} with \code{ARTIFACT_EVIDENCE_SCORE}, \code{TECHNICAL_AUTHENTICITY_SCORE},
\code{BIOLOGICAL_SUPPORT_SCORE}, \code{VARIANT_REVIEW_SCORE},
\code{VARIANT_AUTHENTICITY_SCORE} (alias of the technical score),
\code{VARIANT_AUTHENTICITY_CATEGORY} (LIKELY_TRUE/UNCERTAIN/LIKELY_ARTIFACT) and
\code{VARIANT_AUTHENTICITY_COMPONENTS}.
}
\description{
Keeps two distinct concepts separate: technical authenticity (is the call real
vs an artifact -- technical/quality signals only) and biological support
(relevance -- consequence/predictors/hotspot/driver/COSMIC/clinical only). The
authenticity category depends ONLY on technical evidence, so strong biological
evidence can never mask a strong technical artifact. None of these scores set
\code{filter_status}. Weights/thresholds are heuristic and configurable.
}
\examples{
add_variant_authenticity(data.frame(technical_evidence_score = 0.9, artifact_score = 0.1))
}
