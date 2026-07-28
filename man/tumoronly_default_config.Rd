\name{tumoronly_default_config}
\alias{tumoronly_default_config}
\title{Return the centralized tumoronly v2 defaults}
\usage{
tumoronly_default_config()
}
\value{
A nested list containing input, filtering, classification, report, and plotting
defaults.
}
\description{
Returns the default configuration used by \code{tumoronly init}. Scientific
thresholds are centralized here so user-facing wrappers do not hide defaults in
the CLI or report code.
}
\examples{
cfg <- tumoronly_default_config()
names(cfg)
}
