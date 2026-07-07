\name{check_installation}
\alias{check_installation}
\title{Check the tumoronly installation and report capabilities}
\usage{
check_installation(quiet = FALSE)
}
\arguments{
\item{quiet}{suppress console output.}
}
\value{
(invisibly) a list of \code{packages}, \code{tools}, \code{capabilities} and \code{ok}.
}
\description{
Reports which required/optional R packages and external tools are present and what
that enables (cross-build COSMIC liftover, REF validation, OncoKB API, oncoplots,
reporting). Never fails; backs the \code{tumoronly doctor} command.
}
\examples{
check_installation()
}
