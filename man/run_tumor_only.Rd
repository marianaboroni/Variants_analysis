\name{run_tumor_only}
\alias{run_tumor_only}
\title{Run the tumor-only variant analysis}
\usage{
run_tumor_only(config, run_id = NULL)
}
\arguments{
\item{config}{path to a YAML config, or an already-read config list.}
\item{run_id}{optional explicit run id; otherwise taken from the config, else
derived deterministically from the input checksum.}
}
\value{
(invisibly) a list with the annotated \code{variants}, the resolved
\code{run_dir}, the \code{manifest}, and the filter \code{audit}.
}
\description{
Reads a VEP-annotated VCF (or annotated TSV), resolves and validates the genome
build, applies the documented tumor-only filters, annotates variants against a
locally-prepared COSMIC database by canonical key, classifies them, writes the
standardized run directory (tables, MAF, plots, manifest, report) and returns a
structured result. OncoKB is not part of this flow (see \code{\link{annotate_oncokb}}).
}
\examples{
\dontrun{
run_tumor_only("config/example.yml")
}
}
