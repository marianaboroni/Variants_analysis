\name{prepare_cosmic_db}
\alias{prepare_cosmic_db}
\title{Prepare (or reuse) the local processed COSMIC database}
\usage{
prepare_cosmic_db(config, force = FALSE)
}
\arguments{
\item{config}{path to a YAML config, or an already-read config list.}
\item{force}{rebuild even if a valid cached database exists.}
}
\value{
a list with \code{db_path}, \code{manifest_path}, \code{manifest} and \code{reused}.
}
\description{
Reads the raw COSMIC TSV, validates its build, lifts it over to the target build
if needed (via a pluggable backend), normalizes chromosome names and alleles,
validates REF against the target FASTA when available, builds the canonical key
\code{BUILD|CHROM|POS|REF|ALT}, deduplicates, and writes the processed database
plus a checksummed manifest and provenance tables. Idempotent: a valid,
checksum-matching database is reused without re-running the liftover unless
\code{force = TRUE}. No variant is dropped silently.
}
\examples{
\dontrun{
prepare_cosmic_db("config/example.yml")
}
}
