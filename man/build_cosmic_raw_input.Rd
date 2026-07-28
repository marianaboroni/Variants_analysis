\name{build_cosmic_raw_input}
\alias{build_cosmic_raw_input}
\title{Build a flat COSMIC input table from raw COSMIC files}
\usage{
build_cosmic_raw_input(config, force = FALSE)
}
\arguments{
\item{config}{path to a YAML configuration file, or an already-read list.}
\item{force}{logical; rebuild even if an existing manifest matches the inputs.}
}
\value{
A list with \code{raw_file_path}, \code{manifest_path}, \code{manifest}, and
\code{reused}.
}
\description{
Combines COSMIC mutation-coordinate files, normal VCF coordinate files, and
classification metadata into a flat, checksummed input suitable for
\code{\link{prepare_cosmic_db}}. COSMIC data are user-provided and are not
redistributed by the package.
}
\examples{
\dontrun{
build_cosmic_raw_input("config/example.yml")
}
}
