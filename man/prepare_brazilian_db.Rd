\name{prepare_brazilian_db}
\alias{prepare_brazilian_db}
\title{Prepare a Brazilian population-frequency database}
\usage{
prepare_brazilian_db(config, force = FALSE)
}
\arguments{
\item{config}{path to a YAML configuration file, or an already-read list.}
\item{force}{logical; rebuild even if an existing manifest matches the inputs.}
}
\value{
A list with \code{db_path}, \code{manifest_path}, \code{manifest}, and
\code{reused}.
}
\description{
Prepares a VCF-standard Brazilian population-frequency table from a user-supplied
ANNOVAR-style source and reference FASTA. The resulting database can be used as
separate Brazilian-population evidence during tumor-only classification.
}
\examples{
\dontrun{
prepare_brazilian_db("config/example.yml")
}
}
