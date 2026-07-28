\name{run_tumoronly}
\alias{run_tumoronly}
\title{Run tumoronly with the v2 user-facing API}
\usage{
run_tumoronly(input = NULL, config, output_dir = NULL, run_id = NULL,
              strict = TRUE, validate = TRUE)
}
\arguments{
\item{input}{optional path to a VCF, MAF, TSV, or TXT variant file. When set,
overrides \code{input.path} in \code{config}.}
\item{config}{path to a YAML configuration file, or an already-read list.}
\item{output_dir}{optional output directory. When set, overrides
\code{analysis.output_dir}.}
\item{run_id}{optional explicit run identifier.}
\item{strict}{logical; if \code{TRUE}, validation errors stop execution before
the scientific workflow is run.}
\item{validate}{logical; if \code{TRUE}, run \code{\link{validate_tumoronly_input}}
before analysis.}
}
\value{
A list returned by \code{\link{run_tumor_only}} with an additional
\code{validation} element when validation is enabled.
}
\description{
Thin v2 API wrapper around the audited \code{\link{run_tumor_only}} workflow.
It applies user-facing input/output overrides, optionally validates the input,
then delegates all scientific classification logic to the existing pipeline.
The wrapper does not duplicate classification rules.
}
\examples{
\dontrun{
result <- run_tumoronly(
  input = "variants.tsv",
  config = "config.yaml",
  output_dir = "results"
)
}
}
