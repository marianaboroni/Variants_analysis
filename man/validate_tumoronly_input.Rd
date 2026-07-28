\name{validate_tumoronly_input}
\alias{validate_tumoronly_input}
\title{Validate a tumoronly input file}
\usage{
validate_tumoronly_input(input = NULL, config = NULL, strict = TRUE,
                         output = NULL)
}
\arguments{
\item{input}{optional input file path. Overrides \code{input.path} in
\code{config}.}
\item{config}{optional YAML path or configuration list. If omitted, v2 defaults
are used and an input path must be supplied.}
\item{strict}{logical; strict mode records missing indispensable evidence as
errors instead of permissive warnings where appropriate.}
\item{output}{optional report path. Supported extensions are \code{.json},
\code{.tsv}, \code{.txt}, \code{.html}, and \code{.htm}.}
}
\value{
A \code{tumoronly_validation} list with \code{ok}, \code{issues},
\code{summary}, \code{schema_report}, and \code{column_mapping}.
}
\description{
Reads the configured input in read-only mode, resolves the genome build,
performs canonical schema mapping, derives core metrics where possible, and
returns actionable errors or warnings before a full run is attempted.
}
\examples{
\dontrun{
validation <- validate_tumoronly_input(
  input = "variants.tsv",
  config = "config.yaml",
  output = "validation.html"
)
validation$ok
}
}
