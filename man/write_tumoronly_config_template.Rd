\name{write_tumoronly_config_template}
\alias{write_tumoronly_config_template}
\title{Write a tumoronly v2 configuration template}
\usage{
write_tumoronly_config_template(output, force = FALSE)
}
\arguments{
\item{output}{path to the YAML file to create.}
\item{force}{logical; overwrite an existing file only when \code{TRUE}.}
}
\value{
The normalized output path, invisibly.
}
\description{
Writes the centralized v2 default configuration to YAML. The function refuses
to overwrite an existing file unless \code{force = TRUE}, matching
\code{tumoronly init --force}.
}
\examples{
\dontrun{
write_tumoronly_config_template("config.yaml")
}
}
