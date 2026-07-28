\name{infer_ancestry}
\alias{infer_ancestry}
\title{Infer experimental genetic-ancestry proportions}
\usage{
infer_ancestry(variants, cfg = NULL, build = "GRCh38")
}
\arguments{
\item{variants}{Standardized variant table.}
\item{cfg}{Resolved tumoronly configuration with an \code{ancestry} section.}
\item{build}{Genome build, usually \code{GRCh37} or \code{GRCh38}.}
}
\value{
A list with \code{summary}, \code{snps}, \code{qc}, and \code{projection}
tables.
}
\description{
Runs the experimental ancestry module using a user-supplied AIMs panel. Genetic
ancestry is not race or ethnicity, is not diagnostic, and never filters variants
or changes biological classification.
}

