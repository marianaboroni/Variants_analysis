\name{read_variant_input}
\alias{read_variant_input}
\title{Read a tumor-only variant input (VCF / MAF / TSV) into the canonical schema}
\usage{
read_variant_input(path, format = "auto", sample_metadata = NULL,
                   column_map = NULL, chunk_size = NULL)
}
\arguments{
\item{path}{input file (.vcf/.vcf.gz/.maf/.maf.gz/.tsv/.tsv.gz/.txt/.txt.gz).}
\item{format}{"auto" (default) or one of "vcf"/"maf"/"tsv".}
\item{sample_metadata}{optional path/data.frame with per-sample metadata.}
\item{column_map}{optional named list mapping canonical (lowercase) fields to
source columns (for custom TSVs).}
\item{chunk_size}{optional block size for VCF reading (bounded memory).}
}
\value{
a list with \code{variants} (canonical data.frame with originals preserved),
\code{format}, \code{header_metadata}, \code{schema_report}, \code{column_mapping}.
}
\description{
Single ingestion entry point. Detects the format by content (not extension alone),
interprets \code{#}/\code{##} header lines, maps to the canonical internal schema
(CHROM/POS/REF/ALT/SAMPLE_ID/...), preserves all original columns, and emits a
schema report. Raises an actionable error when required fields cannot be resolved.
}
\examples{
\dontrun{
res <- read_variant_input("data/sample.vep.vcf.gz")
head(res$variants); res$schema_report
}
}
