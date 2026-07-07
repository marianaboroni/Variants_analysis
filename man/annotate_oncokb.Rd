\name{annotate_oncokb}
\alias{annotate_oncokb}
\title{Annotate a completed run with OncoKB (optional, non-filtering)}
\usage{
annotate_oncokb(run_dir, scope = c("retained", "all"), cfg = NULL)
}
\arguments{
\item{run_dir}{a completed run directory (from \code{\link{run_tumor_only}}).}
\item{scope}{"retained" (default) or "all".}
\item{cfg}{optional resolved config; read from the run's \code{config.resolved.yml}
when omitted.}
}
\value{
the annotation data.frame (invisibly).
}
\description{
Reads the frozen \code{tables/variants_all.tsv.gz}, annotates the selected scope
against OncoKB, and writes \code{tables/oncokb_annotations.tsv.gz}. It never
modifies \code{filter_status} or the primary run results, and completes even when
the API is unavailable (marking variants as not annotated). The API token is read
only from the \code{ONCOKB_TOKEN} environment variable and is never logged or
serialized.
}
\examples{
\dontrun{
Sys.setenv(ONCOKB_TOKEN = "...")
annotate_oncokb("results/example_run")
}
}
