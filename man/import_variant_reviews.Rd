\name{import_variant_reviews}
\alias{import_variant_reviews}
\title{Import human variant reviews into the append-only evidence database}
\usage{
import_variant_reviews(file, db_dir = "db/variant_evidence", reviewer = NULL)
}
\arguments{
\item{file}{review TSV (from \code{export-review-template}, filled in).}
\item{db_dir}{evidence database directory (default \code{db/variant_evidence}).}
\item{reviewer}{optional default reviewer if absent in the file.}
}
\value{
(invisibly) a summary list (counts by label, number appended, history path).
}
\description{
Confirmed labels (TRUE_POSITIVE / FALSE_POSITIVE) come ONLY from the review file;
the pipeline never converts PASS/FAIL into labels. Storage is append-only:
history is never overwritten, and \code{reviewed_variants} keeps the latest review
per variant. TRUE_POSITIVE/FALSE_POSITIVE require a reviewer (human/independent
evidence).
}
\examples{
\dontrun{ import_variant_reviews("reviews.tsv") }
}
