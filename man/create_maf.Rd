\name{create_maf}
\alias{create_maf}
\title{Create a maftools-compatible MAF from a filtered variant table}
\usage{
create_maf(variants, path = NULL)
}
\arguments{
\item{variants}{variant data.frame with chrom/pos/ref/alt/gene/consequence and,
optionally, vaf/dp/alt_count/filter_status/COSMIC_*/confidence/OncoKB columns.}
\item{path}{optional output path (\code{.maf}/\code{.maf.gz}); if given, the MAF is
written and validated with \code{maftools::read.maf()}.}
}
\value{
the MAF data.frame (invisibly if written).
}
\description{
Maps VEP consequences to MAF \code{Variant_Classification} using an explicit,
tested table (unknown consequences fall into a controlled bucket and raise a QC
warning), derives \code{Variant_Type} from REF/ALT, and assembles the required
MAF fields plus optional depth/VAF/status/COSMIC/confidence/OncoKB columns.
}
\examples{
v <- data.frame(chrom = "17", pos = 7674220, ref = "C", alt = "T", gene = "TP53",
                consequence = "missense_variant", sample_id = "S1", vaf = 0.4)
create_maf(v)
}
