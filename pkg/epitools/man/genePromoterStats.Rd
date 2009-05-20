\name{genePromoterStats}
\alias{genePromoterStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to calculate statistics for a given region of the promoter}
\description{
For each gene promoter, this routine interrogates probes over a specified region of a promoter (relative to TSS) for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.
}
\usage{
genePromoterStats(cs, geneCoords, design, ind = NULL, lookupT = NULL, upStream = 0, downStream = 2000, probeSpacing = 35, verbose = -20, robust = FALSE, minNRobust = 10, adjustMethod = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, presumably from normalized data}
  \item{geneCoords}{\code{data.frame} giving the gene coordinates}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{ind}{cell indices into the \code{AffymetrixCelSet} object, if \code{NULL} the lookup table is created}
  \item{lookupT}{lookup table of indices, if \code{NULL} the lookup table is created}
  \item{upStream}{number of bases upstream to consider in calculation of statistics}
  \item{downStream}{number of bases upstream to consider in calculation of statistics}
  \item{probeSpacing}{minimum probe spacing used}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
}
\details{
At present, the approach simply calculates the contrast for each probe and tests whether the collection of these summaries are significantly non-zero.

If you call the routine with \code{ind=NULL} and \code{lookupT=NULL}, a lookup table is built from the information in the \code{geneCoords} column.  The columns \code{seqname}, \code{strand}, \code{start} and \code{stop} will be necessary.  Note that the lookup table and indices vector is (silently) saved to disk.

If this routine is called successively over slightly different regions (e.g. because you expect changes upstream of TSS for some data and downstream for other data), you may wish to re-use the lookup tables and give subset tables each time.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{geneCoords}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}, \code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}
