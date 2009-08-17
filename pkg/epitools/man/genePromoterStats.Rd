\name{genePromoterStats}
\alias{genePromoterStats}
\alias{genePromoterStats.AffymetrixCelSet}
\alias{genePromoterStats.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to calculate statistics for a given region of the promoter}
\description{
For each gene promoter, this routine interrogates probes over a specified region of a promoter (relative to TSS) for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.
}
\usage{
\method{genePromoterStats}{AffymetrixCelSet}(cs, geneCoords, design, upStream = 0, downStream = 2000, verbose = -20, robust = FALSE, minNRobust = 10, adjustMethod = "fdr", ...)
\method{genePromoterStats}{default}(cs, ndf, geneCoords, design, upStream = 0, downStream = 2000, verbose = -20, robust = FALSE, minNRobust = 10, adjustMethod = "fdr", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, presumably from normalized data}
  \item{ndf}{A \code{data.frame} with columns \code{chr}, \code{position}, and \code{index}}
  \item{geneCoords}{\code{data.frame} giving the gene coordinates. The column names of this are \code{seqname} for the chromosome, \code{probeset_id} for the gene name, \code{start}, \code{stop}, \code{strand}.}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{upStream}{number of bases upstream to consider in calculation of statistics}
  \item{downStream}{number of bases upstream to consider in calculation of statistics}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
  \item{...}{Arguments passed to internal function for statistics calculations}
}
\details{
At present, the approach simply calculates the contrast for each probe and tests whether the collection of these summaries are significantly non-zero.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{geneCoords}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}, \code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}
