\name{blocksStats}
\alias{blocksStats}
\alias{blocksStats.AffymetrixCelSet}
\alias{blocksStats.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to calculate statistics for regions in the genome}
\description{
For each region of interest, this routine interrogates probes for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.  Regions can be surroundings of TSSs, or can be a table of user - specified regions. The function determines if it should calculate around the TSSs or in the user's regions of interest based on whether the variable \code{regionsOfInterestTable} is NULL.
}
\usage{
\method{blocksStats}{AffymetrixCelSet}(cs, geneCoords, design, upStream=0, downStream=2000, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE, regionsOfInterestTable=NULL, ...)
\method{blocksStats}{default}(cs, ndf, geneCoords, design, upStream=0, downStream=2000, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE,   regionsOfInterestTable=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, or a data.frame with columns called, presumably from normalized data}
  \item{ndf}{A \code{data.frame} with columns \code{chr}, \code{position}, and \code{index}}
  \item{geneCoords}{\code{data.frame} giving the gene coordinates. The column names of this are \code{seqname} for the chromosome, \code{probeset_id}, \code{start}, \code{stop}, \code{strand}.}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{upStream}{number of bases upstream to consider in calculation of statistics}
  \item{downStream}{number of bases upstream to consider in calculation of statistics}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{regionsOfInterestTable}{A \code{data.frame} with elements \code{chr} (character), \code{start} (integer), \code{end} (integer), \code{strand} ("+" or "-"), and \code{name}(character)}
  \item{...}{Arguments passed to internal function for statistics calculations}
}
\details{
At present, the approach simply calculates the contrast for each probe and tests whether the collection of these summaries are significantly non-zero.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{regionsOfInterestTable}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}