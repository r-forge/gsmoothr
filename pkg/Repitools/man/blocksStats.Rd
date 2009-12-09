\name{blocksStats}
\alias{blocksStats}
\alias{blocksStats.AffymetrixCelSet}
\alias{blocksStats.GenomeDataList}
\alias{blocksStats.default}
\title{A routine to calculate statistics for regions in the genome}
\description{
For each region of interest, this routine interrogates probes for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.  Regions can be surroundings of TSSs, or can be a table of user - specified regions. The function determines if it should calculate around the TSSs or in the user's regions of interest based on whether the variable \code{regionsOfInterestTable} is NULL. \code{regionsOfInterestTable} and \code{geneCoords} cannot both be given to the function.
}
\usage{
\method{blocksStats}{AffymetrixCelSet}(cs, coordinatesTable, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE, useAsRegions=FALSE, ...)
\method{blocksStats}{GenomeDataList}(cs, coordinatesTable, design, upStream = 0, downStream = 2000, verbose=TRUE, useAsRegions = FALSE, seqLen = NULL, total.lib.size = TRUE, ...)
\method{blocksStats}{default}(cs, ndf, coordinatesTable, design, upStream=0, downStream=2000, verbose=TRUE, robust=FALSE, minNRobust=10, adjustMethod="fdr", log2adjust = TRUE, useAsRegions=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, or a data.frame with columns called, presumably from normalized data}
  \item{ndf}{A \code{data.frame} with columns \code{chr}, \code{position}, and \code{index}}
  \item{coordinatesTable}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{name}, \code{start}, \code{end}, \code{strand}.}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{upStream}{number of bases upstream to consider in calculation of statistics}
  \item{downStream}{number of bases upstream to consider in calculation of statistics}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
  \item{log2adjust}{whether to take log2 of array intensities.}
  \item{useAsRegions}{If true, the \code{start} and \code{end} columns of \code{coordinatesTable} will define the region to calculate statistics in.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{total.lib.size}{whether to use the total number of sequencing reads as the library size. If FALSE, the library size will be the number of reads counted within the coordinatesTable blocks.}

  \item{...}{Arguments passed to internal function for statistics calculations}
}
\details{
At present, the approach simply calculates the contrast for each region and tests whether the collection of these summaries are significantly non-zero.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{coordinatesTable}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}
