\name{geneBlocksStats}
\alias{geneBlocksStats}
\alias{geneBlocksStats.AffymetrixCelSet}
\alias{geneBlocksStats.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A routine to calculate statistics for a given regions in the genome}
\description{
For each region of interested explicitly specified, this routine interrogates probes for either a high level of absolute signal or a change in signal for some specified contrast (i.e. difference) of interest.
}
\usage{
\method{geneBlocksStats}{AffymetrixCelSet}(cs, design, verbose = -20, robust = FALSE, minNRobust = 10, adjustMethod = "fdr", regionsOfInterestTable, ...)
\method{geneBlocksStats}{default}(cs, ndf, design, verbose = -20, robust = FALSE, minNRobust = 10, adjustMethod = "fdr", regionsOfInterestTable, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{An \code{AffymetrixCelSet} object, or a data.frame with columns called, presumably from normalized data}
  \item{ndf}{A \code{data.frame} with columns \code{chr}, \code{position}, and \code{index}}
  \item{design}{design matrix specifying the contrast to compute (i.e. the samples to use and what differences to take)}
  \item{verbose}{logical, whether to output commments of the processing}
  \item{robust}{logical, whether to use a robust test when fitting}
  \item{minNRobust}{minimum number of probes to use}
  \item{adjustMethod}{method used to p-values for multiple testing}
  \item{regionsOfInterestTable}{A \code{data.frame} with elements \code{chr} (character, \code{position} (integer), \code{strand} ("+" or "-"), and \code{name}(character)}
  \item{...}{Arguments passed to internal function for statistics calculations}
}
\details{
At present, the approach simply calculates the contrast for each probe and tests whether the collection of these summaries are significantly non-zero.
}
\value{
  a \code{data.frame}, with the same number of rows as the input \code{regionsOfInterestTable}, but with additional columns for the statistics calculated at each promoter.
}
\author{Mark Robinson}
\seealso{\code{\link{annotationLookup}}, \code{\link{annotationLookup}}}
\examples{

# See the examples in the manual

}
