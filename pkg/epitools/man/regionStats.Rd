\name{regionStats} 
\Rdversion{1.1}
\alias{regionStats}
\alias{regionStats.AffymetrixCelSet}
\alias{regionStats.default}
\alias{regionStats,AffymetrixCelSet-method}

\title{Find Regions of High MAT Scores}
\description{
The function finds the highest MAT score cutoff for a pre-specified FDR. FDR is calculated as the ratio of the number of significant regions found in a premutation-based test, to the number found in the actual experimental microarray data. The function then finds regions that have a number of probes above a pre-specified number, in regions of the genome, of a certain size.
}
\usage{
\method{regionStats}{default}

\method{regionStats}{AffymetrixCelSet}
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{A set of raw microarray intensity data files}
  \item{design}{A design matrix of how to manipulate }
  \item{fdrLevel}{Cutoff of the maximum acceptable FDR}
  \item{nPermutations}{Number of permutations to use }
  \item{probeWindow}{Size of window, in base pairs, to check for }
  \item{meanTrim}{A number representing the top and bottom fraction of ordered values in a window to be removed, before the window mean is calculated.}
  \item{nProbes}{Minimum number of probes in a window, for the region to qualify as a region of significance.}
  \item{verbose}{Print debugging comments while executing.} 
  \item{ind}{A vector of the positions of the probes on the array}
  \item{fdrProbes}{If TRUE, a count of the number of probes in all significant regions is calculated.}
  \item{ndf}{The Nimblegen Definition File for Nimblegen array data.}
}
\details{}
\value{
	\item{chr}{Chromosome number that the region is on.}
	\item{start}{Number of the base pair where the region starts. The number is from the start of the particular chromosome the region is on.}
	\item{end}{Number of the nucelotide where the region ends. The number is from the start of the particular chromosome the region is on.}
	\item{score}{MAT score of the region.}
}
\references{}
\author{Mark Robinson}
\note{}

\seealso{}
\examples{
function (cs, design, fdrLevel = 0.05, nPermutations = 5, probeWindow = 600, meanTrim = 0.1, nProbes = 10, verbose = TRUE, fdrProbes = FALSE, ...) 
UseMethod("regionStats")
}
