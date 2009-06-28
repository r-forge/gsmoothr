\name{topSplices}
\alias{topSplices}
\title{Highlights the top differentially spliced  genes}
\description{Using the simulated null distribution, the top FIRMAGene (gene-sample/class) scores are retrieved.}
\usage{
topSplices(fgObj, n = 20, plm=NULL, cls=NULL)
}
\arguments{
  \item{fgObj}{an object output from 'FIRMAGene'}
  \item{n}{The number of top differentially spliced genes to display}
  \item{plm}{The corresponding RmaPlm object to gather probe, default \code{NULL}}
  \item{cls}{The corresponding group identifer that was used to build the object output from 'FIRMAGene'}
}
\value{
A table giving the gene identifier, sample (or class) and FIRMAGene score for the top 'n' genes.  If \code{plm} and \code{cls} are specified (and these should match what was used as input to \code{FIRMAGene}), the probe numbers that start and end the maximum MUF score are retrieved.
}
\author{Mark Robinson}
\seealso{\code{\link{FIRMAGene}}}
\examples{
## not run
##	
}
