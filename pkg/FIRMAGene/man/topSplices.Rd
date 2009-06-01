\name{topSplices}
\alias{topSplices}
\title{Highlights the top differentially spliced  genes}
\description{Using the simulated null distribution, the top FIRMAGene (gene-sample/class) scores are retrieved.}
\usage{
topSplices(fgObj, n = 20)
}
\arguments{
  \item{fgObj}{an object output from 'FIRMAGene'}
  \item{n}{The number of top differentially spliced genes to }
}
\value{
A table giving the gene identifier, sample (or class) and FIRMAGene score for the top 'n' genes
}
\author{Mark Robinson}
\seealso{\code{\link{FIRMAGene}}}
\examples{
## not run
##	
}
