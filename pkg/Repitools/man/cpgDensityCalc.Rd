\name{cpgDensityCalc}
\alias{cpgDensityCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate CpG Density in a Window}
\description{Function to calculate CpG density around a position.
}
\usage{
cpgDensityCalc(locationsTable, windowSize=500, wFunction="linear", organism = "Hsapiens")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locationsTable}{A \code{data.frame}, with columns \code{chr} and \code{position}.}
  \item{windowSize}{Bases around the locations that are in the window. Calculation will consider \code{windowSize} bases upstream, and \code{windowSize} bases downstream.}
  \item{wFunction}{Weighting function to use. Can be \code{none}, \code{linear}, \code{log}, or \code{exp}}
  \item{organism}{A character string specifying which genome to use. Possible values are \code{Hsapiens}, \code{Amellifera}, \code{Athaliana}, \code{Btaurus}, \code{Celegans}, \code{Cfamiliaris}, \code{Dmelanogaster}, \code{Drerio}, \code{Ecoli}, \code{Ggallus}, \code{Mmusculus}, \code{Ptroglodytes}, \code{Rnorvegicus}, \code{Scerevisiae}.}
}
\details{
  For weighting scheme \code{none}, this is equivalent to the number of CG matches in the region. For \code{linear} weighting, each match is given a score \code{1/x} where \code{x} is the number of bases from the postition that the match occurred, and the scores are summed. For exponential weighting and logarithmic weighting, the idea is similar, but the scores decay exponentially (\code{exp^-5x/windowSize}) and logarithmically (\code{log2(2 - (distancesForRegion / windowSize)}).
}
\value{
  A vector of CpG densities for each region.
}
\author{Dario Strbenac}
\examples{
TSSTable <- data.frame(chr=paste("chr",c(1,2),sep=""), position=c(100000,200000))
cpgDensityCalc(TSSTable, 600)
}
