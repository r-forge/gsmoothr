\name{cpgDensityCalc}
\alias{cpgDensityCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate CpG Density in a Window}
\description{Function to calculate CpG density around a position.
}
\usage{
cpgDensityCalc(locationsTable, windowSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{locationsTable}{A \code{data.frame}, with columns \code{chr} and \code{position}.}
  \item{windowSize}{Bases around the locations that are in the window. Calculation will consider \code{windowSize} bases upstream, and \code{windowSize} bases downstream.}
}
\details{
  Linear weighting of counts is used, so that CpGs further away from the position will count less.
}
\value{
  A vector of CpG densities for each region.
}
\author{Dario Strbenac}
\examples{
TSSTable <- data.frame(chr=c(1,2), position=c(100000,200000))
cpgDensityCalc(TSSTable, 600)
}