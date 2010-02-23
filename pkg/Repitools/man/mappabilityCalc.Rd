\name{mappabilityCalc}
\alias{mappabilityCalc}
\title{Calculate The Mappability of a Region}
\description{Function to calculate mappability of windows surrounding the given positions.
}
\usage{
mappabilityCalc(locationsTable, windowSize = 500)
}
\arguments{
  \item{locationsTable}{A \code{data.frame}, with columns \code{chr} and \code{position}.}
  \item{windowSize}{Bases around the locations that are in the window. Calculation will consider \code{windowSize} bases upstream, and \code{windowSize} bases downstream.}
}
\details{
The windows considered will be \code{windowSize} bases upstream and downstream of the given position. The value returned for each region is a percentage of bases in that region that are not N (any base in IUPAC nomenclature).
}
\value{
  A vector of mappability percentages, one for each region.
}
\author{Dario Strbenac}
\examples{
TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
mappabilityCalc(TSSTable, 200)
}

