\name{mappabilityCalc}
\alias{mappabilityCalc}
\title{Calculate The Mappability of a Region}
\description{Function to calculate mappability of windows surrounding the given positions.
}
\usage{
mappabilityCalc(locationsTable, windowSize = 500, organism)
}
\arguments{
  \item{locationsTable}{A \code{data.frame}, with columns \code{chr} and \code{position}.}
  \item{windowSize}{Bases around the locations that are in the window. Calculation will consider \code{windowSize / 2 - 1} bases upstream, and \code{windowSize / 2} bases downstream.}
  \item{organism}{The \code{BSgenome} object to calculate mappability upon.}
}
\details{
The windows considered will be \code{windowSize / 2 - 1} bases upstream and \code{windowSize / 2} bases downstream of the given position, for each position. The value returned for each region is a percentage of bases in that region that are not N (any base in IUPAC nomenclature).
}
\value{
  A vector of mappability percentages, one for each region.
}
\author{Dario Strbenac}
\examples{
require(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
TSSTable <- data.frame(chr = paste("chr", c(1,2), sep = ""), position = c(100000, 200000))
mappabilityCalc(TSSTable, 200, organism=Hsapiens36bp)
}

