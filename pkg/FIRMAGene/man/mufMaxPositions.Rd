\name{mufMaxPositions}
\Rdversion{1.1}
\alias{mufMaxPositions}
\title{Returns the maximum location of MUF score for specified identifer/sample group pairs of the given RmaPlm object}
\description{Reads the residuals, calculates the MUF score and locates the probe sequence that gave the maximum score and returns this.}
\usage{
mufMaxPositions(object, cls, idsOfInterest, clsOfInterest, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{RmaPlm object}
  \item{cls}{grouping variable (as used with \code{FIRMAGene}}
  \item{idsOfInterest}{identifiers to retrieve location of MUF score for}
  \item{clsOfInterest}{samples to retrieve (must be a subset of \code{cls} and be the same length as \code{idsOfInterest})}
  \item{verbose}{whether or not to display verbose output, default \code{TRUE}}
}
\details{
}
\value{
Returns a data.frame with the start and stop probe locations that give the maximum MUF score corresponding to the (identifier,sample class) pair that is given as input.
}
\author{
Mark Robinson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{topSplices}}
}
\examples{
##---- no example
}
