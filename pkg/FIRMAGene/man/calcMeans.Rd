\name{calcMeans}
\alias{calcMeans}
\title{Calculates columnwise "means"}
\description{Given a numeric matrix and a class variable, this function calculates the sum divided by the sqrt of the number of observations for each class, returning a gene-by-class matrix of summarized values.}
\usage{
calcMeans(x, cls)
}
\arguments{
  \item{x}{numeric matrix}
  \item{cls}{class variable (must be same length as number of columns of \code{x}}
}
\details{
This function, of course, does not calculate means, per se.  The denominator here is the square root of the number of observations.  This is to ensure the resultant summaries have the same variance.
}
\value{
numeric matrix of dimension \code{nrow(x)} by \code{length(unique(cls))}
}
\author{Mark Robinson}
\examples{
x <- matrix( rnorm(20), nr=5 )
calcMeans(x, cls=c(1,1,2,3))
}
%\keyword{ ~kwd1 }
