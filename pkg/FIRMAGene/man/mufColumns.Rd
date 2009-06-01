\name{mufColumns}
\alias{mufColumns}
\title{Calculates the MUF statistic for each column of a matrix}
\description{Apply the 'mufMax' function to each column of a matrix}
\usage{
mufColumns(d)
}
\arguments{
  \item{d}{matrix of numeric values}
}
\value{ vector of maximum (in magnitude) MUF scores, one for each column }
\author{Mark Robinson}

\seealso{ 
\code{\link{mufMax}}
\code{\link{mufC}}
}
\examples{
x <- matrix( rnorm(20), nr=5 )
mufColumns(x)
}
