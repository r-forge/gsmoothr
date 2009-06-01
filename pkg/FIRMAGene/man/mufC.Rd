\name{mufC}
\alias{mufC}
\title{Calls the C code for the partial sum 'material unaccounted for' (MUF) statistic}
\description{
  Given an input vector, every partial sum is calculated and returned.
}
\usage{
mufC(v)
}
\author{Mark Robinson}
\arguments{
  \item{v}{vector of numeric values}
}
\value{
  A 'list' with the following elements
  \item{v}{original vector}
  \item{x}{vector of all partial sums}
  \item{n}{length of original vector}
}
\seealso{ 
\code{\link{mufMax}}
\code{\link{mufColumns}}
}
\examples{
 x <- rnorm(20)
 m <- mufC(x)
}
