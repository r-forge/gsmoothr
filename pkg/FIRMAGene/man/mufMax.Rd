\name{mufMax}
\alias{mufMax}
\title{Function to calculate MUF score of a vector}
\description{Function to calculate maximum (in magnitude) partial sum of a vector}
\usage{
mufMax(x)
}
\arguments{
  \item{x}{input (numeric) vector}
}
\value{ scalar MUF score }
\author{Mark Robinson}
\seealso{ 
\code{\link{mufC}} 
\code{\link{mufColumns}} 
}
\examples{
mufMax(-4:3)
}
