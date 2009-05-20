\name{scoreIntensity}
\alias{scoreIntensity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreIntensity(lookup, intensities, minProbes = 1, removeZeros = FALSE, returnMatrix = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lookup}{ ~~Describe \code{lookup} here~~ }
  \item{intensities}{ ~~Describe \code{intensities} here~~ }
  \item{minProbes}{ ~~Describe \code{minProbes} here~~ }
  \item{removeZeros}{ ~~Describe \code{removeZeros} here~~ }
  \item{returnMatrix}{ ~~Describe \code{returnMatrix} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lookup, intensities, minProbes = 1, removeZeros = FALSE, 
    returnMatrix = FALSE) 
{
    windowMeans <- function(indexes, intensities) {
        x <- intensities[indexes]
        if (removeZeros) 
            x <- x[x != 0]
        return(mean(x))
    }
    if (class(lookup) == "list") {
        useLookup = sapply(lookup, length) >= minProbes
        return(median(sapply(lookup[useLookup], windowMeans, 
            intensities), na.rm = TRUE))
    }
    else if (class(lookup) == "matrix") {
        d <- lookup
        for (i in 1:ncol(lookup)) {
            d[, i] <- intensities[lookup[, i]]
        }
        if (removeZeros) 
            d[d == 0] <- NA
        if (returnMatrix) 
            return(d)
        return(apply(d, 2, median, na.rm = TRUE))
    }
    else {
        stop("lookup is neither a 'list' or a 'matrix'.")
    }
  }
}
