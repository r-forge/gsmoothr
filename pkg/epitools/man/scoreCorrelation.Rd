\name{scoreCorrelation}
\alias{scoreCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreCorrelation(lookup, intensities, correlateTo, minProbes = 1, cor.method = "pearson")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lookup}{ ~~Describe \code{lookup} here~~ }
  \item{intensities}{ ~~Describe \code{intensities} here~~ }
  \item{correlateTo}{ ~~Describe \code{correlateTo} here~~ }
  \item{minProbes}{ ~~Describe \code{minProbes} here~~ }
  \item{cor.method}{ ~~Describe \code{cor.method} here~~ }
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
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lookup, intensities, correlateTo, minProbes = 1, cor.method = "pearson") 
{
    windowMeans <- function(indexes, intensities) {
        return(mean(intensities[indexes]))
    }
    useLookup = sapply(lookup, length) >= minProbes
    return(cor(sapply(lookup[useLookup], windowMeans, intensities), 
        correlateTo[useLookup], method = cor.method))
  }
}