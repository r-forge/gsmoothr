\name{trimmedMean}
\alias{trimmedMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
trimmedMean(pos, score, probeWindow = 600, meanTrim = 0.1, nProbes = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{score}{ ~~Describe \code{score} here~~ }
  \item{probeWindow}{ ~~Describe \code{probeWindow} here~~ }
  \item{meanTrim}{ ~~Describe \code{meanTrim} here~~ }
  \item{nProbes}{ ~~Describe \code{nProbes} here~~ }
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
function (pos, score, probeWindow = 600, meanTrim = 0.1, nProbes = 10) 
{
    st <- 1
    en <- 1
    n <- length(pos)
    stopifnot(length(score) == n)
    tmean <- rep(0, n)
    for (ii in 1:n) {
        while ((pos[ii] - pos[st]) > probeWindow) st <- st + 
            1
        while ((pos[en] - pos[ii]) < probeWindow & (en < n)) en <- en + 
            1
        if ((en - st + 1) < nProbes) 
            next
        lo <- floor((en - st + 1) * meanTrim) + 1
        hi <- (en - st + 1) + 1 - lo
        tmean[ii] <- mean(score[st:en], trim = meanTrim) * sqrt(hi - 
            lo + 1)
    }
    tmean
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
