\name{tmeanC}
\alias{tmeanC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tmeanC(sp, x, nProbes = 10, probeWindow = 600, trim = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp}{ ~~Describe \code{sp} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{nProbes}{ ~~Describe \code{nProbes} here~~ }
  \item{probeWindow}{ ~~Describe \code{probeWindow} here~~ }
  \item{trim}{ ~~Describe \code{trim} here~~ }
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
function (sp, x, nProbes = 10, probeWindow = 600, trim = 0.1) 
{
    n <- length(x)
    stopifnot(length(x) == length(sp))
    out <- .C("tmean", x = as.double(x), xs = double(n), sp = as.integer(sp), 
        n = as.integer(n), tr = as.double(trim), np = as.integer(nProbes), 
        pw = as.integer(probeWindow), PACKAGE = "epitools")
    out$xs
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
