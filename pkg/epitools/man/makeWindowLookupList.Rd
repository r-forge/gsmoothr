\name{makeWindowLookupList}
\alias{makeWindowLookupList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeWindowLookupList(indexes, offsets, starts, ends)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexes}{ ~~Describe \code{indexes} here~~ }
  \item{offsets}{ ~~Describe \code{offsets} here~~ }
  \item{starts}{ ~~Describe \code{starts} here~~ }
  \item{ends}{ ~~Describe \code{ends} here~~ }
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
function (indexes, offsets, starts, ends) 
{
    getProbesInWindow <- function(indexes, offsets, start, end) {
        if (length(indexes) == 0) 
            return(integer())
        else return(indexes[which(offsets > start & offsets < 
            end)])
    }
    lookupList <- vector(mode = "list", length = length(starts))
    names(lookupList) <- (starts + ends)/2
    for (i in 1:length(starts)) {
        cat("Processing", starts[i], "to", ends[i], "\n")
        lookupList[[i]] = vector(mode = "list", length = length(indexes))
        if (!is.null(names(indexes))) 
            names(lookupList[[i]]) <- names(indexes)
        for (j in 1:length(indexes)) lookupList[[i]][[j]] = getProbesInWindow(indexes[[j]], 
            offsets[[j]], starts[i], ends[i])
    }
    return(lookupList)
  }
}
