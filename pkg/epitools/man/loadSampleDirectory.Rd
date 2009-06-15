\name{loadSampleDirectory}
\alias{loadSampleDirectory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadSampleDirectory(path, ndf, what="Cy3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{ ~~Describe \code{path} here~~ }
  \item{ndf}{ ~~Describe \code{ndf} here~~ }
  \item{what}{ ~~Describe \code{what} here~~ }
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
function (path, ndf) 
{
    files <- dir(path, pattern = ".pair$")
    samples <- unique(gsub("_532.pair|_635.pair", "", files))
    tempData <- matrix(NA, nrow = nrow(ndf), ncol = length(samples), 
        dimnames = list(NULL, samples))
    for (i in 1:length(samples)) {
        tempData[, i] <- loadPairFile(paste(path, "/", samples[i], 
            "_532.pair", sep = ""), ndf)
    }
    return(tempData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
