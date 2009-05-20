\name{loadPairFile}
\alias{loadPairFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadPairFile(filename, ndf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ ~~Describe \code{filename} here~~ }
  \item{ndf}{ ~~Describe \code{ndf} here~~ }
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
function (filename, ndf) 
{
    pairTemp <- read.table(filename, sep = "\t", header = T, 
        stringsAsFactors = F)
    pairTemp$PM[pairTemp$PM == 0] = 1
    return(log2(pairTemp$PM[match(ndf$index, pairTemp$Y * 768 + 
        pairTemp$X)]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
