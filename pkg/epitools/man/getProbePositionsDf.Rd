\name{getProbePositionsDf}
\alias{getProbePositionsDf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getProbePositionsDf(cdf, ..., verbose = -20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{ ~~Describe \code{cdf} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function (cdf, ..., verbose = -20) 
{
    chrs <- paste("chr", c(1:22, "X", "Y", "M"), sep = "")
    names(chrs) <- 1:25
    ind <- getCellIndices(cdf, ..., useNames = FALSE, unlist = TRUE, 
        verbose = verbose)
    acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
    ch <- acp[ind, 1, drop = TRUE]
    sp <- acp[ind, 2, drop = TRUE]
    data.frame(chr = chrs[as.character(ch)], position = sp, index = ind, 
        stringsAsFactors = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
