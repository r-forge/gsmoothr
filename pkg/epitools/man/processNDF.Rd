\name{processNDF}
\alias{processNDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processNDF(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ ~~Describe \code{filename} here~~ }
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
function (filename) 
{
    ndfTemp <- read.table(filename, sep = "\t", header = T, stringsAsFactors = F)
    chr <- gsub("FS.*|_RS.*", "", ndfTemp$PROBE_ID)
    chr <- gsub("CHR", "chr", chr)
    chr <- gsub("chr0", "chr", chr)
    position <- ndfTemp$POSITION + trunc(nchar(ndfTemp$PROBE_SEQUENCE)/2)
    strand <- rep("+", nrow(ndfTemp))
    strand[grep("_RS", ndfTemp$PROBE_ID)] = "-"
    GC <- sapply(ndfTemp$PROBE_SEQUENCE, function(x) {
        sum(strsplit(x, split = "")[[1]] \%in\% c("C", "G"))/nchar(x)
    }, USE.NAMES = F)
    ndf <- data.frame(chr, position, strand, index = ndfTemp$Y * 
        768 + ndfTemp$X, sequence = ndfTemp$PROBE_SEQUENCE, GC, 
        stringsAsFactors = F)
    ndf <- ndf[order(ndf$chr, ndf$position, ndf$strand), ]
    return(ndf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
