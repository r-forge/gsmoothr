\name{processAffyAnnotation}
\alias{processAffyAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
processAffyAnnotation(csvFile, skip = 19, getRefseq = FALSE, ..., verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{csvFile}{ ~~Describe \code{csvFile} here~~ }
  \item{skip}{ ~~Describe \code{skip} here~~ }
  \item{getRefseq}{ ~~Describe \code{getRefseq} here~~ }
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
function (csvFile, skip = 19, ..., verbose = TRUE) 
{
    if (verbose) 
        cat("Reading file:", csvFile, "\n")
    anno <- read.csv(csvFile, sep = ",", skip = skip, header = TRUE, 
        comment.char = "", stringsAsFactors = FALSE)
    if (verbose) 
        cat("Gathering source of annotation.\n")
    mr <- as.list(anno$mrna_assignment)
    names(mr) <- anno$probeset_id
    ms <- sapply(mr, FUN = function(u) {
        s <- strsplit(u, split = " // ")[[1]]
        if (length(s) > 1) {
            if (s[2] == "ENSEMBL") {
                gtype <- strsplit(s[3], " ")[[1]][1]
                return(paste(s[2], gtype, sep = "-"))
            }
            else {
                return(s[2])
            }
        }
        else {
            return(s[1])
        }
    })
    if (verbose) 
        cat("Parsing gene symbol.\n")
    g <- as.list(anno$gene_assignment)
    names(g) <- anno$probeset_id
    h <- sapply(g, FUN = function(u) {
        s <- strsplit(u, split = " // ")[[1]]
        ifelse(length(s) > 1, s[2], s[1])
    })
    anno <- anno[, c("probeset_id", "seqname", "strand", "start", 
        "stop", "total_probes", "category")]
    data.frame(anno, symbol = h, type = ms)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
