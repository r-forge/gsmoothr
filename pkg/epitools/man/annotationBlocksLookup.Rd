\name{annotationBlocksLookup}
\alias{annotationBlocksLookup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
annotationBlocksLookup(probes, annotation, probeIndex = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{ ~~Describe \code{probes} here~~ }
  \item{annotation}{ ~~Describe \code{annotation} here~~ }
  \item{probeIndex}{ ~~Describe \code{probeIndex} here~~ }
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
function (probes, annotation, probeIndex = NULL) 
{
    processChunk <- function(probePositions, annotation) {
        numAnnot = nrow(annotation)
        annotProbes = list(indexes = vector(mode = "list", length = numAnnot), 
            offsets = vector(mode = "list", length = numAnnot))
        if (is.null(probePositions) || (length(probePositions) == 
            0)) 
            return(annotProbes)
        chromosomeSize = max(probePositions)
        chromosomeLookup <- rep.int(NA, chromosomeSize)
        chromosomeLookup[probePositions] = names(probePositions)
        for (i in 1:numAnnot) {
            tempProbes = chromosomeLookup[annotation$start[i]:annotation$end[i]]
            annotProbes$indexes[[i]] = as.integer(na.omit(tempProbes))
            annotProbes$offsets[[i]] = which(!is.na(tempProbes))
        }
        return(annotProbes)
    }
    processChromosome <- function(probePositions, annotation) {
        numAnnot = nrow(annotation)
        annotProbes = list(indexes = vector(mode = "list", length = numAnnot), 
            offsets = vector(mode = "list", length = numAnnot))
        if (length(probePositions) == 0) 
            return(annotProbes)
        chromosomeSize = max(probePositions)
        annotChunks = split(1:nrow(annotation), trunc(annotation$start/1e+07))
        for (i in annotChunks) {
            chunkRange = range(c(annotation$start[i], annotation$end[i]))
            chunkPositions = probePositions[(probePositions >= 
                chunkRange[1]) & (probePositions <= chunkRange[2])] - 
                chunkRange[1] + 1
            chunkAnnot = annotation[i, ]
            chunkAnnot$start = chunkAnnot$start - chunkRange[1] + 
                1
            chunkAnnot$end = chunkAnnot$end - chunkRange[1] + 
                1
            tempAnnot = processChunk(chunkPositions, chunkAnnot)
            annotProbes$indexes[i] = tempAnnot$indexes
            annotProbes$offsets[i] = tempAnnot$offsets
        }
        return(annotProbes)
    }
    if (is.null(annotation$strand)) {
        probesStrandChr <- probes$chr
        annotationStrandChr <- annotation$chr
    }
    else {
        probesStrandChr <- paste(probes$chr, probes$strand, sep = "")
        annotationStrandChr <- paste(annotation$chr, annotation$strand, 
            sep = "")
    }
    annotChr = split(1:nrow(annotation), annotationStrandChr)
    annot = list(indexes = vector(mode = "list", length = nrow(annotation)), 
        offsets = vector(mode = "list", length = nrow(annotation)))
    for (i in annotChr) {
        thisChr = annotationStrandChr[i[1]]
        cat("Processing", thisChr, "\n")
        tempIndex = which(probesStrandChr == thisChr)
        tempProbes = probes$position[tempIndex]
        if (is.null(probeIndex)) 
            names(tempProbes) <- tempIndex
        else names(tempProbes) <- probeIndex[tempIndex]
        tempAnnot = processChromosome(tempProbes, annotation[i, 
            ])
        annot$indexes[i] = tempAnnot$indexes
        annot$offsets[i] = tempAnnot$offsets
    }
    if (!is.null(rownames(annotation))) {
        names(annot$indexes) <- annotation$name
        names(annot$offsets) <- annotation$name
    }
    return(annot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
