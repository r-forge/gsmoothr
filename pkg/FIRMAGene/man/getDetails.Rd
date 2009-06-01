\name{getDetails}
\alias{getDetails}
\title{Gathers details from biomaRt}
\description{
Given a ProbeLevelModel object and information regarding the probes (start and stop locations), this function retrieves public information from Biomart
}
\usage{
getDetails(plm, probesets, id = "7952953", mart = NULL, verbose = FALSE, biomart = "ensembl", dataset = "hsapiens_gene_ensembl", col = "blue", lwd = 1, o = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plm}{a ProbeLevelModel object}
  \item{probesets}{a data.frame object giving the coordinates}
  \item{id}{identifier}
  \item{mart}{a mart returned from \code{useMart} (optional, default is \code{NULL})}
  \item{verbose}{logical, whether to print out progress comments.}
  \item{biomart}{biomart to use}
  \item{dataset}{biomart dataset to use}
  \item{col}{scalar or vector (of the same length as the number of samples) to use for line plot colours}
  \item{lwd}{scalar or vector (of the same length as the number of samples) to use for line plot widths}
  \item{o}{order to reorder plotting of lines}
}
\details{
This function gathers the raw data and residuals as 'ExonArray' objects (see GenomeGraphs) and creates a 'GeneRegion' object for the gene selected.  A list is returned which can be plotted with 'gdPlot'
}
\value{
a 'list' object that can be plotted with 'gdPlot' (from GenomeGraphs)
}
\author{Mark Robinson}
\examples{
## not run
}
