\name{getDetails}
\alias{getDetails}
\title{Gathers details from biomaRt}
\description{
Given a ProbeLevelModel object and information regarding the probes (start and stop locations), this function retrieves public information from Biomart
}
\usage{
getDetails(plm, probesets, id = "7952953", mart = NULL, verticalBars=FALSE, geneSymbolId="external_gene_id", transcriptClusterId="Transcript.Cluster.ID", colours = "blue", lwd = 1, o = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plm}{a ProbeLevelModel object}
  \item{probesets}{a data.frame object giving the coordinates, must have the following columns: \code{start}, \code{stop}, \code{seqname} (chromosome)}
  \item{id}{identifier to gather annotation details for}
  \item{mart}{a mart returned from \code{useMart}}
  \item{verticalBars}{logical, whether to draw vertical bars between probes/probesets}
  \item{geneSymbolId}{the Ensembl/biomaRt column to use as }
  \item{transcriptClusterId}{the column of the probesets dataframe to match the \code{id} on}
  \item{colours}{scalar or vector (of the same length as the number of samples) to use for line plot colours}
  \item{lwd}{scalar or vector (of the same length as the number of samples) to use for line plot widths}
  \item{o}{order to reorder plotting of lines}
  \item{verbose}{logical, whether to print out progress comments.}
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
