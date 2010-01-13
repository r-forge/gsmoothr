\name{binPlots}
\alias{binPlots}
\alias{binPlots.GenomeDataList}
\alias{binPlots.matrix}
\title{Create line plots of averaged signal across a promoter}
\description{
Using a specified ordering of genes, they are split into multiple bins.  In each bin, the signal across is summarized and displayed visually.  The plots are created in an array, one row for each column of the data matrix and one column for each ordering specified.
}
\usage{
\method{binPlots}{GenomeDataList}(rs, coordinatesTable, upStream=7500, downStream=2500, by=100, bw=300, total.lib.size=TRUE, seqLen=NULL, verbose=FALSE, ...)
\method{binPlots}{matrix}(dataMatrix, lookupTable, orderingList, plotType = c("line", "heatmap", "terrain", "boxplot"), nbins = 10, cols = NULL, lwd = 3, lty = 1, sameScale = TRUE, symmScale = FALSE, verbose = FALSE, removeZeros=TRUE, useMean=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{\code{matrix} of data to use}
  \item{rs}{\code{GenomeDataList} of sequencing data to use}
  \item{lookupTable}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{orderingList}{\code{list} of either numeric or factor variables, of the same length as the number of rows in \code{lookupTable}}
  \item{plotType}{type of plot}
  \item{nbins}{number of bins}
  \item{cols}{vector of colours}
  \item{lwd}{line width (either scalar or vector)}
  \item{lty}{line type (either scalar or vector)}
  \item{sameScale}{should the scale on all plots be the same}
  \item{symmScale}{should the scale on plots be symmetrical ie centred around 0}
  \item{verbose}{Print details of processing.}
  \item{removeZeros}{Remove zero values.}
  \item{useMean}{Whether to use mean instead of median within bins (useful for seq data).}
  \item{coordinatesTable}{\code{data.frame} giving the gene coordinates or regions of interest. The column names of this are \code{chr}, \code{position}, \code{strand}}
  \item{upStream}{number of bases upstream to plot}
  \item{downStream}{number of bases upstream to plot}
  \item{by}{how many bases between data points}
  \item{bw}{how many bases up and downstream to be counted}
  \item{total.lib.size}{whether to use the total number of sequencing reads per lane to normalise.}
  \item{seqLen}{If sequencing reads need to be extended, the fragment size to be used}
  \item{...}{further arguments to send to \code{bxp}}
}
\details{
If \code{plotType="line"}, a line is plotted for each bin across the promoter.
  
If \code{plotType="heatmap"}, a series of bins are plotted as a heatmap.  This can be useful to display a larger number of bins.

If \code{plotType="terrain"}, a series of bins are plotted as a 3d-terrain map.  This can be useful to display a larger number of bins.

The option \code{plotType="boxplot"} is really only suitable for 1- or 2-group plots.  In this case, instead of plotting the summarized data, it plots the whole distribution of signal across the whole promoter.
}
\value{
  Either a single- or multiple-panel figure.
}
\author{Mark Robinson}
\examples{
  # See examples in manual.
}
