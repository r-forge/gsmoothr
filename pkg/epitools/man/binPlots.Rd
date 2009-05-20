\name{binPlots}
\alias{binPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create line plots of averaged signal across a promoter}
\description{
Using a specified ordering of genes, they are split into multiple bins.  In each bin, the signal across is summarized and displayed visually.  The plots are created in an array, one row for each column of the data matrix and one column for each ordering specified.
}
\usage{
binPlots(dataMatrix, lookupTable, orderingList, plotType = c("line", "heatmap", "boxplot"), nbins = 10, cols = NULL, lwd = 3, lty = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataMatrix}{\code{matrix} of data to use}
  \item{lookupTable}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{orderingList}{\code{list} of either numeric or factor variables, of the same length as the number of rows in \code{dataTable}}
  \item{plotType}{type of plot}
  \item{nbins}{number of bins}
  \item{cols}{vector of colours}
  \item{lwd}{line width (either scalar or vector)}
  \item{lty}{line type (either scalar or vector)}
  \item{...}{further arguments to send to \code{bxp}}
}
\details{
If \code{plotType="line"}, a line is plotted for each bin across the promoter.
  
If \code{plotType="heatmap"}, a series of bins are plotted as a heatmap.  This can be useful to display a larger number of bins.

The option \code{plotType="boxplot"} is really only suitable for 1- or 2-group plots.  In this case, instead of plotting the summarized data, it plots the whole distribution of signal across the whole promoter.
}
\value{
  Either a single- or multiple-panel figure.
}
\author{Mark Robinson}
\examples{
  # See examples in manual.
}
