\name{plotDecompAroma}
\alias{plotDecompAroma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the RMA decomposition}
\description{
Creates a 4-panel plot with 
}
\usage{
plotDecompAroma(plm, id = "7952953", col = "black", lwd = 1, vlines = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plm}{a ProbeLevelModel object}
  \item{id}{text identifier of the gene of interest}
  \item{col}{scalar or vector (of the same length as the number of samples) to use for line plot colours}
  \item{lwd}{scalar or vector (of the same length as the number of samples) to use for line plot widths}
  \item{vlines}{optional vector of x-values to plot vertical lines}
}
\details{
From the 'ProbeLevelModel' fit, the original data, the probe effects, chip effects and residuals are read and plotted.
}
\value{
Called for the result, a plot.

Invisibly returns a list of all the data that is plotted.
}
\author{Mark Robinson}
\examples{
## Not run.
## plotDecompAroma(plm, id = "7952953", col = "black", lwd = 1, vlines = NULL)
}
