\name{confidenceGraphs}
\alias{confidenceGraphs}
\title{Create line plots of averaged signal across a promoter compared to random sampling}
\description{
Creates a plot where the average signal across a promoter of supplied gene lists is compared to random samplings of all genes, with a shaded confidence area.
}
\usage{
confidenceGraphs(dataMatrix, lookupTable, geneList, titles, nSamples=100, confidence=0.975, legend.plot="topleft", cols=rainbow(length(geneList)), ...)
}
\arguments{
  \item{dataMatrix}{\code{matrix} of data to use}
  \item{lookupTable}{\code{matrix}, e.g. from \code{makeWindowLookupTable}}
  \item{geneList}{Names \code{list} of \code{logical} or \code{integer} vectors, specifying the genes to be averaged and plotted}
  \item{titles}{\code{character} vector of titles for the plots}
  \item{nSamples}{number of times to randomly sample from all genes}
  \item{confidence}{percentage confidence interval to be plotted (must be > 0.5 and < 1.0}
  \item{legend.plot}{Where to plot the legend - directly passed to \code{legend}. NA suppresses the legend}
  \item{cols}{the colour for each of the genelists supplied}
  \item{...}{further arguments sent to \code{matplot}}
}
\details{
For each column of data in \code{dataMatrix} a plot is created showing the average signal of the genes specified in each list element of \code{geneList} compared to \code{nSamples} random samplings of all genes, with \code{confidence} % intervals shaded. If an element of \code{geneList} is a \code{logical} vector, its length must be the same as the number of rows in \code{lookupTable}.
}
\value{
  A plot.
}
\author{Aaron Statham}
\examples{
  # See examples in manual.
}
