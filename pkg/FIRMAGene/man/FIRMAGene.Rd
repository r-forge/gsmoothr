
\name{FIRMAGene}

\alias{FIRMAGene}


\title{FIRMAGene -- Finding Isoforms using Robust Multichip Analysis for (Affymetrix) Gene platforms}

\description{Using auxiliary information from the RMA robust fit, differential splicing can be used }


\usage{FIRMAGene(object, nSamples = 2000, seed = 1976, cls = NULL, verbose = TRUE, idsToUse = NULL, minProbes = 4) }

\arguments{

\item{object}{a ProbeLevelModel object (from package aroma.affymetrix).  In future, support will be added for \code{PLMset}}

\item{nSamples}{number of random samples to take for simulated null distribution.}

\item{seed}{seed for random number generator.}

\item{cls}{variable giving the class (aligned with the columns or sample names of the input object).}

\item{verbose}{logical, whether to display progress.}

\item{idsToUse}{indices of the units to calculate FIRMAGene scores for.}

\item{minProbes}{minimum number of probes to operate on.}

}


\details{

Basically, the following steps take place:

* gather the residual from the object in question
* standardize the residuals
* calculate the MUF (maximum partial sum) score for the observed data
* calculate the MUF (maximum partial sum) score for the simulated null data
* compare the observed scores to the null distibution

}



\value{

'FGset' object

}



\author{Mark Robinson}

\references{

Robinson MD and Speed TP (2009).  ...

}



\seealso{
topSplices
}



\examples{

## not run
# hgnetaffx <- read.csv("HuGene-1_0-st-v1.na25.hg18.transcript.csv",sep=",",skip=19,header=TRUE,comment.char="")
# u <- which(getUnitNames(cdf) %in% hgnetaffx$probeset_id[hgnetaffx$category == "main" & hgnetaffx$total_probes > 7 & hgnetaffx$total_probes < 200])
# fg <- FIRMAGene(plm, idsToUse=u, cls=cls)


}

