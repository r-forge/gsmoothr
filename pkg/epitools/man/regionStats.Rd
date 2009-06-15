\name{regionStats}
\Rdversion{1.1}
\alias{regionStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calculate statistics and find regions}
\description{
Using design matrix, input data, forms a score and finds extreme values of it at a given FDR. 
}
\usage{
regionStats(cs, design, fdrLevel = 0.05, nPermutations = 5, probeWindow = 600, meanTrim = 0.1, nProbes = 10, verbose = TRUE, fdrProbes = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cs}{
%%     ~~Describe \code{cs} here~~
}
  \item{design}{
%%     ~~Describe \code{design} here~~
}
  \item{fdrLevel}{
%%     ~~Describe \code{fdrLevel} here~~
}
  \item{nPermutations}{
%%     ~~Describe \code{nPermutations} here~~
}
  \item{probeWindow}{
%%     ~~Describe \code{probeWindow} here~~
}
  \item{meanTrim}{
%%     ~~Describe \code{meanTrim} here~~
}
  \item{nProbes}{
%%     ~~Describe \code{nProbes} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{fdrProbes}{
%%     ~~Describe \code{fdrProbes} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cs, design, fdrLevel = 0.05, nPermutations = 5, probeWindow = 600, 
    meanTrim = 0.1, nProbes = 10, verbose = TRUE, fdrProbes = FALSE, 
    ...) 
UseMethod("regionStats")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
