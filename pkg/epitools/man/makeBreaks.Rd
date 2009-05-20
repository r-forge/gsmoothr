\name{makeBreaks}
\alias{makeBreaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeBreaks(u, n = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ ~~Describe \code{u} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
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
function (u, n = 100) 
{
    if (ncol(u) == 6) {
        q <- quantile(as.numeric(u), p = c(0.025, 0.975), na.rm = TRUE)
        v <- max(abs(c(q[1], q[2])))
        list(breaks = seq(q[1] - 0.001, q[2] + 0.001, length = n), 
            colors = colorpanel(low = "lightgrey", mid = "green", 
                high = "darkblue", n = n - 1))
    }
    else {
        q <- quantile(as.numeric(u), p = c(0.5, 0.95), na.rm = TRUE)
        v <- max(abs(c(q[1], q[2])))
        list(breaks = seq(-v, v, length = n), colors = colorpanel(low = "lightgrey", 
            mid = "green", high = "darkblue", n = n - 1))
    }
  }
}
