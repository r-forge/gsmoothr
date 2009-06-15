\name{galaxyJoin}
\alias{galaxyJoin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
galaxyJoin(tabOrFile1 = NULL, tabOrFile2 = NULL, cols1 = 1:3, cols2 = 1:3, mincols = 1, fill = c("none", "right", "left", "both"), delete = TRUE, script = "gops_join.py")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tabOrFile1}{ ~~Describe \code{tabOrFile1} here~~ }
  \item{tabOrFile2}{ ~~Describe \code{tabOrFile2} here~~ }
  \item{cols1}{ ~~Describe \code{cols1} here~~ }
  \item{cols2}{ ~~Describe \code{cols2} here~~ }
  \item{mincols}{ ~~Describe \code{mincols} here~~ }
  \item{fill}{ ~~Describe \code{fill} here~~ }
  \item{delete}{ ~~Describe \code{delete} here~~ }
  \item{script}{ ~~Describe \code{script} here~~ }
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
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tabOrFile1 = NULL, tabOrFile2 = NULL, cols1 = 1:3, 
    cols2 = 1:3, mincols = 1, fill = c("none", "right", "left", 
        "both"), delete = TRUE, script = "gops_join.py") 
{
    options(digits.secs = 6)
    label1 <- gsub("[-: .]", "", Sys.time())
    fn1 <- writeFileIfNeeded(tabOrFile1, label1)
    label2 <- gsub("[-: .]", "", Sys.time())
    fn2 <- writeFileIfNeeded(tabOrFile2, label2)
    fill <- match.arg(fill)
    if (is.null(options()$galaxyPath)) 
        stop("Need to set options(galaxyPath=\"path/to/galaxy/dist\")")
    systemCall <- paste(options()$galaxyPath, script, sep = "/")
    systemCall <- paste(systemCall, "--cols1", paste(paste(cols1, 
        collapse = ","), ",", sep = ""))
    systemCall <- paste(systemCall, "--cols2", paste(paste(cols2, 
        collapse = ","), ",", sep = ""))
    systemCall <- paste(systemCall, "--mincols", mincols)
    systemCall <- paste(systemCall, "--fill", fill)
    systemCall <- paste(systemCall, fn1$file)
    systemCall <- paste(systemCall, fn2$file)
    systemCall <- paste(systemCall, label1)
    systemCall <- paste("python", systemCall)
    cat("Calling:", systemCall, "\n")
    system(systemCall)
    if (!file.exists(label1)) 
        stop("python call did not produce output")
    join <- read.table(label1, sep = "\t", comment.char = "", 
        header = FALSE, stringsAsFactors = FALSE)
    if (delete) {
        if (fn1$written) 
            unlink(fn1$file)
        if (fn2$written) 
            unlink(fn2$file)
        unlink(label1)
    }
    join
  }
}
