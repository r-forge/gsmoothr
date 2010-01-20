# TODO: Add comment
# 
# Author: aarsta
###############################################################################


laneCounts <- function(cs) {
	stopifnot("GenomeDataList" %in% class(cs))
	#require(IRanges)
	#require(BSgenome)
	#lib.sizes <- IRanges::as.list(BSgenome::gdapply(cs, function(x) length(unlist(x))))
	#lib.sizes <- lapply(lib.sizes, IRanges::as.list)
	#lib.sizes <- sapply(lib.sizes, function(x) sum(unlist(x)))
	#return(lib.sizes)	
        # modification by Mark R -- Jan 20 2010
        sapply(cs, FUN=function(u) sum( sapply(u, FUN=function(v) length(v[["+"]])+length(v[["-"]])) ) )
}
