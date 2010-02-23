cpgDensityCalc <- function(locationsTable, windowSize=500, wFunction="linear", organism)
{
	if (!any(c("linear", "exp", "log", "none") %in% wFunction)) 
        	stop("wFunction needs to be one of [linear, exp, log or none] ..")
	
	
	cpgDensity <- numeric()
	require(BSgenome)
	sequences <- getSeq(x = organism, names = locationsTable$chr, start = locationsTable$position - windowSize, end = locationsTable$position + windowSize)
	CGfinds <- gregexpr("CG", sequences, fixed = TRUE)
	if(wFunction == "none") {
		cpgDensity <- sapply(CGfinds, length)
	} else {
		distances <- lapply(CGfinds, function(positionsInRegion) {abs(positionsInRegion - windowSize)})
		if(wFunction == "linear") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(1 - (distancesInRegion / windowSize)))
		} else if(wFunction == "log") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(log2(2 - (distancesInRegion / windowSize))))
		} else { # Exponential decay was requested.
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(exp(-5 * distancesInRegion / windowSize)))	
		}
	}
	rm(sequences, CGfinds)
	gc()
	return(cpgDensity)
}
