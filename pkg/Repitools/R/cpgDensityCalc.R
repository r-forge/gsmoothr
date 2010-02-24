cpgDensityCalc <- function(locationsTable, windowSize=500, wFunction=c("linear","exp","log","none"), organism)
{
	gregexpr.2 <- function(...) sapply(gregexpr(...), function(x) {
					if (x[1]=="-1") return(integer(0))
					else return(x)
				}, USE.NAMES = FALSE)

	wFunction <- match.arg(wFunction)
	if(any(locationsTable$position < (windowSize / 2)))
		stop("Not all locations' windows are obtainable. Remove locations that are too close to the edge of the start of chromosomes.")

	locationsByChr <- split(locationsTable$position, locationsTable$chr)

	for(chrIndex in 1:length(locationsByChr))
	{
		endOfChromosome <- length(organism[[as.character(names(locationsByChr[chrIndex]))]])
		offFarEdge <- sapply(locationsByChr[[chrIndex]], function(location) location > endOfChromosome - windowSize / 2)
		if(any(offFarEdge))
			stop("Not all locations' windows are obtainable. Remove locations that are too close to the end of chromosomes.")
	}
	
	cpgDensity <- numeric()
	sequences <- getSeq(x = organism, names = locationsTable$chr, start = locationsTable$position - (windowSize / 2) + 1, end = locationsTable$position + windowSize / 2)
	CGfinds <- gregexpr.2("CG", sequences, fixed = TRUE)
	if(wFunction == "none") {
		cpgDensity <- sapply(CGfinds, length)	
	} else {
		distances <- lapply(CGfinds, function(positionsInRegion) {abs(positionsInRegion - windowSize / 2)})
		if(wFunction == "linear") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(1 - (distancesInRegion / (windowSize / 2))))
		} else if(wFunction == "log") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(log2(2 - (distancesInRegion / (windowSize / 2)))))
		} else { # Exponential decay was requested.
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(exp(-5 * distancesInRegion / (windowSize / 2))))	
		}
	}
	rm(sequences, CGfinds)
	gc()
	return(cpgDensity)
}
