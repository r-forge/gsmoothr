mappabilityCalc <- function(locationsTable, windowSize = 500, organism)
{
	offFarEdge <- logical()
	for (rowIndex in 1:nrow(locationsTable))
		offFarEdge[rowIndex] <- locationsTable$position[rowIndex] > (length(organism[[as.character(locationsTable$chr[rowIndex])]]) - windowSize)
	if(locationsTable$position < windowSize || any(offFarEdge))
		stop("Not all locations' windows are obtainable. Remove locations that are too close to the edge of chromosomes.")

	locationsByChr <- split(locationsTable$position, locationsTable$chr)
	
	mapScores <- list()
	for(chromosomeIndex in 1:length(locationsByChr))
	{
		currentlocations <- locationsByChr[[chromosomeIndex]]
		currentSeq <- organism[[names(locationsByChr)[chromosomeIndex]]]
		mapScores[[chromosomeIndex]] <- sapply(currentlocations, function(currentPosition){indices <- (currentPosition - windowSize) : (currentPosition + windowSize); 1-alphabetFrequency(currentSeq[indices], as.prob=TRUE)[["N"]]})
	}

	return(unsplit(mapScores, locationsTable$chr))
}
