mappabilityCalc <- function(locationsTable, windowSize = 500)
{
	require(BSgenome.Hsapiens.UCSC.hg18)
	require(BSgenome.Hsapiens36bp.UCSC.hg18mappability)
	
	if(any(locationsTable$position < (windowSize / 2)))
		stop("Not all locations' windows are obtainable. Remove locations that are too close to the edge of the start of chromosomes.")

	locationsByChr <- split(locationsTable$position, locationsTable$chr)

	for(chrIndex in 1:length(locationsByChr))
	{
		endOfChromosome <- length(Hsapiens[[as.character(names(locationsByChr[chrIndex]))]])
		offFarEdge <- sapply(locationsByChr[[chrIndex]], function(location) location > endOfChromosome - windowSize / 2)
		if(any(offFarEdge))
			stop("Not all locations' windows are obtainable. Remove locations that are too close to the end of chromosomes.")
	}
		
	mapScores <- list()
	for(chromosomeIndex in 1:length(locationsByChr))
	{
		currentlocations <- locationsByChr[[chromosomeIndex]]
		currentSeq <- Hsapiens36bp[[names(locationsByChr)[chromosomeIndex]]]
		mapScores[[chromosomeIndex]] <- sapply(currentlocations, function(currentPosition){indices <- (currentPosition - windowSize) : (currentPosition + windowSize); 1-alphabetFrequency(currentSeq[indices], as.prob=TRUE)[["N"]]})
	}

	return(unsplit(mapScores, locationsTable$chr))
}
