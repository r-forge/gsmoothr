mappabilityCalc <- function(locations, window = 500, organism, verbose=FALSE) {	
	chr.length <- seqlengths(organism)
	if (any((locations$position<window/2) | (locations$position+window/2>chr.length[locations$chr])))
		warning("Not all locations' windows are within chromosome boundaries.")
	
	Npattern <- DNAString("N")
	mapScores <- numeric(nrow(locations))
	for(chr in unique(locations$chr)) {
		if (verbose) cat(chr, " ")
		thisChr <- which(locations$chr==chr)
		chrseq <- organism[[chr]]
		Nmatches <- start(matchPattern(Npattern, chrseq))
		N.IRanges <- IRanges(start=Nmatches, width=1)
		loc.IRanges <- IRanges(start=locations$position[thisChr]-window/2+1, end=locations$position[thisChr]+window/2)
		loc.overlaps <- findOverlaps(query=N.IRanges, subject=loc.IRanges)
		loc.overlaps <- tapply(loc.overlaps@matchMatrix[,1], loc.overlaps@matchMatrix[,2], list)
		if (length(loc.overlaps)==0) next
		thisChr2 <- thisChr[as.integer(names(loc.overlaps))]
		mapScores[thisChr2] <- sapply(loc.overlaps, length)/window
	}
	if (verbose) cat("\n")
	return(mapScores)
}
