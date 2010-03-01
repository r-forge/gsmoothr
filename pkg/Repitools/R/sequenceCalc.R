sequenceCalc <- function(locations, window=500, organism, pattern, fixed=TRUE, Nmask=FALSE, positions=FALSE, verbose=FALSE) {	
	chr.length <- seqlengths(organism)
	if (any((locations$position<window/2) | (locations$position+window/2>chr.length[locations$chr])))
		warning("Not all locations' windows are within chromosome boundaries.")
	
	scores <- if (positions) vector(mode='list', length=nrow(locations)) else numeric(nrow(locations))
	for(chr in unique(locations$chr)) {
		if (verbose) cat(chr, " ")
		thisChr <- which(locations$chr==chr)
		#Grab the smallest chunk of the chromosome possible
		chrStart <- min(1, min(locations$position[thisChr])-window/2+1)
		chrEnd <- max(chr.length[chr], max(locations$position[thisChr])+window/2)
		chrseq <- subseq(organism[[chr]], chrStart, chrEnd)
		if (Nmask) chrseq <- mask(organism[[chr]], pattern=pattern) 
		matches <- start(matchPattern(pattern, chrseq, fixed=fixed))
		matches.IRanges <- IRanges(start=matches, width=1)
		loc.IRanges <- IRanges(start=locations$position[thisChr]-window/2+1-chrStart, end=locations$position[thisChr]+window/2-chrStart)
		loc.overlaps <- findOverlaps(query=matches.IRanges, subject=loc.IRanges)
		loc.overlaps <- tapply(loc.overlaps@matchMatrix[,1], loc.overlaps@matchMatrix[,2], list)
		if (length(loc.overlaps)==0) next
		thisChr2 <- thisChr[as.integer(names(loc.overlaps))]
		if (positions) {
			temp <- mapply(function(x,y) matches[x]-y+chrStart-1, loc.overlaps, locations$position[thisChr2], SIMPLIFY=FALSE)
			scores[thisChr2] <- temp
		} else scores[thisChr2] <- sapply(loc.overlaps, length)
	}
	if (verbose) cat("\n")
	return(scores)
}
