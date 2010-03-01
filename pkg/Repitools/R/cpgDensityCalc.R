setMethodS3("cpgDensityCalc", "GenomeDataList", function(rs, seqLen, ...) {
	return(lapply(IRanges::as.list(rs), cpgDensityCalc, seqLen, ...))
})


setMethodS3("cpgDensityCalc", "GenomeData", function(rs, seqLen, ...) {
	rs.midpt <- vector(mode='list', length=length(rs))
	names(rs.midpt) <- names(rs)
	
	for (chr in names(rs)) {
		rs.midpt[[chr]] <- data.frame(chr=chr, position=c(rs[[chr]][["+"]]+seqLen/2, rs[[chr]][["-"]]-seqLen/2), stringsAsFactors=FALSE)
	}
	rs.midpt <- do.call(rbind, rs.midpt)
	
	return(cpgDensityCalc(rs.midpt, window=seqLen, ...))
			
})

setMethodS3("cpgDensityCalc", "data.frame", function(locations, window=500, wFunction=c("linear","exp","log","none"), organism, verbose=FALSE, ...) {
	gregexpr.2 <- function(...) sapply(gregexpr(...), function(x) {
					if (x[1]=="-1") return(integer(0))
					else return(x)
				}, USE.NAMES = FALSE)

	wFunction <- match.arg(wFunction)

	chr.length <- seqlengths(organism)
	if (any((locations$position<window/2) | (locations$position+window/2>chr.length[locations$chr])))
		warning("Not all locations' windows are within chromosome boundaries.")

	CGpattern <- DNAString("CG")
	CGfinds <- vector(mode='list', length=nrow(locations))
	for (chr in unique(locations$chr)) {
		if (verbose) cat(chr, " ")
		thisChr <- which(locations$chr==chr)
		chrseq <- organism[[chr]]
		CGmatches <- start(matchPattern(CGpattern, chrseq))
		CG.IRanges <- IRanges(start=CGmatches, width=1)
		loc.IRanges <- IRanges(start=locations$position[thisChr]-window/2+1, end=locations$position[thisChr]+window/2)
		loc.overlaps <- findOverlaps(query=CG.IRanges, subject=loc.IRanges)
		loc.overlaps <- tapply(loc.overlaps@matchMatrix[,1], loc.overlaps@matchMatrix[,2], list)
		if (length(loc.overlaps)==0) next
		thisChr2 <- thisChr[as.integer(names(loc.overlaps))]
		temp <- mapply(function(x,y) CGmatches[x]-y+window/2, loc.overlaps, locations$position[thisChr2], SIMPLIFY=FALSE)
		CGfinds[thisChr2] <- temp
	}
	if (verbose) cat("\n")
	if(wFunction == "none") {
		cpgDensity <- sapply(CGfinds, length)	
	} else {
		distances <- lapply(CGfinds, function(positionsInRegion) {abs(positionsInRegion - window / 2)})
		if(wFunction == "linear") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(1 - (distancesInRegion / (window / 2))))
		} else if(wFunction == "log") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(log2(2 - (distancesInRegion / (window / 2)))))
		} else { # Exponential decay was requested.
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(exp(-5 * distancesInRegion / (window / 2))))	
		}
	}
	rm(CGfinds)
	gc()
	return(cpgDensity)
})
