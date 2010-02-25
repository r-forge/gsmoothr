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

setMethodS3("cpgDensityCalc", "data.frame", function(locations, window=500, wFunction=c("linear","exp","log","none"), organism, shift=TRUE, ...) {
	gregexpr.2 <- function(...) sapply(gregexpr(...), function(x) {
					if (x[1]=="-1") return(integer(0))
					else return(x)
				}, USE.NAMES = FALSE)

	wFunction <- match.arg(wFunction)

	chr.length <- seqlengths(organism)
	if(any((locations$position<window/2) | (locations$position+window/2>chr.length[locations$chr]))) {
		if (!shift) stop("Not all locations' windows are within chromosome boundaries")
		warning("Not all locations' windows are within chromosome boundaries, shifting.")
		locations$position[locations$position<window/2] <- window/2
		locations$position <- ifelse((locations$position+window/2)>chr.length[locations$chr], chr.length[locations$chr]-window/2, locations$position)
	}

	sequences <- getSeq(x = organism, names = locations$chr, start = locations$position - (window / 2) + 1, end = locations$position + window / 2)
	CGfinds <- gregexpr.2("CG", sequences, fixed = TRUE)
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
	rm(sequences, CGfinds)
	gc()
	return(cpgDensity)
})
