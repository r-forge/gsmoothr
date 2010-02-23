cpgDensityCalc <- function(locationsTable, windowSize=500, wFunction=c("linear","exp","log","none"), organism = "Hsapiens")
{

	pkgName <- switch(organism,
		Hsapiens = "BSgenome.Hsapiens.UCSC.hg18",
		Amellifera = "BSgenome.Amellifera.UCSC.apiMel2",
		Athaliana = "BSgenome.Athaliana.TAIR.04232008",
		Btaurus = "BSgenome.Btaurus.UCSC.bosTau4",
		Celegans = "BSgenome.Celegans.UCSC.ce2",
		Cfamiliaris = "BSgenome.Cfamiliaris.UCSC.canFam2",
		Dmelanogaster = "BSgenome.Dmelanogaster.UCSC.dm3",
		Drerio = "BSgenome.Drerio.UCSC.danRer5",
		Ecoli = "BSgenome.Ecoli.NCBI.20080805",
		Ggallus = "BSgenome.Ggallus.UCSC.galGal3",
		Mmusculus = "BSgenome.Mmusculus.UCSC.mm9",
		Ptroglodytes = "BSgenome.Ptroglodytes.UCSC.panTro2",
		Rnorvegicus = "BSgenome.Rnorvegicus.UCSC.rn4",
		Scerevisiae = "BSgenome.Scerevisiae.UCSC.sacCer2",
	)
	if(is.null(pkgName))
		stop("Error : The organism name you provided does not match any of the known sequences. See ?cpgDensityCalc for supported organism names.")

	do.call(require, list(pkgName))

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
	
	cpgDensity <- numeric()
	sequences <- do.call(getSeq, list(x = get(organism), names = locationsTable$chr, start = locationsTable$position - windowSize / 2, end = locationsTable$position + windowSize / 2))
	CGfinds <- gregexpr("CG", sequences, fixed = TRUE)
	if(wFunction == "none") {
		cpgDensity <- sapply(CGfinds, length)
	} else {
		distances <- lapply(CGfinds, function(positionsInRegion) {abs(positionsInRegion - windowSize / 2)})
		if(wFunction == "linear") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(1 - (distancesInRegion / windowSize / 2)))
		} else if(wFunction == "log") {
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(log2(2 - (distancesInRegion / windowSize / 2))))
		} else { # Exponential decay was requested.
			cpgDensity <- sapply(distances, function(distancesInRegion) sum(exp(-5 * distancesInRegion / windowSize / 2)))	
		}
	}
	rm(sequences, CGfinds)
	gc()
	return(cpgDensity)
}
