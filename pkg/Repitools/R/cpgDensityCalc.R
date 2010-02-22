cpgDensityCalc <- function(locationsTable, windowSize=500, wFunction="linear", organism = "Hsapiens")
{
	if (wFunction != "linear" && wFunction != "exp" && wFunction != "log" && wFunction != "none") 
        	stop("wFunction needs to be one of [linear, exp, log or none] ..")
	
	if(organism == "Hsapiens")
	{
		do.call(require, list("BSgenome.Hsapiens.UCSC.hg19"))
	} else if (organism == "Amellifera"){
		do.call(require, list("BSgenome.Amellifera.UCSC.apiMel2"))
	} else if (organism == "Athaliana"){
		do.call(require, list("BSgenome.Athaliana.TAIR.04232008"))
	} else if (organism == "Btaurus"){
		do.call(require, list("BSgenome.Btaurus.UCSC.bosTau4"))
	} else if (organism == "Celegans"){
		do.call(require, list("BSgenome.Celegans.UCSC.ce2"))
	} else if (organism == "Cfamiliaris"){
		do.call(require, list("BSgenome.Cfamiliaris.UCSC.canFam2"))
	} else if (organism == "Dmelanogaster"){
		do.call(require, list("BSgenome.Dmelanogaster.UCSC.dm3"))
	} else if (organism == "Drerio"){
		do.call(require, list("BSgenome.Drerio.UCSC.danRer5"))
	} else if (organism == "Ecoli"){
		do.call(require, list("BSgenome.Ecoli.NCBI.20080805"))
	} else if (organism == "Ggallus"){
		do.call(require, list("BSgenome.Ggallus.UCSC.galGal3"))
	} else if (organism == "Mmusculus"){
		do.call(require, list("BSgenome.Mmusculus.UCSC.mm9"))
	} else if (organism == "Ptroglodytes"){
		do.call(require, list("BSgenome.Ptroglodytes.UCSC.panTro2"))
	} else if (organism == "Rnorvegicus"){
		do.call(require, list("BSgenome.Rnorvegicus.UCSC.rn4"))
	} else if (organism == "Scerevisiae"){
		do.call(require, list("BSgenome.Scerevisiae.UCSC.sacCer2"))
	} else {
		stop("The organism name you provided does not match any of the known sequences. See ?cpgDensityCalc for supported organism names.")
	}
	
	cpgDensity <- numeric()
	sequences <- do.call(getSeq, list(x = get(organism), names = locationsTable$chr, start = locationsTable$position - windowSize, end = locationsTable$position + windowSize))
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
