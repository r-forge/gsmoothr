cpgDensityCalc <- function(locationsTable, windowSize)
{
	require(MEDME)
	require(BSgenome.Hsapiens.UCSC.hg18)

	dummy <- matrix(nrow=dim(locationsTable)[1], ncol=1)
	rownames(dummy) <- 1:nrow(dummy)

	mms <- new("MEDMEset",chr=as.character(locationsTable$chr),pos=locationsTable$position,logR=dummy,organism="hsa")
	cgDensity <- CGcount(data=mms,wsize=windowSize)@CGcounts
	cgDensity
}
