cpgDensityCalc <- function(locationsTable, windowSize=500, wFunction="linear")
{
	require(MEDME)
	require(BSgenome.Hsapiens.UCSC.hg18)

	dummy <- matrix(nrow=dim(locationsTable)[1], ncol=2)
	rownames(dummy) <- 1:nrow(dummy)

	mms <- new("MEDMEset",chr=as.character(locationsTable$chr),pos=locationsTable$position,logR=dummy,organism="hsa")
	cgDensity <- CGcount(data=mms,wsize=windowSize, wFunction=wFunction)@CGcounts
	cgDensity
}
