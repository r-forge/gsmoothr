confidenceGraphs <- function(dataMatrix, lookupTable, geneList, titles, nSamples=100, confidence=0.975, ...) {
	#Test geneList for sanity
	if (class(geneList)=="logical") {
		if (length(geneList)!=nrow(lookupTable)) 
		  stop("boolean geneList length must equal num of rows in lookupTable")
	} else if (class(geneList)=="integer") {
		if(max(geneList)>nrow(lookupTable)) 
		  stop("geneList value greater than num of rows in lookupTable") 
	} else stop("geneList must a be boolean or integer vector")
	stopifnot(confidence>=0, confidence<=1) 
	x.p <- as.numeric(colnames(lookupTable))
	#grab Intensities for all genes
	for (i in 1:ncol(dataMatrix)) {
		sMat <- .scoreIntensity(lookupTable, dataMatrix[,i], removeZeros=TRUE, returnMatrix=TRUE)
		sMat.geneList <- sMat[geneList, ]
		trace.geneList <- apply(sMat.geneList, 2, median, na.rm=TRUE)

		#choose nSamples random genelists
		if (class(geneList)=="logical") sMat.rest <- sMat[!geneList, ] else sMat.rest <- sMat[-geneList]
		inds <- lapply(seq_len(nSamples), FUN=function(u) sample(seq_len(nrow(sMat.rest)), nrow(sMat.geneList)))
		meds <- sapply(inds, FUN=function(u) apply(sMat.rest[u,],2,median,na.rm=TRUE))
		meds.conf <- apply(meds, 1, quantile, p=c(1-confidence, 0.5, confidence))
	
		#plot tiem
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="n", lty=c(2,1,2,1), lwd=c(1,3,1,3), xlab="Position relative to TSS", ylab="Signal", main=titles[i])
		polygon(x=c(x.p, rev(x.p)), y=c(meds.conf[1,], rev(meds.conf[3,])), col="lightblue")
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="l", lty=c(2,1,2,1), lwd=c(1,3,1,3), add=TRUE, col=c("blue","blue","blue","black"))

	}
}

