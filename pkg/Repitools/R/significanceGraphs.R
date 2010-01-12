significanceGraphs <- function(dataMatrix, lookupTable, geneList, titles, nSamples=100, confidence=0.975, legend.plot="topleft", cols=rainbow(length(geneList)), ...) {
	#Test geneList for sanity
	for (i in 1:length(geneList)) if (class(geneList[[i]])=="logical") {
		if (length(geneList[[i]])!=nrow(lookupTable)) 
		  stop("boolean geneList element length must equal num of rows in lookupTable")
	} else if (class(geneList[[i]])=="integer") {
		if(max(geneList[[i]])>nrow(lookupTable)) 
		  stop("geneList element value greater than num of rows in lookupTable") 
	} else stop("geneList elements must a be boolean or integer vector")
	stopifnot(confidence>0.5, confidence<1)
	if (length(legend.plot)!=ncol(dataMatrix)) if (length(legend.plot)!=1) stop("legend.plot must be either same length as columns in dataMatrix or 1") else legend.plot <- rep(legend.plot, ncol(dataMatrix))
	x.p <- as.numeric(colnames(lookupTable))
	#grab Intensities for all genes
	geneList.max <- which.max(sapply(geneList, FUN=function(u) if(class(u)=="logical") sum(u) else length(u)))
	for (i in 1:ncol(dataMatrix)) {
		sMat <- .scoreIntensity(lookupTable, dataMatrix[,i], removeZeros=TRUE, returnMatrix=TRUE)
		sMat.geneList <- lapply(geneList, FUN=function(u) sMat[u, ])
		trace.geneList <- sapply(sMat.geneList, FUN=function(u) apply(u, 2, median, na.rm=TRUE))

		#choose nSamples random genelists
		inds <- lapply(seq_len(nSamples), FUN=function(u) sample(nrow(sMat), nrow(sMat.geneList[[geneList.max]])))
		meds <- sapply(inds, FUN=function(u) apply(sMat[u,], 2, median,na.rm=TRUE))
		meds.conf <- apply(meds, 1, quantile, p=c(1-confidence, 0.5, confidence))
	
		#plot tiem
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="n", lty=c(2,1,2,1), lwd=c(1,3,1,3), xlab="Position relative to TSS", ylab="Signal", main=titles[i])
		polygon(x=c(x.p, rev(x.p)), y=c(meds.conf[1,], rev(meds.conf[3,])), col="lightblue")
		matplot(x.p, cbind(t(meds.conf), trace.geneList), type="l", lty=c(2,1,2,rep(1, length(geneList))), lwd=c(1,3,1,rep(3, length(geneList))), add=TRUE, col=c("blue","blue","blue",cols))
		if (!is.na(legend.plot[i])) legend(legend.plot[i], legend=names(geneList), col=cols, lwd=3)
	}
}

