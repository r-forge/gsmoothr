.geneBlocksStats <- function(cs, design, ind, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", regionsOfInterestTable, annot)
{
  # cut down on the amount of data read, if some rows of the design matrix are all zeros
  w <- which( rowSums(design != 0) > 0 )
  cs <- extract(cs,w, verbose=verbose)
  dmP <- log2(extractMatrix(cs,cells=ind,verbose=verbose))

  diffs <- dmP %*% design[w,]
  means <- tstats <- matrix(NA, nr=nrow(regionsOfInterestTable), nc=ncol(diffs), dimnames=list(NULL,colnames(design)))
  df <- rep(0,nrow(regionsOfInterestTable))

  for(i in 1:nrow(regionsOfInterestTable)) {
    vind <- annot[[1]][[i]]
    vind <- unique(vind[!is.na(vind)])
    if( length(vind) < 2 )
      next
	for(j in 1:ncol(diffs)) {
	  if( robust & length(vind) >= minNRobust ) {
	    require(MASS)
        rf <- summary(rlm(diffs[vind,j]~1))
        tstats[i,j] <- rf$coef[3]
	    means[i,j] <- rf$coef[1]
        df[i] <- rf$df[2]
	  } else {
        tt <- t.test(diffs[vind,j])
        tstats[i,j] <- tt$statistic
	    means[i,j] <- tt$estimate
        df[i] <- tt$parameter
	  }
	}
  }

  pvals <- 2*pt(-abs(tstats),df)
  
  # adjust p-values for multiple testing
  adjpvals <- pvals
  for(i in 1:ncol(pvals))
    adjpvals[,i] <- p.adjust(pvals[,i],method=adjustMethod)
  
  xDf <- data.frame(df=df, meandiff=means, tstats=tstats, pvals=pvals, adjpvals=adjpvals)
  colnames(xDf)[1] <- paste("df",sep=".")
  if( ncol(xDf)==5 )
    colnames(xDf)[2:5] <- paste(c("meandiff","tstats","pvals","adjpvals"), gsub(".[1-9]$","",colnames(xDf)[2:5]), sep=".")
  
    cbind(regionsOfInterestTable, xDf)
}

setMethodS3("geneBlocksStats", "AffymetrixCelSet", function(cs, design, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", regionsOfInterestTable, ...)
{
	require(aroma.affymetrix)

  # cs - AffymetrixCelSet to read probe-level data from
  # dmP - data matrix of probes
  
  if( nrow(design) != nbrOfArrays(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
	
	probePositions <- getProbePositionsDf( getCdf(cs), verbose=verbose )
  	
	if(regionsOfInterestTable != NULL)
		annot <- annotationBlocksLookup(probePositions, regionsOfInterestTable)

	saveObject(annot, file=paste("annot",".Rdata", sep=""))
	
	return(.geneBlocksStats(cs, design, ind = probePositions$index, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", regionsOfInterestTable, annot))

}
)

setMethodS3("geneBlocksStats", "default", function(cs, ndf, design, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", regionsOfInterestTable, ...)
{

	probePositions <- data.frame(chr = ndf$chr, position = ndf$position, index = ndf$index, strand = ndf$strand, stringsAsFactors=FALSE)	
  
  if( nrow(design) != ncol(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
  
	annot <- annotationBlocksLookup(probePositions, regionsOfInterestTable)

    # saving objects to disk
    saveObject(annot, file=paste("annot",".Rdata", sep=""))

	return(.geneBlocksStats(cs, design, ind = probePositions$index, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", regionsOfInterestTable, annot))
  }
)