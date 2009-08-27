.genePromoterStats <- function(diffs, geneCoords, design, ind, upStream, downStream, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", annot)
{
  means <- tstats <- matrix(NA, nr=nrow(geneCoords), nc=ncol(diffs), dimnames=list(NULL,colnames(design)))
  df <- rep(0,nrow(geneCoords))

  for(i in 1:nrow(geneCoords)) {
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
  colnames(xDf)[1] <- paste("df",upStream,downStream,sep=".")
  if( ncol(xDf)==5 )
    colnames(xDf)[2:5] <- paste(c("meandiff","tstats","pvals","adjpvals"), gsub(".[1-9]$","",colnames(xDf)[2:5]), sep=".")
   
  cbind(geneCoords,xDf)
}

setMethodS3("genePromoterStats", "AffymetrixCelSet", function(cs, geneCoords, design, upStream=0, downStream=2000, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", ...)
{
	require(aroma.affymetrix)

  # cs - AffymetrixCelSet to read probe-level data from
  # geneCoords - data frame giving genome coordinates
  # lookupT - (optional) table matching geneCoord giving indices of probe-level data across promoter
  # ind - (optional) 
  # dmP - data matrix of probes
  
  if( upStream < 0 ) {
    warning("Making 'upStream' a positive number.")
	upStream <- abs(upStream)
  }
  
  if( nrow(design) != nbrOfArrays(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
	
	w <- which( rowSums(design != 0) > 0 )
	cs <- extract(cs,w, verbose=verbose)

	
	probePositions <- getProbePositionsDf( getCdf(cs), verbose=verbose )
  	

		upS <- max(min(7500,upStream),2500)
		dnS <- max(min(2500,downStream),2500)
		
		pos <- rep(NA,nrow(geneCoords))
		pos[geneCoords$strand=="+"] <- geneCoords$start[geneCoords$strand=="+"]
		pos[geneCoords$strand=="-"] <- geneCoords$stop[geneCoords$strand=="-"]

		# will have to change this later ... standardize column names of this table ("seqname","probeset_id")
		genePositions <- data.frame(chr=geneCoords$seqname,position=pos, strand=geneCoords$strand, 
									row.names=geneCoords$probeset_id, stringsAsFactors=FALSE)
	
	  dmP <- log2(extractMatrix(cs,cells=ind,verbose=verbose))
	  diffs <- dmP %*% design[w,]									
																												
		# run lookup twice.  first to get a list of smaller list of probes to use
		annot <- annotationLookup(probePositions, genePositions, upS, dnS)
		pb <- unique(unlist(annot$indexes, use.names=FALSE))
		probePositions <- probePositions[pb,]
		annot <- annotationLookup(probePositions, genePositions, upS, dnS)
	
	return(.genePromoterStats(diffs, geneCoords, design, ind = probePositions$index, upStream, downStream, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", annot))
  
}
)

setMethodS3("genePromoterStats", "default", function(cs, ndf, geneCoords, design, upStream=0, downStream=2000, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", ...)
{

	probePositions <- data.frame(chr = ndf$chr, position = ndf$position, index = ndf$index, stringsAsFactors=FALSE)	
  
  if( upStream < 0 ) {
    warning("Making 'upStream' a positive number.")
	upStream <- abs(upStream)
  }
  
  if( nrow(design) != ncol(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")

	upS <- max(min(7500,upStream),2500)
	dnS <- max(min(2500,downStream),2500)
	
	pos <- rep(NA,nrow(geneCoords))
	pos[geneCoords$strand=="+"] <- geneCoords$start[geneCoords$strand=="+"]
	pos[geneCoords$strand=="-"] <- geneCoords$stop[geneCoords$strand=="-"]

	# will have to change this later ... standardize column names of this table ("seqname","probeset_id")
	genePositions <- data.frame(chr=geneCoords$seqname,position=pos, strand=geneCoords$strand, 
								row.names=geneCoords$probeset_id, stringsAsFactors=FALSE)
	
	w <- which( rowSums(design != 0) > 0 )							
	dmP <- log2(cs)
	diffs <- dmP %*% design[w,]		
									
	# run lookup twice.  first to get a list of smaller list of probes to use
	annot <- annotationLookup(probePositions, genePositions, upS, dnS)
	pb <- unique(unlist(annot$indexes, use.names=FALSE))
	probePositions <- probePositions[pb,]
	annot <- annotationLookup(probePositions, genePositions, upS, dnS)


    # saving objects to disk
    saveObject(annot, file=paste("annot.up",upS,".dn",dnS,".Rdata", sep=""))

	return(.genePromoterStats(diffs, geneCoords, design, ind = probePositions$index, upStream, downStream, verbose=-20, robust=FALSE, minNRobust=10, adjustMethod="fdr", annot))
  }
)