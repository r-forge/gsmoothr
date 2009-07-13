#regionStats <- function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, fdrProbes = FALSE, ...) UseMethod("regionStats")

regionStats <- function(...) UseMethod("regionStats")

.regionStats<- function(diffs, design, ch, sp, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, fdrProbes = FALSE) {
  getBed <- function(scoreV, chrV, posV, cut=NULL, nProbes=10, indexExtend=0, ...) {
    if( is.null(cut) )
      stop("Need to specify 'cut'.")
    posInd <- getRegions( scoreV, chrV, nProbes=nProbes, cutoff=cut, count=FALSE , ...)
    if( is.null(posInd) )
      return(list())
    posReg <- data.frame(chr=paste("chr",chrV[posInd$start],sep=""),
                         start=posV[posInd$start-indexExtend],
	                     end=posV[posInd$end+indexExtend], score=0, stringsAsFactors=F)
    for(i in 1:nrow(posInd))
      posReg$score[i] <- round(median(scoreV[ (posInd$start[i]:posInd$end[i]) ]),3)
    posReg
  }

  
  getRegions <- function(score, chrV, nProbes=3, cutoff=2, count=TRUE, doAbsolute = TRUE, fdrProbes = FALSE) {
    getRegionsChr <- function(ind, score, nProbes=3, cutoff=2, count=TRUE, doAbsolute = TRUE, fdrProbes = FALSE) {
      if (doAbsolute)
        probes <- abs(score[ind]) > cutoff 
      else
        probes <- score[ind] > cutoff 
      df <- diff(probes)
      st <- ind[which( df==1 )+1]
      en <- ind[which(df ==-1 )]
      en <- en[ en>=st[1] ]
      st <- st[1:length(en)]
      w <- (en-st+1) >= nProbes
      if (count)
        if (fdrProbes)
          sum(w,en[w]-st[w])
        else
          sum(w,na.rm=TRUE)
      else {
        if (sum(w)==0)
          return(data.frame(start=NULL, end=NULL))
        else
	      data.frame(start=st,end=en)[w,]
      }
    }

    chrInds <- split(1:length(score), chrV)
    if (count) 
      return(sum(sapply(chrInds, getRegionsChr, score, nProbes, cutoff, count, doAbsolute, fdrProbes=fdrProbes)))
    else {
      regTable <- data.frame(start=NULL, end=NULL)
      for (i in 1:length(chrInds)) regTable <- rbind(regTable, getRegionsChr(chrInds[[i]],score, nProbes, cutoff, count, doAbsolute))  
      return(regTable)
    }
  }

  

  fdrTable <- function(realScore, permScore, chrV, minCutoff = .5, maxCutoff=max( abs(permScore), na.rm=TRUE ), cutsLength=100, nProbes = 10, verbose=TRUE, fdrProbes = FALSE, ...) {
    require(gsmoothr)
    cuts <- seq(minCutoff,maxCutoff,length=cutsLength)

    fdr <- matrix(,nr=length(cuts),nc=4)
    colnames(fdr) <- c("cutoff","neg","pos","fdr")
    for(i in 1:length(cuts)) {
      pos <- getRegions( realScore, chrV, nProbes=nProbes, cutoff=cuts[i], fdrProbes=fdrProbes, ... )
      neg <- getRegions( permScore, chrV, nProbes=nProbes, cutoff=cuts[i], fdrProbes=fdrProbes, ... )
      fdr[i,] <- c(cuts[i],neg,pos,min(neg/pos,1))
      cat(".")
    }
    cat("\n")
    as.data.frame(fdr)
  }

  
  uch <- unique(ch)

  tmeanReal <- matrix(,nr=nrow(diffs),nc=ncol(diffs))
  tmeanPerms <- lapply( as.list(colnames(design)), FUN=function(u) {
    matrix(NA,nr=nrow(diffs),nc=nPermutations)
  })
  regions <- fdrTabs <- vector("list", length(tmeanPerms))
  names(regions) <- names(fdrTabs) <- colnames(design)
  
  ifelse( verbose, print(gc()), gc())
  
  
  # calculate smoothed statistics
  for(col in 1:ncol(diffs)) {
    if( verbose )
	  cat("Calculating trimmed means for column", col, "of design matrix:\n")
    for(ii in 1:length(uch)) {
      if( verbose )
	    cat(" ", uch[ii], "-", sep="")
	  w <- which(ch == uch[ii])
      #tmeanReal[w,col] <- trimmedMean(sp[w], diffs[w,col], probeWindow=probeWindow, meanTrim=meanTrim, nProbes=nProbes)
	  
	  tmeanReal[w,col] <- gsmoothr::tmeanC(sp[w], diffs[w,col], probeWindow=probeWindow, trim=meanTrim, nProbes=nProbes)
  	  if( verbose )
	      cat("R")
	  for(j in 1:ncol(tmeanPerms[[col]])) {
	    s <- sample(1:nrow(tmeanReal))
	    #tmeanPerms[[col]][w,j] <- trimmedMean(sp[w], diffs[s,col][w], probeWindow=probeWindow, meanTrim=meanTrim, nProbes=nProbes)
		#return(list(s=s, diffs=diffs, col=col, sp=sp, w=w))
		#save(diffs,s,col,w,sp,j,file="preERROR.Rdata")
	    tmeanPerms[[col]][w,j] <- gsmoothr::tmeanC(sp[w], diffs[s,col][w], probeWindow=probeWindow, trim=meanTrim, nProbes=nProbes)
		if( verbose )
	      cat(".")

	  }
	}
    if( verbose )
      cat("\nCalculating FDR table.\n")
    # calculate FDR table
	mx <- max(abs(tmeanPerms[[col]]),na.rm=TRUE)
	z <- apply(tmeanPerms[[col]], 2, FUN=function(u) fdrTable(realScore=tmeanReal[,col], permScore=u, chrV=ch,  maxCut=mx, nProbes=nProbes, cutsLength=40, fdrProbes=fdrProbes))
	fdrTabs[[col]] <- z[[1]]
	for(i in 2:length(nPermutations))
	  fdrTabs[[col]][,2:3] <- fdrTabs[[col]][,2:3] + z[[i]][,2:3]
	fdrTabs[[col]]$fdr <- pmin(fdrTabs[[col]]$neg/fdrTabs[[col]]$pos,1)  # re-adjust FDR calculation over all permutations
	
	# select lowest cutoff such that FDR is achieved
	w <- which(fdrTabs[[col]]$fdr < fdrLevel )
	cut <- min( fdrTabs[[col]]$cut[w], na.rm=TRUE )
	
    if( verbose )
      cat("Using cutoff of", cut, "for FDR of", fdrLevel,"\n")
	  
	regions[[col]] <- getBed(tmeanReal[,col], chrV=ch, posV=sp,  nProbes=nProbes, cut=cut)
	
	
  }
  
  return(list(regions=regions,tmeanReal=tmeanReal,tmeanPerms=tmeanPerms,fdrTables=fdrTabs))

}

#regionStats.AffymetrixCelSet <- function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, ind=NULL, fdrProbes = FALSE) {
setMethodS3("regionStats","AffymetrixCelSet",function(cs, design, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, ind=NULL, fdrProbes = FALSE) {

  require(aroma.affymetrix)

  # cs - AffymetrixCelSet to read probe-level data from
  # design - design matrix
  # ind - (optional) 

  
  cdf <- getCdf(cs)
    
  if( is.null(ind) )
    ind <- getCellIndices( cdf, useNames=FALSE, unlist=TRUE)

  if( nrow(design) != nbrOfArrays(cs) )
    stop("The number of rows in the design matrix does not equal the number of columns in the probes data matrix")
	
  acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
  ch <- acp[ind,1,drop=TRUE]
  sp <- acp[ind,2,drop=TRUE]
  
  # cut down on the amount of data read, if some rows of the design matrix are all zeros
  w <- which( rowSums(design != 0) > 0 )
  cs <- extract(cs,w, verbose=verbose)
  dmP <- log2(extractMatrix(cs,cells=ind,verbose=verbose))
  
  # compute probe-level score of some contrast
  diffs <- dmP %*% design[w,]

  w <- rowSums( is.na(diffs) )==0
  if( verbose )
    cat("Removing", sum(!w), "rows, due to NAs.\n")
	
  diffs <- diffs[w,,drop=FALSE]
  ch <- ch[w]
  sp <- sp[w]
  
  rm(dmP)
  ifelse( verbose, print(gc()), gc())

  return(.regionStats(diffs, design, ch, sp, fdrLevel, nPermutations, probeWindow, meanTrim, nProbes, verbose, fdrProbes))
})

#regionStats.default <- function(cs, design, ind=NULL, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, fdrProbes=FALSE,  ndf) {
setMethodS3("regionStats","default",function(cs, design, ind=NULL, fdrLevel=0.05, nPermutations=5, probeWindow=600, meanTrim=.1, nProbes=10, verbose=TRUE, fdrProbes=FALSE,  ndf) {
  #nimblegen data

  # cut down on the amount of data read, if some rows of the design matrix are all zeros
  w <- which( rowSums(design != 0) > 0 )
  diffs = cs %*% design

  w <- rowSums( is.na(diffs) )==0
  if( verbose )
    cat("Removing", sum(!w), "rows, due to NAs.\n")
  return(.regionStats(diffs, design, ch=gsub("chr","",ndf$chr), sp=ndf$position, fdrLevel, nPermutations, probeWindow, meanTrim, nProbes, verbose, fdrProbes))
})

