
FIRMAGene <- function(object, nSamples=2000, seed=1976, cls=NULL, verbose=TRUE, idsToUse=NULL, minProbes=4) {

  library("aroma.affymetrix")
  
  if (is.null(cls)) {
    ds <- getDataSet(object)
	getNames <- function(...) UseMethod("getNames")
	nm <- getNames(ds)
    cls <- 1:length(nm)
	rm(ds,nm)
  }
	
  if (is.null(idsToUse))
    idsToUse <- 1:nbrOfUnits( getCdf(object) )
  unitsToUse <- idsToUse
  
  if(verbose)
    cat("Gathering/calculating residuals.\n")
  rs<-calculateResidualSet(object,verbose=1*verbose)
  
  if(verbose)
    cat("Reading units.\n")
  rsu <- readUnits(rs,units=unitsToUse,verbose=1*verbose)

  if(verbose)
    cat("Extracting standardized residuals.\n")
	
  rsu1 <- lapply(rsu,FUN=function(u,cls) {
    r <- log2(u[[1]]$eps)
    m <- mad(r)
    calcMeans(r/m,cls)
  },cls=cls)
  nProbe <- sapply(rsu,FUN=function(u) nrow(u[[1]]$eps))
  
  rm(rsu)
  gco <- gc()
  if(verbose)
    print(gco)

  if(verbose)
    cat("Calculating MUF score (observed data).\n")
	
	cat( length(rsu1), "\n")

  upLimit <- options()$aroma.affymetrix.settings$models$RmaPlm$medianPolishThreshold[1]	
  w <- (nProbe >= minProbes) & (nProbe < upLimit)
  nProbe <- nProbe[w]
	
  mufScores <- t(sapply( rsu1[w], FUN=function(u) mufColumns(u)))
  colnames(mufScores) <- unique(cls)
  
  cat( nrow(mufScores), "\n")

  if(verbose)
    cat("Extracting residual matrix.\n")

  resMatrix <- matrix( unlist(rsu1, use.names=FALSE), byrow=TRUE, ncol=length(unique(cls)) )
  
  rm(rsu1)
  
  resMatrix <- resMatrix[rowSums(is.na(resMatrix) | is.infinite(resMatrix))==0,]
  gco <- gc()
  if(verbose)
    print(gco)
	
  if(verbose)
    cat("Calculating MUF score (residual permutation).\n")	

  v <- sort(unique(nProbe))
  names(v) <- v
  v <- as.list(v)
  set.seed(seed)
  #return(list(r=resMatrix,v=v,n=nProbe))
  z <- lapply(v,FUN=function(u) {
    if(verbose) cat(u," ",sep="")
    abs(mufColumns(matrix(sample(resMatrix, size=u*nSamples),nc=nSamples)))
  })
  if(verbose) cat("\n")
  
  mv <- lapply( z, FUN=function(u) c(mean(u),sd(u)))

  mufScoresZ <- mufScores
  for(i in 1:nrow(mufScoresZ)) {
    nPch <- as.character(nProbe[i])
    mufScoresZ[i,] <- (abs(mufScoresZ[i,])-mv[[nPch]][1])/mv[[nPch]][2]
  }
    
  list(mufScores=mufScores,firmaGeneScores=mufScoresZ,nProbe=nProbe,nullDistributions=z,nProbeNull=v,plmObject=object,cls=cls)

}
