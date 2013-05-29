mufMaxPositions <- function(object, cls, idsOfInterest, clsOfInterest, verbose=TRUE) {

  getStartStop <- function(n, idx) {
    ens <- unlist(sapply(1:n,function(u) seq(u,n)))
    sts <- rep(1:n,n:1)
    c(sts[idx],ens[idx])
  }
  
  # get ResidualSet for PLM fit
  rs <- calculateResidualSet(object,verbose=verbose)
  
  # find unit indices for specified ids
  cdf <- getCdf(object)
  ind <- indexOf(cdf,names=idsOfInterest)

  # extract subset of samples
  sampleInds <- which(cls %in% unique(clsOfInterest))
  rs <- extract(rs, sampleInds)

  # read residuals for specified ids, summarize over replicates
  rsu <- readUnits(rs,units=ind,verbose=verbose)
  rsu <- lapply(rsu,FUN=function(u,cls) {
    r <- log2(u[[1]]$eps)
    m <- mad(r)
    calcMeans(r/m,cls)
  },cls=cls[sampleInds])

  
  starts <- stops <- rep(NA,length(idsOfInterest))
  for(i in 1:length(starts)) {
    resVector <- rsu[[ idsOfInterest[i] ]][,clsOfInterest[i],drop=TRUE]
    mC <- mufC(resVector)
    ss <- getStartStop(mC$n, which.max(abs(mC$x)))
    starts[i] <- ss[1]
    stops[i] <- ss[2]
  }
  data.frame(ids=idsOfInterest, samp=clsOfInterest, probeStart=starts,probeEnd=stops)
}
