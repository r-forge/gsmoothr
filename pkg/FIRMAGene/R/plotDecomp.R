
setMethodS3("plotDecompAroma", "ProbeLevelModel",
  plotDecompAroma<-function(plm,id="7952953",col="black",lwd=1,vlines=NULL) {

  library("aroma.affymetrix")
  
  ces <- getChipEffectSet(plm)
  rs <- getResidualSet(plm)
  ds <- getDataSet(plm)
  paf <- getProbeAffinities(plm)
  cdf <- getCdf(plm)

  if (is.integer(id))
    unit <- id
  else 
    unit <- which( getUnitNames(cdf)== id)
	
  if (length(unit) != 1)
    stop("can only have 1 unit plotted.")
	
  pe <- log2(readUnits(paf,units=unit)[[1]][[1]]$phi)
  r <- log2(readUnits(rs,units=unit)[[1]][[1]]$eps)
  d <- log2(readUnits(ds,units=unit)[[1]][[1]]$intensities)
  ch <- log2(readUnits(ces,units=unit)[[1]][[1]]$intensities)
  
  if(length(col)==ncol(r)) {
    o<-order(col,decreasing=TRUE)
    col<-col[o]
    lwd<-lwd[o]
  } else {
    o<-1:ncol(r)
  }
  par(mfrow=c(2,2))
  xl<-"Probe Index"
  matplot(d,xlab="",ylab="",main="",lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  plot(pe,type="h",main="",xlab="",ylab="")
  matplot(ch,xlab="",ylab="",main="",lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  matplot(r[,o],xlab=xl,ylab="",main="",lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  
  invisible(list(data=d,residuals=r,probeeffects=pe,chipeffects=ch,col=col,lwd=lwd))
}) # plotDecompAroma


.plotDecomp <- function(obj,id="7952953",col="black",lwd=1,vlines=NULL) {
  par(mfrow=c(2,2))
  pe<-obj@probe.coefs[[id]]
  pe<-c(pe,-sum(pe))
  r<-obj@residuals$PM.resid
  r<-r[which(rownames(r)==id),]
  ce<-obj@chip.coefs[id,]

  y<-r+outer(rep(1,nrow(r)),ce)+outer(pe,rep(1,ncol(r)))
  pay<-r+outer(rep(1,nrow(r)),ce)
  
  if(length(col)==ncol(r)) {
    o<-order(col,decreasing=TRUE)
    col<-col[o]
    lwd<-lwd[o]
  } else {
    o<-1:ncol(r)
  }
  par(mfrow=c(2,2))
  xl<-"Probe Index"
  matplot(y[,o],xlab="",ylab="",main=paste("BG Corrected, Normalized data","--",id,sep=""),lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  plot(pe,type="h",main="Probe Effects",xlab="",ylab="")
  matplot(pay[,o],xlab=xl,ylab="",main="Probe-adjusted Expression",lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  matplot(r[,o],xlab=xl,ylab="",main="Residuals",lty=1,type="l",col=col,lwd=lwd)
  abline(v=vlines,lty=3)
  
  invisible(list(data=y,residuals=r,probeadjusted=pay,probeeffects=pe,col=col,lwd=lwd))
}

