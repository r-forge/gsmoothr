

binPlots <- function(dataMatrix, lookupTable, orderingList, plotType=c("line","heatmap","boxplot"), nbins=10, cols=NULL, lwd=3, lty=1, ...) {

  par(mfrow=c(ncol(dataMatrix),length(orderingList)))
  
  plotType <- match.arg(plotType)
  
  if( is.null(cols) ) {
    require(gplots)
	if(plotType=="line") {
	  cols <- colorpanel(nbins,"blue","green","red")
	} else {
	  cols <- colorpanel(64,"blue","white","red")
	}
  }
  
  label <- ""
  
  makeBins <- function(u) {
    if(class(u)=="numeric") {
      br <- quantile(u,p=(0:nbins)/nbins)
      cut(u,breaks=br)
	} else if(class(u)=="factor") {
	  nbins <- length(levels(u))
	  u
	}
  }
  
  if(class(orderingList[[1]])=="numeric")
    label <- " Order:"
  else if (class(orderingList[[1]])=="factor")
    label <- " Factor:"
  
  breakList <- lapply(orderingList,makeBins)
  
  #if( plotType=="boxplot" & length(orderingList) > 1)
  #  stop("plotType='boxplot' only suitable with

  if( plotType %in% c("line","heatmap")) {
    intensScores <- array(NA,dim=c(ncol(dataMatrix), ncol(lookupTable), nbins, length(orderingList)),
	                      dimnames=list(colnames(dataMatrix),colnames(lookupTable),NULL,names(orderingList)))
  } else {
    intensScores <- vector("list",ncol(dataMatrix))
	for(k in 1:length(orderingList))
	  for(i in 1:length(intensScores))
	    intensScores[[i]][[k]] <- vector("list",length(levels(breakList[[k]])))
  }
  
  for(k in 1:length(orderingList)) {
    cat(names(orderingList)[k],": ",sep="")
	
	cutLevels <- levels( breakList[[k]] )
  
    for(j in 1:length(cutLevels)) {
      lev <- cutLevels[j]
	  cat(lev," ")
      for(i in 1:ncol(dataMatrix)) {
        if( plotType %in% c("line","heatmap")) {
          intensScores[i,,j,k] <- .scoreIntensity(lookupTable[breakList[[k]]==lev,], intensities=dataMatrix[,i], 
		                                         minProbes=2, removeZeros=TRUE)
		} else {
          d <- .scoreIntensity(lookupTable[breakList[[k]]==lev,], intensities=dataMatrix[,i], 
		       minProbes=2, returnMatrix=TRUE,removeZeros=TRUE)
		  
          intensScores[[i]][[k]][[j]] <- boxplot(as.data.frame(d), plot=FALSE)
		}
      }
    }
    cat("\n")
  }
  
  
  
  xval <- as.numeric(colnames(lookupTable))
  
  for(i in 1:ncol(dataMatrix)) {
    if( plotType %in% c("line","heatmap"))
      rng <- range(intensScores[i,,,], na.rm=TRUE)
    for(k in 1:length(orderingList)) {
	  if(plotType=="boxplot") {
	    iS <- intensScores[[i]][[k]]
		n <- length(iS)
		df <- diff(xval)[1]
	    for(j in 1:length(iS)) {
		  xvals <- as.numeric(iS[[j]]$names)
		  bxp(iS[[j]], at=xval+(j-1)*df/n,pars=list(boxwex=.7*df/n,medcol=cols[j],boxcol=cols[j],whiskcol=cols[j],outcol=cols[j]),
		      add=(j>1),show.names=(j==1),xlim=c(xvals[1]-df/2,xvals[length(xvals)]+df/2),...)
		}
	  } else {
	    dm <- intensScores[i,,,k]
	    titName <- paste("Signal:", colnames(dataMatrix)[i], label, names(orderingList)[k], sep="")
	    if(plotType=="line")
          matplot(xval,dm,type="l",col=cols,lty=lty,lwd=lwd,main=titName,xlab="Position relative to TSS",ylab="Signal",ylim=rng)
	    else if(plotType=="heatmap")
          image(xval,1:nbins,dm,xlab="Position relative to TSS",main=titName,ylab="Bin",col=cols,zlim=rng)
	  }
	}
  }
  
  invisible(intensScores)
  

}
