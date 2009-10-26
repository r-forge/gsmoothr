binPlots <- function(dataMatrix, lookupTable, orderingList, plotType=c("line","heatmap","boxplot"), nbins=10, cols=NULL, lwd=3, lty=1, ...) {

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
      list(breakpoints=br, intervals=cut(u,breaks=br))
	} else if(class(u)=="factor") {
	  nbins <- length(levels(u))
	  list(breakpoints=u, intervals=u)
	}
  }
  
  if(class(orderingList[[1]])=="numeric")
    label <- " Order:"
  else if (class(orderingList[[1]])=="factor")
    label <- " Factor:"
  
  breaks <- lapply(orderingList,makeBins)

  if( plotType %in% c("line","heatmap")) {
    intensScores <- array(NA,dim=c(ncol(dataMatrix), ncol(lookupTable), nbins, length(orderingList)),
	                      dimnames=list(colnames(dataMatrix),colnames(lookupTable),NULL,names(orderingList)))
  } else {
    intensScores <- vector("list",ncol(dataMatrix))
	for(k in 1:length(orderingList))
	  for(i in 1:length(intensScores))
	    intensScores[[i]][[k]] <- vector("list",length(levels(breaks[[k]][["intervals"]])))
  }
  
  for(k in 1:length(orderingList)) {
    if(class(orderingList[[k]])=="numeric")
		cat(names(orderingList)[k],": ",sep="")
	
	cutLevels <- levels( breaks[[k]][["intervals"]] )
  
    for(j in 1:length(cutLevels)) {
      lev <- cutLevels[j]
	  if(class(orderingList[[k]])=="numeric")
		cat(lev," ")
      for(i in 1:ncol(dataMatrix)) {
        if( plotType %in% c("line","heatmap")) {
          intensScores[i,,j,k] <- .scoreIntensity(lookupTable[breaks[[k]][["intervals"]]==lev,], intensities=dataMatrix[,i], 
		                                         minProbes=2, removeZeros=TRUE)
		} else {
          d <- .scoreIntensity(lookupTable[breaks[[k]][["intervals"]]==lev,], intensities=dataMatrix[,i], 
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
  if(plotType=="boxplot") {
	iS <- intensScores[[i]][[i]]
	n <- length(iS)
	df <- diff(xval)[1]
	for(j in 1:length(iS)) {
	  xvals <- as.numeric(iS[[j]]$names)
	  bxp(iS[[j]], at=xval+(j-1)*df/n,pars=list(boxwex=.7*df/n,medcol=cols[j],boxcol=cols[j],whiskcol=cols[j],outcol=cols[j]),
		  add=(j>1),show.names=(j==1),xlim=c(xvals[1]-df/2,xvals[length(xvals)]+df/2),...)
	}
  } else {
	dm <- intensScores[i,,,i]
	titName <- paste("Signal:", colnames(dataMatrix)[i], label, names(orderingList)[i], sep="")
	if(plotType=="line")
	{
	  layout(rbind(c(1, 2)), widths=c(3,2))
	  par(mai=c(1.02,0.90,0.82,0))
	  matplot(xval,dm,type="l",col=cols,lty=lty,lwd=lwd,xlab="Position relative to TSS",ylab="Signal",ylim=rng)
	  par(mai=c(1.02,0.05,0.82,0))
	  plot.new()
	  legend(x="top", title ="Line Colours", col=cols, lty = 1, legend=cutLevels)
	  print(cols)
	  print(breaks[[i]][["intervals"]])
	  par(oma = c(0, 0, 2, 0))
	  mtext(titName, line = 0.5, outer = TRUE)
	} else if(plotType=="heatmap") {
	  layout(rbind(c(1,2,3)), widths=c(1,3,1))
	  par(mai=c(1.02,0.50,0.82,0.05))
	  par(oma = c(0, 0, 0, 0))
	  image(rbind(1:nbins), col=cols,axes=F, xlab="Signal Intensity")
	  axis(2, at=(0:nbins)/nbins, labels=format(seq(rng[1], rng[2], length.out=nbins+1), digits=1))
	  par(mai=c(1.02,0.05,0.82,0.05))
	  image(xval,1:nbins,dm,xlab="Position relative to TSS",main=titName, yaxt="n", ylab="Bin",col=cols,zlim=rng)
	  par(mai=c(1.02,0.05,0.82,0.50))
	  plot(x=breaks[[i]][["breakpoints"]],y=0:nbins, type="l", yaxt="n", lwd=3,xlab="log2 Expression", yaxs="i")
	}
   }
  }
  
  invisible(intensScores)
}
