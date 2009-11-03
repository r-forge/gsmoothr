binPlots <- function(dataMatrix, lookupTable, orderingList, plotType=c("line","heatmap","terrain","boxplot"), nbins=10, cols=NULL, lwd=3, lty=1, verbose=FALSE, ...) {

  plotType <- match.arg(plotType)
  
  if(!length(orderingList) == ncol(dataMatrix)) {
    if (!length(orderingList) == 1)
      stop("orderingList must be either length 1 or as long as the number of columns in dataMatrix.")
      orderingIndex <- rep(1, ncol(dataMatrix))
  } else orderingIndex <- 1:ncol(dataMatrix)

  
  if( is.null(cols) ) {
    require(gplots)
	if(plotType=="line") {
	  cols <- colorpanel(nbins,"blue","green","red")
	} else {
	  cols <- colorpanel(64,"blue","white","red")
	}
  }
  
  label <- vector("list", length(orderingList))
  .makeBins <- function(u) {
    if(class(u)=="numeric") {
      br <- quantile(u,p=(0:nbins)/nbins)
      list(breakpoints=br, intervals=cut(u,breaks=br))
	} else if(class(u)=="factor") {
	  nbins <- length(levels(u))
	  list(breakpoints=u, intervals=u)
	}
  }
  for(index in 1:length(orderingList))
  {
	if(class(orderingList[[index]])=="numeric")
		label[[index]] <- " Order:"
	else if (class(orderingList[[index]])=="factor")
		label[[index]] <- " Factor:"
  }
  breaks <- lapply(orderingList, .makeBins)
  if( plotType %in% c("line","heatmap","terrain")) {
    intensScores <- array(NA,dim=c(ncol(dataMatrix), ncol(lookupTable), nbins),
	                      dimnames=list(colnames(dataMatrix),colnames(lookupTable),NULL))
  } else {
    intensScores <- vector("list",ncol(dataMatrix))
	for(i in 1:length(intensScores))
		intensScores[[i]] <- vector("list", length(levels(breaks[[orderingIndex[i]]][["intervals"]])))
  }

  xval <- as.numeric(colnames(lookupTable))

  for(i in 1:ncol(dataMatrix)) {
	if (verbose) cat(names(orderingList)[orderingIndex[i]],": ",sep="")
	cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	
	for(j in 1:length(cutLevels)){
		level <- cutLevels[j]
	    lookupTableSubset <- lookupTable[breaks[[orderingIndex[i]]][["intervals"]]==level, ]
	  if( plotType %in% c("line","heatmap","terrain")) {
	    intensScores[i,,j] <- .scoreIntensity(lookupTableSubset, intensities=dataMatrix[,i], minProbes=2, removeZeros=TRUE)	
      } else {
		d <- .scoreIntensity(lookupTableSubset, intensities=dataMatrix[,i], minProbes=2, returnMatrix=TRUE,removeZeros=TRUE)
		intensScores[[i]][[j]] <- boxplot(as.data.frame(d), plot=FALSE)
	  }
	}
  }

    if( plotType %in% c("line","heatmap","terrain"))
      rng <- range(intensScores, na.rm=TRUE)


  for(i in 1:ncol(dataMatrix)) {
  cutLevels <- levels( breaks[[orderingIndex[i]]][["intervals"]] )

	  
  if(plotType=="boxplot") {
	iS <- intensScores[[i]]
	n <- length(iS)
	df <- diff(xval)[1]
	for(j in 1:length(iS)) {
	  xvals <- as.numeric(iS[[j]]$names)
	  bxp(iS[[j]], at=xval+(j-1)*df/n,pars=list(boxwex=.7*df/n,medcol=cols[j],boxcol=cols[j],whiskcol=cols[j],outcol=cols[j]),
		  add=(j>1),show.names=(j==1),xlim=c(xvals[1]-df/2,xvals[length(xvals)]+df/2), ...)
	}
  } else {
	dm <- intensScores[i,,]
	titName <- paste("Signal:", colnames(dataMatrix)[i], label[orderingIndex[i]], names(orderingList)[orderingIndex[i]], sep="")
	if(plotType=="line")
	{
		  layout(rbind(c(1, 2)), widths=c(3,2))
		  par(oma = c(0, 0, 2, 0))
		  par(mai=c(1.02,0.90,0.82,0))
		  matplot(xval,dm,type="l",col=cols,lty=lty,lwd=lwd,xlab="Position relative to TSS",ylab="Signal",ylim=rng)
		  par(mai=c(1.02,0.05,0.82,0))
		  plot.new()
		  legend(x="top", title ="Line Colours", col=cols, lty = 1, legend=cutLevels)
		  if (verbose) print(cols)
		  intervals <- breaks[[orderingIndex[i]]][["intervals"]]
		  if (verbose) print(intervals)
		  mtext(titName, line = 0.5, outer = TRUE)
	} else if(plotType=="heatmap") {
		  layout(rbind(c(1,2,3)), widths=c(1,3,1))
		  par(mai=c(1.02,0.50,0.82,0.05))
		  par(oma = c(0, 0, 0, 0))
		  image(rbind(1:nbins), col=cols,axes=F, xlab="Signal Intensity")
		  axis(2, at=(0:nbins)/nbins, labels=format(seq(rng[1], rng[2], length.out=nbins+1), digits=1))
		  par(mai=c(1.02,0.05,0.82,0.05))
		  image(xval,1:nbins,dm,xlab="Position relative to TSS", yaxt="n", ylab="Bin",col=cols,zlim=rng)
		  par(mai=c(1.02,0.05,0.82,0.50))
		  breakpoints <- breaks[[orderingIndex[i]]][["breakpoints"]]
		  plot(x=breakpoints,y=0:nbins, type="l", yaxt="n", lwd=3,xlab="log2 Expression", yaxs="i")
		  par(oma = c(0, 0, 2, 0))
		  mtext(titName, line = 0, outer = TRUE)
		} else if(plotType=="terrain") {
		  layout(1)
		  par(oma = c(0, 0, 2, 0))
  		  dm.avg <- (dm[-1, -1] + dm[-1, -(ncol(dm) - 1)] +
             		dm[-(nrow(dm) -1), -1] + dm[-(nrow(dm) -1), -(ncol(dm) - 1)]) / 4

  		  this.cols = cols[cut(dm.avg, breaks = seq(rng[1], rng[2], length.out=length(cols)), include.lowest = T)] 
		  persp(xval, 1:nbins, dm, xlab="Position relative to TSS", yaxt="n", ylab="Bin", col=this.cols, zlim=rng, theta=-25, phi=20, d=1.5, border=NA, ticktype="detailed", zlab="Signal")
		  mtext(titName, line = 0, outer = TRUE)
		}

	  }
	  invisible(intensScores)
    }
}
