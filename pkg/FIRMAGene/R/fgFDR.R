.fgFDR <- function(fgObj) {

  quanSeq <- seq(1,0.5,length=100)
  
  for(i in 1:length(quanSeq)) {
    q <- sapply(fgObj$nullDistributions,quantile,p=quanSeq[i])

    tp <- fp <- rep(NA,length(fgObj$nullDistributions))
    for(ii in 1:length(fgObj$nullDistributions)) {
      w <- which(fgObj$nProbe==fgObj$nProbeNull[[ii]])
      tp[ii] <- sum( fgObj$mufScores[w,] > q[[ii]] )
	  nd <- fgObj$nullDistributions[[ii]]
      fp[ii] <- sum( nd >= q[[ii]] )/length(w) * ncol(fgObj$mufScores)
	}
	cat( quanSeq[i], sum(fp)/sum(tp),sum(tp),"\n")
  }

}

