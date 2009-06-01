
calcMeans<-function(x,cls) {
  ucls <- unique(cls) 
  v <- matrix(,nrow=nrow(x),ncol=length(ucls))
  for(i in 1:length(ucls)) {
    ind<-which(cls==ucls[i])
    if( length(ind) > 1 ) {
	  #xx <- matrix(x[,ind],nc=length(ind))
	  xx <- as.matrix(x[,ind])
      #v[,i] <- rowMeans(xx/sqrt(length(ind)))
      v[,i] <- rowSums(xx)/sqrt(length(ind))
    } else {
      v[,i] <- x[,ind]
    }
  }
  colnames(v) <- ucls
  v
}

