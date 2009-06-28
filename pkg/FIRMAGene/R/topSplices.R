
topSplices <- function(fgObj, n=20, plm=NULL, cls=NULL) {

  x <- matrix(rank( -fgObj$firmaGeneScores), nc=ncol(fgObj$firmaGeneScores), dimnames=dimnames(fgObj$mufScores))
  z <- which( x <= n, arr.ind=TRUE)
  rn <- rownames(z)
  zz <- data.frame(ID=rownames(z),Sample=colnames(fgObj$firmaGeneScores)[z[,2]],Score=0)
  for(i in 1:nrow(z))
    zz[i,"Score"] <- fgObj$firmaGeneScores[ rn[i], z[i,"col"]]

  if( !is.null(plm) & !is.null(cls) ) {
    mmP <- mufMaxPositions(plm, cls, as.character(zz$ID), as.character(zz$Sample), verbose=FALSE)
    zz <- cbind(zz, mmP[,-c(1:2)])
  }

  o <- order(zz$Score,decreasing=TRUE)
  zz[o,]
}
