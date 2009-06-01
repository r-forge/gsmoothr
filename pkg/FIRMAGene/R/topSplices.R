
topSplices <- function(fgObj, n=20) {

  x <- matrix(rank( -fgObj$firmaGeneScores), nc=ncol(fgObj$firmaGeneScores), dimnames=dimnames(fgObj$mufScores))
  z <- which( x <= n, arr.ind=TRUE)
  rn <- rownames(z)
  zz <- data.frame(ID=rownames(z),Sample=colnames(fgObj$firmaGeneScores)[z[,2]],Score=0)
  for(i in 1:nrow(z))
    zz[i,"Score"] <- fgObj$firmaGeneScores[ rn[i], z[i,"col"]]

  o <- order(zz$Score,decreasing=TRUE)
  zz[o,]
}
