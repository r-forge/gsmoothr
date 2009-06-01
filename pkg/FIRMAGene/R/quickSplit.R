

.splitQuick<-function(r) {
  rn3<-substr(rownames(r),1,3)
  split.matrix<-split.data.frame
  rr<-split(r,factor(rn3))
  rr<-unlist(lapply(rr,FUN=function(u) split(u,rownames(u))),recursive=FALSE)
  rr
}

