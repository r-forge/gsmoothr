

mufC<-function(v) {
  n<-length(v)
  x<-rep(0,n*(n+1)/2)
  out<-.C("muf",v=as.double(v),x=as.double(x),n=as.integer(n),PACKAGE="FIRMAGene")
  out
}


