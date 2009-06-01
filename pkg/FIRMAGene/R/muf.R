

mufMax<-function(x) {
  m<-mufC(x)
  a<-which.max( abs(m$x) )
  m$x[a]
}

mufColumns<-function(d) {
  apply(d,2,mufMax)
}


