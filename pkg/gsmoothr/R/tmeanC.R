tmeanC<-function(sp, x, nProbes=10, probeWindow=600, trim=0.1) {
  # void tmean(double *x, double *xs, long *sp, long *n_, double *tr, long *_np, long *_pw)

  n <- length(x)
  outx <- rep(NA,n)

  stopifnot( length(x)==length(sp) )

  # screen out NA and NaNs before sending to C
  k <- !is.na(x) & !is.nan(x)
  n <- sum(k)

  out<-.C("tmean", x=as.double(x[k]), xs=double(n), sp=as.integer(sp[k]), n=as.integer(n), 
                   tr=as.double(trim), np=as.integer(nProbes), pw=as.integer(probeWindow), PACKAGE="gsmoothr")

  outx[k] <- out$xs
  outx
}
