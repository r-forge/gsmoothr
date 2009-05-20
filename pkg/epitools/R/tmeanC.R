tmeanC<-function(sp, x, nProbes=10, probeWindow=600, trim=0.1) {
  # void tmean(double *x, double *xs, long *sp, long *n_, double *tr, long *_np, long *_pw)

  n <- length(x)
  stopifnot( length(x)==length(sp) )

  out<-.C("tmean", x=as.double(x), xs=double(n), sp=as.integer(sp), n=as.integer(n), 
                   tr=as.double(trim), np=as.integer(nProbes), pw=as.integer(probeWindow), PACKAGE="epitools")

  out$xs
}
