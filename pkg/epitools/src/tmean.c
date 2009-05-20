#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rmath.h>

// void tmean(double *x, double *xs, long *sp, long *n_, double *_tr, long *_np, long *_pw)
void tmean(double *x, double *xs, int *sp, int *n_, double *_tr, int *_np, int *_pw)
{

/*
  #  x  -- original probe-level score
  # xs  -- smoothed probe-level score
  # sp  -- position
  # n_  -- length of each of the above vectors
  # tr  -- amount of trim (0-0.49)
  # np  -- (min) number of probes
  # _pw -- probe window
*/

  int ii=0, jj=0, kk=0, lo, hi, st=0, en=0, n=*n_, np=*_np, pw=*_pw;
  double dummy[1000], sm=0., tr=*_tr;

/*  
  for (ii=0; ii< n; ii++)
    printf("ii=%d data=%f smoothed=%f pos=%d\n", ii, x[ii], xs[ii], sp[ii]);
  printf("probewindow=%d\n", pw);
  printf("numprobes=%d\n", np);
  printf("trim=%f\n", *tr);
  printf("length=%d\n", n);
*/

  for (ii=0; ii< n; ii++) {
  
	/* find set of probes to use */
    while( (sp[ii]-sp[st]) > pw )
	  st++;
    while( ((sp[en]-sp[ii]) < pw) && (en < (n-1)) )
	  en++;
	if ((en-st+1) < np)
	  continue;
	
	/* assign current vector to dummy */
	kk=0;
	for(jj=st; jj <= en; jj++)
	  dummy[kk++] = x[jj];
	  
	/* sort vector inline */
	/*
    printf("before sort:");
	for (jj=0; jj< kk; jj++)
      printf(" %f", dummy[jj]);
    printf("\nafter sort:");
	*/
	R_rsort(dummy, kk);
	/*
	for (jj=0; jj< kk; jj++)
      printf(" %f", dummy[jj]);
    printf("\n");
	*/

	/* take mean of (internal portion) of dummy */
	sm = 0.;
    lo = floor((float) kk * tr);
    hi = kk-lo-1;
    // printf("lo=%d hi=%d H-L+1=%d\n", lo, hi, (hi-lo+1));
	for(jj=lo; jj <= hi; jj++)
	  sm += dummy[jj];
	
	/* assign smoothed data */
	xs[ii] = sm/sqrt((float)(hi-lo+1));
  }


/*
# ------------------------	
# R code for trimmed mean
# ------------------------	
#trimmedMean <- function(pos, score, probeWindow=600, meanTrim=.1, nProbes=10) {
#  st <- 1
#  en <- 1
#  n <- length(pos)
#  stopifnot( length(score)==n )
#  tmean <- rep(0,n)
#  for(ii in 1:n) {
#    while( (pos[ii]-pos[st]) > probeWindow )
#      st <- st + 1
#    while( (pos[en]-pos[ii]) < probeWindow & en < (n-1))
#      en <- en + 1
#    if ( (en-st+1) < nProbes )
#      next
#    tmean[ii] <- mean( score[st:en], trim=meanTrim )*sqrt(en-st+1)
#  }
#  tmean
#}
# ------------------------	
*/

}
