#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void muf(double *v, double *x, long *n_)
{

  int i, j, k, n = *n_, count=0;
  double sum=0.;
  
  /*
  printf("%d\n",n);
  for (i=0; i< n; i++)
	printf("v[%d]=%f\n",i,v[i]);
  for (i=0; i< n*(n+1)/2; i++)
	printf("x[%d]=%f\n",i,x[i]);
  */
  
  for (i=0; i< n; i++)
    for (j=i; j<n; j++) {
	   sum=0.;
       for (k=i; k<=j; k++) {
         sum+=v[k];
	     /* printf("i=%d, j=%d, k=%d, sum=%f\n",i,j,k,sum); */
	   }
	   x[count]=sum/sqrt(j-i+1.);
	   /* printf("(denom=%d) x[%d]=%f\n",(j-i+1),count,x[count]); */
	   count+=1;
    }
   
}
   

