#include <math.h>
#include "mex.h"

/*
 *
 *    mss_gsum.c:
 *    MEX-file for y=nlid_gsum(x)
 *    for m-by-n double x generates q-by-n double y, where
 *    y(k,:)=x(i(k)+1,:)+x(i(k)+1,:)+...+x(i(k+1),:) for k=1,...,q
 *    where 0=i(1)<i(2)<...<i(q+1)=n are integers such that
 *    x(i,1)<0.5 for i=q(k), k=1,...,q-1, x(i,1)>=0.5 for other i<m
 *
 */


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  double *x,*y;              /* input x, output y      */
  double z;
  mwSize m,n,q,qq,i,j;              /* integer counters */
  
  if(nrhs!=1) 
      mexErrMsgTxt("exactly one input required");
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) 
      mexErrMsgTxt("input must be a double");
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if(m*n==0)
      mexErrMsgTxt("input must be not empty");
  x=mxGetPr(prhs[0]);
  q=0;
  for(i=0;i<m-1;i++)
      if(x[i]<0.5) 
          q++;
  qq=q+1;
  plhs[0] = mxCreateDoubleMatrix(qq,n, mxREAL);
  y = mxGetPr(plhs[0]);
  q=0;
  for(i=0;i<m-1;i++) {
      for(j=0;j<n;j++)
          y[q+qq*j]+=x[i+m*j];
      if(x[i]<0.5) 
          q++;
  }
  for(j=0;j<n;j++)
      y[q+qq*j]+=x[i+m*j];
}
  
  
