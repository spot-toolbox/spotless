#include <math.h>
#include "mex.h"

/*
 *
 *    mss_gsum.c:
 *    MEX-file for y=nlid_gset(x) 
 *    for m-by-2 double x such that 
 *    x(i,2)<=x(i+1,2) whenever x(i,1)>0.5, i=1,...,m-1
 *    generates q-by-2 double y, the rows of which list 
 *    all pairs (-x(i,2),x(k,2)) such that
 *    i<k, x(i,2)<0, x(k,2)>=0, and x(j,1)>0.5 for all j=i,i+1,...,k-1

 *
 */


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  double *x,*y;              /* input x, output y      */
  mwSize m,n,q,qq,d,r,i,j,k;              /* integer counters */
  
  if(nrhs!=1) 
      mexErrMsgTxt("exactly one input required");
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) 
      mexErrMsgTxt("input must be a noncomplex double");
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if((n!=2)||(m==0))
      mexErrMsgTxt("input must be a non-empty two-column matrix");
  x=mxGetPr(prhs[0]);
  n=m-1;
  q=0;
  i=0;
  while(i<n) {
      if((x[i+m]>=0.0)||(x[i]<=0.5))  /* not a promising set begfinning */
          i++;
      else {
          d=i+1;     /* find d: the end of neighbors <0.0 */
          while(d<m) {
              if((x[d+m]>=0.0)||(x[d]<=0.5))
                  break;
              else
                  d++;
          }
          if(x[d+m]<0.0)  /* nothing to catch in this group */
              i=d+1;
          else {
              r=d;  /* find r: the end of neighbors */
              while(r<n) {
                  if(x[r]<=0.5)
                      break;
                  else
                      r++;
              }
              q+=(d-i)*(r-d+1);
              i=r+1;
          }
      }
  }
  qq=q;
  plhs[0] = mxCreateDoubleMatrix(qq,2, mxREAL);
  y = mxGetPr(plhs[0]);
  q=0;
  i=0;
  while(i<n) {
      if((x[i+m]>=0.0)||(x[i]<=0.5))  /* not a promising set begfinning */
          i++;
      else {
          d=i+1;     /* find d: the end of neighbors <0.0 */
          while(d<m) {
              if((x[d+m]>=0.0)||(x[d]<=0.5))
                  break;
              else
                  d++;
          }
          if(x[d+m]<0.0)          /* nothing to catch in this group */
              i=d+1;
          else {
              r=d;                /* find r: the end of neighbors */
              while(r<n) {
                  if(x[r]<=0.5)
                      break;
                  else
                      r++;
              }
              for(j=i;j<d;j++)
                  for(k=d;k<=r;k++) {
                      y[q]=-x[j+m];
                      y[q+qq]=x[k+m];
                      q++;
                  }
              i=r+1;
          }
      }
  }
}