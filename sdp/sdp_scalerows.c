#include <math.h>
#include "mex.h"

/* place in real sparse matrix y the result of scaling sparse matrix s
   by the entries of v  
   Usage:  mss_scalerows(y,v,s)   */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *y, *v, *s;
  mwSize ns, nv, ny, i;
  mwSize *irs;
    
  if (nrhs != 3) mexErrMsgTxt("Three input arguments required.");
  ny = mxGetJc(prhs[0])[mxGetN(prhs[0])];   /* nnz of y */
  ns = mxGetJc(prhs[2])[mxGetN(prhs[2])];   /* nnz of s */
  nv = mxGetNumberOfElements(prhs[1]);      /* length of v */
  irs = mxGetIr(prhs[2]);
  if (ns != ny) mexErrMsgTxt("Incompatible sparse matrices.");
  s = mxGetPr(prhs[2]);
  v = mxGetPr(prhs[1]);
  y = mxGetPr(prhs[0]);
  for (i = 0; i < ns; i++) y[i] = v[irs[i]]*s[i];
}