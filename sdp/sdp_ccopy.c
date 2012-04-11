#include <math.h>
#include "mex.h"
#include "matrix.h"

/* sdp_ccopy(u,v) copies real entries of u into real entries of v*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *u, *v;
  mwSize nu, nv, i;
    
  if (nrhs != 2) mexErrMsgTxt("Two input arguments required.");
  nu = mxGetNumberOfElements(prhs[0]);      /* length of u data */
  if (mxIsSparse(prhs[1]))
      nv = mxGetJc(prhs[1])[mxGetN(prhs[1])];
  else
      nv = mxGetNumberOfElements(prhs[1]);      /* length of v data */
  if (nu != nv) mexErrMsgTxt("Incompatible data sizes.");
  u = mxGetPr(prhs[0]);
  v = mxGetPr(prhs[1]);
  for (i = 0; i < nu; i++) v[i] = u[i];
}