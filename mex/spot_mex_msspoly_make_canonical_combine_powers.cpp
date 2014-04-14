#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

#include "spot_mex_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /*  [pow, var] = spot_mex_msspoly_make_canonical_combine_powers(powIn, varIn)
   *
   *  powIn -- m-by-n array of integers.
   *  varIn -- m-by-n array of integers, each row sorted positive integers,
   *           padded by zeros.
   *
   *  pow
   */
  mxArray *mx_powIn;
  mxArray *mx_varIn;
  int m,n;
  double *powIn;
  double *varIn;

  // Take in arguments.
  getArgSized(&mx_powIn, "double", 0, nrhs, prhs, -1, -1);
  powIn = mxGetPr(mx_powIn);
  m = mxGetM(mx_powIn);
  n = mxGetN(mx_powIn);

  getArgSized(&mx_varIn, "double", 1, nrhs, prhs, m, n);
  varIn = mxGetPr(mx_varIn);

  int i,j;

  mxArray *mx_powOut;
  mxArray *mx_varOut;
  double *powOut;
  double *varOut;

  mx_powOut = mxCreateDoubleMatrix(m,n,mxREAL);
  powOut = mxGetPr(mx_powOut);

  mx_varOut = mxCreateDoubleMatrix(m,n,mxREAL);
  varOut = mxGetPr(mx_varOut);

  //Initialize to all zeros.
  for(i = 0; i < m*n; i++){
    powOut[i] = 0;
    varOut[i] = 0;
  }

  int kmax = 0;


  for (i = 0; i < m; i++) { // For each row.
    int k = 0;
    for (j = 0; j < n; j++) { // Walk through the incoming variables.
      int idxIn  = j*m + i;
      int idxOut = k*m + i;

      if(varIn[idxIn] == 0) // We are at the end of the list.
	break;

      if(powIn[idxIn] == 0) // Skip zero powers.
	continue;

      // First item related to a power.
      if(varOut[idxOut] == 0) {
	varOut[idxOut] = varIn[idxIn];
	powOut[idxOut] = powIn[idxIn];
      } else if(varOut[idxOut] == varIn[idxIn]){
	powOut[idxOut] += powIn[idxIn];
	// With negative powers, we could be set back.
	if(powOut[idxOut] == 0){
	  varOut[idxOut] = 0;
	}
      } else { // Existing variable is there but does not match.
	k++;
	idxOut+=m;

	if(idxOut > m*n)
	  mexErrMsgTxt("Indexing exception.");
	varOut[idxOut] = varIn[idxIn];
	powOut[idxOut] = powIn[idxIn];
      }
    }
    if(k > kmax)
      kmax = k;
  }

  if(nlhs > 0)
    plhs[0] = mx_powOut;
  if(nlhs > 1)
    plhs[1] = mx_varOut;
  if(nlhs > 2) {
    mxArray *mx_kout = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *kout = mxGetPr(mx_kout);
    *kout = kmax + 1;
    plhs[2] = mx_kout;
  }
}
