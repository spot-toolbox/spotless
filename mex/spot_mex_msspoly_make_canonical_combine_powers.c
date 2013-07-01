#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

#include "spot_mex_helpers.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*  [pow,var] = spot_mex_msspoly_make_canonical_helper(powIn,varIn)
   *
   *  powIn -- m-by-n array of integers.
   *  varIn -- m-by-n array of integers, each row sorted positive integers, padded by zeros.
   *
   * 
   *
   */

  mxArray *mx_powIn;
  mxArray *mx_varIn;
  int m,n;
  double *powIn;
  double *varIn;

  // Take in arguments.
  getArgSized(&mx_powIn,"double",0,nrhs,prhs,-1,-1);
  powIn = mxGetPr(mx_powIn);
  m = mxGetM(mx_powIn);
  n = mxGetN(mx_powIn);

  getArgSized(&mx_varIn,"double",1,nrhs,prhs,m,n);
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
  

  
  for(i = 0; i < m*n; i++){
    powOut[i] = 0;
    varOut[i] = 0;
  }

  int kmax = 0;

  /* For each row. */
  for(i = 0; i < m; i++){
    int k = 0;

    for(j = 0; j < n; j++){
      int idxIn  = j*m + i;
      int idxOut = k*m + i;
      /*
	printf("(%d,%d):(%d,%d) [%f,%f] <- [%f,%f]\n",
	     i,j,idxIn,idxOut,
	     varIn[idxIn],powIn[idxIn],
	     varOut[idxOut],powOut[idxOut]);
      */
      if(varIn[idxIn] == 0)
	break;

      if(j == 0) {
	/* printf("k==0: %f\n",varIn[idxIn]); */
	varOut[idxOut] = varIn[idxIn];
      }

      if(varIn[idxIn] == varOut[idxOut]){
	/*printf("Match!\n"); */
	powOut[idxOut] += powIn[idxIn];
      } else {
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
    mxArray *mx_kout = mxCreateDoubleMatrix(1,1,mxREAL);
    double *kout = mxGetPr(mx_kout);
    *kout = kmax+1;
    plhs[2] = mx_kout;
  }
  
  return;
  

}
