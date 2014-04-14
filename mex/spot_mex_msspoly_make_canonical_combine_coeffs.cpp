#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

#include "spot_mex_helpers.h"

inline bool rowMatch(double *a, int ra, double *b, int rb, int m, int n) {
  int j;
  for(j = 0; j < n ; j++){
    if(a[j*m+ra] != b[j*m+rb])
      return 0;
  }
  return 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /*  [k, keyOut, valOut] = spot_mex_msspoly_make_canonical_combine_coeff(key, val)
   *
   *  key   -- m-by-n array, identical rows must be contiguous (e.g. key is sorted)
   *  val   -- m-by-1 array
   *
   *  keyOut -- m-by-n array
   *  valOut -- m-by-1 array
   *  k      -- 1-by-1 non-negative integer.
   *
   *  k is the nubmer of unique rows in key.  keyOut(1:k,:) contains the unique rows
   *  valOut(j) = sum_{ i | key(i,:)=keyOut(j,:)} val(i).
   *
   */

  mxArray *mx_key;
  mxArray *mx_val;
  int m,n;
  double *key;
  double *valR;
  double *valI;

  // Take in arguments.
  getArgSized(&mx_key, "double", 0, nrhs, prhs, -1, -1);
  key = mxGetPr(mx_key);
  m = mxGetM(mx_key);
  n = mxGetN(mx_key);

  getArgSized(&mx_val, "double", 1, nrhs, prhs, m, 1);
  valR = mxGetPr(mx_val);
  valI = mxGetPi(mx_val);

  int i,j;

  mxArray *mx_keyOut;
  mxArray *mx_valOut;
  int k;
  double *keyOut;
  double *valOutR;
  double *valOutI;

  bool flag = mxIsComplex(mx_val);

  mx_keyOut = mxCreateDoubleMatrix(m,n,mxREAL);
  keyOut = mxGetPr(mx_keyOut);
  
  mx_valOut = mxCreateDoubleMatrix(m,1,flag ? mxCOMPLEX : mxREAL);
  valOutR = mxGetPr(mx_valOut);
  valOutI = mxGetPi(mx_valOut);
  
  k=-1;
  for(i = 0; i < m; i++){
    if( i == 0 || !rowMatch(keyOut,k,key,i,m,n)){
      k++;
      for(j = 0; j < n; j++){
	keyOut[j*m+k] = key[j*m+i];
      }
      valOutR[k] = valR[i];
      if(flag == mxCOMPLEX)
	valOutI[k] = valI[i];
    } else {
      valOutR[k] += valR[i];
      if(flag == mxCOMPLEX)
	valOutI[k] += valI[i];
    }
  }
  if(nlhs > 0)
    plhs[0] = mxCreateDoubleScalar(k+1);
  if(nlhs > 1)
    plhs[1] = mx_keyOut;
  if(nlhs > 2)
    plhs[2] = mx_valOut;
  
}
