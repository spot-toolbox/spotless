#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

void getArgSized(mxArray **ptr, const char * type, int idx,int nrhs, const mxArray *prhs[], int m, int n)
{
  if( nrhs <= idx ){
    mexErrMsgTxt("Need more arguments.");
  } else if (!mxIsClass(prhs[idx],type)) {
    mexErrMsgTxt("Argument has wrong type.");
  } else if (m >= 0 && mxGetM(prhs[idx]) != m){
    mexErrMsgTxt("Wrong number of rows.");
  } else if (n >= 0 && mxGetN(prhs[idx]) != n){
    mexErrMsgTxt("Wrong number of cols.");
  } else {
    *ptr = mxDuplicateArray(prhs[idx]);
  }
}

void getScalarDArgDefault(double *dst,int idx,int nrhs, const mxArray *prhs[],double def)
{
  if( nrhs <= idx ){
    *dst = def;
  } else {
    mxArray *ptr;
    getArgSized(&ptr,"double",idx,nrhs,prhs,1,1);
    *dst = *mxGetPr(ptr);
    mxDestroyArray(ptr);
  }
}

double poly_eval(int d, double *coeff, double x)
{
  int i;
  double v = 0;
  for(i = 0; i < d; i++){
    v = x*v + coeff[d-i-1];
  }
  return v;
}
