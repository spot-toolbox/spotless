#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

// Fetch an mxArray from a list of arguments.  Also perform basic size and type
// checks.
//
// Args:
//   ptr: A pointer to an mxArray pointer to write to.
//   type: A string argument two be matched against mxIsClass(..., type).
//   idx: Which argument to be fetched (must satisfy 0 <= idx < nhrs).
//   nrhs: The number of arguments in the argument list.
//   prhs: An array of mxArray pointers (the argument comes form this list).
//   m: The number of rows required (-1 for no constraint).
//   n: The number of columns required (-1 for no constraint).
void getArgSized(mxArray **ptr, const char *type, int idx, int nrhs,
		 const mxArray *prhs[], int m, int n) {
  if ( nrhs <= idx ){
    mexErrMsgTxt("Need more arguments.");
  } else if (!mxIsClass(prhs[idx], type)) {
    mexErrMsgTxt("Argument has wrong type.");
  } else if (m >= 0 && mxGetM(prhs[idx]) != m){
    mexErrMsgTxt("Wrong number of rows.");
  } else if (n >= 0 && mxGetN(prhs[idx]) != n){
    mexErrMsgTxt("Wrong number of cols.");
  } else {
    *ptr = mxDuplicateArray(prhs[idx]);
  }
}

// Fetches a scalar number from an array of mxArray arguments.
//
// Args:
//   dst: Pointer to a double to be written to.
//   idx: Index of the argument to be read (must satisfy 0 <= idx < nhrs).
//   nrhs: Number of elements in prhs.
//   prhs: The list of mxArray pointers to fetch the argument from.
//   def: If nrhs <= idx, the default value to assign to the argument.
void getScalarDArgDefault(double *dst, int idx, int nrhs,
			  const mxArray *prhs[], double def) {
  if ( nrhs <= idx ) {
    *dst = def;
  } else {
    mxArray *ptr;
    getArgSized(&ptr, "double", idx, nrhs, prhs, 1, 1);
    *dst = *mxGetPr(ptr);
    mxDestroyArray(ptr);
  }
}

// Evaluates a polynomial by Horner's rule.
//
// Args:
//   d: Number of coefficients stored (degree + 1).
//   coeff: Array of coefficients.
//   x: Point to evalutate the polynomial st.
double poly_eval(int d, double *coeff, double x)
{
  int i;
  double v = 0;
  for (i = 0; i < d; i++) {
    v = x*v + coeff[d-i-1];
  }
  return v;
}
