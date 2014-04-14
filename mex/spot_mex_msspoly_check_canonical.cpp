#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>

int generateErrorCode(const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*  errno = spot_mex_msspoly_check_canonical(dim, sub, var, pow, coeff)
     *
     *  [errno, v] = ...;
     *
     *  v -- 1 if there is a value error.
     *
     *  errno -- Integer Error Code:
     *
     *    0     |  No error.
     *    [1-5] |  Argument of given number is not a double.
     *   1[1-5] |  Argument of given number has wrong size.
     *   16     |  (sub,pow,var,coeff) do not have matching row counts.
     *   17     |  (pow,var) do not have matching col counts.
     *   2[1-5] |  Argument of given numer is not integer / non-negative integer.
     *   32     |  Argument sub is not in legal range.
     *   33     |  A negative value of pow corresponds to a non-trig. variable.
     *   34     |  A non-zero pow corresponds to a zero variable Id.
     *   35     |  A value in coeff is NaN/Inf/Complex.
     *  102     |  [sub var pow] is not sorted.
     *  103     |  [sub var pow] has a duplicate.
     *  104     |  A row of var is not sorted.
     *  105     |  There are zero coefficients.
     *  106     |  A zero power has a non-zero variable id.
     *  107     |  Row of var has repeated entries.
     *  108     |  Column of pow is all zeros.
     */

  if(nrhs != 5)
    mexErrMsgTxt("spot_mex_msspoly_check_values requires 5 arguments");

  int err = generateErrorCode(prhs);

  if(nlhs > 0)
    plhs[0] = mxCreateDoubleScalar(err);

  if(nlhs > 1){
    int value_error = (err < 100) && (err != 0);
    plhs[1] = mxCreateDoubleScalar(value_error);
  }
}

inline bool isInt(double value) {
  double intpart;
  return modf(value, &intpart) == 0.0;
}

const int TEST_NONNEG = 1;
const int DONT_TEST_NONNEG = 0;

inline bool isTrigId(int v) {
  return v % 2;
}

inline bool isIntArray(const int cnt, const double *values, const bool testPos) {
  int i;
  double intpart;
  for(i = 0; i < cnt; i++){
    if(!isInt(values[i]))
      return 0;
    if (testPos == TEST_NONNEG) {
      if (values[i] < 0.0)
	return 0;
    }
  }
  return 1;
}

int generateErrorCode(const mxArray *prhs[])
{
  int canon_err = 0;
  
  /* Test type of each argument. */
  const char *sdouble = "double";
  int i;
  
  for(i = 0; i < 5; i++){
    if(!mxIsClass(prhs[i],sdouble))
      return i+1;
  }
  
  const mxArray *mx_dim   = prhs[0];
  const mxArray *mx_sub   = prhs[1];
  const mxArray *mx_var   = prhs[2];
  const mxArray *mx_pow   = prhs[3];
  const mxArray *mx_coeff = prhs[4];

  if(mxGetM(mx_dim) != 1 || mxGetN(mx_dim) != 2)
    return 11;

  if(mxGetN(mx_sub) != 2)
    return 12;

  if(mxGetN(mx_coeff) != 1)
    return 15;

  int rows = mxGetM(mx_sub);
  
  if(rows != mxGetM(mx_pow) ||
     rows != mxGetM(mx_var) ||
     rows != mxGetM(mx_coeff))
    return 16;

  int cols = mxGetN(mx_pow);
  
  if(cols != mxGetN(mx_var))
    return 17;

  double *dim = mxGetPr(mx_dim);
  double *sub = mxGetPr(mx_sub);
  double *pow = mxGetPr(mx_pow);  
  double *var = mxGetPr(mx_var);
  double *coeffR = mxGetPr(mx_coeff);
  double *coeffI = NULL;

  bool complexCoeff = mxIsComplex(mx_coeff);

  if(complexCoeff)
    coeffI = mxGetPi(mx_coeff);
  
  if(!isIntArray(2,dim,TEST_NONNEG))
    return 21;

  int m = dim[0];
  int n = dim[1];

  /* The careful reader will at this point notice it was an
     error to store everything "row-wise" in terms of memory-locality. */

  int r;
  for(r = 0; r < rows; r++){
    int subI = sub[r];
    int subJ = sub[rows+r];

    if(!isInt(subI) || !isInt(subJ))
      return 22;

    /* It is important that this test go here
       and not later, see the line inside the for loop below. */
    if(subI <= 0 || subI > m || subJ <= 0 || subJ > n)
      return 32;

  }

  if(!isIntArray(rows*cols,var,TEST_NONNEG))
    return 24;

  int j;

  for(i = 0; i < rows; i++){
    for(j = 0; j < cols; j++){
      int idx = j*rows + i;

      if(!isInt(pow[idx]))
	return 23;

      if(pow[idx] < 0 && !isTrigId((int)var[idx]))
	return 33;

      if(var[idx] == 0.0 && pow[idx] != 0.0)
	return 34;

      if(var[idx] != 0.0 && pow[idx] == 0.0)
	canon_err = 106;


      if(j > 0){
	int prevIdx = (j-1)*rows + i;
	if(var[idx] != 0.0 && var[prevIdx] == var[idx]){
	  canon_err = 107;
	}

	if(var[prevIdx] < var[idx]){
	  canon_err = 104;
	}
      }
    }
  }

  for(i = 0; i < rows; i++){
    if(complexCoeff){
      if(!mxIsFinite(coeffR[i]) || !mxIsFinite(coeffI[i]))
	return 35;
      if(coeffR[i] == 0.0 && coeffI[i] == 0.0)
	canon_err = 105;
    } else {
      if(!mxIsFinite(coeffR[i]))
	return 35;
      if(coeffR[i] == 0.0)
	canon_err = 105;
    }
  }

  /* Return any error found so far, o.w. continue */
  if(canon_err != 0)
    return canon_err;

  /* The code below will exit with an error if it finds one. */
  
  int prev_subI = 0;
  int prev_subJ = 0;
  for(r = 0; r < rows; r++){
    int subI = sub[r];
    int subJ = sub[rows+r];

    /* First case: indices not sorted */
    if(prev_subI > subI ||
       ( prev_subI == subI && prev_subJ > subJ)){
      return 102;
    } else if(prev_subI == subI && prev_subJ == subJ) {
      /* Possible to have a duplicate (103)
	 or be out of order (102) */
      int c;
      int dup = 1;
      int varSorted = 0;
      int powSorted = 0;
      for(c = 0; c < cols; c++){
	/* The loop invariant:
	   varSorted = 0, powSorted = {-1,0,1} */
	int idx = r+c*rows;
	int prevIdx = r-1+c*rows;

	if(var[idx] > var[prevIdx]){
	  varSorted = 1;  /* Leave with success */
	  break;
	} else if(var[idx] < var[prevIdx]){
	  varSorted = -1; /* Leave with failure */
	  break;
	}

	if(powSorted == 0){
	  if(pow[idx] > pow[prevIdx])
	    powSorted = 1;
	  else if(pow[idx] < pow[prevIdx])
	    powSorted = -1;
	}
      }
      if(varSorted == 0 && powSorted == 0){
	return 103;
      }
      if(varSorted == -1 || (varSorted == 0 && powSorted == -1)){
	return 102;
      }
      
      
      
    }

    prev_subI = subI;
    prev_subJ = subJ;
  }

  /* Test for any columns that are all zeros. */
  if(rows == 0){
    return 0;
  }

  for(j = 0; j < cols; j++){ /* For each column */
    int all_zero = 1;
    for(i = 0; i < rows ; i ++){
      int idx = j*rows + i;

      if(pow[idx] != 0){
	all_zero = 0;
	break;
      }
    }
    if(all_zero){
      return 108;
    }
  }

  
  return 0;
}
