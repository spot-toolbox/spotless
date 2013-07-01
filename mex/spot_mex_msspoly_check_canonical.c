#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>


inline bool isInt(double value)
{
  double intpart;
  return modf(value,&intpart) == 0.0;
}

const int TEST_NONNEG = 1;
const int DONT_TEST_NONNEG = 0;

inline bool isTrigId(int v)
{
  return v%2;
}

inline bool isIntArray(const int cnt, const double *values, const bool testPos)
{
  int i;
  double intpart;
  for(i = 0; i < cnt; i++){
    if(!isInt(values[i]))
      return 0;
    if (testPos == TEST_NONNEG){
      if(values[i] < 0.0)
	return 0;
    }
  }
  return 1;
}


int generateErrorCode(const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*  errno = spot_mex_msspoly_check_dimensions(dim,sub,pow,var,coeff)
     *
     *  [errno,v] = ...;
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
     *  102     |  Sub is not sorted.
     *  104     |  A row of var is not sorted or has repeated entries.
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



int generateErrorCode(const mxArray *prhs[])
{
  int canon_err = 0;
  
  /* Test type of each argument. */
  char *sdouble = "double";
  int i;
  
  for(i = 0; i < 5; i++){
    if(!mxIsClass(prhs[i],sdouble))
      return i+1;
  }
  
  const mxArray *mx_dim   = prhs[0];
  const mxArray *mx_sub   = prhs[1];
  const mxArray *mx_pow   = prhs[2];
  const mxArray *mx_var   = prhs[3];
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
  double *coeff = mxGetPr(mx_coeff);
  
  if(!isIntArray(2,dim,TEST_NONNEG))
    return 21;

  int m = dim[0];
  int n = dim[1];

  /* The careful reader will at this point notice it was an
     error to store everything "row-wise" in terms of memory-locality. */

  int prev_i = 0;
  int prev_j = 0;

  int r;
  for(r = 0; r < rows; r++){
    int i = sub[r];
    int j = sub[rows+r];

    if(prev_i > i || (prev_i == i && prev_j > j)){
      /* printf("(%d,%d) (%d,%d)\n",i,j,prev_i,prev_j); */
      canon_err = 102;
    }
    
    if(!isInt(i) || !isInt(j))
      return 22;

    if(i <= 0 || i > m || j <= 0 || j > n)
      return 32;

    prev_i = i;
    prev_j = j;
  }

  if(!isIntArray(rows*cols,var,TEST_NONNEG))
    return 24;

  int j;

  for(i = 0; i < rows; i++){
    for(j = 0; j < cols; j++){
      int idx = j*rows + i;

      if(!isInt(pow[idx]))
	return 23;

      if(pow[i] < 0 && !isTrigId((int)var[idx]))
	return 33;

      if(var[idx] == 0.0 && pow[idx] != 0.0)
	return 34;

      if(j > 0){
	int prevIdx = (j-1)*rows + i;
	if(var[idx] != 0.0 && var[prevIdx] == var[idx]){
	  /* printf("foo: %d %d: %f %f\n",prevIdx,idx,var[prevIdx],var[idx]); */
	  canon_err = 104;
	}

	if(var[prevIdx] < var[idx]){
	  /* printf("%d %d: %f %f\n",prevIdx,idx,var[prevIdx],var[idx]); */
	  canon_err = 104;
	}
      }
    }
  }

  for(i = 0; i < rows; i++){
    if(!mxIsFinite(coeff[i]))
      return 35;
  }

  return canon_err;
}
