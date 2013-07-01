#ifndef __SSID_MEX_HELPERS__
#define __SSID_MEX_HELPERS__
#include <mex.h>
void getArgSized(mxArray **ptr, char * class, int idx,int nrhs, const mxArray *prhs[], int m, int n);

void getScalarDArgDefault(double *dst,int idx,int nrhs, const mxArray *prhs[],double def);


double poly_eval(int d, double *coeff, double x);

#endif /*  __SSID_MEX_HELPERS__ */
