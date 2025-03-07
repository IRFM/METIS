/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu May 19 15:32:55 2005
 * Arguments: "-d" "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op" "-x" "zregular" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __zregular_h
#define __zregular_h 1

#include "matlab.h"

extern mxArray * mlfZregular(mxArray * * a,
                             mxArray * * b,
                             mxArray * * c,
                             mxArray * * d,
                             mxArray * x,
                             mxArray * y_);
extern void mlxZregular(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
