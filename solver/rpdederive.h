/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu May 19 15:32:59 2005
 * Arguments: "-d" "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/solver" "-x"
 * "rpdederive" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __rpdederive_h
#define __rpdederive_h 1

#include "matlab.h"

extern mxArray * mlfRpdederive(mxArray * x,
                               mxArray * y,
                               mxArray * c0,
                               mxArray * c1,
                               mxArray * d,
                               mxArray * o);
extern void mlxRpdederive(int nlhs,
                          mxArray * plhs[],
                          int nrhs,
                          mxArray * prhs[]);

#endif
