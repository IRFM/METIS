/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu Jun 14 15:51:44 2001
 * Arguments: "-x" "pde1dsolver.m" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __squeeze_h
#define __squeeze_h 1

#include "matlab.h"

extern mxArray * mlfSqueeze(mxArray * a);
extern void mlxSqueeze(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#endif
