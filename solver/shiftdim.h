/*
 * MATLAB Compiler: 3.0
 * Date: Mon Mar 21 17:18:54 2005
 * Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
 * "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on" "-O"
 * "array_indexing:on" "-O" "optimize_conditionals:on" "-v" "-d"
 * "/usr/drfc/cgc/matlab5/zineb/v3.0/solver" "-x" "-W" "mex" "-L" "C" "-t" "-T"
 * "link:mexlibrary" "libmatlbmx.mlib" "rpdederive" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __shiftdim_h
#define __shiftdim_h 1

#ifdef __cplusplus
extern "C" {
#endif

#include "libmatlb.h"

extern void InitializeModule_shiftdim(void);
extern void TerminateModule_shiftdim(void);
extern _mexLocalFunctionTable _local_function_table_shiftdim;

extern mxArray * mlfShiftdim(mxArray * * nshifts, mxArray * x, mxArray * n);
extern void mlxShiftdim(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[]);

#ifdef __cplusplus
}
#endif

#endif
