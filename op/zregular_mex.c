/*
 * MATLAB Compiler: 2.0.1
 * Date: Thu May 19 15:32:55 2005
 * Arguments: "-d" "/home/sccp/gsem/cgc/matlab5/zineb/v3.0/op" "-x" "zregular" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "matlab.h"
#include "zregular.h"

static mlfFunctionTableEntry function_table[1]
  = { { "zregular", mlxZregular, 2, 5 } };

/*
 * The function "mexFunction" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxZregular". Finally, it clears the feval table and exits.
 */
void mexFunction(int nlhs, mxArray * * plhs, int nrhs, mxArray * * prhs) {
    mlfTry {
        mlfFunctionTableSetup(1, function_table);
        mclImportGlobal(0, NULL);
        mlxZregular(nlhs, plhs, nrhs, prhs);
        mlfFunctionTableTakedown(1, function_table);
    } mlfCatch {
        mlfFunctionTableTakedown(1, function_table);
        mclMexError();
    } mlfEndCatch
}
