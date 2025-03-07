# include "mex.h"

extern void bspbase();

# define T  prhs[0]
# define K  prhs[1]
# define X  prhs[2]
# define B  plhs[0]

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
 int		k, sizeX, sizeT, nbase;
 int		*d;
 register int	i;

 k = mxGetScalar(K);
 sizeX = mxGetM(X) * mxGetN(X);
 sizeT = mxGetM(T) * mxGetN(T);
 nbase = sizeT - k;

 B = mxCreateDoubleMatrix(nbase, sizeX, mxREAL);
 bspbase(mxGetPr(B), mxGetPr(T), sizeT, k, mxGetPr(X), sizeX);
}
