# include "mex.h"
# include "math.h"
/* Function from J-M Moret toolbox  (EPFL/CRPP)*/
 
double ellipfast(double y, int kind);

# define M	prhs[0]
# define KIND	prhs[1]
# define F	plhs[0]

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
 int kind;
 register int k;
 register double *m, *f;

 F = mxCreateDoubleMatrix(mxGetM(M), mxGetN(M), mxREAL);
 kind = mxGetScalar(KIND);
 m = mxGetPr(M);
 f = mxGetPr(F);
 for (k = mxGetM(M)*mxGetN(M); k--;)
  *f++ = ellipfast(*m++,kind);
}


/* ellipfast.c */
double ellipfast(double y, int kind)

{ y = 1 - y;
  return(
   y < 0 || y > 1 ?
     log(-1.0) 
     /* NaN */
   :
     kind == 1 ?
       (((0.01451196212 * y +
	  0.03742563713) * y +
	  0.03590092383) * y +
	  0.09666344259) * y +
	  1.38629436112 -
      ((((0.00441787012 * y +
	  0.03328355346) * y +
	  0.06880248576) * y +
	  0.12498593597) * y +
	  0.50000000000) * log(y)
     : /* kind == 2 */
       y == 0 ?
	 1
       :
	 (((0.01736506451 * y +
	    0.04757383546) * y +
	    0.06260601220) * y +
	    0.44325141463) * y +
	    1.00000000000 -
	((((0.00526449639 * y +
	    0.04069697526) * y +
	    0.09200180037) * y +
	    0.24998368310) * y) * log(y));
}
