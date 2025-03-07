# include "mex.h"
# include <math.h>

# define max(a,b)	((a) > (b) ? (a) : (b))
# define min(a,b)	((a) < (b) ? (a) : (b))
# define sqr(a)   ((a)*(a))

void inter2d(double *s1,double *s2, double *s3, double *s4,
	double rd, double zd, double pvd, double tvd,
	double r1, double r2, double z1, double z2)
		
{
 double sinpvd, cospvd, costvd, t1, t2, rdc, l1, l2;

	/* default values */
	*s2 = mxGetInf(); /* other are assumed to be 0 already */

	/* non crossing horizontal or vertical chord */
	sinpvd = sin(pvd);
	cospvd = cos(pvd);
	if ((sinpvd == 0 && (zd >= z2 || zd < z1)) ||
	  	 (cospvd == 0 && (rd >= r2 || rd < r1))) {
		*s2 = -1;
 	return;
	}

	/* z boundary */
	if (sinpvd > 0) {
		l1 = (z1 - zd) / sinpvd;
		l2 = (z2 - zd) / sinpvd;
 	*s1 = max(*s1,l1);
		*s2 = min(*s2,l2);
	} else if (sinpvd < 0) {
		l1 = (z2 - zd) / sinpvd;
		l2 = (z1 - zd) / sinpvd;
 	*s1 = max(*s1,l1);
		*s2 = min(*s2,l2);
	}
	if (*s2 <= *s1) return;

	/* r2 boundary */
	costvd = cos(tvd);
	t2 = sqr(rd) * (sqr(costvd) - 1);
	t1 = t2 + sqr(r1);
	t2 = t2 + sqr(r2);
	if (t2 < 0) {
		*s2 = -1;
		return;
	} else if (cospvd != 0) {
 	rdc = rd * costvd / cospvd;
 	t2 = sqrt(t2) / fabs(cospvd);
		l1 = rdc - t2;
		l2 = rdc + t2;
		*s1 = max(*s1,l1);
		*s2 = min(*s2,l2);
	}
	if (*s2 <= *s1) return;

	/* r1 boundary */
	if (t1 >=0 & cospvd != 0) {
		l1 = sqrt(t1) / fabs(cospvd);
		l2 = rdc + l1; l1 = rdc - l1;
		if (l1 <= *s1)
			*s1 = max(*s1,l2);
	 else {
			if (l2 <= *s2) {
 			*s3 = l2;
 	 	*s4 = *s2;
			}
			*s2 = min(*s2,l1);
		}
	}
}

/* s = intersect2dmex(n,rd,zd,pvd,tvd,r1,r2,z1,z2) */

# define S	 plhs[0]
# define RD prhs[1]
# define ZD prhs[2]
# define PVD prhs[3]
# define TVD prhs[4]
# define R1 prhs[5]
# define R2 prhs[6]
# define Z1 prhs[7]
# define Z2 prhs[8]

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	double *s,*rd,*zd,*pvd,*tvd,r1,r2,z1,z2;
	int k,n;
	n = mxGetScalar(prhs[0]);
	S = mxCreateDoubleMatrix(4,n,mxREAL); s = mxGetPr(S);
	rd = mxGetPr(RD);
	zd = mxGetPr(ZD);
	pvd = mxGetPr(PVD);
	tvd = mxGetPr(TVD);
	r1 = mxGetScalar(R1);
	z1 = mxGetScalar(Z1);
	r2 = mxGetScalar(R2);
	z2 = mxGetScalar(Z2);
	for (k = 0; k < n; k++, s += 4)
		inter2d(s,s+1,s+2,s+3,rd[k],zd[k],pvd[k],tvd[k],r1,r2,z1,z2);
}
		


