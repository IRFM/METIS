# include "mex.h"
# ifdef BSPNTH
# include <pthread.h>
# endif

/* void bspsum(V,T,nT,C,nC,X,nX,d)

computes the combination or the derivative of a combination of B-spline
base functions.

double	*V;	the values of the combination for each X
double	*T;	the sorted knot sequence
int	nT;
double	*C;	the coefficient of each base funtion; the spline order is
int	nC;		  k = nT - nC > 0
double	*X;	the values where V must be computed; X must be in the interval
int	nX;		  [T(k),T(end-k+1)]
int	d;	    the derivative order
*/

void bspsum(double	*V, double	*T, int	nT, double	*C, int	nC,
            double	*X, int	nX, int	d)
{
 int i, jl, jh, jn, r, s, k, kk;
 double *ti, *tj, o, *Xi, *w;

 k = nT - nC;
 w = (double *)malloc(k*sizeof(double));
 jl = k - 1; jh = nT - k;
 for (i = 0, Xi = X; i < nX; i++, Xi++) {
  /* find j, the index in T of the larger knot <= X[i]
     note that j == nT-1 produces indexing problem */
  if (*Xi < T[jl]) jl = k - 1;
  if (*Xi > T[jh]) jh = nT - k;
  while (jh - jl > 1) {
   jn = (jl + jh) / 2;
   if (*Xi > T[jn])
    jl = jn;
   else
    jh = jn;
  }
  /* copy c in w */
  for (s=0; s < k; s++) w[s] = C[jl-s];
  /* perform derivation */
  kk = k;
  for (r = 0; r < d; r++) {
   tj = (ti = T + jl) + --kk;
   for (s = 0; s < kk; s++) {
    o = *(tj--) - *(ti--);
    w[s] = o == 0 ? 0 : (w[s] - w[s+1]) * kk / o;
   }
  }
  /* do the evaluation */
  while (kk > 1) {
   tj = (ti = T + jl) + --kk;
   for (s = 0; s < kk; s++) {
    o = *(tj--) - *ti;
    if (o != 0) o = (*Xi - *(ti--)) / o;
    w[s] = w[s] * o + (1 - o) * w[s+1];
   }
  }
  V[i] = w[0];
 }
 free(w);
}

# ifdef BSPSUM4
void bspsum40(double	*V, double	*T, int	nT, double	*C, int dum1, double	*X, int	nX, int d)

{
 int i, jl, jh, k;
 double	c0, c1, c2, c3, t0, t1, t2, t3, t4, t5, o, Xi, Xit3, Xit4;

 jh = 0; jl = nT - 1;
 for (i = 0; i < nX; i++) {
  Xi = X[i];
  /* find j, the index in T of the larger knot <= X[i] */
  if (Xi > T[jh] || Xi < T[jl]) {
   jh = nT - 4;
   jl = 3;
   while (jh - jl > 1) {
    k = (jl + jh) / 2;
    if (Xi > T[k]) jl = k; else jh = k;
   }
   t0 = T[jl+3];
   t1 = T[jl+2];
   t2 = T[jl+1];
   t3 = T[jl];
   t4 = T[jl-1];
   t5 = T[jl-2];
  }
  c0 = C[jl];
  c1 = C[jl-1];
  c2 = C[jl-2];
  c3 = C[jl-3];
  
  Xit3 = Xi - t3;
  Xit4 = Xi - t4;

# ifdef BSPSUM44 
  if (d == 0) {
   o    =  Xit3     / (t0 - t3);
   c0   = (c0 - c1) * o + c1;
   o    =  Xit4     / (t1 - t4);
   c1   = (c1 - c2) * o + c2;
   o    = (Xi - t5) / (t2 - t5);
   c2   = (c2 - c3) * o + c3;
  } else {
   c0 = (c0 - c1) * 3 / (t0 - t3);
   c1 = (c1 - c2) * 3 / (t1 - t4);
   c2 = (c2 - c3) * 3 / (t2 - t5);
  }
  if (d == 2) {
   c0 = (c0 - c1) * 2 / (t1 - t3);
   c1 = (c1 - c2) * 2 / (t2 - t4);
  } else {
   o    =  Xit3     / (t1 - t3);
   c0   = (c0 - c1) * o + c1;
   o    =  Xit4     / (t2 - t4);
   c1   = (c1 - c2) * o + c2;
  }
  o    =  Xit3     / (t2 - t3);
  V[i] = (c0 - c1) * o + c1;
# else

  o    =  Xit3     / (t0 - t3);
  c0   = (c0 - c1) * o + c1;
  o    =  Xit4     / (t1 - t4);
  c1   = (c1 - c2) * o + c2;
  o    = (Xi - t5) / (t2 - t5);
  c2   = (c2 - c3) * o + c3;

  o    =  Xit3     / (t1 - t3);
  c0   = (c0 - c1) * o + c1;
  o    =  Xit4     / (t2 - t4);
  c1   = (c1 - c2) * o + c2;

  o    =  Xit3     / (t2 - t3);
  V[i] = (c0 - c1) * o + c1;
# endif
 }
}

void bspsum41(double	*V, double	*T, int	nT, double	*C, int dum1, double	*X, int	nX, int d)

{
 int i, jl, jh, k;
 double	c0, c1, c2, c3, t0, t1, t2, t3, t4, t5, o, Xi, Xit3;

 jh = 0; jl = nT - 1;
 for (i = 0; i < nX; i++) {
  Xi = X[i];
  /* find j, the index in T of the larger knot <= X[i] */
  if (Xi > T[jh] || Xi < T[jl]) {
   jh = nT - 4;
   jl = 3;
   while (jh - jl > 1) {
    k = (jl + jh) / 2;
    if (Xi > T[k]) jl = k; else jh = k;
   }
   t0 = T[jl+3];
   t1 = T[jl+2];
   t2 = T[jl+1];
   t3 = T[jl];
   t4 = T[jl-1];
   t5 = T[jl-2];
  }
  c0 = C[jl];
  c1 = C[jl-1];
  c2 = C[jl-2];
  c3 = C[jl-3];
  
  Xit3 = Xi - t3;

  c0 = (c0 - c1) * 3 / (t0 - t3);
  c1 = (c1 - c2) * 3 / (t1 - t4);
  c2 = (c2 - c3) * 3 / (t2 - t5);
  
  o    =  Xit3     / (t1 - t3);
  c0   = (c0 - c1) * o + c1;
  o    = (Xi - t4) / (t2 - t4);
  c1   = (c1 - c2) * o + c2;
  
  o    =  Xit3     / (t2 - t3);
  V[i] = (c0 - c1) * o + c1;
 }
}

void bspsum42(double	*V, double	*T, int	nT, double	*C, int dum1, double	*X, int	nX, int d)

{
 int i, jl, jh, k;
 double	c0, c1, c2, c3, t0, t1, t2, t3, t4, t5, o, Xi;

 jh = 0; jl = nT - 1;
 for (i = 0; i < nX; i++) {
  Xi = X[i];
  /* find j, the index in T of the larger knot <= X[i] */
  if (Xi > T[jh] || Xi < T[jl]) {
   jh = nT - 4;
   jl = 3;
   while (jh - jl > 1) {
    k = (jl + jh) / 2;
    if (Xi > T[k]) jl = k; else jh = k;
   }
   t0 = T[jl+3];
   t1 = T[jl+2];
   t2 = T[jl+1];
   t3 = T[jl];
   t4 = T[jl-1];
   t5 = T[jl-2];
  }
  c0 = C[jl];
  c1 = C[jl-1];
  c2 = C[jl-2];
  c3 = C[jl-3];

  c0 = (c0 - c1) * 3 / (t0 - t3);
  c1 = (c1 - c2) * 3 / (t1 - t4);
  c2 = (c2 - c3) * 3 / (t2 - t5);
  
  c0 = (c0 - c1) * 2 / (t1 - t3);
  c1 = (c1 - c2) * 2 / (t2 - t4);
  
  o    = (Xi - t3) / (t2 - t3);
  V[i] = (c0 - c1) * o + c1;
 }
}
#endif

# ifdef BSPNTH
struct thstruct {
 pthread_t th;
 void (*f)();
 int kth,i1,i2,p,lT,mC,lX,d,Vinc,Xinc;
 double *V,*T,*C,*X;
};

void *thdo(void *thdum)
{
 struct thstruct *th;
 int i;
 th = (struct thstruct *)thdum;
 for (i = th->i1; i < th->i2; i++)
  th->f(th->V + i*th->Vinc, th->T, th->lT, th->C + i*th->mC, th->mC, 
        th->X + (i*th->Xinc)%th->lX, th->Vinc, th->d);
 printf("thdo #%d : %d to %d\n",th->kth,th->i1,th->i2);
 pthread_exit(NULL);
} 
# endif

# define MXT  prhs[0]
# define MXC  prhs[1]
# define MXNC prhs[2]
# define MXX  prhs[3]
# define MXLX prhs[4]
# define MXD  prhs[5]
# define MXP  prhs[6]
# define MXV  plhs[0]

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
 int		lT, mC, nC, lX, k, p, d, i, Vinc, Xinc;
# ifdef BSPNTH
 struct thstruct th[BSPNTH];
 int kth;
# endif
 void (*f)();

 lT = mxGetM(MXT)*mxGetN(MXT);
 mC = mxGetM(MXC);
 nC = mxGetScalar(MXNC);
 lX = mxGetScalar(MXLX);
 d = mxGetScalar(MXD);
 p = mxGetScalar(MXP);
 k = lT - mC;

 /* memory allocation */
 MXV = mxCreateDoubleMatrix(lX, p ? nC/lX : nC, mxREAL);
 
 /* do the work */
 if (d >= k) return;
 
 Vinc = p ? 1 : lX;
 Xinc = p ? 1 : 0;
 f = &bspsum;
# ifdef BSPSUM4
 if (k == 4 && d <= 2) f = d ? (d == 2 ? &bspsum42 : &bspsum41) : &bspsum40;
# endif
# ifndef BSPNTH
# ifdef BSPPAR
# pragma ibm parallel_loop
# endif
  for (i = 0; i < nC; i++)
   f(mxGetPr(MXV) + i*Vinc, mxGetPr(MXT), lT, mxGetPr(MXC) + i*mC, mC,
     mxGetPr(MXX) + (i*Xinc)%lX, Vinc, d);
# else
  th[0].f    = f;
  th[0].p    = p;
  th[0].lT   = lT;
  th[0].mC   = mC;
  th[0].lX   = lX;
  th[0].d    = d;
  th[0].Vinc = Vinc;
  th[0].Xinc = Xinc;
  th[0].V    = mxGetPr(MXV);
  th[0].T    = mxGetPr(MXT);
  th[0].C    = mxGetPr(MXC);
  th[0].X    = mxGetPr(MXX);
  for (kth = 0; kth < BSPNTH; kth++) {
   th[kth] = th[0];
   th[kth].i1 = (kth*nC)/BSPNTH;
   th[kth].i2 = ((kth+1)*nC)/BSPNTH;
   th[kth].kth = kth;
   pthread_create(&th[kth].th,NULL,thdo,th+kth);
  }
  for (kth = 0; kth < BSPNTH; kth++)
   pthread_join(th[kth].th,NULL);
# endif
}
