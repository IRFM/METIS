/*
BSPSUM   computes the combination or the derivative of a combination of B-spline
         base functions.

void bspsum(V,T,nT,C,nC,X,nX,d,w)

double	*V;	the values of the combination for each X
double	*T;	the sorted knot sequence
int	nT;
double	*C;	the coefficient of each base funtion; the spline order is
		k = nT - nC > 0
int	nC;
double	*X;	the values where V must be computed; X must be in the interval
		[T(k),T(last-k+1)]
int	nX;
int	d;	the derivative order
double	*w;	work area of n k

*/

void bspsum(V,T,nT,C,nC,X,nX,d,w)

double	*V;
double	*T;
int	nT;
double	*C;
int	nC;
double	*X;
int	nX;
int	d;
double	*w;

{
 register int i, jl, jh, jn, r, s, k, kk;
 register double *ti, *tj, o, *Xi;

 k = nT - nC;
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
}
