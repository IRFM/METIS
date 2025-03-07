/*
bspbase

computes the B-spline base functions at given coordinates for a
given knot sequence.

void bspbase(B,T,sizeT,k,X,sizeX)

double *B; computed B-spline base functions, stored consecutively
  for each value of X; there are sizeT - k base functions
double *T; the sorted knot sequence
int sizeT;
int k; the spline order (degree + 1); sizeT - k > 0
double *X; the values where to compute B; X must be in the interval
  [T(k),T(last-k+1)] 
*/

void bspbase(B,T,sizeT,k,X,sizeX)

double *B;
double *T;
int sizeT;
int k;
double *X;
int sizeX;

{
 register int i, j, r, s;
 register double *b, *ti, *tj, o;
 int  nbase;

 nbase = sizeT - k;
 j = k - 1; /* the index in T of the larger knot <= X[i] */
 for (i = 0; i < sizeX; i++) {
  /* find j, note that j == sizeT-1 produces indexing problem */
  for (; j < sizeT-k; j++) if (X[i] >= T[j] && X[i] <= T[j+1]) break;
  /* where in B the last of the k non zero base functions must be stored */
  b = B + j + nbase*i;
  /* recurrence */
  b[0] = 1;
  for (r = -1; r > -k; r--) {
   ti = (tj = T + j) + r;
   for (s = r; s < 0; s++) {
    o = *(++tj) - *(++ti);
    o = o == 0 ? 0 : (X[i] - *ti) / o;
    b[s] += (1 - o) * b[s+1];
    b[s+1] *= o;
   }
  }
 }
}
