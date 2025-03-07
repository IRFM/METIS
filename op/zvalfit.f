#ifdef MATLAB73
#include "fintrf.h"
#endif

C-----------------------------------------------------------------------
C syntaxe de la fonction :
C  [F,FAIL] = zvalfit(NX,NY,LAMDA,MU,C,X,Y)
C
C entrees :
C  NX,NY,LAMDA,MU,C  -> matrices des bsplines
C  (X,Y)   -> coordonnees de la grille [M,N] ou il faut evaluer la fonction
C  
C sorties :
C  F                -> valeur des donnees aux points (X,Y) [M,N]
C  FAIl              -> code d'erreur [1,1] (0 = ok)                     
C
C ordre de compilation :
C   mex -v -lnag zvalfit.f
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
      SUBROUTINE mexFunction(NLHS, PLHS, NRHS, PRHS)
      CRONOSINT PLHS(*), PRHS(*)
      INTEGER*4 NLHS, NRHS
	CRONOSINT  mxGetM, mxGetN, mxCreateFull, mxGetPr
C---------------------------------------------------------------------
C  variables intermediaire
C-----------------------------------------------------------------------
C
      INTEGER*4 MIN, NIN, NN
      REAL*8 R

C---------------------------------------------------------------------
C  Constante de dimensionnements des variables
C  MMAX est la valeur maximale de M
C  NESTMAX est la valeur maximale de NXEST et NYEST
C-----------------------------------------------------------------------
C
      INTEGER*4 MMAX, NESTMAX
      PARAMETER (MMAX=5000000,NESTMAX = 10000)

C---------------------------------------------------------------------
C  variable de la fonction de fit
C  MMAX est la valeur maximale de M
C-----------------------------------------------------------------------
C
      INTEGER M,MX,MY, IFAIL ,NX, NY ,i,j
      REAL*8 LAMDA(NESTMAX), MU(NESTMAX)
cWRK(NESTMAX)! ATTENTION MMAX * MMAX trop grand
      REAL*8 X(MMAX),Y(MMAX),F(MMAX),C(MMAX),z(MMAX)
      integer lwrk,kwrk,kx,ky
      parameter (kwrk=36000)
      parameter (lwrk=20000)
      integer iwrk(kwrk)
      real*8 wrk(lwrk)

C---------------------------------------------------------------------
C   gestion des arguments d'entrees
C-----------------------------------------------------------------------
C
      IF (NRHS.NE.7) THEN
        CALL mexErrMsgTxt('Il faut 7 arguments en entree')
      ENDIF
      IF (NLHS.NE.2) THEN
        CALL mexErrMsgTxt('Il faut 2 arguments en sortie')
      ENDIF

      MIN = mxGetM(PRHS(6))
      NIN = mxGetN(PRHS(6))

      M   = MIN * NIN
      IF (M.GT.MMAX) THEN
        CALL mexErrMsgTxt('Trop de points en entree')
      ENDIF

C  securirte avant recopie

      IF (mxGetM(PRHS(7)).NE.MIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatibles entre X et Y')
      ENDIF
      IF (mxGetN(PRHS(7)).NE.NIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatibles entre X et Y')
      ENDIF

C  recopie des donnees matlab dans des tableaux fortran
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(1)),R,1)
      NX = R
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(2)),R,1)
      NY = R
C      WRITE(*,*) NX, NY

      NN  =  mxGetM(PRHS(3))*mxGetN(PRHS(3))
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(3)),LAMDA,NN)

      NN  =  mxGetM(PRHS(4))*mxGetN(PRHS(4))
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(4)),MU,NN)

      NN  =  mxGetM(PRHS(5))*mxGetN(PRHS(5))
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(5)),C,NN)

      CALL mxCopyPtrToReal8(mxGetPr(PRHS(6)),X,M)

      CALL mxCopyPtrToReal8(mxGetPr(PRHS(7)),Y,M)

C---------------------------------------------------------------------
C   appel de la fonction de fit
C-----------------------------------------------------------------------
C
      IFAIL = -1
 
      
      
c      CALL E02DEF(M,NX,NY,X,Y,LAMDA,MU,C,F,WRK,IWRK,IFAIL)
           
c--- remplacement de E02DEF par bispev
      kx = 4
      ky = 3
      MX = MIN
      MY = NIN

      call bispevi(LAMDA,NX,MU,NY,C,kx,ky,X,MX,Y,MY,z,wrk,lwrk,iwrk,
     &             kwrk,ifail)     
c-- reconstruction de F (transposition de la matrice z)
      do 66 i=1,MX
        do 67 j=1,MY
	    F(MX*(j-1)+i)=z(MY*(i-1)+j)
   67 continue     
   66 continue
C---------------------------------------------------------------------
C   creation des donnees de sortie
C-----------------------------------------------------------------------
C
c      WRITE(*,*) IFAIL

      PLHS(1) = mxCreateFull(MIN,NIN,0)
      CALL mxCopyReal8ToPtr(F, mxGetPr(PLHS(1)),M)

      PLHS(2) = mxCreateFull(1,1,0)
      R = 1.0*IFAIL
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(2)), 1)

      RETURN
      END
c
c
c
      SUBROUTINE bispevi(LAMDA,NX,MU,NY,C,kx,ky,X,MX,Y,MY,z,wrk,lwrk,
     &           iwrk,kwrk,ifail) 
c cette routine sert à prendre uniquement la premiere colonne de X et la permiere ligne de Y. 
c En effet, la fonction NAG demandait des matrices MX*MY avec répétition des colonnes pour X et répétition des lignes pour Y
c bispev ne demande qu'un vecteur de taille MX(respectivement MY) pour X(respectivement pour Y)
      
          implicit none
	  integer NX,NY,KX,KY,MX,MY,MMAX,i,j,lwrk,kwrk,ifail,incr
	  parameter (MMAX = 5000000)
	  integer iwrk(kwrk)
c	  real*8 X(MX,MY),Y(MMAX),z(MX,MY)
          real*8 X(MX,MY),Y(MX,MY),z(MX,MY)
	  real*8 LAMDA(NX),MU(NY),C(*),wrk(lwrk)
  

         call bispev(LAMDA,NX,MU,NY,C,kx,ky,X(:,1),MX,Y(1,:),MY,z,wrk,
     &                lwrk,iwrk,kwrk,ifail) 
      END
      
c--- definition de bispev      
      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
     * iwrk,kwrk,ier)
c  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
c  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
c  b-spline representation.
c
c  calling sequence:
c     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
c    * iwrk,kwrk,ier)
c
c  input parameters:
c   tx    : real array, length nx, which contains the position of the
c           knots in the x-direction.
c   nx    : integer, giving the total number of knots in the x-direction
c   ty    : real array, length ny, which contains the position of the
c           knots in the y-direction.
c   ny    : integer, giving the total number of knots in the y-direction
c   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
c           b-spline coefficients.
c   kx,ky : integer values, giving the degrees of the spline.
c   x     : real array of dimension (mx).
c           before entry x(i) must be set to the x co-ordinate of the
c           i-th grid point along the x-axis.
c           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
c   mx    : on entry mx must specify the number of grid points along
c           the x-axis. mx >=1.
c   y     : real array of dimension (my).
c           before entry y(j) must be set to the y co-ordinate of the
c           j-th grid point along the y-axis.
c           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
c   my    : on entry my must specify the number of grid points along
c           the y-axis. my >=1.
c   wrk   : real array of dimension lwrk. used as workspace.
c   lwrk  : integer, specifying the dimension of wrk.
c           lwrk >= mx*(kx+1)+my*(ky+1)
c   iwrk  : integer array of dimension kwrk. used as workspace.
c   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
c
c  output parameters:
c   z     : real array of dimension (mx*my).
c           on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
c           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
c   ier   : integer error flag
c    ier=0 : normal return
c    ier=10: invalid input data (see restrictions)
c
c  restrictions:
c   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
c   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
c   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
c
c  other subroutines required:
c    fpbisp,fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
c  ..array arguments..
      integer iwrk(kwrk)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
c  ..local scalars..
      integer i,iw,lwest
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
  10  do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  if(my-1) 100,60,40
  40  do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
  50  continue
  60  ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),
     * iwrk(1),iwrk(mx+1))
 100  return
      end
      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
c  ..scalar arguments..
      integer nx,ny,kx,ky,mx,my
c  ..array arguments..
      integer lx(mx),ly(my)
      real tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wx(mx,kx+1),wy(my,ky+1)
c  ..local scalars..
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real arg,sp,tb,te
c  ..local arrays..
      real h(6)
c  ..subroutine references..
c    fpbspl
c  ..
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  10    if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
  20    call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
  30    continue
  40  continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
  50    if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
  60    call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
  70    continue
  80  continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
  90    continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100        continue
            l1 = l1+nky1
 110      continue
          m = m+1
          z(m) = sp
 120    continue
 130  continue
      return
      end

      subroutine fpbspl(t,n,k,x,l,h)
c  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
c  degree k at t(l) <= x < t(l+1) using the stable recurrence
c  relation of de boor and cox.
c  ..
c  ..scalar arguments..
      real x
      integer n,k,l
c  ..array arguments..
      real t(n),h(6)
c  ..local scalars..
      real f,one
      integer i,j,li,lj
c  ..local arrays..
      real hh(5)
c  ..
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
  10    continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
  20  continue
      return
      end
