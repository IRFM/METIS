#ifdef MATLAB73
#include "fintrf.h"
#endif

C-----------------------------------------------------------------------
C syntaxe de la fonction :
C  [NX,NY,LAMDA,MU,C,FP,DEFRANK,FAIL] = zfitdata1(X,Y,F,W,S,MINMAX,ndimtab)
C
C entrees :
C  (X,Y)   -> coordonnees de la grille [M,N] 
C  F       -> valeurs sur la grilles [M,N] 
C  W       -> poids de chaque donnee [M,N] 
C  S       -> lissage  > 0 [1,1]
C  MINMAX =[XMIN,XMAX,YMIN,YMAX] -> zone de validite du fit 
C  ndimtab -> dimension du tableau d'entree (ndimtab(1)*ndimtab(2) = M*N)
c  si l'on donne en entree un vecteur
C sorties :
C  NX,NY,LAMDA,MU,C  -> matrices des bsplines
C  FP                -> residu [1,1]
C  DEFRANK           -> deficite du rang du systeme [1,1]
C                       ((NX-4)*(NY-4) -RANK )
C  FAIl              -> code d'erreur [1,1] (0 = ok)                     
C
C ordre de compilation :
C   mex -v -lnag zfitdata1.f
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
      SUBROUTINE mexFunction(NLHS, PLHS, NRHS, PRHS)
      CRONOSINT PLHS(*), PRHS(*)
		CRONOSINT  mxGetM, mxGetN, mxCreateFull, mxGetPr
      INTEGER*4 NLHS, NRHS
C---------------------------------------------------------------------
C  variables intermediaire
C-----------------------------------------------------------------------
C
      INTEGER*4 MIN, NIN, MV
      REAL*8 R, MINMAX(4), ndimtab(2)

C---------------------------------------------------------------------
C  Constante de dimensionnements des variables
C  MMAX est la valeur maximale de M
C  NESTMAX est la valeur maximale de NXEST et NYEST
C-----------------------------------------------------------------------
C
      INTEGER*4 MMAX, NESTMAX, LWRK, LIWRK 
      PARAMETER (MMAX=100000,NESTMAX = 8)
      PARAMETER (LWRK=MMAX*NESTMAX)
      PARAMETER (LIWRK=MMAX+2*(NESTMAX-7)*(NESTMAX-7))
C---------------------------------------------------------------------
C  variable de la fonction de fit
C  MMAX est la valeur maximale de M
C-----------------------------------------------------------------------
C
      CHARACTER START
      INTEGER*4 M, NXEST, NYEST, NX, NY, RANK, IFAIL, L , IWRK(LIWRK)
      REAL*8 LAMDA(NESTMAX), MU(NESTMAX)
      REAL*8 X(MMAX), Y(MMAX), F(MMAX), W(MMAX), C(MMAX), WRK(LWRK)
      REAL*8 S, FP

C---------------------------------------------------------------------
C   gestion des arguments d'entrees
C-----------------------------------------------------------------------
C
      IF (NRHS.LT.6) THEN
        CALL mexErrMsgTxt('Il faut 6 arguments en entree')
      ENDIF
      IF (NLHS.NE.8) THEN
        CALL mexErrMsgTxt('Il faut 8 arguments en sortie')
      ENDIF
      
      MV     = mxGetM(PRHS(7)) * mxGetN(PRHS(7))
      write(6,*) 'MV=',MV
      IF (MV.NE.2) THEN
             CALL mexErrMsgTxt('ndimtab doit avoir 2 elements')
      ENDIF
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(7)),ndimtab,2)
      NXEST  = ndimtab(1)
      NYEST  = ndimtab(2)
      write(6,*) 'ndimtab',ndimtab

      MIN = mxGetM(PRHS(1))
      NIN = mxGetN(PRHS(1))

      write(6,*) MIN,NIN,NESTMAX
      M   = MIN * NIN
      MV  = M + 2
      IF (MV.GT.MMAX) THEN
        CALL mexErrMsgTxt('Trop de points en entree')
      ENDIF

C  securirte avant recopie

      IF (mxGetM(PRHS(2)).NE.MIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et Y')
      ENDIF
      IF (mxGetN(PRHS(2)).NE.NIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et Y')
      ENDIF
      IF (mxGetM(PRHS(3)).NE.MIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et F')
      ENDIF
      IF (mxGetN(PRHS(3)).NE.NIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et F')
      ENDIF
      IF (mxGetM(PRHS(4)).NE.MIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et W')
      ENDIF
      IF (mxGetN(PRHS(4)).NE.NIN) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre X et W')
      ENDIF
      IF (mxGetM(PRHS(5)).LT.1) THEN
        CALL mexErrMsgTxt('Il faut donner S')
      ENDIF
      IF (mxGetN(PRHS(5)).LT.1) THEN
        CALL mexErrMsgTxt('Il faut donner S')
      ENDIF

C  recopie des donnees matlab dans des tableaux fortran
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(1)),X,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(2)),Y,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(3)),F,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(4)),W,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(5)),S,1)

      IF (NRHS.EQ.6) THEN
          MV = mxGetM(PRHS(6)) * mxGetN(PRHS(6))
          IF (MV.NE.4) THEN
             CALL mexErrMsgTxt('MINMAX doit avoir 4 elements')
          ENDIF
          CALL mxCopyPtrToReal8(mxGetPr(PRHS(6)),MINMAX,4)
          X(M+1) = MINMAX(1)
          X(M+2) = MINMAX(2)
          Y(M+1) = MINMAX(3)
          Y(M+2) = MINMAX(4)
          F(M+1) = 0.0
          F(M+2) = 0.0
          W(M+1) = 0.0
          W(M+2) = 0.0
          MV     = M + 2
      ELSE
          MV = M
      ENDIF

      
C---------------------------------------------------------------------
C   appel de la fonction de fit
C-----------------------------------------------------------------------
C
 
      IFAIL = -1
      START  = 'Cold'
      CALL E02DDF(START,MV,X,Y,F,W,S,NXEST,NYEST,NX,LAMDA,NY,MU,C,FP,
     #             RANK,WRK,LWRK,IWRK,LIWRK,IFAIL)

C      WRITE(*,*) 'Nx et Ny :', NX, NY
C---------------------------------------------------------------------
C   creation des donnees de sortie
C-----------------------------------------------------------------------
C
      PLHS(1) = mxCreateFull(1,1,0)
      R = 1.0*NX
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(1)), 1)

      PLHS(2) = mxCreateFull(1,1,0)
      R = 1.0*NY
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(2)), 1)

      PLHS(3) = mxCreateFull(NXEST,1,0)
      CALL mxCopyReal8ToPtr(LAMDA, mxGetPr(PLHS(3)), NXEST)

      PLHS(4) = mxCreateFull(NYEST,1,0)
      CALL mxCopyReal8ToPtr(MU, mxGetPr(PLHS(4)), NYEST)

      PLHS(5) = mxCreateFull(NXEST-4,NYEST-4,0)
      CALL mxCopyReal8ToPtr(C, mxGetPr(PLHS(5)), (NXEST-4)*(NYEST-4))

      PLHS(6) = mxCreateFull(1,1,0)
      CALL mxCopyReal8ToPtr(FP, mxGetPr(PLHS(6)),1)

      PLHS(7) = mxCreateFull(1,1,0)
      R = 1.0*((NX-4)*(NY-4)-RANK)
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(7)),1)

      PLHS(8) = mxCreateFull(1,1,0)
      R = 1.0*IFAIL
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(8)), 1)

      RETURN
      END
