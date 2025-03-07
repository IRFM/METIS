C-----------------------------------------------------------------------
C syntaxe de la fonction :
C  [FF,DFFDX,DFFDY,FAIL] = zinterp2d(X,Y,F,XX,YY)
C
C entrees :
C  (X,Y)   -> coordonnees des donnees d'entrees [M,N] 
C  F       -> valeurs des donnees d'entrees [M,N] 
C  (XX,YY)   -> coordonnees des points de sorties [MM,NN] 
C  
C sorties :
C  FF                -> donnees interpolees
C  DFFDX             -> derivee /a X des donnees interpolees
C  DFFDY             -> derivee /a Y des donnees interpolees
C  FAIl              -> code d'erreur [1,1] (0 = ok)                     
C
C  remrque : si FF, DFFDX ou DFFDY > 1e308 -> NaN
C
C ordre de compilation :
C   mex -v -lnag zinterp2d_hel.f
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
      SUBROUTINE mexFunction(NLHS, PLHS, NRHS, PRHS)
      INTEGER PLHS(*), PRHS(*)
		integer  mxGetM, mxGetN, mxCreateFull, mxGetPr
      INTEGER*4 NLHS, NRHS
C---------------------------------------------------------------------
C  variables intermediaire
C-----------------------------------------------------------------------
C
      INTEGER*4 MIN, NIN, MOUT, NOUT
      REAL*8 R

C---------------------------------------------------------------------
C  Constante de dimensionnements des variables
C  MMAX est la valeur maximale de M
C-----------------------------------------------------------------------
C
      INTEGER*4 MMAX, LIQ, LRQ
      PARAMETER (MMAX=30000)
      PARAMETER (LIQ = 2 * MMAX + 1)
      PARAMETER (LRQ = 6  *MMAX + 1)
C---------------------------------------------------------------------
C  variable de la fonction de fit
C  MMAX est la valeur maximale de M
C-----------------------------------------------------------------------
C
      INTEGER*4 M, N, NW, NQ, IQ(LIQ), IFAIL
      REAL*8 X(MMAX), Y(MMAX), F(MMAX), XX(MMAX), YY(MMAX), FF(MMAX) 
      REAL*8 RQ(LRQ), DFFDX(MMAX) , DFFDY(MMAX)

C---------------------------------------------------------------------
C   gestion des arguments d'entrees
C-----------------------------------------------------------------------
C
      IF (NRHS.NE.5) THEN
        CALL mexErrMsgTxt('Il faut 5 arguments en entree')
      ENDIF
      IF (NLHS.NE.4) THEN
        CALL mexErrMsgTxt('Il faut 4 arguments en sortie')
      ENDIF

      MIN = mxGetM(PRHS(1))
      NIN = mxGetN(PRHS(1))

      M   = MIN * NIN
      IF (M.GT.MMAX) THEN
        CALL mexErrMsgTxt('Trop de points en entree')
      ENDIF

      MOUT = mxGetM(PRHS(4))
      NOUT = mxGetN(PRHS(4))

      N   = MOUT * NOUT
      IF (N.GT.MMAX) THEN
        CALL mexErrMsgTxt('Trop de points en sortie')
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
      IF (mxGetM(PRHS(5)).NE.MOUT) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre XX et YY')
      ENDIF
      IF (mxGetN(PRHS(5)).NE.NOUT) THEN
        CALL mexErrMsgTxt('Dimensions incompatible entre XX et YY')
      ENDIF

C  recopie des donnees matlab dans des tableaux fortran
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(1)),X,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(2)),Y,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(3)),F,M)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(4)),XX,N)
      CALL mxCopyPtrToReal8(mxGetPr(PRHS(5)),YY,N)

C---------------------------------------------------------------------
C   appel de la fonction de fit
C-----------------------------------------------------------------------
C
 
      IFAIL = -1
      NQ    = 40
      NW    = 40 
      CALL E01SGF(M,X,Y,F,NW,NQ,IQ,LIQ,RQ,LRQ,IFAIL)

      IF (IFAIL.NE.0) THEN
         PLHS(1) = mxCreateFull(0,0,0)
         PLHS(2) = mxCreateFull(0,0,0)
         PLHS(3) = mxCreateFull(0,0,0)
         PLHS(4) = mxCreateFull(1,1,0)
         R = 1.0*IFAIL
         CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(4)), 1)
         RETURN
      ENDIF
C---------------------------------------------------------------------
C   appel de la fonction d'evaluation
C-----------------------------------------------------------------------
C
 
      IFAIL = 1
      CALL E01SHF(M,X,Y,F,IQ,LIQ,RQ,LRQ,N,XX,YY,FF,DFFDX,DFFDY,IFAIL)

C---------------------------------------------------------------------
C   creation des donnees de sortie
C-----------------------------------------------------------------------
C
      PLHS(1) = mxCreateFull(MOUT,NOUT,0)
      CALL mxCopyReal8ToPtr(FF, mxGetPr(PLHS(1)), N)

      PLHS(2) = mxCreateFull(MOUT,NOUT,0)
      CALL mxCopyReal8ToPtr(DFFDX, mxGetPr(PLHS(2)), N)

      PLHS(3) = mxCreateFull(MOUT,NOUT,0)
      CALL mxCopyReal8ToPtr(DFFDY, mxGetPr(PLHS(3)), N)

      PLHS(4) = mxCreateFull(1,1,0)
      R = 100.0*IFAIL
      CALL mxCopyReal8ToPtr(R, mxGetPr(PLHS(4)), 1)

      RETURN
      END
