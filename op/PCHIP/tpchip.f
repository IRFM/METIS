#ifdef MATLAB73
#include "fintrf.h"
#endif
      subroutine mexFunction(nlhs,plhs,nrhs,prhs)
C
C--------------------------------------------------------------

      IMPLICIT none

C-----------------------------------------------------------------------                               
C DECLARATION DES VARIABLES
c pointeurs mex   
      CRONOSINT plhs(*), prhs(*)
      integer nlhs, nrhs
C taille des tableau
      integer NMAX
      PARAMETER (NMAX = 1000000)
c variables de sortie
      real*8 Y_OUT(NMAX)
c taille des tableaux d'entree fixee a 10**6 elements
      real*8 X_IN(NMAX),Y_IN(NMAX),X_OUT(NMAX)
c pointeurs et tailles      
      CRONOSINT YOUT_pr,XIN_pr,YIN_pr,XOUT_pr,d_pr,FAIL_pr
      integer nb_in,nb_out,m ,n
      integer mxGetM, mxGetN
      CRONOSINT mxCreateDoubleMatrix, mxGetPr

C variables pour le calcul
      real*8 d(NMAX), FAIL
      integer INCFD
      PARAMETER  (INCFD=1)
      integer    IERR
      logical    SKIP

      SKIP = .TRUE.

C      write(*,*) 'debut mexfile',nlhs,nrhs, NMAX
C      write(*,*) X_IN(100000)
C---------------------------------------------------------------------
C GESTION DES ENTREES     
C taille en entree
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      nb_in = m*n
C      write(*,*) 'nbpoints_in=',nb_in
C taille en entree
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      nb_out = m*n
C      write(*,*) 'nbpoints_out=',nb_out
C lecture des donnees en entree
      XIN_pr = mxGetPr(prhs(1))
C      write(*,*) '0',XIN_pr
      call mxCopyPtrToReal8(XIN_pr,X_IN,nb_in)
C      write(*,*) '1'
      YIN_pr = mxGetPr(prhs(2))
C      write(*,*) '2',YIN_pr
      call mxCopyPtrToReal8(YIN_pr,Y_IN,nb_in)
C      write(*,*) '3'
      XOUT_pr = mxGetPr(prhs(3))
C      write(*,*) '4',XOUT_pr
      call mxCopyPtrToReal8(XOUT_pr,X_OUT,nb_out)
C      write(*,*) 'apres lecture'

C calcul des pentes
      call DPCHIM (nb_in, X_IN, Y_IN, d, INCFD, IERR)   
      if (ierr.lt.0) THEN
 		write(*,*) 'PCHIP: ERROR in DPCHIM'
      ELSE
C evaluation  
      	 call DPCHFE(nb_in,X_IN,Y_IN,d,INCFD,SKIP,
     &               nb_out,X_OUT,Y_OUT,IERR)
      ENDIF
C      write(*,*) 'apres DPCHFE avant ierr',ierr
      if (ierr.lt.0) THEN
 		write(*,*) 'PCHIP: ERROR in DPCHFE'
      ENDIF
C      write(*,*) 'apres DPCHFE'
      FAIL = 1.0 * ierr

C---------------------------------------------------------------------
C GESTION DES SORTIES
C  YOUT
      plhs(1) = mxCreateDoubleMatrix(nb_out,1,0)
      YOUT_pr = mxGetPr(plhs(1))
      call mxCopyReal8ToPtr(Y_OUT,YOUT_pr,nb_out)
      plhs(2) = mxCreateDoubleMatrix(nb_in,1,0)
      d_pr = mxGetPr(plhs(2))
      call mxCopyReal8ToPtr(d,d_pr,nb_in)
      plhs(3) = mxCreateDoubleMatrix(1,1,0)
      FAIL_pr = mxGetPr(plhs(3))
      call mxCopyReal8ToPtr(FAIL,FAIL_pr,1)
C      write(*,*) 'fin mexfile'

c FIN DE LA FONCTION
      return
      END
