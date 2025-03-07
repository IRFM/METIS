C
C   THIS FILE IS PART OF MUMPS VERSION 4.7.3
C   This Version was built on Fri May  4 15:54:01 2007
C
C
C  This version of MUMPS is provided to you free of charge. It is public
C  domain, based on public domain software developed during the Esprit IV
C  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
C  Since this first public domain version in 1999, the developments are
C  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
C  INRIA.
C
C  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
C  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
C
C  Up-to-date copies of the MUMPS package can be obtained
C  from the Web pages http://mumps.enseeiht.fr/
C  or http://graal.ens-lyon.fr/MUMPS
C
C
C   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
C
C
C  User documentation of any code that uses this software can
C  include this complete notice. You can acknowledge (using
C  references [1], [2], and [3] the contribution of this package
C  in any scientific publication dependent upon the use of the
C  package. You shall use reasonable endeavours to notify
C  the authors of the package of this publication.
C
C   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
C   Multifrontal parallel distributed symmetric and unsymmetric solvers,
C   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
C
C   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
C   S. Pralet, Hybrid scheduling for the parallel solution of linear
C   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
C
      INTEGER FUNCTION CMUMPS_50 
     *         ( SLAVEF, K48, K10, K50,
     *         NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K10, K50, NFRONT, NCB
      INTEGER NSLAVESMIN, NASS, KMAX
      REAL Wmaster, Wtotal, Wmax
      INTEGER I,ACC,X
      REAL CMUMPS_45
      INTEGER CMUMPS_497    
      EXTERNAL  CMUMPS_45, CMUMPS_497
      KMAX = CMUMPS_497( K10, NCB )
      NASS = NFRONT - NCB
      NSLAVESMIN  = 1
      IF ( K48 .EQ.0 .OR. (K48.EQ.5 .AND.K50.EQ.0)) THEN 
         NSLAVESMIN = MAX(NCB/MAX(1,KMAX),1)
      ELSE IF (K48 .EQ. 3 .OR.(K48.EQ.5 .AND.K50.NE.0) ) THEN
        Wmax    = CMUMPS_45(KMAX,NFRONT,NASS)
        Wtotal  = CMUMPS_45(NCB,NFRONT,NASS)
        Wmaster = REAL(NASS*NASS)*REAL(NASS)/(3.0)
        IF ( Wmaster .GT.  Wmax ) THEN 
         NSLAVESMIN = MAX ( NINT ( Wtotal / Wmaster ), 1 )
        ELSE
         NSLAVESMIN = MAX ( NINT ( Wtotal / Wmax ), 1 )
        ENDIF
        IF (K48 .EQ. 5) THEN
          NSLAVESMIN = MAX ( NSLAVESMIN/2, 1 )
        END IF
      ELSE IF (K48 .EQ. 4 ) THEN
         IF ( K10 > 0 ) THEN
           WRITE(*,*) 'Internal Error 1 in CMUMPS_50'
           STOP
         ENDIF
         KMAX=ABS(K10)
         IF(K50.EQ.0)THEN
            NSLAVESMIN = MAX(INT((NCB*NCB)/KMAX),1)
         ELSE
            ACC=0
            NSLAVESMIN=0
            DO WHILE (ACC.NE.NCB)
               x=INT((-REAL(NFRONT-NCB+ACC)
     $              +SQRT(((REAL(NFRONT-NCB+ACC)*
     $              REAL(NFRONT-NCB+ACC))+REAL(4)*
     *              REAL(KMAX))))/
     $              REAL(2))
               ACC=ACC+x
               NSLAVESMIN=NSLAVESMIN+1
               IF (((NCB-ACC)*NCB).LT.KMAX)THEN
                  ACC=NCB 
                  NSLAVESMIN=NSLAVESMIN+1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      NSLAVESMIN = MIN ( NSLAVESMIN,(SLAVEF-1) )
      CMUMPS_50 = 
     *               MIN ( NSLAVESMIN, NCB )
      RETURN
      END FUNCTION CMUMPS_50 
      INTEGER FUNCTION CMUMPS_52 
     *        ( SLAVEF, K48, K10, K50,
     *          NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K10, K50,NFRONT, NCB
      INTEGER NSLAVESMAX, KMAX, KMIN
      INTEGER NSLAVESMIN
      INTEGER CMUMPS_497,CMUMPS_442,
     *        CMUMPS_50,
     *        CMUMPS_46
      EXTERNAL CMUMPS_497,CMUMPS_442,
     *        CMUMPS_50,
     *        CMUMPS_46
      IF (K48 .eq. 0 .OR. K48.eq.3.OR.K48.EQ.5) THEN
         KMAX = CMUMPS_497( K10, NCB )
         KMIN = CMUMPS_442( K10, K50, KMAX, NCB)
         NSLAVESMAX = CMUMPS_46(
     &                SLAVEF, K48, K50, KMIN, NFRONT, NCB )
      ELSE
         NSLAVESMAX = SLAVEF-1
      ENDIF
      NSLAVESMIN = CMUMPS_50(
     &     SLAVEF, K48, K10, K50, NFRONT, NCB )
      NSLAVESMAX = MAX ( NSLAVESMAX, NSLAVESMIN )
      CMUMPS_52 = 
     &               MIN ( NSLAVESMAX, NCB )
      RETURN
      END FUNCTION CMUMPS_52
      SUBROUTINE CMUMPS_503( WHAT, KEEP,KEEP8,
     *           NCB, NFR, SLAVEF, NBROWMAX, MAXSURFCB )
      IMPLICIT NONE
      INTEGER, intent(in) :: WHAT, NCB, NFR, SLAVEF
      INTEGER, intent(in) :: KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER, intent(out) :: NBROWMAX, MAXSURFCB
      INTEGER KMAX, KMIN, NSLAVES, SIZEDUMMY, TABDUMMY
      EXTERNAL CMUMPS_497, CMUMPS_442,
     *         CMUMPS_50
      INTEGER CMUMPS_497, CMUMPS_442,
     *        CMUMPS_50
      IF ( WHAT .NE. 1 .and. WHAT .NE. 2 ) THEN
        IF (WHAT .NE. 4 .and. WHAT .NE. 5 .AND.
     *       KEEP(48).NE.5 ) THEN
        WRITE(*,*) "Internal error 1 in CMUMPS_503"
        CALL CMUMPS_ABORT()
        END IF
      ENDIF
      KMAX    = CMUMPS_497( KEEP(10), NCB )
      IF (WHAT .EQ.1.OR.WHAT.EQ.2) THEN
        NSLAVES = CMUMPS_50( SLAVEF, KEEP(48),
     *            KEEP(10), KEEP(50),
     *            NFR, NCB )
      ELSE
        NSLAVES=SLAVEF
      ENDIF
      IF ( KEEP(48) == 0 .OR. (KEEP(48).EQ.5.AND.KEEP(50).EQ.0)) THEN
        NBROWMAX = NCB / NSLAVES + MOD( NCB, NSLAVES )
        IF ( WHAT == 2 .OR. WHAT == 5 ) MAXSURFCB = NBROWMAX * NCB
      ELSE IF (KEEP(48) == 3.OR.(KEEP(48).EQ.5.AND.KEEP(50).NE.0))THEN
        KMIN = CMUMPS_442( KEEP(10), KEEP(50), KMAX, NCB )  
        SIZEDUMMY        = 1
        IF(WHAT.GT.3)THEN
           CALL  CMUMPS_440(
     *          WHAT-3, NSLAVES, NFR, NCB, 
     *          KMIN, KMAX, SLAVEF,
     *          NBROWMAX, MAXSURFCB, TABDUMMY, SIZEDUMMY) 
        ELSE
           CALL  CMUMPS_440(
     *          WHAT, NSLAVES, NFR, NCB, 
     *          KMIN, KMAX, SLAVEF,
     *          NBROWMAX, MAXSURFCB, TABDUMMY, SIZEDUMMY) 
        ENDIF
      ELSE IF ( KEEP(48) == 4 ) THEN
         IF (KEEP(10) > 0) THEN
            WRITE(*,*) "Internal error 2 in CMUMPS_503"
            CALL CMUMPS_ABORT()
         END IF
         IF(KEEP(50).EQ.0)THEN
            IF ( ABS(KEEP(10)) * ( SLAVEF - 1 ) > NCB *NFR ) THEN
              NBROWMAX = (NCB + SLAVEF -2 ) / ( SLAVEF - 1 )
              IF ( WHAT == 2 ) MAXSURFCB = NBROWMAX *NCB
            ELSE
              NBROWMAX=INT (( ABS(KEEP(10)) + NFR - 1 )/ NFR)
              IF ( WHAT == 2 ) MAXSURFCB = ABS(KEEP(10))
            ENDIF
         ELSE
            NBROWMAX=INT((-(NFR-NCB)
     $              +SQRT(REAL(((NFR-NCB)*
     $              (NFR-NCB))+REAL(4)*
     *              (ABS(KEEP(10))))))/
     $              REAL(2))
            IF ( WHAT == 2 ) MAXSURFCB = ABS(KEEP(10))
         ENDIF
      ELSE
        NBROWMAX = NCB
        IF (WHAT == 2) MAXSURFCB = NCB * NCB
      ENDIF 
      NBROWMAX = MIN ( MAX(NBROWMAX, 1), NCB)
      RETURN
      END SUBROUTINE CMUMPS_503
      INTEGER FUNCTION CMUMPS_46( SLAVEF, K48, K50,
     *         BLSIZE, NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K50, BLSIZE, NFRONT, NCB
      INTEGER NSLAVES, NASS
      REAL Wtotal, Wblsize
      REAL CMUMPS_45
      EXTERNAL          CMUMPS_45
      NASS = NFRONT - NCB
      NSLAVES  = SLAVEF-1
      IF ( K48 .EQ.0 .OR. (K48.EQ.5 .AND. K50.EQ.0)) THEN 
         NSLAVES = MAX(NCB/MAX(1,BLSIZE),1)
      ELSE IF (K48.EQ.3 .OR. (K48.EQ.5 .AND. K50.NE.0))THEN
        Wblsize = CMUMPS_45(BLSIZE,NFRONT,NASS)
        Wtotal  = CMUMPS_45(NCB,NFRONT,NASS)
        NSLAVES = MAX(NINT ( Wtotal / Wblsize ), 1)
      ENDIF
      CMUMPS_46 = 
     *               MIN ( NSLAVES,(SLAVEF-1) )
      RETURN
      END FUNCTION CMUMPS_46
      SUBROUTINE  CMUMPS_440( 
     *    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWMAX, MAXSURFCB, TABPOS, SIZETABPOS) 
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: GETPOSITIONS, 
     *    NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF, SIZETABPOS 
      INTEGER, INTENT (OUT) :: NBROWMAX, MAXSURFCB
      INTEGER, INTENT (OUT) :: TABPOS(SIZETABPOS) 
      REAL Wtotal, W, COSTni
      REAL delta
      INTEGER  SumNi, NCOLim1, I, BLSIZE, NASS
      LOGICAL GETROW, GETSURF, GETPOS, GET_AVGROW, GET_AVGSURF
      REAL CMUMPS_45
      EXTERNAL          CMUMPS_45
      GETROW = (GETPOSITIONS.EQ.1)
      GETSURF= (GETPOSITIONS.EQ.2)
      GETPOS = (GETPOSITIONS.EQ.3)
      GET_AVGROW = (GETPOSITIONS.EQ.4)
      GET_AVGSURF = (GETPOSITIONS.EQ.5)
      NBROWMAX  = 0
      MAXSURFCB = 0
      IF (GETPOS) THEN 
        TABPOS (1) = 1
        TABPOS (NSLAVES+1)= NCB+1
        TABPOS (SLAVEF+2) = NSLAVES
      ENDIF
      IF (NSLAVES.EQ.1) THEN
       IF ( GETSURF ) THEN
         NBROWMAX  = NCB
         MAXSURFCB = NCB * NCB
       ELSEIF ( GETROW ) THEN
         NBROWMAX  = NCB
       ENDIF
      ELSE
        NASS    = NFRONT - NCB
        W       = CMUMPS_45(NCB,NFRONT,NASS)
        SumNi   = 0
        NCOLim1 = NASS
        DO I = 1, NSLAVES-1
          delta   = REAL(2*NCOLim1-NASS+1)**2 + 
     &                  (REAL(4)*W)/REAL(NASS*(NSLAVES-I+1))
          delta   = SQRT(delta)
          delta   = (REAL(-2*NCOLim1+NASS-1) + delta )/REAL(2)
          BLSIZE  = MAX(INT(delta), 1)
          IF ( (NFRONT-NCOLim1-BLSIZE) .LE. NSLAVES-I ) THEN
            BLSIZE = 1
          ENDIF
          NCOLim1 = NCOLim1+BLSIZE
          COSTni  = CMUMPS_45(BLSIZE,NCOLim1,NASS)
          W       = W - COSTni
          IF (GETPOS) TABPOS(I) = SumNi + 1
          IF (GETSURF) THEN
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
            MAXSURFCB = MAX ( MAXSURFCB, 
     *       BLSIZE* (SumNi+BLSIZE) )
          ELSEIF ( GETROW ) THEN         
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
             RETURN
          ELSEIF (GET_AVGSURF) THEN
            NBROWMAX = NBROWMAX + BLSIZE
            MAXSURFCB = MAXSURFCB + BLSIZE*(SumNi+BLSIZE) 
          ELSEIF (GET_AVGROW) THEN
             NBROWMAX = NBROWMAX + BLSIZE
          ENDIF
          SumNi   = SumNi + BLSIZE
        ENDDO
        BLSIZE = NCB - SumNi
        IF (BLSIZE.LE.0) THEN
          write(*,*) ' Error in CMUMPS_440: ', 
     &     ' size lastbloc ', BLSIZE
          CALL CMUMPS_ABORT()
        ENDIF
        if (NCOLim1+BLSIZE.NE.NFRONT) then
          write(*,*) ' Error in CMUMPS_440: ', 
     &     ' NCOLim1, BLSIZE, NFRONT=', 
     *       NCOLim1, BLSIZE, NFRONT
          CALL CMUMPS_ABORT()
        endif
        IF (GETPOS) TABPOS(NSLAVES) = SumNi + 1
        IF (GETSURF) THEN
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
            MAXSURFCB = MAX ( MAXSURFCB, 
     *       BLSIZE* (SumNi+BLSIZE) )
        ELSEIF ( GETROW ) THEN         
            NBROWMAX  = MAX ( NBROWMAX, 
     *       BLSIZE )
        ELSEIF (GET_AVGSURF) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          MAXSURFCB = MAXSURFCB + BLSIZE*(SumNi+BLSIZE)
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
          MAXSURFCB=(MAXSURFCB+NSLAVES-1)/NSLAVES
        ELSEIF (GET_AVGROW) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_440
      SUBROUTINE CMUMPS_441(      
     *            KEEP,KEEP8, SLAVEF,
     *            TAB_POS_IN_PERE,
     *            NSLAVES, NFRONT, NCB
     *             )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NSLAVES, SLAVEF, NFRONT,
     *                         KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER TAB_POS_IN_PERE(SLAVEF+2)
      INTEGER :: I, BLSIZE
      INTEGER KMIN, KMAX, NBROWDUMMY, MAXSURFDUMMY,
     *        GETPOSITIONS, SIZECOLTAB
      INTEGER CMUMPS_442, CMUMPS_497 
      EXTERNAL CMUMPS_442, CMUMPS_497, 
     *        CMUMPS_440
       IF (KEEP(48).EQ.0) THEN
        BLSIZE = NCB / NSLAVES
        TAB_POS_IN_PERE( 1 ) = 1
        DO I = 1, NSLAVES-1
          TAB_POS_IN_PERE( I+1 ) = TAB_POS_IN_PERE(I) +
     *    BLSIZE
        ENDDO
        TAB_POS_IN_PERE(NSLAVES+1) = NCB+1
        TAB_POS_IN_PERE(SLAVEF+2)  = NSLAVES
        RETURN
      ELSE IF (KEEP(48).EQ.3 ) THEN
        KMAX = CMUMPS_497(KEEP(10), NCB)
        KMIN = CMUMPS_442(KEEP(10), KEEP(50), KMAX, NCB)
        GETPOSITIONS = 3
        SIZECOLTAB       = SLAVEF+2
        CALL  CMUMPS_440(
     *    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     *    KMIN, KMAX, SLAVEF,
     *    NBROWDUMMY, MAXSURFDUMMY, 
     *    TAB_POS_IN_PERE(1), SIZECOLTAB)
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_441
      SUBROUTINE CMUMPS_49(
     *            KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     *            ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *            ISLAVE, NCB, NSLAVES, SIZE, FIRST_INDEX )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: ISLAVE, NCB, NSLAVES, SLAVEF,
     *                         KEEP(500), INODE, N 
      INTEGER*8 KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER, INTENT( OUT ):: SIZE, FIRST_INDEX
      INTEGER BLSIZE, J
      IF (KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       IF ( ISLAVE .NE. NSLAVES ) THEN
        SIZE = BLSIZE
       ELSE
        SIZE = BLSIZE + MOD( NCB, NSLAVES )
       END IF
       FIRST_INDEX = ( ISLAVE - 1 ) * BLSIZE + 1
      ELSEIF (KEEP(48).EQ.3) THEN
       J = ISTEP_TO_INIV2 ( STEP(INODE) )
       FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
       SIZE        = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.5) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSE
       WRITE(*,*) 'Error in CMUMPS_BLOC2 undef strat'  
       CALL CMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_49
      REAL FUNCTION CMUMPS_45(NROW,NCOL,NASS)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: NROW,NCOL,NASS
      CMUMPS_45 = REAL(NASS*NROW)*
     &                 REAL(2*NCOL - NASS - NROW + 1)
      RETURN
      END FUNCTION CMUMPS_45
      INTEGER FUNCTION CMUMPS_12 
     &      (K10, K48, K50, SLAVEF, 
     &      NCB, NFRONT, NSLAVES_less, NMB_OF_CAND )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NFRONT, NSLAVES_less, 
     *                  K10, K48, K50, SLAVEF, NMB_OF_CAND
      INTEGER J, NSLAVES, NSLAVES_REAL, ISLAVE
      INTEGER KMIN, KMAX, BLSIZE, NPIV,
     *        NSLAVES_ref, KMINSURF, NSLAVES_max
      REAL WK_MASTER, WK_SLAVE
      INTEGER  CMUMPS_497, CMUMPS_50, 
     *         CMUMPS_52
      REAL  CMUMPS_45
      EXTERNAL CMUMPS_497, CMUMPS_50, 
     *         CMUMPS_52
      EXTERNAL CMUMPS_45
      IF (NMB_OF_CAND.LE.0) THEN
      ENDIF
      IF ( (K48.EQ.0).OR. (K48.EQ.3) ) THEN
         KMAX = CMUMPS_497( K10, NCB )
         NSLAVES_ref = CMUMPS_50(
     &     SLAVEF, K48, K10, K50, NFRONT, NCB )
         NSLAVES = NSLAVES_ref
         IF ( NSLAVES_ref.LT.SLAVEF ) THEN
           NSLAVES_max = CMUMPS_52(
     *       SLAVEF, K48, K10, K50, NFRONT, NCB )
           IF ( NSLAVES_max .LT. NSLAVES_less ) THEN
            NSLAVES =  NSLAVES_max
           ELSE 
            NSLAVES =  NSLAVES_less
           ENDIF
           NSLAVES = MAX(NSLAVES_ref,NSLAVES)
         ENDIF
         NSLAVES = MIN (NSLAVES, NMB_OF_CAND)
         IF ( NSLAVES.GT.NSLAVES_ref) THEN
          NPIV = NFRONT - NCB
          IF ( K50.EQ.0 ) THEN
           WK_SLAVE = REAL( NPIV ) * REAL( NCB ) *
     *         ( 2.0E0 * NFRONT - NPIV ) / NSLAVES
           WK_MASTER = 0.6667E0 * NPIV*NPIV*NPIV+
     *                NPIV*NPIV*NCB
          ELSE
           WK_SLAVE = CMUMPS_45(NCB,NFRONT,NPIV)
     *         / NSLAVES
           WK_MASTER =  REAL(NPIV*NPIV)*REAL(NPIV)/(3.0)
          ENDIF
          IF ( (WK_MASTER.GT.WK_SLAVE).AND.
     *     (WK_SLAVE.GT.1.0) ) THEN
           NSLAVES = 
     *         INT( REAL(NSLAVES) * (WK_SLAVE/WK_MASTER))
           NSLAVES = MAX(NSLAVES_ref, NSLAVES)
          ENDIF
         ENDIF
      ELSE
       NSLAVES = NSLAVES_LESS
      ENDIF
      NSLAVES = MIN (NSLAVES, NCB)
      NSLAVES = MIN (NSLAVES, NMB_OF_CAND)
      CMUMPS_12 = NSLAVES
      RETURN
      END FUNCTION CMUMPS_12
      SUBROUTINE CMUMPS_47( 
     *   KEEP,KEEP8, INODE, STEP, N, SLAVEF, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *   NASS, NCB, 
     *   NSLAVES, POSITION, ISLAVE, IPOSSLAVE )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: KEEP(500),INODE,N,SLAVEF 
      INTEGER*8 KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER, INTENT( IN  ) :: NASS, NCB, 
     &                          NSLAVES, POSITION
      INTEGER, INTENT( OUT ) :: ISLAVE, IPOSSLAVE
      INTEGER BLSIZE, J, ISHIFT
      IF ((NSLAVES.LE.0).OR.(POSITION.LE.NASS)) THEN
       ISLAVE = 0
       IPOSSLAVE = POSITION
       RETURN
      ENDIF
      IF ( KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       ISLAVE    = MIN( NSLAVES,
     *               ( POSITION - NASS - 1 ) / BLSIZE + 1 )
       IPOSSLAVE = POSITION - NASS - ( ISLAVE - 1 ) * BLSIZE
      ELSEIF (KEEP(48).EQ.3) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSEIF (KEEP(48).EQ.5) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         ISHIFT = POSITION - NASS
         DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
           IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
           EXIT
          END IF
         END DO
      ELSE
       WRITE(*,*) 'Error in CMUMPS_47: undef strat'  
       CALL CMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_47
      INTEGER FUNCTION CMUMPS_442( K10, K50, KMAX, NCB )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: K10, KMAX, NCB, K50
      INTEGER KMIN, MINGRAN, KMINSURF, MINSURF
      IF ( ( NCB .LE.0 ).OR. (KMAX.LE.0) ) THEN 
        CMUMPS_442 = 1
        RETURN
      ENDIF
      IF (K50.EQ.0) THEN
       KMINSURF = 60000
#if defined(t3e) || defined(sgi)
       MINGRAN = 40 
#else
       MINGRAN = 50 
#endif
      ELSE
       KMINSURF = 30000
#if defined(t3e) || defined(sgi)
       MINGRAN = 10 
#else
       MINGRAN = 20 
#endif
      ENDIF
      IF (K10.GT.0) THEN
#if defined(t3e) || defined(sgi)
           KMIN = MAX(MINGRAN,KMAX/10)
#else
           KMIN = MAX(MINGRAN,KMAX/20)
#endif
      ELSE
           KMINSURF = MAX ( 
     &                INT(IABS(K10)/500), KMINSURF)
           KMIN     = MAX(INT(KMINSURF/MAX(NCB,1)), 1)
      ENDIF
      KMIN = min(KMIN,KMAX)
      KMIN = max(KMIN,1)  
      CMUMPS_442 = KMIN
      RETURN
      END FUNCTION CMUMPS_442
      INTEGER FUNCTION CMUMPS_497( KEEP10, NCB )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: KEEP10, NCB
      INTEGER KMAX 
      IF ( NCB .LE.0 ) THEN 
        CMUMPS_497 = 1
        RETURN
      ENDIF
      IF ( KEEP10.GT.0 ) THEN
       KMAX = KEEP10
      ELSE
       KMAX =  -KEEP10/NCB 
      ENDIF
      KMAX = MIN (NCB, KMAX)
      CMUMPS_497 = MAX ( KMAX, 1 )
      RETURN
      END FUNCTION CMUMPS_497
