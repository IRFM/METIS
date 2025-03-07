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
      SUBROUTINE ZMUMPS_324(A, NFRONT, NPIV, NBROW, K59)
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NBROW, K59
      COMPLEX*16 A(NFRONT*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J
      IF ((NPIV.EQ.0).OR.(NFRONT.EQ.NPIV)) GOTO 500
      IOLD = NFRONT + 1
      INEW = NPIV + 1
      IF ( NPIV < K59 ) THEN
        DO I = 1, NBROW+NPIV-1
          DO J = 0, NPIV-1
            A(INEW+J) = A(IOLD+J)
          ENDDO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      ELSE
        DO I = 1, NPIV - 2
          DO J = 0, I+1
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        END DO
        I = NPIV - 1
        DO J = 0, I
            A( INEW + J ) = A( IOLD + J )
        END DO
        INEW = INEW + NPIV
        IOLD = IOLD + NFRONT
        DO I = NPIV, NBROW+NPIV-1
          DO J = 0, NPIV - 1
            A( INEW + J ) = A( IOLD + J )
          END DO
          INEW = INEW + NPIV
          IOLD = IOLD + NFRONT
        ENDDO
      END IF
 500  RETURN
      END SUBROUTINE ZMUMPS_324
      SUBROUTINE ZMUMPS_325(A, NPIV, NBROW)
      IMPLICIT NONE
      INTEGER NPIV, NBROW
      COMPLEX*16 A((NBROW+NPIV)*(NBROW+NPIV))
      INTEGER IOLD, INEW, I , J, LDA
      IF ((NPIV.EQ.0).OR.(NBROW.EQ.0)) GOTO 500
      LDA  = NBROW+NPIV
      IOLD = LDA + 1
      INEW = NPIV + 1
      DO I = 1, NBROW+NPIV-1
        DO J = 0, NPIV-1
         A(INEW+J) = A(IOLD+J)
        ENDDO
        INEW = INEW + NPIV
        IOLD = IOLD + LDA
      ENDDO
 500  RETURN
      END SUBROUTINE ZMUMPS_325
      SUBROUTINE ZMUMPS_308( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_308
      SUBROUTINE ZMUMPS_309( N, NASS1, PERM, ITLOC,
     *                             SMALL, LSMALL,
     *                             LARGE, LLARGE,
     *                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N )
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE ZMUMPS_309
      SUBROUTINE ZMUMPS_88( MYID, N, PERM,
     *                             ITLOC,
     *                             IW, LIW, NASS1 )
      IMPLICIT NONE
      INTEGER MYID, N, LIW, NASS1
      INTEGER PERM( N ), ITLOC( N ), IW( LIW )
      INTEGER I, J
      IF ( LIW .eq. 0 ) RETURN
      IF ( NASS1 .eq. 0 ) THEN
      END IF
      DO I = 1, LIW - 1
        IF ( PERM( IW( I ) ) .GE. PERM( IW( I + 1 ) ) ) THEN
          WRITE(*,*) MYID,':Error in CHECK_ORDER_LIST: not sorted'
          WRITE(*,*) MYID,':List = ',IW(1:LIW)
          WRITE(*,*) MYID,':Permuted list =', (PERM(IW(J)),J=1,LIW)
          CALL ZMUMPS_ABORT()
        END IF
        IF (ITLOC( IW( I + 1 ) ) .NE. ITLOC( IW( I ) ) + 1 ) THEN
          WRITE(*,*) MYID,':Error in MUMP_CHECK_ORDER_LIST: bad ITLOC'
          WRITE(*,*) MYID,':ITLOC( IW ) =',(ITLOC(IW(J)),J=1,LIW)
          CALL ZMUMPS_ABORT() 
        END IF
      END DO
      IF ( ITLOC( IW( 1 ) )  .NE. NASS1 + 1 ) THEN
        WRITE(*,*) MYID,': Bad First ITLOC ITLOC(IW(1)),NASS1+1=',
     *              ITLOC( IW( 1 ) ), NASS1 + 1
        CALL ZMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_88
      SUBROUTINE ZMUMPS_181( N, NBROOT, NBLEAF, POOL, 
     &          LPOOL, PERM, FILS, FRERE, NSTK )
      IMPLICIT NONE
      INTEGER N, NBROOT, NBLEAF, LPOOL
      INTEGER POOL( LPOOL ), PERM( N ), FILS( N ),
     *        FRERE( N ), NSTK( N )
      INTEGER IPERM, INODE, IN
      IPERM = 1
      DO WHILE ( NBLEAF .NE. 0 )
        INODE = POOL( NBLEAF )
        NBLEAF = NBLEAF - 1
        IN = INODE
        DO WHILE ( IN .GT. 0 )
          PERM ( IN ) = IPERM
          IPERM = IPERM + 1
          IN = FILS( IN )
        END DO
        IN = FRERE( INODE )
        DO WHILE( IN .GT. 0 )
          IN = FRERE( IN )
        END DO
        IF ( IN .eq. 0 ) THEN
          NBROOT = NBROOT - 1
        ELSE
          IN = - IN
          NSTK( IN ) = NSTK( IN ) - 1
          IF ( NSTK( IN ) .eq. 0 ) THEN
            NBLEAF = NBLEAF + 1
            POOL( NBLEAF ) = IN
          END IF
        END IF
      END DO
      END SUBROUTINE ZMUMPS_181
      SUBROUTINE ZMUMPS_81(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF,
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG
      INTEGER KEEP(500), NFRONT_EFF
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE(KEEP(28)),
     *        NBPROCFILS(KEEP(28))
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K
      LOGICAL LEVEL1
      INCLUDE 'mumps_headers.h'
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1   ! HF includes XSIZE shift
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF 
      INEW1 = 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      IF (NUMSTK .NE. 0) THEN
        NTOTFS = NUMORG
        ISON = IFSON
        ICT11 = IOLDP2 + NFRONT
        DO 100 IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+XSIZE)
          NELIM = IW(J2 + 1+XSIZE)
          NPIVS = IW(J2 + 3+XSIZE)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .EQ. 0) GOTO 70
          DO 60 JJ = J1, J3
            NTOTFS = NTOTFS + 1
            JT1 = IW(JJ)
            IW(ICT11 + NTOTFS) = JT1
            ITLOC(JT1) = NTOTFS
            IW(JJ) = NTOTFS
            IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
   60     CONTINUE
   70     J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO 80 JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN
                NEWEL = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW(NEWEL) = J
                IW(JJ) = NFRONT_EFF
                ITLOC(J) = NFRONT_EFF
              ELSE
                IW(JJ) = ITLOC(J)
              ENDIF
   80       CONTINUE
          ELSE
            DO 90 JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
   90       CONTINUE
          ENDIF
          ISON = FRERE(STEP(ISON))
  100   CONTINUE
      ENDIF
      IBROT = INODE
      DO 120 IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO 110 JJ = J1, J2
            J = INTARR(JJ)
            IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL) = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J) = NFRONT_EFF
            ELSE
              INTARR(JJ) = ITLOC(J)
            ENDIF
  110     CONTINUE
        ENDIF
  120 CONTINUE
      IF (NFRONT.NE.NFRONT_EFF) THEN
       IF (NUMORG.EQ.NASS1) THEN
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NFRONT_EFF - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ELSE
        IP1 = IOLDPS + NFRONT + HF + NUMORG 
        IP2 = IOLDPS + HF + NFRONT_EFF + NUMORG 
        IW(IP2:IP2+NASS-1) = IW(IP1:IP1+NASS-1)
        IP1 = IOLDPS + NASS1 + HF 
        IP2 = IOLDPS + HF + NFRONT - 1 
        DO I = IP1, IP2
         IW(I + NFRONT_EFF) = IW(I)
        ENDDO
        IP1 = IOLDPS + HF 
        IP2 = IOLDPS + HF + NUMORG - 1 
        DO I = IP1, IP2
          IW(I + NFRONT_EFF) = IW(I)
        ENDDO
       ENDIF
      ELSE
       IP1 = IOLDPS + NASS1 + HF 
       IP2 = IOLDPS + HF + NFRONT - 1 
       DO I = IP1, IP2
        IW(I + NFRONT) = IW(I)
       ENDDO
       IP1 = IOLDPS + HF 
       IP2 = IOLDPS + HF + NUMORG - 1 
       DO I = IP1, IP2
         IW(I + NFRONT) = IW(I)
       ENDDO
      ENDIF
      K1 = IOLDPS + HF + NUMORG 
      K2 = K1 + NFRONT_EFF - 1 + NASS
      DO 150 K = K1, K2
        I = IW(K)
        ITLOC(I) = 0
  150 CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_81
      SUBROUTINE ZMUMPS_86(MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF, PERM, 
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, IW, LIW, 
     *           INTARR, ITLOC, FILS, FRERE_STEPS, 
     *           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG)
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG, 
     *        NFRONT_EFF
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PTRAIW(N), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP(28)),
     *        NBPROCFILS(KEEP(28)), PERM(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     *        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1
      INCLUDE 'mumps_headers.h'
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      ALLOCATE(PTTRI(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      ALLOCATE(PTLAST(NUMSTK+1), stat=allocok)
      IF (allocok .GT. 0) THEN
       IFLAG = -13
       GOTO 800
      ENDIF
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      JDEBROW = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = JDEBROW
      PTLAST(NUMSTK+1) = JDEBROW + INTARR(JDEBROW-3) - 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW (INEW+NFRONT) = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK   = IW(J2    +XSIZE)
          NELIM  = IW(J2 + 1+XSIZE)
          NPIVS  = IW(J2 + 3+XSIZE)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
          ELSE
           IF (LEVEL1) THEN
            NBPROCFILS(STEP(ISON)) = 1
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     *                               NBPROCFILS(STEP(ISON))
          ENDIF
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT) THEN
            DO JJ = J1, J2
              J = IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              IW(JJ) = ITLOC(IW(JJ))
            ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT.EQ.NASS1) GOTO 600
 199  CONTINUE
      IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     *        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT) THEN 
       IBROT = FILS( INODE )
       DO IORG = 2, NUMORG
         J1    = PTRAIW(IBROT) + 2
         IBROT = FILS( IBROT )
         J2    = J1 + INTARR(J1 - 2)
         J1    = J1 + 1
         DO JJ = J1, J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
       IW( IOLDP2+NASS1+1 : IOLDP2+NFRONT_EFF ) = 
     *  IW( ICT11+NASS1+1 : ICT11+NFRONT_EFF ) 
      ELSE
        CALL ZMUMPS_308( N, PERM, 
     *           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL ZMUMPS_309( N, NASS1, PERM, ITLOC,
     *               IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     *               IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     *               IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
      END IF
      IF ( NFRONT_EFF .NE. NEWEL1_SAVE ) THEN
        IP1 = IOLDPS + NASS1 + HF    ! XSIZE is included in HF
        IP2 = IOLDPS + HF + NFRONT - 1
        DO I = IP1, IP2
          IW(I + NFRONT) = IW(I)
        ENDDO
      END IF
      IF (NFRONT.NE.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF 
        IP2 = IOLDPS + NFRONT_EFF + HF
        IW (IP2:IP2+NFRONT_EFF) = IW (IP1:IP1+NFRONT_EFF)
      ENDIF
  600 CONTINUE
      IF ((NUMSTK .NE. 0).AND.(NFRONT.GT.NASS1)) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          LSTK = IW(J2+XSIZE)
          NELIM = IW(J2 + 1+XSIZE)
          NPIVS = IW(J2 + 3+XSIZE)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5+XSIZE)
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.IWPOSCB) THEN
           NROWS = IW(J2 + 2+XSIZE)
          ENDIF
          HS = NSLSON + 6 +XSIZE
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          DO JJ = J1, J2
              J = IW(JJ)
                IW(JJ) = ITLOC(J)
          ENDDO
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        IF (J1 .LE. J2) THEN
          DO JJ = J1, J2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
        ENDIF
      ENDDO
        K1 = IOLDPS + HF + NUMORG   ! XSIZE included in HF
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (ALLOCATED(PTTRI)) DEALLOCATE(PTTRI)
      IF (ALLOCATED(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE ZMUMPS_86
      SUBROUTINE ZMUMPS_140( N, INODE, IW, LIW, A, LA,
     *                           IOLDPS, POSELT, IFLAG,
     *                           UU, NNEG, NPVW,
     *                           KEEP,KEEP8,
     *                           MYID, SEUIL, AVOID_DELAYED, ETATASS,
     *     DKEEP,PIVNUL_LIST,LPN_LIST )
      IMPLICIT NONE
      INTEGER N, INODE, LIW, LA, IFLAG, NNEG, NPVW
      INTEGER MYID, SLAVEF, IOLDPS, POSELT
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION UU, SEUIL
      INTEGER POSFAC, NIRBDU, COMP
      COMPLEX*16 A( LA )
      INTEGER IW( LIW )
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, LAELL, NBTLKJ,IBEG_BLOCK
      INTEGER NASS, NEL1, LDA
      DOUBLE PRECISION UUTEMP
      INCLUDE 'mumps_headers.h'
      EXTERNAL ZMUMPS_330, ZMUMPS_222, ZMUMPS_234, 
     *         ZMUMPS_230, ZMUMPS_226, 
     *         ZMUMPS_237
      INTEGER  ZMUMPS_330
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER PIVSIZ,IWPOS
      INOPV = 0
      SEUIL_LOC = SEUIL
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = MAX(SEUIL,EPSILON(SEUIL))
      ELSE
        UUTEMP=UU
      ENDIF
      IBEG_BLOCK = 1
      NFRONT = IW(IOLDPS+XSIZE)
      LDA    = NFRONT
      NASS   = IABS(IW(IOLDPS+2+XSIZE))
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      IW(IOLDPS+3+XSIZE) = MIN0(NASS,NBTLKJ)
      UUTEMP = UU
 50   CONTINUE
      CALL ZMUMPS_222(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     *                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     *     DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      IF (IFLAG.LT.0) GOTO 500
      IF(KEEP(109).GT. 0) THEN
         IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
            IWPOS = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE
     *              +IW(IOLDPS+5+XSIZE)
            PIVNUL_LIST(KEEP(109)) = IW(IWPOS)
         ENDIF
      ENDIF
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
         CALL ZMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *        LDA, IOLDPS,POSELT, KEEP,KEEP8, UUTEMP, ETATASS )
         GOTO 500
      END IF
      IF (INOPV.EQ.2) THEN
         CALL ZMUMPS_234(IBEG_BLOCK,
     *            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
       CALL ZMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     *                 IOLDPS,POSELT)
       IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
       GO TO 500
      ENDIF
       CALL ZMUMPS_226(IBEG_BLOCK,
     *             NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             LDA, UUTEMP, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,PIVSIZ)
       IF(PIVSIZ .EQ. 2) THEN
          IWPOS = IOLDPS+IW(IOLDPS+1+XSIZE)+6
          IW(IWPOS+NFRONT+XSIZE) = -IW(IWPOS+NFRONT+XSIZE)
       ENDIF
       IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       NPIV   = IW(IOLDPS+1+XSIZE)
       NEL1   = NASS - NPIV
       CALL ZMUMPS_234(IBEG_BLOCK,
     *            NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            LDA, IOLDPS,POSELT,NBTLKJ,KEEP(4),UUTEMP,
     *            KEEP,KEEP8)
       IF (IFINB.EQ.-1) THEN 
         CALL ZMUMPS_237(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *         LDA, IOLDPS,POSELT, KEEP,KEEP8, UUTEMP, ETATASS)
         GOTO 500
       ENDIF
      GO TO 50
 500  CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_140
      SUBROUTINE ZMUMPS_222 
     *   (NFRONT,NASS,N,INODE,IW,LIW,
     *    A,LA, INOPV,
     *    NNEG,
     *    IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     *     DKEEP,PIVNUL_LIST,LPN_LIST)
#if defined (PROFILE_BLAS_ASS_G)
      USE ZMUMPS_LOAD
#endif
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,
     *        IOLDPS, POSELT, NNEG
      INTEGER PIVSIZ,LPIV
      COMPLEX*16 A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      include 'mpif.h'
      INTEGER POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      DOUBLE PRECISION RMAX,AMAX,TMAX,SWOP,TOL
      DOUBLE PRECISION DELTA,MAXPIV
      DOUBLE PRECISION PIVNUL,FIXA
      COMPLEX*16 PIVOT,DETPIV
      PARAMETER(TOL = 1.0D-20)
      INCLUDE 'mumps_headers.h'
      INTEGER APOS, LDA
      DOUBLE PRECISION ZERO,ONE
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,K
      INTRINSIC MAX
      DATA ZERO /0.0D0/
      DATA ONE /1.0D0/
      PIVNUL = DKEEP(1)
      FIXA = DKEEP(2)
      LDA   = NFRONT
      UULOC = UU
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      NASSW   = IABS(IW(IOLDPS+3+XSIZE))
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + LDA*NPIV + NPIV
         POSPV1 = APOS
         IF(ABS(A(APOS)).LT.SEUIL) THEN
            IF(DBLE(A(APOS)) .GE. ZERO) THEN
               A(APOS) = SEUIL
            ELSE
               A(APOS) = -SEUIL
            ENDIF
            KEEP(98) = KEEP(98)+1
         ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      DO 460 IPIV=NPIVP1,NASSW
         APOS = POSELT + LDA*(IPIV-1) + NPIV
         POSPV1 = APOS + IPIV - NPIVP1
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.ZERO) THEN 
            IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
            GO TO 420
         ENDIF
         AMAX = ZERO
         JMAX = 0
         J1 = APOS
         J2 = POSPV1 - 1
         DO JJ=J1,J2
            IF(ABS(A(JJ)) .GT. AMAX) THEN
               AMAX = ABS(A(JJ))
               JMAX = IPIV - (POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA
         J2 = NASSW - IPIV
         DO JJ=1,J2
            IF(ABS(A(J1)) .GT. AMAX) THEN
               AMAX = ABS(A(J1))
               JMAX = IPIV + JJ
            ENDIF
            J1 = J1 + LDA
         ENDDO
         RMAX = ZERO
         J2 = NFRONT - NASSW 
         DO JJ=1,J2
            RMAX = MAX(ABS(A(J1)),RMAX)
            J1 = J1 + LDA
         ENDDO
         IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(FIXA.GT.ZERO) THEN
               IF(DBLE(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA
               J2 = NASSW - IPIV
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDA
               ENDDO
               J2 = NFRONT - NASSW 
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDA
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
         IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
           IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. EPSILON(SEUIL)) THEN
               IF(DBLE(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = SEUIL
               ELSE
                  A(POSPV1) = -SEUIL
               ENDIF
               PIVOT = A(POSPV1)
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
           ENDIF
         ENDIF
         IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 460
         IF (ABS(PIVOT).GT.MAX(UULOC*MAX(RMAX,AMAX),SEUIL)) THEN
               A(POSPV1) = PIVOT
               GO TO 415
         END IF
         IF (AMAX.LE.TOL) GO TO 460
         IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  IF(POSPV1-JJ .NE. IPIV-JMAX) THEN
                     RMAX = MAX(RMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDA
               J2 = NASS - IPIV
               DO JJ=1,J2
                  IF(IPIV+JJ .NE. JMAX) THEN
                     RMAX = MAX(ABS(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDA
               ENDDO
           ENDIF
           APOSJ = POSELT + (JMAX-1)*LDA + NPIV
           POSPV2 = APOSJ + JMAX - NPIVP1
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + IPIV - NPIVP1
           ELSE
              OFFDAG = APOS + JMAX - NPIVP1
           END IF
           TMAX = ZERO
           IF(JMAX .LT. IPIV) THEN
              JJ = POSPV2
              DO K = 1, NFRONT-JMAX
                 JJ = JJ+NFRONT
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=MAX(TMAX,ABS(A(JJ)))
                 ENDIF
              ENDDO
              DO K =  APOSJ, POSPV2-1
                 TMAX = MAX(TMAX,ABS(A(K)))
              ENDDO
           ELSE
              JJ = POSPV2
              DO K = 1, NFRONT-JMAX
                 JJ = JJ+NFRONT
                 TMAX=MAX(TMAX,ABS(A(JJ)))
              ENDDO
              DO K =  APOSJ, POSPV2-1
                 IF (K.NE.OFFDAG) THEN
                    TMAX = MAX(TMAX,ABS(A(K)))
                 ENDIF
              ENDDO
           ENDIF
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           IF (SEUIL.GT.ZERO) THEN
                IF (SQRT(ABS(DETPIV)) .LE. SEUIL ) GOTO 460
           ENDIF
           MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
           IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
           IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 460
           IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     +          ABS(DETPIV)) GO TO 460
           IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     +          ABS(DETPIV)) GO TO 460
           PIVSIZ = 2
           KEEP(103) = KEEP(103)+1
 415       CONTINUE
           IF (PIVSIZ.EQ.2) THEN
              LPIV = MIN(IPIV,JMAX)
           ELSE
              LPIV = IPIV              
           ENDIF
           DO K=1,PIVSIZ
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL ZMUMPS_319( A, LA, IW, LIW,
     *             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     *             LDA, NFRONT, 1, KEEP(219), KEEP(50) )
 416          LPIV = MAX(IPIV,JMAX)
              NPIVP1 = NPIVP1+1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA+1)*NPIV+1) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_222
      SUBROUTINE ZMUMPS_226(IBEG_BLOCK,
     *     NFRONT,NASS,N,INODE,IW,LIW,
     *     A,LA,LDA, UU,
     *     IOLDPS,POSELT,IFINB,LKJIB,PIVSIZ)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDA,
     &        NPBEG, IBEG_BLOCK
      COMPLEX*16    A(LA)
      DOUBLE PRECISION   UU
      INTEGER IW(LIW)
      COMPLEX*16    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER NPIV,JROW2
      INTEGER NEL2,LPOS,NEL
      COMPLEX*16 ONE, ALPHA
      DOUBLE PRECISION  ZERO
      INTEGER PIVSIZ,NPIV_NEW,OFFDAG_OLD,LPOS1,LPOS2,JJ
      INTEGER POSPV1,POSPV2,OFFDAG,IBEG,IEND,J2,K1,K2,IROW
      COMPLEX*16 SWOP,DETPIV,MULT1,MULT2
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      NEL    = NFRONT - NPIV_NEW
      IFINB  = 0
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + NPIV*(NFRONT + 1)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDA
         IF ( UU.NE.ZERO) THEN
            CALL ZCOPY(NFRONT-NPIV_NEW, A(LPOS), LDA, A(APOS+1), 1)
         ELSE
            CALL ZCOPY(NASS-NPIV_NEW, A(LPOS), LDA, A(APOS+1), 1)
         END IF
         CALL ZMUMPS_ZSYR('U', NEL2, -VALPIV, A(LPOS), LDA, A(LPOS+1), 
     &        LDA)
         IF (UU.NE.ZERO) THEN
            CALL ZSCAL(NFRONT-NPIV_NEW, VALPIV, A(LPOS), LDA)
         ELSE
            CALL ZSCAL(NASS-NPIV_NEW, VALPIV, A(LPOS), LDA)
         ENDIF
         IF (NEL2.GT.0) THEN
            K1POS = LPOS + NEL2*LDA
            IF ( UU.NE.ZERO ) THEN
               NCB1  = NFRONT - JROW2
            ELSE
               NCB1  = NASS   - JROW2
            END IF
            CALL ZGERU(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &           A(K1POS), LDA, A(K1POS+1), LDA)
         ENDIF
      ELSE
         POSPV1 = POSELT + NPIV*(NFRONT + 1)
         POSPV2 = POSPV1+NFRONT+1
         OFFDAG_OLD = POSPV2 - 1
         OFFDAG = POSPV1+1
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA - 1
         LPOS2   = LPOS1+1
         CALL ZCOPY(NFRONT-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2), 1)
         CALL ZCOPY(NFRONT-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1), 1)
         JJ = POSPV2 + NFRONT-1
         IBEG = JJ + 2
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NFRONT
            IEND = IEND + NFRONT + 1
            JJ = JJ+NFRONT
         ENDDO
         IEND = IEND-1
         DO J2 = JROW2+1,NFRONT
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NFRONT
            IEND = IEND + NFRONT
            JJ = JJ+NFRONT
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_226
      SUBROUTINE ZMUMPS_230(NFRONT,N,INODE,IW,LIW,A,LA,
     &           IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      COMPLEX*16    A(LA)
      INTEGER IW(LIW)
      COMPLEX*16    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NEL
      INTEGER LPOS,JROW
      INTEGER LDA
      COMPLEX*16 ONE
      PARAMETER(ONE=1.0D0)
        APOS   = POSELT 
        VALPIV = ONE/A(APOS)
        A(APOS) = VALPIV
        NEL    = NFRONT - 1
        IF (NEL.EQ.0) GO TO 500
        LDA    = NFRONT
        LPOS   = APOS + LDA
        CALL ZMUMPS_ZSYR('U',NEL, -VALPIV, 
     &             A(LPOS), NFRONT, A(LPOS+1), NFRONT)
          DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
          END DO
  500   CONTINUE
        RETURN
        END SUBROUTINE ZMUMPS_230
      SUBROUTINE ZMUMPS_234(IBEG_BLOCK,
     *    NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,LKJIB,LKJIT,UUTEMP,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LA,LIW, IBEG_BLOCK
      COMPLEX*16    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER I, Block
      INTEGER BLSIZE
      DOUBLE PRECISION UUTEMP
      COMPLEX*16 ONE, ALPHA
      DOUBLE PRECISION ZERO
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0, ZERO = 0.0D0)
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+XSIZE)
      JROW2  = IABS(IW(IOLDPS+3+XSIZE))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      NEL11  = NFRONT - NPIV
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + NONEL
          IW(IOLDPS+3+XSIZE)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
        IBEG_BLOCK = NPIV + 1
      ELSEIF (JROW2.LT.NASS) THEN
          IBEG_BLOCK   = NPIV + 1
          IW(IOLDPS+3+XSIZE) = MIN0(JROW2+LKJIB,NASS)
          LKJIB  = MIN0(LKJIB,NASS-NPIV)
      ENDIF
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
#if defined(SAK_BYROW)
         DO IROW = JROW2+1, NASS, BLSIZE
           Block = MIN( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + (IROW  - 1) * LDA + NPBEG - 1
           UPOS = POSELT + (NPBEG - 1) * LDA + IROW - 1
           APOS = POSELT + (IROW  - 1) * LDA + JROW2
           CALL ZGEMM( 'N','N', IROW + Block - JROW2 - 1, Block, LKJIW,
     *                ALPHA, A( UPOS ), LDA,
     *                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ENDDO
#else
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + ( IROW - 1 ) * LDA + NPBEG - 1
           UPOS = POSELT + ( NPBEG - 1 ) * LDA + IROW - 1
           APOS = POSELT + ( IROW - 1 ) * LDA + IROW - 1
           CALL ZGEMM( 'N','N', Block, NASS - IROW + 1, LKJIW,
     *                ALPHA, A( UPOS ), LDA,
     *                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         END DO
#endif
        END IF
       LPOS = POSELT + NASS*LDA + NPBEG - 1
       UPOS = POSELT + (NPBEG-1) * LDA + JROW2
       APOS = POSELT + NASS*LDA + JROW2
       IF ( UUTEMP. NE. ZERO ) THEN
         CALL ZGEMM('N', 'N', NEL1, NFRONT-NASS, LKJIW, ALPHA, 
     *              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     *              A(APOS), LDA)
       END IF
      ENDIF
      LBP = JROW2 - NPIV
      IF ( LBP  .GT. 0     .AND.
     *     NEL1 .GE. LKJIB .AND.
     *     NEL1 .GE. LBP + LBP ) THEN
      DO I = 1, LBP
        CALL ZMUMPS_319( A, LA, IW, LIW,
     *                 IOLDPS, NPIV+I, NASS-LBP+I, POSELT, NASS,
     *                 LDA, NFRONT, 1, KEEP(219), KEEP(50))
      END DO
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_234
      SUBROUTINE ZMUMPS_321(NFRONT,N,INODE,IW,LIW,A,LA,
     *                          IOLDPS,POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      COMPLEX*16    A(LA)
      INTEGER IW(LIW)
      COMPLEX*16    VALPIV
      INTEGER APOS, POSELT
      INTEGER IOLDPS,NPIV,NEL
      INTEGER LPOS,JROW
      INCLUDE 'mumps_headers.h'
      COMPLEX*16 ONE
      DATA ONE /1.0D0/
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 500
        VALPIV = A(APOS)
        LPOS   = APOS + NFRONT
        CALL ZMUMPS_ZSYR('U', NEL, -VALPIV, A(LPOS), NFRONT, A(LPOS+1), 
     &             NFRONT)
        DO JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
        ENDDO
  500   RETURN
        END SUBROUTINE ZMUMPS_321
        SUBROUTINE ZMUMPS_319( A, LA, IW, LIW,
     *                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     *                       LDA, NFRONT, LEVEL, K219, K50 )
        IMPLICIT NONE
      INTEGER LA, LIW, IOLDPS, NPIVP1, POSELT, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL, K219, K50
      COMPLEX*16 A( LA )
      INTEGER IW( LIW )
      INCLUDE 'mumps_headers.h'
      INTEGER ISW, ISWPS1, ISWPS2, APOS, IDIAG, HF
      COMPLEX*16 SWOP
            APOS = POSELT + LDA*(IPIV-1) + NPIVP1-1
            IDIAG = APOS + IPIV - NPIVP1
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              CALL ZSWAP( NPIVP1 - 1,
     *            A( POSELT + NPIVP1 - 1 ), LDA,
     *            A( POSELT + IPIV - 1 ), LDA )
            END IF
            CALL ZSWAP( NPIVP1-1, A( POSELT+(NPIVP1-1) * LDA ), 1,
     *                          A( POSELT + (IPIV-1) * LDA ), 1 )
             CALL ZSWAP( IPIV - NPIVP1 - 1,
     *                  A( POSELT + NPIVP1 * LDA + NPIVP1 - 1),
     *                  LDA, A( APOS + 1 ), 1 )
            SWOP = A( IDIAG )
            A( IDIAG  ) = A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1 )
            A( POSELT+(NPIVP1-1)*LDA+NPIVP1-1) = SWOP
            CALL ZSWAP( NASS - IPIV, A( APOS + LDA ), LDA,
     *                  A( IDIAG + LDA ), LDA )
            IF ( LEVEL .eq. 1 ) THEN
              CALL ZSWAP( NFRONT - NASS,
     *        A( APOS  + ( NASS - IPIV + 1 ) * LDA ), LDA,
     *        A( IDIAG + ( NASS - IPIV + 1 ) * LDA ), LDA )
            END IF
            IF (K219.NE.0 .AND.K50.EQ.2) THEN
             IF ( LEVEL .eq. 2) THEN
              APOS=POSELT+LDA*LDA-1
              SWOP=A(APOS+NPIVP1)
              A(APOS+NPIVP1)=A(APOS+IPIV)
              A(APOS+IPIV)=SWOP
             ENDIF
            ENDIF
        RETURN
        END SUBROUTINE ZMUMPS_319
      SUBROUTINE ZMUMPS_237(NFRONT,NASS,N,INODE,
     *    IW,LIW,A,LA,
     *    LDA,
     *    IOLDPS,POSELT,KEEP,KEEP8,UU, ETATASS)
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,INODE,LA,LIW
      COMPLEX*16    A(LA)
      DOUBLE PRECISION   UU
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDA
      INTEGER IOLDPS, ETATASS
      INCLUDE 'mumps_headers.h'
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, I, IROWEND
      INTEGER UPOS, APOS, LPOS
      INTEGER I2, I2END, Block2
      COMPLEX*16  ONE, ALPHA, BETA
      DOUBLE PRECISION ZERO
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0,ZERO=0.0D0)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(57) ) THEN
        BLSIZE = KEEP(58)
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = KEEP(218)
      NPIV = IW( IOLDPS + 1 + XSIZE)
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( UU.eq.ZERO ) THEN
         CALL ZTRSM( 'L', 'U', 'T', 'U',
     *               NPIV, NFRONT-NPIV, ONE,
     *               A( POSELT ), LDA,
     *               A( POSELT + LDA * NPIV ), LDA )
       ENDIF
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = MIN( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + NASS*LDA + (IROW-1) * LDA
        APOS = POSELT + NASS*LDA + (IROW-1) * LDA + NASS + IROW - 1
        UPOS = POSELT + NASS
        IF (UU.NE.ZERO) THEN
          UPOS = POSELT + NASS + IROW - 1
        ENDIF
        IF (UU.EQ.ZERO) THEN
         DO I = 1, NPIV
          CALL ZCOPY( Block, A(LPOS+I-1), LDA, A(UPOS+(I-1)*LDA),1)
          CALL ZSCAL( Block, A(POSELT+(LDA+1)*(I-1)),
     *                A( LPOS + I - 1 ), LDA )
         ENDDO
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = MIN(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL ZGEMM('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     *               A(UPOS+I2-1), LDA, A(LPOS+(I2-1)*LDA),LDA,
     *               BETA,
     *               A(APOS + I2-1 + (I2-1)*LDA), LDA)
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL ZGEMM( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     *              ALPHA,  A( UPOS ), LDA,
     *              A( LPOS + LDA * Block ), LDA,
     *              BETA,
     *              A( APOS + LDA * Block ), LDA )
        ENDIF
       END DO
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_237
      SUBROUTINE ZMUMPS_320( BUF, BLOCK_SIZE,
     *                           MYROW, MYCOL, NPROW, NPCOL,
     *                           A, LOCAL_M, LOCAL_N, N, MYID, COMM )
      IMPLICIT NONE
      INTEGER BLOCK_SIZE, NPROW, NPCOL, LOCAL_M, LOCAL_N, N, COMM
      INTEGER MYROW, MYCOL, MYID
      COMPLEX*16 BUF( BLOCK_SIZE * BLOCK_SIZE )
      COMPLEX*16 A( LOCAL_M, LOCAL_N )
      INTEGER NBLOCK, IBLOCK, JBLOCK, IBLOCK_SIZE, JBLOCK_SIZE
      INTEGER ROW_SOURCE, ROW_DEST, COL_SOURCE, COL_DEST
      INTEGER IGLOB, JGLOB
      INTEGER IROW_LOC_SOURCE, JCOL_LOC_SOURCE
      INTEGER IROW_LOC_DEST, JCOL_LOC_DEST
      INTEGER PROC_SOURCE, PROC_DEST
      NBLOCK = ( N - 1 ) / BLOCK_SIZE + 1
      DO IBLOCK = 1, NBLOCK
        IF ( IBLOCK .NE. NBLOCK
     *    ) THEN
          IBLOCK_SIZE = BLOCK_SIZE
        ELSE
          IBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
        END IF
        ROW_SOURCE = MOD( IBLOCK - 1, NPROW ) 
        COL_DEST   = MOD( IBLOCK - 1, NPCOL )
        IGLOB = ( IBLOCK - 1 ) * BLOCK_SIZE + 1
        IROW_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        JCOL_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( IGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( IGLOB - 1, BLOCK_SIZE ) + 1
        DO JBLOCK = 1, IBLOCK
          IF ( JBLOCK .NE. NBLOCK
     *      ) THEN
            JBLOCK_SIZE = BLOCK_SIZE
          ELSE
            JBLOCK_SIZE = N - ( NBLOCK - 1 ) * BLOCK_SIZE
          END IF
          COL_SOURCE = MOD( JBLOCK - 1, NPCOL )
          ROW_DEST   = MOD( JBLOCK - 1, NPROW )
          PROC_SOURCE = ROW_SOURCE * NPCOL + COL_SOURCE
          PROC_DEST   = ROW_DEST   * NPCOL + COL_DEST
          IF ( PROC_SOURCE .eq. PROC_DEST ) THEN
           IF ( MYID .eq. PROC_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                  ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            IF ( IBLOCK .eq. JBLOCK ) THEN
              IF ( IBLOCK_SIZE .ne. JBLOCK_SIZE ) THEN
                WRITE(*,*) MYID,': Error in calling transdiag:unsym'
                CALL ZMUMPS_ABORT()
              END IF
              CALL ZMUMPS_327( A( IROW_LOC_SOURCE,
     *                 JCOL_LOC_SOURCE),
     *                 IBLOCK_SIZE, LOCAL_M )
            ELSE
              CALL ZMUMPS_326(
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ),
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ),
     *           IBLOCK_SIZE, JBLOCK_SIZE, LOCAL_M )
            END IF
           END IF
          ELSE IF (  MYROW .eq. ROW_SOURCE 
     *    .AND. MYCOL .eq. COL_SOURCE ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            JCOL_LOC_SOURCE = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPCOL) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL ZMUMPS_293( BUF,
     *           A( IROW_LOC_SOURCE, JCOL_LOC_SOURCE ), LOCAL_M,
     *           IBLOCK_SIZE, JBLOCK_SIZE, COMM, PROC_DEST )
          ELSE IF ( MYROW .eq. ROW_DEST 
     *    .AND.     MYCOL .eq. COL_DEST ) THEN
            JGLOB = ( JBLOCK - 1 ) * BLOCK_SIZE + 1
            IROW_LOC_DEST   = BLOCK_SIZE *
     *                    ( ( JGLOB - 1 ) / (BLOCK_SIZE*NPROW) )
     *                  + MOD( JGLOB - 1, BLOCK_SIZE ) + 1
            CALL ZMUMPS_281( BUF,
     *           A( IROW_LOC_DEST, JCOL_LOC_DEST ), LOCAL_M,
     *           JBLOCK_SIZE, IBLOCK_SIZE, COMM, PROC_SOURCE )
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_320
      SUBROUTINE ZMUMPS_293( BUF, A, LDA, M, N, COMM, DEST )
      IMPLICIT NONE
      INTEGER M, N, LDA, DEST, COMM
      COMPLEX*16 BUF(*), A(LDA,*)
      INTEGER I, IBUF, IERR
      INTEGER J
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      IBUF = 1
      DO J = 1, N
        BUF( IBUF: IBUF + M - 1 ) = A( 1 : M, J )
        DO I = 1, M
        END DO
        IBUF = IBUF + M
      END DO
      CALL MPI_SEND( BUF, M * N, MPI_DOUBLE_COMPLEX,
     *     DEST, SYMMETRIZE, COMM, IERR )
      RETURN
      END SUBROUTINE ZMUMPS_293
      SUBROUTINE ZMUMPS_281( BUF, A, LDA, M, N, COMM, SOURCE )
      IMPLICIT NONE
      INTEGER LDA, M, N, COMM, SOURCE
      COMPLEX*16 BUF(*), A( LDA, *)
      INTEGER I, IBUF, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      CALL MPI_RECV( BUF(1), M * N, MPI_DOUBLE_COMPLEX, SOURCE,
     *               SYMMETRIZE, COMM, STATUS, IERR )
      IBUF = 1
      DO I = 1, M
        CALL ZCOPY( N, BUF(IBUF), 1, A(I,1), LDA )
        IBUF = IBUF + N
      END DO
      RETURN
      END
      SUBROUTINE ZMUMPS_327( A, N, LDA )
      IMPLICIT NONE
      INTEGER N,LDA
      COMPLEX*16 A( LDA, * )
      INTEGER I, J
      DO I = 2, N
        DO J = 1, I - 1
          A( J, I ) = A( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_327
      SUBROUTINE ZMUMPS_326( A1, A2, M, N, LD )
      IMPLICIT NONE
      INTEGER M,N,LD
      COMPLEX*16 A1( LD,* ), A2( LD, * )
      INTEGER I, J
      DO J = 1, N
        DO I = 1, M
          A2( J, I ) = A1( I, J )
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_326
      RECURSIVE SUBROUTINE ZMUMPS_274( 
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE ZMUMPS_BUFFER
      USE ZMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      COMPLEX*16 A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        ITLOC(N), FILS(N), ND(KEEP(28))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      COMPLEX*16 DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      LOGICAL FLAG
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER PIVI
      INTEGER POSPV1,POSPV2,OFFDAG,J2,K1,K2,JJ,LPOS1
      COMPLEX*16 MULT1,MULT2
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR, IERR_MPI
      INTEGER LAELL, NCOL, POSBLOCFACTO, NROW
      INTEGER MEM_GAIN ! for S_NOL* bands
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11, LPOS, LPOS2, DPOS, UPOS
      INTEGER I, IPOS, KPOS, IPIV, FPERE, NSLAVES_TOT,
     *        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER LCONT,NELIM,NASS, LDA, NCOL_TO_SEND,
     *        SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      COMPLEX*16, DIMENSION(:),ALLOCATABLE :: UIP21K
      INTEGER ITYPE2
      PARAMETER (ITYPE2=2)
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      COMPLEX*16 ONE,ALPHA
      COMPLEX*16 VALPIV
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      FPERE    = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     *                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     *                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOL
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          IF (LP > 0 ) WRITE(LP,*) MYID,
     *": FAILURE IN ZMUMPS_274,
     * REAL WORKSPACE TOO SMALL"
          GOTO 700
        END IF
        CALL ZMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *       NIRBDU, LRLU, IPTRLU,
     *       IWPOS, IWPOSCB, PTRIST, PTRAST,
     *       STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     *": FAILURE IN ZMUMPS_274,
     * INTEGER WORKSPACE TOO SMALL"
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
      ENDIF
      KEEP(67) = MIN(LRLUS, KEEP(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0,LAELL,KEEP,KEEP8,LRLU)
      IF ( NPIV.GT.0 ) THEN
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IPIV ), NPIV,
     *                 MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *              A(POSBLOCFACTO), NPIV*NCOL, MPI_DOUBLE_COMPLEX,
     *              COMM, IERR )
      ENDIF
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        WRITE(*,*) MYID,
     *   ': Internal ERROR 1  in ZMUMPS_274',
     *   ' INODE =', INODE,
     *   ' MAITRE_DESC_BANDE not yet received from ', MSGSOU
        CALL ZMUMPS_ABORT()
      ENDIF
      DO WHILE ( NBPROCFILS(STEP(INODE)) .NE. 0 ) 
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL ZMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, CONTSIPERENIV2,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL ZMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC, 
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS + XSIZE)
      NASS1  = IW( IOLDPS + 1 + XSIZE)
      NROW1  = IW( IOLDPS + 2 + XSIZE)
      NPIV1  = IW( IOLDPS + 3 + XSIZE)
      NSLAV1 = IW( IOLDPS + 5 + XSIZE)
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      HS     = 6 + NSLAV1 + XSIZE
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     *                       NB_BLOC_FAC
      END IF
      IF (NPIV.EQ.0) GOTO 200
      IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
        WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     *             NPIV1,NCOL,NASS1
        CALL ZMUMPS_ABORT()
      END IF
      ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
      DO 100 I = 1, NPIV
         PIVI = ABS(IW(IPIV+I-1))
         IF (PIVI.EQ.I) GOTO 100
        ISW = IW(ICT11+I)
        IW(ICT11+I) = IW(ICT11+PIVI)
        IW(ICT11+PIVI) = ISW
        IPOS = POSELT + NPIV1 + I - 1
        KPOS = POSELT + NPIV1 + PIVI - 1
        CALL ZSWAP(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
 100  CONTINUE
      ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     *": ALLOCATION FAILURE FOR UIP21K IN ZMUMPS_274"
        IFLAG = -13
        IERROR = NPIV * NROW1
        GOTO 700
      END IF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
        ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     *            stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     *": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     * IN ZMUMPS_274"
          IFLAG = -13
          IERROR = NSLAVES_FOLLOW
          GOTO 700
        END IF
        LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     *  IW(IOLDPS+6+XTRA_SLAVES_SYM+XSIZE:
     *   IOLDPS+5+XTRA_SLAVES_SYM+XSIZE+NSLAVES_FOLLOW)
      END IF
      CALL ZTRSM( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     *             A( POSBLOCFACTO ), NCOL, A(POSELT+NPIV1), NCOL1 )
      LPOS = POSELT + NPIV1
      UPOS = 1
      DO I = 1, NROW1
        UIP21K( UPOS: UPOS + NPIV-1 ) = A( LPOS: LPOS+NPIV-1)
        LPOS = LPOS + NCOL1
        UPOS = UPOS + NPIV
      END DO
      LPOS = POSELT + NPIV1
      DPOS = POSBLOCFACTO
      I = 1
      DO
         IF(I .GT. NPIV) EXIT
         IF(IW(IPIV+I-1) .GT. 0) THEN
            CALL ZSCAL( NROW1, A(DPOS), A(LPOS), NCOL1 )
            LPOS = LPOS + 1
            DPOS = DPOS + NCOL + 1
            I = I+1
         ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS+ NCOL + 1
            OFFDAG = POSPV1+1
            LPOS1 = LPOS
            DO J2 = 1,NROW1
               K1 = JJ
               K2 = JJ+1
               MULT1 = A(POSPV1)*A(LPOS1)+A(OFFDAG)*A(LPOS1+1)
               MULT2 = A(OFFDAG)*A(LPOS1)+A(POSPV2)*A(LPOS1+1)
               A(LPOS1) = MULT1
               A(LPOS1+1) = MULT2
               LPOS1 = LPOS1 + NCOL1
            ENDDO
            LPOS = LPOS + 2
            DPOS = POSPV2 + NCOL + 1
            I = I+2
         ENDIF
      ENDDO
      LPOS2 = POSELT + NPIV1
      UPOS = POSBLOCFACTO+NPIV
      LPOS  = LPOS2 + NPIV
      CALL ZGEMM('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     *           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
      DPOS = POSELT + NCOL1 - NROW1
      IF ( NROW1 .GT. KEEP(7) ) THEN
        BLSIZE = KEEP(8)
      ELSE
        BLSIZE = NROW1
      ENDIF
      IF ( NROW1 .GT. 0 ) THEN
      DO IROW = 1, NROW1, BLSIZE
        Block = MIN( BLSIZE, NROW1 - IROW + 1 )
        DPOS  = POSELT + NCOL1 - NROW1 + ( IROW - 1 ) * ( NCOL1 + 1 )
        LPOS2 = POSELT + NPIV1         + ( IROW - 1 ) * NCOL1
        UPOS  = ( IROW - 1 ) * NPIV + 1
        DO I = 1, Block
          CALL ZGEMV( 'T', NPIV, Block-I+1, ALPHA,
     *                A( LPOS2 + (I - 1) * NCOL1 ), NCOL1,
     *                UIP21K( UPOS + NPIV * ( I - 1 ) ), 1, ONE,
     *                A(DPOS+(NCOL1+1)*(I-1)),NCOL1 )
        END DO
       IF ( NROW1-IROW+1-Block .ne. 0 )
     * CALL ZGEMM( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     *             UIP21K( UPOS ), NPIV,
     *             A( LPOS2 + Block * NCOL1 ), NCOL1, ONE,
     *             A( DPOS + Block * NCOL1 ), NCOL1 )
      ENDDO
      ENDIF
      FLOP1 = DBLE(NROW1) * DBLE(NPIV) *
     *           ( 2 * NCOL  - NPIV + NROW1 +1 )
      FLOP1 = -FLOP1
      CALL ZMUMPS_190( 1, .FALSE., FLOP1, KEEP,KEEP8 )
 200  CONTINUE
      IW(IOLDPS+XSIZE) = IW(IOLDPS+XSIZE) - NPIV
      IW(IOLDPS + 3+XSIZE) = IW(IOLDPS+3+XSIZE) + NPIV
      IF (LASTBL) IW(IOLDPS+1+XSIZE) = IW(IOLDPS + 3+XSIZE)
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
          IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL ZMUMPS_64(
     *                    INODE, NPIVSENT, FPERE,
     *                    IPOSK, JPOSK,
     *                    UIP21K, NROW1,
     *                    NSLAVES_FOLLOW,
     *                    LIST_SLAVES_FOLLOW(1),
     *                    COMM, IERR )
            IF (IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV= .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL ZMUMPS_329( COMM_LOAD, ASS_IRECV,
     *         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *         MPI_ANY_SOURCE, MPI_ANY_TAG,
     *         STATUS, 
     *         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *         IWPOS, IWPOSCB, IPTRLU,
     *         LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *         PTLUST_S, PTRFAC,
     *         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, COMM,
     *         NBPROCFILS,
     *         IPOOL, LPOOL, LEAF,
     *         NBFIN, MYID, SLAVEF,
     *         root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *         LPTRAR, NELT, FRTPTR, FRTELT, 
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
             IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
           END DO
           IF ( IERR .eq. -2 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     *": FAILURE, SEND BUFFER TOO SMALL DURING
     * ZMUMPS_274"
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     *": FAILURE, RECV BUFFER TOO SMALL DURING
     * ZMUMPS_274"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( NPIV .NE. 0 ) DEALLOCATE( UIP21K )
      IOLDPS = PTRIST(STEP(INODE))
      IF (LASTBL) THEN
         IW(IOLDPS+6+XSIZE) =  IW(IOLDPS+6+XSIZE) -
     &                            TO_UPDATE_CPT_END 
         IF ( IW(IOLDPS+6+XSIZE ) .eq. 0
     *   .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 )
     *   THEN
         DEST = ZMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL ZMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
           IFLAG = -99
           GOTO 700
         END IF
         ENDIF
      END IF
      IF (LASTBL) THEN ! (1)
      IF (IW(IOLDPS+6+XSIZE) .eq. 0 ) THEN ! (2)
        IW(IOLDPS+XXS)=S_ALL
         IF (KEEP(214).EQ.1) THEN
          CALL ZMUMPS_314( N, INODE,
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     $     )
          IOLDPS = PTRIST(STEP(INODE))
          IF (KEEP(38).NE.FPERE) THEN
            IW(IOLDPS+XXS)=S_NOLCBNOCONTIG
            IF (KEEP(216).NE.3) THEN
             MEM_GAIN=IW( IOLDPS + 2 + XSIZE )*
     *                IW( IOLDPS + 3 + XSIZE )
             LRLUS = LRLUS+MEM_GAIN
             CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *              LA-LRLUS,0,-MEM_GAIN,KEEP,KEEP8,LRLU)
            ENDIF
          ENDIF
          IF (KEEP(216).EQ.2) THEN
           IF (FPERE.NE.KEEP(38)) THEN
           CALL ZMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     *         IW( IOLDPS + 2 + XSIZE ),
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 3 + XSIZE )+
     *         IW( IOLDPS + XSIZE ), 0,
     *         IW( IOLDPS + XXS ), 0 )
           IW(IOLDPS+XXS)=S_NOLCBCONTIG
           IW(IOLDPS+XXS)=S_NOLCBCONTIG
           ENDIF
          ENDIF ! End: compact type 2 CB systematically ?
         ENDIF ! KEEP(214).EQ.1
      IF ( KEEP(38).EQ.FPERE) THEN
       LCONT  = IW(IOLDPS+XSIZE)
       NROW   = IW(IOLDPS+2+XSIZE)
       NPIV   = IW(IOLDPS+3+XSIZE)
       NASS   = IW(IOLDPS+4+XSIZE)
       NELIM  = NASS-NPIV
       NCOL_TO_SEND =  LCONT-NELIM
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NASS
       SHIFT_VAL_SON      = NASS
       LDA                = LCONT + NPIV
      IF (IW(IOLDPS+8+XSIZE).EQ.0) THEN
        IW(IOLDPS+8+XSIZE) = 1
      ELSE
      ENDIF
       CALL ZMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     *    PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM.EQ.0) THEN
         IF (KEEP(214).EQ.2) THEN
          CALL ZMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     *    )
         ENDIF
         CALL ZMUMPS_626( N, INODE, PTRIST, PTRAST, IW, LIW,
     *        A, LA, LRLU, LRLUS, IWPOSCB, NIRBDU, IPTRLU, STEP,
     *        MYID, KEEP
     $         )
       ELSE
         IOLDPS = PTRIST(STEP(INODE))
         IF (IW(IOLDPS+8+XSIZE).EQ.S_ROOT2SON_CALLED) THEN
           CALL ZMUMPS_626( N, INODE, PTRIST, PTRAST, IW, LIW,
     *        A, LA, LRLU, LRLUS, IWPOSCB, NIRBDU, IPTRLU, STEP,
     *        MYID, KEEP
     $         )
         ELSE
          IW(IOLDPS+8+XSIZE) = 0
          IF (KEEP(214).EQ.1.AND.KEEP(216).NE.3) THEN
           IW(IOLDPS+XXS)=S_NOLCBNOCONTIG38
           CALL ZMUMPS_628( IW(IOLDPS),
     *                     LIW-IOLDPS+1,
     *                     MEM_GAIN )
           LRLUS = LRLUS + MEM_GAIN
              CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                LA-LRLUS,0,-MEM_GAIN,KEEP,KEEP8,LRLU)
            IF (KEEP(216).EQ.2) THEN
              CALL ZMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     *         IW( IOLDPS + 2 + XSIZE ),
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 3 + XSIZE )+
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 4 + XSIZE ) -
     *         IW( IOLDPS + 3 + XSIZE ),!NELIM
     *         IW( IOLDPS + XXS ),0)
              IW(IOLDPS+XXS)=S_NOLCBCONTIG38
            ENDIF
          ENDIF
         ENDIF ! Test on ROOT2SON already called
       ENDIF ! NELIM = 0 or not
      ENDIF ! Test if father is KEEP(38)
      ENDIF ! Test if factorization finished
      ENDIF ! Test if factorization finished
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL ZMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_274
      SUBROUTINE ZMUMPS_141( COMM_LOAD, ASS_IRECV, 
     *           N, INODE, FPERE, IW, LIW, A, LA,
     *           UU, NOFFW,
     *           NPVW,
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     *            DKEEP,PIVNUL_LIST,LPN_LIST )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, LA, NOFFW, NPVW
      COMPLEX*16 A( LA )
      DOUBLE PRECISION UU, SEUIL
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL, SLAVEF, 
     *        POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU,
     *        COMP 
      INTEGER NB_BLOC_FAC
      INTEGER ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW( LIW )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL), ITLOC(N)
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND( KEEP(28) )
      INTEGER FRERE(KEEP(28)), FILS(N) 
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     * PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28)), STEP(N)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX*16  DBLARR(MAX(1,KEEP(13)))
      LOGICAL AVOID_DELAYED 
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBTLKJ
      INTEGER NASS, IEND, IOLDPS, POSELT, LDAFS,allocok, IBEG_BLOCK
      LOGICAL LASTBL 
      LOGICAL RESET_TO_ONE, TO_UPDATE
      INTEGER K109_ON_ENTRY
      INTEGER I,J,JJ,K,IDEB
      DOUBLE PRECISION UUTEMP
      INCLUDE 'mumps_headers.h'
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      EXTERNAL ZMUMPS_223, ZMUMPS_235,
     *         ZMUMPS_227, ZMUMPS_294,
     *         ZMUMPS_44
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      INTEGER PIVSIZ,IWPOSPIV
      COMPLEX*16 ONE
      PARAMETER(ONE = 1.0D0)
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = MAX(SEUIL,EPSILON(SEUIL))
      ELSE
        SEUIL_LOC=SEUIL
        UUTEMP=UU
      ENDIF
      RESET_TO_ONE = ((KEEP(110).GT.0).AND.(DKEEP(2).LE.0))
      IF (RESET_TO_ONE) THEN
        K109_ON_ENTRY = KEEP(109)
      ENDIF
      IBEG_BLOCK=1
      NB_BLOC_FAC = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST( STEP( INODE ))
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = IABS(IW(IOLDPS+2+XSIZE))
      LDAFS  = NASS
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5), NASS )
      ENDIF
      IW(IOLDPS+3+XSIZE) = MIN0(NASS,NBTLKJ)
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     * ' integers'
        IFLAG=-13
        IERROR=NASS
        GO TO 490
      END IF
 50   CONTINUE
      IBEGKJI = IBEG_BLOCK
      CALL ZMUMPS_223(
     *                NFRONT,NASS,IBEGKJI, NASS, IPIV,
     *                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     *                IFLAG,IOLDPS,POSELT,UU, SEUIL_LOC,
     *                KEEP,KEEP8,PIVSIZ,
     *           DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      IF (IFLAG.LT.0) GOTO 490
      IF(KEEP(109).GT. 0) THEN
         IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
            IWPOSPIV = IOLDPS+IW(IOLDPS+1+XSIZE)+6
     *              +IW(IOLDPS+5+XSIZE)
            PIVNUL_LIST(KEEP(109)) = IW(IWPOSPIV+XSIZE)
         ENDIF
      ENDIF
         IF(INOPV.EQ. 1 .AND. STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1+XSIZE)
          CALL ZMUMPS_294( COMM_LOAD, ASS_IRECV,
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS,
     *             IBEGKJI, IEND, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
      ENDIF
      IF (INOPV.EQ.1) GO TO 500
      IF (INOPV.EQ.2) THEN
         CALL ZMUMPS_235(IBEG_BLOCK,
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP,KEEP8)
         GOTO 50
      ENDIF
      NPVW = NPVW + PIVSIZ
      IF (NASS.LE.1) THEN
        IFINB = -1
        IF (NASS == 1) A(POSELT)=ONE/A(POSELT)
      ELSE
         CALL ZMUMPS_227(IBEG_BLOCK,
     *             NASS, N,INODE,IW,LIW,A,LA,
     *             LDAFS, IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4),PIVSIZ)
         IF(PIVSIZ .EQ. 2) THEN
            IWPOSPIV = IOLDPS+XSIZE+IW(IOLDPS+1+XSIZE)+6+
     &                 IW(IOLDPS+5+XSIZE)
            IW(IWPOSPIV+NFRONT) = -IW(IWPOSPIV+NFRONT)
         ENDIF
      ENDIF
      IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1+XSIZE)
          CALL ZMUMPS_294(COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDAFS, 
     *             IBEGKJI, IEND, IPIV, NASS, LASTBL,NB_BLOC_FAC,
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     *             PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1+XSIZE)
       CALL ZMUMPS_235(IBEG_BLOCK,
     *            NASS,N,INODE,IW,LIW,A,LA,
     *            LDAFS, 
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4),KEEP,KEEP8)
      GO TO 50
 490  CONTINUE
      CALL ZMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      IF (RESET_TO_ONE.AND.K109_ON_ENTRY.LT.KEEP(109)) THEN
       IDEB =  IOLDPS+XSIZE+IW(IOLDPS+5+XSIZE)+6
       JJ= IDEB
       TO_UPDATE=.FALSE. ! for safety only
       DO K = K109_ON_ENTRY+1, KEEP(109)
        I = PIVNUL_LIST(K)  ! row index
        DO J=JJ,JJ+NASS
         IF (IW(J).EQ.I) THEN
           TO_UPDATE=.TRUE. ! for safety only
           EXIT
         ENDIF
        ENDDO
        IF (TO_UPDATE) THEN
            JJ= J  ! for starting next loop from JJ 
           J =  J-IDEB+1
           A(POSELT+J-1+LDAFS*(J-1))= ONE
           TO_UPDATE=.FALSE. ! for safety only
        ELSE
!       -- should never be reached
        IF (ICNTL(1).GT.0) THEN
         write(ICNTL(1),*) ' Internal error related ', 
     &                'to null pivot row detection'
        ENDIF
         EXIT
        ENDIF
       ENDDO
      ENDIF
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE ZMUMPS_141
      SUBROUTINE ZMUMPS_223( NFRONT, NASS,
     *                   IBEGKJI, NASS2, TIPIV,
     *                   N, INODE, IW, LIW,
     *                   A, LA, NNEG, 
     *                   INOPV, IFLAG,
     *                   IOLDPS, POSELT, UU, 
     *                   SEUIL,KEEP,KEEP8,PIVSIZ,
     *                   DKEEP,PIVNUL_LIST,LPN_LIST)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV
      INTEGER NASS2, IBEGKJI, NNEG
      INTEGER TIPIV( NASS2 )
      INTEGER PIVSIZ,LPIV
      COMPLEX*16 A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW) 
      INTEGER   IOLDPS, POSELT
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      include 'mpif.h'
      INTEGER POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      DOUBLE PRECISION RMAX,AMAX,TMAX,SWOP,TOL
      DOUBLE PRECISION DELTA,MAXPIV
      COMPLEX*16 PIVOT,DETPIV
      PARAMETER(TOL = 1.0D-20)
      INCLUDE 'mumps_headers.h'
      INTEGER APOSMAX
      INTEGER APOS, LDAFS
      DOUBLE PRECISION ZERO,ONE
      DOUBLE PRECISION PIVNUL,FIXA
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,J1,JJ,J2,ILOC,K
      INTRINSIC MAX
      DATA ZERO /0.0D0/
      DATA ONE /1.0D0/
        PIVNUL = DKEEP(1)
        FIXA = DKEEP(2)
        LDAFS = NASS
        UULOC = UU
        PIVSIZ = 1
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        ILOC = NPIVP1 - IBEGKJI + 1
        TIPIV( ILOC ) = ILOC
        NASSW   = IABS(IW(IOLDPS+3+XSIZE))
        APOSMAX = POSELT+LDAFS**2-1
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + LDAFS*(NPIVP1-1) + NPIV
           POSPV1 = APOS
           IF(ABS(A(APOS)).LT.SEUIL) THEN
              IF(DBLE(A(APOS)) .GE. ZERO) THEN
                 A(APOS) = SEUIL
              ELSE
                 A(APOS) = -SEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
        DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + LDAFS*(IPIV-1) + NPIV
            POSPV1 = APOS + IPIV - NPIVP1
            PIVOT = A(POSPV1)
            IF (UULOC.EQ.ZERO) THEN 
              IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
              GO TO 420
            ENDIF
            AMAX = ZERO
            JMAX = 0
            J1 = APOS
            J2 = POSPV1 - 1
            DO JJ=J1,J2
               IF(ABS(A(JJ)) .GT. AMAX) THEN
                  AMAX = ABS(A(JJ))
                  JMAX = IPIV - (POSPV1-JJ)
               ENDIF
            ENDDO
            J1 = POSPV1 + LDAFS
            J2 = NASSW - IPIV
            DO JJ=1,J2
               IF(ABS(A(J1)) .GT. AMAX) THEN
                  AMAX = MAX(ABS(A(J1)),AMAX)
                  JMAX = IPIV + JJ
               ENDIF
               J1 = J1 + LDAFS
            ENDDO
            IF (KEEP(219).NE.0) THEN
             RMAX = DBLE(A(APOSMAX+IPIV))
            ELSE
             RMAX = ZERO
            ENDIF
            J2 = NASS - NASSW
            DO JJ=1,J2
               RMAX = MAX(ABS(A(J1)),RMAX)
               J1 = J1 + LDAFS
            ENDDO
         IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(FIXA.GT.ZERO) THEN
               IF(DBLE(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDAFS
               J2 = NASSW - IPIV
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDAFS
               ENDDO
               J2 = NASS - NASSW 
               DO JJ=1,J2
                  A(J1) = ZERO
                  J1 = J1 + LDAFS
               ENDDO
                A(POSPV1) = MAX(1.0D10*RMAX, SQRT(HUGE(RMAX))/1.0D8)
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
        IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
         IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. EPSILON(SEUIL)) THEN
               IF(DBLE(PIVOT) .GE. ZERO) THEN
                  A(POSPV1) = SEUIL
               ELSE
                  A(POSPV1) = -SEUIL
               ENDIF
               PIVOT = A(POSPV1)
               WRITE(*,*) 'WARNING matrix may be singular'
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
         ENDIF
        ENDIF
        IF (MAX(AMAX,ABS(PIVOT)) .LE. TOL) GOTO 460
        IF (ABS(PIVOT).GT.MAX(UULOC*MAX(RMAX,AMAX),SEUIL)) THEN
               GO TO 415
           END IF
            IF (AMAX.LE.TOL) GO TO 460
            IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1
               DO JJ=J1,J2
                  IF(POSPV1-JJ .NE. IPIV-JMAX) THEN
                     RMAX = MAX(RMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDAFS
               J2 = NASS - IPIV
               DO JJ=1,J2
                  IF(IPIV+JJ .NE. JMAX) THEN
                     RMAX = MAX(ABS(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDAFS
               ENDDO
            ENDIF            
            APOSJ = POSELT + (JMAX-1)*LDAFS + NPIV
            POSPV2 = APOSJ + JMAX - NPIVP1
            IF (IPIV.LT.JMAX) THEN
               OFFDAG = APOSJ + IPIV - NPIVP1
            ELSE
               OFFDAG = APOS + JMAX - NPIVP1
            END IF
            IF (KEEP(219).NE.0) THEN
             TMAX = MAX(SEUIL/UULOC,DBLE(A(APOSMAX+JMAX)))
            ELSE
             TMAX = SEUIL/UULOC
            ENDIF
            IF(JMAX .LT. IPIV) THEN
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+NASS
                  IF (JMAX+K.NE.IPIV) THEN
                     TMAX=MAX(TMAX,ABS(A(JJ)))
                  ENDIF
               ENDDO
               DO K =  APOSJ, POSPV2-1
                  TMAX = MAX(TMAX,ABS(A(K)))
               ENDDO
            ELSE
               JJ = POSPV2
               DO K = 1, NASS-JMAX
                  JJ = JJ+NASS
                  TMAX=MAX(TMAX,ABS(A(JJ)))
               ENDDO
               DO K =  APOSJ, POSPV2-1
                  IF (K.NE.OFFDAG) THEN
                     TMAX = MAX(TMAX,ABS(A(K)))
                  ENDIF
               ENDDO
            ENDIF
            DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
            IF (SEUIL.GT.ZERO) THEN
               IF (SQRT(ABS(DETPIV)) .LE. SEUIL ) GOTO 460
            ENDIF
            MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
            IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
            IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 460
            IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     +           ABS(DETPIV)) GO TO 460
            IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     +           ABS(DETPIV)) GO TO 460
           PIVSIZ = 2
           KEEP(105) = KEEP(105)+1
 415       CONTINUE
           LPIV = IPIV
           IF (PIVSIZ.EQ.2) THEN
              LPIV = MIN(IPIV,JMAX)
              TIPIV(ILOC) = -(LPIV - IBEGKJI + 1)
              TIPIV(ILOC+1) = -(MAX(IPIV,JMAX) - IBEGKJI + 1)
           ELSE
              TIPIV(ILOC) = IPIV - IBEGKJI + 1
           ENDIF
           DO K=1,PIVSIZ
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL ZMUMPS_319( A, LA, IW, LIW,
     *             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     *             LDAFS, NFRONT, 2, KEEP(219), KEEP(50))
 416          LPIV = MAX(IPIV,JMAX)
              NPIVP1 = NPIVP1+1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+LDAFS*NPIV+NPIV+1) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
  420 CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_223
      SUBROUTINE ZMUMPS_235(
     *                 IBEG_BLOCK,
     *                 NASS, N, INODE,
     *                 IW, LIW, A, LA,
     *                 LDAFS,
     *                 IOLDPS, POSELT,
     *                 LKJIB, LKJIT, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW
      COMPLEX*16    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSELT, LDAFS, IBEG_BLOCK
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1
      INTEGER HF
      INTEGER LPOS,UPOS,APOS
      INTEGER LKJIT
      INTEGER LKJIBOLD, IROW
      INTEGER J, Block
      INTEGER BLSIZE
      COMPLEX*16 ONE, ALPHA
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      LKJIBOLD = LKJIB
      NPIV   = IW(IOLDPS+1+XSIZE)
      JROW2  = IABS(IW(IOLDPS+3+XSIZE))
      NPBEG  = IBEG_BLOCK
      HF     = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
      NEL1   = NASS - JROW2
      LKJIW  = NPIV - NPBEG + 1
      IF ( LKJIW .NE. LKJIB ) THEN
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB + 2 * NONEL
          IW(IOLDPS+3+XSIZE)= MIN0(NPIV+LKJIB,NASS)
          LKJIB       = MIN0(LKJIB, NASS - NPIV)
        ELSE
          LKJIB = NASS - NPIV
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
      ELSEIF (JROW2.LT.NASS) THEN
          IW(IOLDPS+3+XSIZE) = MIN0(JROW2+LKJIB,NASS)
      ENDIF
      IBEG_BLOCK = NPIV + 1
      IF (LKJIW.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( NASS - JROW2 > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NASS - JROW2
        END IF
        IF ( NASS - JROW2 .GT. 0 ) THEN
         DO IROW = JROW2+1, NASS, BLSIZE
          Block = MIN( BLSIZE, NASS - IROW + 1 )
          LPOS = POSELT + ( IROW - 1 ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS =  POSELT + (IROW-1) * LDAFS + IROW - 1
          DO J=1, Block
            CALL ZGEMV( 'T', LKJIW, Block - J + 1, ALPHA,
     &                  A( LPOS ), LDAFS, A( UPOS ), LDAFS,
     &                  ONE, A( APOS ), LDAFS )
            LPOS = LPOS + LDAFS
            APOS = APOS + LDAFS + 1
            UPOS = UPOS + 1
          END DO
          LPOS = POSELT + ( IROW - 1 + Block ) * LDAFS + NPBEG - 1
          UPOS = POSELT + ( NPBEG - 1 ) * LDAFS + IROW - 1
          APOS = POSELT + ( IROW - 1 + Block ) * LDAFS + IROW - 1
          CALL ZGEMM( 'N','N', Block, NASS - IROW + 1 - Block, LKJIW,
     *                ALPHA, A( UPOS ), LDAFS,
     *                A( LPOS ), LDAFS, ONE, A( APOS ), LDAFS )
         END DO
        END IF
      END IF
  500 CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_235
      SUBROUTINE ZMUMPS_227
     *     ( IBEG_BLOCK, NASS, N, INODE, IW, LIW,
     *     A, LA, LDAFS, 
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,PIVSIZ)
      IMPLICIT NONE
      INTEGER NASS,N,LA,LIW,INODE,IFINB,LKJIB,LDAFS,
     &        NPBEG
      COMPLEX*16    A(LA)
      INTEGER IW(LIW)
      COMPLEX*16    VALPIV
      INTEGER APOS, IOLDPS, POSELT, K1POS, NCB1
      INTEGER LKJIT, IBEG_BLOCK
      INTEGER NPIV,JROW2
      INTEGER NEL2,LPOS
      COMPLEX*16 ONE, ALPHA
      DOUBLE PRECISION  ZERO
      INTEGER PIVSIZ,NPIV_NEW,OFFDAG_OLD,LPOS1,LPOS2,JJ
      INTEGER POSPV1,POSPV2,OFFDAG,IBEG,IEND,J2,K1,K2,IROW
      COMPLEX*16 SWOP,DETPIV,MULT1,MULT2
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      INCLUDE 'mumps_headers.h'
      NPIV   = IW(IOLDPS+1+XSIZE)
      NPIV_NEW = NPIV + PIVSIZ
      IFINB  = 0
      IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
         IW(IOLDPS+3+XSIZE) = MIN0(NASS,LKJIB)
      ENDIF
      JROW2 = IW(IOLDPS+3+XSIZE)
      NPBEG = IBEG_BLOCK
      NEL2   = JROW2 - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + NPIV*(LDAFS + 1)
         VALPIV = ONE/A(APOS)
         A(APOS) = VALPIV
         LPOS   = APOS + LDAFS
         CALL ZCOPY(NASS-NPIV_NEW, A(LPOS), LDAFS, A(APOS+1), 1)
         CALL ZMUMPS_ZSYR('U', NEL2, -VALPIV, A(LPOS), LDAFS, A(LPOS+1), 
     &        LDAFS)
         CALL ZSCAL(NASS-NPIV_NEW, VALPIV, A(LPOS), LDAFS)
         IF (NEL2.GT.0) THEN
            K1POS = LPOS + NEL2*LDAFS
            NCB1  = NASS - JROW2
            CALL ZGERU(NEL2, NCB1 , ALPHA, A(APOS+1), 1, 
     &           A(K1POS), LDAFS, A(K1POS+1), LDAFS)
         ENDIF
      ELSE
         POSPV1 = POSELT + NPIV*(LDAFS + 1)
         POSPV2 = POSPV1+LDAFS+1
         OFFDAG_OLD = POSPV2 - 1
         OFFDAG = POSPV1+1
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
         A(POSPV2) = A(POSPV1)/DETPIV
         A(POSPV1) = SWOP/DETPIV
         A(OFFDAG) = -A(OFFDAG_OLD)/DETPIV
         A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDAFS - 1
         LPOS2   = LPOS1+1
         CALL ZCOPY(NASS-NPIV_NEW, A(LPOS1), LDAFS, A(POSPV1+2), 1)
         CALL ZCOPY(NASS-NPIV_NEW, A(LPOS2), LDAFS, A(POSPV2+1), 1)
         JJ = POSPV2 + NASS-1
         IBEG = JJ + 2
         IEND = IBEG
         DO J2 = 1,NEL2
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(OFFDAG)*A(K2))
            MULT2 = - (A(OFFDAG)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NASS
            IEND = IEND + NASS + 1
            JJ = JJ+NASS
         ENDDO
         IEND = IEND-1
         DO J2 = JROW2+1,NASS
            K1 = JJ
            K2 = JJ+1
            MULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
            MULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
            K1 = POSPV1+2
            K2 = POSPV2+1
            DO IROW = IBEG,IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1
               K2 = K2 + 1
            ENDDO
            A(JJ) = -MULT1
            A(JJ+1) = -MULT2
            IBEG = IBEG + NASS
            IEND = IEND + NASS
            JJ = JJ+NASS
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_227
      RECURSIVE SUBROUTINE ZMUMPS_263(
     *   COMM_LOAD, ASS_IRECV,
     *   BUFR, LBUFR,
     *   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     *   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     *   A, LA, NIRBDU, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     *   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     *   MYID, COMM, IFLAG, IERROR, NBFIN,
     *
     *    PTLUST_S, PTRFAC, root, OPASSW, OPELIW, ITLOC, FILS,  
     *    PTRARW, PTRAIW, INTARR, DBLARR,
     *    ICNTL, KEEP,KEEP8, IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *     )
      USE ZMUMPS_BUFFER
      USE ZMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, IPTRLU, 
     &        LRLU, LRLUS, LIW, LA
      INTEGER COMP
      INTEGER NIRBDU, IFLAG, IERROR, POSFAC, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     *        PTRAST(KEEP(28)),
     *        NSTK_S(KEEP(28))
      INTEGER NBPROCFILS(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER IW( LIW )
      COMPLEX*16 A( LA )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER ITLOC( N ), FILS( N )
      INTEGER ND( KEEP(28) ), FRERE_STEPS( KEEP(28) )
      INTEGER INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      COMPLEX*16  DBLARR( MAX(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ITYPE2
      PARAMETER(ITYPE2=2)
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      INTEGER INODE, IPOSK, JPOSK, NCOLU, NPIV, POSITION, IERR
      INTEGER LAELL, POSBLOCFACTO
      INTEGER MEM_GAIN ! for S_NOL* bands
      INTEGER IOLDPS, POSELT, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NELIM1, NCOL_TO_SEND
      INTEGER LDA
      INTEGER NSLAV1, HS, DEST, NSLAVES_FOLLOW
      INTEGER FPERE
      INTEGER CPOS, LPOS
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON, SHIFT_VAL_SON
      LOGICAL DYNAMIC
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER allocok
      COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: UDYNAMIC
      COMPLEX*16 ONE,ALPHA
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      DYNAMIC = .FALSE.
      POSITION  = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, JPOSK, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     *                 MPI_INTEGER, COMM, IERR )
      IF ( NPIV .LE. 0 ) THEN
      NPIV = - NPIV
        WRITE(*,*) MYID,':error, received negative NPIV in BLFAC'
        CALL ZMUMPS_ABORT()
      END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     *                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOLU, 1,
     *                 MPI_INTEGER, COMM, IERR )
      LAELL = NPIV * NCOLU
      IF ( LRLU .LT. LAELL ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          IERROR = LAELL - LRLU
          GOTO 700
        END IF
        CALL ZMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress ass..blocfacto: LRLU,LRLUS='
     *       ,LRLU,LRLUS
             IFLAG = -9
             IERROR = LAELL -LRLU
             GOTO 700
        END IF
      END IF
      LRLU  = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP(67) = MIN(LRLUS, KEEP(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0, LAELL,KEEP,KEEP8,LRLU)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(POSBLOCFACTO), NPIV*NCOLU,
     *                 MPI_DOUBLE_COMPLEX,
     *                 COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) DYNAMIC = .TRUE.
      IF ( (PTRIST(STEP( INODE )).NE.0) .AND.
     &  (IPOSK + NPIV -1 .GT. IW(PTRIST(STEP(INODE))+3+XSIZE)) )THEN
        DYNAMIC = .TRUE.
      ENDIF
      IF (DYNAMIC)  THEN
        ALLOCATE(UDYNAMIC(LAELL), stat=allocok)
        if (allocok .GT. 0) THEN
          write(*,*) myid, ' : PB allocation U in blfac_slave '
     *     , LAELL
          IFLAG = -13 
          IERROR = LAELL
          GOTO 700
        endif
        UDYNAMIC(1:LAELL) = A(POSBLOCFACTO:POSBLOCFACTO+LAELL-1)
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *          LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      DO WHILE ( PTRIST(STEP(INODE)) .EQ. 0 )
        MSGSOU = ZMUMPS_275( STEP(INODE),
     *           PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL ZMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MAITRE_DESC_BANDE,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      ENDDO
      DO WHILE ( IPOSK + NPIV -1 .GT.
     *            IW( PTRIST(STEP( INODE )) + 3 +XSIZE) )
        MSGSOU = ZMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
        SET_IRECV = .FALSE.
        BLOCKING  = .TRUE.
        MESSAGE_RECEIVED = .FALSE.
        CALL ZMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, BLOC_FACTO_SYM,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END DO
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL ZMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IOLDPS  = PTRIST(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      LCONT1 = IW( IOLDPS + XSIZE )
      NROW1  = IW( IOLDPS + 2  + XSIZE)
      NPIV1  = IW( IOLDPS + 3  + XSIZE)
      NSLAV1 = IW( IOLDPS + 5  + XSIZE)
      HS     = 6 + NSLAV1 + XSIZE
      NCOL1  = LCONT1 + NPIV1
      CPOS = POSELT + JPOSK - 1
      LPOS = POSELT + IPOSK - 1
      IF ( NPIV .GT. 0 ) THEN
       IF (DYNAMIC) THEN
        CALL ZGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            UDYNAMIC(1), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ELSE
        CALL ZGEMM( 'T', 'N', NCOLU, NROW1, NPIV, ALPHA,
     *            A( POSBLOCFACTO ), NPIV,
     *            A( LPOS ), NCOL1, ONE,
     *            A( CPOS ), NCOL1 )
       ENDIF
       FLOP1 = DBLE(NCOLU*NPIV)*DBLE(2*NROW1)
       FLOP1 = -FLOP1
       CALL ZMUMPS_190(1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF
           IW( IOLDPS + 6 + XSIZE ) = IW( IOLDPS + 6+ XSIZE ) + 1
      IF (DYNAMIC) THEN
       DEALLOCATE(UDYNAMIC)
      ELSE
        LRLU  = LRLU + LAELL
        LRLUS = LRLUS + LAELL
        POSFAC = POSFAC - LAELL
      CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                      LA-LRLUS,0,-LAELL,KEEP,KEEP8,LRLU)
      ENDIF
      NSLAVES_FOLLOW = IW( IOLDPS + 5 +XSIZE ) - XTRA_SLAVES_SYM
      IF ( IW( IOLDPS + 6  +XSIZE) .eq. 0 .and.
     *  KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 1 ) THEN
         DEST = ZMUMPS_275( STEP(INODE), PROCNODE_STEPS, SLAVEF )
         CALL ZMUMPS_62( INODE, DEST, END_NIV2_LDLT,
     *                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_BLFAC_SLAVE.'
           IFLAG = -99
           GOTO 700
         END IF
      END IF
      IF (IW(PTRIST(STEP(INODE)) + 6+XSIZE ) .eq. 0) THEN
        IW(PTRIST(STEP(INODE))+XXS)=S_ALL
        IF (KEEP(214).EQ.1) THEN
          CALL ZMUMPS_314( N, INODE,
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     $     )
          IOLDPS = PTRIST(STEP(INODE))
          IF (KEEP(38).NE.FPERE) THEN
            IW(PTRIST(STEP(INODE))+XXS)=S_NOLCBNOCONTIG
            IF (KEEP(216).NE.3) THEN
              MEM_GAIN=IW( IOLDPS + 2 + XSIZE )*
     *                 IW( IOLDPS + 3 + XSIZE )
              LRLUS = LRLUS+MEM_GAIN
              CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                LA-LRLUS,0,-MEM_GAIN,KEEP,KEEP8,LRLU)
            ENDIF
          ENDIF
          IF (KEEP(216).EQ.2)THEN
           IF (FPERE.NE.KEEP(38)) THEN
            CALL ZMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     *         IW( IOLDPS + 2 + XSIZE ),
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 3 + XSIZE )+
     *         IW( IOLDPS + XSIZE ), 0,
     *         IW( IOLDPS + XXS ), 0 )
             IW(IOLDPS+XXS) = S_NOLCBCONTIG
           ENDIF
          ENDIF ! End: compact type 2 CB systematically ?
        ENDIF ! KEEP(214).EQ.1
       IF (KEEP(38).EQ.FPERE) THEN
       NPIV1  = IW( IOLDPS + 3 +XSIZE)
       NASS1  = IW(IOLDPS+4+XSIZE)
       NELIM1  = NASS1-NPIV1
       NCOL_TO_SEND =  LCONT1-NELIM1
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+XSIZE)  + XSIZE
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW1 + NASS1
       SHIFT_VAL_SON      = NASS1
       LDA                = LCONT1 + NPIV1
      IF (IW(IOLDPS+8+XSIZE).EQ.0) THEN
        IW(IOLDPS+8+XSIZE) = 1
      ELSE
      ENDIF
       CALL ZMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE, 
     *    PTRIST, PTRAST, 
     *    root, NROW1, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA, 
     *    ROOT_CONT_STATIC, MYID, COMM,
     *    
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     *    .FALSE., ND, FRERE_STEPS,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF ( IFLAG < 0 ) GOTO 600
       IF (NELIM1.EQ.0) THEN
         IF (KEEP(214).EQ.2) THEN
          CALL ZMUMPS_314( N, INODE,  
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE2
     $       )
         ENDIF
         CALL ZMUMPS_626( N, INODE,
     *         PTRIST, PTRAST, IW, LIW, A, LA,
     *         LRLU, LRLUS, IWPOSCB,
     *         NIRBDU, IPTRLU, STEP,
     *         MYID, KEEP
     $         )
       ELSE
         IOLDPS = PTRIST(STEP(INODE))
         IF (IW(IOLDPS+8+XSIZE).EQ.S_ROOT2SON_CALLED) THEN
           CALL ZMUMPS_626( N, INODE,
     *         PTRIST, PTRAST, IW, LIW, A, LA,
     *         LRLU, LRLUS, IWPOSCB,
     *         NIRBDU, IPTRLU, STEP,
     *         MYID, KEEP
     $         )
         ELSE
          IW(IOLDPS+8+XSIZE) = 0
          IF (KEEP(214).EQ.1.AND.KEEP(216).NE.3) THEN
            IW(IOLDPS+XXS)=S_NOLCBNOCONTIG38
            CALL ZMUMPS_628( IW(IOLDPS),
     *                     LIW-IOLDPS+1,
     *                     MEM_GAIN )
            LRLUS = LRLUS + MEM_GAIN
            CALL ZMUMPS_471(.FALSE.,.FALSE.,
     *                LA-LRLUS,0,-MEM_GAIN,KEEP,KEEP8,LRLU)
            IF (KEEP(216).EQ.2) THEN
              CALL ZMUMPS_627(A,LA,PTRAST(STEP(INODE)),
     *         IW( IOLDPS + 2 + XSIZE ),
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 3 + XSIZE )+
     *         IW( IOLDPS + XSIZE ),
     *         IW( IOLDPS + 4 + XSIZE ) -
     *         IW( IOLDPS + 3 + XSIZE ),!NELIM
     *         IW( IOLDPS + XXS ),0)
              IW(IOLDPS + XXS)=S_NOLCBCONTIG38
            ENDIF
          ENDIF
         ENDIF ! Test on ROOT2SON already called
       ENDIF ! NELIM1 = 0 or not
       ENDIF ! Father = root node
       ENDIF ! Test if factorization finished
 600  CONTINUE
      RETURN
 700  CONTINUE
      CALL ZMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_263
      SUBROUTINE ZMUMPS_137( INODE, N, PROCNODE_STEPS,
     *           SLAVEF,
     &           ND, FILS, FRERE_STEPS, STEP, PIMASTER,
     *           KEEP28, KEEP50,
     *           FLOP1,
     &           IW, LIW )
      IMPLICIT NONE
      INTEGER INODE, N, KEEP50, LIW, SLAVEF, KEEP28
      INTEGER PROCNODE_STEPS(KEEP28), ND(KEEP28),
     &        FILS(N), FRERE_STEPS(KEEP28),
     *        STEP(N), 
     * PIMASTER(KEEP28),
     *  IW( LIW )
      DOUBLE PRECISION FLOP1
      INTEGER NUMORG, IN, NASS, IFSON, NUMSTK, NFRONT, NPIV, NCB, 
     &        LEVEL, ISON
      LOGICAL ZMUMPS_170
      INTEGER ZMUMPS_330
      EXTERNAL ZMUMPS_170, ZMUMPS_330
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INCLUDE 'mumps_headers.h'
      FLOP1 = 0.0D0
      IF (ZMUMPS_170(STEP(INODE),
     *    PROCNODE_STEPS, SLAVEF) ) RETURN
      IN     = INODE
      NUMORG = 0
   10 NUMORG = NUMORG + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 10
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .EQ. 0) GOTO 30
   20 NUMSTK = NUMSTK + 1
      NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 +XSIZE)
      ISON = FRERE_STEPS(STEP(ISON))
      IF (ISON .GT. 0) GOTO 20
   30 NFRONT = ND(STEP(INODE)) + NASS
      NPIV  = NASS + NUMORG
      NCB   = NFRONT - NPIV
      LEVEL = ZMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      CALL ZMUMPS_511(NFRONT,NPIV,NPIV,KEEP50,LEVEL,FLOP1)
      RETURN
      END SUBROUTINE ZMUMPS_137
      SUBROUTINE ZMUMPS_511(NFRONT,NPIV,NASS,
     *                                 KEEP50,LEVEL,COST)
      IMPLICIT NONE
      INTEGER NFRONT,NPIV,KEEP50,LEVEL, NASS
      DOUBLE PRECISION COST
      IF (KEEP50.EQ.0) THEN
        IF (LEVEL.EQ.1) THEN
          COST = DBLE(2) * DBLE(NFRONT) * DBLE(NPIV) *
     *      DBLE(NFRONT - NPIV - 1) +
     *      DBLE(NPIV) * DBLE(NPIV + 1) * DBLE(2 * NPIV + 1)
     *          / DBLE(3)
          COST = COST + DBLE(2 * NFRONT - NPIV - 1)
     *      * DBLE(NPIV) /DBLE(2)
        ELSEIF (LEVEL.EQ.2) THEN
          COST = DBLE(2*NASS)*DBLE(NFRONT) -
     *          DBLE(NASS+NFRONT)*DBLE(NPIV+1)
          COST = DBLE(NPIV)*COST +
     *     DBLE(2 * NASS - NPIV - 1) * DBLE(NPIV) / DBLE(2) +
     *     DBLE(NPIV) * DBLE(NPIV + 1) *
     *     DBLE(2 * NPIV + 1) /DBLE(3)
        ENDIF
      ELSE
        IF (LEVEL.EQ.1) THEN
          COST = DBLE(NPIV) * (
     *          DBLE( NFRONT ) * DBLE( NFRONT ) +
     *          DBLE( NFRONT ) - (
     *          DBLE( NFRONT)*DBLE(NPIV) + DBLE(NPIV+1)
     *          )) +( DBLE(NPIV)*DBLE(NPIV+1)
     *          *DBLE(2*NPIV+1))/ DBLE(6)
        ELSE IF (LEVEL.EQ.3.AND.KEEP50.EQ.2) THEN
          COST = DBLE(2) * DBLE(NFRONT) * DBLE(NPIV) *
     *      DBLE(NFRONT - NPIV - 1) +
     *      DBLE(NPIV) * DBLE(NPIV + 1) *
     *      DBLE(2 * NPIV + 1) / DBLE(3)
          COST = COST + DBLE(2 * NFRONT - NPIV - 1)
     *         * DBLE(NPIV) / DBLE(2)
        ELSE
          COST = DBLE(NPIV) * (
     *          DBLE( NASS ) * DBLE( NASS ) + DBLE( NASS )
     *        - DBLE( NASS*NPIV +  NPIV + 1 )    )
     *        + ( DBLE(NPIV)*DBLE(NPIV+1)*DBLE(2*NPIV+1) )
     *        / DBLE( 6 )
        ENDIF
      ENDIF
      END SUBROUTINE ZMUMPS_511
      SUBROUTINE ZMUMPS_38( NROW_SON, NCOL_SON, INDROW_SON,
     *                         INDCOL_SON, VAL_SON, VAL_ROOT,
     *                         LOCAL_M, LOCAL_N )
      IMPLICIT NONE
      INTEGER NCOL_SON, NROW_SON
      INTEGER INDROW_SON( NROW_SON ), INDCOL_SON( NCOL_SON )
      INTEGER LOCAL_M, LOCAL_N
      COMPLEX*16 VAL_SON( NCOL_SON, NROW_SON )
      COMPLEX*16 VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER I, J
      DO I = 1, NROW_SON
        DO J = 1, NCOL_SON
          VAL_ROOT( INDROW_SON( I ), INDCOL_SON( J ) ) =
     *    VAL_ROOT( INDROW_SON( I ), INDCOL_SON( J ) ) + VAL_SON(J,I)
        END DO
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_38
      RECURSIVE SUBROUTINE ZMUMPS_80
     *  ( COMM_LOAD, ASS_IRECV, N, ISON, IROOT,
     *    PTRI, PTRR,
     *    root,
     *    NBROW, NBCOL, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON,
     *    SHIFT_VAL_SON, LDA, TAG,
     *    MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE ZMUMPS_OOC        
      USE ZMUMPS_BUFFER
      USE ZMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, ISON, IROOT, TAG
      INTEGER PTRI( KEEP(28) ), PTRR( KEEP(28) )
      INTEGER NBROW, NBCOL, LDA
      INTEGER SHIFT_VAL_SON
      INTEGER SHIFT_LIST_ROW_SON, SHIFT_LIST_COL_SON
      INTEGER MYID, COMM
      LOGICAL INVERT
      INCLUDE 'mpif.h'
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS
      INTEGER LIW, LA
      INTEGER IW( LIW )
      COMPLEX*16 A( LA )
      INTEGER NIRBDU
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTRIST( KEEP(28) ), PTRAST(KEEP(28)), PTRFAC(KEEP(28)),
     *        PTLUST_S(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), NSTK( N )
      INTEGER COMP, IFLAG, IERROR
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NBFIN, SLAVEF
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER ITLOC( N ), FILS( N ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER INTARR( MAX(1,KEEP(14)) )
      COMPLEX*16 DBLARR( MAX(1,KEEP(13)) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTRROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTRCOL
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ROW_INDEX_LIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: COL_INDEX_LIST
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER I, POS_IN_ROOT, IROW, JCOL
      INTEGER PDEST, IERR, IERR_MPI
      INTEGER LOCAL_M, LOCAL_N, POSROOT
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER NRLOCAL, NCLOCAL
      LOGICAL SET_IRECV, BLOCKING, MESSAGE_RECEIVED
      INTEGER SIZE_MSG
      INTEGER LP
      INCLUDE 'mumps_headers.h'
      LOGICAL FLAG
      LP = ICNTL(1)
      IF ( ICNTL(4) .LE. 0 ) LP = -1
      ALLOCATE(PTRROW(root%NPROW + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPROW + 1
      endif
      ALLOCATE(PTRCOL(root%NPCOL + 1 ),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = root%NPCOL + 1
      endif
      ALLOCATE(ROW_INDEX_LIST(NBROW+1),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = NBROW + 1
      endif
      ALLOCATE(COL_INDEX_LIST(NBCOL+1),  stat=allocok)
      if (allocok .GT. 0) THEN
       IFLAG  =-13
       IERROR = NBCOL + 1
      endif
      IF (IFLAG.LT.0) THEN
         IF (LP > 0) write(6,*) myid, ' : MEMORY ALLOCATION ',
     &     'FAILURE in ZMUMPS_80'
         CALL ZMUMPS_44( MYID, SLAVEF, COMM )
         RETURN
      ENDIF
      PTRROW = 0
      PTRCOL = 0
      DO I = 1, NBROW
        POS_IN_ROOT = root%RG2L_ROW( IW( PTRI(STEP(ISON)) +
     *                          SHIFT_LIST_ROW_SON + I - 1 ) )
        IF ( .NOT. INVERT ) THEN
          IROW  = MOD( ( POS_IN_ROOT - 1 ) / root%MBLOCK, root%NPROW )
          PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
        ELSE
          JCOL =  MOD( ( POS_IN_ROOT - 1 ) / root%NBLOCK, root%NPCOL ) 
          PTRCOL( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
        END IF
      END DO
      DO I = 1, NBCOL 
        POS_IN_ROOT = root%RG2L_COL( IW( PTRI(STEP(ISON)) +
     *                SHIFT_LIST_COL_SON + I - 1 ) )
        IF ( .NOT. INVERT ) THEN
          JCOL        = MOD( ( POS_IN_ROOT - 1 ) /
     *                  root%NBLOCK, root%NPCOL )
          PTRCOL ( JCOL + 2 ) = PTRCOL( JCOL + 2 ) + 1
        ELSE
          IROW        = MOD( ( POS_IN_ROOT - 1 ) /
     *                  root%MBLOCK, root%NPROW )
          PTRROW ( IROW + 2 ) = PTRROW( IROW + 2 ) + 1
        END IF
      END DO
      PTRROW( 1 ) = 1
      DO IROW = 2, root%NPROW + 1
        PTRROW( IROW ) = PTRROW( IROW ) + PTRROW( IROW - 1 )
      END DO
      PTRCOL( 1 ) = 1
      DO JCOL = 2, root%NPCOL + 1
        PTRCOL( JCOL ) = PTRCOL( JCOL ) + PTRCOL( JCOL - 1 )
      END DO
      DO I = 1, NBROW
        POS_IN_ROOT = root%RG2L_ROW( IW( PTRI(STEP(ISON)) +
     *                          SHIFT_LIST_ROW_SON + I - 1 ) )
        IF ( .NOT. INVERT ) THEN
          IROW        = MOD( ( POS_IN_ROOT - 1 ) / root%MBLOCK,
     *                       root%NPROW )
          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ) = I
          PTRROW ( IROW + 1 ) = PTRROW( IROW + 1 ) + 1
        ELSE
          JCOL        = MOD( ( POS_IN_ROOT - 1 ) / root%NBLOCK,
     *                       root%NPCOL )
          ROW_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        END IF
      END DO
      DO I = 1, NBCOL 
        POS_IN_ROOT = root%RG2L_COL( IW( PTRI(STEP(ISON)) +
     *                SHIFT_LIST_COL_SON + I - 1 ) )
        IF ( .NOT. INVERT ) THEN
          JCOL        = MOD( ( POS_IN_ROOT - 1 ) /
     *                root%NBLOCK, root%NPCOL )
          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ) = I 
          PTRCOL ( JCOL + 1 ) = PTRCOL( JCOL + 1 ) + 1
        ELSE
          IROW        = MOD( ( POS_IN_ROOT - 1 ) /
     *                root%MBLOCK, root%NPROW )
          COL_INDEX_LIST( PTRROW( IROW + 1 ) ) = I
          PTRROW( IROW + 1 ) = PTRROW( IROW + 1 ) + 1
        END IF
      END DO
      DO IROW = root%NPROW, 2, -1
        PTRROW( IROW ) = PTRROW( IROW - 1 )
      END DO
      PTRROW( 1 ) = 1
      DO JCOL = root%NPCOL, 2, -1
        PTRCOL( JCOL ) = PTRCOL( JCOL - 1 )
      END DO
      PTRCOL( 1 ) = 1
      JCOL  = root%MYCOL
      IROW  = root%MYROW
      IF ( root%yes ) THEN
        if (irow .ne. root%myrow .or. jcol.ne.root%mycol) then
        write(*,*) ' error in grid position buildandsendcbroot'
        stop
        end if
        IF ( PTRIST(STEP(IROOT)).EQ.0.AND.
     *       PTLUST_S(STEP(IROOT)).EQ.0) THEN
           NBPROCFILS( STEP(IROOT) ) = -1
           IF (KEEP(60) .EQ. 0) THEN
             CALL ZMUMPS_284(root, IROOT, N, IW, LIW,
     *                     A, LA,
     *                     NIRBDU,
     *                     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *                     LRLU, IPTRLU,
     *                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     *                     STEP, PIMASTER, PAMASTER, ITLOC,
     *                     COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
              IF (IFLAG.LT.0) THEN
                CALL ZMUMPS_44( MYID, SLAVEF, COMM )
                RETURN
              ENDIF
           ELSE
              PTRIST(STEP(IROOT))=-66666
           ENDIF
        ELSE
           NBPROCFILS(STEP(IROOT)) = NBPROCFILS(STEP(IROOT)) - 1
           IF ( NBPROCFILS( STEP(IROOT) ) .eq. 0 ) THEN
              IF (KEEP(201).NE.0) THEN
                CALL ZMUMPS_580(IERR)              
              ENDIF
              CALL ZMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     *        SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     *        STEP, IROOT+N )
              IF (KEEP(47) .GE. 3) THEN
                 CALL ZMUMPS_500(
     $                IPOOL, LPOOL, 
     *                PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *                MYID, STEP, N, ND, FILS )
              ENDIF
          END IF
        END IF
       IF (KEEP(60) .NE. 0 ) THEN
         LOCAL_M = root%SCHUR_LLD
         LOCAL_N = root%SCHUR_NLOC
          IF ( .NOT. INVERT ) THEN
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            CALL ZMUMPS_285( N,
     *        root%SCHUR_POINTER(1),
     *        LOCAL_M, LOCAL_N,
     *        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     *        NBCOL, NBROW,
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *        LDA, A( PTRR( STEP(ISON)) + SHIFT_VAL_SON ),
     *        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *        NRLOCAL,
     *        NCLOCAL,
     *        root%RG2L_ROW, root%RG2L_COL, INVERT )
          ELSE
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            CALL ZMUMPS_285( N,
     *        root%SCHUR_POINTER(1),
     *        LOCAL_M, LOCAL_N,
     *        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     *        NBCOL, NBROW,
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *        LDA, A(PTRR(STEP(ISON)) + SHIFT_VAL_SON ),
     *        ROW_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *        COL_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *        NCLOCAL,
     *        NRLOCAL,
     *        root%RG2L_ROW, root%RG2L_COL, INVERT )
          END IF
       ELSE
        IF ( PTRIST(STEP( IROOT )) .GE. 0 ) THEN
          IF ( PTRIST(STEP( IROOT )) .EQ. 0 ) THEN
            LOCAL_N = IW( PTLUST_S(STEP(IROOT)) + 1 + XSIZE)
            LOCAL_M = IW( PTLUST_S(STEP(IROOT)) + 2 + XSIZE)
            POSROOT = PTRFAC(IW( PTLUST_S(STEP(IROOT)) +4+XSIZE ))
          ELSE
            LOCAL_N = - IW( PTRIST(STEP(IROOT)) +XSIZE)
            LOCAL_M = IW( PTRIST(STEP(IROOT)) + 1 +XSIZE)
            POSROOT = PAMASTER(STEP( IROOT ))
          ENDIF
          IF ( .NOT. INVERT ) THEN
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            CALL ZMUMPS_285( N, A( POSROOT ),
     *        LOCAL_M, LOCAL_N,
     *        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     *        NBCOL, NBROW,
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *        LDA, A( PTRR( STEP(ISON)) + SHIFT_VAL_SON ),
     *        ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *        COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *        NRLOCAL,
     *        NCLOCAL,
     *        root%RG2L_ROW, root%RG2L_COL, INVERT )
          ELSE
            NCLOCAL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
            NRLOCAL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
            CALL ZMUMPS_285( N, A( POSROOT ),
     *        LOCAL_M, LOCAL_N,
     *        root%NPCOL, root%NPROW, root%MBLOCK, root%NBLOCK,
     *        NBCOL, NBROW,
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *        IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *        LDA, A(PTRR(STEP(ISON)) + SHIFT_VAL_SON ),
     *        ROW_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *        COL_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *        NCLOCAL,
     *        NRLOCAL,
     *        root%RG2L_ROW, root%RG2L_COL, INVERT )
          END IF
        END IF
       ENDIF
      END IF
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
          PDEST = IROW * root%NPCOL + JCOL
          IF ( (root%MYROW.eq.IROW.and.root%MYCOL.eq.JCOL) .and.
     *      MYID.ne.PDEST) THEN
            write(*,*) 'error: myrow,mycol=',root%MYROW,root%MYCOL
            write(*,*) ' MYID,PDEST=',MYID,pDEST
            stop
          END IF
          IF ( root%MYROW .NE. IROW .OR. root%MYCOL .NE. JCOL) THEN
            IERR = -1
            DO WHILE ( IERR .EQ. -1 )
              IF ( .NOT. INVERT ) THEN
              NSUBSET_ROW = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
              NSUBSET_COL = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
              ELSE
              NSUBSET_ROW = PTRCOL( JCOL + 2 ) - PTRCOL( JCOL + 1 )
              NSUBSET_COL = PTRROW( IROW + 2 ) - PTRROW( IROW + 1 )
              END IF
              IF ( LRLU .LT. NSUBSET_ROW * NSUBSET_COL
     *         .and. LRLUS .GT. NSUBSET_ROW * NSUBSET_COL ) THEN
                CALL ZMUMPS_94(N, KEEP(28),
     *          IW, LIW, A, LA, NIRBDU,
     *          LRLU, IPTRLU,
     *          IWPOS, IWPOSCB, PTRIST, PTRAST,
     *          STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
                COMP = COMP + 1
                IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) MYID,': Error in b&scbroot: pb compress'
                  WRITE(*,*) MYID,': LRLU, LRLUS=',LRLU,LRLUS
                  CALL ZMUMPS_ABORT()
                END IF
              END IF
              IF ( .NOT. INVERT ) THEN
                CALL ZMUMPS_75( N, ISON,
     *          NBCOL, NBROW,
     *          IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *          IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *          LDA, A( PTRR(STEP(ISON)) + SHIFT_VAL_SON ),
     *          TAG,
     *          ROW_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *          COL_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *          NSUBSET_ROW, NSUBSET_COL,
     *          root%NPROW, root%NPCOL, root%MBLOCK,
     *          root%RG2L_ROW, root%RG2L_COL,
     *          root%NBLOCK, PDEST,
     *          COMM, IERR, A( POSFAC ), LRLU, INVERT,
     *          SIZE_MSG )
              ELSE
                CALL ZMUMPS_75( N, ISON,
     *          NBCOL, NBROW,
     *          IW( PTRI(STEP(ISON)) + SHIFT_LIST_COL_SON ),
     *          IW( PTRI(STEP(ISON)) + SHIFT_LIST_ROW_SON ),
     *          LDA, A( PTRR(STEP(ISON)) + SHIFT_VAL_SON ),
     *          TAG,
     *          ROW_INDEX_LIST( PTRCOL( JCOL + 1 ) ),
     *          COL_INDEX_LIST( PTRROW( IROW + 1 ) ),
     *          NSUBSET_ROW, NSUBSET_COL,
     *          root%NPROW, root%NPCOL, root%MBLOCK,
     *          root%RG2L_ROW, root%RG2L_COL,
     *          root%NBLOCK, PDEST,
     *          COMM, IERR, A( POSFAC ), LRLU, INVERT,
     *          SIZE_MSG )
              END IF
              IF ( IERR .EQ. -1 ) THEN
                  BLOCKING  = .FALSE.
                  SET_IRECV = .TRUE.
                  MESSAGE_RECEIVED = .FALSE.
                  CALL ZMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *            BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *            MPI_ANY_SOURCE, MPI_ANY_TAG, 
     *            STATUS, BUFR, LBUFR,
     *            LBUFR_BYTES, PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB,
     *            IPTRLU, LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *            PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     *            PIMASTER, PAMASTER, NSTK,
     *            COMP, IFLAG, IERROR, COMM, NBPROCFILS, IPOOL, LPOOL,
     *            LEAF, NBFIN, MYID, SLAVEF, root,
     *            OPASSW, OPELIW, ITLOC, FILS,
     *            PTRARW, PTRAIW, INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     *            ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     *            ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
                  IF ( IFLAG .LT. 0 ) GOTO 500
              END IF
            END DO
            IF ( IERR == -2 ) THEN
              IFLAG  = -17
              IERROR = SIZE_MSG
              IF (LP > 0) WRITE(LP, *) "FAILURE, SEND BUFFER TOO
     * SMALL DURING ZMUMPS_80"
              CALL ZMUMPS_44( MYID, SLAVEF, COMM )
              GOTO 500
            ENDIF
            IF ( IERR == -3 ) THEN
              IF (LP > 0) WRITE(LP, *) "FAILURE, RECV BUFFER TOO
     * SMALL DURING ZMUMPS_80"
              IFLAG  = -20
              IERROR = SIZE_MSG
              CALL ZMUMPS_44( MYID, SLAVEF, COMM )
              GOTO 500
            ENDIF
          END IF
        END DO
      END DO
 500  CONTINUE
      DEALLOCATE(PTRROW)
      DEALLOCATE(PTRCOL)
      DEALLOCATE(ROW_INDEX_LIST)
      DEALLOCATE(COL_INDEX_LIST)
      RETURN
      END SUBROUTINE ZMUMPS_80
      SUBROUTINE ZMUMPS_285( N, VAL_ROOT,
     *   LOCAL_M, LOCAL_N,
     *   NPCOL, NPROW, MBLOCK, NBLOCK, NBCOL_SON, NBROW_SON, INDCOL_SON,
     *   INDROW_SON, LD_SON, VAL_SON, SUBSET_ROW, SUBSET_COL,
     *   NSUBSET_ROW, NSUBSET_COL, RG2L_ROW, RG2L_COL, INVERT )
      IMPLICIT NONE
      INTEGER N, LOCAL_M, LOCAL_N
      COMPLEX*16 VAL_ROOT( LOCAL_M, LOCAL_N )
      INTEGER NPCOL, NPROW, MBLOCK, NBLOCK
      INTEGER NBCOL_SON, NBROW_SON
      INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
      INTEGER LD_SON
      COMPLEX*16 VAL_SON( LD_SON, NBROW_SON )
      INTEGER NSUBSET_ROW, NSUBSET_COL
      INTEGER SUBSET_ROW( NSUBSET_ROW ), SUBSET_COL( NSUBSET_COL )
      INTEGER RG2L_ROW( N ), RG2L_COL( N )
      LOGICAL INVERT
      INTEGER ISUB, JSUB, I, J, IPOS_ROOT, JPOS_ROOT
      INTEGER ILOC_ROOT, JLOC_ROOT, IGLOB, JGLOB
      IF ( .NOT. INVERT ) THEN
        DO ISUB = 1, NSUBSET_ROW
          I         = SUBSET_ROW( ISUB )
          IGLOB     = INDROW_SON( I )
          IPOS_ROOT = RG2L_ROW( IGLOB )
          ILOC_ROOT = MBLOCK
     *            * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *            + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
          DO JSUB = 1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
            JPOS_ROOT = RG2L_COL( JGLOB )
            JLOC_ROOT = NBLOCK
     *                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     *            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
          END DO
        END DO
      ELSE
        DO ISUB = 1, NSUBSET_ROW
          I         = SUBSET_ROW( ISUB )
          IGLOB     = INDROW_SON( I )
          JPOS_ROOT = RG2L_ROW( IGLOB )
          JLOC_ROOT = NBLOCK
     *                * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
          DO JSUB = 1, NSUBSET_COL
            J         = SUBSET_COL( JSUB )
            JGLOB     = INDCOL_SON( J )
            IPOS_ROOT = RG2L_COL( JGLOB )
            ILOC_ROOT = MBLOCK
     *                * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *                + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) =
     *            VAL_ROOT( ILOC_ROOT, JLOC_ROOT ) + VAL_SON( J, I )
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_285
      SUBROUTINE ZMUMPS_164
     *( MYID, NPROCS, N, root, COMM_ROOT, IROOT, FILS,
     *  K50, K46, K51
     *     , K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
     * )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INTEGER MYID, MYID_ROOT
      TYPE (ZMUMPS_ROOT_STRUC)::root
      INTEGER COMM_ROOT
      INTEGER N, IROOT, NPROCS, K50, K46, K51
      INTEGER FILS( N )
      INTEGER K60, IDNPROW, IDNPCOL, IDMBLOCK, IDNBLOCK
      INTEGER INODE, NPROWtemp, NPCOLtemp
      LOGICAL SLAVE
      root%ROOT_SIZE     = 0
      root%TOT_ROOT_SIZE = 0
      NULLIFY( root%RG2L_ROW )
      NULLIFY( root%RG2L_COL )
      SLAVE = ( MYID .ne. 0 .or.
     *        ( MYID .eq. 0 .and. K46 .eq. 1 ) )
      IF ( SLAVE ) THEN
        INODE = IROOT
        DO WHILE ( INODE .GT. 0 )
          INODE = FILS( INODE )
          root%ROOT_SIZE = root%ROOT_SIZE + 1
        END DO
      END IF
      IF ( ( K60 .NE. 2 .AND. K60 .NE. 3 ) .OR.
     *       IDNPROW .LE. 0 .OR. IDNPCOL .LE. 0
     *      .OR. IDMBLOCK .LE.0 .OR. IDNBLOCK.LE.0
     *      .OR. IDNPROW * IDNPCOL .GT. NPROCS ) THEN
      root%MBLOCK = K51
      root%NBLOCK = K51
      CALL ZMUMPS_99( NPROCS, root%NPROW, root%NPCOL,
     *                         root%ROOT_SIZE, K50 )
      IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        IDNPROW = root%NPROW
        IDNPCOL = root%NPCOL
        IDMBLOCK = root%MBLOCK
        IDNBLOCK = root%NBLOCK
      ENDIF
      ELSE IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        root%NPROW = IDNPROW
        root%NPCOL = IDNPCOL
        root%MBLOCK = IDMBLOCK
        root%NBLOCK = IDNBLOCK
      ENDIF
      IF  ( K60 .EQ. 2 .OR. K60 .EQ. 3 ) THEN
        IF (SLAVE) THEN
          root%LPIV = 0
          IF (K46.EQ.0) THEN
            MYID_ROOT=MYID-1
          ELSE
            MYID_ROOT=MYID
          ENDIF
          IF (MYID_ROOT < root%NPROW*root%NPCOL) THEN
            root%MYROW = MYID_ROOT / root%NPCOL
            root%MYCOL = MOD(MYID_ROOT, root%NPCOL)
            root%yes  = .true.
          ELSE
            root%MYROW = -1
            root%MYCOL = -1
            root%yes  = .FALSE.
          ENDIF
        ELSE
          root%yes  = .FALSE.
        ENDIF
      ELSE
      IF ( SLAVE ) THEN
        IF ( root%gridinit_done) THEN
           CALL BLACS_GRIDEXIT( root%CNTXT_BLACS )
           root%gridinit_done = .FALSE.
        END IF
        root%CNTXT_BLACS = COMM_ROOT
        CALL BLACS_GRIDINIT( root%CNTXT_BLACS, 'R',
     *                       root%NPROW, root%NPCOL )
        root%gridinit_done = .TRUE.
        CALL BLACS_GRIDINFO( root%CNTXT_BLACS,
     *                       NPROWtemp, NPCOLtemp,
     *                       root%MYROW, root%MYCOL )
        IF ( root%MYROW .NE. -1 ) THEN
          root%yes = .true.
        ELSE
          root%yes = .false.
        END IF
        root%LPIV = 0
      ELSE
        root%yes = .FALSE.
      END IF
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_164
      SUBROUTINE ZMUMPS_165( N, root, FILS, IROOT, INFO )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      TYPE ( ZMUMPS_ROOT_STRUC ):: root
      INTEGER N, IROOT, INFO(40)
      INTEGER FILS( N )
      INTEGER INODE, I, allocok
      IF ( ASSOCIATED( root%RG2L_ROW ) ) DEALLOCATE( root%RG2L_ROW )
      IF ( ASSOCIATED( root%RG2L_COL ) ) DEALLOCATE( root%RG2L_COL )
      ALLOCATE( root%RG2L_ROW( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      ALLOCATE( root%RG2L_COL( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        RETURN
      ENDIF
      INODE = IROOT
      I = 1
      DO WHILE ( INODE .GT. 0 )
        root%RG2L_ROW( INODE ) = I
        root%RG2L_COL( INODE ) = I
        I = I + 1
        INODE = FILS( INODE )
      END DO
      RETURN
      END SUBROUTINE ZMUMPS_165
      SUBROUTINE ZMUMPS_99( NPROCS, NPROW, NPCOL, SIZE, K50 )
      IMPLICIT NONE
      INTEGER NPROCS, NPROW, NPCOL, SIZE, K50
      INTEGER NPROWtemp, NPCOLtemp, NPROCSUSED, FLATNESS
      COMPLEX*16 TWO
      PARAMETER(TWO=2.0D0)
      LOGICAL KEEPIT
      IF ( K50 .EQ. 1 ) THEN
        FLATNESS = 2
      ELSE
        FLATNESS = 3
      ENDIF
      NPROW  = INT(SQRT(DBLE(NPROCS)))
      NPROWtemp = NPROW
      NPCOL  = INT(NPROCS / NPROW)
      NPCOLtemp = NPCOL
      NPROCSused = NPROWtemp * NPCOLtemp
 10   CONTINUE
      IF ( NPROWtemp >= NPCOLtemp/FLATNESS .AND. NPROWtemp > 1) THEN
        NPROWtemp = NPROWtemp - 1
        NPCOLtemp = INT(NPROCS / NPROWtemp)
        KEEPIT=.FALSE.
        IF ( NPROWtemp * NPCOLtemp .GE. NPROCSused ) THEN
          IF ( ( K50 .NE. 1 .AND. NPROWtemp >= NPCOLtemp/FLATNESS)
     *         .OR. NPROWtemp * NPCOLtemp .GT. NPROCSused )
     *         KEEPIT=.TRUE.
        END IF
        IF ( KEEPIT ) THEN
          NPROW = NPROWtemp
          NPCOL = NPCOLtemp
          NPROCSused = NPROW * NPCOL
        END IF
        GO TO 10
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_99
      SUBROUTINE ZMUMPS_290(MYID, M, N, ASEQ,
     *                    LOCAL_M, LOCAL_N,
     *                    MBLOCK, NBLOCK,
     *                    APAR,
     *                    MASTER_ROOT,
     *                    NPROW, NPCOL,
     *                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      COMPLEX*16 APAR( LOCAL_M, LOCAL_N )
      COMPLEX*16 ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, IDEST, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = MOD ( IBLOCK, NPROW )
            ICOL = MOD ( JBLOCK, NPCOL )
            IDEST = IROW * NPCOL + ICOL
            IF ( IDEST .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  WK(KK)=ASEQ(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK, SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         IDEST, 128, COMM, IERR )
              ELSE IF ( MYID .EQ. IDEST ) THEN
                CALL MPI_RECV( WK(1),
     *                         SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         MASTER_ROOT,128,COMM,STATUS,IERR)
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  APAR(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     *              JAPAR:JAPAR+SIZE_JBLOCK-1 )
     *        = ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END SUBROUTINE ZMUMPS_290
      SUBROUTINE ZMUMPS_156(MYID, M, N, ASEQ,
     *                    LOCAL_M, LOCAL_N,
     *                    MBLOCK, NBLOCK,
     *                    APAR,
     *                    MASTER_ROOT,
     *                    NPROW, NPCOL,
     *                    COMM)
      IMPLICIT NONE
      INTEGER MYID, MASTER_ROOT, COMM
      INTEGER M, N
      INTEGER NPROW, NPCOL
      INTEGER LOCAL_M, LOCAL_N
      INTEGER MBLOCK, NBLOCK
      COMPLEX*16 APAR( LOCAL_M, LOCAL_N )
      COMPLEX*16 ASEQ( M, N )
      INCLUDE 'mpif.h'
      INTEGER I, J, SIZE_IBLOCK, SIZE_JBLOCK, ISOUR, IROW, ICOL
      INTEGER IBLOCK, JBLOCK, II, JJ, KK
      INTEGER IAPAR, JAPAR, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMPLEX*16 WK( MBLOCK * NBLOCK )
      LOGICAL JUPDATE
        IAPAR = 1
        JAPAR = 1
        DO J = 1, N, NBLOCK
          SIZE_JBLOCK = NBLOCK
          IF ( J + NBLOCK > N ) THEN
            SIZE_JBLOCK = N - J + 1
          END IF
          JUPDATE = .FALSE.
          DO I = 1, M, MBLOCK
            SIZE_IBLOCK = MBLOCK
            IF ( I + MBLOCK > M ) THEN
              SIZE_IBLOCK = M - I + 1
            END IF
            IBLOCK = I / MBLOCK
            JBLOCK = J / NBLOCK
            IROW = MOD ( IBLOCK, NPROW )
            ICOL = MOD ( JBLOCK, NPCOL )
            ISOUR = IROW * NPCOL + ICOL
            IF ( ISOUR .NE. MASTER_ROOT ) THEN
              IF ( MYID .EQ. MASTER_ROOT ) THEN
                CALL MPI_RECV( WK(1), SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         ISOUR, 128, COMM, STATUS, IERR )
                KK=1
                DO JJ=J,J+SIZE_JBLOCK-1
                DO II=I,I+SIZE_IBLOCK-1
                  ASEQ(II,JJ)=WK(KK)
                  KK=KK+1
                END DO
                END DO
              ELSE IF ( MYID .EQ. ISOUR ) THEN
                KK=1
                DO JJ=JAPAR,JAPAR+SIZE_JBLOCK-1
                DO II=IAPAR,IAPAR+SIZE_IBLOCK-1
                  WK(KK)=APAR(II,JJ)
                  KK=KK+1
                END DO
                END DO
                CALL MPI_SSEND( WK( 1 ),
     *                         SIZE_IBLOCK*SIZE_JBLOCK,
     *                         MPI_DOUBLE_COMPLEX,
     *                         MASTER_ROOT,128,COMM,IERR)
                JUPDATE = .TRUE.
                IAPAR = IAPAR + SIZE_IBLOCK
              END IF
            ELSE IF ( MYID.EQ. MASTER_ROOT ) THEN
              ASEQ(I:I+SIZE_IBLOCK-1,J:J+SIZE_JBLOCK-1)
     *        = APAR( IAPAR:IAPAR+SIZE_IBLOCK-1,
     *                JAPAR:JAPAR+SIZE_JBLOCK-1 )
              JUPDATE = .TRUE.
              IAPAR = IAPAR + SIZE_IBLOCK
            END IF
          END DO
          IF ( JUPDATE ) THEN
            IAPAR = 1
            JAPAR = JAPAR + SIZE_JBLOCK
          END IF
        END DO
      RETURN
      END SUBROUTINE ZMUMPS_156
      SUBROUTINE ZMUMPS_284(root, IROOT, N,
     *                  IW, LIW, A, LA,
     *                  NIRBDU,
     *                  FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *                  LRLU, IPTRLU,
     *                  IWPOS, IWPOSCB, PTRIST, PTRAST,
     *                  STEP, PIMASTER, PAMASTER, ITLOC,
     *                  COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      TYPE (ZMUMPS_ROOT_STRUC ) :: ROOT
      INTEGER IROOT, LIW, N, LA, NIRBDU, LRLU, IPTRLU, IWPOS, IWPOSCB
      INTEGER IW( LIW )
      COMPLEX*16 A( LA )
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), STEP(N)
      INTEGER 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER ITLOC( N )
      INTEGER COMP, LRLUS, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER FILS( N ), PTRAIW(N), PTRARW( N )
      INTEGER INTARR(MAX(1,KEEP(14)))
      COMPLEX*16 DBLARR(MAX(1,KEEP(13)))
      INTEGER NUMROC
      EXTERNAL NUMROC
      INTEGER LREQI_ROOT, LREQA_ROOT, LOCAL_M, LOCAL_N
            LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *                root%MYROW, 0, root%NPROW )
            LOCAL_M = MAX( 1, LOCAL_M )
            LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *                root%MYCOL, 0, root%NPCOL )
            LREQI_ROOT = 2 + XSIZE
            LREQA_ROOT = LOCAL_M * LOCAL_N
            IF (LREQA_ROOT.EQ.0) THEN
              PTRIST(STEP(IROOT)) = -9999999
              RETURN
            ENDIF
            CALL ZMUMPS_22(.FALSE.,.FALSE.,
     *                     MYID,N,KEEP,KEEP8,IW,LIW,A,LA, NIRBDU,
     *                     LRLU, IPTRLU,
     *                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     *                     STEP, PIMASTER, PAMASTER, ITLOC, LREQI_ROOT,
     *                     LREQA_ROOT, IROOT, S_NOTFREE, .TRUE., COMP,
     *                     LRLUS, IFLAG, IERROR )
            IF ( IFLAG .LT. 0 ) RETURN
            PTRIST  ( STEP(IROOT) ) = IWPOSCB + 1
            PAMASTER( STEP(IROOT) ) = IPTRLU  + 1
            IW( IWPOSCB + 1 + XSIZE) = - LOCAL_N
            IW( IWPOSCB + 2 + XSIZE) =   LOCAL_M
      RETURN
      END SUBROUTINE ZMUMPS_284
      SUBROUTINE ZMUMPS_35( N, root, IROOT,
     *   VALROOT, LOCAL_M, LOCAL_N, FILS,
     *       PTRAIW, PTRARW,
     *       INTARR, DBLARR,
     *       KEEP,KEEP8,
     *       MYID)
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      TYPE (ZMUMPS_ROOT_STRUC) :: ROOT
      INTEGER N, MYID, IROOT, LOCAL_M, LOCAL_N, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER FILS( N ), PTRARW( N ), PTRAIW( N )
      INTEGER INTARR(MAX(1,KEEP(14)))
      COMPLEX*16 VALROOT(LOCAL_M,LOCAL_N)
      COMPLEX*16 DBLARR(MAX(1,KEEP(13)))
      COMPLEX*16 VAL
      INTEGER IORG, IBROT, NUMORG, JJ, J1,JK, J2,J3, J4,
     *        IROW, JCOL, AINPUT
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INTEGER ILOCROOT, JLOCROOT
      NUMORG = root%ROOT_SIZE
      IBROT  = IROOT
      DO IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        IBROT = FILS(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        JCOL = INTARR(J1)
        DO JJ = J1, J2
         IROW = INTARR(JJ)
         VAL  = DBLARR(AINPUT)
         AINPUT = AINPUT + 1
         IPOSROOT = root%RG2L_ROW( IROW )
         JPOSROOT = root%RG2L_COL( JCOL )
         IROW_GRID  = MOD( ( IPOSROOT - 1 ) / root%MBLOCK, root%NPROW )
         JCOL_GRID  = MOD( ( JPOSROOT - 1 ) / root%NBLOCK, root%NPCOL )
         IF ( IROW_GRID .EQ. root%MYROW .AND.
     *        JCOL_GRID .EQ. root%MYCOL ) THEN
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
            VALROOT( ILOCROOT, JLOCROOT ) = VAL
         END IF
        END DO
        IF (J3 .LE. J4) THEN
         IROW =  INTARR(J1)
         DO JJ= J3,J4
          JCOL = INTARR(JJ)
          VAL  = DBLARR(AINPUT)
          AINPUT = AINPUT + 1
          IPOSROOT = root%RG2L_ROW( IROW )
          JPOSROOT = root%RG2L_COL( JCOL )
          IROW_GRID= MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW)
          JCOL_GRID= MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL)
          IF ( IROW_GRID .EQ. root%MYROW .AND.
     *        JCOL_GRID .EQ. root%MYCOL ) THEN
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
            VALROOT( ILOCROOT, JLOCROOT ) = VAL
          END IF
         END DO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_35
