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
      SUBROUTINE DMUMPS_250(N, NZ, NSTEPS, A, LA,
     * IW, LIW, PERM, NA, LNA, NE_STEPS,
     * NSTK_STEPS, NBPROCFILS,
     * MAXFRT, ND, FILS, STEP, FRERE_STEPS, DAD, CAND,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * PTRARW, PTRAIW,
     * PTRAST, PTRIST, PTLUST_S, PTRFAC, PIMASTER, PAMASTER, ITLOC, 
     * NIRBDU, IPOOL, LPOOL,
     * UULOC, ICNTL, INFO, RINFO, KEEP,KEEP8,PROCNODE_STEPS,
     * SLAVEF,COMM_NODES,MYID,MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,
     * INTARR,DBLARR,root, NELT, FRTPTR, 
     * FRTELT, LPTRAR, COMM_LOAD, ASS_IRECV, SEUIL,
     * SEUIL_LDLT_NIV2, MEM_DISTRIB,
     * DKEEP,PIVNUL_LIST,LPN_LIST)
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER N, LIW, LA, LPTRAR, LNA
      INTEGER NZ, NSTEPS, LPOOL, KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION A(LA)
      INTEGER MAXFRT, INFO(40), NIRBDU
      INTEGER SLAVEF,COMM_NODES,MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IW(LIW)
      INTEGER PERM(N), NA(LNA),
     *        NE_STEPS(KEEP(28)), NSTK_STEPS(KEEP(28)),
     *        DAD(KEEP(28)), STEP(N),
     *        ND(KEEP(28)), FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR),PTLUST_S(KEEP(28)),
     *        PTRFAC(KEEP(28)),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        PTRAST(KEEP(28)), PTRIST(KEEP(28)),
     *        ITLOC(N), NBPROCFILS( KEEP(28) )
      INTEGER   CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER IPOOL(LPOOL), PROCNODE_STEPS(KEEP(28))
      INTEGER NELT, COMM_LOAD, ASS_IRECV
      INTEGER FRTPTR(N+1), FRTELT(NELT)
      DOUBLE PRECISION    RINFO(20)
      INTEGER ICNTL(40)
      DOUBLE PRECISION UULOC
      DOUBLE PRECISION SEUIL, SEUIL_LDLT_NIV2
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER INTARR(MAX(1,KEEP(14)))
      DOUBLE PRECISION DBLARR(MAX(1,KEEP(13)))
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER MPRINT,LP
      INTEGER I, NBLEAF
      INTEGER INODE, IERR, NTOTPVTOT
      INTEGER NTOTPV
      INTEGER LRLU, IPTRLU, LRLUS,
     *    LEAF, NBROOT, NROOT
      INTRINSIC INT,REAL,LOG
      MPRINT = ICNTL(2)
      LP     = ICNTL(1)
      IF (KEEP(4).LE.0) KEEP(4)=32
      IF (KEEP(5).LE.0) KEEP(5)=16
      IF (KEEP(5).GT.KEEP(4)) KEEP(5) = KEEP(4)
      IF (KEEP(6).LE.0) KEEP(6)=24
      IF (KEEP(3).LE.KEEP(4)) KEEP(3)=KEEP(4)*2
      IF (KEEP(6).GT.KEEP(3)) KEEP(6) = KEEP(3)
      INFO(9)   = 1
      INFO(10)  = 1
      MAXFRT = 0
        LRLU     = LA
        LRLUS    = LRLU
        KEEP(67) = LRLUS
        IPTRLU   = LRLU
      KEEP(18) = LRLU
      NTOTPV = 0
      NSTK_STEPS(1:KEEP(28)) = NE_STEPS(1:KEEP(28))
      CALL DMUMPS_362(N, LEAF, NBROOT, NROOT,
     *                     MYID_NODES,
     *                     SLAVEF, NA, LNA,
     *                     KEEP,KEEP8, STEP,
     *                     PROCNODE_STEPS,
     *                     IPOOL, LPOOL)
      CALL DMUMPS_506(IPOOL, LPOOL, LEAF)
      CALL DMUMPS_555(IPOOL, LPOOL,KEEP,KEEP8)
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        NBROOT = NBROOT + root%NPROW * root%NPCOL - 1
      END IF
      IF ( root%yes )  THEN 
         IF ( DMUMPS_275( STEP(KEEP(38)),PROCNODE_STEPS, SLAVEF)
     *         .NE. MYID_NODES ) THEN
             NROOT = NROOT + 1
         END IF
      END IF
      CALL DMUMPS_251(N,IW,LIW,A,LA,NSTK_STEPS,NBPROCFILS,
     *         INFO(1),ND,FILS,STEP,FRERE_STEPS, DAD, CAND,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *         MAXFRT, NTOTPV,PTRIST,PTRAST,
     *         PIMASTER, PAMASTER, PTRARW,PTRAIW,
     *         ITLOC, INFO(2), IPOOL, LPOOL, 
     *         NIRBDU, RINFO,
     *         INFO(9),INFO(10),LRLU,IPTRLU, 
     *         LRLUS, LEAF, NROOT, NBROOT,
     *         UULOC,ICNTL,PTLUST_S,PTRFAC,NSTEPS,INFO,
     *         KEEP,KEEP8,
     *         PROCNODE_STEPS,SLAVEF,MYID,COMM_NODES,
     *         MYID_NODES,
     *         BUFR,LBUFR, LBUFR_BYTES,
     *         INTARR, DBLARR, root, PERM,
     *         NELT, FRTPTR, FRTELT, LPTRAR, 
     *         COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     *         MEM_DISTRIB(0),NE_STEPS,
     *     DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      INFO(9)  = INFO(9)  - 1
      INFO(10) = INFO(10) - 1
      IF (KEEP(201).EQ.0) THEN
        KEEP(31) = INFO(9)
      ENDIF
      KEEP(32) = INFO(10)
      KEEP(67) = LA - KEEP(67)
      KEEP(89) = NTOTPV
      INFO(23) = KEEP(89)
      CALL MPI_ALLREDUCE(NTOTPV, NTOTPVTOT, 1, MPI_INTEGER, MPI_SUM, 
     *                COMM_NODES, IERR)
      IF ( ( (INFO(1).EQ.-10) .and. (NTOTPVTOT.EQ.N)) 
     *              .OR. ( NTOTPVTOT.GT.N ) ) THEN
       write(*,*) ' Error 1 in mc51d NTOTPVTOT=', NTOTPVTOT 
       CALL DMUMPS_ABORT()
      endif
      IF ( (KEEP(19).NE.0 ) .and. (NTOTPVTOT.NE.N))  THEN
       write(*,*) ' Error 2 in mc51d NTOTPVTOT=', NTOTPVTOT 
       CALL DMUMPS_ABORT()
      endif
      IF ( (INFO(1) .GE. 0 ) 
     *      .AND. (NTOTPVTOT.NE.N) ) THEN
         INFO(1) = -10
         INFO(2) = NTOTPVTOT
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_250
!      
      SUBROUTINE DMUMPS_26(id)
      USE DMUMPS_LOAD
      USE DMUMPS_CV
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      TYPE(DMUMPS_STRUC), TARGET :: id
      INTEGER LIW, IKEEP, FILS, FRERE, PTRAR, NFSIZ
      INTEGER NE, NA
      INTEGER I, ALLOCOK
      INTEGER MAXIS1_CHECK
      INTEGER INDX, NRECV, NB_NIV2, LAST_ALLOWED, IDEST
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER LOCAL_M, LOCAL_N, K44MAX
      INTEGER NUMROC
      EXTERNAL NUMROC
      INTEGER IRANK
      INTEGER MP, LP, MPG
      LOGICAL PROK, PROKG, LISTVAR_SCHUR_2BE_FREED
      INTEGER SIZE_SCHUR_PASSED
      INTEGER SBUF_SEND, SBUF_REC, TOTAL_MBYTES
      INTEGER NSLAVES_SAVE
      INTEGER LEAF, INODE, ISTEP, INN, LPTRAR
      INTEGER NBLEAF, NBROOT, MYROW_CHECK, INIV2
      INTEGER K16TMP
      INTEGER MM_WRITE
      INTEGER MM_WRITE_CHECK
      CHARACTER(LEN=20) MM_IDSTR
      CHARACTER (LEN=10) :: symm
      CHARACTER (LEN=8) :: prec
      DOUBLE PRECISION PEAK
      INTEGER PERM
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: REQPTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SSARBR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PAR2_NODES
      INTEGER, DIMENSION(:), POINTER :: IWtemp 
      INTEGER, DIMENSION(:), POINTER :: XNODEL, NODEL
      INTEGER, POINTER ::  NELT, LELTVAR
      INTEGER, DIMENSION(:), POINTER :: KEEP,INFO, INFOG
      INTEGER*8, DIMENSION(:), POINTER :: KEEP8
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFO
      DOUBLE PRECISION, DIMENSION(:), POINTER :: RINFOG
      INTEGER, DIMENSION(:), POINTER :: ICNTL
      LOGICAL I_AM_SLAVE, OOC_STAT, PERLU_ON
      INTEGER DMUMPS_330, DMUMPS_275
      EXTERNAL DMUMPS_330, DMUMPS_275
      INTEGER K,J, IORD, IFS
      INTEGER SIZE_TEMP_MEM,SIZE_DEPTH_FIRST,SIZE_COST_TRAV
      LOGICAL IS_REORDERTREE_CALLED
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_MEM
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_ROOT
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_LEAF
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: TEMP_SIZE
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH_FIRST
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: COST_TRAV_TMP
      INTEGER*8 :: TOTAL_BYTES
      INTEGER LOC_DIM_BUF_IO,MAX_SIZE_FACTOR
      LOC_DIM_BUF_IO=0
      IS_REORDERTREE_CALLED=.FALSE.
      KEEP   => id%KEEP
      KEEP8  => id%KEEP8
      INFO   => id%INFO
      RINFO  => id%RINFO
      INFOG  => id%INFOG
      RINFOG => id%RINFOG
      ICNTL  => id%ICNTL
      NELT    => id%NELT
      LELTVAR => id%LELTVAR
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               id%KEEP(46) .eq. 1 ) )
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IORD = ICNTL(7)
      IF (PROK) WRITE( MP, 220 )
      IF (PROKG.AND.(MPG .NE. MP)) WRITE( MPG, 220 ) id%VERSION_NUMBER
 220  FORMAT( /' DMUMPS Version ',A6,' -- Jan 2007' )
      IF ( PROK ) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          WRITE(MP, '(A)') 'L U Solver for unsymmetric matrices'
        ELSE IF ( KEEP(50) .eq. 1 ) THEN
          WRITE(MP, '(A)') 
     * 'L D L^T Solver for symmetric positive definite matrices'
        ELSE
          WRITE(MP, '(A)') 
     * 'L D L^T Solver for general symmetric matrices'
        END IF
        IF ( KEEP(46) .eq. 1 ) THEN
          WRITE(MP, '(A)') 'Type of parallelism: Working host'
        ELSE
          WRITE(MP, '(A)') 'Type of parallelism: Host not working'
        END IF
      END IF
      IF ( PROKG .AND. (MP.NE.MPG)) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          WRITE(MPG, '(A)') 'L U Solver for unsymmetric matrices'
        ELSE IF ( KEEP(50) .eq. 1 ) THEN
          WRITE(MPG, '(A)') 
     * 'L D L^T Solver for symmetric positive definite matrices'
        ELSE
          WRITE(MPG, '(A)') 
     * 'L D L^T Solver for general symmetric matrices'
        END IF
        IF ( KEEP(46) .eq. 1 ) THEN
          WRITE(MPG, '(A)') 'Type of parallelism: Working host'
        ELSE
          WRITE(MPG, '(A)') 'Type of parallelism: Host not working'
        END IF
      END IF
      IF (PROK) WRITE( MP, 110 )
      IF (PROKG .AND. (MPG.NE.MP)) WRITE( MPG, 110 )
      IF ( (id%KEEP(24).NE.0) .AND.
     *     id%NSLAVES.eq.1 ) THEN
        id%KEEP(24) = 0
        IF ( PROKG ) THEN
          WRITE(MPG, '(A)')
     *    ' Resetting candidate strategy to 0 because NSLAVES=1'
          WRITE(MPG, '(A)') ' '
        END IF
      END IF
      IF ( (id%KEEP(24).EQ.0) .AND.
     *     id%NSLAVES.GT.1 ) THEN
         id%KEEP(24) = 8
      ENDIF
      IF ( (id%KEEP(24).NE.0) .AND.(id%KEEP(24).NE.1).AND.
     *     (id%KEEP(24).NE.8).AND.(id%KEEP(24).NE.10)
     *  .AND.(id%KEEP(24).NE.12).AND.(id%KEEP(24).NE.14)
     *  .AND.(id%KEEP(24).NE.16).AND.(id%KEEP(24).NE.18)) THEN
       id%KEEP(24) = 8
        IF ( PROKG ) THEN
          WRITE(MPG, '(A)')
     *    ' Resetting candidate strategy to 8 '
          WRITE(MPG, '(A)') ' '
        END IF
      END IF
      KEEP(10)=KEEP(85)
      IF ( id%MYID.EQ.MASTER ) THEN
        KEEP(201)=ICNTL(22)
        IF ( IORD .EQ.1 ) THEN
          IF ( .not. ASSOCIATED( id%PERM_IN ) ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 3
            GO TO 333
          ELSE IF ( SIZE( id%PERM_IN ) < id%N ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 3
            GO TO 333
          END IF
        ENDIF
        IF (KEEP(9) .LE. 1 ) KEEP(9) = 500
        IF ( KEEP(10) .GT. 0 ) THEN 
         IF ((KEEP(10).LE.1) .OR. (KEEP(10).GT.KEEP(9))) 
     *      KEEP(10) = MIN(KEEP(9),100)
        ENDIF
        IF (KEEP(48). EQ. 1 ) KEEP(48) = -12345
        IF ( (KEEP(48).LT.0) .OR. (KEEP(48).GT.5) ) THEN
          KEEP(48)=5
        ENDIF
        KEEP(60) = ICNTL(19)
        IF ( KEEP( 60 ) .LE. 0 ) KEEP( 60 ) = 0
        IF ( KEEP( 60 ) .GT. 3 ) KEEP( 60 ) = 0
        IF ( KEEP(60) .NE.0 ) THEN
         IF ( .not. ASSOCIATED( id%LISTVAR_SCHUR ) ) THEN
            id%INFO(1) = -22
            id%INFO(2) = 8
            GO TO 333
         ELSE IF (SIZE(id%LISTVAR_SCHUR)<id%SIZE_SCHUR) THEN
            id%INFO(1) = -22
            id%INFO(2) = 8
            GO TO 333
         END IF
        ENDIF
        IF (KEEP(60) == 3 .AND. KEEP(50).NE.0) THEN
          IF (id%MBLOCK > 0 .and. id%NBLOCK > 0 .and.
     *        id%NPROW > 0 .and. id%NPCOL > 0 ) THEN
            IF (id%NPROW *id%NPCOL .le. id%NSLAVES) THEN
              IF (id%MBLOCK .NE. id%NBLOCK ) THEN
                id%INFO(1)=-31
                id%INFO(2)=id%MBLOCK - id%NBLOCK
                GOTO 333
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF ( KEEP(60) .ne. 0 .and. IORD .eq. 1 ) THEN
          DO I = 1,id%SIZE_SCHUR
            IF (id%PERM_IN(id%LISTVAR_SCHUR(I))
     *        ==id%N-id%SIZE_SCHUR+I)
     *        CYCLE
            IORD = 0
            IF (MPG.GT.0) THEN
              WRITE(MPG,'(A)')
     *    ' ** Ignoring user-ordering, because incompatible with Schur.'
              WRITE(MPG,'(A)') ' ** ICNTL(7) treated as 0.'
            END IF
            EXIT
          ENDDO
        END IF
        KEEP(23) = ICNTL(6)
        KEEP(95) = ICNTL(12)
        IF (KEEP(50).NE.2) KEEP(95) = 1
        IF ((KEEP(95).GT.3).OR.(KEEP(95).LT.0)) KEEP(95) = 0
        KEEP(54) = ICNTL(18)
        IF ( KEEP(54) .LT. 0 .or. KEEP(54).GT.3 ) THEN
          IF ( PROKG ) THEN
             WRITE(MPG, *) ' Out-of-range value for ICNTL(18).'
             WRITE(MPG, *) ' Used 0 ie matrix not distributed'
          END IF
          KEEP(54) = 0
        END IF
        KEEP(55) = ICNTL(5)
        IF ( KEEP(55) .LT. 0 .OR. KEEP(55) .GT. 1 ) THEN
           IF ( PROKG ) THEN
              WRITE(MPG, *) ' Out-of-range value for ICNTL(5).'
              WRITE(MPG, *) ' Used 0 ie matrix is assembled'
           END IF
           KEEP(55) = 0
        END IF
        IF (KEEP(23).LT.0.OR.KEEP(23).GT.7) KEEP(23) = 7
        IF ( KEEP(50) .EQ. 1 ) THEN
           IF (KEEP(23) .NE. 0) THEN
              IF (MPG.GT.0) THEN
                 WRITE(MPG,'(A)')
     *             ' ** Max-trans not compatible with LLT factorization'
              END IF
              KEEP(23) = 0
           ENDIF
           KEEP(95) = 1
        END IF
        IF  (KEEP(60) .GT. 0) THEN
          IF (KEEP(23) .NE. 0) THEN
             IF (MPG.GT.0) THEN
                WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because of Schur'
             END IF
             KEEP(23) = 0
          ENDIF
          IF (MPG.GT.0) THEN
             WRITE(MPG,'(A)')
     *  ' ** Scaling not allowed because of Schur'
             KEEP(52) = 0
          ENDIF
          KEEP(95) = 1
        END IF
        IF ( (KEEP(23) .NE. 0) .AND. (IORD .EQ. 1) ) THEN
          KEEP(23) = 0
          KEEP(95) = 1
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because ordering is given'
          END IF
        END IF
        IF (KEEP(54) .NE. 0) THEN
           IF( KEEP(23) .NE. 0 ) THEN
              IF (MPG.GT.0) THEN
                 WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed because matrix is distributed'
              END IF
              KEEP(23) = 0
           ENDIF
           IF (MPG.GT.0) THEN
              WRITE(MPG,'(A)')
     *  ' ** Scaling not allowed because matrix is distributed'
              KEEP(52) = 0
           ENDIF
           KEEP(95) = 1
        END IF
        IF ( KEEP(55) .NE. 0 ) THEN
           IF( KEEP(23) .NE. 0 ) THEN
              IF (MPG.GT.0) THEN
                 WRITE(MPG,'(A)')
     *  ' ** Max-trans not allowed for element matrix'
              END IF
              KEEP(23) = 0
           ENDIF
           IF (MPG.GT.0) THEN
              WRITE(MPG,'(A)')
     *  ' ** Scaling not allowed for element matrix'
              KEEP(52) = 0
           ENDIF
           KEEP(95) = 1
        ENDIF
        IF ( KEEP(54) .ne. 0 .AND. KEEP(55) .ne. 0 ) THEN
          KEEP(54) = 0
          IF (MPG.GT.0) THEN
            WRITE(MPG,'(A)')
     *  ' ** Distributed entry not available for element matrix'
          END IF
        ENDIF
      END IF
      IF(KEEP(50) .EQ. 2) THEN
         IF( .NOT. ASSOCIATED(id%A) ) THEN
            IF(KEEP(95) .EQ. 3) THEN
               KEEP(95) = 2
            ENDIF
         ENDIF
         IF(KEEP(95) .EQ. 3 .AND. IORD .NE. 2) THEN
            IF (PROK) WRITE(MP,*)
     *        'WARNING: ANAL_0 constrained ordering not available ', 
     *        'with selected ordering'
            KEEP(95) = 2
         ENDIF 
         IF(KEEP(95) .EQ. 2 .AND. IORD .EQ. 0) THEN 
            IORD = 6
            ICNTL(7) = 6
         ENDIF
         IF(KEEP(95) .EQ. 3) THEN
            KEEP(23) = 5
            KEEP(52) = -2
         ELSE IF(KEEP(95) .EQ. 2 .AND. 
     *         (KEEP(23) .EQ. 0 .OR. KEEP(23) .EQ. 7) ) THEN
            IF( ASSOCIATED(id%A) ) THEN
               KEEP(23) = 5
            ELSE
               KEEP(23) = 1
            ENDIF
         ELSE IF(KEEP(95) .EQ. 1) THEN
            KEEP(23) = 0
         ELSE IF(KEEP(95) .EQ. 0 .AND. KEEP(23) .EQ. 0) THEN
            KEEP(95) = 1
         ENDIF
      ELSE
         KEEP(95) = 1
      ENDIF
 333  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                  id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( KEEP(60), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF (id%KEEP(60) .EQ. 2 .or. id%KEEP(60). EQ. 3) THEN
        CALL MPI_BCAST( id%NPROW, 1,
     *                 MPI_INTEGER, MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%NPCOL, 1,
     *                 MPI_INTEGER, MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%MBLOCK, 1,
     *                 MPI_INTEGER, MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%NBLOCK, 1,
     *                 MPI_INTEGER, MASTER, id%COMM, IERR )
      ENDIF
      CALL DMUMPS_276( ICNTL, INFO,
     *                  id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( KEEP(54), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(69), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER, MASTER, id%COMM, IERR )
      IF ( ASSOCIATED(id%MEM_DIST) ) DEALLOCATE( id%MEM_DIST )
      ALLOCATE( id%MEM_DIST( 0:id%NSLAVES-1 ), STAT=IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO(1) = -7
         INFO(2) = id%NSLAVES
         IF ( LP .GT. 0 ) THEN
           WRITE(LP, 150) 'MEM_DIST'
         END IF
         RETURN
      END IF
      id%MEM_DIST(0:id%NSLAVES-1) = 0
      CALL DMUMPS_427(
     *     id%COMM,id%COMM_NODES,KEEP(69),KEEP(46),
     *     id%NSLAVES,id%MEM_DIST,INFO)
      IF ( KEEP(54) .eq. 3 ) THEN
          IF ( KEEP(46) .eq. 0 .and. id%MYID .eq. MASTER ) THEN
            id%NZ_loc = 0
          END IF
          CALL MPI_REDUCE( id%NZ_loc, id%NZ, 1, MPI_INTEGER, MPI_SUM,
     *                     MASTER, id%COMM, IERR )
          IF ( id%MYID .eq. MASTER ) THEN
            ALLOCATE( REQPTR( id%NPROCS, 3 ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -7
              id%INFO(2) = 3 * id%NPROCS
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'REQPTR'
              END IF
              GOTO 13
            END IF
            ALLOCATE( id%IRN( id%NZ ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -7
              id%INFO(2) = id%NZ
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'IRN'
              END IF
              GOTO 13
            END IF
            ALLOCATE( id%JCN( id%NZ ), stat = allocok )
            IF ( allocok .GT. 0 ) THEN
              id%INFO(1) = -7
              id%INFO(2) = id%NZ
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'JCN'
              END IF
              GOTO 13
            END IF
          END IF
 13       CONTINUE
          CALL DMUMPS_276( ICNTL, INFO,
     *                              id%COMM, id%MYID )
          IF ( INFO(1) < 0 ) RETURN
          IF ( id%MYID .eq. MASTER ) THEN
            DO I = 1, id%NPROCS - 1
              CALL MPI_RECV( REQPTR( I+1, 1 ), 1, 
     *        MPI_INTEGER, I,
     *        COLLECT_NZ, id%COMM, STATUS, IERR )
            END DO
            IF ( KEEP(46) .eq. 0 ) THEN
              REQPTR( 1, 1 ) = 1
            ELSE
              REQPTR( 1, 1 ) = id%NZ_loc + 1
            END IF
            DO I = 2, id%NPROCS
              REQPTR( I, 1 ) = REQPTR( I, 1 ) + REQPTR( I-1, 1 )
            END DO
          ELSE
            CALL MPI_SEND( id%NZ_loc, 1, MPI_INTEGER, MASTER,
     *        COLLECT_NZ, id%COMM, IERR )
          END IF
          IF ( id%MYID .eq. MASTER ) THEN
            NRECV = 0
            DO I = 1, id%NPROCS - 1
              IF ( REQPTR( I + 1, 1 ) - REQPTR( I, 1 ) .NE. 0 ) THEN
                NRECV = NRECV + 2
                CALL MPI_IRECV( id%IRN( REQPTR( I, 1 ) ),
     *          REQPTR( I + 1, 1 ) - REQPTR( I, 1 ), 
     *          MPI_INTEGER,
     *          I, COLLECT_IRN, id%COMM, REQPTR(I, 2), IERR )
                CALL MPI_IRECV( id%JCN( REQPTR( I, 1 ) ),
     *          REQPTR( I + 1, 1 ) - REQPTR( I, 1 ),   
     *          MPI_INTEGER,
     *          I, COLLECT_JCN, id%COMM, REQPTR(I, 3), IERR )
              ELSE
                REQPTR(I, 2) = MPI_REQUEST_NULL
                REQPTR(I, 3) = MPI_REQUEST_NULL
              END IF
            END DO
          ELSE
            IF ( id%NZ_loc .NE. 0 ) THEN
              CALL MPI_SEND( id%IRN_loc(1), id%NZ_loc, 
     *          MPI_INTEGER, MASTER,
     *          COLLECT_IRN, id%COMM, IERR )
              CALL MPI_SEND( id%JCN_loc(1), id%NZ_loc, 
     *          MPI_INTEGER, MASTER,
     *          COLLECT_JCN, id%COMM, IERR )
            END IF
          END IF
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( id%NZ_loc .NE. 0 ) THEN
              DO I=1,id%NZ_loc
                id%IRN(i) = id%IRN_loc(i)
                id%JCN(i) = id%JCN_loc(i)
              ENDDO
            END IF
            REQPTR( id%NPROCS, 2 ) = MPI_REQUEST_NULL
            REQPTR( id%NPROCS, 3 ) = MPI_REQUEST_NULL
            DO I = 1, NRECV
              CALL MPI_WAITANY
     *        ( 2 * id%NPROCS, REQPTR( 1, 2 ), INDX, STATUS, IERR )
            END DO
            DEALLOCATE( REQPTR )
          END IF
      END IF
      IF ( id%MYID .eq. MASTER .AND. KEEP(54).NE.3) THEN
        IF (id%WRITE_PROBLEM(1:20) .NE. "NAME_NOT_INITIALIZED")THEN
          IF (ASSOCIATED(id%A)) THEN
      prec='real '
          ELSE
            prec='pattern '
          ENDIF
          IF (id%SYM .eq. 0) THEN
            symm="general"
          ELSE
            symm="symmetric"
          END IF
          OPEN(69,FILE=TRIM(id%WRITE_PROBLEM))
          WRITE(69,FMT=*)'%%MatrixMarket matrix coordinate ',TRIM(prec),
     &          " ",TRIM(symm)
          WRITE(69,*) id%N, id%N, id%NZ
          IF (ASSOCIATED(id%A)) THEN
            DO I=1,id%NZ
          WRITE(69,*) id%IRN(I), id%JCN(I), id%A(I)  
            ENDDO
          ELSE
            DO I=1,id%NZ
              WRITE(69,*) id%IRN(I), id%JCN(I)
            ENDDO
          ENDIF
          CLOSE(69)
          OPEN(69,FILE=TRIM(id%WRITE_PROBLEM) //".rhs")
          WRITE(69,FMT=*)'%%MatrixMarket matrix array ',TRIM(prec),
     &          ' general'
          WRITE(69,*) id%N, 1
          IF (ASSOCIATED(id%RHS)) THEN
            DO I=1,id%N
              WRITE(69,*) id%RHS(I)
            ENDDO
          ENDIF
          CLOSE(69)
        END IF
      ENDIF
      IF ( KEEP(54).EQ.3) THEN
        IF (id%WRITE_PROBLEM(1:20) .EQ. "NAME_NOT_INITIALIZED"
     *    .OR. .NOT. I_AM_SLAVE )THEN
          MM_WRITE = 0
        ELSE
          MM_WRITE = 1
        ENDIF
        CALL MPI_ALLREDUCE(MM_WRITE, MM_WRITE_CHECK, 1,
     *         MPI_INTEGER, MPI_SUM, id%COMM, IERR)
        IF (MM_WRITE_CHECK.EQ.id%NSLAVES .AND. I_AM_SLAVE) THEN
          IF (ASSOCIATED(id%A_loc)) THEN
        prec='real '
          ELSE
              prec='pattern '
          ENDIF
          IF (id%SYM .eq. 0) THEN
            symm="general"
          ELSE
            symm="symmetric"
          END IF
          WRITE(MM_IDSTR,'(I7)') id%MYID
          OPEN(69,FILE=TRIM(id%WRITE_PROBLEM)//TRIM(ADJUSTL(MM_IDSTR)))
          WRITE(69,FMT=*)'%%MatrixMarket matrix coordinate ',TRIM(prec),
     &          " ",TRIM(symm)
          WRITE(69,*) id%N, id%N, id%NZ_loc
          IF (ASSOCIATED(id%A_loc)) THEN
            DO I=1,id%NZ
          WRITE(69,*) id%IRN_loc(I), id%JCN_loc(I), id%A_loc(I)  
            ENDDO
          ELSE
            DO I=1,id%NZ
              WRITE(69,*) id%IRN(I), id%JCN(I)
            ENDDO
          ENDIF
          CLOSE(69)
        ENDIF
      ENDIF
      IF ( id%MYID .eq. MASTER ) THEN
 1234   CONTINUE
        IF ( ( (KEEP(23) .NE. 0) .AND.
     *       ( (KEEP(23).NE.7) .OR. KEEP(50).EQ. 2 ) )
     *       .OR.
     *       (
     *       ASSOCIATED(id%A) .AND. KEEP(52) .EQ. 7 .AND.
     *       (KEEP(50).EQ.2)
     *       )
     *       .OR.
     *       KEEP(52) .EQ. -2
     *       ) THEN
          IF (.not.ASSOCIATED(id%A)) THEN
            IF (KEEP(23).GT.2) KEEP(23) = 1
          ENDIF
          CALL DMUMPS_203(id%N, id%NZ, KEEP(23), id%IS1(1), id,
     *        ICNTL, INFO)
          IF (INFO(1) .LT. 0) THEN
            KEEP(23) = 0
            GOTO 12
          END IF
        END IF
        IF (ICNTL(39).NE.1 .and. ICNTL(39).NE.2) THEN
          KEEP(106)=1
        ELSE
          KEEP(106)=ICNTL(39)
        ENDIF
        IF (KEEP(55) .EQ. 0) THEN
          IF ( IORD .eq. 1 ) THEN
             LIW = 2 * id%NZ + 3 * id%N + 2
          ELSE
            LIW = 2 * id%NZ + 3 * id%N + 2
          ENDIF
          IF (LIW.LT.3*id%N) LIW = 3*id%N
        ELSE
#if defined(metis)
         IF (KEEP(60).NE.0 .OR. IORD.EQ.5) THEN
#else
         IF (KEEP(60).NE.0) THEN
#endif
           LIW = id%N + id%N + 1
         ELSE
           LIW =  id%N + id%N + id%N+3 + id%N+1
         ENDIF
        ENDIF
        IF (LIW.LT.3*id%N) LIW = 3*id%N
        IF (KEEP(23) .NE. 0) THEN
          IKEEP = id%N + 1
        ELSE
          IKEEP = 1
        END IF
        NA      = IKEEP +     id%N
        NE      = IKEEP + 2 * id%N
        FILS    = IKEEP + 3 * id%N
        FRERE   = FILS  +     id%N
        PTRAR   = FRERE +     id%N
        IF (KEEP(55) .EQ. 0) THEN
          NFSIZ  = PTRAR + 4 * id%N
          MAXIS1_CHECK = NFSIZ + id%N - 1
        ELSE
          NFSIZ  = PTRAR + 2 * (NELT + 1)
          MAXIS1_CHECK = NFSIZ + id%N -1
        ENDIF
        IF ( id%MAXIS1 .LT. MAXIS1_CHECK ) THEN
          IF (LP.GE.0) THEN
          WRITE(LP,*) '***********************************'
          WRITE(LP,*) 'MAXIS1 and MAXIS1_CHECK are different !!'
          WRITE(LP,*) 'MAXIS1, MAXIS1_CHECK=',id%MAXIS1,
     *                 MAXIS1_CHECK
          WRITE(LP,*) 'This might cause problems ...'
          WRITE(LP,*) '***********************************'
          ENDIF
        END IF
        IF ( IORD .EQ. 1 ) THEN
          DO I = 1, id%N
            id%IS1( IKEEP + I - 1 ) = id%PERM_IN( I )
          END DO
        END IF
        CALL DMUMPS_633(KEEP(12),ICNTL(14),
     *       KEEP(50),KEEP(54),ICNTL(6),ICNTL(8))
#if defined(try_null_space)
        KEEP(53) = ICNTL(16)
        IF ( KEEP(53) .LT. 0 .OR.
     *       KEEP(53) .GT. 10 ) THEN
           KEEP(53) = 0
        END IF
        IF ( KEEP(53) > 2 .and. KEEP(53) < 9 ) THEN
          IF ( MPG > 0 ) THEN
            WRITE( MPG, '(A)' )
     *    ' Null space option currently disabled.'
            WRITE( MPG, '(A)' )
     *    ' Contact authors to experiment it.'
          END IF
          KEEP(53) = 0
        END IF
#else
        KEEP(53)=0
#endif
        IF(KEEP(86).EQ.1)THEN
           IF(KEEP(47).LT.2) KEEP(47)=2
        ENDIF
#ifndef FORCE5_
        IF(KEEP(48).EQ.5)THEN
#endif
           IF(KEEP(50).EQ.0)THEN
              KEEP(87)=50
              KEEP(88)=50
           ELSE
              KEEP(87)=70
              KEEP(88)=70
           ENDIF
#ifndef FORCE5_
        ENDIF
#endif
        IF((id%NSLAVES.EQ.1).AND.(id%KEEP(76).GT.3))THEN
           id%KEEP(76)=2
        ENDIF
        IF(KEEP(81).GT.0)THEN
           IF(KEEP(47).LT.2) KEEP(47)=2
        ENDIF
        INFOG(1) = 0
        INFOG(2) = 0
        INFOG(8) = -1
        IF ( .not. ASSOCIATED( id%LISTVAR_SCHUR ) ) THEN
          SIZE_SCHUR_PASSED = 1
          LISTVAR_SCHUR_2BE_FREED=.TRUE.
          ALLOCATE( id%LISTVAR_SCHUR( 1 ), stat=allocok )
          IF ( allocok .GT. 0 ) THEN
            WRITE(*,*) 'PB allocating an array of size 1 in Schur '
            CALL DMUMPS_ABORT()
          END IF
        ELSE
          SIZE_SCHUR_PASSED=id%SIZE_SCHUR
          LISTVAR_SCHUR_2BE_FREED=.FALSE.
        END IF
        IF (KEEP(55) .EQ. 0) THEN
          CALL DMUMPS_195(id%N, id%NZ, id%IRN(1), id%JCN(1),
     *        LIW, id%IS1(IKEEP),
     *        id%IS1(PTRAR), IORD, id%IS1(NFSIZ),
     *        id%IS1(FILS), id%IS1(FRERE),
     *        id%LISTVAR_SCHUR(1), SIZE_SCHUR_PASSED,
     *        ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),id%NSLAVES, 
     *        id%SYM,id%IS1(1),id)
          IF ( (KEEP(23).LE.-1).AND.(KEEP(23).GE.-6) ) THEN
            KEEP(23) = -KEEP(23)
            IF (.not.ASSOCIATED(id%A)) KEEP(23) = 1
            GOTO 1234
          ENDIF
          INFOG(7)=IORD
        ELSE
          ALLOCATE( IWtemp ( 3*id%N ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = 3*id%N
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'IWtemp'
              END IF
            GOTO 12
          ENDIF
          ALLOCATE( XNODEL ( id%N+1 ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = id%N + 1
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'XNODEL'
              END IF
            GOTO 12
          ENDIF
          IF (LELTVAR.ne.id%ELTPTR(NELT+1)-1)  THEN
           INFO(1) = -2002
           INFO(2) = id%ELTPTR(NELT+1)-1
           GOTO 12
          ENDIF
          ALLOCATE( NODEL ( LELTVAR ), stat = ierr )
          IF ( ierr .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = LELTVAR
              IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'NODEL'
              END IF
            GOTO 12
          ENDIF
          CALL DMUMPS_128(id%N, NELT,
     *        id%ELTPTR(1), id%ELTVAR(1), LIW,
     *        id%IS1(IKEEP),
     *        IWtemp, IORD, id%IS1(NFSIZ), id%IS1(FILS),
     *        id%IS1(FRERE), id%LISTVAR_SCHUR, SIZE_SCHUR_PASSED, 
     *        ICNTL(1), INFOG(1), KEEP(1),KEEP8(1),
     *        id%ELTPROC(1), id%NSLAVES, 
     *        XNODEL(1), NODEL(1))
          INFOG(7)=IORD
          IF ( KEEP(24).EQ.0 ) DEALLOCATE( IWtemp )
        ENDIF
       IF ( KEEP(52).EQ.7 ) THEN
          IF (id%SYM.EQ.0) THEN
             KEEP(52) = 4
          ELSE
             KEEP(52) = 0
          ENDIF
       ENDIF
        IF ( LISTVAR_SCHUR_2BE_FREED ) THEN
          DEALLOCATE( id%LISTVAR_SCHUR )
          NULLIFY   ( id%LISTVAR_SCHUR )
        ENDIF
        IF ( KEEP(60).NE.0) THEN
         I = id%SIZE_SCHUR
         IF (id%SYM.EQ.0) THEN
          IF (INFOG(20).GT.0) 
     &          INFOG(20) = INFOG(20) - I*I
          KEEP8(1)  = KEEP8(1) - I*I
         ELSE
          IF (INFOG(20).GT.0) 
     &          INFOG(20) = INFOG(20) - (I*(I+1))/2
          KEEP8(1)  = KEEP8(1) -  (I*(I+1))/2
         ENDIF
        END IF
        INFO(1)=INFOG(1)
        INFO(2)=INFOG(2)
        KEEP(28) = INFOG(6)
        IF ( INFO(1) .LT. 0 ) THEN
          GO TO 12
        ENDIF
        ALLOCATE(SSARBR(id%N), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          id%INFO(1) = -7
          id%INFO(2) = id%N
          IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'SSARBR'
          END IF
          GOTO 12
        END IF
        IF (KEEP(24).EQ.0) THEN
          id%NBSA = 0
          IF ( (id%KEEP(60).EQ.0).
     &            AND.(id%KEEP(53).EQ.0))  THEN 
               id%KEEP(20)=0
               id%KEEP(38)=0
          ENDIF
          id%KEEP(56)=0
          id%PROCNODE = 0
          IF (id%KEEP(60) .EQ. 2 .OR. id%KEEP(60).EQ.3) THEN
            CALL DMUMPS_564(id%KEEP(38), id%PROCNODE(1),
     *           1+2*id%NSLAVES, id%IS1(FILS),id%N)
          ENDIF
      ELSE
      PEAK=0.0D0
      PERM=0
         CALL DMUMPS_534(id%N, id%IS1(FRERE),
     *      id%IS1(FILS), id%IS1(NA),
     *      id%IS1(NE),id%IS1(NFSIZ),
     *      PERM,
     *      id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47),id%KEEP(81)
     *    ,id%KEEP(215)
     *    ,id%PROCNODE(1),id%NSLAVES, PEAK
     $     )
      IF (INFO(1) < 0 ) GOTO 12
      IF ( PROKG ) THEN
        IF (PERM == 0) THEN
        WRITE( MPG,'(A,F10.0) ')
     *' ** Peak of sequential stack size (number of real entries)   :',
     * PEAK
        ELSE
        WRITE( MPG,'(A,F10.0) ')
     *' ** Peak of sequential total memory (number of real entries) :',
     * PEAK
        ENDIF
       END IF
        call DMUMPS_537(id%N,id%NSLAVES,ICNTL,INFOG,
     *              id%IS1(NE),
     *              id%IS1(NFSIZ),
     *              id%IS1(FRERE),
     *              id%IS1(FILS),
     *              KEEP,KEEP8,id%PROCNODE(1),
     *              SSARBR(1),id%NBSA,PEAK,ierr)
          if(ierr.eq.-999) then 
              write(6,*) ' Internal error in DMUMPS_369'
             INFO(1) = ierr
             GOTO 12
          endif
          if(ierr.ne.0) then 
             INFO(1) = -135
             INFO(2) = ierr
             GOTO 12
          endif
        IF (KEEP(55) .EQ. 0) THEN
         CALL DMUMPS_202(id%N, id%NZ, id%IS1(IKEEP),
     *      id%IS1(FILS), id%IS1(FRERE), id%IS1(NFSIZ),
     *      id%IS1(IKEEP+2*id%N), 
     *      id%IS1(IKEEP+id%N),
     *      id%IRN(1), id%JCN(1), id%IS1(PTRAR+2*id%N) , 
     *      id%IS1(PTRAR+3*id%N), id%IS1(PTRAR),
     *      ICNTL(1), INFOG(1), KEEP(1), KEEP8(1))
        ELSE
         CALL DMUMPS_353(id%N, id%IS1(IKEEP),
     *      id%IS1(FILS), 
     *      id%IS1(FRERE), id%IS1(NFSIZ),
     *      id%IS1(IKEEP+2*id%N), 
     *      id%IS1(IKEEP+id%N), IWtemp,
     *      ICNTL(1), INFOG(1), KEEP(1),KEEP8(1))
          DEALLOCATE( IWtemp )
        ENDIF
      ENDIF
       DEALLOCATE( SSARBR)
        IF ( KEEP(54) .eq. 3 ) THEN
          DEALLOCATE( id%IRN )
          DEALLOCATE( id%JCN )
        END IF
       IF (KEEP(55).NE.0) THEN
          IF (ASSOCIATED(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
          IF (ASSOCIATED(id%FRTELT)) DEALLOCATE(id%FRTELT)
          ALLOCATE( id%FRTPTR(id%N+1), id%FRTELT(id%NELT),
     *              stat=allocok )
          IF ( allocok .ne. 0) THEN
            INFO(1) = -7
            INFO(2) = id%NELT+id%N+1
            IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'FRTPTR'
            END IF
            GOTO 12
          ENDIF
         CALL DMUMPS_153(
     *   id%N, NELT, id%ELTPTR(NELT+1)-1, id%IS1(FRERE),
     *   id%IS1(FILS),
     *   id%IS1(IKEEP+id%N), id%IS1(IKEEP+2*id%N), XNODEL, NODEL,
     *   id%FRTPTR(1), id%FRTELT(1), id%ELTPROC(1))
         DEALLOCATE(XNODEL)
         DEALLOCATE(NODEL)
       ENDIF
        IF ( INFO( 1 ) .LT. 0 ) GOTO 12
        IF ( KEEP(55) .ne. 0 ) THEN
           CALL DMUMPS_120(id%N, NELT, id%ELTPROC(1),id%NSLAVES,
     *                        id%PROCNODE)
        END IF
        NB_NIV2 = KEEP(56)
        IF ( NB_NIV2.GT.0 ) THEN
          ALLOCATE(PAR2_NODES(NB_NIV2),
     *             stat=allocok)
          if (allocok .gt.0) then
             info(1)= -7
             info(2)= NB_NIV2
            IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'PAR2_NODES'
            END IF
             goto 12
          end if
        ENDIF
        IF ((NB_NIV2.GT.0) .AND. (KEEP(24).EQ.0)) THEN
          INIV2 = 0
          DO 777 INODE = 1, id%N
            IF ( ( id%IS1(FRERE+INODE-1) .NE. id%N+1 ) .AND.
     *           ( DMUMPS_330(INODE,id%PROCNODE(1),id%NSLAVES)
     *             .eq. 2) ) THEN
              INIV2 = INIV2 + 1
              PAR2_NODES(INIV2) = INODE
            END IF
 777      CONTINUE
          IF ( INIV2 .NE. NB_NIV2 ) THEN
            WRITE(*,*) "Internal Error 2 in DMUMPS_26",
     *                  INIV2, NB_NIV2
            CALL DMUMPS_ABORT()
          ENDIF
        ENDIF
        IF ( (KEEP(24) .NE. 0) .AND. (NB_NIV2.GT.0) ) THEN
          IF ( ASSOCIATED(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
          ALLOCATE( id%CANDIDATES(id%nslaves+1,NB_NIV2),
     *             stat=allocok)
          if (allocok .gt.0) then
             info(1)= -7
             info(2)= NB_NIV2*(id%nslaves+1)
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'CANDIDATES'
             END IF
             goto 12
          end if
          call DMUMPS_393
     *         (PAR2_NODES,id%CANDIDATES,ierr)
          if(ierr.ne.0)  then
             INFO(1) = -2002
             GOTO 12
          endif
          call DMUMPS_494()
          if(ierr.ne.0)  then
             INFO(1) = -2002
             GOTO 12
          endif
        ELSE
          ALLOCATE(id%CANDIDATES(1,1), stat=allocok)
          IF (allocok .NE. 0) THEN
             info(1)= -7
             info(2)= 1
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'CANDIDATES'
             END IF
             goto 12
          ENDIF
        ENDIF
        KEEP(129)=INFO(7)
        KEEP(130)=INFO(8)
 12     CONTINUE
      END IF
      CALL DMUMPS_276( ICNTL, INFO, id%COMM, id%MYID )
      IF ( INFO(1) < 0 ) RETURN
      CALL MPI_BCAST( id%KEEP, 110, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      CALL MPI_BCAST( id%KEEP(205), 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      IF (id%MYID==MASTER) KEEP(127)=INFOG(5)
      CALL MPI_BCAST( id%KEEP(127), 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      CALL MPI_BCAST( id%N, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      IF ( id%KEEP(55) .EQ. 0) THEN
        CALL MPI_BCAST( id%NZ, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      ELSE
        CALL MPI_BCAST( id%NA_ELT, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      ENDIF
      CALL MPI_BCAST( id%NBSA, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      IF ( id%KEEP(55) .NE. 0) THEN
        CALL MPI_BCAST( id%NELT, 1, MPI_INTEGER, MASTER,
     *                                   id%COMM, IERR )
      ENDIF
      IF ( ASSOCIATED( id%STEP ) ) DEALLOCATE(id%STEP)
      ALLOCATE  (id%STEP(id%N),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
           WRITE(LP, 150) 'id%STEP'
        END IF
        INFO(1) = -7
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%PROCNODE_STEPS ))
     *    DEALLOCATE(id%PROCNODE_STEPS)
      ALLOCATE(id%PROCNODE_STEPS( id%KEEP(28)), stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%PROCNODE_STEPS'
        END IF
        INFO(1) = -7
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%NE_STEPS ) ) DEALLOCATE(id%NE_STEPS)
      ALLOCATE (id%NE_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%NE_STEPS'
        END IF
        INFO(1) = -7
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%ND_STEPS ) ) DEALLOCATE(id%ND_STEPS)
      ALLOCATE (id%ND_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%ND_STEPS'
        END IF
        INFO(1) = -7
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%FRERE_STEPS ) ) DEALLOCATE(id%FRERE_STEPS)
      ALLOCATE (id%FRERE_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%FRERE_STEPS'
        END IF
        INFO(1) = -7
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%DAD_STEPS ) ) DEALLOCATE(id%DAD_STEPS)
      ALLOCATE (id%DAD_STEPS(id%KEEP(28)),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%DAD_STEPS'
        END IF
        INFO(1) = -7
        INFO(2) = id%KEEP(28)
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%SYM_PERM ) ) DEALLOCATE(id%SYM_PERM)
      ALLOCATE(id%SYM_PERM(id%N),stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%SYM_PERM'
        END IF
        INFO(1) = -7
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%UNS_PERM ) ) DEALLOCATE(id%UNS_PERM)
      IF ( ASSOCIATED( id%FILS ) ) DEALLOCATE(id%FILS)
      ALLOCATE(id%FILS( id%N ), stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%FILS'
        END IF
        INFO(1) = -7
        INFO(2) = id%N
        GOTO 94
      ENDIF
      IF ( ASSOCIATED( id%PTRAR ) ) DEALLOCATE(id%PTRAR)
      IF ( id%KEEP(55).eq.0 ) THEN
        LPTRAR = id%N+id%N
      ELSE
        LPTRAR = id%NELT+id%NELT+2
      ENDIF
      ALLOCATE(id%PTRAR( LPTRAR ), stat=allocok)
      IF ( allocok .ne. 0) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%PTRAR'
        END IF
        INFO(1) = -7
        INFO(2) = LPTRAR
        GOTO 94
      ENDIF
      IF ( id%KEEP(55) /= 0 ) THEN
        IF ( id%MYID .NE. MASTER ) THEN
          IF (ASSOCIATED(id%FRTPTR)) DEALLOCATE(id%FRTPTR)
          IF (ASSOCIATED(id%FRTELT)) DEALLOCATE(id%FRTELT)
          ALLOCATE( id%FRTPTR(id%N+1), id%FRTELT(id%NELT),
     *              stat=allocok )
          IF ( allocok .ne. 0) THEN
            INFO(1) = -7
            INFO(2) = id%NELT+id%N+1
            IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'id%FRTPTR/FRTELT'
            END IF
            GOTO 94
          ENDIF
        END IF
      ELSE
        ALLOCATE( id%FRTPTR(1), id%FRTELT(1),
     *            stat=allocok )
        IF ( allocok .ne. 0) THEN
          INFO(1) = -7
          INFO(2) = 2
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'id%FRTPTR/FRTELT'
          END IF
          GOTO 94
        ENDIF
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
        DO I=1,id%N
          id%SYM_PERM(I) = id%IS1(IKEEP+I-1)
        END DO
      END IF
      IF ( id%MYID == MASTER .AND. id%KEEP(23) .NE. 0 ) THEN
        ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
        IF ( allocok .ne. 0) THEN
          INFO(1) = -7
          INFO(2) = id%N
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'id%UNS_PERM'
          END IF
          GOTO 94
        ENDIF
        DO I=1,id%N
          id%UNS_PERM(I) = id%IS1(I)
        END DO
      ENDIF
  94  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
      IF ( id%MYID .EQ. MASTER ) THEN
        DO I=1,id%N
          id%FILS(I) = id%IS1(FILS+I-1)
        ENDDO
        IF ( id%KEEP(55) == 0 ) THEN
          DO I= 1, LPTRAR
            id%PTRAR(I) = id%IS1(PTRAR+I-1)
          END DO
        ELSE
          DO I=1, id%NELT+1
            id%PTRAR(id%NELT+I+1)=id%ELTPTR(I)
          ENDDO
        ENDIF
      END IF
      IF (id%MYID .EQ. MASTER ) THEN
        IF (id%N.eq.1) THEN
          NBROOT = 1
          NBLEAF = 1
        ELSE IF (id%IS1(NA+id%N-1) .LT.0) THEN
          NBLEAF = id%N
          NBROOT = id%N
        ELSE IF (id%IS1(NA+id%N-2) .LT.0) THEN
          NBLEAF = id%N-1
          NBROOT = id%IS1(NA+id%N-1)
        ELSE
          NBLEAF = id%IS1(NA+id%N-2)
          NBROOT = id%IS1(NA+id%N-1)
        ENDIF
        id%LNA = 2+NBLEAF+NBROOT
      ENDIF
      CALL MPI_BCAST( id%LNA, 1, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF ( ASSOCIATED( id%NA ) ) DEALLOCATE(id%NA)
      ALLOCATE(id%NA(id%LNA),stat=allocok)
      IF ( allocok .ne. 0) THEN
        INFO(1) = -7
        INFO(2) = LPTRAR
        IF ( LP .GT. 0 ) THEN
          WRITE(LP, 150) 'id%NA'
        END IF
        GOTO 96
      ENDIF
      IF (id%MYID .EQ.MASTER ) THEN
        id%NA(1) = NBLEAF
        id%NA(2) = NBROOT
        LEAF = 3
        IF ( id%N == 1 ) THEN
          id%NA(LEAF) = 1
          LEAF = LEAF + 1
        ELSE IF (id%IS1(NA+id%N-1) < 0) THEN
          id%NA(LEAF) = - id%IS1(NA+id%N-1)-1
          LEAF = LEAF + 1
          DO I = 1, NBLEAF - 1
            id%NA(LEAF) = id%IS1(NA+I-1)
            LEAF = LEAF + 1
          ENDDO
        ELSE IF (id%IS1(NA+id%N-2) < 0 ) THEN
          INODE = - id%IS1(NA+id%N-2) - 1
          id%NA(LEAF) = INODE
          LEAF =LEAF + 1
          IF ( NBLEAF > 1 ) THEN
            DO I = 1, NBLEAF - 1
              id%NA(LEAF) = id%IS1(NA+I-1)
              LEAF = LEAF + 1
            ENDDO
          ENDIF
        ELSE
          DO I = 1, NBLEAF
            id%NA(LEAF) = id%IS1(NA+I-1)
            LEAF = LEAF + 1
          ENDDO
        END IF
      END IF
 96   CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( id%MYID .EQ. MASTER ) THEN
        ISTEP = 0
        DO I = 1, id%N
          IF ( id%IS1(FRERE+I-1) .ne. id%N + 1 ) THEN
            ISTEP = ISTEP + 1
            id%STEP(I)=ISTEP
            INN = id%IS1(FILS+I-1)
            DO WHILE ( INN .GT. 0 )
              id%STEP(INN) = - ISTEP
              INN = id%IS1(FILS + INN -1)
            END DO
            IF (id%IS1(FRERE+I-1) .eq. 0) THEN
              id%NA(LEAF) = I
              LEAF = LEAF + 1
            ENDIF
          ENDIF
        END DO
        IF ( LEAF - 1 .NE. 2+NBROOT + NBLEAF ) THEN
          WRITE(*,*) 'Internal error 2 in DMUMPS_26'
          CALL DMUMPS_ABORT()
        ENDIF
        IF ( ISTEP .NE. id%KEEP(28) ) THEN
          write(*,*) 'Internal error 3 in DMUMPS_26'
          CALL DMUMPS_ABORT()
        ENDIF
        DO I = 1, id%N
          IF (id%IS1(FRERE+I-1) .NE. id%N+1) THEN
            id%PROCNODE_STEPS(id%STEP(I)) = id%PROCNODE( I )
            id%FRERE_STEPS(id%STEP(I))    = id%IS1(FRERE+I-1)
            id%NE_STEPS(id%STEP(I))    = id%IS1(NE+I-1)
            id%ND_STEPS(id%STEP(I))    = id%IS1(NFSIZ+I-1)
          ENDIF
        ENDDO
        DO I = 1, id%N
          IF ( id%STEP(I) .LE. 0) CYCLE
          IF (id%IS1(FRERE+I-1) .eq. 0) THEN
            id%DAD_STEPS(id%STEP(I)) = 0
          ENDIF
          IFS = id%IS1(FILS+I-1)
          DO WHILE ( IFS .GT. 0 )
              IFS= id%IS1(FILS + IFS -1)
          END DO
          IFS = -IFS
          DO WHILE (IFS.GT.0) 
              id%DAD_STEPS(id%STEP(IFS)) = I
              IFS   = id%IS1(FRERE+IFS-1)
          ENDDO
        END DO
        DEALLOCATE(id%PROCNODE)
        NULLIFY(id%PROCNODE)
        DEALLOCATE(id%IS1)
        NULLIFY(id%IS1)
       IF (( id%KEEP(70) .NE. -1 ).OR.(id%KEEP(47).EQ.4).OR.
     *     (( id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2)))
     *      THEN
          IS_REORDERTREE_CALLED=.TRUE.
         IF ((id%KEEP(47) .EQ. 4).OR.
     *       (( id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
           IF(id%NSLAVES.GT.1) THEN
              SIZE_TEMP_MEM = id%NBSA
           ELSE
              SIZE_TEMP_MEM = id%NA(2)
           ENDIF
         ELSE
           SIZE_TEMP_MEM = 1
         ENDIF
         IF(id%KEEP(76).EQ.4)THEN
            SIZE_DEPTH_FIRST=id%KEEP(28)
         ELSE
            SIZE_DEPTH_FIRST=1
         ENDIF
         ALLOCATE(TEMP_MEM(SIZE_TEMP_MEM,id%NSLAVES),stat=allocok) 
         IF (allocok .ne.0) then
           info(1)= -7
           info(2)= SIZE_TEMP_MEM*id%NSLAVES
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'TEMP_MEM'
           END IF
           GOTO 80
         end if
         ALLOCATE(TEMP_LEAF(SIZE_TEMP_MEM,id%NSLAVES),stat=allocok) 
         IF (allocok .ne.0) then
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'TEMP_LEAF'
           END IF
           info(1)= -7
           info(2)= SIZE_TEMP_MEM*id%NSLAVES
           GOTO 80
         end if
         ALLOCATE(TEMP_SIZE(SIZE_TEMP_MEM,id%NSLAVES),stat=allocok) 
         IF (allocok .ne.0) then
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'TEMP_SIZE'
           END IF
           info(1)= -7
           info(2)= SIZE_TEMP_MEM*id%NSLAVES
           GOTO 80
         end if
         ALLOCATE(TEMP_ROOT(SIZE_TEMP_MEM,id%NSLAVES),stat=allocok) 
         IF (allocok .ne.0) then
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'TEMP_ROOT'
           END IF
           info(1)= -7
           info(2)= SIZE_TEMP_MEM*id%NSLAVES
           GOTO 80
         end if
         ALLOCATE(DEPTH_FIRST(SIZE_DEPTH_FIRST),stat=allocok) 
         IF (allocok .ne.0) then
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'DEPTH_FIRST'
           END IF
           info(1)= -7
           info(2)= SIZE_DEPTH_FIRST
           GOTO 80
         end if
         IF(id%KEEP(76).EQ.5)THEN
            SIZE_COST_TRAV=id%KEEP(28)
         ELSE
            SIZE_COST_TRAV=1
         ENDIF
         ALLOCATE(COST_TRAV_TMP(SIZE_COST_TRAV),stat=allocok) 
         IF (allocok .ne.0) then
           IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'COST_TRAV_TMP'
           END IF
           info(1)= -7
           info(2)= SIZE_COST_TRAV
           GOTO 80
         end if
         IF(id%KEEP(76).EQ.5)THEN
            IF(id%KEEP(70).EQ.0)THEN
               id%KEEP(70)=5
            ENDIF
            IF(id%KEEP(70).EQ.1)THEN
               id%KEEP(70)=6
            ENDIF
         ENDIF
         IF(id%KEEP(76).EQ.4)THEN
            IF(id%KEEP(70).EQ.0)THEN
               id%KEEP(70)=3
            ENDIF
            IF(id%KEEP(70).EQ.1)THEN
               id%KEEP(70)=4
            ENDIF
         ENDIF
         CALL DMUMPS_363(id%N, id%FRERE_STEPS(1), id %STEP(1),
     *      id%FILS(1), id%NA(1), id%LNA, id%NE_STEPS(1),
     *      id%ND_STEPS(1), id%DAD_STEPS, id%KEEP(28), .TRUE.,
     *      id%KEEP(28), id%KEEP(70),
     *      id%KEEP(50), id%INFO(1), id%ICNTL(1),id%KEEP(47),id%KEEP(81)
     *      ,id%KEEP(76),id%KEEP(215),
     *      id%PROCNODE_STEPS(1),TEMP_MEM,id%NSLAVES, 
     *      SIZE_TEMP_MEM, PEAK,id%KEEP(90),SIZE_DEPTH_FIRST,
     $      SIZE_COST_TRAV,DEPTH_FIRST(1),COST_TRAV_TMP(1),
     $      TEMP_LEAF,TEMP_SIZE,TEMP_ROOT
     $        )
       END IF
 80     CONTINUE
      ENDIF
      CALL DMUMPS_276( ICNTL, INFO,
     *                      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      CALL MPI_BCAST( id%SYM_PERM(1), id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FILS(1), id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF ( id%KEEP(55) .eq. 0 ) THEN
        CALL MPI_BCAST( id%PTRAR(1), LPTRAR, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
      ELSE
        CALL MPI_BCAST( id%PTRAR(id%NELT+2), id%NELT+1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
      END IF
      CALL MPI_BCAST( id%NA(1), id%LNA, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%STEP, id%N, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%PROCNODE_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%DAD_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%FRERE_STEPS, id% KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR)
      CALL MPI_BCAST( id%NE_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      CALL MPI_BCAST( id%ND_STEPS, id%KEEP(28), MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      IF(id%KEEP(76).EQ.4)THEN
         IF(ASSOCIATED(id%DEPTH_FIRST))
     *        DEALLOCATE(id%DEPTH_FIRST)
         ALLOCATE(id%DEPTH_FIRST(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            info(1)= -7
            info(2)= id%KEEP(28)
            IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%DEPTH_FIRST(1:id%KEEP(28))=DEPTH_FIRST(1:id%KEEP(28))
         ENDIF
         CALL MPI_BCAST( id%DEPTH_FIRST, id%KEEP(28), MPI_INTEGER,
     *        MASTER, id%COMM, IERR )         
      ELSE
         IF(ASSOCIATED(id%DEPTH_FIRST))
     *        DEALLOCATE(id%DEPTH_FIRST)
         ALLOCATE(id%DEPTH_FIRST(1),stat=allocok)
         IF (allocok .ne.0) then
            info(1)= -7
            info(2)= 1
            IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'id%DEPTH_FIRST'
            END IF
            GOTO 87
         END IF
         id%DEPTH_FIRST(1)=0
      ENDIF
      IF(id%KEEP(76).EQ.5)THEN
         IF(ASSOCIATED(id%COST_TRAV))
     *        DEALLOCATE(id%COST_TRAV)
         ALLOCATE(id%COST_TRAV(id%KEEP(28)),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'id%COST_TRAV'
            END IF
            info(1)= -7
            info(2)= id%KEEP(28)
            GOTO 87
         END IF
         IF(id%MYID.EQ.MASTER)THEN
            id%COST_TRAV(1:id%KEEP(28))=COST_TRAV_TMP(1:id%KEEP(28))
         ENDIF
         CALL MPI_BCAST( id%COST_TRAV, id%KEEP(28),
     $        MPI_DOUBLE_PRECISION,MASTER, id%COMM, IERR )         
      ELSE
         IF(ASSOCIATED(id%COST_TRAV))
     *        DEALLOCATE(id%COST_TRAV)
         ALLOCATE(id%COST_TRAV(1),stat=allocok)
         IF (allocok .ne.0) then
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'id%COST_TRAV(1)'
             END IF
            info(1)= -7
            info(2)= 1
            GOTO 87
         END IF
         id%COST_TRAV(1)=0
      ENDIF
      IF (id%KEEP(47) .EQ. 4 .OR.
     *    ((id%KEEP(81) .GT. 0).AND.(id%KEEP(47).GE.2))) THEN
       IF(id%MYID .EQ. MASTER)THEN
         DO K=1,id%NSLAVES
           DO J=1,SIZE_TEMP_MEM
             IF(TEMP_MEM(J,K).EQ.-1) GOTO 666 
           ENDDO
 666       CONTINUE
           J=J-1
           IF (id%KEEP(46) == 1) THEN
             IDEST = K - 1
           ELSE
             IDEST = K
           ENDIF
           IF (IDEST .NE. MASTER) THEN
             CALL MPI_SEND(J,1,MPI_INTEGER,IDEST,0,
     *         id%COMM,IERR)
             CALL MPI_SEND(TEMP_MEM(1,K),J,MPI_INTEGER,
     *         IDEST, 0, id%COMM,IERR)
             CALL MPI_SEND(TEMP_LEAF(1,K),J,MPI_INTEGER,
     *         IDEST, 0, id%COMM,IERR)
             CALL MPI_SEND(TEMP_SIZE(1,K),J,MPI_INTEGER,
     *         IDEST, 0, id%COMM,IERR)
             CALL MPI_SEND(TEMP_ROOT(1,K),J,MPI_INTEGER,
     *         IDEST, 0, id%COMM,IERR)             
           ELSE
             IF(ASSOCIATED(id%MEM_SUBTREE))
     *         DEALLOCATE(id%MEM_SUBTREE)
             ALLOCATE(id%MEM_SUBTREE(J),stat=allocok)
             IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MEM_SUBTREE'
               END IF
               info(1)= -7
               info(2)= J
               GOTO 87
             END IF
             id%NBSA_LOCAL = J
             id%MEM_SUBTREE(1:J)=TEMP_MEM(1:J,1)
             IF(ASSOCIATED(id%MY_ROOT_SBTR))
     *         DEALLOCATE(id%MY_ROOT_SBTR)
             ALLOCATE(id%MY_ROOT_SBTR(J),stat=allocok)
             IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_ROOT_SBTR'
               END IF
               info(1)= -7
               info(2)= J
               GOTO 87
             END IF
             id%MY_ROOT_SBTR(1:J)=TEMP_ROOT(1:J,1)
             IF(ASSOCIATED(id%MY_FIRST_LEAF))
     *         DEALLOCATE(id%MY_FIRST_LEAF)
             ALLOCATE(id%MY_FIRST_LEAF(J),stat=allocok)
             IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_FIRST_LEAF'
               END IF
               info(1)= -7
               info(2)= J
               GOTO 87
             END IF
             id%MY_FIRST_LEAF(1:J)=TEMP_LEAF(1:J,1)
             IF(ASSOCIATED(id%MY_NB_LEAF))
     *         DEALLOCATE(id%MY_NB_LEAF)
             ALLOCATE(id%MY_NB_LEAF(J),stat=allocok)
             IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_NB_LEAF'
               END IF
               info(1)= -7
               info(2)= J
               GOTO 87
             END IF
             id%MY_NB_LEAF(1:J)=TEMP_SIZE(1:J,1)
           ENDIF
         ENDDO
       ELSE
         CALL MPI_RECV(id%NBSA_LOCAL,1,MPI_INTEGER,
     *                 MASTER,0,id%COMM,STATUS, IERR)
         IF(ASSOCIATED(id%MEM_SUBTREE))
     *     DEALLOCATE(id%MEM_SUBTREE)
         ALLOCATE(id%MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
         IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MEM_SUBTREE'
               END IF
               info(1)= -7
               info(2)= id%NBSA_LOCAL
               GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_ROOT_SBTR))
     *        DEALLOCATE(id%MY_ROOT_SBTR)
         ALLOCATE(id%MY_ROOT_SBTR(id%NBSA_LOCAL),stat=allocok)
         IF (allocok .ne.0) then
            IF ( LP .GT. 0 ) THEN
               WRITE(LP, 150) 'id%MY_ROOT_SBTR'
            END IF
            info(1)= -7
            info(2)= J
            GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_FIRST_LEAF))
     *        DEALLOCATE(id%MY_FIRST_LEAF)
         ALLOCATE(id%MY_FIRST_LEAF(id%NBSA_LOCAL),stat=allocok)
         IF (allocok .ne.0) then
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'MY_FIRST_LEAF'
             END IF
            info(1)= -7
            info(2)= J
            GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_NB_LEAF))
     *        DEALLOCATE(id%MY_NB_LEAF)
         ALLOCATE(id%MY_NB_LEAF(id%NBSA_LOCAL),stat=allocok)
         IF (allocok .ne.0) then
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'MY_NB_LEAF'
             END IF
            info(1)= -7
            info(2)= J
            GOTO 87
         END IF
         CALL MPI_RECV(id%MEM_SUBTREE(1),id%NBSA_LOCAL,
     *       MPI_INTEGER,MASTER,0,
     *       id%COMM,STATUS,IERR)
         CALL MPI_RECV(id%MY_FIRST_LEAF(1),id%NBSA_LOCAL,
     *       MPI_INTEGER,MASTER,0,
     *       id%COMM,STATUS,IERR)
         CALL MPI_RECV(id%MY_NB_LEAF(1),id%NBSA_LOCAL,
     *       MPI_INTEGER,MASTER,0,
     *       id%COMM,STATUS,IERR)
         CALL MPI_RECV(id%MY_ROOT_SBTR(1),id%NBSA_LOCAL,
     *       MPI_INTEGER,MASTER,0,
     *       id%COMM,STATUS,IERR)
       ENDIF
      ELSE
         id%NBSA_LOCAL = -999999
         IF(ASSOCIATED(id%MEM_SUBTREE))
     *     DEALLOCATE(id%MEM_SUBTREE)
         ALLOCATE(id%MEM_SUBTREE(1),stat=allocok)
         IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MEM_SUBTREE(1)'
               END IF
               info(1)= -7
               info(2)= 1
               GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_ROOT_SBTR))
     *        DEALLOCATE(id%MY_ROOT_SBTR)
         ALLOCATE(id%MY_ROOT_SBTR(1),stat=allocok)
         IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_ROOT_SBTR(1)'
               END IF
            info(1)= -7
            info(2)= 1
            GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_FIRST_LEAF))
     *        DEALLOCATE(id%MY_FIRST_LEAF)
         ALLOCATE(id%MY_FIRST_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_FIRST_LEAF(1)'
               END IF
            info(1)= -7
            info(2)= 1
            GOTO 87
         END IF
         IF(ASSOCIATED(id%MY_NB_LEAF))
     *        DEALLOCATE(id%MY_NB_LEAF)
         ALLOCATE(id%MY_NB_LEAF(1),stat=allocok)
         IF (allocok .ne.0) then
               IF ( LP .GT. 0 ) THEN
                 WRITE(LP, 150) 'id%MY_NB_LEAF(1)'
               END IF
            info(1)= -7
            info(2)= 1
            GOTO 87
         END IF
      ENDIF
      IF(id%MYID.EQ.MASTER)THEN
         IF(IS_REORDERTREE_CALLED)THEN 
            DEALLOCATE(TEMP_MEM)
            DEALLOCATE(TEMP_SIZE)
            DEALLOCATE(TEMP_ROOT)
            DEALLOCATE(TEMP_LEAF)
            DEALLOCATE(COST_TRAV_TMP)
            DEALLOCATE(DEPTH_FIRST)
         ENDIF
      ENDIF
87    CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                      id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( id%KEEP(55) /= 0 ) THEN
        CALL MPI_BCAST( id%FRTPTR(1), id%N+1, MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
        CALL MPI_BCAST( id%FRTELT(1), id%NELT,  MPI_INTEGER,
     *               MASTER, id%COMM, IERR )
      END IF
      NB_NIV2 = KEEP(56) ! KEEP(1:100) was broadcast in DMUMPS_BDCAST
      IF (  NB_NIV2.GT.0  ) THEN
        if (id%myid.ne.master) then
        IF (ASSOCIATED(id%CANDIDATES)) DEALLOCATE(id%CANDIDATES)
          allocate(PAR2_NODES(NB_NIV2),
     *             id%CANDIDATES(id%NSLAVES+1,NB_NIV2),
     *             STAT=allocok)
          IF (allocok .ne.0) then
             info(1)= -7
             info(2)= NB_NIV2*(id%NSLAVES+1)
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'PAR2_NODES/id%CANDIDATES'
             END IF
          end if
        end if
        CALL DMUMPS_276( ICNTL, INFO,
     *                       id%COMM, id%MYID )
        IF ( INFO(1).LT.0 ) RETURN
        CALL MPI_BCAST(PAR2_NODES(1),NB_NIV2,
     *         MPI_INTEGER, MASTER, id%COMM, IERR )
        IF (KEEP(24) .NE.0 ) THEN
          CALL MPI_BCAST(id%CANDIDATES(1,1),(NB_NIV2*(id%NSLAVES+1)),
     *         MPI_INTEGER, MASTER, id%COMM, IERR )
        ENDIF
      ENDIF
      IF ( ASSOCIATED(id%ISTEP_TO_INIV2)) THEN
         DEALLOCATE(id%ISTEP_TO_INIV2)
         NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF ( ASSOCIATED(id%I_AM_CAND)) THEN
        DEALLOCATE(id%I_AM_CAND)
        NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (NB_NIV2.EQ.0) THEN 
         id%KEEP(71) = 1
      ELSE
         id%KEEP(71) = id%KEEP(28)
      ENDIF
      ALLOCATE(id%ISTEP_TO_INIV2(id%KEEP(71)),
     *         id%I_AM_CAND(MAX(NB_NIV2,1)),
     *         stat=allocok)
      IF (allocok .gt.0) THEN
           IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
              WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
           END IF
           info(1)= -7
           IF (NB_NIV2.EQ.0) THEN
            info(2)= 2
           ELSE
            info(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
           END IF
           GOTO 321
      ENDIF
      IF ( NB_NIV2 .GT.0 ) THEN
        DO INIV2 = 1, NB_NIV2
         INN = PAR2_NODES(INIV2)
         id%ISTEP_TO_INIV2(id%STEP(INN)) = INIV2
        END DO 
         DO INIV2=1, NB_NIV2
          id%I_AM_CAND(INIV2)=.FALSE.
          DO I=1, id%CANDIDATES(id%NSLAVES+1,INIV2)
           IF (id%CANDIDATES(I,INIV2).EQ.id%MYID_NODES) THEN
             id%I_AM_CAND(INIV2)=.TRUE.
             EXIT
           ENDIF
          ENDDO
         END DO
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (ASSOCIATED(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
      ENDIF
      ALLOCATE(id%FUTURE_NIV2(id%NSLAVES), stat=allocok)
      IF (allocok .gt.0) THEN
           IF ( LP .GT. 0 ) THEN
              WRITE(LP, 150) 'FUTURE_NIV2'
           END IF
           info(1)= -7
           info(2)= id%NSLAVES
           GOTO 321
      ENDIF
      id%FUTURE_NIV2=0
      DO INIV2 = 1, NB_NIV2
        IDEST = DMUMPS_275(id%STEP(PAR2_NODES(INIV2)),
     *                       id%PROCNODE_STEPS,
     *                       id%NSLAVES)
        id%FUTURE_NIV2(IDEST+1)=id%FUTURE_NIV2(IDEST+1)+1
      ENDDO
#endif
      IF ( I_AM_SLAVE ) THEN
        IF ( ASSOCIATED(id%TAB_POS_IN_PERE)) THEN
          DEALLOCATE(id%TAB_POS_IN_PERE)
          NULLIFY(id%TAB_POS_IN_PERE)
        ENDIF
        ALLOCATE(id%TAB_POS_IN_PERE(id%NSLAVES+2,MAX(NB_NIV2,1)),
     *      stat=allocok)
        IF (allocok .gt.0) THEN
             IF ( LP .GT. 0 ) THEN
                WRITE(LP, 150) 'id%ISTEP_TO_INIV2'
                WRITE(LP, 150) 'id%TAB_POS_IN_PERE'
             END IF
             info(1)= -7
             IF (NB_NIV2.EQ.0) THEN
              info(2)= 2
             ELSE
              info(2)= id%KEEP(28)+NB_NIV2*(id%NSLAVES+2)
             END IF
             GOTO 321
        ENDIF
      END IF
      IF (NB_NIV2.GT.0) DEALLOCATE (PAR2_NODES)
  321 CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *     id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( KEEP(23).NE.0 .and. id%MYID .EQ. MASTER ) THEN
        IKEEP = id%N + 1
      ELSE
        IKEEP = 1
      END IF
      FILS   = IKEEP + 3 * id%N
        NE   = IKEEP + 2 * id%N
        NA   = IKEEP +     id%N
      FRERE  = FILS  + id%N
      PTRAR  = FRERE + id%N
      IF (KEEP(55) .EQ. 0) THEN
        IF ( id%MYID.EQ.MASTER ) THEN
          NFSIZ   = PTRAR  + 4 * id%N
        ELSE
          NFSIZ   = PTRAR  + 2 * id%N
        ENDIF
      ELSE
        NFSIZ   = PTRAR  + 2 * (NELT + 1)
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
          CALL DMUMPS_164( id%MYID,
     *    id%NSLAVES, id%N, id%root,
     *    id%COMM_NODES, KEEP( 38 ), id%FILS(1),
     *    id%KEEP(50), id%KEEP(46),
     *    id%KEEP(51)
     *     , id%KEEP(60), id%NPROW, id%NPCOL, id%MBLOCK, id%NBLOCK
     *    )
      ELSE
          id%root%yes = .FALSE.
      END IF
      IF ( KEEP(38) .NE. 0 .and. I_AM_SLAVE ) THEN
          CALL MPI_ALLREDUCE(id%root%MYROW, MYROW_CHECK, 1,
     *         MPI_INTEGER, MPI_MAX, id%COMM_NODES, IERR)
          IF ( MYROW_CHECK .eq. -1) THEN
            INFO(1) = -25
            INFO(2) = 0
          END IF
          IF ( id%root%MYROW .LT. -1 .OR.
     *         id%root%MYCOL .LT. -1 ) THEN
            INFO(1) = -25
            INFO(2) = 0
          END IF
          IF ( LP > 0 .AND. INFO(1) == -25 ) THEN
            WRITE(LP, '(A)') 'Problem with your version of the BLACS.'
            WRITE(LP, '(A)') 'Try using a BLACS version from netlib.'
          ENDIF
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) RETURN
      IF ( I_AM_SLAVE ) THEN
        IF (KEEP(55) .EQ. 0) THEN
          CALL DMUMPS_24( id%MYID,
     *        id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     *        id%STEP(1), id%PTRAR(1),
     *        id%PTRAR(id%N +1),
     *        id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     *        id%KEEP(1),id%KEEP8(1), ICNTL, id )
        ELSE
          CALL DMUMPS_25( id%MYID,
     *        id%NSLAVES, id%N, id%PROCNODE_STEPS(1),
     *        id%STEP(1),
     *        id%PTRAR(1),
     *        id%PTRAR(id%NELT+2 ),
     *        id%NELT, 
     *        id%FRTPTR(1), id%FRTELT(1),
     *        id%KEEP(1), id%KEEP8(1), ICNTL, id%SYM )
        ENDIF
        IF ( id%root%yes ) THEN
          LOCAL_M = NUMROC( id%ND_STEPS(id%STEP(KEEP(38))),
     *              id%root%MBLOCK, id%root%MYROW, 0,
     *              id%root%NPROW )
          LOCAL_M = MAX(1, LOCAL_M)
          LOCAL_N = NUMROC( id%ND_STEPS(id%STEP(KEEP(38))),
     *              id%root%NBLOCK, id%root%MYCOL, 0,
     *              id%root%NPCOL )
        ELSE
          LOCAL_M = 0
          LOCAL_N = 0
        END IF
        IF  ( KEEP(60) .EQ. 2 .OR. KEEP(60) .EQ. 3 ) THEN
          id%SCHUR_MLOC=LOCAL_M
          id%SCHUR_NLOC=LOCAL_N
          id%root%SCHUR_MLOC=LOCAL_M
          id%root%SCHUR_NLOC=LOCAL_N
        ENDIF
        IF ( .NOT. ASSOCIATED(id%CANDIDATES)) THEN
          ALLOCATE(id%CANDIDATES(id%NSLAVES+1,1))
        ENDIF
        NB_NIV2 = MAX(KEEP(56), 1) 
        CALL DMUMPS_246( id%MYID_NODES, id%N,
     *    id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     *    id%NA, id%LNA, id%NE_STEPS(1), id%DAD_STEPS(1),
     *    id%ND_STEPS(1), id%PROCNODE_STEPS(1),
     *    id%NSLAVES,
     *    INFO( 3 ), INFO( 4 ), KEEP(15),
     *    KEEP(16),
     *    id%KEEP(200),
     *    INFO(5), RINFO(1),
     *    KEEP(1), KEEP8, LOCAL_M, LOCAL_N,
     *    SBUF_SEND, SBUF_REC, id%COST_SUBTREES, KEEP(28),
     *    id%I_AM_CAND(1), NB_NIV2, id%ISTEP_TO_INIV2,
     *    id%CANDIDATES, 
     *    INFO(1), INFO(2)
     $    ,id%MAX_SURF_MASTER
     *    ,MAX_SIZE_FACTOR
     $     )
      KEEP(202)=MAX_SIZE_FACTOR
      CALL MPI_ALLREDUCE (KEEP(202), KEEP(212), 1, 
     *     MPI_INTEGER, MPI_MAX, 
     *     id%COMM_NODES, IERR)
      IF (KEEP(205) .GT. 0) THEN
        KEEP(100) = KEEP(205)
      ELSE
         LOC_DIM_BUF_IO=KEEP(212)+INT(MAX(KEEP(12),0)*
     $     (DBLE(KEEP(212))/DBLE(100)+DBLE(1)))
         KEEP(100)= MIN(LOC_DIM_BUF_IO*2, 12000000)
      ENDIF
      IF (KEEP(99) .LT.3) KEEP(100)=0
        K44MAX = 0
        CALL MPI_ALLREDUCE (SBUF_REC, K44MAX, 1, 
     *    MPI_INTEGER, MPI_MAX, 
     *    id%COMM_NODES, IERR)
        SBUF_REC  = K44MAX
        IF (KEEP(48)==5) SBUF_SEND=SBUF_REC
        KEEP(66) = INT(RINFO(1)/1000000)
        KEEP( 43 ) = SBUF_SEND
        KEEP( 44 ) = SBUF_REC
        KEEP( 25 ) = INFO( 3 )
        KEEP( 26 ) = INFO( 4 )
        KEEP( 27 ) = INFO( 5 )
        KEEP( 29 ) = KEEP(15) + 2* MAX(KEEP(12),10)
     *                          * ( KEEP(15) / 100 + 1)
        INFO( 7 )  = KEEP(29)
        KEEP(30)   = KEEP(16) + KEEP(12) * ( KEEP(16) / 100 + 1 )
        INFO( 8 )  = KEEP(30)
        KEEP(206) = KEEP(200) + KEEP(12) * ( KEEP(200) /100 +1)
        IF ( MP .GT. 0 ) THEN
          WRITE(MP,'(A,I10) ') 
     *    ' estimated INTEGER space for factors         :',
     *                KEEP(26)
          WRITE(MP,'(A,I10) ') 
     *    ' estimated REAL    space for factors         :',
     *                KEEP(25)
          WRITE(MP,'(A,I10) ') 
     *    ' Minimum value of MAXIS for factorization    :',
     *                KEEP(15)
          WRITE(MP,'(A,I10) ') 
     *    ' Minimum value of MAXS  for IC factorization :',
     *                KEEP(16)
          IF (KEEP(207) == 1) THEN
          WRITE(MP,'(A,I10) ') 
     *    ' Minimum value of MAXS  for OOC factorization:',
     *                KEEP(200)
          ENDIF
          WRITE(MP,'(A,I10) ') 
     *    ' Current value of space relaxation parameter :',
     *                KEEP(12)
          WRITE(MP,'(A,I10) ') 
     *    ' Relaxed Value of MAXIS for factorization    :',
     *                KEEP(29)
          WRITE(MP,'(A,I10) ') 
     *    ' Relaxed Value of MAXS for IC factorization  :',
     *                KEEP(30)
          IF (KEEP(207)==1) THEN
          WRITE(MP,'(A,I10) ') 
     *    ' Relaxed Value of MAXS for OOC factorization :',
     *                KEEP(206)
          ENDIF
          WRITE(MP,'(A,I10) ') 
     *    ' Estimated maximum SENT message size (bytes) :',
     *                KEEP(43) * KEEP(35)
          WRITE(MP,'(A,I10) ') 
     *    ' Estimated maximum RECV message size (bytes) :',
     *                KEEP(44) * KEEP(35)
        END IF
      ELSE
       KEEP(30) = 0
       KEEP(29) = 0
       KEEP(206)= 0
      END IF
      CALL DMUMPS_100(id%MYID, id%COMM, KEEP, KEEP8,
     *     INFO, INFOG, RINFO, RINFOG, ICNTL)
      IF ( id%MYID. eq. MASTER .AND. KEEP(54) .eq. 1 ) THEN
        IF (ASSOCIATED( id%MAPPING))
     *        DEALLOCATE( id%MAPPING)
        ALLOCATE( id%MAPPING(id%NZ), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -7
          INFO(2) = id%NZ
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'id%MAPPING'
          END IF
          GOTO 92
        END IF
        ALLOCATE(IWtemp( id%N ), stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-7
          INFO(2)=id%N
          IF ( LP .GT. 0 ) THEN
            WRITE(LP, 150) 'IWtemp(N)'
          END IF
          GOTO 92
        END IF
        CALL DMUMPS_83(
     *    id%N, id%MAPPING(1),
     *    id%NZ, id%IRN(1),id%JCN(1), id%PROCNODE_STEPS(1),
     *    id%STEP(1),
     *    id%NSLAVES, id%SYM_PERM(1),
     *    id%FILS, IWtemp, id%KEEP(1),id%KEEP8(1),
     *    id%root%MBLOCK, id%root%NBLOCK,
     *    id%root%NPROW, id%root%NPCOL )
        DEALLOCATE( IWtemp )
 92     CONTINUE
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                  id%COMM, id%MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      OOC_STAT = .FALSE.
      PERLU_ON = .FALSE.  ! switch off PERLU to compute KEEP8(2)
      CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *                    id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *                    id%NA_ELT,
     *                    id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     *                    OOC_STAT, PERLU_ON, TOTAL_BYTES)
      id%KEEP8(2) = TOTAL_BYTES
      PERLU_ON  = .TRUE.
      CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *                    id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *                    id%NA_ELT,
     *                    id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     *                    OOC_STAT, PERLU_ON, TOTAL_BYTES)
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     *' Estimated space in MBYTES for IC factorization            :',
     *                TOTAL_MBYTES
      END IF
      id%INFO(15) = TOTAL_MBYTES
      CALL DMUMPS_243( id%MYID, id%COMM,
     *                           id%INFO(15), id%INFOG(16), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     *' ** Rank of proc needing largest memory in IC facto        :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *' ** Estimated corresponding MBYTES for IC facto            :',
     *  id%INFOG(16)
        IF ( KEEP(46) .eq. 0 ) THEN
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     *  ,(id%INFOG(17)-id%INFO(15))/id%NSLAVES
        ELSE
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. MBYTES per work. proc at facto (IC)     :'
     *  ,id%INFOG(17)/id%NSLAVES
        END IF
        WRITE(MPG,'(A,I10) ')
     *' ** TOTAL     space in MBYTES for IC factorization         :'
     * ,id%INFOG(17)
      END IF
      OOC_STAT = .TRUE.
      PERLU_ON = .FALSE.  ! PERLU NOT taken into account
      CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *                    id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *                    id%NA_ELT,
     *                    id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     *                    OOC_STAT, PERLU_ON, TOTAL_BYTES)
      id%KEEP8(3) = TOTAL_BYTES
      PERLU_ON  = .TRUE.  ! PERLU taken into account
      CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *                    id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *                    id%NA_ELT,
     *                    id%NSLAVES, TOTAL_MBYTES, .FALSE.,
     *                    OOC_STAT, PERLU_ON, TOTAL_BYTES)
      id%INFO(17) = TOTAL_MBYTES
      CALL DMUMPS_243( id%MYID, id%COMM,
     *                           id%INFO(17), id%INFOG(26), IRANK )
      IF ( PROKG .AND. KEEP(207)==1 ) THEN
        WRITE( MPG,'(A,I10) ')
     *' ** Rank of proc needing largest memory for OOC facto      :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *' ** Estimated corresponding MBYTES for OOC facto           :',
     *  id%INFOG(26)
        IF ( KEEP(46) .eq. 0 ) THEN
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     *  ,(id%INFOG(27)-id%INFO(15))/id%NSLAVES
        ELSE
         WRITE( MPG,'(A,I10) ')
     *' ** Estimated avg. MBYTES per work. proc at facto (OOC)    :'
     *  ,id%INFOG(27)/id%NSLAVES
        END IF
        WRITE(MPG,'(A,I10) ')
     *' ** TOTAL     space in MBYTES for OOC factorization        :'
     * ,id%INFOG(27)
      END IF
      RETURN
 110  FORMAT(/' ****** ANALYSIS STEP ********'/)
 150  FORMAT(
     */' ** FAILURE DURING DMUMPS_26, DYNAMIC ALLOCATION OF',
     *A30)
      END SUBROUTINE DMUMPS_26
        SUBROUTINE DMUMPS_537(N,NSLAVES,
     *              ICNTL,INFOG, NE, NFSIZ,
     *              FRERE, FILS,
     *              KEEP,KEEP8,PROCNODE,
     *              SSARBR,NBSA,PEAK,ierr)
        USE DMUMPS_CV
        IMPLICIT NONE
        INTEGER N, NSLAVES, NBSA, ierr
        INTEGER ICNTL(40),INFOG(40),KEEP(500)
      INTEGER*8 KEEP8(150)
        INTEGER NE(N),NFSIZ(N),FRERE(N),FILS(N),PROCNODE(N)
        INTEGER SSARBR(NBSA)
        DOUBLE PRECISION PEAK
        CALL DMUMPS_369(N,NSLAVES,
     *              ICNTL,INFOG, NE, NFSIZ,
     *              FRERE, FILS,
     *              KEEP,KEEP8,PROCNODE,
     *              SSARBR,NBSA,PEAK,ierr)
        RETURN
        END SUBROUTINE DMUMPS_537
        SUBROUTINE DMUMPS_564(INODE, PROCNODE, VALUE, FILS, N)
        INTEGER, intent(in) :: INODE, N, VALUE
        INTEGER, intent(in) :: FILS(N)
        INTEGER, intent(inout) :: PROCNODE(N)
        INTEGER IN
        IN=INODE
        DO WHILE ( IN > 0 )
           PROCNODE( IN ) = VALUE
           IN=FILS( IN )
        ENDDO
        RETURN
        END SUBROUTINE DMUMPS_564
      SUBROUTINE DMUMPS_251(N,IW,LIW,A,LA,
     *             NSTK_STEPS, NBPROCFILS,IFLAG,ND,FILS,STEP,
     *             FRERE, DAD, CAND, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *             MAXFRT, NTOTPV,PTRIST,PTRAST,
     *             PIMASTER, PAMASTER, PTRARW,PTRAIW,
     *             ITLOC, IERROR,IPOOL, LPOOL, NIRBDU,
     *             RINFO, POSFAC ,IWPOS, LRLU, IPTRLU, 
     *             LRLUS, LEAF, NBROOT, NBRTOT,
     *             UU, ICNTL, PTLUST_S, PTRFAC, NSTEPS, INFO,
     *             KEEP,KEEP8,
     *             PROCNODE_STEPS,SLAVEF,MYID, COMM_NODES,
     *             MYID_NODES,
     *             BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,root,
     *             PERM, NELT, FRTPTR, FRTELT, LPTRAR, 
     *             COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     *             MEM_DISTRIB, NE ,
     *     DKEEP,PIVNUL_LIST,LPN_LIST)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER N,IFLAG,NTOTPV,MAXFRT,LA,LIW, LPTRAR,
     *        IERROR,NIRBDU, NSTEPS, INFO(40)
      DOUBLE PRECISION A(LA)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      INTEGER LPOOL 
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N)
      INTEGER IW(LIW), NSTK_STEPS(KEEP(28)), NBPROCFILS(KEEP(28))
      INTEGER PTRARW(LPTRAR), PTRAIW(LPTRAR), ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28)),PTRAST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER PIMASTER(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)), PERM(N)
      INTEGER CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER IPOOL(LPOOL)
      INTEGER NE(KEEP(28))
      DOUBLE PRECISION RINFO(20)
      INTEGER POSFAC,IWPOS,LRLU,
     *    IPTRLU, LRLUS,
     *    LEAF, NBROOT 
      INTEGER COMM_LOAD, ASS_IRECV
      DOUBLE PRECISION UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER        INTARR( MAX(1,KEEP(14)) )
      DOUBLE PRECISION DBLARR( MAX(1,KEEP(13)) )
      LOGICAL IS_ISOLATED_NODE
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(30)
      INTEGER  LIWK_RR, LWK_RR, LOCAL_M, LOCAL_N, PHASE,
     &         MBLOCK, NBLOCK
      LOGICAL ROOT_OWNER
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      LOGICAL FLAG
      DOUBLE PRECISION ZERO
      INTEGER INODE
      INTEGER IWPOSCB
      INTEGER FPERE, TYPEF
      INTEGER MP, LP, DUMMY(1)
      INTEGER NBFIN, NBRTOT, NBROOT_TRAITEES
      INTEGER NFRONT, IOLDPS, POSELT
      INTEGER IPOSROOT
      INTEGER GLOBK109
      INTEGER LBUFRX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUFRX
      DOUBLE PRECISION FLOP1
      INTEGER TYPE
      LOGICAL SON_LEVEL2, SET_IRECV, BLOCKING, 
     &        MESSAGE_RECEIVED
      LOGICAL AVOID_DELAYED
      EXTERNAL DMUMPS_330, DMUMPS_275 
      INTEGER DMUMPS_330, DMUMPS_275
      LOGICAL DMUMPS_167,DMUMPS_283
      EXTERNAL DMUMPS_167,DMUMPS_283
      LOGICAL DMUMPS_508
      EXTERNAL DMUMPS_508, DMUMPS_509
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, COMP,
     &        JOBASS, ETATASS
      LOGICAL CONCERNED
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION    OPASSW, OPELIW
      DATA ZERO /0.0D0/
       ASS_IRECV = MPI_REQUEST_NULL
       ITLOC(1:N) =0
       PTRIST  (1:KEEP(28))=0
       PTLUST_S(1:KEEP(28))=0
       PTRAST(1:KEEP(28))=0
       PTRFAC(1:KEEP(28))=-99999
      MP = ICNTL(2)
      LP = ICNTL(1)
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      COMP  = 0
      OPASSW = ZERO
      OPELIW = ZERO
      IWPOSCB = NIRBDU    
      CALL DMUMPS_22(
     *     .FALSE., .FALSE., MYID_NODES, N, KEEP, KEEP8,
     *     IW, LIW, A, LA, NIRBDU, LRLU, IPTRLU, IWPOS, IWPOSCB,
     *     PTRIST, PTRAST, STEP, PIMASTER,
     *     PAMASTER, ITLOC, XSIZE, 0, -444, -444, .true.,
     *     COMP, LRLUS,
     *     IFLAG, IERROR)
      JOBASS  = 0
      ETATASS = 0
      NBFIN = NBRTOT
      NBROOT_TRAITEES = 0
      NBPROCFILS(1:KEEP(28)) = 0
      IF ( KEEP(38).NE.0 ) THEN
        IF (root%yes) THEN
          IF (KEEP(60).eq.0) THEN
            CALL DMUMPS_284( 
     *        root, KEEP(38), N, IW, LIW,
     *        A, LA,
     *        NIRBDU, FILS, MYID_NODES, PTRAIW, PTRARW,
     *        INTARR, DBLARR,
     *        LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,
     *        COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
          ELSE
             PTRIST(STEP(KEEP(38)))=-77777
          ENDIF
        ENDIF
        IF ( IFLAG .LT. 0 ) GOTO 635
      END IF
 20   CONTINUE
      NIV1_FLAG=0
      SET_IRECV = .TRUE.
      BLOCKING = .FALSE.
      MESSAGE_RECEIVED = .FALSE.
      CALL DMUMPS_329(
     *      COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV, 
     *      MESSAGE_RECEIVED,
     *      MPI_ANY_SOURCE, MPI_ANY_TAG,
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC,
     *      IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, 
     *      NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS,
     *      COMP, IFLAG,
     *      IERROR, COMM_NODES,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID_NODES, SLAVEF,
     *
     *      root, OPASSW, OPELIW, ITLOC, FILS, 
     *      PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      CALL DMUMPS_467(COMM_LOAD, KEEP)
      IF (MESSAGE_RECEIVED) THEN
          IF ( IFLAG .LT. 0 ) GO TO 640
          IF ( NBFIN .eq. 0 ) GOTO 640
      ELSE
          IF ( .NOT. DMUMPS_508( IPOOL, LPOOL) )THEN
            CALL DMUMPS_509( N, IPOOL, LPOOL,
     *      PROCNODE_STEPS,
     *      SLAVEF, STEP, INODE, KEEP,KEEP8, MYID_NODES, ND )
            IF (KEEP(47) .GE. 3) THEN
              CALL DMUMPS_500(
     $              IPOOL, LPOOL, 
     *              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *              MYID_NODES, STEP, N, ND, FILS )
            ENDIF
            IF (KEEP(47).EQ.4) THEN
               IF(INODE.GT.0.AND.INODE.LE.N)THEN
                  IF((NE(STEP(INODE)).EQ.0).AND.
     $                 (FRERE(STEP(INODE)).EQ.0))THEN
                     IS_ISOLATED_NODE=.TRUE.
                  ELSE
                     IS_ISOLATED_NODE=.FALSE.
                  ENDIF
               ENDIF
               CALL DMUMPS_501(
     *              IS_ISOLATED_NODE,INODE,IPOOL,LPOOL,
     *              MYID_NODES,slavef,COMM_LOAD,KEEP,KEEP8)
            ENDIF
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     *           ( KEEP(47) == 4 )).OR.
     $           (KEEP(80) == 1 .AND. KEEP(47) .GE. 1)) THEN
               CALL DMUMPS_512(INODE,STEP,KEEP(28),
     *         PROCNODE_STEPS,FRERE,ND,COMM_LOAD,SLAVEF,
     *         MYID_NODES,KEEP,KEEP8,N)
            END IF
            GOTO 30
          ENDIF
      ENDIF
      GO TO 20
 30   CONTINUE
      IF ( INODE .LT. 0 ) THEN
        INODE = -INODE
        FPERE = DAD(STEP(INODE))
        GOTO 130
      ELSE IF (INODE.GT.N) THEN
       INODE = INODE - N
       IF (INODE.EQ.KEEP(38)) THEN
         NBROOT_TRAITEES = NBROOT_TRAITEES + 1
         IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID_NODES,
     *          COMM_NODES, RACINE, SLAVEF)
            END IF
         ENDIF
         IF (NBFIN.EQ.0) GOTO 640
         GOTO 20
       ENDIF
       TYPE = DMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
       IF (TYPE.EQ.1) GOTO 100
       FPERE = DAD(STEP(INODE))
       AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     *                   .AND. KEEP(60).ne.0 )
       IF ( KEEP(50) .eq. 0 ) THEN
         CALL  DMUMPS_144( COMM_LOAD, ASS_IRECV,
     *             N, INODE, FPERE, IW, LIW, A, LA,
     *             UU, NOFFW,
     *             NPVW,
     *             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     *             NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     *             STEP, PIMASTER, PAMASTER,
     *             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     *             root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED, 
     *              DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
          IF ( IFLAG .LT. 0 ) GOTO 640
       ELSE
         CALL  DMUMPS_141( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, A, LA,
     *             UU, NOFFW,
     *             NPVW,
     *             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     *             NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     *             STEP, PIMASTER, PAMASTER,
     *             NSTK_STEPS,NBPROCFILS,PROCNODE_STEPS,
     *             root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL_LDLT_NIV2, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     *              DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
        IF ( IFLAG .LT. 0 ) GOTO 640
        IF ( IW( PTLUST_S(STEP(INODE)) + XSIZE + 5 ) .GT. 1 ) THEN
             GOTO 20
        END IF
       END IF
       GOTO 130
      ENDIF
      IF (INODE.EQ.KEEP(38)) THEN
         CALL  DMUMPS_176( COMM_LOAD, ASS_IRECV, 
     *    ROOT, FRERE,
     *    INODE,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS, COMP,
     *    IFLAG, IERROR, COMM_NODES,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID_NODES, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
         IF ( IFLAG .LT. 0 ) GOTO 640
        GOTO 20
      ENDIF
      TYPE = DMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (TYPE.EQ.1) THEN
        IF (KEEP(55).NE.0) THEN
         CALL DMUMPS_36( COMM_LOAD, ASS_IRECV,
     *        NELT, FRTPTR, FRTELT,
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *     PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE )
        ELSE
         JOBASS = 0
         CALL DMUMPS_252(COMM_LOAD, ASS_IRECV, 
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    JOBASS,ETATASS )
        ENDIF
       IF ( IFLAG .LT. 0 ) GOTO 640
        IF ((NBPROCFILS(STEP(INODE)).GT.0).OR.(SON_LEVEL2)) GOTO 20
      ELSE
        IF ( KEEP(55) .eq. 0 ) THEN
          CALL DMUMPS_253(COMM_LOAD, ASS_IRECV, 
     *    N, INODE, IW, LIW, A, LA,
     *    IFLAG, IERROR,
     *    ND, NIRBDU, FILS, FRERE, CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW,
     *    root, OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     *    PTRAIW, ITLOC, NSTEPS,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     *    ICNTL, KEEP,KEEP8,INTARR,DBLARR,
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     *    MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     *    MEM_DISTRIB(0)
     *    )
        ELSE
          CALL DMUMPS_37( COMM_LOAD, ASS_IRECV, 
     *    NELT, FRTPTR, FRTELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG, IERROR,
     *    ND, NIRBDU, FILS, FRERE, CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW,
     *    root, OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     *    PTRAIW, ITLOC, NSTEPS,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     *    ICNTL, KEEP,KEEP8,INTARR,DBLARR,
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     *    MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES,
     *    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     *    MEM_DISTRIB(0))
        END IF
        IF (IFLAG.LT.0) GOTO 640
        GOTO 20
      ENDIF
 100  CONTINUE
       FPERE = DAD(STEP(INODE))
      IF ( INODE .eq. KEEP(20) ) THEN
        POSELT = PTRAST(STEP(INODE))
        IF (PTRFAC(STEP(INODE)).NE.POSELT) THEN
          WRITE(*,*) "ERROR 2 in DMUMPS_251", POSELT
          CALL DMUMPS_ABORT()
        ENDIF
        CALL DMUMPS_87( IW(PTLUST_S(STEP(INODE))+XSIZE) )
        GOTO 200
      END IF
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST_S(STEP(INODE))
      AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     *                   .AND. KEEP(60).ne.0 )
      IF (KEEP(50).EQ.0) THEN
         CALL DMUMPS_143( N, INODE, IW, LIW, A, LA,
     *               IOLDPS, POSELT,
     *               IFLAG, UU, NOFFW, NPVW,
     *               KEEP,KEEP8,
     *               STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     *               SEUIL, AVOID_DELAYED, ETATASS,
     *              DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
         JOBASS = ETATASS
         IF (JOBASS.EQ.1) THEN
           CALL DMUMPS_252(COMM_LOAD, ASS_IRECV, 
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    JOBASS,ETATASS )
         ENDIF
      ELSE
         IW( IOLDPS+4+XSIZE ) = 1
         CALL DMUMPS_140( N, INODE, IW, LIW, A, LA,
     *           IOLDPS, POSELT,
     *           IFLAG, UU, NOFFW, NPVW,
     *           KEEP,KEEP8, MYID_NODES, SEUIL, AVOID_DELAYED, 
     *           ETATASS,
     *        DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
         IW( IOLDPS+4+XSIZE ) = STEP(INODE)
         JOBASS = ETATASS  ! ETATASS=1 indicates
         IF (JOBASS.EQ.1) THEN
           CALL DMUMPS_252(COMM_LOAD, ASS_IRECV, 
     *        N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,NIRBDU,
     *        FILS,FRERE,MAXFRW,root,OPASSW, OPELIW,
     *      PTRIST,PTLUST_S,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     *        PTRARW,PTRAIW,
     *        ITLOC, NSTEPS, SON_LEVEL2,
     *        COMP, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, POSFAC, LRLUS,
     *        ICNTL, KEEP,KEEP8, INTARR, DBLARR,
     *    NSTK_STEPS,NBPROCFILS, PROCNODE_STEPS, SLAVEF,
     *    COMM_NODES, MYID_NODES,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    JOBASS,ETATASS )
         ENDIF
      ENDIF
      IF (IFLAG.LT.0) GOTO 635
 130  CONTINUE
      TYPE  = DMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF ( FPERE .NE. 0 ) THEN
        TYPEF = DMUMPS_330(STEP(FPERE),PROCNODE_STEPS,SLAVEF)
      ELSE
        TYPEF = -9999
      END IF
      CALL DMUMPS_254( COMM_LOAD, ASS_IRECV, 
     *       N,INODE,TYPE,TYPEF,LA,IW,LIW,A,
     *       IFLAG,IERROR,OPELIW,NELVAW,
     *       PTRIST,PTLUST_S,PTRFAC,
     *       PTRAST, STEP, PIMASTER, PAMASTER,
     *       NE, POSFAC,LRLU,
     *       LRLUS,IPTRLU,ICNTL,KEEP,KEEP8,NIRBDU,COMP,IWPOS,IWPOSCB,
     *       PROCNODE_STEPS,SLAVEF,FPERE,COMM_NODES,MYID_NODES,
     *       IPOOL, LPOOL, LEAF,
     *       NSTK_STEPS, NBPROCFILS, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     * root, OPASSW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     * ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IF (IFLAG.LT.0) GOTO 640
 200  CONTINUE
      IF ( INODE .eq. KEEP(38) ) THEN
        WRITE(*,*) 'Error .. in DMUMPS_251: ',
     &             ' INODE == KEEP(38)'
        Stop
      END IF
      IF ( FPERE.EQ.0 ) THEN
        NBROOT_TRAITEES = NBROOT_TRAITEES + 1
        IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
           IF (KEEP(201).NE.0) THEN
            CALL DMUMPS_580(IERR)              
           ENDIF
            NBFIN = NBFIN - NBROOT
            IF ( NBFIN .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_251: ', 
     &                   ' NBFIN=', NBFIN 
              CALL DMUMPS_ABORT()
            END IF
            IF ( NBROOT .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_251: ', 
     &                   ' NBROOT=', NBROOT
              CALL DMUMPS_ABORT()
            END IF
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_242( DUMMY(1), 1, MPI_INTEGER,
     *          MYID_NODES, COMM_NODES, RACINE, SLAVEF)
            END IF
        ENDIF
        IF (NBFIN.EQ.0)THEN
           GOTO 640
        ENDIF
      ELSEIF ( FPERE.NE.KEEP(38) .AND.
     *         DMUMPS_275(STEP(FPERE),PROCNODE_STEPS,SLAVEF).EQ.
     *         MYID_NODES ) THEN
        NSTK_STEPS(STEP(FPERE)) = NSTK_STEPS(STEP(FPERE))-1
        IF ( NSTK_STEPS( STEP( FPERE )).EQ.0) THEN
          CALL DMUMPS_507(N, IPOOL, LPOOL,
     *         PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *         KEEP(80), KEEP(47), STEP, FPERE )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *            MYID_NODES, STEP, N, ND, FILS )
          ENDIF
          CALL DMUMPS_137( FPERE, N, PROCNODE_STEPS,SLAVEF,
     &           ND, FILS, FRERE, STEP, PIMASTER, KEEP(28),
     *           KEEP(50), FLOP1,
     &           IW, LIW )
          IF (FPERE.NE.KEEP(20))
     *    CALL DMUMPS_190(1,.FALSE.,FLOP1,KEEP,KEEP8)
        ENDIF
      ENDIF
      GO TO 20
 635  CONTINUE
      CALL DMUMPS_44( MYID_NODES, SLAVEF, COMM_NODES )
 640  CONTINUE
        CALL DMUMPS_255( INFO(1),
     *       ASS_IRECV, BUFR, LBUFR,
     *       LBUFR_BYTES,
     *       COMM_NODES,
     *       MYID_NODES, SLAVEF)
       CALL DMUMPS_180( INFO(1),
     *      BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      COMM_NODES, COMM_LOAD, SLAVEF, MP)
      CALL MPI_BARRIER( COMM_NODES, IERR )
       IF ( INFO(1) .GE. 0 ) THEN
          IF( KEEP(38) .NE. 0 .OR. KEEP(20).NE.0) THEN
             IF ( KEEP(38) .NE. 0 )THEN
                IF ( LRLU .GT. LBUFR ) THEN
#if defined (null_space_old)
                   IF ( KEEP(19) .ne. 0 ) THEN
                      CALL DMUMPS_147( MYID_NODES, 
     *                     root, N, 
     *                     KEEP(38), COMM_NODES, IW, LIW, IWPOS + 1,
     *                     A, LA, PTRAST, PTLUST_S, PTRFAC, STEP,
     *                     INFO(1), KEEP(50), KEEP(19), A(POSFAC),
     *                     LRLU, KEEP,KEEP8,A(POSFAC),LRLU,
     *                     A(POSFAC+3*root%tot_root_size+3),
     *                     A(POSFAC+2*root%tot_root_size+2))
                   ELSE
#endif
                      CALL DMUMPS_146( MYID_NODES, 
     *                     root, N, 
     *                     KEEP(38), COMM_NODES, IW, LIW, IWPOS + 1,
     *                     A, LA, PTRAST, PTLUST_S, PTRFAC, STEP,
     *                     INFO(1), KEEP(50), KEEP(19), A(POSFAC),
     *                     LRLU, KEEP,KEEP8 )
#if defined(null_space_old)
                   END IF
#endif
                ELSE
#if defined (null_space_old)
                   IF ( KEEP(19) .ne. 0 ) THEN
                      NFRONT = IW(PTLUST_S(STEP(KEEP(38)))+XSIZE)
                      OPELIW = OPELIW + 
     *                     ( DBLE(4)*DBLE(NFRONT)*DBLE(NFRONT)*
     *                     DBLE(NFRONT)/DBLE(3) ) /
     *                     DBLE(SLAVEF)
                      CALL DMUMPS_147( MYID_NODES,
     *                     root, N, KEEP(38),
     *                     COMM_NODES, IW, LIW, IWPOS + 1,
     *                     A, LA, PTRAST, PTLUST_S, PTRFAC, STEP,
     *                     INFO(1), KEEP(50), KEEP(19),
     *                     BUFR, MAX(LBUFR/2,1), KEEP,KEEP8,
     *                     BUFR, MAX(LBUFR/2,1),
     *                     BUFR(2*(3*root%tot_root_size+4)),
     *                     BUFR(2*(2*root%tot_root_size+3)))
                   ELSE
#endif
                      LBUFRX = MIN(root%MBLOCK * root%NBLOCK,
     *                     root%TOT_ROOT_SIZE*root%TOT_ROOT_SIZE )
                      ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                      IF (IERR.gt.0) THEN
                         INFO(1) = -9
                         INFO(2) =  LBUFRX
                         write(LP,*) ' Error allocating, real array ',
     1                        'of size before DMUMPS_146',  LBUFRX
                         stop
                      ENDIF
                      CALL DMUMPS_146( MYID_NODES,
     *                     root, N, KEEP(38),
     *                     COMM_NODES, IW, LIW, IWPOS + 1,
     *                     A, LA, PTRAST, PTLUST_S, PTRFAC, STEP,
     *                     INFO(1), KEEP(50), KEEP(19),
     *                     BUFRX, LBUFRX, KEEP,KEEP8 )
                      DEALLOCATE ( BUFRX )
#if defined(null_space_old)
                   END IF
#endif
                ENDIF
                IF ( MYID_NODES .eq. DMUMPS_275(STEP(KEEP(38)),
     *               PROCNODE_STEPS,SLAVEF)) THEN
                   IF ( INFO(1) .eq. -10 ) THEN
                      NPVW = NPVW + INFO(2)
                   ELSE
                      NPVW = NPVW + root%tot_root_size
                   END IF
                END IF
                IF (root%yes.AND.KEEP(60).EQ.0) THEN
                   CALL DMUMPS_471(.FALSE.,.FALSE.,
     *                  LA-LRLUS
     *                  ,IW(PTLUST_S(STEP(KEEP(38))) + 2 +XSIZE)*
     *                  IW( PTLUST_S(STEP(KEEP(38))) + 1 +XSIZE),0,
     *                  KEEP,KEEP8,LRLU)
                  IF (KEEP(201) .NE.0) THEN
                    KEEP(31)=KEEP(31)+
     $              (IW(PTLUST_S(STEP(KEEP(38)))+2+XSIZE)*
     $              IW(PTLUST_S(STEP(KEEP(38)))+1+XSIZE))
                    CALL DMUMPS_576(KEEP(38),PTRFAC,
     $              KEEP,KEEP8,A,LA,
     $              (IW(PTLUST_S(STEP(KEEP(38)))+2+XSIZE)*
     $              IW(PTLUST_S(STEP(KEEP(38)))+1+XSIZE)),IERR)
                    IF(IERR.LT.0)THEN
               WRITE(*,*)MYID,': Internal error in DMUMPS_576'
                     CALL DMUMPS_ABORT()
                    ENDIF
                  ENDIF
                ENDIF
             ELSE
                IF (KEEP(19).NE.0) THEN
                  CALL MPI_REDUCE(KEEP(109), GLOBK109, 1,
     *                 MPI_INTEGER, MPI_SUM,
     *                 DMUMPS_275(STEP(KEEP(20)),
     *                           PROCNODE_STEPS,SLAVEF),
     *                 COMM_NODES, IERR)
                ENDIF
                CONCERNED =  (MYID_NODES .EQ. 
     *                        DMUMPS_275(STEP(KEEP(20)),
     *                              PROCNODE_STEPS,SLAVEF))
                IF (CONCERNED ) THEN
                   IPOSROOT = PTLUST_S(STEP(KEEP(20)))
                   NFRONT   = IW(IPOSROOT+XSIZE)
                   NPVW = NPVW + NFRONT
                   CALL DMUMPS_471(.FALSE.,.FALSE.,LA-LRLUS,
     *                  NFRONT*NFRONT,0,KEEP,KEEP8,LRLU)
                END IF
#if defined(try_null_space)
                IF ( KEEP(60) .EQ. 0 ) THEN
                   IF ( CONCERNED ) THEN
                    LOCAL_M = 0
                    LOCAL_N = 0
                    MBLOCK  = 0
                    NBLOCK  = 0
                    PHASE   = 1  ! compute effective size
                    ROOT_OWNER = .TRUE.
                    CALL DMUMPS_279( PHASE, 
     *               MBLOCK, NBLOCK, NFRONT, LOCAL_M, LOCAL_N,
     *               ROOT_OWNER, KEEP,KEEP8,
     *               LIWK_RR, LWK_RR )
                    IF (LIWK_RR.NE.0) THEN
                 WRITE(*,*) "Error after DMUMPS_279",LIWK_RR
                      CALL DMUMPS_ABORT()
                    ENDIF
                      IF ( LRLU .GT. LWK_RR ) THEN
                         CALL DMUMPS_295(
     1                        NFRONT,A(PTRAST(STEP(KEEP(20)))),
     2                        root,
     3                        A(POSFAC), LRLU, 
     4                        KEEP,KEEP8,INFO, LP, DKEEP,
     5                        GLOBK109, OPELIW, NOFFW)
                      ELSE
                         LBUFRX = LWK_RR
                         ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                         IF (IERR.gt.0) THEN
                            INFO(1) = -13
                            INFO(2) =  LBUFRX
                            write(LP,*) ' Error allocating, real
     $                           array ','of size ', LBUFRX
                            CALL DMUMPS_ABORT()
                         ENDIF
                         CALL DMUMPS_295(
     1                        NFRONT,A(PTRAST(STEP(KEEP(20)))),
     2                        root,
     3                        BUFRX(1), LBUFRX, 
     4                        KEEP,KEEP8, INFO, LP, DKEEP,
     5                        GLOBK109, OPELIW, NOFFW)
                         DEALLOCATE ( BUFRX )
                      ENDIF
                      IF (KEEP(201).NE.0) THEN
                        KEEP(31)=KEEP(31)+NFRONT*NFRONT
                        CALL DMUMPS_576(KEEP(20),PTRFAC,
     $                  KEEP,KEEP8,A,LA, NFRONT*NFRONT, IERR)
                        IF(IERR.LT.0)THEN
              WRITE(*,*)MYID,': Internal error in DMUMPS_576'
                        CALL DMUMPS_ABORT()
                        ENDIF
                      ENDIF
                   ENDIF
                   CALL MPI_BCAST( KEEP(17), 1, MPI_INTEGER,
     *                  DMUMPS_275(STEP(KEEP(20)),
     $                  PROCNODE_STEPS,SLAVEF),
     *                  COMM_NODES, IERR )
                END IF
#endif
             END IF
          END IF
       END IF
       IF ( KEEP(38) .NE. 0 ) THEN
         IF (MYID_NODES.EQ. DMUMPS_275(STEP(KEEP(38)),
     *           PROCNODE_STEPS,SLAVEF)) THEN
           MAXFRW = MAX ( MAXFRW, root%TOT_ROOT_SIZE)
         END IF
       END IF
       MAXFRT       = MAXFRW
       NTOTPV       = NPVW
       INFO(12)     = NOFFW
       RINFO(2)     = DBLE(OPASSW)
       RINFO(3)     = DBLE(OPELIW)
       INFO(13)     = NELVAW
       INFO(14)     = COMP
      RETURN
      END SUBROUTINE DMUMPS_251
      SUBROUTINE DMUMPS_87( HEADER )
        INTEGER HEADER( 6 )
        INTEGER NFRONT, NASS
        NFRONT = HEADER(1)
        IF ( HEADER(2) .ne. 0 ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 1 :',HEADER(2)
          CALL DMUMPS_ABORT()
        END IF
        NASS   = IABS( HEADER( 3 ) )
        IF ( NASS .NE. IABS( HEADER( 4 ) ) ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 2 :',HEADER(3:4)
          CALL DMUMPS_ABORT()
        END IF
        IF ( NASS .NE. NFRONT ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 3 : not root'
          CALL DMUMPS_ABORT()
        END IF
        HEADER( 1 ) = NFRONT
        HEADER( 2 ) = 0
        HEADER( 3 ) = NFRONT
        HEADER( 4 ) = NFRONT
        RETURN
      END SUBROUTINE DMUMPS_87
      SUBROUTINE DMUMPS_136( id )
      USE DMUMPS_OOC
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_BUFFER
      IMPLICIT NONE
      include 'mpif.h'
      TYPE( DMUMPS_STRUC ) :: id
      LOGICAL I_AM_SLAVE
      INTEGER IERR, MASTER, STATUS( MPI_STATUS_SIZE )
      INTEGER MAXS_BYTES
      PARAMETER ( MASTER = 0 )
      I_AM_SLAVE = ( id%MYID .ne. MASTER .OR. id%KEEP(46) .NE. 0 )
      IF (id%KEEP(201).NE.0 .AND. I_AM_SLAVE) THEN
        CALL DMUMPS_587(id,IERR)
      END IF
      IF (id%root%gridinit_done) THEN
        IF ( id%KEEP(38).NE.0 .and. id%root%yes ) THEN
          CALL BLACS_GRIDEXIT( id%root%CNTXT_BLACS )
          id%root%gridinit_done = .FALSE.
        END IF
      END IF
      IF ( id%MYID .NE. MASTER .OR. id%KEEP(46) .ne. 0 ) THEN
        CALL MPI_COMM_FREE( id%COMM_NODES, IERR )
        CALL MPI_COMM_FREE( id%COMM_LOAD, IERR )
      END IF
      IF (ASSOCIATED(id%MEM_DIST))  THEN
         DEALLOCATE(id%MEM_DIST)
         NULLIFY(id%MEM_DIST)
      ENDIF
      IF (ASSOCIATED(id%MAPPING)) THEN
        DEALLOCATE(id%MAPPING)
        NULLIFY(id%MAPPING)
      END IF
      IF (ASSOCIATED(id%NULL_SPACE)) THEN
        DEALLOCATE(id%NULL_SPACE)
        NULLIFY(id%NULL_SPACE)
      ENDIF
       NULLIFY(id%SCHUR_CINTERFACE)
      IF ( id%KEEP(52) .NE. -1 .or. id%MYID .ne. MASTER ) THEN
        IF (ASSOCIATED(id%COLSCA)) THEN
          DEALLOCATE(id%COLSCA)
          NULLIFY(id%COLSCA)
        ENDIF
        IF (ASSOCIATED(id%ROWSCA)) THEN
          DEALLOCATE(id%ROWSCA)
          NULLIFY(id%ROWSCA)
        ENDIF
      END IF
      IF (ASSOCIATED(id%PTLUST_S)) THEN
        DEALLOCATE(id%PTLUST_S)
        NULLIFY(id%PTLUST_S)
      END IF
      IF (ASSOCIATED(id%PTRFAC)) THEN
        DEALLOCATE(id%PTRFAC)
        NULLIFY(id%PTRFAC)
      END IF
      IF (ASSOCIATED(id%POIDS)) THEN
        DEALLOCATE(id%POIDS)
        NULLIFY(id%POIDS)
      ENDIF
      IF (ASSOCIATED(id%IS)) THEN
        DEALLOCATE(id%IS)
        NULLIFY(id%IS)
      ENDIF
      IF (ASSOCIATED(id%IS1)) THEN
        DEALLOCATE(id%IS1)
        NULLIFY(id%IS1)
      ENDIF
      IF (ASSOCIATED(id%STEP))      THEN
        DEALLOCATE(id%STEP)
        NULLIFY(id%STEP)
      ENDIF
      IF (ASSOCIATED(id%NE_STEPS))  THEN
        DEALLOCATE(id%NE_STEPS)
        NULLIFY(id%NE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%ND_STEPS))  THEN
        DEALLOCATE(id%ND_STEPS)
        NULLIFY(id%ND_STEPS)
      ENDIF
      IF (ASSOCIATED(id%FRERE_STEPS))  THEN
        DEALLOCATE(id%FRERE_STEPS)
        NULLIFY(id%FRERE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%DAD_STEPS))  THEN
        DEALLOCATE(id%DAD_STEPS)
        NULLIFY(id%DAD_STEPS)
      ENDIF
      IF (ASSOCIATED(id%SYM_PERM))  THEN
        DEALLOCATE(id%SYM_PERM)
        NULLIFY(id%SYM_PERM)
      ENDIF
      IF (ASSOCIATED(id%UNS_PERM))  THEN
        DEALLOCATE(id%UNS_PERM)
        NULLIFY(id%UNS_PERM)
      ENDIF
      IF (ASSOCIATED(id%PIVNUL_LIST))  THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
      IF (ASSOCIATED(id%FILS))      THEN
        DEALLOCATE(id%FILS)
        NULLIFY(id%FILS)
      ENDIF
      IF (ASSOCIATED(id%PTRAR))     THEN
        DEALLOCATE(id%PTRAR)
        NULLIFY(id%PTRAR)
      ENDIF
      IF (ASSOCIATED(id%FRTPTR))    THEN
        DEALLOCATE(id%FRTPTR)
        NULLIFY(id%FRTPTR)
      ENDIF
      IF (ASSOCIATED(id%FRTELT))    THEN
        DEALLOCATE(id%FRTELT)
        NULLIFY(id%FRTELT)
      ENDIF
      IF (ASSOCIATED(id%NA))        THEN
        DEALLOCATE(id%NA)
        NULLIFY(id%NA)
      ENDIF
      IF (ASSOCIATED(id%PROCNODE_STEPS)) THEN
        DEALLOCATE(id%PROCNODE_STEPS)
        NULLIFY(id%PROCNODE_STEPS)
      ENDIF
      IF (ASSOCIATED(id%PROCNODE)) THEN
        DEALLOCATE(id%PROCNODE)
        NULLIFY(id%PROCNODE)
      ENDIF
      IF (ASSOCIATED(id%RHSCOMP)) THEN
        DEALLOCATE(id%RHSCOMP)
        NULLIFY(id%RHSCOMP)
      ENDIF
      IF (ASSOCIATED(id%POSINRHSCOMP)) THEN
        DEALLOCATE(id%POSINRHSCOMP)
        NULLIFY(id%POSINRHSCOMP)
      ENDIF
      IF (id%KEEP(46).eq.1 .and.
     *    id%KEEP(55).ne.0 .and.
     *    id%MYID .eq. MASTER .and.
     *    id%KEEP(52) .eq. 0 ) THEN
        NULLIFY(id%DBLARR)
      ELSE
        IF (ASSOCIATED(id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
      END IF
      IF (ASSOCIATED(id%INTARR))       THEN
        DEALLOCATE(id%INTARR)
        NULLIFY(id%INTARR)
      ENDIF
      IF (ASSOCIATED(id%root%RG2L_ROW))THEN
        DEALLOCATE(id%root%RG2L_ROW)
        NULLIFY(id%root%RG2L_ROW)
      ENDIF
      IF (ASSOCIATED(id%root%RG2L_COL))THEN
        DEALLOCATE(id%root%RG2L_COL)
        NULLIFY(id%root%RG2L_COL)
      ENDIF
      IF (ASSOCIATED(id%root%IPIV))    THEN
        DEALLOCATE(id%root%IPIV)
        NULLIFY(id%root%IPIV)
      ENDIF
      CALL DMUMPS_636(id)
      IF (ASSOCIATED(id%ELTPROC))     THEN
        DEALLOCATE(id%ELTPROC)
        NULLIFY(id%ELTPROC)
      ENDIF
      IF (ASSOCIATED(id%CANDIDATES)) THEN
        DEALLOCATE(id%CANDIDATES)
        NULLIFY(id%CANDIDATES)
      ENDIF
      IF (ASSOCIATED(id%I_AM_CAND)) THEN
        DEALLOCATE(id%I_AM_CAND)
        NULLIFY(id%I_AM_CAND)
      ENDIF
      IF (ASSOCIATED(id%ISTEP_TO_INIV2)) THEN
        DEALLOCATE(id%ISTEP_TO_INIV2)
        NULLIFY(id%ISTEP_TO_INIV2)
      ENDIF
      IF (I_AM_SLAVE) THEN
       IF (ASSOCIATED(id%TAB_POS_IN_PERE)) THEN
        DEALLOCATE(id%TAB_POS_IN_PERE)
        NULLIFY(id%TAB_POS_IN_PERE)
       ENDIF
       IF (ASSOCIATED(id%FUTURE_NIV2)) THEN
        DEALLOCATE(id%FUTURE_NIV2)
        NULLIFY(id%FUTURE_NIV2)
       ENDIF
      ENDIF
      IF(ASSOCIATED(id%DEPTH_FIRST))THEN
        DEALLOCATE(id%DEPTH_FIRST)
        NULLIFY(id%DEPTH_FIRST)
      ENDIF
      IF (ASSOCIATED(id%MEM_SUBTREE)) THEN
        DEALLOCATE(id%MEM_SUBTREE)
        NULLIFY(id%MEM_SUBTREE)
      ENDIF
      IF (ASSOCIATED(id%MY_ROOT_SBTR)) THEN
        DEALLOCATE(id%MY_ROOT_SBTR)
        NULLIFY(id%MY_ROOT_SBTR)
      ENDIF
      IF (ASSOCIATED(id%MY_FIRST_LEAF)) THEN
        DEALLOCATE(id%MY_FIRST_LEAF)
        NULLIFY(id%MY_FIRST_LEAF)
      ENDIF
      IF (ASSOCIATED(id%MY_NB_LEAF)) THEN
        DEALLOCATE(id%MY_NB_LEAF)
        NULLIFY(id%MY_NB_LEAF)
      ENDIF
      IF (ASSOCIATED(id%COST_TRAV)) THEN
        DEALLOCATE(id%COST_TRAV)
        NULLIFY(id%COST_TRAV)
      ENDIF
        IF(ASSOCIATED (id%OOC_INODE_SEQUENCE))THEN
         DEALLOCATE(id%OOC_INODE_SEQUENCE)
         NULLIFY(id%OOC_INODE_SEQUENCE)
        ENDIF
        IF(ASSOCIATED (id%OOC_NUM_FILE))THEN
         DEALLOCATE(id%OOC_NUM_FILE)
         NULLIFY(id%OOC_NUM_FILE)
        ENDIF
        IF(ASSOCIATED (id%OOC_POS_IN_FILE))THEN
         DEALLOCATE(id%OOC_POS_IN_FILE)
         NULLIFY(id%OOC_POS_IN_FILE)
        ENDIF
        IF(ASSOCIATED (id%OOC_SIZE_OF_BLOCK))THEN
         DEALLOCATE(id%OOC_SIZE_OF_BLOCK)
         NULLIFY(id%OOC_SIZE_OF_BLOCK)
        ENDIF
      IF (ASSOCIATED(id%S))        DEALLOCATE(id%S)
      NULLIFY(id%S)
      IF (I_AM_SLAVE) THEN
        CALL DMUMPS_57( IERR )
        CALL DMUMPS_59( IERR )
      END IF
      IF ( ASSOCIATED( id%BUFR ) ) DEALLOCATE( id%BUFR )
      NULLIFY( id%BUFR )
      RETURN
      END SUBROUTINE DMUMPS_136
      SUBROUTINE DMUMPS_150(MYID,COMM,S,MAXS,MAXS_BYTES)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER COMM, MYID, MAXS, MAXS_BYTES
      INTEGER S( MAXS )
      INTEGER MSGTAG, MSGSOU, MSGLEN
      LOGICAL FLAG
      FLAG = .TRUE.
      DO WHILE ( FLAG )
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
        IF (FLAG) THEN
          MSGTAG=STATUS(MPI_TAG)
          MSGSOU=STATUS(MPI_SOURCE)
          CALL MPI_GET_COUNT(STATUS,MPI_PACKED,MSGLEN,IERR)
          IF (MSGLEN <= MAXS_BYTES) THEN
            CALL MPI_RECV(S(1),MAXS_BYTES,MPI_PACKED,
     *      MSGSOU, MSGTAG, COMM, STATUS, IERR)
          ELSE
            EXIT
          ENDIF
        END IF
      END DO
      CALL MPI_BARRIER( COMM, IERR )
      RETURN
      END SUBROUTINE DMUMPS_150
      SUBROUTINE DMUMPS_254(COMM_LOAD, ASS_IRECV, 
     *    N, INODE, TYPE, TYPEF, 
     *    LA, IW, LIW, A,
     *    IFLAG, IERROR, OPELIW, NELVAW, PTRIST, PTLUST_S,
     *    PTRFAC, PTRAST,
     *    STEP, PIMASTER, PAMASTER, NE,
     *    POSFAC, LRLU, LRLUS, IPTRLU, ICNTL, KEEP,KEEP8,
     *    NIRBDU, COMP, IWPOS, IWPOSCB, PROCNODE_STEPS, SLAVEF,
     *    FPERE, COMM, MYID,
     *    IPOOL, LPOOL, LEAF, NSTK_S,
     *    NBPROCFILS,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, root,
     *    OPASSW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *    ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE DMUMPS_BUFFER
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER COMM, MYID, TYPE, TYPEF
      INTEGER N, LA, LIW, INODE,IFLAG,IERROR
      INTEGER ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POSFAC, LRLU,LRLUS, IPTRLU,NIRBDU,IWPOSCB,IWPOS,
     &        FPERE, SLAVEF, NELVAW
      INTEGER IW(LIW),PROCNODE_STEPS(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), NE(KEEP(28))
      DOUBLE PRECISION    A(LA)
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(MAX(1,KEEP(13)))
      INTEGER INTARR(MAX(1,KEEP(14))) 
      INTEGER ITLOC( N ), FILS( N ),
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER LPOOL, LEAF, COMP
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NBFIN
      INTEGER IWRITEFAC, IREADFAC,NFRONT_ESTIM,NELIM_ESTIM
      INTEGER DMUMPS_275
      EXTERNAL DMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER APOS, POSELT, OPSFAC, NSLAVES, NCB,
     *        H_INODE, IERR, IERR_MPI, NBCOL, NBROW, NBROW_SEND,
     *        NBROW_STACK, NELIM
      INTEGER SHIFT_LIST_ROW_SON,
     *        SHIFT_LIST_COL_SON, SHIFT_VAL_SON,
     *        LIST_ROW_SON, LIST_COL_SON, LIST_SLAVES
      INTEGER IOLDPS,NFRONT,NPIV,NASS,IOLDP1,PTROWEND,
     *        LREQI, LREQCB_SAMEPROC, LREQCB_OTHERPROC, LCONT
      INTEGER I,NPOS,J,J1,LDA, INIV2
      INTEGER MSGDEST, MSGTAG, CHK_LOAD
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mumps_headers.h'
      LOGICAL  COMPRESSCB
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL SSARBR, SSARBR_ROOT, DMUMPS_167,
     *DMUMPS_170
      EXTERNAL DMUMPS_167, DMUMPS_170
      IOLDPS = PTLUST_S(STEP(INODE))
      NFRONT = IW(IOLDPS+XSIZE)
      NPIV   = IW(IOLDPS + 1+XSIZE)
      NASS   = IABS(IW(IOLDPS + 2+XSIZE))
      NSLAVES=  IW(IOLDPS+5+XSIZE)
      H_INODE= 6 + NSLAVES + XSIZE
      LCONT = NFRONT - NPIV
      NBCOL = LCONT
      SSARBR = DMUMPS_167(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      SSARBR_ROOT = DMUMPS_170
     *              (STEP(INODE),PROCNODE_STEPS,SLAVEF)
      COMPRESSCB= ((KEEP(215).EQ.0)
     *             .AND.(KEEP(50).NE.0)
     *             .AND.(TYPEF.EQ.1
     *             .OR.TYPEF.EQ.2
     *              )
     *             .AND.(TYPE.EQ.1))
      IF ((FPERE.EQ.0).AND.(NASS.NE.NPIV)) THEN 
        IFLAG = -10
        GOTO 600
      ENDIF
      NBROW      = LCONT
      IF (TYPE.EQ.2) NBROW = NASS - NPIV
      IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
        LDA = NASS
      ELSE
        LDA = NFRONT
      ENDIF
      NBROW_SEND = NBROW
      IF (TYPEF.EQ.2) NBROW_SEND = NASS - NPIV  ! only uneliminated rows
      IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
        LREQCB_SAMEPROC  = NBROW * NBROW 
      ELSE
        LREQCB_SAMEPROC  = NBROW * NBCOL
      ENDIF
      IF (COMPRESSCB) THEN
         LREQCB_SAMEPROC  = (LCONT*(LCONT+1))/2
      ENDIF
      LREQCB_OTHERPROC = (NBROW-NBROW_SEND) * NBCOL
      IF (COMPRESSCB) THEN
        LREQCB_OTHERPROC = (LCONT*(LCONT+1))/2 -
     &                     (NBROW_SEND*(NBROW_SEND+1))/2
      ENDIF
      NELIM = NASS-NPIV
      POSELT = PTRAST(STEP(INODE))
      IF (POSELT .ne. PTRFAC(STEP(INODE))) THEN
        WRITE(*,*) "Error 1 in G"
        CALL DMUMPS_ABORT()
      END IF
      NELVAW = NELVAW + NASS - NPIV
      CALL DMUMPS_511( NFRONT, NPIV, NASS,
     *     KEEP(50), TYPE,FLOP1 )
      IF ( (.NOT. SSARBR_ROOT) .and. TYPE == 1) THEN
        IF (NE(STEP(INODE))==0) THEN
          CHK_LOAD=0
        ELSE
          CHK_LOAD=1
        ENDIF
        CALL DMUMPS_190(CHK_LOAD, .FALSE., -FLOP1,
     *                      KEEP,KEEP8)
      ENDIF
      FLOP1_EFFECTIVE = FLOP1
      OPELIW = OPELIW + FLOP1
      IF ( NPIV .NE. NASS ) THEN
        CALL DMUMPS_511( NFRONT, NASS, NASS,
     *       KEEP(50), TYPE,FLOP1 )
         IF (.NOT. SSARBR_ROOT ) THEN
            IF (NE(STEP(INODE))==0) THEN
              CHK_LOAD=0
            ELSE
              CHK_LOAD=1
            ENDIF
            CALL DMUMPS_190(CHK_LOAD, .FALSE.,
     *                         FLOP1_EFFECTIVE-FLOP1,
     *                         KEEP,KEEP8)
         ENDIF
      END IF
      IF ( SSARBR_ROOT ) THEN
        NFRONT_ESTIM=ND(STEP(INODE))
        NELIM_ESTIM=NASS-(NFRONT-NFRONT_ESTIM)
        CALL DMUMPS_511(NFRONT_ESTIM,NELIM_ESTIM,NELIM_ESTIM,
     *       KEEP(50),1,FLOP1)
      END IF
      FLOP1=-FLOP1
      IF (SSARBR_ROOT) THEN
        CALL DMUMPS_190(0,.FALSE.,FLOP1,KEEP,KEEP8)
      ELSE
        CALL DMUMPS_190(2,.FALSE.,FLOP1,KEEP,KEEP8)
      ENDIF
      IF ( FPERE .EQ. 0 ) GOTO 200
      IF ( FPERE.EQ.KEEP(38) ) THEN
       NCB   = NFRONT - NASS
       SHIFT_LIST_ROW_SON = H_INODE + NASS
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NASS
       SHIFT_VAL_SON      = NASS*NFRONT + NASS
       IF (TYPE.EQ.1) THEN
         CALL DMUMPS_80(
     *    COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE,
     *    PTLUST_S, PTRAST,
     *    ROOT, NCB, NCB, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     *    ROOT_CONT_STATIC, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF (IFLAG < 0 ) GOTO 500
       ENDIF
       MSGDEST=  DMUMPS_275(STEP(FPERE),PROCNODE_STEPS,SLAVEF)
       IOLDPS = PTLUST_S(STEP(INODE))
       LIST_ROW_SON = IOLDPS + H_INODE + NPIV
       LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
       LIST_SLAVES  = IOLDPS + 6 + XSIZE
       IF (MSGDEST.EQ.MYID) THEN
         CALL DMUMPS_273( root, 
     *      INODE, NELIM, NSLAVES, IW(LIST_ROW_SON), 
     *      IW(LIST_COL_SON), IW(LIST_SLAVES),
     *
     *      PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *      PTLUST_S, PTRFAC,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, ITLOC, COMP,
     *      IFLAG, IERROR, 
     *      IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,KEEP8,
     *      COMM, COMM_LOAD, FILS, ND)
         IF (IFLAG.LT.0) GOTO 600
       ELSE
        IERR = -1
        DO WHILE (IERR.EQ.-1)
         CALL DMUMPS_76( INODE, NELIM, 
     *     IW(LIST_ROW_SON), IW(LIST_COL_SON), NSLAVES, 
     *     IW(LIST_SLAVES), MSGDEST, COMM, IERR)
         IF ( IERR .EQ. -1 ) THEN
            BLOCKING   =.FALSE.
            SET_IRECV  =.TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL DMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *       MPI_ANY_SOURCE, MPI_ANY_TAG,
     *       STATUS,
     *       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *       IWPOS, IWPOSCB, IPTRLU,
     *       LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *       PTLUST_S, PTRFAC,
     *       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *       IFLAG, IERROR, COMM,
     *       NBPROCFILS,
     *       IPOOL, LPOOL, LEAF,
     *       NBFIN, MYID, SLAVEF,
     *
     *       root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP,KEEP8,
     *       ND, FRERE, LPTRAR, NELT,
     *       FRTPTR, FRTELT, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
            IF ( IFLAG .LT. 0 ) GOTO 500
            IOLDPS = PTLUST_S(STEP(INODE))
            LIST_ROW_SON = IOLDPS + H_INODE + NPIV
            LIST_COL_SON = IOLDPS + H_INODE + NFRONT + NPIV
            LIST_SLAVES  = IOLDPS + 6 + XSIZE
         ENDIF
        ENDDO
        IF ( IERR .EQ. -2 ) THEN 
            IERROR = ( 3 + NSLAVES + 2 * NELIM ) * KEEP( 34 ) 
            IFLAG  = - 17
            GOTO 600
        ELSE IF ( IERR .EQ. -3 ) THEN
            IERROR = ( 3 + NSLAVES + 2 * NELIM ) * KEEP( 34 ) 
            IFLAG  = -20
            GOTO 600
        ENDIF
       ENDIF
       IF (NELIM.EQ.0) THEN 
          POSELT = PTRAST(STEP(INODE))
          OPSFAC = POSELT + NPIV * NFRONT + NPIV 
          GOTO 190
       ELSE
          GOTO 500
       ENDIF
      ENDIF
      OPSFAC = POSELT + NPIV * LDA + NPIV
      IF ( DMUMPS_275(STEP(FPERE),PROCNODE_STEPS,
     *     SLAVEF) .EQ. MYID ) THEN
        LREQI = 2 + XSIZE
        CALL DMUMPS_22(SSARBR, .FALSE.,
     *   MYID,N,KEEP,KEEP8,IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQI, LREQCB_SAMEPROC, INODE, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR )
        IF (IFLAG.LT.0) GOTO 600
        NPOS = IPTRLU + 1
        PAMASTER(STEP(INODE)) = NPOS
        PTRAST(STEP(INODE)) = -99999999
        IF ((KEEP(50).NE.0).AND.(TYPE.EQ.2)) THEN
          IF (NBROW.GT.0) THEN
            APOS = OPSFAC
            DO I = 1, NBROW
             J1 = APOS
             DO J = 1, NBROW
               A(NPOS) = A(J1)
               NPOS = NPOS + 1
               J1 = J1 + 1
             ENDDO
             APOS = APOS + LDA
            ENDDO
          ENDIF
          PIMASTER (STEP(INODE)) = PTLUST_S(STEP(INODE))
          PTRIST(STEP(INODE)) = IWPOSCB+1
          IW(IWPOSCB+1+XSIZE) = MIN(-NBROW,-1)
          IW(IWPOSCB+2+XSIZE) = NBROW
        ELSE
          IF ((KEEP(50).NE.0).AND.(NPIV.GT.0).AND.(NFRONT.NE.NPIV)) THEN
            IWRITEFAC = POSELT + NPIV
            IREADFAC  = POSELT + NFRONT
            IF ( NPIV < KEEP(59) ) THEN
             DO I = 1, NPIV - 1
              DO J = 0, NPIV-1
                A(IWRITEFAC+J)=A(IREADFAC+J)
              ENDDO
              IREADFAC = IREADFAC + NFRONT
              IWRITEFAC = IWRITEFAC + NPIV
             ENDDO
            ELSE
             DO I = 1, NPIV - 2
                DO J = 0, I+1
                   A(IWRITEFAC+J)=A(IREADFAC+J)
                ENDDO
                IREADFAC = IREADFAC + NFRONT
                IWRITEFAC = IWRITEFAC + NPIV
             ENDDO
             I = NPIV - 1
             DO J = 0, I
                A(IWRITEFAC+J)=A(IREADFAC+J)
             ENDDO
             IREADFAC = IREADFAC + NFRONT
             IWRITEFAC = IWRITEFAC + NPIV
            ENDIF
          ENDIF
          IF (NBROW.GT.0) THEN
            APOS  = OPSFAC
            DO I = 1, NBROW
              IF ( KEEP(50) .eq. 0 ) THEN
                CALL DCOPY( NBCOL, A(APOS), 1, A(NPOS), 1 )
                IF  ( I == 1 ) THEN
                   IWRITEFAC  = OPSFAC
                   IREADFAC   = POSELT + (NPIV + 1) * NFRONT
                ELSE IF (NPIV.GT.0) THEN
                   DO J=0,NPIV-1
                      A(IWRITEFAC + J) = A(IREADFAC + J)
                   ENDDO
                   IWRITEFAC=IWRITEFAC+NPIV
                   IREADFAC=IREADFAC+NPIV+NBCOL
                  END IF
              ELSE
                CALL DCOPY( I    , A(APOS), 1, A(NPOS), 1 )
                IF (NPIV.GT.0) THEN
                  DO J = 0, NPIV - 1
                    A(IWRITEFAC+J) = A(IREADFAC+J)
                  ENDDO
                  IWRITEFAC = IWRITEFAC + NPIV
                  IREADFAC = IREADFAC + NFRONT
                ENDIF
#if ! defined(ALLOW_NON_INIT)
                IF (.NOT.COMPRESSCB) THEN
                  A(NPOS + I: NPOS + NBCOL - 1) = DBLE(ZERO)
                ENDIF
#endif
              END IF
              IF (COMPRESSCB) THEN
               NPOS = NPOS + I
              ELSE
               NPOS = NPOS + NBCOL
              ENDIF
              APOS = APOS + LDA
            ENDDO
          ENDIF
          PIMASTER (STEP(INODE)) = PTLUST_S(STEP(INODE))
          PTRIST(STEP(INODE)) = IWPOSCB+1
          IF (COMPRESSCB) IW(IWPOSCB+1+XXS) = S_CB1COMP
          IW(IWPOSCB+1+XSIZE) = MIN( -NBCOL, -1 )
          IW(IWPOSCB+2+XSIZE) = NBROW
          IF ( NBCOL .eq. 0 .and. NBROW .ne. 0 ) THEN
            write(*,*) 'error: more rows than columns',NBCOL,NBROW
            stop
          END IF
          GOTO 200
        ENDIF
      ELSE
        MSGTAG =NOEUD
        MSGDEST=DMUMPS_275( STEP(FPERE), PROCNODE_STEPS, SLAVEF )
        IERR = -1 
        DO WHILE (IERR.EQ.-1)
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
             CALL DMUMPS_66( INODE, FPERE, NFRONT, 
     *         LCONT, NASS, NPIV, IW( IOLDPS +  H_INODE + NPIV ),
     *         IW( IOLDPS +  H_INODE + NPIV + NFRONT ),
     *         A( OPSFAC ), COMPRESSCB, LREQCB_SAMEPROC, 
     *         MSGDEST, MSGTAG, COMM, IERR )
          ELSE
             IF ( (TYPE.EQ.2).AND.(KEEP(48).NE.0)) THEN
              INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
             ELSE
              INIV2 = -9999
             ENDIF
             CALL DMUMPS_70( FPERE, INODE, 
     *             NBROW_SEND, IW(IOLDPS +  H_INODE + NPIV ), 
     *             NBCOL, IW(IOLDPS +  H_INODE + NPIV + NFRONT ),
     *             A(OPSFAC), LDA, NELIM, TYPE,
     *             NSLAVES, IW(IOLDPS+6+XSIZE), MSGDEST, 
     *             COMM, IERR, 
     * 
     *             SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
          END IF
          IF ( IERR .EQ. -1 ) THEN
            BLOCKING  = .FALSE.
            SET_IRECV = .TRUE.
            MESSAGE_RECEIVED = .FALSE.
            CALL DMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *       BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *       MPI_ANY_SOURCE, MPI_ANY_TAG,
     *       STATUS,
     *       BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *       IWPOS, IWPOSCB, IPTRLU,
     *       LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *       PTRIST, PTLUST_S, PTRFAC,
     *       PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *       IFLAG, IERROR, COMM,
     *       NBPROCFILS,
     *       IPOOL, LPOOL, LEAF,
     *       NBFIN, MYID, SLAVEF,
     *
     *       root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *       LPTRAR, NELT, FRTPTR, FRTELT,
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
             IF ( IFLAG .LT. 0 ) GOTO 500
          ENDIF
          IOLDPS = PTLUST_S(STEP( INODE ))
        END DO
        IF ( IERR .EQ. -2 .OR. IERR .EQ. -3 ) THEN 
          IF ( (TYPE.EQ.1) .AND. (TYPEF.EQ.1) ) THEN
            IERROR = ( 2*LCONT + 9 ) * KEEP( 34 ) + 
     *          LCONT*LCONT * KEEP( 35 )
          ELSE IF (KEEP(50).ne.0 .AND. TYPE .eq. 2 ) THEN
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) 
     *                 * KEEP( 34 ) + 
     *          NBROW_SEND*NBROW_SEND*KEEP( 35 )
          ELSE
            IERROR = ( NBROW_SEND + NBCOL+ 5 + NSLAVES) * KEEP( 34 ) + 
     *          NBROW_SEND*NBCOL*KEEP( 35 )
          ENDIF
          IF (IERR .EQ. -2) IFLAG = -17
          IF (IERR .EQ. -3) IFLAG = -20
          GOTO 600
        ENDIF
        IF ((FPERE.EQ.0).OR.(LREQCB_OTHERPROC.EQ.0)) GOTO 190
        NBROW_STACK = NBROW-NBROW_SEND
        LREQI       = 6 + NBROW_STACK + NBCOL + XSIZE
        CALL DMUMPS_22(SSARBR, .FALSE.,
     *   MYID,N,KEEP,KEEP8,IW,LIW,A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQI, LREQCB_OTHERPROC, INODE, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR )
        IF (IFLAG.LT.0) GOTO 600
        PTRIST(STEP(INODE)) = IWPOSCB+1
        IF (COMPRESSCB) IW(IWPOSCB+1+XXS)=S_CB1COMP
        IW(IWPOSCB+1+XSIZE) = NBCOL
        IW(IWPOSCB+2+XSIZE) = 0
        IW(IWPOSCB+3+XSIZE) = NBROW_STACK
        IW(IWPOSCB+4+XSIZE) = 0
        IW(IWPOSCB+5+XSIZE) = 1
        IW(IWPOSCB+6+XSIZE) = 0
        IOLDP1   = PTLUST_S(STEP(INODE))+H_INODE
        PTROWEND = IWPOSCB+6+NBROW_STACK+XSIZE
        IW(IWPOSCB+7+XSIZE:PTROWEND) =
     *    IW(IOLDP1+NFRONT-NBROW_STACK:IOLDP1+NFRONT-1)
        IW(PTROWEND+1:PTROWEND+NBCOL) =
     *    IW(IOLDP1+NFRONT+NPIV: IOLDP1+2*NFRONT-1)
        PTRAST(STEP(INODE)) = IPTRLU+1
        NPOS = IPTRLU+1
        APOS = OPSFAC + NBROW_SEND*NFRONT
        IF ( NELIM .ne. NBCOL-NBROW_STACK ) THEN
          write(*,*) 'Internal error in DMUMPS_254',
     *  NBCOL, NBROW_STACK, NELIM
          CALL DMUMPS_ABORT()
        ENDIF
        DO I = 1, NBROW_STACK
          IF ( KEEP(50) .eq. 0 ) THEN
            CALL DCOPY( NBCOL, A(APOS), 1, A(NPOS), 1 )
          ELSE
            CALL DCOPY( I+NELIM, A(APOS), 1, A(NPOS), 1 )
#if ! defined(ALLOW_NON_INIT)
            IF (.NOT. COMPRESSCB) THEN
              A(NPOS+I+NELIM:NPOS+NBCOL-1)= DBLE(ZERO)
            ENDIF
#endif
          END IF
          IF (COMPRESSCB) THEN
            NPOS=NPOS+I+NELIM
          ELSE
            NPOS = NPOS + NBCOL
          ENDIF
          APOS = APOS + NFRONT
        ENDDO
      END IF
  190 CONTINUE
      IF (KEEP(50).EQ.0) THEN
        IF ((NPIV .NE. 0) .AND. (NBROW .GT. 1)) THEN
          IF ( IOLDPS .NE. PTLUST_S( STEP(INODE)) ) THEN
            WRITE(*,*)
     *      ' ERROR 1 in 51G : IOLDPS,PTLUST_S(STEP(INODE))=',
     *      IOLDPS,PTLUST_S(STEP(INODE))
            CALL DMUMPS_ABORT()
          END IF
          IOLDPS = PTLUST_S(STEP(INODE))
          POSELT = PTRFAC(IW(IOLDPS + 4+XSIZE))
          APOS     = POSELT + (NPIV + 1) * NFRONT
          DO 180 I = 1, NBROW - 1
            DO J = 0, NPIV-1
              A(OPSFAC + J) = A(APOS + J)
            END DO
            OPSFAC = OPSFAC + NPIV
            APOS   = APOS + NPIV + NBCOL
  180     CONTINUE
        ENDIF
      ELSE
        IOLDPS = PTLUST_S(STEP(INODE))
        POSELT = PTRFAC(IW(IOLDPS + 4+XSIZE))
        IF ( TYPE .eq. 1 ) THEN
          CALL DMUMPS_324(A(POSELT), NFRONT, NPIV, NBROW,
     *                               KEEP(59))
        ELSE
          CALL DMUMPS_325(A(POSELT), NPIV, NBROW)
        ENDIF
      END IF
  200 CONTINUE
      IOLDPS = PTLUST_S(STEP(INODE))
      IW(IOLDPS+XSIZE)     = NBCOL
      IW(IOLDPS + 1+XSIZE) = NASS - NPIV
      IF (TYPE.EQ.2) THEN
        IW(IOLDPS + 2+XSIZE) = NASS
      ELSE
        IW(IOLDPS + 2+XSIZE) = NFRONT
      ENDIF
      IW(IOLDPS + 3+XSIZE) = NPIV
      CALL DMUMPS_93(MYID,N,IOLDPS,TYPE, IW, LIW, 
     *    A, LA, POSFAC, LRLU, LRLUS, 
     *    IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8, SSARBR,INODE,IERR)
      IF(IERR.LT.0)THEN
         IFLAG=IERR
         IERROR=0
         GOTO 600
      ENDIF
 500  CONTINUE
      RETURN
  600 CONTINUE
      IF (IFLAG .NE. -1) CALL DMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_254
      SUBROUTINE DMUMPS_142( id)
      USE DMUMPS_BUFFER
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
#ifndef SUN_
      INTERFACE
      SUBROUTINE DMUMPS_84( id )
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_84
      SUBROUTINE DMUMPS_27(id, ANORMINF, LSCAL)
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      DOUBLE PRECISION, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      END SUBROUTINE DMUMPS_27
      END INTERFACE
#endif
      TYPE(DMUMPS_STRUC), TARGET :: id
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, MASTER
      PARAMETER( MASTER = 0 )
      INCLUDE 'mumps_headers.h'
      INTEGER NSEND, NSEND_TOT, LDPTRAR, NELT
      INTEGER NLOCAL, NLOCAL_TOT, KEEP13_SAVE, K67, ITMP
      INTEGER  DMUMPS_275
      EXTERNAL DMUMPS_275
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL
      INTEGER DMUMPS_LBUF, DMUMPS_LBUFR_BYTES, DMUMPS_LBUF_INT
      INTEGER PTRIST, PTRWB, MAXELT_SIZE, LOC_DIM_BUF_IO,
     *        ITLOC, IPOOL, NSTEPS, K28, LPOOL, LIW
      INTEGER IRANK, ID_SCHUR, ID_ROOT, SIZE_SCHUR
      INTEGER KKKK
      INTEGER MAXS_BASE, MAXS_BASE_RELAXED, MEMORY_MD_ARG
      DOUBLE PRECISION CNTL4
      INTEGER MIN_PERLU
      DOUBLE PRECISION TIME
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER PERLU, TOTAL_MBYTES
      INTEGER ITOT
      LOGICAL I_AM_SLAVE, OOC_STAT, PERLU_ON
      DOUBLE PRECISION :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2
      DOUBLE PRECISION :: CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER DUMMY,I
      INTEGER, DIMENSION(:), POINTER :: IWK
      DOUBLE PRECISION, DIMENSION(:), POINTER :: WK
      INTEGER, POINTER :: MAXS, JOB, NZ, MAXIS
      DOUBLE PRECISION,DIMENSION(:),POINTER::RINFO, RINFOG
      DOUBLE PRECISION,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER::INFO, INFOG, KEEP
      INTEGER*8,DIMENSION(:),POINTER::KEEP8
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL DMUMPS_505
      INTEGER DMUMPS_505
      INTEGER SIZE_DEPTH_FIRST
      INTEGER*8 TOTAL_BYTES
#if defined(LARGEMATRICES)
      INTEGER NUMROC
      EXTERNAL NUMROC
#endif
      JOB=>id%JOB
      MAXS=>id%MAXS
      NZ=>id%NZ
      MAXIS=>id%MAXIS
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFO=>id%INFO
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      KEEP8=>id%KEEP8
      ICNTL=>id%ICNTL
      N = id%N
      EPS = EPSILON ( ZERO )
      IF (KEEP(219).NE.0) THEN
       CALL DMUMPS_617(MAX(KEEP(108),1),IERR)
       IF (IERR .NE. 0) THEN
          INFO(1) = -13
          INFO(2) = MAX(KEEP(108),1)
          GOTO 137
       END IF
      ENDIF
      KEEP13_SAVE = id%KEEP(13)
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = ( MP  .GT. 0 )
      PROKG = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
      IF ( PROKG .and. KEEP(53).GT.0 ) THEN
        WRITE(MPG,'(/A,I3)') ' Null space option :', KEEP(19)
        IF ( KEEP(21) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Max deficiency    : ', KEEP(21)
        END IF
        IF ( KEEP(22) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Min deficiency    : ', KEEP(22)
        END IF
      END IF
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
        IF (id%MYID .EQ. MASTER) THEN
          KEEP(201)=id%ICNTL(22)
        ENDIF
        CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        IF (id%MYID.EQ.MASTER) THEN
          IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
            KEEP(217)=0
          ENDIF
          KEEP(214)=KEEP(217)
          IF (KEEP(214).EQ.0) THEN
            IF (KEEP(201).NE.0) THEN
              KEEP(214)=1
            ELSE
              KEEP(214)=2
            ENDIF
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        IF (KEEP(201).NE.0) THEN
          CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
        ENDIF
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (id%CNTL(1) .ne. ZERO ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)')
     *' ** Warning : SPD solver called, resetting CNTL(1) to 0.0D0'
            END IF
          END IF
          id%CNTL(1) = ZERO
        END IF
      CALL DMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF (INFO(1).LT.0) GOTO 550
      IF ( PROKG ) THEN
          WRITE( MPG, 172 ) id%NSLAVES,
     *    KEEP8(5), KEEP(126), KEEP(127), KEEP(28)
      ENDIF
      IF ( id%MYID.EQ.MASTER ) THEN
        LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 6))
        IF ( LSCAL ) THEN
          IF ( KEEP(52) .eq. 5 .or. 
     &         KEEP(52) .eq. 6 ) THEN
            ITOT = 5 * id%N + NZ
          ELSE
            ITOT = 5 * id%N
          END IF
          ALLOCATE( WK( ITOT ), stat = IERR )
          IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=ITOT
            GOTO 137
          END IF
          CALL DMUMPS_217(id%N, NZ, KEEP(52), id%A(1),
     *             id%IRN(1), id%JCN(1),
     *             id%COLSCA(1), id%ROWSCA(1),
     *             WK, ITOT, ICNTL, INFO )
          DEALLOCATE( WK )
        END IF
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
      END IF
 137  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 550
      KEEP(110)=ICNTL(24)
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     *               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_DOUBLE_PRECISION,
     *               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_DOUBLE_PRECISION,
     *               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_DOUBLE_PRECISION,
     *               MASTER, id%COMM, IERR)
      ANORMINF = ZERO
      IF (KEEP(19).EQ.0) THEN 
         SEUIL = ZERO
      ELSE
         CALL DMUMPS_27(  id , ANORMINF, LSCAL )
         IF (CNTL6 .LT. ZERO) THEN
           SEUIL = EPS*ANORMINF
         ELSE
           SEUIL = CNTL6*ANORMINF
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     *   ' ABSOLUTE PIVOT THRESHOLD for rank revealing =',SEUIL
      ENDIF
      SEUIL_LDLT_NIV2 = SEUIL
      IF (KEEP(110).EQ.0) THEN 
         id%DKEEP(1) = -1.0D0
         id%DKEEP(2) = ZERO
      ELSE
         IF (ANORMINF.EQ.ZERO) 
     &       CALL DMUMPS_27(  id , ANORMINF, LSCAL )
         IF (CNTL3 .LE. ZERO) THEN
           id%DKEEP(1)  = 1.0E-5*EPS*ANORMINF
         ELSE
           id%DKEEP(1)  = CNTL3*ANORMINF
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     *    ' ZERO PIVOT DETECTION ON, THRESHOLD          =',id%DKEEP(1)
         IF (CNTL5.GT.ZERO) THEN
            id%DKEEP(2) = CNTL5 * ANORMINF
            IF (PROKG) WRITE(MPG,*) 
     *    ' FIXATION FOR NULL PIVOTS                    =',id%DKEEP(2)
         ELSE
            IF (PROKG) WRITE(MPG,*) 'INFINITE FIXATION '
            id%DKEEP(2) = ZERO
         ENDIF
      ENDIF
      IF (KEEP(53).NE.0) THEN
        ID_ROOT =DMUMPS_275(id%STEP(id%KEEP(20)),
     *            id%PROCNODE_STEPS(1),id%NSLAVES)
        IF ( id%KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
      IF ( ASSOCIATED( id%PIVNUL_LIST) ) DEALLOCATE(id%PIVNUL_LIST)
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = N
      ELSE
         LPN_LIST = 1
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     *   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = N
      ENDIF
      ALLOCATE( id%PIVNUL_LIST(LPN_LIST),stat = IERR )
      IF ( IERR .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=LPN_LIST
      END IF
      id%PIVNUL_LIST(1:LPN_LIST) = 0
      KEEP(109) = 0
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 550
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_DOUBLE_PRECISION,
     *                MASTER, id%COMM, IERR )
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
            IF(ANORMINF .EQ. ZERO) THEN
               CALL DMUMPS_27(  id , ANORMINF, LSCAL )
            ENDIF
            SEUIL = SQRT(EPS) * ANORMINF
         ELSE
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
      KEEP(98) = 0
      KEEP(103) = 0
      KEEP(105) = 0
      IF ( id%MYID.EQ.MASTER ) THEN
        ITMP = ICNTL(23)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     *                MASTER, id%COMM, IERR )
      KEEP8(4) = ITMP * 1000000
      PERLU = KEEP(12)
      IF (KEEP(201) .EQ. 0) THEN
        MAXS_BASE=KEEP(16)
      ELSE
        MAXS_BASE=KEEP(200)
      ENDIF
      IF ( MAXS_BASE .gt. 0 ) THEN
          MAXS_BASE_RELAXED =
     *         MAXS_BASE + PERLU * ( MAXS_BASE / 100 + 1)
          id%MAXS = MAX( id%MAXS, MAXS_BASE_RELAXED )
      ELSE
        id%MAXS = 1
        MAXS_BASE_RELAXED = 1
      END IF
      IF (I_AM_SLAVE) THEN
#if ! defined(SIMULATE_OOC_FAC)
        IF (KEEP(201) .NE. 0) THEN
#endif
          IF(id%KEEP(96).GT.0)THEN
            id%MAXS=id%KEEP(96)
          ELSE
            IF (KEEP8(4) .NE. 0) THEN
             PERLU_ON = .TRUE.
             IF (KEEP(201).NE.0) THEN 
               OOC_STAT = .TRUE.
             ELSE
               OOC_STAT = .FALSE.
             ENDIF
             CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *       id%MYID, id%N, id%NELT, id%LNA, id%NZ, id%NA_ELT,
     *       id%NSLAVES, TOTAL_MBYTES, .FALSE., OOC_STAT, 
     *       PERLU_ON, TOTAL_BYTES)
               id%MAXS=MAXS_BASE_RELAXED +
     *         NINT( DBLE(KEEP8(4)-TOTAL_BYTES)/DBLE(KEEP(35)) )
               IF (id%MAXS < 0) THEN
                 id%INFO(1)=-9
                 id%INFO(2)=id%KEEP(200)
               ENDIF
            ENDIF
          ENDIF
#if ! defined(SIMULATE_OOC_FAC)
        ENDIF
#endif
      ENDIF
      CALL DMUMPS_536(PROKG, MPG, id%MAXS, id%NSLAVES,
     * id%COMM, "effective relaxed value of MAXS          =")
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 550
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        CALL DMUMPS_188( DBLE(id%COST_SUBTREES),
     $        KEEP(64), KEEP(66),id%MAXS )
        K28=KEEP(28)
        MEMORY_MD_ARG = MIN(PERLU * ( MAXS_BASE / 100 + 1 ),
     $                      MAX(0, id%MAXS-MAXS_BASE))
        CALL DMUMPS_185( id, MEMORY_MD_ARG )
        CALL DMUMPS_587(id, IERR)
        IF (KEEP(205) .GT. 0) THEN
          KEEP(100) = KEEP(205)
        ELSE
           LOC_DIM_BUF_IO=KEEP(212)+INT(MAX(KEEP(12),0)*
     $     (DBLE(KEEP(212))/DBLE(100)+DBLE(1)))
           KEEP(100)=MIN(LOC_DIM_BUF_IO*2, 12000000)
        ENDIF
        IF (KEEP(99) .LT.3) KEEP(100)=0
        IF((DBLE(KEEP(100))*DBLE(KEEP(35))/DBLE(2)).GT.
     $     (DBLE(1999999999)))THEN
           WRITE(*,*)id%MYID,': Warning: DIM_BUF_IO might be
     $  too big for Filesystem'
        ENDIF
        IF (KEEP(201) .NE. 0) THEN
           ALLOCATE (id%OOC_INODE_SEQUENCE(id%KEEP(28)),
     $          stat=ierr)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 111
           ENDIF
           ALLOCATE (id%OOC_NUM_FILE(id%KEEP(28)),
     $          stat=ierr)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%KEEP(28)
              NULLIFY(id%OOC_NUM_FILE)
              GOTO 111
           ENDIF
           ALLOCATE (id%OOC_POS_IN_FILE(id%KEEP(28)),
     $          stat=ierr)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%KEEP(28)
              NULLIFY(id%OOC_POS_IN_FILE)
              GOTO 111
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(id%KEEP(28)),
     $          stat=ierr)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 111
           ENDIF
           CALL DMUMPS_575(id)
           IF(INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
        CALL DMUMPS_190(0,.FALSE.,DBLE(id%COST_SUBTREES),
     *          KEEP,KEEP8)
#endif
        IF (INFO(1).LT.0) GOTO 111
#if defined(stephinfo)
        write(*,*) 'proc ',id%MYID,' array of dist : ',
     *       id%MEM_DIST(0:id%NSLAVES - 1)
#endif
      END IF
      IF ( ASSOCIATED (id%S) ) THEN
        DEALLOCATE(id%S)
        NULLIFY(id%S)
      ENDIF
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      ALLOCATE (id%S(id%MAXS),stat=ierr)
      IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = id%MAXS
          NULLIFY(id%S)
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
 111  CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP(55) .eq. 0 ) THEN
        IF (ASSOCIATED( id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
        IF ( I_AM_SLAVE .and. KEEP(13) .ne. 0 ) THEN
          ALLOCATE( id%DBLARR( KEEP(13) ), stat = IERR )
        ELSE
          ALLOCATE( id%DBLARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     *       ':Error allocating DBLARR : IERR = ', IERR
          INFO(1)=-13
          INFO(2)=KEEP(13)
          NULLIFY(id%DBLARR)
          GOTO 100
        END IF
      ELSE
         IF ( ASSOCIATED( id%INTARR ) ) THEN
           DEALLOCATE( id%INTARR )
           NULLIFY( id%INTARR )
         END IF
         IF ( I_AM_SLAVE .and. KEEP(14) .ne. 0 ) THEN
           ALLOCATE( id%INTARR( KEEP(14) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = KEEP(14)
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         ELSE
           ALLOCATE( id%INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         END IF
         IF (ASSOCIATED( id%DBLARR)) THEN
           DEALLOCATE(id%DBLARR)
           NULLIFY(id%DBLARR)
         ENDIF
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     *       .AND.   KEEP(46)   .eq. 1
     *       .AND.   KEEP(52)   .eq. 0 ) THEN
             id%DBLARR => id%A_ELT
           ELSE
             IF ( KEEP(13) .ne. 0 ) THEN
               ALLOCATE( id%DBLARR( KEEP(13) ), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = KEEP(13)
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF
             ELSE
               ALLOCATE( id%DBLARR(1), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF 
             END IF
           END IF
         ELSE
           ALLOCATE( id%DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%DBLARR)
             GOTO 100
           END IF
         END IF
      END IF
      IF ( KEEP(38).NE.0 .AND.  I_AM_SLAVE ) THEN
         CALL DMUMPS_165( id%N,
     *   id%root, id%FILS(1), KEEP(38), id%INFO )
      END IF
 100  CONTINUE
      CALL DMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP( 55 ) .eq. 0 ) THEN
      IF ( KEEP(54) .eq. 0 ) THEN
      IF ( id%MYID .eq. MASTER ) THEN
        ALLOCATE(IWK(id%N), stat=allocok)
        IF ( allocok .NE. 0 ) THEN
          INFO(1)=-13
          INFO(2)=id%N
        END IF
      ENDIF
      CALL DMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( id%MYID .eq. MASTER ) THEN
        IF (PROKG ) THEN
          CALL DMUMPS_291(TIME)
        END IF
        IF ( .not. associated( id%INTARR ) ) THEN
          ALLOCATE( id%INTARR( 1 ) )
        ENDIF
#if defined(LARGEMATRICES)
         IF (KEEP(38).NE.0) THEN
          ITOT = NUMROC( id%root%ROOT_SIZE, id%root%MBLOCK,
     *               id%root%MYROW, 0, id%root%NPROW ) 
          ITOT = MAX( 1, ITOT )
          ITOT = ITOT*
     *          NUMROC( id%root%ROOT_SIZE, id%root%NBLOCK,
     *          id%root%MYCOL, 0, id%root%NPCOL )
          ITOT = MAX( 1, ITOT )
         ELSE
          ITOT = 1
         ENDIF
         IF ( ASSOCIATED (id%S) ) THEN
           DEALLOCATE(id%S)
           NULLIFY(id%S)
         ENDIF
         ALLOCATE (WK(ITOT),stat=ierr)
         IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ITOT
          write(6,*) ' PB1 ALLOC LARGEMAT'
          CALL DMUMPS_ABORT()
         ENDIF
           CALL DMUMPS_148(id%N, NZ, id%A(1),
     *      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     *      LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     *      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     *      MIN(KEEP(39),id%NZ),
     *      LP, id%COMM, id%root, KEEP,KEEP8,
     *      id%FILS(1), IWK(1), ! workspace of size N
     *
     *      id%INTARR(1), id%DBLARR(1),
     *      id%PTRAR(1), id%PTRAR(id%N+1),
     *      id%FRERE_STEPS(1), id%STEP(1), WK(1), ITOT,
     *      id%ISTEP_TO_INIV2, id%I_AM_CAND,
     *      id%CANDIDATES) 
        write(6,*) '!!! A,IRN,JCN are freed during facto '
        DEALLOCATE (id%A)
        NULLIFY(id%A)
        DEALLOCATE (id%IRN)
        NULLIFY (id%IRN)
        DEALLOCATE (id%JCN)
        NULLIFY (id%JCN)
        ALLOCATE (id%S(id%MAXS),stat=ierr)
        IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = id%MAXS
          NULLIFY(id%S)
          write(6,*) ' PB2 ALLOC LARGEMAT'
          CALL DMUMPS_ABORT()
        ENDIF
        id%S(id%MAXS-ITOT+1:id%MAXS) = WK(1:ITOT)
        DEALLOCATE (WK)
#else
        CALL DMUMPS_148(id%N, NZ, id%A(1),
     *   id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     *   LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     *   id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     *   MIN(KEEP(39),id%NZ),
     *   LP, id%COMM, id%root, KEEP,KEEP8,
     *   id%FILS(1), IWK(1),
     *
     *   id%INTARR(1), id%DBLARR(1),
     *   id%PTRAR(1), id%PTRAR(id%N+1),
     *   id%FRERE_STEPS(1), id%STEP(1), id%S(1), id%MAXS,
     *   id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     *   id%CANDIDATES(1,1) ) 
#endif
      DEALLOCATE(IWK)
        IF ( PROKG ) THEN
          CALL DMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL DMUMPS_291(TIME)
        END IF
      ELSE
        CALL DMUMPS_145( id%N,
     *       id%DBLARR( 1 ), MAX(1,KEEP( 13 )),
     *       id%INTARR( 1 ), MAX(1,KEEP( 14 )),
     *       id%PTRAR( 1 ),
     *       id%PTRAR(id%N+1),
     *       KEEP( 1 ), KEEP8(1), id%MYID, id%COMM,
     *       MIN(id%KEEP(39),id%NZ),
     *
     *       id%S(1), MAXS,
     *       id%root,
     *       id%PROCNODE_STEPS(1), id%NSLAVES,
     *       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     *       id%INFO(1), id%INFO(2) )
      ENDIF
      ELSE
      IF (PROKG ) THEN
        CALL DMUMPS_291(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
        CALL DMUMPS_282( id%N,
     * id%NZ_loc,
     * id,
     * id%DBLARR(1), KEEP(13), id%INTARR(1),
     * KEEP(14), id%PTRAR(1), id%PTRAR(id%N+1),
     * KEEP(1), KEEP8(1), id%MYID_NODES,
     * id%COMM_NODES, MIN(id%KEEP(39),id%NZ),
     * id%S(1), MAXS, id%root, id%PROCNODE_STEPS,
     * id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     * id%ICNTL(1), id%INFO(1), NSEND, NLOCAL,
     * id%ISTEP_TO_INIV2,
     * id%CANDIDATES )
#if defined(LARGEMATRICES)
         IF (ASSOCIATED(id%IRN_loc)) THEN
            DEALLOCATE(id%IRN_loc)
            NULLIFY(id%IRN_loc)
         ENDIF
         IF (ASSOCIATED(id%JCN_loc)) THEN
            DEALLOCATE(id%JCN_loc)
            NULLIFY(id%JCN_loc)
         ENDIF
         IF (ASSOCIATED(id%A_loc)) THEN
            DEALLOCATE(id%A_loc)
            NULLIFY(id%A_loc)
         ENDIF
       write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
      IF (PROK) THEN
        WRITE(MP,120) NLOCAL, NSEND
      END IF
      END IF
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
        NSEND  = 0
        NLOCAL = 0
      END IF
      CALL MPI_REDUCE( NSEND, NSEND_TOT, 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL, NLOCAL_TOT, 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT, NSEND_TOT
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL DMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL DMUMPS_291(TIME)
      END IF
      END IF
      ELSE
        IF (PROKG ) THEN
          CALL DMUMPS_291(TIME)
        END IF
      IF ( id%MYID.eq.MASTER)
     *CALL DMUMPS_213( id%ELTPTR(1),
     *                        id%NELT,
     *                        MAXELT_SIZE )
      CALL DMUMPS_126( id%N, id%NELT, id%NA_ELT,
     *     id%COMM, id%MYID,
     *     id%NSLAVES, id%PTRAR(1),
     *     id%PTRAR(id%NELT+2),
     *     id%INTARR(1), id%DBLARR(1),
     *     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     *     id%FRTPTR(1), id%FRTELT(1),
     *     id%S(1), id%MAXS, id%FILS(1),
     *     id, id%root )
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( PROKG ) THEN
          CALL DMUMPS_292(TIME)
          WRITE(MPG,160) TIME
          CALL DMUMPS_291(TIME)
      END IF
      END IF
      IF ( I_AM_SLAVE )  THEN
        CALL DMUMPS_528(id%MYID_NODES)
        DMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        DMUMPS_LBUFR_BYTES = MAX( DMUMPS_LBUFR_BYTES,
     *                      100000 )
        PERLU = KEEP( 12 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        DMUMPS_LBUFR_BYTES = DMUMPS_LBUFR_BYTES
     *        + INT( DBLE(MAX(PERLU,MIN_PERLU))*
     *        DBLE(DMUMPS_LBUFR_BYTES)/100D0)
        IF (KEEP(48)==5) THEN
           KEEP(10) = INT( DMUMPS_LBUFR_BYTES / KEEP(35) )
        ENDIF
        IF (KEEP(50).NE.1) THEN
           DMUMPS_LBUFR_BYTES = DMUMPS_LBUFR_BYTES
     *            + INT( DBLE(MAX(PERLU,0))*
     *            DBLE(DMUMPS_LBUFR_BYTES)/100D0)
        ENDIF
        IF (KEEP(48).EQ.5) THEN
          DMUMPS_LBUFR_BYTES=DMUMPS_LBUFR_BYTES+5*KEEP(34)
        ENDIF
        DMUMPS_LBUFR_BYTES = DMUMPS_LBUFR_BYTES +
     *        (6 + id%NSLAVES - 2 + 2 * KEEP(127))
     *        * KEEP(34)
        IF(KEEP(219).NE.0.AND.
     *  KEEP(50) .EQ. 2)
     *  DMUMPS_LBUFR_BYTES=DMUMPS_LBUFR_BYTES+
     *              KEEP(108)*KEEP(35)+KEEP(34)
        IF (id%NSLAVES.EQ.2) THEN
         DMUMPS_LBUF = INT( 1.01D0 * KEEP( 43 ) *KEEP( 35 )  )
        ELSE
         DMUMPS_LBUF = INT( DBLE(KEEP(213)) / 100.0D0 *
     *                     KEEP( 43 ) *KEEP( 35 )  )
        ENDIF
        DMUMPS_LBUF = MAX( DMUMPS_LBUF, 100000 )
        DMUMPS_LBUF = DMUMPS_LBUF
     *                 + INT( 2.0D0 * DBLE(MAX(PERLU,MIN_PERLU))*
     *                   DBLE(DMUMPS_LBUF)/100D0)
        IF(id%KEEP(48).EQ.4)THEN
           DMUMPS_LBUFR_BYTES=DMUMPS_LBUFR_BYTES*5
           DMUMPS_LBUF=DMUMPS_LBUF*5
        ENDIF
        DMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     *               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
          KKKK = DMUMPS_275( id%STEP(KEEP(38)),
     *         id%PROCNODE_STEPS, id%NSLAVES )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             DMUMPS_LBUF_INT = DMUMPS_LBUF_INT + 5 * id%NSLAVES
     *                      * KEEP( 34 )
          END IF
        END IF
        IF ( MP .GT. 0 ) THEN
          WRITE( MP, 9999 ) DMUMPS_LBUFR_BYTES,
     *                      DMUMPS_LBUF, DMUMPS_LBUF_INT
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     *  ' Size of reception buffer in bytes ...... = ', I10,
     *  /,
     *  ' Size of async. emission buffer (bytes).. = ', I10,/,
     *  ' Small emission buffer (bytes) .......... = ', I10)
        CALL DMUMPS_55( DMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     *   ':Error allocating small Send buffer:IERR='
     *   ,IERR
          INFO(1)= -13
          INFO(2)= (DMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
        CALL DMUMPS_53( DMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating Send buffer:IERR='
     *   ,IERR
          INFO(1)= -13
          INFO(2)= (DMUMPS_LBUF+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
        id%LBUFR_BYTES = DMUMPS_LBUFR_BYTES
        id%LBUFR = (DMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
        IF (ASSOCIATED(id%BUFR)) DEALLOCATE(id%BUFR)
        ALLOCATE( id%BUFR( id%LBUFR ),stat=ierr )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating BUFR:IERR='
     *   ,IERR
          INFO(1)=-13
          INFO(2)=id%LBUFR
          NULLIFY(id%BUFR)
          GO TO 110
        END IF
        PERLU          = KEEP( 12 )
        id%MAXIS = MAX( id%MAXIS, 1,
     *       KEEP(15) + 2 * MAX(PERLU,10) * ( KEEP(15) / 100 + 1 )
     *  )
        IF (ASSOCIATED(id%IS)) DEALLOCATE( id%IS )
        ALLOCATE( id%IS( id%MAXIS  ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) id%MYID,':Error allocating IS:IERR=',IERR
         INFO(1)=-13
         INFO(2)=id%MAXIS
         NULLIFY(id%IS)
         GO TO 110
        END IF
        LIW = id%MAXIS
        IF (ASSOCIATED( id%PTLUST_S )) DEALLOCATE(id%PTLUST_S)
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTLUST:IERR = ',
     *    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTLUST_S)
          GOTO 100
        END IF
        IF (ASSOCIATED( id%PTRFAC )) DEALLOCATE(id%PTRFAC)
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTRFAC:IERR = ',
     *    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          GOTO 100
        END IF
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 5 * id%KEEP(28)
        IPOOL  = ITLOC  + id%N
        LPOOL  = DMUMPS_505(id%KEEP(1),id%KEEP8(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     *    IERR
          INFO(1)=-13
          INFO(2)=IPOOL + LPOOL - 1
          GOTO 110
        END IF
      ENDIF
 110  CONTINUE
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( I_AM_SLAVE )  THEN
        CALL DMUMPS_60( id%LBUFR_BYTES )
        IF (MP .GT. 0) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, KEEP(16), KEEP(15), KEEP(13),
     *    KEEP(14), KEEP(25), KEEP(26), KEEP(27)
        ENDIF
      END IF
      PERLU_ON = .TRUE.
      IF (KEEP(201).NE.0) THEN 
        OOC_STAT = .TRUE.
      ELSE
        OOC_STAT = .FALSE.
      ENDIF
      CALL DMUMPS_214( id%KEEP(1), id%KEEP8(1),
     *     id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *     id%NA_ELT,
     *     id%NSLAVES, TOTAL_MBYTES, .FALSE., OOC_STAT, 
     *     PERLU_ON, TOTAL_BYTES)
      id%INFO(16) = TOTAL_MBYTES
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     *    ' ** Space in MBYTES used during factorization  :',
     *                id%INFO(16)
      END IF
      CALL DMUMPS_243( id%MYID, id%COMM,
     *                           id%INFO(16), id%INFOG(18), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** Memory relaxation parameter ( ICNTL(14)  )            :',
     *  KEEP(12)
        WRITE( MPG,'(A,I10) ')
     *  ' ** Rank of processor needing largest memory in facto     :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *  ' ** Space in MBYTES used by this processor for facto      :',
     *  id%INFOG(18)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** Avg. Space in MBYTES per working proc during facto    :',
     *  ( id%INFOG(19)-id%INFO(16) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     *  ' ** Avg. Space in MBYTES per working proc during facto    :',
     *  id%INFOG(19) / id%NSLAVES
        END IF
      END IF
      IF ( I_AM_SLAVE ) THEN
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = id%N
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT = id%NELT
        ELSE
          NELT = 1
        END IF
        CALL DMUMPS_244( id%N, NZ, NSTEPS, id%S(1),
     *      id%MAXS, id%IS( 1 ), LIW,
     *      id%SYM_PERM(1), id%NA(1), id%LNA, id%NE_STEPS(1),
     *      id%ND_STEPS(1), id%FILS(1), id%STEP(1),
     *      id%FRERE_STEPS(1), id%DAD_STEPS(1), id%CANDIDATES(1,1), 
     *      id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     *      id%PTRAR(1), LDPTRAR, IWK( PTRIST ),
     *      id%PTLUST_S(1), id%PTRFAC(1), IWK( PTRWB ),
     *      IWK( ITLOC ), IWK( IPOOL ), LPOOL,
     *      id%CNTL(1), ICNTL, INFO, RINFO, KEEP,KEEP8,
     *      id%PROCNODE_STEPS(1),
     *      id%NSLAVES, id%COMM_NODES,
     *      id%MYID, id%MYID_NODES,
     *      id%BUFR(1),id%LBUFR,id%LBUFR_BYTES,
     *      id%INTARR(1), id%DBLARR(1), id%root,
     *      NELT, id%FRTPTR(1), 
     *      id%FRTELT(1), id%COMM_LOAD, id%ASS_IRECV, SEUIL,
     *      SEUIL_LDLT_NIV2, id%MEM_DIST(0),
     *       id%DKEEP(1),id%PIVNUL_LIST(1),LPN_LIST)
        IF ( MP . GT. 0 .and. KEEP(38) .ne. 0 ) THEN
          WRITE( MP, 175 ) KEEP(49)
        END IF
        DEALLOCATE( IWK )
      ENDIF
        IF ( KEEP(55) .eq. 0 ) THEN
          IF (ASSOCIATED( id%DBLARR)) THEN
            DEALLOCATE(id%DBLARR)
            NULLIFY(id%DBLARR)
          ENDIF
        ELSE
          DEALLOCATE( id%INTARR)
          NULLIFY( id%INTARR )
          IF (      id%MYID_NODES .eq. MASTER
     *      .AND.   KEEP(46)   .eq. 1
     *      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( id%DBLARR )
          ELSE
            IF (ASSOCIATED( id%DBLARR)) THEN
              DEALLOCATE(id%DBLARR)
              NULLIFY(id%DBLARR)
            ENDIF
          END IF
        END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     *                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     *                   id%COMM, IERR )
          END IF
        END IF
      END IF
#if defined(null_space_old)
      IF ( KEEP(19) .EQ. 2 .OR. KEEP(19) .EQ. 4
     *.OR. KEEP(19) .EQ. 6 .OR. KEEP(19) .EQ. 8
     *.OR. KEEP(19) .EQ. 10
     * ) THEN
        CALL DMUMPS_84( id )
      END IF
#endif
      IF (ASSOCIATED(id%BUFR)) THEN
        DEALLOCATE(id%BUFR)
        NULLIFY(id%BUFR)
      END IF
      CALL DMUMPS_57( IERR )
      CALL DMUMPS_59( IERR )
      IF (KEEP(219).NE.0) THEN
      CALL DMUMPS_620()
      ENDIF
      IF ( KEEP(60) .EQ. 1 ) THEN
        ID_SCHUR =DMUMPS_275(id%STEP(id%KEEP(20)),
     *            id%PROCNODE_STEPS(1),id%NSLAVES)
        IF ( id%KEEP( 46 )  .NE. 1 ) THEN
          ID_SCHUR = ID_SCHUR + 1
        END IF
        CALL DMUMPS_276( id%ICNTL, id%INFO,
     *                        id%COMM, id%MYID )
        IF ( id%INFO(1) < 0 ) GOTO 500
        IF ( id%MYID .eq. MASTER .and. ID_SCHUR .eq. MASTER ) THEN
          CALL DCOPY( id%SIZE_SCHUR * id%SIZE_SCHUR,
     *    id%S(id%PTRFAC(id%IS(
     *      id%PTLUST_S(id%STEP(id%KEEP(20)))+4+XSIZE))),
     *    1, id%SCHUR(1), 1 )
        ELSE
          IF ( id%MYID .eq. ID_SCHUR ) THEN
            SIZE_SCHUR =
     *        id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+0+XSIZE)
            CALL MPI_SEND( id%S(
     *          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     *                    +4+XSIZE))),
     *          SIZE_SCHUR*SIZE_SCHUR,
     *          MPI_DOUBLE_PRECISION,
     *          MASTER, TAG_SCHUR,
     *          id%COMM, IERR )
          ELSE IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( id%SCHUR(1), id%SIZE_SCHUR*id%SIZE_SCHUR,
     *                     MPI_DOUBLE_PRECISION, ID_SCHUR, TAG_SCHUR,
     *                     id%COMM, STATUS, IERR )
          END IF
        END IF
      END IF
      IF (KEEP(201). NE. 0) THEN
       IF ( I_AM_SLAVE ) THEN
         CALL DMUMPS_591(id,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 500
         ENDIF
       END IF
       DEALLOCATE(id%S)
       NULLIFY(id%S)
      END IF
      IF ( PROKG ) THEN
        CALL DMUMPS_292(TIME)
        WRITE(MPG,180) TIME
      END IF
      IF (INFO(1).LT.0) GOTO 500
      PERLU_ON = .TRUE.
      IF (KEEP(201).NE.0) THEN 
        OOC_STAT = .TRUE.
      ELSE
        OOC_STAT = .FALSE.
      ENDIF
      CALL DMUMPS_214( id%KEEP(1),id%KEEP8(1),
     *     id%MYID, id%N, id%NELT, id%LNA, id%NZ,
     *     id%NA_ELT,
     *     id%NSLAVES, TOTAL_MBYTES, .TRUE., OOC_STAT, 
     *     PERLU_ON, TOTAL_BYTES)
      IF ( MP .gt. 0 ) THEN
          WRITE(MP,'(A,I10) ')
     *    ' ** Effective minimum Space in BYTES for facto  :',
     *                TOTAL_MBYTES
      ENDIF
      IF (I_AM_SLAVE) THEN
       K67 = KEEP(67)
      ELSE
       K67 = 0
      ENDIF
      CALL DMUMPS_536(PROKG, MPG, K67, id%NSLAVES,
     * id%COMM, "effective space used in S    (KEEP(67)   =")
      CALL DMUMPS_243( id%MYID, id%COMM,
     *                    TOTAL_MBYTES, id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** EFF Min: Rank of processor needing largest memory :',
     *  IRANK
        WRITE( MPG,'(A,I10) ')
     *  ' ** EFF Min: Space in MBYTES used by this processor   :',
     *  id%INFOG(21)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     *  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     *  ( id%INFOG(22)-TOTAL_MBYTES ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     *  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     *  id%INFOG(22) / id%NSLAVES
        END IF
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        IF ( id%KEEP(46) .eq. 0 ) THEN
          RINFO(2)=ZERO
          RINFO(3)=ZERO
          INFO(9:14) = 0
        END IF
      ENDIF
      KEEP(33) = INFO(11)
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     *                 MPI_DOUBLE_PRECISION,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     *                 MPI_DOUBLE_PRECISION,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL DMUMPS_621( INFO(9), INFOG(9), KEEP8(6), id%COMM)
      CALL MPI_REDUCE( INFO(10), INFOG(10), 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_ALLREDUCE( INFO(11), INFOG(11), 1, MPI_INTEGER,
     *                 MPI_MAX, id%COMM, IERR)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( INFO(12), INFOG(12), 3, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(103) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(105) = INFOG(25)
      INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( INFO(25), INFOG(25), 1, MPI_INTEGER,
     *                 MPI_SUM, id%COMM, IERR)
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )
         IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%NPROCS
         END IF
         CALL DMUMPS_276( ICNTL, INFO,
     *     id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 550
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     *        ITMP2(1), 1, MPI_INTEGER, 
     *        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
               CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     *              MPI_INTEGER,I, 
     *              ZERO_PIV, id%COMM, STATUS, IERR)
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     *              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
            ENDDO
         ELSE
            CALL MPI_SEND( id%PIVNUL_LIST(1), KEEP(109), MPI_INTEGER,
     *           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     *           id%COMM, STATUS, IERR )
         ENDIF
         DEALLOCATE(ITMP2)
         INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     *        MPI_SUM, id%COMM, IERR)
      ELSE
         INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
      INFOG(28)=KEEP(112)+KEEP(17)
      IF ( PROKG ) THEN
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP8(6),INFOG(10),
     *                    INFOG(11)
          IF (id%SYM == 1 .OR. id%SYM == 2) THEN
            WRITE(MPG, 99987) INFOG(12)
          END IF
          IF (id%SYM == 0) THEN
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%SYM .NE. 1) THEN
            WRITE(MPG, 99982) INFOG(13)
          END IF
          WRITE(MPG, 99986) KEEP(98)
          IF (id%SYM == 2) THEN
             WRITE(MPG, 99988) KEEP(103)
             WRITE(MPG, 99989) KEEP(105)
          ENDIF
          IF (KEEP(110) .NE.0) THEN
            IF (KEEP(19).EQ.0) THEN
              WRITE(MPG, 99990) KEEP(112)
            ELSE
              WRITE(MPG, 99991) KEEP(112)
            ENDIF
          ENDIF
          IF ( KEEP(17) .ne. 0 )
     *    WRITE(MPG, 99983) KEEP(17)
          IF (KEEP(110).NE.0.OR.KEEP(17).NE.0)
     *    WRITE(MPG, 99992) KEEP(17)+KEEP(112)
          WRITE(MPG, 99981) INFOG(14)
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
      END IF
99981 FORMAT(
     1  ' INFOG(14)  NUMBER OF MEMORY COMPRESS           =',I12)
99982 FORMAT( ' INFOG(12)  NUMBER OF DELAYED PIVOTS      =',I12)
99983 FORMAT( '            DEFICIENCY ON LAST ROOT NODE  =',I12)
99990 FORMAT( '           NUMBER OF ZERO PIVOTS         =',I12)
99991 FORMAT( '           NUMBER OF ZERO PIVOTS WT ROOT =',I12)
99992 FORMAT( ' INFOG(28)  ESTIMATED DEFICIENCY          =',I12)
99984 FORMAT(/' GLOBAL STATISTICS '/
     1  ' RINFOG(2)  OPERATIONS DURING NODE ASSEMBLY     =',1PD10.3/
     2  ' ------(3)  OPERATIONS DURING NODE ELIMINATION  =',1PD10.3/
     3  ' INFOG (9)  REAL SPACE FOR FACTORS              =',I12/
     4  ' INFOG(10)  INTEGER SPACE FOR FACTORS           =',I12/
     5  ' INFOG(11)  MAXIMUM FRONT SIZE                  =',I12)
99985 FORMAT( ' INFOG(13)  NUMBER OF OFF DIAGONAL PIVOTS =',I12)
99986 FORMAT( ' INFOG(25)  NUMBER OF TINY PIVOTS         =',I12)
99987 FORMAT( ' INFOG(12)  NEGATIVE PIVOTS               =',I12)
99988 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 1 nodes     =',I12)
99989 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 2 nodes     =',I12)
 500  CONTINUE
      IF ( I_AM_SLAVE ) THEN
        CALL DMUMPS_183( INFO(1), IERR )
        IF (KEEP(201).NE.0) THEN
          CALL DMUMPS_592(id,IERR)
        ENDIF
      END IF
      CALL DMUMPS_276( ICNTL, INFO,
     *     id%COMM, id%MYID )
 550  CONTINUE
      id%KEEP(13) = KEEP13_SAVE
      RETURN
 120  FORMAT(/' LOCAL REDISTRIB: DATA LOCAL/SENT     =',I12,I12)
 125  FORMAT(/' REDISTRIB: TOTAL DATA LOCAL/SENT     =',I12,I12)
 130  FORMAT(/' ****** FACTORIZATION STEP ********'/)
 160  FORMAT(' GLOBAL TIME FOR MATRIX DISTRIBUTION  =',F12.4)
 170  FORMAT(/' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     *        ' MAXS                                 =',I12/
     *        ' MAXIS                                =',I12/
     *        ' MINIMUM (ICNTL(14)=0) VALUE OF MAXS  =',I12/
     *        ' MINIMUM (ICNTL(14)=0) VALUE OF MAXIS =',I12/
     *        ' REAL SPACE FOR ORIGINAL MATRIX       =',I12/
     *        ' INTEGER SPACE FOR ORIGINAL MATRIX    =',I12/
     *        ' REAL SPACE FOR FACTORS               =',I12/
     *        ' INTEGER SPACE FOR FACTORS            =',I12/
     *        ' MAXIMUM FRONTAL SIZE (ESTIMATED)     =',I12)
 172  FORMAT(/' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     *        ' NUMBER OF WORKING PROCESSES          =',I12/
     *        ' REAL SPACE FOR FACTORS               =',I12/
     *        ' INTEGER SPACE FOR FACTORS            =',I12/
     *        ' MAXIMUM FRONTAL SIZE (ESTIMATED)     =',I12/
     *        ' NUMBER OF NODES IN THE TREE          =',I12)
 175  FORMAT(/' NUMBER OF ENTRIES FOR // ROOT                 =',I12)
 180  FORMAT(/' ELAPSED TIME FOR FACTORIZATION       =',F12.4)
      END SUBROUTINE DMUMPS_142
      SUBROUTINE DMUMPS_536(PROKG, MPG, VAL, NSLAVES,
     *     COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL PROKG
      INTEGER MPG
      INTEGER VAL
      INTEGER NSLAVES
      INTEGER COMM
      CHARACTER*42 MSG 
      INTEGER MAX_VAL
      INTEGER IERR, MASTER
      DOUBLE PRECISION LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
      CALL MPI_REDUCE( VAL, MAX_VAL, 1, MPI_INTEGER,
     *                 MPI_MAX, MASTER, COMM, IERR )
      LOC_VAL = DBLE(VAL)/DBLE(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_DOUBLE_PRECISION,
     *                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        WRITE(MPG,100) " Maximum ", MSG, MAX_VAL
        WRITE(MPG,100) " Average ", MSG, INT(AVG_VAL)
      ENDIF
      RETURN
 100  FORMAT(A9,A42,I12)
      END SUBROUTINE DMUMPS_536
      SUBROUTINE DMUMPS_83
     * ( N, MAPPING, NZ, IRN, JCN, PROCNODE, STEP,
     *   SLAVEF, PERM, FILS,
     *   RG2L, KEEP,KEEP8, MBLOCK, NBLOCK, NPROW, NPCOL )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NZ, SLAVEF, MBLOCK, NBLOCK, NPROW, NPCOL
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN( NZ ), JCN( NZ ) 
      INTEGER MAPPING( NZ ), STEP( N )
      INTEGER PROCNODE( KEEP(28) ), PERM( N ), FILS( N ), RG2L( N )
      INTEGER DMUMPS_275, DMUMPS_330
      EXTERNAL DMUMPS_275, DMUMPS_330
      INTEGER K, IOLD, JOLD, INEW, JNEW, ISEND, JSEND, IARR, INODE
      INTEGER TYPE_NODE, DEST
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INODE = KEEP(38)
      K = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = K
        INODE = FILS( INODE )
        K = K + 1
      END DO
      DO K = 1, NZ
        IOLD = IRN( K )
        JOLD = JCN( K )
        IF ( IOLD .GT. N .OR. IOLD .LT. 1 .OR.
     *       JOLD .GT. N .OR. JOLD .LT. 1 ) THEN
           MAPPING( K ) = -1
           CYCLE
        END IF
        IF ( IOLD .eq. JOLD ) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM( IOLD )
          JNEW = PERM( JOLD )
          IF ( INEW .LT. JNEW ) THEN
            ISEND = IOLD
            IF ( KEEP(50) .ne. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          END IF
        END IF
        IARR = ABS( ISEND )
        TYPE_NODE = DMUMPS_330( ABS(STEP(IARR)),
     *                              PROCNODE, SLAVEF )
        IF ( TYPE_NODE .eq. 1 .or. TYPE_NODE .eq. 2 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = DMUMPS_275( ABS(STEP(IARR)),
     *                             PROCNODE, SLAVEF ) + 1
          ELSE
            DEST = DMUMPS_275( ABS(STEP(IARR)),
     *                             PROCNODE, SLAVEF )
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L( JSEND )
            JPOSROOT = RG2L( IARR  )
          ELSE
            IPOSROOT = RG2L( IARR  )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/MBLOCK, NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/NBLOCK, NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * NPCOL + JCOL_GRID
          END IF
        END IF
        MAPPING( K ) = DEST
      END DO
      RETURN
      END SUBROUTINE DMUMPS_83
      SUBROUTINE DMUMPS_282(
     * N, NZ_loc, id,
     * DBLARR, LDBLARR, INTARR, LINTARR,
     * PTRAIW, PTRARW, KEEP,KEEP8, MYID, COMM, NBRECORDS,
     *
     * A, LA, root, PROCNODE_STEPS, SLAVEF, PERM, STEP,
     * ICNTL, INFO, NSEND, NLOCAL,
     * ISTEP_TO_INIV2, CANDIDATES
     * )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NZ_loc
      TYPE (DMUMPS_STRUC) :: id
      INTEGER LDBLARR, LINTARR
      DOUBLE PRECISION DBLARR( LDBLARR )
      INTEGER INTARR( LINTARR )
      INTEGER PTRAIW( N ), PTRARW( N )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER MYID, COMM, NBRECORDS
      INTEGER LA, SLAVEF
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      INTEGER CANDIDATES(SLAVEF+1, MAX(1,KEEP(56)))
      DOUBLE PRECISION A( LA )
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER PROCNODE_STEPS(KEEP(28)), PERM( N ), STEP( N )
      INTEGER INFO( 40 ), ICNTL(40)
      INTEGER DMUMPS_275, DMUMPS_330, NUMROC
      EXTERNAL DMUMPS_275, DMUMPS_330, NUMROC
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE ), MSGSOU
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      INTEGER END_MSG_2_RECV
      INTEGER I, K, I1, IA
      INTEGER TYPE_NODE, DEST
      INTEGER IOLD, JOLD, IARR, ISEND, JSEND, INEW, JNEW
      INTEGER allocok 
      DOUBLE PRECISION VAL
      INTEGER PTR_ROOT, LOCAL_M, LOCAL_N, ARROW_ROOT
      INTEGER IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT
      INTEGER MP,LP
      INTEGER KPROBE, FREQPROBE
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      INTEGER IACT( SLAVEF ), IREQI( SLAVEF ), IREQR( SLAVEF )
      LOGICAL SEND_ACTIVE( SLAVEF )
      LOGICAL FLAG
      INTEGER NSEND, NLOCAL
      INTEGER MASTER_NODE, ISTEP
      NSEND = 0
      NLOCAL = 0
      LP = ICNTL(1)
      MP = ICNTL(2)
      END_MSG_2_RECV = SLAVEF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
      END IF
      ALLOCATE( BUFR( NBRECORDS, 2, SLAVEF), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *     '** Error allocating real buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( BUFRECR( NBRECORDS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     *    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS
        GOTO 20
      END IF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        INFO(2) = N * 2
      END IF
 20   CONTINUE
      CALL DMUMPS_276( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      ARROW_ROOT = 0
      DO I = 1, N
          I1 = PTRAIW( I )
          IA = PTRARW( I )
          IF ( IA .GT. 0 ) THEN
            DBLARR( IA ) = DBLE(ZERO)
            IW4( I, 1 ) = INTARR( I1 )
            IW4( I, 2 ) = -INTARR( I1 + 1 )
            INTARR( I1 + 2 ) = I
          END IF
      END DO
      IF ( KEEP(38) .NE. 0 ) THEN
          IF (KEEP(60)==0) THEN
          LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *               root%MYROW, 0, root%NPROW )
          LOCAL_M = MAX( 1, LOCAL_M )
          LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
          PTR_ROOT = LA - LOCAL_M * LOCAL_N + 1
          IF ( PTR_ROOT .LE. LA ) THEN
            A( PTR_ROOT:LA ) = DBLE(ZERO)
          END IF
          ELSE
            DO I = 1, root%SCHUR_NLOC
              root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     *        (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC)=DBLE(ZERO)
            ENDDO
          ENDIF
      END IF
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      KPROBE = 0
      FREQPROBE = MAX(1,NBRECORDS/10)
      DO K = 1, NZ_loc
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     *                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     *                 MPI_INTEGER,
     *                 MSGSOU, ARR_INT, COMM, STATUS, IERR )
            CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     *                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
            CALL DMUMPS_102(
     *             BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *             A, LA,
     *             END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *             ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *             INTARR, LINTARR, DBLARR, LDBLARR
     *             )
          END IF
        END IF
        IOLD = id%IRN_loc(K)
        JOLD = id%JCN_loc(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) CYCLE
        VAL = id%A_loc(K)
        IF (IOLD.EQ.JOLD) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          IF (INEW.LT.JNEW) THEN
            ISEND = IOLD
            IF ( KEEP(50) .NE. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
        ENDIF
        IARR = ABS( ISEND )
        ISTEP = ABS(STEP(IARR))
        TYPE_NODE = DMUMPS_330( ISTEP,
     *              PROCNODE_STEPS, SLAVEF )
        MASTER_NODE= DMUMPS_275( ISTEP,
     *           PROCNODE_STEPS, SLAVEF )
        IF ( TYPE_NODE .eq. 1 ) THEN
          DEST = MASTER_NODE
        ELSE IF ( TYPE_NODE .eq. 2 ) THEN
          IF ( ISEND .LT. 0 ) THEN
            DEST = -1
          ELSE
            DEST = MASTER_NODE
          END IF
        ELSE
          IF ( ISEND < 0 ) THEN
            IPOSROOT = root%RG2L_ROW(JSEND)
            JPOSROOT = root%RG2L_ROW(IARR )
          ELSE
            IPOSROOT = root%RG2L_ROW(IARR )
            JPOSROOT = root%RG2L_ROW(JSEND)
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          DEST = IROW_GRID * root%NPCOL + JCOL_GRID
        END IF
        if (dest .eq. -1) then
          nlocal = nlocal + 1
          nsend = nsend + slavef -1
        else
          if (dest .eq.myid ) then
            nlocal = nlocal + 1
          else
            nsend = nsend + 1
          endif
        end if
        IF ( DEST.EQ.-1) THEN
         DO I=1, CANDIDATES(SLAVEF+1,ISTEP_TO_INIV2(ISTEP))
            DEST=CANDIDATES(I,ISTEP_TO_INIV2(ISTEP))
            CALL DMUMPS_101( DEST, ISEND, JSEND, VAL,
     *   BUFI, BUFR, BUFRECI, BUFRECR,
     *   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     *   root, KEEP,KEEP8 )
         ENDDO
         DEST=MASTER_NODE
         CALL DMUMPS_101( DEST, ISEND, JSEND, VAL,
     *   BUFI, BUFR, BUFRECI, BUFRECR,
     *   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     *   root, KEEP,KEEP8 )
        ELSE
         CALL DMUMPS_101( DEST, ISEND, JSEND, VAL,
     *   BUFI, BUFR, BUFRECI, BUFRECR,
     *   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     *   root, KEEP,KEEP8 )
        ENDIF
      END DO
      DEST = -2
        CALL DMUMPS_101( DEST, ISEND, JSEND, VAL,
     *  BUFI, BUFR, BUFRECI, BUFRECR,
     *  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     *  N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, 
     *  IW4(1,1), root, KEEP,KEEP8 )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     *                 MPI_ANY_SOURCE, ARR_INT, COMM, STATUS, IERR )
        MSGSOU = STATUS( MPI_SOURCE )
        CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     *                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
        CALL DMUMPS_102(
     *           BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *           KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *           A, LA,
     *           END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *           ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *           INTARR, LINTARR, DBLARR, LDBLARR
     *           )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
          CALL MPI_WAIT( IREQR( I ), STATUS, IERR )
        END IF
      END DO
      KEEP(49) = ARROW_ROOT
      DEALLOCATE( IW4 )
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( BUFRECI )
      DEALLOCATE( BUFRECR )
      RETURN
      END SUBROUTINE DMUMPS_282
      SUBROUTINE DMUMPS_192
     *( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     *  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      DOUBLE PRECISION A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      Y_loc = DBLE(ZERO)
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *          (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     *        .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *        (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_192
      SUBROUTINE DMUMPS_101( DEST, ISEND, JSEND, VAL,
     *  BUFI, BUFR, BUFRECI, BUFRECR,
     *  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     *  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR, N,
     *  PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     *  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4, root,
     *  KEEP,KEEP8 )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ISEND, JSEND, DEST, NBRECORDS, SLAVEF, COMM, MYID, N
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER ARROW_ROOT, END_MSG_2_RECV, LOCAL_M, LOCAL_N
      INTEGER LINTARR, LDBLARR, LA, PTR_ROOT
      INTEGER BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER BUFRECI( NBRECORDS * 2 + 1 )
      INTEGER IREQI(SLAVEF), IREQR(SLAVEF), IACT(SLAVEF)
      INTEGER IW4( N, 2 )
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR( LINTARR )
      DOUBLE PRECISION DBLARR( LDBLARR ), A( LA )
      LOGICAL SEND_ACTIVE(SLAVEF)
      DOUBLE PRECISION BUFR( NBRECORDS, 2, SLAVEF )
      DOUBLE PRECISION BUFRECR( NBRECORDS )
      DOUBLE PRECISION VAL
      INTEGER ISLAVE, IBEG, IEND, NBREC, IREQ
      INTEGER TAILLE_SEND_I, TAILLE_SEND_R, MSGSOU
      LOGICAL FLAG, SEND_LOCAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS(MPI_STATUS_SIZE)
      IF ( DEST .eq. -2 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      SEND_LOCAL = .FALSE.
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -2 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -2 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQR( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     *                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     *                  MPI_INTEGER, MSGSOU, ARR_INT, COMM,
     *                  STATUS, IERR )
                  CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     *                  MPI_DOUBLE_PRECISION, MSGSOU,
     *                  ARR_REAL, COMM, STATUS, IERR )
                  CALL DMUMPS_102(
     *              BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     *              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *              A, LA,
     *              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *              INTARR, LINTARR, DBLARR, LDBLARR
     *              )
                END IF
            ELSE
                CALL MPI_WAIT( IREQI( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            TAILLE_SEND_R = NBREC
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     *          TAILLE_SEND_I,
     *          MPI_INTEGER, ISLAVE - 1, ARR_INT, COMM,
     *          IREQI( ISLAVE ), IERR )
            CALL MPI_ISEND( BUFR(1, IACT(ISLAVE), ISLAVE ),
     *          TAILLE_SEND_R,
     *          MPI_DOUBLE_PRECISION, ISLAVE - 1, ARR_REAL, COMM,
     *          IREQR( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            SEND_LOCAL = .TRUE.
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
          BUFR(IREQ,IACT(ISLAVE),ISLAVE )    = VAL
        END IF
      END DO
      IF ( SEND_LOCAL ) THEN
            ISLAVE = MYID + 1
            CALL DMUMPS_102(
     *              BUFI(1,3-IACT(ISLAVE),ISLAVE),
     *              BUFR(1,3-IACT(ISLAVE),ISLAVE),
     *              NBRECORDS, N, IW4(1,1),
     *              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     *              A, LA,
     *              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     *              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     *              INTARR, LINTARR, DBLARR, LDBLARR
     *              )
      END IF
      RETURN
      END SUBROUTINE DMUMPS_101
      SUBROUTINE DMUMPS_102
     *           ( BUFI, BUFR, NBRECORDS, N, IW4,
     *             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT, A, LA,
     *             END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     *             SLAVEF, ARROW_ROOT,
     *             PTRAIW, PTRARW, PERM, STEP,
     *             INTARR, LINTARR, DBLARR, LDBLARR )
      IMPLICIT NONE
      INCLUDE 'dmumps_root.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER NBRECORDS, N, ARROW_ROOT, MYID, SLAVEF
      INTEGER BUFI( NBRECORDS * 2 + 1 )
      DOUBLE PRECISION BUFR( NBRECORDS )
      INTEGER IW4( N, 2 )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER END_MSG_2_RECV
      INTEGER PTRAIW( N ), PTRARW( N ), PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LINTARR, LDBLARR
      INTEGER INTARR( LINTARR )
      INTEGER LOCAL_M, LOCAL_N, PTR_ROOT, LA
      DOUBLE PRECISION A( LA ), DBLARR( LDBLARR )
      INTEGER DMUMPS_330, DMUMPS_275
      EXTERNAL DMUMPS_330, DMUMPS_275
      INTEGER IREC, NB_REC, NODE_TYPE, IPROC
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID,
     *        ILOCROOT, JLOCROOT
      INTEGER IA, IS1, ISHIFT, IIW, IS, IAS, IARR, JARR
      INTEGER TAILLE
      DOUBLE PRECISION VAL
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        NODE_TYPE = DMUMPS_330( ABS(STEP(ABS( IARR ))),
     *              PROCNODE_STEPS, SLAVEF )
        IF ( NODE_TYPE .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
          IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( IROW_GRID .NE. root%MYROW .OR.
     *       JCOL_GRID .NE. root%MYCOL ) THEN
            WRITE(*,*) MYID,':INTERNAL Error: recvd root arrowhead '
            WRITE(*,*) MYID,':not belonging to me. IARR,JARR=',IARR,JARR
            WRITE(*,*) MYID,':IROW_GRID,JCOL_GRID=',IROW_GRID,JCOL_GRID
            WRITE(*,*) MYID,':MYROW, MYCOL=', root%MYROW, root%MYCOL
            WRITE(*,*) MYID,':IPOSROOT,JPOSROOT=', IPOSROOT, JPOSROOT
            CALL DMUMPS_ABORT()
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
          A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *      + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *      * LOCAL_M + ILOCROOT - 1 ) + VAL
          ELSE
          root%SCHUR_POINTER(( JLOCROOT - 1 ) * root%SCHUR_LLD
     *      + ILOCROOT ) = root%SCHUR_POINTER( (JLOCROOT - 1)
     *      * root%SCHUR_LLD + ILOCROOT) + VAL
          ENDIF
        ELSE IF (IARR.GE.0) THEN
         IF (IARR.EQ.JARR) THEN
          IA = PTRARW(IARR)
          DBLARR(IA) = DBLARR(IA) + VAL
         ELSE
          IS1 =  PTRAIW(IARR)
          ISHIFT      = INTARR(IS1) + IW4(IARR,2)
          IW4(IARR,2) = IW4(IARR,2) - 1
          IIW         = IS1 + ISHIFT + 2
          INTARR(IIW)     = JARR
          IS          = PTRARW(IARR)
          IAS         = IS + ISHIFT
          DBLARR(IAS) = VAL
         ENDIF
        ELSE
           IARR = -IARR
           ISHIFT      = PTRAIW(IARR)+IW4(IARR,1)+2
           INTARR(ISHIFT)  = JARR
           IAS         = PTRARW(IARR)+IW4(IARR,1)
           IW4(IARR,1) = IW4(IARR,1) - 1
           DBLARR(IAS)      = VAL
           IPROC = DMUMPS_275( ABS(STEP(IARR)),
     *             PROCNODE_STEPS,SLAVEF )
           IF ( KEEP(50) .NE. 0 .AND.
     *          IW4(IARR,1) .EQ. 0 .AND.
     *          IPROC .EQ. MYID
     *          .AND. STEP(IARR) > 0 ) THEN
             TAILLE = INTARR( PTRAIW(IARR) )
             CALL DMUMPS_310( N, PERM,
     *            INTARR( PTRAIW(IARR) + 3 ),
     *            DBLARR( PTRARW(IARR) + 1 ),
     *            TAILLE, MYID )
           END IF
        ENDIF
      ENDDO
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_102
      SUBROUTINE DMUMPS_193
     *( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     *  LDLT, MTYPE)
      IMPLICIT NONE
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      DOUBLE PRECISION A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
      INTEGER I, J, K
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      Y_loc = DBLE(ZERO)
      IF ( LDLT .eq. 0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *          (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + ABS( A_loc(K) * X(J) )
        ENDDO
        ELSE
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     *        .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(J) = Y_loc(J) + ABS( A_loc(K) * X(I) )
        ENDDO
        END IF
      ELSE
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     *        (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y_loc(I) = Y_loc(I) + ABS( A_loc(K) * X(J) )
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + ABS( A_loc(K) * X(I) )
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_193
      SUBROUTINE DMUMPS_151( NRHS, N, KEEP28, IWCB, LIWW,
     *       W, LWC,
     *       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER NRHS, N,LIWW,LWC,POSWCB,IWPOSCB, KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      DOUBLE PRECISION W(LWC)
      INTEGER SIZFI, SIZFR
      IF ( IWPOSCB .eq. LIWW ) RETURN
      DO WHILE ( IWCB( IWPOSCB + 2 ) .eq. 0 )
        SIZFR = IWCB( IWPOSCB + 1 )
        SIZFI = SIZFR + 2
        SIZFR = SIZFR * NRHS
        IWPOSCB = IWPOSCB + SIZFI
        POSWCB  = POSWCB  + SIZFR
        IF ( IWPOSCB .eq. LIWW ) RETURN
      END DO
      RETURN
      END SUBROUTINE DMUMPS_151
      SUBROUTINE DMUMPS_95(NRHS,N,KEEP28,IWCB,LIWW,W,LWC,
     *       POSWCB,IWPOSCB,PTRICB,PTRACB)
      IMPLICIT NONE
      INTEGER NRHS,N,LIWW,LWC,POSWCB,IWPOSCB,KEEP28
      INTEGER IWCB(LIWW),PTRICB(KEEP28),PTRACB(KEEP28)
      DOUBLE PRECISION W(LWC)
      INTEGER IPTIW,IPTA,SIZFI,SIZFR,LONGI,LONGR
      INTEGER I
      IPTIW = IWPOSCB
      IPTA  = POSWCB
      LONGI = 0
      LONGR = 0
      IF ( IPTIW .EQ. LIWW ) RETURN
10    CONTINUE
      IF (IWCB(IPTIW+2).EQ.0) THEN
        SIZFR  = IWCB(IPTIW+1)
        SIZFI  = SIZFR + 2
        SIZFR  = SIZFR * NRHS
        IF (LONGI.NE.0) THEN
          DO 20 I=0,LONGI-1
            IWCB(IPTIW + SIZFI - I) = IWCB (IPTIW - I )
 20       CONTINUE 
          DO 30 I=0,LONGR-1
            W(IPTA + SIZFR - I)   = W(IPTA - I )
 30       CONTINUE
        ENDIF
        DO 40 I=1,KEEP28
          IF ((PTRICB(I).LE.(IPTIW+1)).AND.
     &        (PTRICB(I).GT.IWPOSCB) ) THEN
            PTRICB(I) = PTRICB(I) + SIZFI
            PTRACB(I) = PTRACB(I) + SIZFR
          ENDIF 
40      CONTINUE 
        IWPOSCB = IWPOSCB + SIZFI
        IPTIW   = IPTIW + SIZFI
        POSWCB = POSWCB + SIZFR
        IPTA   = IPTA + SIZFR     
       ELSE
        SIZFR  = IWCB(IPTIW+1)
        SIZFI  = SIZFR + 2
        SIZFR  = SIZFR * NRHS
        IPTIW = IPTIW + SIZFI
        LONGI = LONGI + SIZFI
        IPTA  = IPTA + SIZFR
        LONGR = LONGR + SIZFR
       ENDIF
       IF (IPTIW.NE.LIWW) GOTO 10
       RETURN
       END SUBROUTINE DMUMPS_95
      SUBROUTINE DMUMPS_205(MTYPE, IFLAG, N, NZ,
     *    LHS, WRHS, W, RHS, GIVSOL, SOL, ANORM, XNORM, SCLNRM,
     *    MPRINT, ICNTL, KEEP,KEEP8)
      INTEGER MTYPE,N,NZ,IFLAG,ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION RHS(N),LHS(N)
      DOUBLE PRECISION WRHS(N),SOL(*)
      DOUBLE PRECISION ERL2,ERREL,SCLNRM,W(N)
      DOUBLE PRECISION RESMAX,RESL2,XNORM, ERMAX,MAXSOL,
     *     COMAX
      DOUBLE PRECISION ANORM,DZERO,EPSI
      LOGICAL GIVSOL,PROK
      INTEGER MPRINT, MP
      INTEGER I,J,K
      INTRINSIC ABS, MAX, SQRT
      MP = ICNTL(2)
      PROK = (MPRINT .GT. 0)
      DZERO = 0.0D0
      EPSI = 0.1D-9
      ANORM = DZERO
      RESMAX = DZERO
      RESL2 = DZERO
      DO 40 K = 1, N
        RESMAX = MAX(RESMAX, ABS(RHS(K)))
        RESL2 = RESL2 + ABS(RHS(K)) * ABS(RHS(K))
        ANORM = MAX(ANORM, ABS(W(K)))
   40 CONTINUE
      XNORM = DZERO
      DO 50 K = 1, N
        XNORM = MAX(XNORM, ABS(LHS(K)))
   50 CONTINUE
      IF (XNORM .GT. EPSI) THEN
        SCLNRM = RESMAX / (ANORM * XNORM)
      ELSE
        IFLAG = IFLAG + 2
        IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     *' MAX-NORM of computed solut. is zero'
        SCLNRM = RESMAX / ANORM
      ENDIF
      RESL2 = SQRT(RESL2)
      ERMAX = DZERO
      COMAX = DZERO
      ERL2 = DZERO
      IF (.NOT.GIVSOL) THEN
        IF (PROK) WRITE( MPRINT, 90 ) RESMAX, RESL2, ANORM, XNORM, 
     *      SCLNRM
      ELSE
        MAXSOL = DZERO
        DO 60 K = 1, N
          MAXSOL = MAX(MAXSOL, ABS(SOL(K)))
   60   CONTINUE
        DO 70 K = 1, N
          ERL2 = (LHS(K) - SOL(K)) ** 2 + ERL2
          ERMAX = MAX(ERMAX, ABS(LHS(K) - SOL(K)))
   70   CONTINUE
        DO 80 K = 1, N
          IF (ABS(SOL(K)) .GT. EPSI) THEN
            COMAX = MAX(COMAX, (ABS(LHS(K) - SOL(K)) / ABS(SOL(K))))
          ENDIF
   80   CONTINUE
        ERL2 = SQRT(ERL2)
        IF (MAXSOL .GT. EPSI) THEN
          ERREL = ERMAX / MAXSOL
        ELSE
          IFLAG = IFLAG + 2
          IF ((MP .GT. 0) .AND. (ICNTL(4) .GE. 2)) WRITE( MP, * ) 
     *' MAX-NORM of exact solution is zero'
          ERREL = ERMAX
        ENDIF
        IF (PROK) WRITE( MPRINT, 100 ) ERMAX, ERL2, ERREL, COMAX, RESMAX
     *      , RESL2, ANORM, XNORM, SCLNRM
      ENDIF
   90  FORMAT (/' RESIDUAL IS ............ (MAX-NORM)        =',1PD9.2/
     *       '                       .. (2-NORM)          =',1PD9.2/
     *       ' RINFOG(4):NORM OF input  Matrix  (MAX-NORM)=',1PD9.2/
     *       ' RINFOG(5):NORM OF Computed SOLUT (MAX-NORM)=',1PD9.2/
     *       ' RINFOG(6):SCALED RESIDUAL ...... (MAX-NORM)=',1PD9.2)
      RETURN
  100  FORMAT (/' ERROR IS     ............ (MAX-NORM)       =',1PD9.2/
     *       '              ............ (2-NORM)         =',1PD9.2/
     *       ' RELATIVE ERROR........... (MAX-NORM)       =',1PD9.2/
     *       ' Comp. Wise ERROR......... (MAX-NORM)       =',1PD9.2/
     *       ' AND RESIDUAL IS ......... (MAX-NORM)       =',1PD9.2/
     *       '                        .. (2-NORM)         =',1PD9.2/
     *       ' NORM OF input  MATRIX ... (MAX-NORM)       =',1PD9.2/
     *       ' NORM of computed SOLUT... (MAX-NORM)       =',1PD9.2/
     *       ' SCALED RESIDUAL ......... (MAX-NORM)       =',1PD9.2)
      END SUBROUTINE DMUMPS_205
      SUBROUTINE DMUMPS_206(NZ, N, RHS,
     *    X, Y, D, W, IW, KASE,
     *    OMEGA, ERX, JOB, COND, MAXIT, NOITER, LP, KEEP,KEEP8,
     *    ARRET )
      INTEGER NZ, N, KASE, KEEP(500), JOB
      INTEGER*8 KEEP8(150)
      INTEGER IW(N,2)
      DOUBLE PRECISION RHS(N)
      DOUBLE PRECISION X(N), Y(N), D(N), W(N,3)
      INTEGER LP, MAXIT, NOITER
      DOUBLE PRECISION COND(2),OMEGA(2)
      DOUBLE PRECISION ARRET
      DOUBLE PRECISION CGCE, CTAU
      DATA  CTAU /1.0D3/, CGCE /0.2D0/
      LOGICAL LCOND1, LCOND2
      INTEGER IFLAG, JUMP, I, IMAX
      DOUBLE PRECISION ERX, DXMAX
      DOUBLE PRECISION CONVER, OM1, OM2, DXIMAX
      DOUBLE PRECISION  ZERO, ONE,TAU, DD
      DOUBLE PRECISION OLDOMG(2)
      INTEGER IDAMAX
      INTRINSIC     ABS, MAX
      SAVE LCOND1, LCOND2, JUMP,  DXIMAX, DXMAX, CONVER,
     *     OM1, OLDOMG, IFLAG
      DATA ZERO /0.0D0/, ONE /1.0D0/
      IF (KASE .EQ. 0) THEN
        LCOND1 = .FALSE.
        LCOND2 = .FALSE.
        COND(1) = ONE
        COND(2) = ONE
        ERX = DBLE(ZERO)
        OM1 = ZERO
        IFLAG = 0
        NOITER = 0
        JUMP = 1
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 30
      CASE(2)
        GOTO 10
      CASE(3)
        GOTO 110
      CASE(4)
        GOTO 150
      CASE(5)
        GOTO 35
      CASE DEFAULT
      END SELECT
   10 CONTINUE
      DO 20 I = 1, N
        X(I) = X(I) + Y(I)
   20 CONTINUE
      IF (NOITER .GT. MAXIT) THEN
        IFLAG = IFLAG + 8
        GOTO 70
      ENDIF
   30 CONTINUE
      KASE = 14
      JUMP = 5
      RETURN
   35 CONTINUE
      IMAX = IDAMAX(N, X, 1)
      DXMAX = ABS(X(IMAX))
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      DO 40 I = 1, N
        TAU = (W(I, 2) * DXMAX + ABS(RHS(I))) * N * CTAU
        DD = W(I, 1) + ABS(RHS(I))
        IF ((DD + TAU) .GT. TAU) THEN
          OMEGA(1) = MAX(OMEGA(1), ABS(Y(I) / DD))
          IW(I, 1) = 1
        ELSE
          IF (TAU .GT. ZERO) THEN
            OMEGA(2) = MAX(OMEGA(2), ABS(Y(I) / (DD + W(I, 2) * DXMAX)))
          ENDIF
          IW(I, 1) = 2
        ENDIF
   40 CONTINUE
      OM2 = OMEGA(1) + OMEGA(2)
      IF (OM2 .LT. ARRET ) GOTO 70
      IF (MAXIT .EQ. 0) GOTO 70
      IF (NOITER .GT. 1 .AND. OM2 .GT. OM1 * CGCE) THEN
        CONVER = OM2 / OM1
        IF (OM2 .GT. OM1) THEN
          OMEGA(1) = OLDOMG(1)
          OMEGA(2) = OLDOMG(2)
          DO 50 I = 1, N
            X(I) = W(I, 3)
   50     CONTINUE
        ENDIF
        GOTO 70
      ENDIF
      DO 60 I = 1, N
        W(I, 3) = X(I)
   60 CONTINUE
      OLDOMG(1) = OMEGA(1)
      OLDOMG(2) = OMEGA(2)
      OM1 = OM2
      NOITER = NOITER + 1
      KASE = 2
      JUMP = 2
      RETURN
   70 KASE = 0
      IF (JOB .LE. 0) GOTO 170
      DO 80 I = 1, N
        IF (IW(I, 1) .EQ. 1) THEN
          W(I, 1) = W(I, 1) + ABS(RHS(I))
          W(I, 2) = DBLE(ZERO)
          LCOND1 = .TRUE.
        ELSE
          W(I, 2) = W(I, 2) * DXMAX + W(I, 1)
          W(I, 1) = DBLE(ZERO)
          LCOND2 = .TRUE.
        ENDIF
   80 CONTINUE
      DO 90 I = 1, N
        W(I, 3) = X(I) * D(I)
   90 CONTINUE
      IMAX = IDAMAX(N, W(1, 3), 1)
      DXIMAX = ABS(W(IMAX, 3))
      IF (.NOT.LCOND1) GOTO 130
  100 CALL DMUMPS_218(N, KASE, Y, COND(1), W(1, 3), IW(1, 2))
      IF (KASE .EQ. 0) GOTO 120
      IF (KASE .EQ. 1) CALL DMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL DMUMPS_204(N, Y, W)
      JUMP = 3
      RETURN
  110 CONTINUE
      IF (KASE .EQ. 1) CALL DMUMPS_204(N, Y, W)
      IF (KASE .EQ. 2) CALL DMUMPS_204(N, Y, D)
      GOTO 100
  120 IF (DXIMAX .GT. ZERO) COND(1) = COND(1) / DXIMAX
      ERX = OMEGA(1) * COND(1)
  130 IF (.NOT.LCOND2) GOTO 170
      KASE = 0
  140 CALL DMUMPS_218(N, KASE, Y, COND(2), W(1, 3), IW(1, 2))
      IF (KASE .EQ. 0) GOTO 160
      IF (KASE .EQ. 1) CALL DMUMPS_204(N, Y, D)
      IF (KASE .EQ. 2) CALL DMUMPS_204(N, Y, W(1, 2))
      JUMP = 4
      RETURN
  150 CONTINUE
      IF (KASE .EQ. 1) CALL DMUMPS_204(N, Y, W(1, 2))
      IF (KASE .EQ. 2) CALL DMUMPS_204(N, Y, D)
      GOTO 140
  160 IF (DXIMAX .GT. ZERO) THEN
        COND(2) = COND(2) / DXIMAX
      ENDIF
      ERX = ERX + OMEGA(2) * COND(2)
  170 KASE = -IFLAG
      RETURN
      END SUBROUTINE DMUMPS_206
      SUBROUTINE DMUMPS_207(A, NZ, N, IRN, ICN, Z, KEEP,KEEP8)
      INTEGER NZ, N, I, J, K, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION A(NZ), Z(N)
      DOUBLE PRECISION ZERO
      INTRINSIC     ABS
      DATA ZERO /0.0D0/
      DO 10 I = 1, N
        Z(I) = DBLE(ZERO)
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K))
        IF (J.NE.I) THEN 
          Z(J) = Z(J) + ABS(A(K))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_207
      SUBROUTINE DMUMPS_289(A, NZ, N, IRN, ICN, Z,
     *            KEEP, KEEP8, COLSCA)
      INTEGER NZ, N, I, J, K, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION A(NZ), Z(N), COLSCA(N)
      DOUBLE PRECISION  ZERO
      INTRINSIC     ABS
      DATA ZERO /0.0D0/
      DO 10 I = 1, N
        Z(I) = DBLE(ZERO)
   10 CONTINUE
      IF (KEEP(50) .EQ.0) THEN
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K)*COLSCA(J))
       ENDDO
      ELSE
       DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .LT. 1) .OR. (I .GT. N)) CYCLE
        IF ((J .LT. 1) .OR. (J .GT. N)) CYCLE
        Z(I) = Z(I) + ABS(A(K)*COLSCA(J))
        IF (J.NE.I) THEN
          Z(J) = Z(J) + ABS(A(K)*COLSCA(I))
        ENDIF
       ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_289
      SUBROUTINE DMUMPS_208(A, NZ, N, IRN, ICN, RHS, X, R, W,
     *           KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER NZ, N, I, K, J, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION A(NZ), RHS(N), X(N), R(N), W(N),
     *     D
      DOUBLE PRECISION ZERO
      INTRINSIC    ABS
      DATA ZERO /0.0D0/
      DO I = 1, N
        R(I) = RHS(I)
        W(I) = DBLE(ZERO)
      ENDDO
      DO K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I .GT. N) .OR. (J .GT. N) .OR. (I .LT. 1) .OR. (J .LT. 1))
     *      CYCLE
        D = A(K) * X(J)
        R(I) = R(I) - D
        W(I) = W(I) + ABS(D)
        IF ((I.NE.J) .AND. (KEEP(50).NE.0) ) THEN
          D = A(K) * X(I)
          R(J) = R(J) - D
          W(J) = W(J) + ABS(D)
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_208
      SUBROUTINE DMUMPS_204(N, R, W)
      INTEGER N, I
      DOUBLE PRECISION  R(N), W(N)
      DO 10 I = 1, N
        R(I) = R(I) * W(I)
   10 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_204
      SUBROUTINE DMUMPS_218(N, KASE, X, EST, W, IW)
      INTEGER N, IW(N), KASE
      DOUBLE PRECISION W(N), X(N)
      DOUBLE PRECISION EST
      INTRINSIC ABS, NINT, REAL, SIGN
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN
      DOUBLE PRECISION ZERO, ONE, TEMP
      SAVE ITER, J, JLAST, JUMP
      DATA ZERO /0.0D0/
      DATA ONE  /1.0D0/
      IF (KASE .EQ. 0) THEN
        DO 10 I = 1, N
          X(I) = ONE / DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        RETURN
      ENDIF
      SELECT CASE (JUMP)
      CASE (1)
        GOTO 20
      CASE(2)
        GOTO 40
      CASE(3)
        GOTO 70
      CASE(4)
        GOTO 120
      CASE(5)
        GOTO 160
      CASE DEFAULT
      END SELECT
   20 CONTINUE
      IF (N .EQ. 1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GOTO 190
      ENDIF
      DO 30 I = 1, N
        X(I) = SIGN(ONE,ABS( X(I)))
        IW(I) = NINT(ABS(X(I)))
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
   40 CONTINUE
      J = IDAMAX(N, X, 1)
      ITER = 2
   50 CONTINUE
      DO 60 I = 1, N
        X(I) = DBLE(ZERO)
   60 CONTINUE
      X(J) = DBLE(ONE)
      KASE = 1
      JUMP = 3
      RETURN
   70 CONTINUE
      DO 80 I = 1, N
        W(I) = X(I)
   80 CONTINUE
      DO 90 I = 1, N
        IF (NINT(SIGN(ONE, ABS(X(I)))) .NE. IW(I)) GOTO 100
   90 CONTINUE
      GOTO 130
  100 CONTINUE
      DO 110 I = 1, N
        X(I) = SIGN(ONE, ABS(X(I)))
        IW(I) = NINT(ABS(X(I)))
  110 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  120 CONTINUE
      JLAST = J
      J = IDAMAX(N, X, 1)
      IF ((ABS(X(JLAST)) .NE. ABS(X(J))) .AND. (ITER .LT. ITMAX)) THEN
        ITER = ITER + 1
        GOTO 50
      ENDIF
  130 CONTINUE
      EST = DBLE(ZERO)
      DO 140 I = 1, N
        EST = EST + ABS(W(I))
  140 CONTINUE
      ALTSGN = ONE
      DO 150 I = 1, N
        X(I) = ALTSGN * (ONE + DBLE(I - 1) / DBLE(N - 1))
        ALTSGN = -ALTSGN
  150 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  160 CONTINUE
      TEMP = ZERO
      DO 170 I = 1, N
        TEMP = TEMP + ABS(X(I))
  170 CONTINUE
      TEMP = 2.0 * TEMP / DBLE(3 * N)
      IF (TEMP .GT. EST) THEN
        DO 180 I = 1, N
          W(I) = X(I)
  180   CONTINUE
        EST = TEMP
      ENDIF
  190 KASE = 0
      RETURN
      END SUBROUTINE DMUMPS_218
      SUBROUTINE DMUMPS_278( MTYPE, N, NZ, ASPK, IRN, ICN,
     *    LHS, WRHS, W, RHS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MTYPE, N, NZ
      INTEGER IRN( NZ ), ICN( NZ )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION ASPK( NZ )
      DOUBLE PRECISION LHS( N ), WRHS( N ), W( N ), RHS( N )
      INTEGER K, I, J
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO = 0.0D0)
      DO 10 K = 1, N
        W(K) = DZERO
        RHS(K) = WRHS(K)
   10 CONTINUE
      IF ( KEEP(50) .EQ. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + ABS(ASPK(K))
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(J) = RHS(J) - ASPK(K) * LHS(I)
          W(J) = W(J) + ABS(ASPK(K))
        ENDDO
       ENDIF
      ELSE
       DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          RHS(I) = RHS(I) - ASPK(K) * LHS(J)
          W(I) = W(I) + ABS(ASPK(K))
          IF (J.NE.I) THEN
            RHS(J) = RHS(J) - ASPK(K) * LHS(I)
            W(J) = W(J) + ABS(ASPK(K))
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_278
      SUBROUTINE DMUMPS_121( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    LHS, WRHS, W, RHS, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION A_ELT(NA_ELT)
      DOUBLE PRECISION LHS( N ), WRHS( N ), W( N ), RHS( N )
      CALL DMUMPS_257(N, NELT, ELTPTR, ELTVAR, A_ELT,
     *                         LHS, RHS, KEEP(50), MTYPE )
      RHS = WRHS - RHS
      CALL DMUMPS_119( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP,KEEP8 )
      RETURN
      END SUBROUTINE DMUMPS_121
      SUBROUTINE DMUMPS_119( MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION A_ELT(NA_ELT), W(N), TEMP
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO = 0.0D0)
      W = DBLE(DZERO)
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) = 
     &           W( ELTVAR( IELPTR + I) )
     &           + ABS(A_ELT( K ))
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
               TEMP = TEMP + ABS( A_ELT(K))
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) = 
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) = 
     &        W(ELTVAR( IELPTR + J )) + ABS(A_ELT( K ))
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) = 
     &           W(ELTVAR( IELPTR + J )) + ABS(A_ELT( K ))
              W(ELTVAR( IELPTR + I ) ) = 
     &           W(ELTVAR( IELPTR + I )) + ABS(A_ELT( K ))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_119
      SUBROUTINE DMUMPS_135(MTYPE, N, 
     *    NELT, ELTPTR, LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *    W, KEEP,KEEP8, COLSCA )
      IMPLICIT NONE
      INTEGER MTYPE, N, NELT, LELTVAR, NA_ELT
      INTEGER ELTPTR(NELT+1), ELTVAR(LELTVAR)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION A_ELT(NA_ELT), W(N), COLSCA(N), 
     *                 TEMP, TEMP2
      INTEGER K, I, J, IEL, SIZEI, IELPTR
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO = 0.0D0)
      W = DZERO
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( KEEP(50).EQ.0 ) THEN
         IF (MTYPE.EQ.1) THEN
           DO J = 1, SIZEI
              TEMP2 = COLSCA(ELTVAR( IELPTR + J) )
              DO I = 1, SIZEI
               W( ELTVAR( IELPTR + I) ) =
     &           W( ELTVAR( IELPTR + I) )
     &           + ABS(A_ELT( K )* TEMP2 )
               K = K + 1
              END DO
            END DO
         ELSE
           DO J = 1, SIZEI
              TEMP = W( ELTVAR( IELPTR + J ) )
              TEMP2= COLSCA(ELTVAR( IELPTR + J) )
              DO I = 1, SIZEI
               TEMP = TEMP + ABS( A_ELT(K)*TEMP2)
               K = K + 1
              END DO
              W(ELTVAR( IELPTR + J )) =
     &          W(ELTVAR( IELPTR + J )) + TEMP
            END DO
         ENDIF
        ELSE
         DO J = 1, SIZEI
          W(ELTVAR( IELPTR + J )) =
     &        W(ELTVAR( IELPTR + J )) + 
     &        ABS( A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)) )
          K = K + 1
          DO I = J+1, SIZEI
              W(ELTVAR( IELPTR + J )) =
     &           W(ELTVAR( IELPTR + J )) + 
     &           ABS(A_ELT( K )*COLSCA(ELTVAR( IELPTR + J)))
              W(ELTVAR( IELPTR + I ) ) =
     &           W(ELTVAR( IELPTR + I )) + 
     &           ABS(A_ELT( K )*COLSCA(ELTVAR( IELPTR + I)))
              K = K + 1
          END DO
         ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_135
      SUBROUTINE DMUMPS_122( MTYPE, N, NELT, ELTPTR, 
     *                     LELTVAR, ELTVAR, NA_ELT, A_ELT,
     *                     SAVERHS, X, Y, W, K50 )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE, LELTVAR, NA_ELT
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( LELTVAR )
      DOUBLE PRECISION A_ELT( NA_ELT ), X( N ), Y( N ), 
     *                 SAVERHS(N), W(N)
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION TEMP, TEMP2
      PARAMETER( ZERO = 0.0D0 )
      Y = SAVERHS
      W = DBLE(ZERO)
      K = 1
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
          IF ( MTYPE .eq. 1 ) THEN
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     *          Y( ELTVAR( IELPTR + I ) ) -
     *             A_ELT( K ) * TEMP
                W( ELTVAR( IELPTR + I ) ) =
     *          W( ELTVAR( IELPTR + I ) ) +
     *             ABS( A_ELT( K ) * TEMP )
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              TEMP2 = W( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP - 
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                TEMP2 = TEMP2 +  ABS(
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
              W( ELTVAR( IELPTR + J ) ) = TEMP2
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     *      Y( ELTVAR( IELPTR + J ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            W( ELTVAR( IELPTR + J ) ) =
     *      W( ELTVAR( IELPTR + J ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     *        Y( ELTVAR( IELPTR + I ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     *        Y( ELTVAR( IELPTR + J ) ) -
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              W( ELTVAR( IELPTR + I ) ) =
     *        W( ELTVAR( IELPTR + I ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) ) )
              W( ELTVAR( IELPTR + J ) ) =
     *        W( ELTVAR( IELPTR + J ) ) + ABS(
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_122
