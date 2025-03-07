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
      SUBROUTINE ZMUMPS_301( id)
      USE ZMUMPS_STRUC_DEF
      USE ZMUMPS_BUFFER
      USE ZMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_headers.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE (ZMUMPS_STRUC), TARGET :: id
      INTEGER MP,LP, MPG
      LOGICAL PROK, PROKG
      INTEGER MTYPE, ICNTL20, ICNTL21
      LOGICAL LSCAL, ERANAL, GIVSOL
      INTEGER ICNTL10, ICNTL11
      INTEGER I,K,JPERM, J, II
      INTEGER IZ, NZ_THIS_BLOCK, IRHS_PTR_BEG, SHIFT_PTR
      INTEGER LIW,LIWW
      INTEGER LA, LWC, ISRW2, LIW_PASSED, LA_PASSED
      INTEGER ZMUMPS_LBUF, ZMUMPS_LBUF_INT
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF,
     *        IBEG_GLOB_DEF, IEND_GLOB_DEF,
     *        IROOT_DEF_RHS_COL1
      INTEGER NITREF, NOITER, SOLVET, KASE, JOBIREF
      COMPLEX*16 RSOL(1)
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16, DIMENSION(:), POINTER :: RHS_MUMPS
      COMPLEX*16, ALLOCATABLE :: SAVERHS(:), RW1(:),
     *                                 RW2(:), SRW1(:), SRW3(:),
     *                                 Y(:), W(:), D(:)
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV, 
     *                                      POSINRHSCOMP_N
      INTEGER LIWK_SOLVE
      INTEGER, ALLOCATABLE :: IW1(:), IWK_SOLVE(:)
      INTEGER, POINTER :: N, NZ, MAXIS, MAXS
      INTEGER, POINTER :: NRHS, LRHS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CNTL
      INTEGER, DIMENSION (:), POINTER :: KEEP,ICNTL,INFO
      INTEGER*8, DIMENSION (:), POINTER :: KEEP8
      INTEGER, DIMENSION (:), POINTER :: IS
      DOUBLE PRECISION, DIMENSION(:),POINTER::   RINFOG
      type scaling_data_t
        SEQUENCE
        COMPLEX*16, dimension(:), pointer :: SCALING
        COMPLEX*16, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      DOUBLE PRECISION ARRET
      COMPLEX*16 DUMMY(1)
      INTEGER IDUMMY(1), JDUMMY(1), KDUMMY(1), JJ, WHAT
      INTEGER allocok, PERLU
      INTEGER NBRHS, NBRHS_EFF, BEG_RHS, END_RHS, 
     &        IBEG, LD_RHS, KDEC, MINLWC, NBRHSMIN, 
     &        MASTER_ROOT, MASTER_ROOT_IN_COMM
      INTEGER IPT_RHS_ROOT, SIZE_ROOT, LD_REDRHS
      INTEGER IBEG_REDRHS, IBEG_RHSCOMP, LD_RHSCOMP, LENRHSCOMP
      INTEGER NB_K133
      LOGICAL WORKSPACE_MINIMAL_PREFERRED
      INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE_COPY
      COMPLEX*16, DIMENSION(:), POINTER :: RHS_SPARSE_COPY
#if defined(V_T)
      INTEGER soln_drive_class, glob_comm_ini, perm_scal_ini, soln_dist,
     &        soln_assem, perm_scal_post
#endif
      LOGICAL I_AM_SLAVE, MAXSUPDATED, BUILD_POSINRHSCOMP
      INTEGER  MTYPE_LOC, MAXSSAVED
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      N    =>id%N
      NZ   =>id%NZ
      NRHS =>id%NRHS
      LRHS =>id%LRHS
      MAXIS=>id%MAXIS
      MAXS =>id%MAXS
      CNTL =>id%CNTL
      KEEP =>id%KEEP
      KEEP8=>id%KEEP8
      IS   =>id%IS
      ICNTL=>id%ICNTL
      INFO =>id%INFO
      RINFOG =>id%RINFOG
      MAXSSAVED          = MAXS
      MAXSUPDATED        = .FALSE.
      BUILD_POSINRHSCOMP = .TRUE.
      NB_K133     = 2     
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln driver',soln_drive_class,ierr)
      CALL VTFUNCDEF( 'glob_comm_ini',soln_drive_class,
     &     glob_comm_ini,ierr)
      CALL VTFUNCDEF( 'perm_scal_ini',soln_drive_class,
     &     perm_scal_ini,ierr)
      CALL VTFUNCDEF( 'soln_dist',soln_drive_class,soln_dist,ierr)
      CALL VTFUNCDEF( 'soln_assem',soln_drive_class,soln_assem,ierr)
      CALL VTFUNCDEF( 'perm_scal_post',soln_drive_class,
     &     perm_scal_post,ierr)
#endif
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
      SIZE_ROOT   = -33333
      IF ( KEEP( 38 ) .ne. 0 ) THEN
            MASTER_ROOT = ZMUMPS_275(id%STEP( KEEP(38)),
     *                    id%PROCNODE_STEPS, id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%root%TOT_ROOT_SIZE
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN 
              SIZE_ROOT=id%SIZE_SCHUR
            ENDIF
      ELSE IF (KEEP( 20 ) .ne. 0 ) THEN
            MASTER_ROOT = ZMUMPS_275(id%STEP(KEEP(20)),
     *                    id%PROCNODE_STEPS, id%NSLAVES )
            IF (id%MYID_NODES .eq. MASTER_ROOT) THEN
              SIZE_ROOT = id%IS(
     *               id%PTLUST_S(id%STEP(KEEP(20)))+XSIZE)
            ELSE IF ((id%MYID.EQ.MASTER).AND.KEEP(60).NE.0) THEN 
              SIZE_ROOT=id%SIZE_SCHUR
            ENDIF
      ELSE
            MASTER_ROOT = -44444
      END IF
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LP  = ICNTL( 1 )
      PROK  = (MP.GT.0)
      PROKG = (MPG.GT.0 .and. id%MYID.eq.MASTER)
      IF ( PROK  ) WRITE(MP,100)
      IF ( PROKG ) WRITE(MPG,100)
      IF (id%MYID.EQ.MASTER) THEN
       CALL ZMUMPS_634(KEEP,ICNTL,MPG)
       IF (KEEP(111).eq.-1 .AND. NRHS .NE. KEEP(112)+KEEP(17))THEN
        INFO(1)=-32
        INFO(2)=NRHS
       ENDIF
       IF (KEEP(111).gt.0 .AND. NRHS .NE. 1) THEN
        INFO(1)=-32
        INFO(2)=1
       ENDIF
       IF (( KEEP(111) .LT. -1 ) .OR.
     *   (KEEP(111).GT.KEEP(122)+KEEP(17)) .OR.
     *   (KEEP(111) .EQ.-1 .AND. KEEP(122)+KEEP(17).EQ.0))
     *   THEN
        INFO(1)=-36
        INFO(2)=KEEP(111)
       ENDIF
      ENDIF
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) < 0) RETURN
      CALL MPI_BCAST(KEEP(111),1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      IF (id%MYID .eq. MASTER) THEN
        IF (KEEP(201) .EQ. 0 .OR. KEEP(84) .GT. 0) THEN
          NBRHS = ABS(KEEP(84))
        ELSE
          NBRHS = -2*KEEP(84)
        END IF
        IF (NBRHS .GT. NRHS ) NBRHS = NRHS
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(glob_comm_ini,ierr)
#endif
      CALL MPI_BCAST(NRHS,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(NBRHS,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
        IF (id%KEEP(201).NE.0) THEN
          WORKSPACE_MINIMAL_PREFERRED = .FALSE.
          IF (id%MYID .eq. MASTER) THEN
             KEEP(107) = MAX(0,KEEP(107))
             IF ((KEEP(107).EQ.0).AND.
     &            (KEEP(204).EQ.0).AND.(KEEP(211).NE.1) ) THEN
              WORKSPACE_MINIMAL_PREFERRED=.TRUE.
             ENDIF
          ENDIF
          CALL MPI_BCAST( KEEP(107), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(204), 1, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(208), 2, MPI_INTEGER,
     *                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( WORKSPACE_MINIMAL_PREFERRED, 1,
     *                  MPI_LOGICAL,
     *                  MASTER, id%COMM, IERR )
        ENDIF
      IF ( I_AM_SLAVE ) THEN
        MINLWC = NB_K133*KEEP(133)
        IF (id%KEEP(201).NE.0) THEN
           IF ( KEEP(209).EQ.-1 .AND. WORKSPACE_MINIMAL_PREFERRED)
     *     THEN 
              MAXSUPDATED = .TRUE.
              MAXS = MINLWC*NBRHS + KEEP(203)+1
           ELSE IF ( KEEP(209) .GT.0 ) THEN
              MAXSUPDATED = .TRUE.
              MAXS = MAX(KEEP(209), MINLWC*NBRHS + KEEP(203)+1)
           ENDIF
           ALLOCATE (id%S(MAXS), stat = allocok)
           IF ( allocok .GT. 0 ) THEN
             WRITE(*,*) ' Problem reallocation of MAXS '
             INFO(1) = -13
             INFO(2) = MAXS
           ENDIF
        ENDIF
        IF ( KEEP( 38 ) .NE. 0 .OR. KEEP( 20 ) .NE. 0 ) THEN 
          IF ( MASTER_ROOT .eq. id%MYID_NODES ) THEN
              NB_K133 = NB_K133 + 1
          END IF
        ENDIF
        IF(KEEP(201).EQ.0)THEN
           LA  = KEEP(31)
           LWC = MAXS - LA
           IF (LWC.LT.MINLWC) THEN
             INFO(1) = -11
             INFO(2) = MINLWC + KEEP(31) - MAXS
           ENDIF
        ELSE
           IF ( MAXS .LT. MINLWC + KEEP(203) ) THEN
             INFO(1) = -11
             INFO(2) = MINLWC + KEEP(203) - MAXS
           ENDIF
           LWC = MIN(MAXS - KEEP(203),
     *               MAX( INT(DBLE(MAXS*0.05)), MINLWC ))
           LA  = MAXS-LWC
           IF(LA.GT.(KEEP(31)+KEEP(203)*(KEEP(107)+1)))THEN
              LA=KEEP(31)+KEEP(203)*(KEEP(107)+1)
              LWC=MAXS-LA
           ENDIF
        ENDIF
        IF ( LWC .LT. NBRHS*MINLWC ) THEN
           NBRHS =  LWC / MINLWC
        ENDIF
      ENDIF
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) < 0) RETURN
      CALL MPI_ALLREDUCE(NBRHS, NBRHSMIN, 1, MPI_INTEGER, MPI_MIN,
     &                   id%COMM, IERR)
      NBRHS = NBRHSMIN
      IF ( NBRHS .LT. 1) THEN
         INFO(1) = -11
         INFO(2) = 3*KEEP(133)
         RETURN
      ENDIF
      IF ( I_AM_SLAVE ) THEN
         IF (KEEP(201).NE.0) THEN
          CALL ZMUMPS_590(LA)
          CALL ZMUMPS_586(id)
          IF (INFO(1).LT.0) THEN
             GOTO 111
          ENDIF
        ENDIF
      ENDIF
      IF (id%MYID .eq. MASTER) THEN
        MTYPE = ICNTL(  9 )
        IF ( PROKG )  THEN 
           WRITE( MPG, 150 )
     *             NRHS, NBRHS, ICNTL(9), ICNTL(10), ICNTL(11),
     *             ICNTL(20), ICNTL(21)
           IF (KEEP(111).NE.0) THEN 
            WRITE (MPG, 151) KEEP(111)
           ENDIF
           IF (KEEP(221).NE.0) THEN 
            WRITE (MPG, 152) KEEP(221)
           ENDIF
        ENDIF
        ICNTL20 = ICNTL(20)
        ICNTL21 = ICNTL(21)
        IF (ICNTL20.ne.0.and.ICNTL20.ne.1) ICNTL20=0
        IF (ICNTL20 .NE.0.AND.KEEP(111).NE.0) THEN
          IF (PROKG) WRITE(MPG,'(A)')
     *    ' WARNING: ICNTL(20) treated as if set to 0 (null space)'
          ICNTL20 = 0
        ENDIF
        IF (ICNTL21.ne.0.and.ICNTL21.ne.1) ICNTL21=0
        LSCAL = (((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 6)) .OR. (
     *    KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        ERANAL = ((ICNTL(11) .GT. 0) .OR. (ICNTL(10) .GT. 0))
        IF ( (KEEP(55).eq.0) .AND. KEEP(54).eq.0 .AND. 
     &      .NOT.ASSOCIATED(id%A) ) THEN
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        ELSE
          ICNTL10 = ICNTL(10)
          ICNTL11 = ICNTL(11)
        ENDIF
        IF (KEEP(111).NE.0) THEN
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *    ' WARNING: ICNTL(10) treated as if set to 0 (null space)'
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *    ' WARNING: ICNTL(11) treated as if set to 0 (null space)'
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        END IF
        IF (KEEP(221).NE.0) THEN
          IF (ICNTL10 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *    ' WARNING: ICNTL(10) treated as if set to 0 (reduced RHS))'
          ENDIF
          IF (ICNTL11 .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *    ' WARNING: ICNTL(11) treated as if set to 0 (reduced RHS)'
          ENDIF
          ICNTL10 = 0
          ICNTL11 = 0
          ERANAL = .FALSE.
        END IF
        IF ((ERANAL .AND. NBRHS > 1) .OR. ICNTL(21) > 0) THEN
          IF (ICNTL11 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *     ' WARNING: ICNTL(11) treated as if set to zero'
            ICNTL11=0
          ENDIF
          IF (ICNTL10 > 0) THEN
            IF (PROKG) WRITE(MPG,'(A)')
     *     ' WARNING: ICNTL(10) treated as if set to zero'
            ICNTL10=0
          ENDIF
          ERANAL = .FALSE.     
        ENDIF
        IF (ERANAL) THEN
          ALLOCATE(SAVERHS(N*NBRHS),stat = allocok)
            IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating SAVERHS'
              INFO(1) = -13
              INFO(2) = N*NBRHS
              GOTO 111
            END IF
        ENDIF
      END IF
      CALL MPI_BCAST(ICNTL10,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL11,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL20,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ICNTL21,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(ERANAL,1,MPI_LOGICAL,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(LSCAL,1,MPI_LOGICAL,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(MTYPE,1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
      CALL MPI_BCAST(KEEP(111),1,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
        id%LBUFR_BYTES = ( 4 + KEEP(133) ) * KEEP(34) +
     *                      KEEP(133) * NBRHS * KEEP(35)
        id%LBUFR = ( id%LBUFR_BYTES * KEEP(34) - 1 ) / KEEP(34)
        IF ( ASSOCIATED (id%BUFR) ) DEALLOCATE(id%BUFR)
        ALLOCATE (id%BUFR(id%LBUFR),stat=allocok)
        IF ( allocok .GT. 0 ) THEN
              WRITE(*,*) ' Problem in solve: error allocating BUFR'
              INFO(1) = -13
              INFO(2) = id%LBUFR
              GOTO 111
        ENDIF
      IF ( I_AM_SLAVE ) THEN
        ZMUMPS_LBUF_INT = ( 20 + id%NSLAVES * id%NSLAVES  * 4 )
     *                 * KEEP(34)
        ZMUMPS_LBUF = id%LBUFR_BYTES * id%NSLAVES + 3 * KEEP(34)
        CALL ZMUMPS_55( ZMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ZMUMPS_LBUF_INT
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating small Send buffer:IERR=',IERR
          END IF
          GOTO 111
        END IF
        CALL ZMUMPS_53( ZMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          INFO(1) = -13
          INFO(2) = ZMUMPS_LBUF
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,*) id%MYID,
     *      ':Error allocating Send buffer:IERR=', IERR
          END IF
          GOTO 111
        END IF
      ENDIF
      IF ( id%MYID .NE. MASTER .or.
     *     (id%MYID .EQ. MASTER .AND. ICNTL21 .NE.0 .AND.
     *       ( icntl20.ne.0 .OR. KEEP(111).NE.0 )
     *  )) THEN
        ALLOCATE(RHS_MUMPS(N*NBRHS),stat=ierr)
        IF ( ierr .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N*NBRHS
          IF (LP > 0)
     *      WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
          GOTO 111
        END IF
      ELSE
        RHS_MUMPS=>id%RHS
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        LD_RHSCOMP = MAX(KEEP(89),1)
        IF (KEEP(221).EQ.2) THEN
           IF (.NOT.ASSOCIATED(id%RHSCOMP)) THEN
             INFO(1) = -35
             INFO(2) = 1
             GOTO 111
           ENDIF
           IF (.NOT.ASSOCIATED(id%POSINRHSCOMP)) THEN
             INFO(1) = -35
             INFO(2) = 2
             GOTO 111
           ENDIF
        ELSE IF (KEEP(221).EQ.1) THEN
          IF (ASSOCIATED(id%RHSCOMP)) DEALLOCATE(id%RHSCOMP)
          LENRHSCOMP = LD_RHSCOMP*NRHS
          ALLOCATE (id%RHSCOMP(LENRHSCOMP))
          IF (ASSOCIATED(id%POSINRHSCOMP)) DEALLOCATE(id%POSINRHSCOMP)
          ALLOCATE (id%POSINRHSCOMP(KEEP(28)) )
        ELSE
          LENRHSCOMP = LD_RHSCOMP*NBRHS
          IF (ASSOCIATED(id%RHSCOMP)) DEALLOCATE(id%RHSCOMP)
          ALLOCATE (id%RHSCOMP(LENRHSCOMP))
          IF (ASSOCIATED(id%POSINRHSCOMP)) DEALLOCATE(id%POSINRHSCOMP)
          ALLOCATE (id%POSINRHSCOMP(KEEP(28)) )
        ENDIF
        LIWK_SOLVE = 4 * KEEP(28) + 1
        ALLOCATE ( IWK_SOLVE( LIWK_SOLVE), stat = allocok )
        IF (allocok .GT. 0 ) THEN
         INFO(1)=-13
         INFO(2)=LIWK_SOLVE
         GOTO 111
        END IF
        LIW = KEEP(32)
        LIWW  = MAXIS - LIW
        IF ( LIWW .LT. 0 ) THEN
          IF (LP>0) THEN
            WRITE(LP,*) 'ERROR in ZMUMPS_301: LIWW<0'
          ENDIF
          CALL ZMUMPS_ABORT()
        END IF
        ISRW2   = LA + 1
        ALLOCATE(SRW3(KEEP(133)), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=KEEP(133)
          GOTO 111
        END IF
        IF ( (KEEP(111).NE.0) .OR. (ICNTL20.NE.0) ) THEN
          ALLOCATE(POSINRHSCOMP_N(id%N), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
              IF (LP.GT.0) WRITE(LP,*)
     &       ' ERROR in ZMUMPS_301: allocating POSINRHSCOMP_N'
              INFO(1) = -13
              INFO(2) = N
              GOTO 111
          END IF
        END IF
      ELSE
        LIW=0
      END IF
 111  CONTINUE
#if defined(V_T)
      CALL VTEND(glob_comm_ini,ierr)
#endif
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF ( ICNTL21==1 ) THEN
        IF (LSCAL) THEN
          IF (id%MYID.NE.MASTER) THEN
            IF (MTYPE == 1) THEN
              ALLOCATE(id%COLSCA(id%N),stat=allocok)
            ELSE
              ALLOCATE(id%ROWSCA(id%N),stat=allocok)
            ENDIF
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
                WRITE(LP,*) 'Error allocating temporary scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
          IF (MTYPE == 1) THEN
              CALL MPI_BCAST(id%COLSCA,id%N,
     *                       MPI_DOUBLE_COMPLEX,MASTER,
     *                       id%COMM,IERR)
              scaling_data%SCALING=>id%COLSCA
          ELSE
              CALL MPI_BCAST(id%ROWSCA,id%N,
     *                       MPI_DOUBLE_COMPLEX,MASTER,
     *                       id%COMM,IERR)
              scaling_data%SCALING=>id%ROWSCA
          ENDIF
          IF (I_AM_SLAVE) THEN
            ALLOCATE(scaling_data%SCALING_LOC(id%KEEP(89)),
     *               stat=allocok)
            IF (allocok > 0) THEN
              IF (LP > 0) THEN
                WRITE(LP,*) 'Error allocating local scaling array'
              ENDIF
              INFO(1)=-13
              INFO(2)=id%KEEP(89)
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED=MAX(1,LIW)
          IF (KEEP(89) .GT. 0) THEN
            CALL ZMUMPS_535( MTYPE, id%ISOL_LOC(1),
     *               id%LSOL_LOC, id%PTLUST_S(1),
     *               id%KEEP(1),id%KEEP8(1),
     *               id%IS(1), LIW_PASSED,id%MYID_NODES,
     *               id%N, id%STEP(1), id%PROCNODE_STEPS(1),
     *               id%NSLAVES, scaling_data, LSCAL )
          ENDIF
          IF (id%MYID.NE.MASTER .AND. LSCAL) THEN
            IF (MTYPE == 1) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ELSE
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          IF (id%MYID.NE.MASTER) THEN
            ALLOCATE(id%UNS_PERM(id%N),stat=allocok)
            IF (allocok > 0) THEN
              INFO(1)=-13
              INFO(2)=id%N
              GOTO 40
            ENDIF
          ENDIF
        ENDIF
 40     CONTINUE
        CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
        IF (INFO(1) .LT.0 ) GOTO 90
        IF (KEEP(23) .NE. 0 .AND. MTYPE==1) THEN
          CALL MPI_BCAST(id%UNS_PERM,id%N,MPI_INTEGER,MASTER,
     *               id%COMM,IERR)
          IF (I_AM_SLAVE) THEN
            DO I=1, KEEP(89)
              id%ISOL_LOC(I) = id%UNS_PERM(id%ISOL_LOC(I))
            ENDDO
          ENDIF
          IF (id%MYID.NE.MASTER) THEN
            DEALLOCATE(id%UNS_PERM)
            NULLIFY(id%UNS_PERM)
          ENDIF
        ENDIF
      ENDIF
      IF ( ( KEEP(221) .EQ. 1 ) .OR.
     &     ( KEEP(221) .EQ. 2 ) 
     &   ) THEN
         IF (KEEP(46).EQ.1) THEN
             MASTER_ROOT_IN_COMM=MASTER_ROOT
         ELSE
             MASTER_ROOT_IN_COMM =MASTER_ROOT+1
         ENDIF
         IF ( id%MYID .EQ. MASTER ) THEN
             IF (NRHS.EQ.1) THEN
               LD_REDRHS = id%SIZE_SCHUR
             ELSE
               LD_REDRHS = id%LREDRHS
             ENDIF
         ENDIF
         IF (MASTER.NE.MASTER_ROOT_IN_COMM) THEN
            IF ( id%MYID .EQ. MASTER ) THEN
             CALL MPI_SEND(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
            ELSEIF ( id%MYID.EQ.MASTER_ROOT_IN_COMM) THEN
             CALL MPI_RECV(LD_REDRHS,1,MPI_INTEGER,
     &       MASTER, 0, id%COMM,STATUS,IERR)
            ENDIF
         ENDIF
      ENDIF
      DO BEG_RHS=1, NRHS, NBRHS
        NBRHS_EFF    = MIN(NRHS-BEG_RHS+1, NBRHS)
        END_RHS      = BEG_RHS + NBRHS_EFF - 1
        IF ( (KEEP(111).EQ.0.AND.ICNTL20.eq.0)
     *        .or. ICNTL21.eq.0 ) THEN
          IF (id%MYID .eq. MASTER) THEN
            IF (NRHS.GT.1) THEN 
              LD_RHS    = LRHS
            ELSE
              LD_RHS    = N
            ENDIF
            IBEG      = (BEG_RHS-1) * LD_RHS + 1
          ELSE 
            LD_RHS    = N
            IBEG      = 1
          END IF
        ELSE
          LD_RHS = N
          IBEG   = 1
        ENDIF
        IF (id%MYID.EQ.MASTER .AND. KEEP(221).NE.0) THEN
          IBEG_REDRHS= (BEG_RHS-1)*LD_REDRHS + 1
        ELSE
          IBEG_REDRHS=-142424  ! Should not be used
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          IF ( KEEP(221).EQ.0 ) THEN
             IBEG_RHSCOMP= 1
          ELSE
             IBEG_RHSCOMP= (BEG_RHS-1)*LD_RHSCOMP + 1
          ENDIF
        ELSE
          IBEG_RHSCOMP=-152525  ! Should not be used
        ENDIF
#if defined(V_T)
      CALL VTBEGIN(perm_scal_ini,ierr)
#endif
      IF (id%MYID .eq. MASTER) THEN
        IF (ICNTL20==1) THEN
          NZ_THIS_BLOCK=id%IRHS_PTR(BEG_RHS+NBRHS_EFF)-
     *                    id%IRHS_PTR(BEG_RHS)
          IF (KEEP(23) .NE. 0 .and. MTYPE .NE. 1) THEN
            ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
          ELSE
            IRHS_SPARSE_COPY
     *      =>
     *            id%IRHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     *                        id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
          ENDIF
          IF (LSCAL) THEN
            ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK), stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=NZ_THIS_BLOCK
              GOTO 30
            endif
          ELSE
            RHS_SPARSE_COPY
     *         => id%RHS_SPARSE(id%IRHS_PTR(BEG_RHS):
     *                       id%IRHS_PTR(BEG_RHS)+NZ_THIS_BLOCK-1)
          ENDIF
        ENDIF
        IF (KEEP(23) .NE. 0) THEN
          IF (MTYPE .NE. 1) THEN
            IF (ICNTL20==0) THEN
              ALLOCATE( RW1( N ),stat =allocok )
              IF ( allocok .GT. 0 ) THEN
                INFO(1)=-13
                INFO(2)=N
                IF ( LP .GT. 0 ) THEN
                  WRITE(LP,*) id%MYID,
     *            ':Error allocating RW1 in ZMUMPS_SOLVE_DRIVE'
                END IF
                GOTO 30
              END IF
              DO K = 1, NBRHS_EFF
               KDEC = IBEG+(K-1)*LD_RHS
               DO I = 1, N
                RW1(I)=RHS_MUMPS(I-1+KDEC)
               END DO
               DO I = 1, N
                JPERM = id%UNS_PERM(I)
                RHS_MUMPS(I-1+KDEC) = RW1(JPERM)
               END DO
              END DO
              DEALLOCATE(RW1)
            ELSE
              ALLOCATE(UNS_PERM_INV(N),stat=allocok) !FIXME: build it outside...
              if (allocok .GT.0 ) THEN
                INFO(1)=-13
                INFO(2)=N
                GOTO 30
              endif
              DO I = 1, N
                UNS_PERM_INV(id%UNS_PERM(I))=I
              ENDDO
              DO I = id%IRHS_PTR(BEG_RHS),
     *               id%IRHS_PTR(BEG_RHS+NBRHS_EFF)-1
                JPERM = UNS_PERM_INV(id%IRHS_SPARSE(I))
                IRHS_SPARSE_COPY(I-id%IRHS_PTR(BEG_RHS)+1)=JPERM
              ENDDO
              DEALLOCATE(UNS_PERM_INV) !FIXME: free it oustide
            ENDIF
          ENDIF
        ENDIF
        IF (ERANAL) THEN
         IF ( ICNTL20 == 0 ) THEN
          DO K = 1, NBRHS_EFF
            KDEC = IBEG+(K-1)*LD_RHS
            DO I = 1, N
              SAVERHS(I+(K-1)*N) = RHS_MUMPS(KDEC+I-1)
            END DO
          ENDDO
         ENDIF
        ENDIF
        IF (LSCAL) THEN
         IF (ICNTL20==0) THEN
          IF (MTYPE .EQ. 1) THEN
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%ROWSCA(I)
             END DO
            ENDDO
          ELSE
            DO K =1, NBRHS_EFF 
             KDEC = (K-1) * LD_RHS + IBEG - 1
             DO I = 1, N
              RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%COLSCA(I)
             END DO
            ENDDO
          ENDIF
         ELSE
          KDEC=id%IRHS_PTR(BEG_RHS)
          IF (MTYPE .eq. 1) THEN
            DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*id%ROWSCA(I)
            ENDDO
          ELSE
            DO IZ=1,NZ_THIS_BLOCK
              I=IRHS_SPARSE_COPY(IZ)
              RHS_SPARSE_COPY(IZ)=id%RHS_SPARSE(KDEC+IZ-1)*id%COLSCA(I)
            ENDDO
          ENDIF
         ENDIF
        END IF
      ENDIF
#if defined(V_T)
      CALL VTEND(perm_scal_ini,ierr)
#endif
 30   CONTINUE
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1) .LT.0 ) GOTO 90
      IF ( I_AM_SLAVE ) THEN
       IF ( (KEEP(111).NE.0) .OR. (ICNTL20.NE.0) ) THEN
         IF (BUILD_POSINRHSCOMP) THEN
           IF (KEEP(111).NE.0) THEN
             WHAT      = 2
             MTYPE_LOC = 1
           ELSE
             WHAT      = 1
             MTYPE_LOC = MTYPE
           ENDIF
           LIW_PASSED=MAX(1,LIW)
           CALL ZMUMPS_639(id%NSLAVES,id%N,
     *           id%MYID_NODES, id%PTLUST_S(1),
     *           id%KEEP(1),id%KEEP8(1), 
     *           id%PROCNODE_STEPS(1), id%IS(1), LIW_PASSED, 
     *           id%STEP(1), 
     *           id%POSINRHSCOMP(1), POSINRHSCOMP_N(1), 
     *           id%N, MTYPE_LOC,
     *           WHAT )
           BUILD_POSINRHSCOMP = .FALSE.
         ENDIF
       ENDIF
      ENDIF
#if defined(V_T)
      CALL VTBEGIN(soln_dist,ierr)
#endif
      IF (KEEP(111).eq.0) THEN
        IF (ICNTL20 == 0) THEN
          IF ( .NOT.I_AM_SLAVE ) THEN
            CALL ZMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     *          IDUMMY, 1,
     *          id%STEP(1), KDUMMY, 1, BUILD_POSINRHSCOMP,
     *          id%ICNTL(1),id%INFO(1))
            BUILD_POSINRHSCOMP=.FALSE.
          ELSE
            LIW_PASSED = MAX( LIW, 1 )
            CALL ZMUMPS_638(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     *          id%PROCNODE_STEPS(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), id%POSINRHSCOMP(1), KEEP(28), 
     *          BUILD_POSINRHSCOMP,
     *          id%ICNTL(1),id%INFO(1))
            BUILD_POSINRHSCOMP=.FALSE.
          ENDIF
          IF (INFO(1).LT.0) GOTO 90
        ELSE
         CALL MPI_BCAST( NZ_THIS_BLOCK,1, MPI_INTEGER,
     *                   MASTER, id%COMM,IERR)
         IF (id%MYID==MASTER) THEN
           IRHS_PTR_BEG=BEG_RHS
         ELSE
           ALLOCATE(IRHS_SPARSE_COPY(NZ_THIS_BLOCK))
           ALLOCATE(RHS_SPARSE_COPY(NZ_THIS_BLOCK))
           ALLOCATE(id%IRHS_PTR(NBRHS_EFF+1))
           IRHS_PTR_BEG=1
         ENDIF
         CALL MPI_BCAST(IRHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(RHS_SPARSE_COPY,
     *                NZ_THIS_BLOCK,
     *                MPI_DOUBLE_COMPLEX,
     *                MASTER, id%COMM,IERR)
         CALL MPI_BCAST(id%IRHS_PTR(IRHS_PTR_BEG),
     *                NBRHS_EFF+1,
     *                MPI_INTEGER,
     *                MASTER, id%COMM,IERR)
         SHIFT_PTR=id%IRHS_PTR(IRHS_PTR_BEG)-1
         IF ( I_AM_SLAVE ) THEN
           DO K = 1, NBRHS_EFF
            KDEC = (K-1) * LD_RHS + IBEG - 1
            RHS_MUMPS(KDEC+1:KDEC+id%N)=DCMPLX(ZERO)
            DO IZ=id%IRHS_PTR(IRHS_PTR_BEG+K-1)-SHIFT_PTR,
     *          id%IRHS_PTR(IRHS_PTR_BEG+K)-1-SHIFT_PTR
              I=IRHS_SPARSE_COPY(IZ)
              IF (POSINRHSCOMP_N(I).NE.0) THEN
               RHS_MUMPS(KDEC+I)= RHS_SPARSE_COPY(IZ)
              ENDIF
            ENDDO
           ENDDO
         END IF
         IF (id%MYID .ne. MASTER) THEN
           DEALLOCATE(IRHS_SPARSE_COPY)
           DEALLOCATE(RHS_SPARSE_COPY)
           DEALLOCATE(id%IRHS_PTR)
           NULLIFY(IRHS_SPARSE_COPY)
           NULLIFY(RHS_SPARSE_COPY)
           NULLIFY(id%IRHS_PTR)
         ELSE
           IF (KEEP(23).ne.0 .and. MTYPE.ne.1) THEN
             DEALLOCATE(IRHS_SPARSE_COPY)
           ENDIF
           NULLIFY(IRHS_SPARSE_COPY)
           IF (LSCAL) THEN
             DEALLOCATE(RHS_SPARSE_COPY)
           ENDIF
           NULLIFY(RHS_SPARSE_COPY)
         ENDIF
        ENDIF
#if defined(try_null_space)
      ELSE
        IF (KEEP(111).GT.0) THEN
          IBEG_GLOB_DEF = KEEP(111)
          IEND_GLOB_DEF = KEEP(111)
        ELSE
          IBEG_GLOB_DEF = BEG_RHS
          IEND_GLOB_DEF = BEG_RHS+NBRHS_EFF-1
        ENDIF
        DO K=1, NBRHS_EFF
          KDEC = (K-1) *LD_RHSCOMP
          id%RHSCOMP(KDEC+1:KDEC+LD_RHSCOMP)=DCMPLX(ZERO)
        END DO
        DO I=MAX(IBEG_GLOB_DEF,KEEP(220)),
     *       MIN(IEND_GLOB_DEF,KEEP(220)+KEEP(109)-1)
          JJ= POSINRHSCOMP_N(id%PIVNUL_LIST(I-KEEP(220)+1))
          IF (JJ.GT.0) 
     *     id%RHSCOMP(IBEG_RHSCOMP -1+ (I-IBEG_GLOB_DEF)*LD_RHSCOMP 
     *                + JJ) =  DCMPLX(ONE)
        ENDDO
        IF ( KEEP(17).NE.0 .and. id%MYID_NODES.EQ.MASTER_ROOT) THEN
            IBEG_ROOT_DEF  = MAX(IBEG_GLOB_DEF,KEEP(112)+1)
            IEND_ROOT_DEF  = MIN(IEND_GLOB_DEF,KEEP(112)+KEEP(17))
            IROOT_DEF_RHS_COL1 = IBEG_ROOT_DEF-IBEG_GLOB_DEF + 1
            IBEG_ROOT_DEF = IBEG_ROOT_DEF-KEEP(112)
            IEND_ROOT_DEF = IEND_ROOT_DEF-KEEP(112)
        ELSE
            IBEG_ROOT_DEF = -90999
            IEND_ROOT_DEF = -90999
        ENDIF
#endif
      ENDIF
      IF ( id%MYID_NODES .EQ. MASTER_ROOT ) THEN
        IPT_RHS_ROOT = LWC - NBRHS_EFF * SIZE_ROOT + 1
      ELSE
        IPT_RHS_ROOT = LWC ! Will be passed, but not accessed
      ENDIF
      IF (KEEP(221) .EQ. 2 ) THEN
         IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &        ( id%MYID .EQ. MASTER ) ) THEN
            II = ISRW2+IPT_RHS_ROOT-2
            DO K=1, NBRHS_EFF
             KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
             DO I = 1, SIZE_ROOT
              id%S(II+I) = id%REDRHS(KDEC+I)
             ENDDO
             II = II+SIZE_ROOT
            ENDDO
         ELSE
          IF ( id%MYID .EQ. MASTER) THEN
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               KDEC = IBEG_REDRHS
               CALL MPI_SEND(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
            ELSE
              DO K=1, NBRHS_EFF
                KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
                CALL MPI_SEND(id%REDRHS(KDEC),SIZE_ROOT,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,IERR)
              ENDDO
            ENDIF
          ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
            II = ISRW2+IPT_RHS_ROOT-1
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               CALL MPI_RECV(id%S(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER, 0, id%COMM,STATUS,IERR)
            ELSE
             DO K=1, NBRHS_EFF
              CALL MPI_RECV(id%S(II),SIZE_ROOT,
     &           MPI_DOUBLE_COMPLEX,
     &           MASTER, 0, id%COMM,STATUS,IERR)
              II = II + SIZE_ROOT
             ENDDO
            ENDIF
          ENDIF
         ENDIF
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        LIW_PASSED = MAX( LIW, 1 )
        LA_PASSED  = MAX( LA, 1 )
        CALL ZMUMPS_245(id%root, N, id%S(1), LA_PASSED,
     *    IS(1), LIW_PASSED,
     *    id%S(ISRW2), LWC, IS( LIW + 1 ),
     *    LIWW, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF, 
     *    id%NA,id%LNA,id%NE_STEPS(1), SRW3, MTYPE,
     *    ICNTL, id%STEP(1), id%FRERE_STEPS(1), 
     *    id%DAD_STEPS(1), id%FILS(1),
     *    id%PTLUST_S(1), id%PTRFAC(1),
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS,
     *    id%NSLAVES, INFO, KEEP,KEEP8,
     *    id%COMM, id%COMM_NODES, id%MYID,
     *    id%MYID_NODES,
     *    id%BUFR,
     *    id%LBUFR, id%LBUFR_BYTES, 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     *    IBEG_ROOT_DEF, IEND_ROOT_DEF,
     *    IROOT_DEF_RHS_COL1,
     *    IPT_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,  
     *    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 
     *    id%POSINRHSCOMP(1), BUILD_POSINRHSCOMP )
      END IF
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF (INFO(1).eq.-2) then
        INFO(1)=-11
        write(*,*) ' WARNING : -11 error code obtained in solve'
      END IF
      IF (INFO(1).eq.-3) then
        INFO(1)=-14
        write(*,*) ' WARNING : -14 error code obtained in solve'
      END IF
      IF (INFO(1).LT.0) GO TO 90
      IF ( KEEP(221) .EQ. 1 ) THEN
         IF ( ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) .AND.
     &        ( id%MYID .EQ. MASTER ) ) THEN
            II = ISRW2+IPT_RHS_ROOT-2
            DO K=1, NBRHS_EFF
             KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS -1
             DO I = 1, SIZE_ROOT
              id%REDRHS(KDEC+I) = id%S(II+I) 
             ENDDO
             II = II+SIZE_ROOT
            ENDDO
         ELSE
          IF ( id%MYID .EQ. MASTER ) THEN
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               KDEC = IBEG_REDRHS
               CALL MPI_RECV(id%REDRHS(KDEC),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
            ELSE
             DO K=1, NBRHS_EFF
               KDEC = IBEG_REDRHS+(K-1)*LD_REDRHS
               CALL MPI_RECV(id%REDRHS(KDEC),SIZE_ROOT,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER_ROOT_IN_COMM, 0, id%COMM,
     &              STATUS,IERR)
             ENDDO
            ENDIF
          ELSE IF ( id%MYID .EQ. MASTER_ROOT_IN_COMM ) THEN
            II = ISRW2+IPT_RHS_ROOT-1
            IF (LD_REDRHS.EQ.SIZE_ROOT) THEN
               CALL MPI_SEND(id%S(II),
     &              SIZE_ROOT*NBRHS_EFF,
     &              MPI_DOUBLE_COMPLEX,
     &              MASTER, 0, id%COMM,IERR)
            ELSE
             DO K=1, NBRHS_EFF
              CALL MPI_SEND(id%S(II),SIZE_ROOT,
     &           MPI_DOUBLE_COMPLEX,
     &           MASTER, 0, id%COMM,IERR)
              II = II + SIZE_ROOT
             ENDDO
            ENDIF
          ENDIF
         ENDIF
      ENDIF
      IF ( KEEP(221) .NE. 1 ) THEN
       IF (ICNTL21 == 0) THEN
        LIW_PASSED = MAX( LIW, 1 )
        IF ( .NOT.I_AM_SLAVE ) THEN
        CALL ZMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          JDUMMY, id%KEEP(1), id%KEEP8(1), id%PROCNODE_STEPS(1),
     *          IDUMMY, 1,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ELSE
        CALL ZMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          id%PTLUST_S(1), id%KEEP(1), id%KEEP8(1),
     *          id%PROCNODE_STEPS(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ENDIF
        IF ( id%MYID.eq.MASTER .AND. LSCAL ) THEN
          IF (MTYPE .EQ. 1) THEN
             DO K= 1, NBRHS_EFF
              KDEC = (K-1) * LD_RHS + IBEG - 1
              DO I = 1, N
                RHS_MUMPS(KDEC+ I) = RHS_MUMPS(KDEC+ I) * id%COLSCA(I)
              END DO
             END DO
          ELSE
             DO K= 1, NBRHS_EFF
              KDEC = (K-1) * LD_RHS + IBEG - 1
              DO I = 1, N
                RHS_MUMPS(KDEC+I) = RHS_MUMPS(KDEC+I) * id%ROWSCA(I)
              END DO
             END DO
          ENDIF
        END IF
       ELSE ! Case of distributed solution
        IF ( I_AM_SLAVE ) THEN
         LIW_PASSED = MAX( LIW, 1 )
         IF ( KEEP(89) .GT. 0 ) THEN
           CALL ZMUMPS_532(id%NSLAVES,
     *          id%N, id%MYID,
     *          MTYPE, RHS_MUMPS(IBEG), LD_RHS, NBRHS_EFF,
     *          id%ISOL_LOC(1),
     *          id%SOL_LOC(1), BEG_RHS, id%LSOL_LOC,
     *          id%PTLUST_S(1), id%PROCNODE_STEPS(1),
     *          id%KEEP(1),id%KEEP8(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), scaling_data, LSCAL )
         ENDIF
        ENDIF
       ENDIF
      ENDIF
      IF ( ICNTL10 > 0 .AND. NBRHS_EFF > 1 ) THEN
        DO I = 1, ICNTL10
          write(*,*) 'FIXME: to be implemented'
        END DO
      END IF
      IF (ERANAL) THEN
        IF ((ICNTL10 .GT. 0) .AND. (ICNTL11 .GT. 0)) THEN
          IF (id%MYID .EQ. MASTER) THEN
            GIVSOL = .FALSE.
            IF (MP .GT. 0) WRITE( MP, 170 )
            ALLOCATE(RW1(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
            ALLOCATE(RW2(N),stat=allocok)
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
              GOTO 776
            endif
          END IF
          IF ( KEEP(54) .ne. 0 ) THEN
            ALLOCATE( SRW1( N ), stat =allocok )
            if (allocok .GT.0 ) THEN
              INFO(1)=-13
              INFO(2)=N
            endif
          END IF
 776      CONTINUE
          CALL ZMUMPS_276( ICNTL, INFO,
     *                  id%COMM,id%MYID)
          IF ( INFO(1) .LT. 0 ) GOTO 90
          IF ( KEEP(54) .eq. 0 ) THEN
            IF (id%MYID .EQ. MASTER) THEN
              IF (KEEP(55).EQ.0) THEN
                CALL ZMUMPS_278( ICNTL(9), N, NZ, id%A(1),
     *             id%IRN(1), id%JCN(1),
     *             RHS_MUMPS(IBEG), SAVERHS, RW1, RW2, KEEP,KEEP8 )
              ELSE
                CALL ZMUMPS_121( ICNTL(9), N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          RHS_MUMPS(IBEG), SAVERHS, RW1, RW2, KEEP,KEEP8 )
              ENDIF
            END IF
          ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_COMPLEX, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW2,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              RW2 = SAVERHS - RW2
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, RW1,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            DEALLOCATE( SRW1 )
          END IF
          IF ( id%MYID .EQ. MASTER )  THEN
            CALL ZMUMPS_205(ICNTL(9),INFO(1),N,NZ,
     *        RHS_MUMPS(IBEG), SAVERHS,RW1,RW2,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MP,ICNTL,
     *        KEEP,KEEP8)
            DEALLOCATE(RW1)
            DEALLOCATE(RW2)
          END IF
        END IF
      IF ( PROK  .AND. ICNTL10 .GT. 0 ) WRITE( MP, 270 )
      IF ( PROKG .AND. ICNTL10 .GT. 0 ) WRITE( MPG, 270 )
      ALLOCATE(Y(N), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      IF ( id%MYID .EQ. MASTER ) THEN
        ALLOCATE( IW1( 2 * N ),stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=2 * N
          GOTO 777
        ENDIF
        ALLOCATE( D(N),stat =allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        ALLOCATE( W(3*N), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1)=-13
          INFO(2)=N
          GOTO 777
        ENDIF
        NITREF = ICNTL10
        JOBIREF= ICNTL11
        IF ( PROKG .AND. ICNTL10 .GT. 0 )
     *    WRITE( MPG, 240) 'MAXIMUM NUMBER OF STEPS =', NITREF
        DO I = 1, N
          D( I ) = DCMPLX(ONE)
        END DO
      END IF
      ALLOCATE(SRW1(N),stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=N
        GOTO 777
      ENDIF
      KASE = 0
 777  CONTINUE
      CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
      IF ( INFO(1) .LT. 0 ) GOTO 90
  22    CONTINUE
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 0 ) THEN
              IF (KEEP(55).NE.0) THEN 
               CALL ZMUMPS_119(MTYPE, N, 
     *           id%NELT, id%ELTPTR(1), 
     *           id%LELTVAR, id%ELTVAR(1),
     *           id%NA_ELT, id%A_ELT,
     *           W(N+1), KEEP,KEEP8 )
              ELSE
               IF ( MTYPE .eq. 1 ) THEN
                 CALL ZMUMPS_207
     *       ( id%A(1), NZ, N, id%IRN(1), id%JCN(1), W(N+1), KEEP,KEEP8)
               ELSE
                 CALL ZMUMPS_207
     *       ( id%A(1), NZ, N, id%JCN(1), id%IRN(1), W(N+1), KEEP,KEEP8)
               END IF
              ENDIF
            ENDIF
          END IF
        ELSE
          IF ( KASE .eq. 0 ) THEN
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              IF ( MTYPE .eq. 1 ) THEN
              CALL ZMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
              ELSE
              CALL ZMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%JCN_loc, id%IRN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
              END IF
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W( N + 1 ),
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
          END IF
        END IF
        IF ( id%MYID .eq. MASTER ) THEN
            ARRET = CNTL(2)
            IF (ARRET .LT. 0.0D0) THEN
              ARRET = SQRT(EPSILON(0.0D0))
            END IF
            CALL ZMUMPS_206(NZ,N,SAVERHS,RHS_MUMPS(IBEG),
     *      Y, D, W,
     *      IW1, KASE,RINFOG(7),
     *      RINFOG(9), JOBIREF, RINFOG(10), NITREF, NOITER, MP,
     *      KEEP,KEEP8, ARRET )
        END IF
        IF ( KEEP(54) .ne. 0 ) THEN
          CALL MPI_BCAST( KASE, 1, MPI_INTEGER, MASTER,
     *    id%COMM, IERR )
        END IF
        IF ( KEEP(54) .eq. 0 ) THEN
          IF ( id%MYID .eq. MASTER ) THEN
            IF ( KASE .eq. 14 ) THEN
              IF (KEEP(55).NE.0) THEN
               CALL ZMUMPS_122( MTYPE, N, 
     *            id%NELT, id%ELTPTR, id%LELTVAR,
     *            id%ELTVAR, id%NA_ELT, id%A_ELT,
     *            SAVERHS, RHS_MUMPS(IBEG),
     *            Y, W, KEEP(50))
              ELSE
                 IF ( MTYPE .eq. 1 ) THEN
                   CALL ZMUMPS_208
     *    (id%A(1), NZ, N, id%IRN(1), id%JCN(1), SAVERHS,
     *    RHS_MUMPS(IBEG), Y, W, KEEP,KEEP8)
                 ELSE
                   CALL ZMUMPS_208
     *    (id%A(1), NZ, N, id%JCN(1), id%IRN(1), SAVERHS,
     *    RHS_MUMPS(IBEG), Y, W, KEEP,KEEP8)
                 END IF
              ENDIF
              GOTO 22
            END IF
          END IF
        ELSE
          IF ( KASE.eq.14 ) THEN
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_COMPLEX, MASTER,
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
              Y = SAVERHS - Y
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and. id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_193( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE )
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, W,
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY, 
     *          id%N, MPI_DOUBLE_COMPLEX,
     *          MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            GOTO 22
          END IF
        END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF ( KASE .GT. 0 ) THEN
          IF ( MTYPE .EQ. 1 ) THEN
            SOLVET = KASE - 1
          ELSE
            SOLVET = KASE
          END IF
          IF ( LSCAL ) THEN
            IF ( SOLVET .EQ. 1 ) THEN
              DO K = 1, N
                Y( K ) = Y( K ) * id%ROWSCA( K )
              END DO
            ELSE
              DO K = 1, N
                Y( K ) = Y( K ) * id%COLSCA( K )
              END DO
            END IF
          END IF
        END IF
      END IF
      CALL MPI_BCAST( KASE , 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      CALL MPI_BCAST( SOLVET, 1, MPI_INTEGER, MASTER,
     *                id%COMM, IERR)
      IF ( KASE .GT. 0 ) THEN
        CALL MPI_BCAST( Y, N, MPI_DOUBLE_COMPLEX, MASTER,
     *                id%COMM, IERR )
        IF ( id%MYID_NODES .EQ. MASTER_ROOT ) THEN
          IPT_RHS_ROOT = LWC - NBRHS_EFF * SIZE_ROOT + 1
        ELSE
          IPT_RHS_ROOT = LWC
        ENDIF
        IF ( I_AM_SLAVE ) THEN
          LIW_PASSED = MAX( LIW, 1 )
          LA_PASSED = MAX( LA, 1 )
          CALL ZMUMPS_245( id%root, N,
     *    id%S(1), LA_PASSED, id%IS( 1 ),
     *    LIW_PASSED, id%S( ISRW2 ), LWC, id%IS( LIW + 1 ),
     *    LIWW,
     *    Y, N, NBRHS_EFF, id%NA, id%LNA, id%NE_STEPS,
     *    SRW1, SOLVET, ICNTL,
     *    id%STEP(1), id%FRERE_STEPS(1), id%DAD_STEPS(1), 
     *    id%FILS(1),
     *    id%PTLUST_S(1), id%PTRFAC(1),
     *    IWK_SOLVE, LIWK_SOLVE,
     *    id%PROCNODE_STEPS, id%NSLAVES, INFO, KEEP,KEEP8,
     *    id%COMM,
     *    id%COMM_NODES,
     *    id%MYID, id%MYID_NODES,
     *    id%BUFR, id%LBUFR, id%LBUFR_BYTES , 
     *
     *    id%ISTEP_TO_INIV2(1), id%TAB_POS_IN_PERE(1,1),
     *    IBEG_ROOT_DEF, IEND_ROOT_DEF,
     *    IROOT_DEF_RHS_COL1,
     *    IPT_RHS_ROOT, SIZE_ROOT, MASTER_ROOT, 
     *    id%RHSCOMP(IBEG_RHSCOMP), LD_RHSCOMP, 
     *    id%POSINRHSCOMP(1), BUILD_POSINRHSCOMP )
        END IF
        CALL ZMUMPS_276( ICNTL, INFO,
     *                   id%COMM,id%MYID)
        IF (INFO(1).eq.-2) INFO(1)=-12
        IF (INFO(1).eq.-3) INFO(1)=-15
        IF (INFO(1).LT.0) GO TO 90
        LIW_PASSED = MAX( LIW, 1 )
        IF ( .NOT. I_AM_SLAVE ) THEN
        CALL ZMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, Y, LD_RHS, NBRHS_EFF,
     *          JDUMMY, id%KEEP(1),id%KEEP8(1), id%PROCNODE_STEPS(1),
     *          IDUMMY, 1,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ELSE
        CALL ZMUMPS_521(id%NSLAVES,id%N, id%MYID, id%COMM,
     *          MTYPE, Y, LD_RHS, NBRHS_EFF,
     *          id%PTLUST_S(1), id%KEEP(1),id%KEEP8(1),
     *          id%PROCNODE_STEPS(1),
     *          IS(1), LIW_PASSED,
     *          id%STEP(1), id%BUFR(1), id%LBUFR, id%LBUFR_BYTES)
        ENDIF
        IF ( id%MYID.eq.MASTER) THEN
          IF (LSCAL) THEN
            IF (SOLVET .EQ. 1) THEN
               DO K = 1, N
                 Y(K) = Y(K) * id%COLSCA(K)
               END DO
            ELSE
               DO K = 1, N
                 Y(K) = Y(K) * id%ROWSCA(K)
               END DO
            ENDIF
          END IF
        END IF
        GO TO 22
      ELSEIF ( KASE .LT. 0 ) THEN
           INFO( 1 ) = INFO( 1 ) + 8
      END IF
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(W)
      IF ( id%MYID .eq. MASTER ) DEALLOCATE(IW1)
      IF ( PROKG .AND. NITREF .GT. 0 .AND.
     *id%MYID .EQ. MASTER ) THEN
        WRITE( MPG, 81 ) 
        WRITE( MPG, 141 ) 'NUMBER OF STEPS OF ITERATIVE REFINEMENTS  
     *=', NOITER
      ENDIF
      IF ( id%MYID .EQ. MASTER .AND. NITREF .GT. 0 ) THEN
        id%INFOG(15) = NOITER
      END IF
      IF ( PROK .AND. NITREF .GT.0 ) WRITE( MP, 131 ) 
      IF (ICNTL11 .GT. 0) THEN
        IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
              CALL ZMUMPS_278( MTYPE, N, NZ, id%A(1),
     *          id%IRN(1), id%JCN(1),
     *          RHS_MUMPS(IBEG), SAVERHS, Y, D, KEEP,KEEP8 )
            ELSE
              CALL ZMUMPS_121( MTYPE, N, 
     *          id%NELT, id%ELTPTR, 
     *          id%LELTVAR, id%ELTVAR,
     *          id%NA_ELT, id%A_ELT,
     *          RHS_MUMPS(IBEG), SAVERHS, Y, D, KEEP,KEEP8 )
            ENDIF
          END IF
        ELSE
            CALL MPI_BCAST( RHS_MUMPS(IBEG), id%N,
     *              MPI_DOUBLE_COMPLEX, MASTER, 
     *              id%COMM, IERR )
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_192( id%N, id%NZ_loc,
     *        id%IRN_loc, id%JCN_loc, id%A_loc,
     *        RHS_MUMPS(IBEG), SRW1, KEEP(50), MTYPE ) 
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, D,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
              D = SAVERHS - D
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
            IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
              CALL ZMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc,
     *          SRW1, id%KEEP,id%KEEP8 )
            ELSE
              SRW1 = DCMPLX(ZERO)
            END IF
            IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SRW1, Y,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            ELSE
              CALL MPI_REDUCE( SRW1, DUMMY,
     *        id%N, MPI_DOUBLE_COMPLEX,
     *        MPI_SUM,MASTER,id%COMM, IERR)
            END IF
        END IF
        IF (id%MYID .EQ. MASTER) THEN
         IF ((MPG .GT. 0) .AND. (NITREF .GT. 0)) WRITE( MPG, 65 ) 
         IF ((MPG .GT. 0) .AND. (NITREF .LE. 0)) WRITE( MPG, 170 ) 
         GIVSOL = .FALSE.
         CALL ZMUMPS_205(MTYPE,INFO(1),N,NZ,RHS_MUMPS(IBEG),
     *        SAVERHS,Y,D,GIVSOL,
     *        RSOL,RINFOG(4),RINFOG(5),RINFOG(6),MPG,ICNTL,
     *        KEEP,KEEP8)
         IF ( MPG .GT. 0 ) THEN
          WRITE( MPG, 115 )
     *'RINFOG(7):COMPONENTWISE SCALED RESIDUAL(W1)=', RINFOG(7)
          WRITE( MPG, 115 )
     *'------(8):---------------------------- (W2)=', RINFOG(8)
          WRITE( MPG, 115 )
     *'------(9):Upper bound ERROR ...............=', RINFOG(9)
          WRITE( MPG, 115 )
     *'-----(10):CONDITION NUMBER (1) ............=', RINFOG(10)
          WRITE( MPG, 115 )
     *'-----(11):CONDITION NUMBER (2) ............=', RINFOG(11)
         END IF
        END IF ! MASTER
      END IF ! ICNTL11 > 0
      IF (id%myid == MASTER) DEALLOCATE(D)
      DEALLOCATE(Y)
      DEALLOCATE(SRW1)
      END IF
      IF ( id%MYID .EQ. MASTER .AND. ICNTL21==0
     *     .AND. KEEP(23) .NE. 0) THEN
        IF ((KEEP(221).NE.1 .AND. ICNTL(9) .EQ. 1)
     *     .OR. KEEP(111) .NE.0) THEN
          ALLOCATE( RW1( N ),stat =allocok )
          IF ( allocok .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=N
            WRITE(*,*) 'could not allocate ', N, 'integers.'
            CALL ZMUMPS_ABORT()
          END IF
          DO K = 1, NBRHS_EFF
           KDEC = (K-1)*LD_RHS+IBEG-1
           DO 70 I = 1, N
            RW1(I) = RHS_MUMPS(KDEC+I)
 70        CONTINUE
           DO 80 I = 1, N
            JPERM = id%UNS_PERM(I)
            RHS_MUMPS( KDEC+JPERM ) = RW1( I )
 80        CONTINUE
          END DO
          DEALLOCATE( RW1 )
        END IF
      END IF
      IF (PROKG) WRITE( MPG, 120 ) id%INFOG(1), id%INFOG(2)
      IF (id%MYID.EQ.MASTER .and.ICNTL21==0.and.KEEP(221).NE.1) THEN
        IF ( INFO(1) .GE. 0 .AND. ICNTL(4).GE.3 .AND. ICNTL(3).GT.0)
     *    THEN
          K = min0(10, N)
          IF (ICNTL(4) .eq. 4 ) K = N
          J = min0(10,NBRHS_EFF)
          IF (ICNTL(4) .eq. 4 ) J = NBRHS_EFF
          DO II=1, J
            WRITE(ICNTL(3),110) BEG_RHS+II-1
            WRITE(ICNTL(3),160)
     *    (RHS_MUMPS(IBEG+(II-1)*LRHS+I-1),I=1,K)
          ENDDO
        END IF
      END IF
      ENDDO
 90   CONTINUE
      IF (KEEP(201).NE.0)THEN
        IF (I_AM_SLAVE) THEN
          CALL ZMUMPS_582(IERR)
          IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
        ENDIF
        CALL ZMUMPS_276( ICNTL, INFO,
     *         id%COMM,id%MYID)
      ENDIF
      IF (associated(id%BUFR)) THEN
          DEALLOCATE(id%BUFR)
          NULLIFY(id%BUFR)
      ENDIF
      IF ( I_AM_SLAVE ) THEN
        IF (allocated(IWK_SOLVE)) DEALLOCATE( IWK_SOLVE )
        CALL ZMUMPS_57( IERR )
        CALL ZMUMPS_59( IERR )
      END IF
      IF ( id%MYID .eq. MASTER ) THEN
        IF (allocated(SAVERHS)) DEALLOCATE( SAVERHS)
        IF ((ICNTL20 .ne. 0 .OR.KEEP(111).NE.0)
     *      .and. ICNTL21.ne.0) THEN
          DEALLOCATE(RHS_MUMPS)
        ENDIF
        NULLIFY(RHS_MUMPS)
      ELSE
        IF (ASSOCIATED(RHS_MUMPS)) THEN
          DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
        END IF
      END IF
      IF (I_AM_SLAVE) THEN
        IF (ALLOCATED(SRW3)) DEALLOCATE(SRW3)
        IF (ALLOCATED(POSINRHSCOMP_N)) DEALLOCATE(POSINRHSCOMP_N)
        IF (LSCAL .AND. ICNTL21==1) THEN
          DEALLOCATE(scaling_data%SCALING_LOC)
        ENDIF
        IF (MAXSUPDATED)
     &    MAXS= MAXSSAVED ! restore original value of MAXS
        IF (ASSOCIATED(id%S).AND.KEEP(201).NE.0) THEN
          DEALLOCATE(id%S)
          NULLIFY(id%S)
        ENDIF
        IF (KEEP(221).NE.1) THEN
         IF (ASSOCIATED(id%RHSCOMP)) DEALLOCATE(id%RHSCOMP)
         IF (ASSOCIATED(id%POSINRHSCOMP))
     &                 DEALLOCATE(id%POSINRHSCOMP)
        ENDIF
      ENDIF
      RETURN
 65   FORMAT (//' ERROR ANALYSIS AFTER ITERATIVE REFINEMENT')
 100  FORMAT(//' ****** SOLVE & CHECK STEP ********'/)
 110  FORMAT (//' VECTOR SOLUTION FOR COLUMN ',I12)
 115  FORMAT(1X, A44,1P,D9.2)
 120  FORMAT(//' LEAVING SOLVER WITH:  INFO(1) ............ =',I12/
     *         '                       INFO(2) ............ =',I12)
 150  FORMAT (/' STATISTICS PRIOR SOLVE PHASE     ...........'/
     *        ' NUMBER OF RIGHT-HAND-SIDES                    =',I12/
     *        ' BLOCKING FACTOR FOR MULTIPLE RHS              =',I12/
     *        ' ICNTL (9)                                     =',I12/
     *        '  --- (10)                                     =',I12/
     *        '  --- (11)                                     =',I12/
     *        '  --- (20)                                     =',I12/
     *        '  --- (21)                                     =',I12)
 151  FORMAT ('  --- (25)                                     =',I12)
 152  FORMAT ('  --- (26)                                     =',I12)
 160  FORMAT (' RHS'/(1X,1P,5D14.6))
 170  FORMAT (//' ERROR ANALYSIS' )
 240  FORMAT (1X, A42,I4)
 270  FORMAT (//' BEGIN ITERATIVE REFINEMENT' )
  81  FORMAT (/' STATISTICS AFTER ITERATIVE REFINEMENT ')
 131  FORMAT (/' END   ITERATIVE REFINEMENT ')
 141  FORMAT(1X, A42,I4)
      END SUBROUTINE ZMUMPS_301
      SUBROUTINE ZMUMPS_245(root, N, A, LA, IW, LIW, W, LWC, 
     * IWCB,LIWW,RHS,LRHS,NRHS,NA,LNA,NE_STEPS, W2,
     * MTYPE, ICNTL,
     * STEP, FRERE, DAD, FILS, PTRIST, PTRFAC, IW1,LIW1,
     * PROCNODE_STEPS, SLAVEF,
     * INFO, KEEP,KEEP8, COMM, COMM_NODES, MYID,
     * MYID_NODES,
     * BUFR, LBUFR, LBUFR_BYTES,
     * 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * IBEG_ROOT_DEF, IEND_ROOT_DEF,
     * IROOT_DEF_RHS_COL1, IPT_RHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     * RHSCOMP, LRHSCOMP, POSINRHSCOMP, BUILD_POSINRHSCOMP
     * )
      USE ZMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      TYPE ( ZMUMPS_ROOT_STRUC ) :: root
      INTEGER LA,LWC,N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(40),INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28)),
     * PTRFAC(KEEP(28)), DAD(KEEP(28))
      INTEGER LRHS, NRHS, LRHSCOMP
      COMPLEX*16    A(LA), W(LWC), RHS(LRHS,NRHS),
     *        W2(KEEP(133)), 
     *        RHSCOMP(LRHSCOMP,NRHS)
      INTEGER SLAVEF, COMM, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28)), POSINRHSCOMP(KEEP(28))
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1
      INTEGER NRHS_LOC
      INTEGER IPT_RHS_ROOT, SIZE_ROOT, MASTER_ROOT
      LOGICAL BUILD_POSINRHSCOMP
      INTEGER MP, LP, LDIAG, LWC_LOC
      INTEGER K,I
      INTEGER LPOOL,MYLEAF
      INTEGER NSTK_S,IPOOL,PTRICB,PTRACB
      INTEGER POOLSS, MTYPE_LOC
      INTEGER IDUMMY(1), LIDUMMY, DUMMY, WHAT
      INTEGER IPT_RHS_ROOT_LOC
      INTEGER IERR
      INTEGER       IOLDPS, IAPOS,
     *              LOCAL_M     ,
     *              LOCAL_N
      LOGICAL IS_NODE_IN_MEM
#if defined(V_T)
      INTEGER soln_c_class, forw_soln, back_soln, root_soln
#endif
      INTEGER INODE, IPOS, LIELL, NPIV,J1,JJ
      INTEGER IZERO
      COMPLEX*16 ZERO
      LOGICAL DOFORWARD, DOROOT, DOBACKWARD
      PARAMETER (IZERO = 0 )
      PARAMETER (ZERO = 0.0D0)
      INCLUDE 'mumps_headers.h'
      EXTERNAL ZMUMPS_248, ZMUMPS_249
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      MYLEAF = -1
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln_c',soln_c_class,ierr)
      CALL VTFUNCDEF( 'forw_soln',soln_c_class,forw_soln,ierr)
      CALL VTFUNCDEF( 'back_soln',soln_c_class,back_soln,ierr)
      CALL VTFUNCDEF( 'root_soln',soln_c_class,root_soln,ierr)
#endif
      NSTK_S   = 1
      LIDUMMY  = 1   ! size of IDUMMY array
      DUMMY    = -9998 ! should not be accessed
      PTRICB = NSTK_S + KEEP(28)
      PTRACB = PTRICB + KEEP(28)
      IPOOL  = PTRACB + KEEP(28)
      LPOOL  = LIW1 - IPOOL + 1
      IF ( (LPOOL .ne. KEEP(28) + 1).AND.(MP.GT.0) ) THEN
        WRITE(MP,*) ' warning solve: lpool.NE.nsteps+1 :',
     &       lpool,KEEP(28)+1
      ENDIF
      LWC_LOC = LWC
      IF ( MASTER_ROOT .EQ. MYID_NODES ) THEN
        LWC_LOC = IPT_RHS_ROOT - 1
      ELSE
        LWC_LOC = LWC
      ENDIF
      IF ( KEEP(60).EQ.0 .AND.
     &    ( 
     &      (KEEP(38).NE.0 .AND.  root%yes) 
     &  .OR.
     &      (KEEP(20).NE.0 .AND. MYID_NODES.EQ.MASTER_ROOT)) 
     &   ) 
     &THEN
        DOROOT    = .TRUE.
      ELSE
        DOROOT = .FALSE.
      ENDIF
      DOFORWARD = .TRUE.
      DOBACKWARD= .TRUE.
      IF ( KEEP(111).NE.0 ) THEN
        DOFORWARD = .FALSE.
      ENDIF
      IF (KEEP(221).eq.1) DOBACKWARD = .FALSE.
      IF (KEEP(221).eq.2) DOFORWARD  = .FALSE.
      IF (KEEP(201).NE.0) THEN
        IF (DOFORWARD .OR. DOROOT) THEN
           CALL ZMUMPS_583(PTRFAC,KEEP(28),MTYPE)
          CALL ZMUMPS_585(A,LA,PTRFAC,KEEP(28),IERR)
          IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            RETURN
          ENDIF
        ENDIF
      ENDIF
      IF (DOFORWARD) THEN
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 1
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(forw_soln,ierr)
#endif
        CALL ZMUMPS_248(N, A(1), LA, IW(1), LIW, W(1),
     *           LWC_LOC, RHS, LRHS, NRHS,
     *           IW1(PTRICB), IWCB, LIWW,
     *           RHSCOMP,LRHSCOMP,POSINRHSCOMP,BUILD_POSINRHSCOMP,
     *           NE_STEPS, NA, LNA, STEP, FRERE,DAD,FILS,
     *           IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     *           MYLEAF,INFO,
     *           KEEP,KEEP8,
     *           PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           W( IPT_RHS_ROOT ), MTYPE_LOC, 
     * 
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
         BUILD_POSINRHSCOMP = .FALSE.
#if defined(V_T)
        CALL VTEND(forw_soln,ierr)
#endif
      ENDIF
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     *    ': ** ERROR RETURN FROM ZMUMPS_248,INFO(1:2)=',
     *    INFO(1:2)
        END IF
        RETURN
      END IF
      CALL MPI_BARRIER( COMM_NODES, IERR )
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        IF ( KEEP(60) == 0 ) THEN
          IF ( root%yes ) THEN
        IOLDPS = PTRIST(STEP(KEEP(38)))
        LOCAL_M = IW( IOLDPS + 2 + XSIZE)
        LOCAL_N = IW( IOLDPS + 1 + XSIZE)
         IF (KEEP(201).NE.0) THEN
          IS_NODE_IN_MEM=
     $       ZMUMPS_SOLVE_IS_INODE_IN_MEM(KEEP(38),
     $       PTRFAC,KEEP(28),A,LA,IERR)
          IF(IERR.LT.0)THEN
             INFO(1)=IERR
             INFO(2)=0
             RETURN
          ENDIF
          IF(.NOT.IS_NODE_IN_MEM) THEN
           CALL ZMUMPS_578(KEEP(38),PTRFAC,
     $          KEEP,KEEP8,A,IERR)
            IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
            ENDIF
            CALL ZMUMPS_577(
     &           A(PTRFAC(IW(IOLDPS+4+ XSIZE))),KEEP(38),IERR
     $     )
            IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              RETURN
            ENDIF
          ENDIF
         ENDIF
         IAPOS   = PTRFAC(IW( IOLDPS + 4 + XSIZE))
#if defined(V_T)
      CALL VTBEGIN(root_soln,ierr)
#endif
          CALL DESCINIT( root%DESCB, root%TOT_ROOT_SIZE,
     *      NRHS, root%MBLOCK, root%NBLOCK, 0, 0,
     *      root%CNTXT_BLACS, LOCAL_M, IERR )
          IF (IERR.NE.0) THEN
            WRITE(*,*) 'After DESCINIT, IERR = ', IERR
            CALL ZMUMPS_ABORT()
          END IF
#if defined(null_space_old)
          CALL ZMUMPS_352( NRHS, root%DESCRIPTOR,
     *       root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50), KEEP(19),
     *       root%QR_TAU, W(1), LWC_LOC, KEEP(17),
     *       root%MAXG, root%GIND, root%GROW, root%GCOS, root%GSIN )
#else
          CALL ZMUMPS_286( NRHS, root%DESCRIPTOR, root%DESCB,
     *       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     *       root%MBLOCK, root%NBLOCK,
     *       root%IPIV, root%LPIV, MASTER_ROOT, MYID_NODES,
     *       COMM_NODES,
     *       W( IPT_RHS_ROOT ),
     *       root%TOT_ROOT_SIZE, A( IAPOS ),
     *       INFO(1), MTYPE, KEEP(50))
#endif
          IF(KEEP(201).NE.0)THEN
             CALL ZMUMPS_598(KEEP(38),
     $             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
             IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 RETURN
             ENDIF
          ENDIF
        ENDIF  ! of root%yes
      ELSE     ! of KEEP(60).eq.0
        IF ( 
     *       (KEEP(221).EQ.0) .AND. 
     *       ( MYID_NODES .eq.  ZMUMPS_275( STEP(KEEP(38)),
     *         PROCNODE_STEPS, SLAVEF ) ) 
     *     )  THEN
           W( IPT_RHS_ROOT:IPT_RHS_ROOT+NRHS*SIZE_ROOT- 1) = 0.0D0
        ENDIF
      ENDIF
      ELSE IF ( KEEP(20) .NE. 0 ) THEN
        IF ( MYID_NODES .eq.  ZMUMPS_275( STEP(KEEP(20)),
     *        PROCNODE_STEPS, SLAVEF ) ) THEN
#if defined(try_null_space)
         IF ( KEEP(60) .eq. 0 ) THEN
            IF (KEEP(201).NE.0) THEN
              IS_NODE_IN_MEM=ZMUMPS_SOLVE_IS_INODE_IN_MEM(KEEP(20),
     *        PTRFAC,KEEP(28),A,LA,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                RETURN
              ENDIF
              IF ( .NOT. IS_NODE_IN_MEM ) THEN
                CALL ZMUMPS_578(KEEP(20),PTRFAC,
     $               KEEP,KEEP8,A,IERR)
                IF(IERR.LT.0)THEN
                  INFO(1)=IERR
                  INFO(2)=0
                  RETURN
                ENDIF
                CALL ZMUMPS_577(
     &            A(PTRFAC(STEP(KEEP(20)))),KEEP(20),IERR
     $           )
                IF(IERR.LT.0)THEN
                  INFO(1)=IERR
                  INFO(2)=0
                  RETURN
                ENDIF
              ENDIF
            END IF
            NRHS_LOC         = NRHS
            IPT_RHS_ROOT_LOC = IPT_RHS_ROOT
            IF ( KEEP(111).NE.0 ) THEN
              W( IPT_RHS_ROOT:IPT_RHS_ROOT+NRHS*SIZE_ROOT- 1) = 0.0D0
              NRHS_LOC = IEND_ROOT_DEF - IBEG_ROOT_DEF + 1
              IPT_RHS_ROOT_LOC = IPT_RHS_ROOT_LOC +
     *                           (IROOT_DEF_RHS_COL1-1)*SIZE_ROOT
            ENDIF
           IF (NRHS_LOC .GT. 0) THEN
             CALL ZMUMPS_296(NRHS_LOC,
     *            SIZE_ROOT,A( PTRFAC(
     *            IW( PTRIST(STEP(KEEP(20)))+4+XSIZE))),
     *            root, IBEG_ROOT_DEF, IEND_ROOT_DEF,
     *            W( IPT_RHS_ROOT_LOC ), W(1),LWC_LOC,KEEP,KEEP8,
     *            MTYPE,INFO)
           ENDIF
           IF(KEEP(201).NE.0)THEN
             CALL ZMUMPS_598(KEEP(20),
     $             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
             IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 RETURN
             ENDIF
           ENDIF
         ELSE
#endif
           IF (KEEP(221).EQ.0) 
     &     W( IPT_RHS_ROOT:IPT_RHS_ROOT+NRHS*SIZE_ROOT- 1) = 0.0D0
#if defined(try_null_space)
         END IF
#endif
        END IF
      END IF
#if defined(V_T)
      CALL VTEND(root_soln,ierr)
#endif
      IF ( INFO(1) .LT. 0 ) RETURN
      IF (DOBACKWARD) THEN
        IF (BUILD_POSINRHSCOMP) THEN
          WHAT = 0   ! only compute POSINRHSCOMP
          CALL ZMUMPS_639
     *           (SLAVEF, N, MYID_NODES,
     *           PTRIST,
     *           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     *           POSINRHSCOMP, IDUMMY, LIDUMMY, DUMMY, WHAT)
          BUILD_POSINRHSCOMP=.FALSE.  
        ENDIF
        IF(KEEP(201).NE.0)THEN
          CALL ZMUMPS_584(PTRFAC,KEEP(28),MTYPE)
          IF (KEEP(38).NE.0.AND.root%yes)THEN
             IF(KEEP(60).EQ.0)THEN ! Schur is off              
                CALL ZMUMPS_598(KEEP(38),
     $             PTRFAC,KEEP(28),A,LA,.TRUE.,IERR)
             ENDIF
          ENDIF
          IF (KEEP(20).NE.0 .and. MYID_NODES.EQ.MASTER_ROOT.AND.
     $        KEEP(19).NE.0 ) THEN
            CALL ZMUMPS_598(KEEP(20),
     $             PTRFAC,KEEP(28),A,LA,.TRUE.,IERR)
          ENDIF
        ENDIF
        IF (KEEP(201).NE.0) THEN
          IF (.NOT.DOFORWARD) THEN
            CALL ZMUMPS_585(A,LA,PTRFAC,KEEP(28),IERR)
          ENDIF
        ENDIF
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = IZERO
        ENDIF
        IF (.NOT. DOFORWARD) THEN
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(back_soln,ierr)
#endif
        CALL ZMUMPS_249( N, A, LA, IW, LIW, W(1), LWC_LOC,
     *              RHS, LRHS, NRHS,
     *              RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     *              IW1(PTRICB),IW1(PTRACB),IWCB,LIWW,
     *              W2, NE_STEPS, NA, LNA, STEP, FRERE,FILS,
     *              IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,INFO,
     *              PROCNODE_STEPS, SLAVEF, COMM_NODES,MYID_NODES,
     *              BUFR, LBUFR, LBUFR_BYTES, KEEP,KEEP8,
     *              W( IPT_RHS_ROOT ),
     *              MTYPE_LOC, 
     *              ISTEP_TO_INIV2, TAB_POS_IN_PERE )
#if defined(V_T)
      CALL VTEND(back_soln,ierr)
#endif
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        IF (DOFORWARD) THEN
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        WRITE (MP,99992)
        IF (N.GT.0) WRITE (MP,99993) (RHS(I,1),I=1,K)
        IF (N.GT.0.and.NRHS>1) 
     *              WRITE (MP,99994) (RHS(I,2),I=1,K)
        ENDIF
      ENDIF
      RETURN
99993 FORMAT (' RHS    (first column)'/(1X,1P,5D14.6))
99994 FORMAT (' RHS    (2 nd  column)'/(1X,1P,5D14.6))
99992 FORMAT (//' LEAVING SOLVE (MPI41C) WITH')
      END SUBROUTINE ZMUMPS_245
      SUBROUTINE ZMUMPS_521(NSLAVES, N, MYID, COMM,
     *           MTYPE, RHS, LRHS, NRHS, PTRIST,
     *           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, BUFFER,
     *           SIZE_BUF, SIZE_BUF_BYTES )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS
      COMPLEX*16 RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF_BYTES)
      INTEGER I, J, JJ, J1, ISTEP, MASTER,
     *        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV
      INTEGER MSGSOU, STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND, IROW
      INTEGER SK38, SK20
      COMPLEX*16 ONE_ROW(NRHS)
      INCLUDE 'mumps_headers.h'
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
       IF (NSLAVES.EQ.1 .AND. KEEP(46).EQ.1) RETURN
       IF (NSLAVES.EQ.1 .AND. KEEP(46).EQ.0) THEN
         DO J=1, NRHS
           IF ( MYID .EQ. 1 ) THEN
             CALL MPI_SEND(RHS(1, J), N, MPI_DOUBLE_COMPLEX, MASTER,
     *                 GatherSol, COMM, IERR)
     * 
           ELSE
             CALL MPI_RECV(RHS(1, J), N, MPI_DOUBLE_COMPLEX,
     *                 1,
     *                 GatherSol, COMM, STATUS, IERR )
           ENDIF
         ENDDO
         RETURN
       ENDIF
      N2SEND=0
      N2RECV=N
      POS_BUF=0
      TYPE_PARAL = KEEP(46)
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      IF (NSLAVES > 1 .OR. TYPE_PARAL == 0) THEN
        CALL MPI_PACK_SIZE(2,MPI_INTEGER, COMM, SIZE1, IERR)
        CALL MPI_PACK_SIZE(NRHS,MPI_DOUBLE_COMPLEX, COMM,
     *                   SIZE2, IERR)
        RECORD_SIZE_P_1= SIZE1+SIZE2
      ELSE
        RECORD_SIZE_P_1 = -9999
      ENDIF
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      IF (I_AM_SLAVE) THEN
        POS_BUF = 0
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == ZMUMPS_275(ISTEP,
     *          PROCNODE_STEPS,NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                    IPOS = PTRIST(ISTEP) 
                    LIELL = IW(IPOS+3+XSIZE)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              IF (MYID .EQ. MASTER) THEN
               N2RECV=N2RECV-NPIV
              ELSE
               DO JJ=J1,J1+NPIV-1
                CALL ZMUMPS_522( IW(JJ), RHS(IW(JJ),1:NRHS) )
               ENDDO
              ENDIF
          ENDIF
        ENDDO
        CALL ZMUMPS_523()   ! Send remaining rows
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     *                 MPI_ANY_SOURCE,
     *                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     *                   IROW, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (IROW.NE.0)
          CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     *                   ONE_ROW, NRHS, MPI_DOUBLE_COMPLEX,
     *                   COMM, IERR)
          RHS(IROW,1:NRHS)=ONE_ROW
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     *                   IROW, 1, MPI_INTEGER, COMM, IERR)
          N2RECV=N2RECV-1
        ENDDO
       ENDDO
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE ZMUMPS_522(IROW, RHS_VAL)
        INTEGER IROW
        COMPLEX*16 RHS_VAL(NRHS)
        CALL MPI_PACK(IROW, 1, MPI_INTEGER, BUFFER,
     *                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(RHS_VAL, NRHS, MPI_DOUBLE_COMPLEX,
     *                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     *                IERR)
        N2SEND=N2SEND+1
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL ZMUMPS_523()
        END IF
        RETURN
        END SUBROUTINE ZMUMPS_522
        SUBROUTINE ZMUMPS_523()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(0, 1, MPI_INTEGER, BUFFER,
     *                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     *                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE ZMUMPS_523
      END SUBROUTINE ZMUMPS_521
      SUBROUTINE ZMUMPS_535(MTYPE, ISOL_LOC,
     *             LSOL_LOC, PTRIST, KEEP,KEEP8,
     *             IW, LIW_PASSED, MYID_NODES, N, STEP,
     *             PROCNODE, NSLAVES, scaling_data, LSCAL)
      IMPLICIT NONE
      INTEGER MTYPE, LSOL_LOC, MYID_NODES, N, NSLAVES
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE(KEEP(28))
      INTEGER ISOL_LOC(LSOL_LOC)
      INTEGER LIW_PASSED
      INTEGER IW(LIW_PASSED)
      INTEGER STEP(N)
      LOGICAL LSCAL
      type scaling_data_t
        SEQUENCE
        COMPLEX*16, dimension(:), pointer :: SCALING
        COMPLEX*16, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      INTEGER ISTEP, K
      INTEGER J1, IPOS, LIELL, NPIV, JJ
      LOGICAL ROOT
      INTEGER SK38,SK20
      INCLUDE 'mumps_headers.h'
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      K=0
      DO ISTEP=1, KEEP(28)
          IF ( MYID_NODES == ZMUMPS_275( ISTEP,
     *         PROCNODE, NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                    IPOS = PTRIST(ISTEP)+XSIZE
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 + XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 + XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                  K=K+1
                  ISOL_LOC(K)=IW(JJ)
                  IF (LSCAL) THEN
                    scaling_data%SCALING_LOC(K)=
     *              scaling_data%SCALING(IW(JJ))
                  ENDIF
              ENDDO
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_535
      SUBROUTINE ZMUMPS_532(
     *           SLAVEF, N, MYID_NODES,
     *           MTYPE, RHS, LD_RHS, NRHS,
     *           ISOL_LOC, SOL_LOC, BEG_RHS, LSOL_LOC,
     *           PTRIST,
     *           PROCNODE_STEPS, KEEP,KEEP8, IW, LIW, STEP,
     *           scaling_data, LSCAL)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      type scaling_data_t
        SEQUENCE
        COMPLEX*16, dimension(:), pointer :: SCALING
        COMPLEX*16, dimension(:), pointer :: SCALING_LOC
      end type scaling_data_t
      TYPE (scaling_data_t) :: scaling_data
      LOGICAL LSCAL
      INTEGER SLAVEF, N, MYID_NODES, LIW, MTYPE, NRHS, LD_RHS
      INTEGER LSOL_LOC, BEG_RHS
      INTEGER ISOL_LOC(LSOL_LOC)
      COMPLEX*16 SOL_LOC( LSOL_LOC, *)
      COMPLEX*16 RHS(  LD_RHS , NRHS)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER JJ, J1, ISTEP, K
      INTEGER IPOS, LIELL, NPIV
      LOGICAL ROOT
      INCLUDE 'mumps_headers.h'
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      K=0
        DO ISTEP = 1, KEEP(28)
            IF (MYID_NODES == ZMUMPS_275(ISTEP,
     *          PROCNODE_STEPS,SLAVEF)) THEN
              ROOT=.false.
              IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
              IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
              IF ( ROOT ) THEN
                    IPOS = PTRIST(ISTEP) + XSIZE
                    LIELL = IW(IPOS+3)
                    NPIV = LIELL
                    IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2 +XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +XSIZE)
              END IF
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                K=K+1
                IF (LSCAL) THEN
                  SOL_LOC(K,BEG_RHS:BEG_RHS+NRHS-1) =
     *            scaling_data%SCALING_LOC(K)*RHS(IW(JJ),1:NRHS)
                ELSE
                  SOL_LOC(K,BEG_RHS:BEG_RHS+NRHS-1) =
     *            RHS(IW(JJ),1:NRHS)
                ENDIF
              ENDDO
            ENDIF
        ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_532
      SUBROUTINE ZMUMPS_638
     *           (NSLAVES, N, MYID, COMM,
     *           MTYPE, RHS, LRHS, NRHS, PTRIST,
     *           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     *           POSINRHSCOMP, LENPOSINRHSCOMP,
     *           BUILD_POSINRHSCOMP, ICNTL, INFO)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS, LENPOSINRHSCOMP
      INTEGER ICNTL(40), INFO(40)
      COMPLEX*16 RHS   (LRHS, NRHS)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N), POSINRHSCOMP(LENPOSINRHSCOMP)
      LOGICAL BUILD_POSINRHSCOMP
      INTEGER BUF_MAXSIZE
      PARAMETER (BUF_MAXSIZE=2000)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF_INDX
      COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: BUF_RHS
      INTEGER ENTRIES_2_PROCESS, PROC_WHO_ASKS, BUF_EFFSIZE
      INTEGER INDX ! Temporary position in 1..N of a row in RHS
      INTEGER allocok
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0.0D0)
      INTEGER I, J, K, JJ, J1, ISTEP, MASTER,
     *        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV
      INTEGER MSGSOU, STATUS(MPI_STATUS_SIZE), IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER SK38, SK20, IPOSINRHSCOMP
      INCLUDE 'mumps_headers.h'
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      TYPE_PARAL = KEEP(46)
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      BUF_EFFSIZE = 0
      ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     *          BUF_RHS(NRHS,BUF_MAXSIZE),
     *          stat=allocok)
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)=BUF_MAXSIZE*(NRHS+1)
      ENDIF
      CALL ZMUMPS_276(ICNTL, INFO, COMM, MYID )
      IF (INFO(1).LT.0) RETURN
      IF (MYID.EQ.MASTER) THEN
        ENTRIES_2_PROCESS = N - KEEP(89)
        DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
          CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     *                 MPI_ANY_SOURCE,
     *                 ScatterRhsI, COMM, STATUS, IERR )
          CALL MPI_GET_COUNT( STATUS, MPI_INTEGER, BUF_EFFSIZE, IERR )
          PROC_WHO_ASKS = STATUS(MPI_SOURCE)
          DO I = 1, BUF_EFFSIZE
              INDX = BUF_INDX( I )
            DO K = 1, NRHS
              BUF_RHS( K, I ) = RHS( INDX, K )
              RHS( BUF_INDX(I), K ) = DCMPLX( ZERO )
            ENDDO
          ENDDO
          CALL MPI_SEND( BUF_RHS, NRHS*BUF_EFFSIZE,
     *                   MPI_DOUBLE_COMPLEX, PROC_WHO_ASKS,
     *                   ScatterRhsR, COMM, IERR)
          ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
        ENDDO
        BUF_EFFSIZE= 0  
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (BUILD_POSINRHSCOMP) THEN
           IPOSINRHSCOMP = 1     ! to update POSINRHSCOMP on the fly
           POSINRHSCOMP = -9678  ! to debug  FIXME: SUPPRESS
        ENDIF
        IF (MYID.NE.MASTER) RHS = DCMPLX(ZERO)
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == ZMUMPS_275(ISTEP,
     *          PROCNODE_STEPS,NSLAVES)) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                  IPOS = PTRIST(ISTEP) 
                  LIELL = IW(IPOS+3+XSIZE)
                  NPIV = LIELL
                  IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +XSIZE)
              END IF
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              IF (BUILD_POSINRHSCOMP) THEN
                 POSINRHSCOMP(ISTEP) = IPOSINRHSCOMP
                 IPOSINRHSCOMP       = IPOSINRHSCOMP + NPIV
              ENDIF
              IF (MYID.NE.MASTER) THEN
                DO JJ=J1,J1+NPIV-1
                  BUF_EFFSIZE = BUF_EFFSIZE + 1
                  BUF_INDX(BUF_EFFSIZE) = IW(JJ)
                  IF (BUF_EFFSIZE + 1 .GT. BUF_MAXSIZE) THEN
                   CALL ZMUMPS_640()
                  ENDIF
                ENDDO
              ENDIF
          ENDIF
        ENDDO
        IF ( BUF_EFFSIZE .NE. 0 .AND. MYID.NE.MASTER ) 
     &              CALL ZMUMPS_640()
      ENDIF
      DEALLOCATE (BUF_INDX, BUF_RHS)
      RETURN
      CONTAINS
                  SUBROUTINE ZMUMPS_640()
                  CALL MPI_SEND(BUF_INDX, BUF_EFFSIZE, MPI_INTEGER,
     *            MASTER, ScatterRhsI, COMM, IERR )
                  CALL MPI_RECV(BUF_RHS, BUF_EFFSIZE*NRHS,
     *                 MPI_DOUBLE_COMPLEX,
     *                 MASTER,
     *                 ScatterRhsR, COMM, STATUS, IERR )
                  DO I = 1, BUF_EFFSIZE
                    INDX = BUF_INDX(I)
                    DO K = 1, NRHS
                      RHS( INDX, K ) = BUF_RHS( K, I )
                    ENDDO
                  ENDDO
                  BUF_EFFSIZE = 0
                  RETURN
                  END SUBROUTINE ZMUMPS_640
      END SUBROUTINE ZMUMPS_638
      SUBROUTINE ZMUMPS_639
     *           (NSLAVES, N, MYID_NODES,
     *           PTRIST,
     *           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     *           POSINRHSCOMP, POSINRHSCOMP_N, LPIRC_N, MTYPE,
     *           WHAT )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID_NODES, COMM, LIW
      INTEGER LENPOSINRHSCOMP
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N), POSINRHSCOMP(KEEP(28))
      INTEGER LPIRC_N, WHAT, MTYPE
      INTEGER POSINRHSCOMP_N(LPIRC_N)
      INTEGER ISTEP
      INTEGER NPIV
      INTEGER SK38, SK20, IPOS, LIELL
      INTEGER JJ, J1
      INTEGER IPOSINRHSCOMP
      INCLUDE 'mumps_headers.h'
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      IF (WHAT .NE. 0.AND. WHAT.NE.1.AND.WHAT.NE.2) THEN
        WRITE(*,*) "Internal error in ZMUMPS_639"
        CALL ZMUMPS_ABORT()
      ENDIF
      IF (KEEP(38).NE.0) THEN
        SK38=STEP(KEEP(38))
      ELSE
        SK38=0
      ENDIF
      IF (KEEP(20).NE.0) THEN 
        SK20=STEP(KEEP(20))
      ELSE
        SK20=0
      ENDIF
      IPOSINRHSCOMP   = 1     ! to update POSINRHSCOMP on the fly
      POSINRHSCOMP = -9678  ! to debug  FIXME:SUppress
      IF (WHAT .NE. 0) THEN
        POSINRHSCOMP_N = 0 ! initialized to 0
      ENDIF
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == ZMUMPS_275(ISTEP,
     *     PROCNODE_STEPS,NSLAVES)) THEN
           IPOS = PTRIST(ISTEP)
           NPIV = IW(IPOS+3+XSIZE)
           POSINRHSCOMP(ISTEP) = IPOSINRHSCOMP
           IF (WHAT .NE. 0) THEN
              IF ( ISTEP.EQ.SK38 .OR. ISTEP.EQ.SK20 ) THEN
                  IPOS = PTRIST(ISTEP)
                  LIELL = IW(IPOS+3+XSIZE)
                  NPIV = LIELL
                  IPOS= PTRIST(ISTEP)+5+XSIZE
              ELSE
                  IPOS = PTRIST(ISTEP) + 2+ XSIZE
                  LIELL = IW(IPOS-2)+IW(IPOS+1)
                  IPOS= IPOS+1
                  NPIV = IW(IPOS)
                  IPOS= IPOS+1
                  IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +XSIZE)
              ENDIF
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              DO JJ = J1, J1+NPIV-1
                POSINRHSCOMP_N(IW(JJ)) = IPOSINRHSCOMP+JJ-J1
              END DO
           ENDIF
           IPOSINRHSCOMP       = IPOSINRHSCOMP + NPIV
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ZMUMPS_639
      SUBROUTINE ZMUMPS_248(N, A, LA, IW, LIW, WCB, LWCB,
     *    RHS, LRHS, NRHS,
     *    PTRICB, IWCB, LIWCB, 
     *    RHSCOMP, LRHSCOMP, POSINRHSCOMP, BUILD_POSINRHSCOMP,
     *    NE_STEPS, NA, LNA, STEP,
     *    FRERE, DAD, FILS,
     *    NSTK_S, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO,
     *    KEEP,KEEP8,
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE ZMUMPS_OOC
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N, LA, LIW, LWCB, LPOOL, LIWCB, LNA
      INTEGER SLAVEF, MYLEAF, COMM, MYID
      INTEGER INFO( 40 ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER LRHS, NRHS
      COMPLEX*16 A( LA ), RHS( LRHS, NRHS ), WCB( LWCB )
      COMPLEX*16 RHS_ROOT( * )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NA( LNA ), NE_STEPS( KEEP(28) )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N ),
     *        DAD( KEEP(28) )
      INTEGER NSTK_S(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER PTRICB( KEEP(28) ) 
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER POSINRHSCOMP(KEEP(28)), LRHSCOMP 
      LOGICAL BUILD_POSINRHSCOMP
      COMPLEX*16 RHSCOMP( LRHSCOMP, NRHS )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGTAG, MSGSOU, DUMMY(1)
      LOGICAL FLAG
      INTEGER NBFIN, MYROOT
      INTEGER POSIWCB,POSWCB,PLEFTWCB
      INTEGER INODE
      INTEGER RHSCOMPFREEPOS
      INTEGER I, K
      INTEGER III, NBROOT,NBLEAF,LEAF
      LOGICAL BLOQ
      EXTERNAL ZMUMPS_275
      INTEGER ZMUMPS_275
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1
      IF (BUILD_POSINRHSCOMP) RHSCOMPFREEPOS= 1
      DO I = 1, KEEP(28)
        NSTK_S(I)   = NE_STEPS(I)
      ENDDO
      PTRICB = 0
      CALL ZMUMPS_362(N, LEAF, NBROOT, MYROOT, MYID,
     *     SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     *     PROCNODE_STEPS, IPOOL, LPOOL)
      NBFIN = SLAVEF
      IF ( MYROOT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        DUMMY(1) = 1
        CALL ZMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *       RACINE_SOLVE, SLAVEF)
      END IF
      MYLEAF = LEAF - 1
      III    = 1
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
         CALL ZMUMPS_574
     &        ( IPOOL(1), LPOOL, III, LEAF, INODE,
     &          KEEP(208), MYID )
        GOTO 60
      ENDIF
      BLOQ = ( ( III .EQ. LEAF )
     $     )
      CALL ZMUMPS_303( BLOQ, FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, NRHS, IPOOL, LPOOL, III, LEAF,
     *     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *     IWCB, LIWCB,
     *     WCB, LWCB, POSWCB,
     *     PLEFTWCB, POSIWCB,
     *     PTRICB, INFO, KEEP,KEEP8, STEP,
     *     PROCNODE_STEPS,
     *     RHS, LRHS)
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
         IF (III .NE. LEAF) THEN
            CALL ZMUMPS_574
     &           (IPOOL(1), LPOOL, III, LEAF, INODE,
     &           KEEP(208), MYID )
            GOTO 60
         ENDIF                  ! POOL not empty
      ENDIF                     !   .not.(FLAG) = not_msg_received
      GOTO 50
 60   CONTINUE
      CALL ZMUMPS_302( INODE, BUFR, LBUFR, LBUFR_BYTES,
     *        MSGTAG, MSGSOU, MYID, SLAVEF, COMM,  N,
     *        IPOOL, LPOOL, III, LEAF, NBFIN, NSTK_S,
     *        IWCB, LIWCB, WCB, LWCB, A, LA,
     *        IW, LIW, RHS, LRHS, NRHS, 
     *        POSWCB, PLEFTWCB, POSIWCB,
     *        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     *        FILS, STEP, FRERE, DAD,
     *        MYROOT, INFO, KEEP,KEEP8, RHS_ROOT, MTYPE, 
     *        RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     *        RHSCOMPFREEPOS, BUILD_POSINRHSCOMP,
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *      )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      CALL ZMUMPS_150( MYID,COMM,BUFR,
     *                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE ZMUMPS_248
      RECURSIVE SUBROUTINE ZMUMPS_323
     *         ( BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, NRHS, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     *           PTRFAC, IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB,
     *           INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS, 
     *           RHS, LRHS)
      USE ZMUMPS_OOC 
      USE ZMUMPS_BUFFER 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW, LA
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28)),PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      COMPLEX*16 WCB( LWCB ), A( LA )
      INTEGER LRHS
      COMPLEX*16 RHS(LRHS, NRHS)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, K, JJ
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PTRX, PTRY, APOS, PDEST, I
      LOGICAL FLAG
      EXTERNAL ZMUMPS_275
      INTEGER  ZMUMPS_275
      COMPLEX*16 ALPHA, ONE
      PARAMETER( ONE = 1.0D0, ALPHA = -1.0D0 )
      INCLUDE 'mumps_headers.h'
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
          NBFIN = NBFIN - 1
          IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         LONG, 1, MPI_INTEGER, COMM, IERR )
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
               IPOOL( LEAF ) = FPERE
               LEAF = LEAF + 1
               IF ( LEAF > LPOOL ) THEN
                 WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                 CALL ZMUMPS_ABORT()
               END IF
             END IF
          ELSE
          IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
             PTRICB(STEP(FINODE)) = NCB + 1
          END IF
          IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = LONG
            GOTO 260
          END IF
          IF ( POSWCB - PLEFTWCB + 1 .LT. LONG * NRHS) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = PLEFTWCB - POSWCB - 1 + LONG * NRHS
            GOTO 260
          END IF
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          IWCB( 1 ),
     *          LONG, MPI_INTEGER, COMM, IERR )
           DO K = 1, NRHS
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *               WCB( PLEFTWCB ),
     *               LONG, MPI_DOUBLE_COMPLEX, COMM, IERR )
            DO I = 1, LONG
              RHS(IWCB(I),K) = RHS(IWCB(I),K) + WCB(PLEFTWCB+I-1)
            ENDDO
           END DO
           PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
          ENDIF
          IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          END IF
          IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
            IF ( LEAF > LPOOL ) THEN
              WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
              CALL ZMUMPS_ABORT()
            END IF
          ENDIF
        END IF
        ELSEIF ( MSGTAG .EQ. MASTER2SLAVE ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         NPIV, 1, MPI_INTEGER, COMM, IERR )
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + NCV * NRHS
          PLEFTWCB = PLEFTWCB + (NPIV + NCV) * NRHS
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
            INFO(1) = -11
            INFO(2) = -POSWCB + PLEFTWCB -1
            GO TO 260
          END IF
          DO K=1, NRHS
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *         WCB( PTRY + (K-1) * NCV ), NCV,
     *         MPI_DOUBLE_COMPLEX, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
            DO K=1, NRHS
              CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          WCB( PTRX + (K-1)*NPIV ), NPIV,
     *          MPI_DOUBLE_COMPLEX, COMM, IERR )
            END DO
          END IF
          IF (KEEP(201).NE.0) THEN
           IF(.NOT.ZMUMPS_SOLVE_IS_INODE_IN_MEM(FINODE,PTRFAC,
     $         KEEP(28),A,LA,IERR))THEN
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL ZMUMPS_578(FINODE,PTRFAC,
     $            KEEP,KEEP8,A,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL ZMUMPS_577(
     $            A(PTRFAC(IW( PTRIST( STEP(FINODE)) + 3 +XSIZE))),
     $            FINODE,IERR
     $   )
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ELSE
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ENDIF
          ENDIF
          APOS = PTRFAC(IW( PTRIST( STEP(FINODE)) + 3 + XSIZE ))
          IF ( NRHS == 1 ) THEN
            CALL ZGEMV( 'T', NPIV, NCV, ALPHA, A(APOS), NPIV,
     *                WCB( PTRX ), 1, ONE,
     *                WCB( PTRY ), 1 )
          ELSE
            CALL ZGEMM( 'T', 'N', NCV, NRHS, NPIV, ALPHA, A(APOS), NPIV,
     *                WCB( PTRX), NPIV, ONE,
     *                WCB( PTRY), NCV )
          ENDIF
          IF (KEEP(201).NE.0) THEN
             CALL ZMUMPS_598(FINODE,PTRFAC,
     $            KEEP(28),A,LA,.TRUE.,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - NPIV * NRHS
          PDEST = ZMUMPS_275( STEP(FPERE),
     *            PROCNODE_STEPS, SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
            IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
              NCB = IW( PTRIST(STEP(FINODE)) + 2 + XSIZE )
              PTRICB(STEP(FINODE)) = NCB + 1
            END IF
            DO I = 1, NCV
              JJ=IW(PTRIST(STEP(FINODE))+3+I+ XSIZE )
              DO K=1, NRHS
                RHS(JJ,K)= RHS(JJ,K) + WCB(PTRY+I-1+(K-1)*NCV)
              ENDDO
            END DO
            PTRICB(STEP(FINODE)) =
     *      PTRICB(STEP(FINODE)) - NCV
            IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                 NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            END IF
            IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
              IF ( LEAF > LPOOL ) THEN
                WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                CALL ZMUMPS_ABORT()
              END IF
            ENDIF
          ELSE
 210       CONTINUE
           CALL ZMUMPS_78( NRHS, FINODE, FPERE,
     *       IW(PTRIST(STEP( FINODE )) + 2 + XSIZE ), NCV, NCV,
     *       IW(PTRIST(STEP(FINODE))+4+ XSIZE ),
     *       WCB( PTRY ), PDEST, ContVec, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
              CALL ZMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, NRHS, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP,KEEP8, STEP,
     *               PROCNODE_STEPS, 
     *               RHS, LRHS)
              IF ( INFO( 1 )  .LT. 0 )  GOTO 270
              GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     *              NCV * KEEP( 35 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     *              NCV * KEEP( 35 )
            END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV * NRHS
        ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
        ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
        ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
        ENDIF
        GO TO 270
 260    CONTINUE
        CALL ZMUMPS_44( MYID, SLAVEF, COMM )
 270    CONTINUE
        RETURN
        END SUBROUTINE ZMUMPS_323
      SUBROUTINE ZMUMPS_302( INODE,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *           N, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S,
     *           IWCB, LIWCB,
     *           WCB, LWCB, A, LA, IW, LIW,
     *           RHS, LRHS, NRHS, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     *           FILS, STEP, FRERE, DAD,
     *           MYROOT,
     *           INFO, KEEP,KEEP8, RHS_ROOT, MTYPE,
     *           RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     *           RHSCOMPFREEPOS, BUILD_POSINRHSCOMP,
     *           
     *           ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *           )
      USE ZMUMPS_OOC
      USE ZMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LWCB, LIW, LA, POSWCB, PLEFTWCB, POSIWCB
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER LRHS, NRHS
      COMPLEX*16 WCB( LWCB ), A( LA )
      COMPLEX*16 RHS(LRHS, NRHS ), RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER POSINRHSCOMP(KEEP(28)), LRHSCOMP, RHSCOMPFREEPOS
      COMPLEX*16 RHSCOMP(LRHSCOMP, NRHS)
      LOGICAL BUILD_POSINRHSCOMP
      EXTERNAL ZGEMV, ZTRSV, ZGEMM, ZTRSM, ZMUMPS_275
      INTEGER ZMUMPS_275
      COMPLEX*16 ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB,
     *        IERR,
     *        APOS, APOS1, IF, IFR, IPOSCB, APOSCB, LIELL, IN, JJ,
     *        NELIM, PLEFT, PCB_COURANT, PPIV_COURANT
       INTEGER IPOSINRHSCOMP
      INTEGER Effective_CB_SIZE, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER APOS2,POSWCB1,POSWCB2
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      IF ( INODE .eq. KEEP(38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
        LIELL = IW( PTRIST( STEP(INODE)) + 3 + XSIZE)
        NPIV  = LIELL
        NELIM = 0
        NSLAVES = 0
        IPOS = PTRIST( STEP(INODE)) + 5 + XSIZE
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + XSIZE
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF (KEEP(201).NE.0) THEN
         IF(.NOT.ZMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,PTRFAC,
     $       KEEP(28),A,LA,IERR))THEN
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL ZMUMPS_578(INODE,PTRFAC,
     $          KEEP,KEEP8,A,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL ZMUMPS_577(A(PTRFAC(IW(IPOS))),INODE,IERR
     $   )
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ELSE
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ENDIF                    ! END INODE not in memory
        ENDIF                     ! END OOC_RUN
        APOS = PTRFAC(IW(IPOS))
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE)
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + LIELL
        J3 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + 2 * LIELL
        J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
        IFR = 0
        DO JJ = J1, J3
          J = IW( JJ )
          IFR = IFR + 1
          DO K=1,NRHS
            RHS_ROOT(IFR+NPIV*(K-1)) = RHS(J,K) 
          END DO
        END DO
        IF ( NPIV .LT. LIELL ) THEN
          WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
          CALL ZMUMPS_ABORT()
        END IF
        MYROOT = MYROOT - 1
        IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL ZMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
        END IF
        GO TO 270
      END IF
      PLEFT    = PLEFTWCB
      PPIV_COURANT = PLEFTWCB
      PCB_COURANT = PLEFT + NPIV*NRHS
      PLEFTWCB = PLEFTWCB + LIELL * NRHS
      IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
        INFO(1) = -11
        INFO(2) = PLEFTWCB - POSWCB - 1
        GO TO 260
      END IF
      IFR = PPIV_COURANT - 1
      DO 130 JJ = J1, J3
        J = IW(JJ)
        IFR = IFR + 1
        DO K=1, NRHS
          WCB(IFR+(K-1)*NPIV) = RHS(J,K) 
        END DO
  130 CONTINUE
      IFR = PCB_COURANT - 1
      IF (NPIV .LT. LIELL) THEN
        DO 140 JJ = J3 + 1, J2
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            WCB(IFR+(K-1)*NCB) = RHS(J,K)
            RHS(J,K)=ZERO
          ENDDO
  140   CONTINUE
      ENDIF
      IF ( NPIV .NE. 0 ) THEN
        IF ( KEEP(50).NE.0) THEN
          IF ( NRHS == 1 ) THEN
            CALL ZTRSV( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     *                   WCB(PPIV_COURANT), 1 )
          ELSE
            CALL ZTRSM( 'L','U','T','U', NPIV, NRHS, ONE,
     *                   A(APOS), NPIV, WCB(PPIV_COURANT),
     *                   NPIV )
          ENDIF
        ELSE
          IF ( MTYPE .eq. 1 ) THEN
            IF ( NRHS == 1)  THEN
              CALL ZTRSV( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     *        WCB(PPIV_COURANT), 1 )
            ELSE
              CALL ZTRSM( 'L','U','T','U', NPIV, NRHS, ONE,
     *                   A(APOS), LIELL, WCB(PPIV_COURANT),
     *                   NPIV )
            ENDIF
          ELSE
            IF (NRHS == 1) THEN
               CALL ZTRSV( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     *         WCB(PPIV_COURANT), 1 )
            ELSE
              CALL ZTRSM('L','L','N','N',NPIV, NRHS, ONE,
     *                   A(APOS), LIELL, WCB(PPIV_COURANT),
     *                   NPIV)
            ENDIF
          END IF
        END IF
      END IF
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
        IF ( KEEP(50) .eq. 0 ) THEN
          APOS1 = APOS  + NPIV * LIELL
        ELSE
          APOS1 = APOS + NPIV * NPIV
        END IF
        IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
          NUPDATE = NCB
        ELSE
          NUPDATE = NELIM
        END IF
      ELSE
        APOS1 = APOS + NPIV
        NUPDATE = NCB
      END IF
      IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          IF ( NRHS == 1 ) THEN
                CALL ZGEMV('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     *            NPIV,  WCB(PPIV_COURANT), 1, ONE,
     *            WCB(PCB_COURANT), 1)
          ELSE
             CALL ZGEMM('T', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     *            A(APOS1), NPIV, WCB(PPIV_COURANT), NPIV, ONE,
     *            WCB(PCB_COURANT), NCB)
          END IF
        ELSE
          IF ( NRHS == 1 ) THEN
                CALL ZGEMV('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     *            LIELL, WCB(PPIV_COURANT), 1,
     *            ONE, WCB(PCB_COURANT), 1 )
          ELSE
             CALL ZGEMM('N', 'N', NUPDATE, NRHS, NPIV, ALPHA,
     *            A(APOS1), LIELL, WCB(PPIV_COURANT), NPIV, ONE,
     *            WCB(PCB_COURANT), NCB)
          END IF
        END IF
      END IF
      IFR = PPIV_COURANT - 1
      IF (BUILD_POSINRHSCOMP) THEN
         POSINRHSCOMP(STEP(INODE)) =  RHSCOMPFREEPOS
         RHSCOMPFREEPOS            = RHSCOMPFREEPOS + NPIV
      ENDIF
      IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
      IF ( KEEP(50) .eq. 0 ) THEN
        DO K=1,NRHS
           IFR =  PPIV_COURANT + (K-1)*NPIV
           RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1, K) =
     &            WCB(IFR:IFR+NPIV-1)
        ENDDO
      ELSE
        APOS1 = APOS
        JJ = J1
        DO 
           IF(JJ .GT. J3) EXIT
           IFR = IFR + 1
           IF(IW(JJ+LIELL) .GT. 0) THEN
              DO K=1, NRHS
                 RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &              WCB( IFR+(K-1)*NPIV ) * A( APOS1 )
              END DO
              APOS1 = APOS1 + NPIV + 1
              JJ = JJ+1
           ELSE
              APOS2 = APOS1+NPIV+1
              DO K=1, NRHS
                 POSWCB1 = IFR+(K-1)*NPIV
                 POSWCB2 = POSWCB1+1
                 RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = WCB(POSWCB1)*A(APOS1)
     *                + WCB(POSWCB2)*A(APOS1+1)
                 RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     *                 WCB(POSWCB1)*A(APOS1+1)
     *                + WCB(POSWCB2)*A(APOS2)
              END DO
              APOS1 = APOS2 + NPIV + 1
              JJ = JJ+2
              IFR = IFR+1
           ENDIF
        ENDDO
      END IF
      IF (KEEP(201).NE.0) THEN
        CALL ZMUMPS_598(INODE,PTRFAC,KEEP(28),
     $        A,LA,.TRUE.,IERR)
        IF(IERR.LT.0)THEN
           INFO(1)=IERR
           INFO(2)=0
           GOTO 260
        ENDIF
      END IF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .EQ. 0 ) THEN
          MYROOT = MYROOT - 1
          PLEFTWCB = PLEFTWCB - LIELL *NRHS
          IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
              DUMMY (1) = 1
              CALL ZMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID,
     *             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
          END IF
          GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
        IF (ZMUMPS_275(STEP(FPERE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
         IF ( NCB .ne. 0 ) THEN
          PTRICB(STEP(INODE)) = NCB + 1
          DO 190 I = 1, NUPDATE
            DO K=1, NRHS
             RHS( IW(J3 + I), K ) = RHS( IW(J3 + I), K )
     *       + WCB(PCB_COURANT + I-1 +(K-1)*NCB)
            ENDDO
  190     CONTINUE
          PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
          IF ( PTRICB(STEP(INODE)) == 1 ) THEN
            NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
            IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
              IPOOL( LEAF ) = FPERE
              LEAF = LEAF + 1
            ENDIF
          END IF
         ELSE
          PTRICB(STEP( INODE )) = -1
          NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
          IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
            IPOOL( LEAF ) = FPERE
            LEAF = LEAF + 1
          ENDIF
         ENDIF
        ELSE
 210      CONTINUE
          CALL ZMUMPS_78( NRHS, INODE, FPERE, NCB, NCB,
     *                          NUPDATE,
     *                          IW( J3 + 1 ), WCB( PCB_COURANT ),
     *                          ZMUMPS_275(STEP(FPERE),
     *                          PROCNODE_STEPS,SLAVEF),
     *                          ContVec,
     *                          COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL ZMUMPS_303( .FALSE., FLAG,
     *             BUFR, LBUFR, LBUFR_BYTES,
     *             MYID, SLAVEF, COMM,
     *             N, NRHS, IPOOL, LPOOL, III, LEAF,
     *             NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *             IWCB, LIWCB,
     *             WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *             PTRICB, INFO, KEEP,KEEP8, STEP,
     *             PROCNODE_STEPS, 
     *             RHS, LRHS)
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 210
          ELSE IF ( IERR .EQ. -2 ) THEN
             INFO( 1 ) = -17
             INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     *            ( NUPDATE + 3 ) * KEEP( 34 )
             GOTO 260
          ELSE IF ( IERR .EQ. -3 ) THEN
             INFO( 1 ) = -20
             INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     *            ( NUPDATE + 3 ) * KEEP( 34 )
             GOTO 260
          END IF
        ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     *     .and. NPIV .NE. 0 ) THEN
        DO ISLAVE = 1, NSLAVES
          PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +XSIZE)
          CALL ZMUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB - NELIM, 
     &                NSLAVES, 
     &                Effective_CB_Size, FirstIndex )
 222      CALL ZMUMPS_72( NRHS,
     *    INODE, FPERE,
     *    Effective_CB_Size, NCB, NPIV,
     *    WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     *    WCB( PPIV_COURANT ),
     *    PDEST, COMM, IERR )
          IF ( IERR .EQ. -1 ) THEN
            CALL ZMUMPS_303( .FALSE., FLAG,
     *               BUFR, LBUFR, LBUFR_BYTES,
     *               MYID, SLAVEF, COMM,
     *               N, NRHS, IPOOL, LPOOL, III, LEAF,
     *               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     *               IWCB, LIWCB,
     *               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     *               PTRICB, INFO, KEEP,KEEP8, STEP,
     *               PROCNODE_STEPS, 
     *               RHS, LRHS)
            IF ( INFO( 1 )  .LT. 0 )  GOTO 270
            GOTO 222
          ELSE IF ( IERR .EQ. -2 ) THEN
            INFO( 1 ) = -17
            INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS*KEEP(35) +
     *      ( Effective_CB_Size + 4 ) * KEEP( 34 )
            GOTO 260
          ELSE IF ( IERR .EQ. -3 ) THEN
            INFO( 1 ) = -20
            INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS*KEEP(35) +
     *      ( Effective_CB_Size + 4 ) * KEEP( 34 )
            GOTO 260
          END IF
        END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL*NRHS
  270 CONTINUE
      RETURN
  260 CONTINUE
      CALL ZMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_302
      RECURSIVE SUBROUTINE ZMUMPS_303( BLOQ, FLAG,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           MYID, SLAVEF, COMM,
     *           N, NRHS, IPOOL, LPOOL, III, LEAF,
     *           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     *           IWCB, LIWCB,
     *           WCB, LWCB, POSWCB,
     *           PLEFTWCB, POSIWCB,
     *           PTRICB, INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS,
     *           RHS, LRHS)
      IMPLICIT NONE
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER LIW, LA
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      COMPLEX*16 WCB( LWCB ), A( LA )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LRHS
      COMPLEX*16 RHS(LRHS, NRHS)
      LOGICAL FLAG
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL ZMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     *                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL ZMUMPS_323( BUFR, LBUFR, LBUFR_BYTES,
     *          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     *          N, NRHS, IPOOL, LPOOL, III, LEAF,
     *          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     *          IWCB, LIWCB,
     *          WCB, LWCB, POSWCB,
     *          PLEFTWCB, POSIWCB,
     *          PTRICB, INFO, KEEP,KEEP8, STEP,
     *          PROCNODE_STEPS, 
     *          RHS, LRHS)
         END IF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_303
      SUBROUTINE ZMUMPS_249(N, A, LA, IW, LIW, W, LWC,
     *    RHS, LRHS, NRHS, 
     *    RHSCOMP, LRHSCOMP, POSINRHSCOMP,
     *    PTRICB, PTRACB, IWCB, LIWW, W2, 
     *    NE_STEPS, NA, LNA, STEP,
     *    FRERE, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO, 
     *    PROCNODE_STEPS,
     *    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     *    KEEP,KEEP8, RHS_ROOT, MTYPE, 
     *
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE ZMUMPS_OOC
      USE ZMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER N,LA,LIW,LIWW,LWC,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL)
      INTEGER INFO(40)
      INTEGER PTRIST(KEEP(28)),PTRFAC(KEEP(28)),
     *        PTRICB(KEEP(28)),PTRACB(KEEP(28))
      INTEGER LRHS, NRHS
      COMPLEX*16 A(LA), RHS(LRHS,NRHS), W(LWC)
      COMPLEX*16 W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      COMPLEX*16 RHSCOMP(LRHSCOMP,NRHS)
      COMPLEX*16 RHS_ROOT( * )
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      LOGICAL FLAG
      INTEGER POSIWCB,POSWCB,K
      INTEGER APOS,APOSCB,NPIV
      INTEGER IPOS,IPOSCB,LIELL,NELIM,IFR,JJ,I
      INTEGER J1,J2,J,IST,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,POOL_FIRST_POS,TMP
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER NSLAVES
      COMPLEX*16 ALPHA,ONE,ZERO
      PARAMETER(ONE=1.0D0, ALPHA=-1.0D0, ZERO=0.0D0)
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER SSII,POSII, POSINDICES, IPOSINRHSCOMP
      INTEGER DUMMY(1)
      INTEGER PLEFTW, LDA, PTWCB
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      LOGICAL LTLEVEL2, IN_SUBTREE
      INCLUDE 'mumps_headers.h'
      INTEGER TMPNODE
      LOGICAL BLOCK_SEQUENCE
      LOGICAL SKIP
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      LOGICAL ZMUMPS_283, ZMUMPS_170
      INTEGER ZMUMPS_330
      EXTERNAL ZGEMV, ZTRSV, ZTRSM, ZGEMM,
     *         ZMUMPS_283, ZMUMPS_330, 
     *         ZMUMPS_170
      PLEFTW = 1
      POSIWCB = LIWW
      POSWCB = LWC
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (ZMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1
      IF (MYLEAF .EQ. -1) THEN
        MYLEAF = 0
        DO I=1, NBLEAF
          INODE=NA(I+2)
          IF (ZMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *         SLAVEF) .EQ. MYID) THEN
            MYLEAF = MYLEAF + 1
          ENDIF
        ENDDO
      ENDIF
      MYLEAFE=MYLEAF
      NBFINF = SLAVEF
      IF (MYLEAFE .EQ. 0) THEN
        CALL ZMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     *                  SLAVEF)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) THEN
          GOTO 340
        ENDIF
      ENDIF
      BLOCK_SEQUENCE = .FALSE.
   50 CONTINUE
      BLOQ = ( (  III .EQ. IIPOOL  )
     $     )
      CALL ZMUMPS_41( BLOQ, FLAG, BUFR, LBUFR,
     *     LBUFR_BYTES, MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *     IPOOL, LPOOL, STEP,  FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *     RHS, LRHS, NRHS, MTYPE, 
     *     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *     )
      IF ( INFO(1) .LT. 0 ) GOTO 340
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(IIPOOL-1)
          IIPOOL = IIPOOL - 1
          GO TO 60
        ENDIF
      END IF                    !   .not.(FLAG)
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
   60 CONTINUE
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
        IPOS = PTRIST(STEP(INODE))+XSIZE
        LIELL = IW(IPOS+3)
        NPIV  = LIELL
        IPOS =  PTRIST(STEP(INODE)) + 5 + XSIZE
        IF ( MTYPE .EQ. 1 .AND. KEEP(50) .EQ. 0) THEN
          J1   = IPOS + LIELL + 1
          J2   = IPOS + LIELL + NPIV
        ELSE
          J1   = IPOS + 1
          J2   = IPOS + NPIV
        END IF
        IFR  = 0
        DO JJ = J1, J2
          J  = IW( JJ )
          IFR = IFR + 1
          DO K=1,NRHS
            RHS(J,K) = RHS_ROOT(IFR+NPIV*(K-1))
          ENDDO
        END DO 
        IN = INODE
  270   IN = FILS(IN)
        IF (IN .GT. 0) GOTO 270
        IF (IN .EQ. 0) THEN
          MYLEAFE = MYLEAFE - 1
          IF (MYLEAFE .EQ. 0) THEN
            CALL ZMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
            NBFINF = NBFINF - 1
            IF (NBFINF .EQ. 0) GOTO 340
          ENDIF
          GOTO 50
        ENDIF
        IF   = -IN
        LONG = NPIV
        NBFILS = NE_STEPS(STEP(INODE))
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO I = 1, NBFILS
          IF (ZMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
     *       .EQ. MYID) THEN
             IPOOL(IIPOOL) = IF
             IIPOOL = IIPOOL + 1
          ELSE
            PROCDEST = ZMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *                 SLAVEF)
            IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600          CALL ZMUMPS_78( NRHS, IF, 0, 0,
     *           LONG, LONG, IW( J1 ),
     *           RHS_ROOT( 1 ), PROCDEST,
     *           NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL ZMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *          PLEFTW, KEEP,KEEP8,
     *          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHS, LRHS, NRHS, MTYPE,
     *          RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *           )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 600
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = LONG * KEEP(35) +
     *                      ( LONG + 2 ) * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = LONG * KEEP(35) +
     *                      ( LONG + 2 ) * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF ( IERR .NE. 0 ) CALL ZMUMPS_ABORT()
          ENDIF
          IF = FRERE(STEP(IF))
#if defined(Mila_Print)
#endif
        ENDDO
        IF (IIPOOL.NE.POOL_FIRST_POS) THEN
            DO I=1,(IIPOOL-POOL_FIRST_POS)/2
               TMP=IPOOL(POOL_FIRST_POS+I-1)
               IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
               IPOOL(IIPOOL-I)=TMP
            ENDDO
         ENDIF
        GOTO 50
      END IF
      IN_SUBTREE = ZMUMPS_170( 
     *          STEP (INODE), 
     *          PROCNODE_STEPS, SLAVEF ) 
      LTLEVEL2= ( 
     *   (ZMUMPS_330(STEP(INODE),PROCNODE_STEPS,
     *         SLAVEF).eq.2 ) .AND.
     *   (MTYPE.NE.1)   )
      NPIV = IW(PTRIST(STEP(INODE))+2+XSIZE+1)
      IF ((NPIV.NE.0).AND.(LTLEVEL2)) THEN
            IPOS  = PTRIST(STEP(INODE)) + 2 + XSIZE
            LIELL = IW(IPOS-2)+IW(IPOS+1)
            NELIM = IW(IPOS-1)
            IPOS  = IPOS + 1
            NPIV  = IW(IPOS)
            NCB   = LIELL - NPIV - NELIM
            IPOS  = IPOS + 2
            NSLAVES = IW( IPOS )
            Offset = 0  !!       = NPIV + NELIM
            IPOS = IPOS + NSLAVES   
            IW(PTRIST(STEP(INODE))+XXS)= C_FINI+NSLAVES
           IF ( POSIWCB - NCB - 2 .LT. 0 .or.
     *          POSWCB - NCB*NRHS .LT. PLEFTW - 1 ) THEN
             CALL ZMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
             IF ( POSWCB - NCB*NRHS .LT. PLEFTW - 1 ) THEN
               INFO( 1 ) = -11
               INFO( 2 ) = NCB * NRHS - POSWCB - PLEFTW + 1
               GOTO 330
             END IF
             IF ( POSIWCB - NCB - 2 .LT. 0 ) THEN
               INFO( 1 ) = -14
               INFO( 2 ) = NCB + 2 - POSIWCB
               GO TO 330
             END IF
           END IF
           POSIWCB = POSIWCB - NCB - 2
           POSWCB  = POSWCB - NCB*NRHS
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1
           IWCB( PTRICB(STEP( INODE ))     ) = NCB
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
              POSINDICES = IPOS + LIELL + 1
           ELSE
              POSINDICES = IPOS + 1
           END IF
           IWCB( PTRICB(STEP( INODE )) + 2 :
     *      PTRICB(STEP( INODE )) + 1 + NCB )
     *        = IW( POSINDICES + NPIV + NELIM : POSINDICES + LIELL - 1 )
           IF ( NCB.EQ.0 ) THEN
             write(6,*) ' Internal Error type 2 node with no CB '
             CALL ZMUMPS_ABORT()
           ENDIF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
               J1 = IPOS + LIELL + NPIV + NELIM +1
               J2 = IPOS + 2 * LIELL
           ELSE
               J1 = IPOS + NPIV + NELIM +1
               J2 = IPOS + LIELL
           END IF
           IFR = PTRACB(STEP( INODE )) - 1
           DO JJ = J1, J2
               J = IW(JJ)
               IFR = IFR + 1
               DO K=1, NRHS
                 W(IFR+(K-1)*NCB) = RHS(J,K)
               ENDDO
           ENDDO
           DO ISLAVE = 1, NSLAVES
              CALL ZMUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE+XSIZE)
             CALL ZMUMPS_63(NRHS, INODE,
     *             W(Offset+PTRACB(STEP(INODE))), EffectiveSize, 
     *             NCB, DEST,
     *             BACKSLV_MASTER2SLAVE,
     *             COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL ZMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS,
     *          PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHS, LRHS, NRHS, MTYPE,
     *          RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *     )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     *                            2 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = EffectiveSize * KEEP(35) +
     *                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
           END DO
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
           CALL ZMUMPS_151(NRHS,N, KEEP(28), IWCB, LIWW, W, LWC,
     *             POSWCB,POSIWCB,PTRICB,PTRACB)
           GOTO 50
      ENDIF   ! of LTLEVEL2
      IPOS = PTRIST(STEP(INODE)) + 2 + XSIZE
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      IF (KEEP(201).NE.0) THEN
       IF(.NOT.ZMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $     PTRFAC,KEEP(28),A,LA,IERR))THEN
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
         CALL ZMUMPS_578(INODE,PTRFAC,
     $        KEEP,KEEP8,A,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
         CALL ZMUMPS_577(A(PTRFAC(IW(IPOS))),INODE,IERR
     $   )
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ELSE
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
       ENDIF
      ENDIF                     ! For OOC runs
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + XSIZE )
      IPOS = IPOS + 1 + NSLAVES
      LONG = 0
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ELSE
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      END IF
      IF (IN_SUBTREE) THEN
        PTWCB = PLEFTW
        IF ( POSWCB .LT. LIELL*NRHS ) THEN
          CALL ZMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     *                 POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. LIELL*NRHS ) THEN
            INFO(1) = -11
            INFO(2) = LIELL*NRHS - POSWCB
            GOTO 330
          END IF
        END IF
      ELSE
        IF ( POSIWCB - LIELL - 2 .LT. 0 .or.
     *     POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
          CALL ZMUMPS_95( NRHS, N, KEEP(28), IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = LIELL * NRHS - POSWCB - PLEFTW + 1
            GOTO 330
          END IF
          IF ( POSIWCB - LIELL - 2 .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = LIELL + 2 - POSIWCB
            GO TO 330
          END IF
        END IF
        POSIWCB = POSIWCB - LIELL - 2
        POSWCB  = POSWCB - LIELL*NRHS
        PTRICB(STEP( INODE )) = POSIWCB + 1
        PTRACB(STEP( INODE )) = POSWCB  + 1
        IWCB( PTRICB(STEP( INODE ))     ) = LIELL
        IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
        IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
           POSINDICES = IPOS + LIELL + 1
        ELSE
           POSINDICES = IPOS + 1
        END IF
        IWCB( PTRICB(STEP( INODE )) + 2 :
     *   PTRICB(STEP( INODE )) + 1 + LIELL )
     *     = IW( POSINDICES : POSINDICES + LIELL - 1 )
        PTWCB = PTRACB(STEP( INODE )) 
      ENDIF
      IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
      DO K=1, NRHS
        DO JJ = J1, J2
          W(PTWCB+JJ-J1+(K-1)*LIELL) = RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
        ENDDO
      END DO
      IFR   = PTWCB + NPIV - 1
      IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO JJ = J1, J2
          J = IW(JJ)
          IFR = IFR + 1
          DO K=1, NRHS
            W(IFR+(K-1)*LIELL) = RHS(J,K)
          ENDDO
        ENDDO
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + NPIV
          IF (NRHS == 1) THEN
            CALL ZGEMV( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     *              W(NPIV + PTWCB), 1,
     *              ONE,
     *              W(PTWCB), 1 )
          ELSE
            CALL ZGEMM('T','N', NPIV, NRHS, NCB, ALPHA, A(IST), LIELL,
     *              W(NPIV+PTWCB), LIELL, ONE,
     *              W(PTWCB), LIELL)
          ENDIF
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + NPIV * LIELL
          ELSE
            IST = APOS + NPIV * NPIV
          END IF
            IF ( NRHS == 1 ) THEN
              CALL ZGEMV( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     *                W( NPIV + PTWCB ),
     *                1, ONE,
     *                W(PTWCB), 1 )
            ELSE
                CALL ZGEMM( 'N', 'N', NPIV, NRHS, NCB, ALPHA,
     *                A(IST), NPIV, W(NPIV+PTWCB),LIELL,
     *                ONE, W(PTWCB),LIELL)
            END IF
        END IF ! MTYPE == 1
      ENDIF  ! LIELL > NPIV
      IF ( MTYPE .eq. 1 ) THEN
        IF ( NRHS == 1 ) THEN
          CALL ZTRSV('L', 'T', 'N', NPIV, A(APOS), LIELL,
     *              W(PTWCB), 1)
        ELSE
          CALL ZTRSM('L','L','T','N', NPIV, NRHS, ONE, A(APOS),
     *              LIELL, W(PTWCB), LIELL)
        ENDIF
      ELSE
        IF ( KEEP(50) .EQ. 0 ) THEN
          IF ( NRHS == 1 ) THEN
            CALL ZTRSV('U','N','U', NPIV, A(APOS), LIELL,
     *              W(PTWCB), 1)
          ELSE
            CALL ZTRSM('L','U','N','U', NPIV, NRHS, ONE, A(APOS),
     *                 LIELL,W(PTWCB),LIELL)
          END IF
        ELSE
          IF ( NRHS == 1 ) THEN
            CALL ZTRSV('U','N','U', NPIV, A(APOS), NPIV,
     *              W(PTWCB), 1)
          ELSE
            CALL ZTRSM('L','U','N','U',NPIV, NRHS, ONE, A(APOS),
     *           NPIV, W(PTWCB), LIELL)
          END IF
        END IF
      END IF
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      DO 150 I = 1, NPIV
        JJ = IW(J1 + I - 1)
        DO K=1, NRHS
          RHS(JJ, K) = W(PTWCB+I-1+(K-1)*LIELL)
        ENDDO
  150 CONTINUE
  160 CONTINUE
      IF (KEEP(201).NE.0) THEN
         CALL ZMUMPS_598(INODE,PTRFAC,KEEP(28),
     $        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL ZMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                     FEUILLE, SLAVEF )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      IF (IN_SUBTREE) THEN
        DO I = 1, NBFILS
          IPOOL((IIPOOL-I+1)+NBFILS-I) = IF
          IIPOOL = IIPOOL + 1
          IF = FRERE(STEP(IF))
        ENDDO
      ELSE
        DEBUT = .TRUE.
        DO i = 0, SLAVEF - 1
          DEJA_SEND( i ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO 190 I = 1, NBFILS
          IF (ZMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *      SLAVEF) .EQ. MYID) THEN
            IPOOL(IIPOOL) = IF
            IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
          ELSE
            PROCDEST = ZMUMPS_275(STEP(IF),PROCNODE_STEPS,SLAVEF)
            IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400          CONTINUE
              CALL ZMUMPS_78( NRHS, IF, 0, 0, LIELL,
     *         LIELL,
     *         IWCB( PTRICB(STEP( INODE )) + 2),
     *         W   ( PTRACB(STEP( INODE ))), PROCDEST,
     *         NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL ZMUMPS_41(
     *          .FALSE., FLAG,
     *          BUFR, LBUFR, LBUFR_BYTES,
     *          MYID, SLAVEF, COMM,
     *          N, IWCB, LIWW, POSIWCB,
     *          W, LWC, POSWCB,
     *          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *          IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *          PLEFTW, KEEP,KEEP8,
     *          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *          RHS, LRHS, NRHS, MTYPE, 
     *          RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *     )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 400
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF = FRERE(STEP(IF))
          ENDIF
  190   CONTINUE
         DO I=1,(IIPOOL-POOL_FIRST_POS)/2
            TMP=IPOOL(POOL_FIRST_POS+I-1)
            IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
            IPOOL(IIPOOL-I)=TMP
         ENDDO
        IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
        CALL ZMUMPS_151(NRHS,N, KEEP(28), IWCB, LIWW, 
     *     W, LWC,
     *     POSWCB,POSIWCB,PTRICB,PTRACB)
      ENDIF
      GOTO 50
  330 CONTINUE
      CALL ZMUMPS_242(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     * SLAVEF)
  340 CONTINUE
      CALL ZMUMPS_150( MYID,COMM,BUFR,
     *                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE ZMUMPS_249
      RECURSIVE SUBROUTINE ZMUMPS_41(
     *     BLOQ, FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *     IPOOL, LPOOL, STEP, FRERE, FILS, PROCNODE_STEPS,
     *     PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, RHS,
     *     LRHS, NRHS, MTYPE,
     *     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *        )
      IMPLICIT NONE
      LOGICAL BLOQ, FLAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      COMPLEX*16 W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, LPOOL
      INTEGER IPOOL( LPOOL )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER LIW, LA
      INTEGER PTRIST(KEEP(28)), PTRFAC(KEEP(28)), IW( LIW )
      COMPLEX*16 A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      COMPLEX*16 RHS(LRHS, NRHS)
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      COMPLEX*16 RHSCOMP(LRHSCOMP,NRHS)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     *                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *                   FLAG, STATUS, IERR )
      END IF
      IF (FLAG) THEN
         MSGSOU=STATUS(MPI_SOURCE)
         MSGTAG=STATUS(MPI_TAG)
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL ZMUMPS_44( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV(BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     *                   MSGTAG, COMM, STATUS, IERR)
           CALL ZMUMPS_42( MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     *                KEEP,KEEP8,
     *                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHS, LRHS, NRHS, MTYPE, 
     *                RHSCOMP, LRHSCOMP, POSINRHSCOMP 
     *         )
         END IF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_41
      RECURSIVE SUBROUTINE ZMUMPS_42(
     *                MSGTAG, MSGSOU,
     *                BUFR, LBUFR, LBUFR_BYTES,
     *                MYID, SLAVEF, COMM,
     *                N, IWCB, LIWW, POSIWCB,
     *                W, LWC, POSWCB,
     *                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *                IPOOL, LPOOL, STEP,
     *                FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *                RHS, LRHS, NRHS, MTYPE, 
     *                RHSCOMP, LRHSCOMP, POSINRHSCOMP 
     *           )
      USE ZMUMPS_OOC
      USE ZMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER MSGTAG, MSGSOU
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER LWC
      COMPLEX*16 W( LWC )
      INTEGER POSIWCB, POSWCB
      INTEGER IIPOOL, LPOOL
      INTEGER IPOOL( LPOOL )
      INTEGER NBFINF, INFO(40)
      INTEGER PLEFTW, KEEP( 500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRICB(KEEP(28)), PTRACB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER FRERE(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LIW, LA
      INTEGER IW( LIW ), PTRIST( KEEP(28) ), PTRFAC(KEEP(28))
      COMPLEX*16 A( LA ), W2( KEEP(133) )
      INTEGER LRHS, NRHS
      COMPLEX*16  RHS(LRHS, NRHS)
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP(KEEP(28))
      COMPLEX*16 RHSCOMP(LRHSCOMP,NRHS)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER POSITION, IF, INODE, IERR, LONG, DUMMY(1)
      INTEGER P_UPDATE, P_SOL_MAS, LIELL, K
      INTEGER NPIV, NPIV_CHECK, NROW_L, APOS, IPOS, NROW_RECU
      INTEGER I, JJ, IN, PROCDEST, J1, J2, IFR, IST, LDA
      INTEGER NSLAVES, NELIM, J, POSINDICES, INODEPOS,
     &        IPOSINRHSCOMP
      LOGICAL FLAG
      LOGICAL OOCON
      COMPLEX*16 ZERO, ALPHA, ONE
      PARAMETER( ZERO = 0.0D0, ALPHA = -1.0D0, ONE = 1.0D0)
      INCLUDE 'mumps_headers.h'
       INTEGER POOL_FIRST_POS, TMP
      LOGICAL DEJA_SEND( 0:SLAVEF-1 )
      INTEGER ZMUMPS_275
      EXTERNAL ZMUMPS_275, ZTRSV, ZTRSM, ZGEMV, ZGEMM
      IF (MSGTAG .EQ. FEUILLE) THEN
          NBFINF = NBFINF - 1
      ELSE IF (MSGTAG .EQ. NOEUD) THEN
          POSITION = 0
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        INODE, 1, MPI_INTEGER,
     *        COMM, IERR)
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        LONG, 1, MPI_INTEGER,
     *        COMM, IERR)
          IF (   POSIWCB - LONG - 2 .LT. 0
     *      .OR. POSWCB - PLEFTW + 1 .LT. LONG ) THEN
            CALL ZMUMPS_95(NRHS, N, KEEP(28), IWCB,
     *      LIWW, W, LWC,
     *      POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ((POSIWCB - LONG - 2 ) .LT. 0) THEN
              INFO(1)=-14
              INFO(2)=-POSIWCB + LONG + 2
              GOTO 260
            END IF
            IF ( POSWCB - PLEFTW + 1 .LT. LONG ) THEN
              INFO(1) = -11
              INFO(2) = LONG + PLEFTW - POSWCB - 1
              GOTO 260
            END IF
          ENDIF
          POSIWCB = POSIWCB - LONG
          POSWCB = POSWCB - LONG
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          IWCB(POSIWCB + 1), 
     *          LONG, MPI_INTEGER, COMM, IERR)
            DO K=1,NRHS
             CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          W(POSWCB + 1), LONG, 
     *          MPI_DOUBLE_COMPLEX, COMM, IERR)
             DO JJ=0, LONG-1
               RHS(IWCB(POSIWCB+1+JJ),K) = W(POSWCB+1+JJ)
             ENDDO
            ENDDO
            POSIWCB = POSIWCB + LONG
            POSWCB = POSWCB + LONG
          ENDIF
          POOL_FIRST_POS = IIPOOL
          IPOOL( IIPOOL ) = INODE
          IIPOOL = IIPOOL + 1
          IF = FRERE( STEP(INODE) )
          DO WHILE ( IF .GT. 0 )
            IF ( ZMUMPS_275(STEP(IF),PROCNODE_STEPS,
     *           SLAVEF) .eq. MYID ) THEN
              IPOOL( IIPOOL ) = IF
              IIPOOL = IIPOOL + 1
            END IF
            IF = FRERE( STEP( IF ) )
          END DO
        DO I=1,(IIPOOL-POOL_FIRST_POS)/2
           TMP=IPOOL(POOL_FIRST_POS+I-1)
           IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
           IPOOL(IIPOOL-I)=TMP
        ENDDO
      ELSE IF ( MSGTAG .EQ. BACKSLV_MASTER2SLAVE ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NROW_RECU, 1, MPI_INTEGER, COMM, IERR )
        IPOS   = PTRIST( STEP(INODE) ) + XSIZE
        NPIV   = - IW( IPOS     )
        NROW_L =   IW( IPOS + 1 )
        IF (KEEP(201).NE.0) THEN
         IF(.NOT.ZMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $       PTRFAC,KEEP(28),A,LA,IERR))THEN
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL ZMUMPS_578(INODE,PTRFAC,
     $          KEEP,KEEP8,A,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           CALL ZMUMPS_577(A(PTRFAC(IW( IPOS + 3 ))),INODE,IERR
     $   )
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ELSE
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
         ENDIF
        ENDIF                     ! For OOC run
        APOS   =   PTRFAC(IW( IPOS + 3 ))
        IF ( NROW_L .NE. NROW_RECU ) THEN
          WRITE(*,*) 'Error1 in 41S : NROW L/RECU=',NROW_L, NROW_RECU
          CALL ZMUMPS_ABORT()
        END IF
        LONG = NROW_L + NPIV
        IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
           CALL ZMUMPS_95(NRHS, N, KEEP(28), IWCB,
     *          LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
           IF ( POSWCB - LONG*NRHS .LT. PLEFTW - 1 ) THEN
             INFO(1) = -11
             INFO(2) = LONG * NRHS- POSWCB
             GOTO 260
           END IF
        END IF
        P_UPDATE  = PLEFTW
        P_SOL_MAS = PLEFTW + NPIV * NRHS
        PLEFTW    = P_SOL_MAS + NROW_L * NRHS
        DO K=1, NRHS
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W( P_SOL_MAS+(K-1)*NROW_L),NROW_L,
     *                   MPI_DOUBLE_COMPLEX,
     *                   COMM, IERR )
        ENDDO
        IF ( NRHS == 1 ) THEN
          CALL ZGEMV( 'N', NPIV, NROW_L, ALPHA, A( APOS ), NPIV,
     *              W( P_SOL_MAS ), 1, ZERO,
     *              W( P_UPDATE ), 1 )
        ELSE
          CALL ZGEMM( 'N', 'N', NPIV, NRHS, NROW_L, ALPHA, A(APOS),
     *                NPIV, W( P_SOL_MAS ), NROW_L, ZERO, W( P_UPDATE ),
     *                NPIV )
        END IF
        IF (KEEP(201).NE.0) THEN
         CALL ZMUMPS_598(INODE,PTRFAC,KEEP(28),
     $          A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
        ENDIF
        PLEFTW = PLEFTW - NROW_L * NRHS
 100    CONTINUE
        CALL ZMUMPS_63( NRHS, INODE, W(P_UPDATE),
     *                               NPIV, NPIV,
     *                                MSGSOU, 
     *                                BACKSLV_UPDATERHS,
     *                                COMM, IERR )
        IF ( IERR .EQ. -1 ) THEN
          CALL ZMUMPS_41(
     *     .FALSE., FLAG,
     *     BUFR, LBUFR, LBUFR_BYTES,
     *     MYID, SLAVEF, COMM,
     *     N, IWCB, LIWW, POSIWCB,
     *     W, LWC, POSWCB,
     *     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *     IPOOL, LPOOL, STEP,
     *     FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE,
     *     RHS, LRHS, NRHS, MTYPE,
     *     RHSCOMP, LRHSCOMP, POSINRHSCOMP
     * )
          IF ( INFO( 1 ) .LT. 0 ) GOTO 270
          GOTO 100
        ELSE IF ( IERR .EQ. -2 ) THEN
          INFO( 1 ) = -17
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        ELSE IF ( IERR .EQ. -3 ) THEN
          INFO( 1 ) = -20
          INFO( 2 ) = NPIV * KEEP(35) + 2 * KEEP(34)
          GOTO 260
        END IF
        PLEFTW = PLEFTW - NPIV * NRHS
      ELSE IF ( MSGTAG .EQ. BACKSLV_UPDATERHS ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   INODE, 1, MPI_INTEGER, COMM, IERR )
        IPOS  = PTRIST(STEP(INODE)) + 2 + XSIZE
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NPIV, 1, MPI_INTEGER, COMM, IERR )
          NELIM = IW(IPOS-1)
          IPOS = IPOS + 1
          NPIV = IW(IPOS)
          IPOS = IPOS + 1
          NSLAVES = IW( IPOS + 1 )
          IPOS = IPOS + 1 + NSLAVES
          OOCON  = (KEEP(201).NE.0)  ! OOC case
          IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
             J1 = IPOS + 1
             J2 = IPOS + NPIV
          ELSE
             J1 = IPOS + LIELL + 1
             J2 = IPOS + NPIV + LIELL
          END IF
        DO K=1, NRHS
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   W2, NPIV, MPI_DOUBLE_COMPLEX,
     *                   COMM, IERR )
         IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
         I = 1
         DO JJ = J1,J2   ! NPIV terms
            RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = 
     &      RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) + W2(I)
            I = I+1
         ENDDO
        ENDDO  ! Loop on NRHS
        IW(PTRIST(STEP(INODE))+XXS) = 
     *      IW(PTRIST(STEP(INODE))+XXS) - 1
        IF ( IW(PTRIST(STEP(INODE))+XXS).EQ.C_FINI ) THEN
          IF (OOCON) THEN
           INODEPOS = PTRIST(STEP(INODE)) + XSIZE + 4
           IF(.NOT.ZMUMPS_SOLVE_IS_INODE_IN_MEM(INODE,
     $         PTRFAC,KEEP(28),A,LA,IERR))THEN
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL ZMUMPS_578(INODE,PTRFAC,
     $            KEEP,KEEP8,A,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             CALL ZMUMPS_577(A(PTRFAC(IW( INODEPOS ))),INODE,IERR
     $       )
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ELSE
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
           ENDIF
           IF ( POSIWCB - LIELL - 2 .LT. 0 .or.
     *         POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
            CALL ZMUMPS_95( NRHS, N, KEEP(28), IWCB, 
     *          LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ( POSWCB - LIELL*NRHS .LT. PLEFTW - 1 ) THEN
              INFO( 1 ) = -11
              INFO( 2 ) = LIELL * NRHS - POSWCB - PLEFTW + 1
              GOTO 260
            END IF
            IF ( POSIWCB - LIELL - 2 .LT. 0 ) THEN
              INFO( 1 ) = -14
              INFO( 2 ) = LIELL + 2 - POSIWCB
              GO TO 260
            END IF
           END IF
          ENDIF    ! if OOCON
           POSIWCB = POSIWCB - LIELL - 2
           POSWCB  = POSWCB - LIELL*NRHS
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1
           IWCB( PTRICB(STEP( INODE ))     ) = LIELL
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IPOS = PTRIST(STEP(INODE)) + XSIZE + 5 + NSLAVES
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
             POSINDICES = IPOS + LIELL + 1
           ELSE
             POSINDICES = IPOS + 1
           END IF
           IWCB( PTRICB(STEP( INODE )) + 2 :
     *      PTRICB(STEP( INODE )) + 1 + LIELL )
     *       = IW( POSINDICES : POSINDICES + LIELL - 1 )
           IPOSINRHSCOMP =  POSINRHSCOMP(STEP(INODE))
           IFR = PTRACB(STEP( INODE ))
           DO K=1, NRHS
             DO JJ = J1, J2
               W(IFR+JJ-J1+(K-1)*LIELL) = 
     &           RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
             ENDDO
           END DO
           IFR = PTRACB(STEP(INODE))-1+NPIV
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
             J1 = IPOS + LIELL + NPIV + 1
             J2 = IPOS + 2 * LIELL
           ELSE
             J1 = IPOS + NPIV + 1
             J2 = IPOS + LIELL
           END IF
           DO JJ = J1, J2   ! LIELL - NPIV terms
              J = IW(JJ)
              IFR = IFR + 1
              DO K=1, NRHS
                W(IFR+(K-1)*LIELL) = RHS(J,K)
              ENDDO
           ENDDO
           IF (NELIM .GT.0) THEN
            INODEPOS = PTRIST(STEP(INODE)) + 4 + XSIZE
            APOS = PTRFAC(IW(INODEPOS))
            IF ( KEEP(50) .eq. 0 ) THEN
                IST = APOS + NPIV * LIELL
            ELSE
                IST = APOS + NPIV * NPIV
            END IF
            IF ( NRHS == 1 ) THEN
                CALL ZGEMV( 'N', NPIV, NELIM, ALPHA,
     *                A( IST ), NPIV,
     *                W( NPIV + PTRACB(STEP(INODE)) ),
     *                1, ONE,
     *                W(PTRACB(STEP(INODE))), 1 )
             ELSE
                CALL ZGEMM( 'N', 'N', NPIV, NRHS, NELIM, ALPHA,
     *                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     *                ONE, W(PTRACB(STEP(INODE))),LIELL)
             END IF
           ENDIF ! update due to nelim var
          IPOS  = PTRIST(STEP(INODE)) + 2 + XSIZE
          LIELL = IW(IPOS-2)+IW(IPOS+1)
          NPIV_CHECK = IW( IPOS + 1 )
          IF ( NPIV .NE. NPIV_CHECK ) THEN
            WRITE(*,*) 'Error 2 in usedbympima41s.F'
            CALL ZMUMPS_ABORT()
          END IF
          APOS  = PTRFAC(IW( IPOS +2 ))
          NSLAVES = IW (IPOS+3)
          IF ( KEEP(50) .eq. 0 ) THEN
           LDA = LIELL
          ELSE
           LDA = NPIV
          ENDIF
          IF ( NRHS == 1 ) THEN
              CALL ZTRSV( 'U', 'N', 'U', NPIV, A(APOS), LDA,
     *                  W(PTRACB(STEP(INODE))),1)
          ELSE
             CALL ZTRSM( 'L','U', 'N', 'U', NPIV, NRHS, ONE,
     *                   A(APOS), LDA,
     *                   W(PTRACB(STEP(INODE))),LIELL)
          END IF
          IF (OOCON) THEN
           CALL ZMUMPS_598(INODE,PTRFAC,KEEP(28),
     $          A,LA,.TRUE.,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
          ENDIF
          IPOS =   PTRIST(STEP(INODE)) +  XSIZE + 6 + NSLAVES   ! pointer to first index in the list
          DO I = 1, NPIV
            JJ = IW( IPOS + I - 1 )
            DO K=1,NRHS
              RHS( JJ, K ) = W( PTRACB(STEP(INODE))+I-1
     *         + (K-1)*LIELL )
            ENDDO
          END DO
          IN = INODE
  200     IN = FILS(IN)
          IF (IN .GT. 0) GOTO 200
          IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
              CALL ZMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     *                       FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
            ENDIF
            IWCB( PTRICB(STEP(INODE)) + 1 ) = 0
            CALL ZMUMPS_151(NRHS, N, KEEP(28),
     *          IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
            GOTO 270
          ENDIF  ! IN is a leave node
          DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
          END DO
          IN = -IN
 300      CONTINUE
           POOL_FIRST_POS  = IIPOOL
            IF (ZMUMPS_275(STEP(IN),PROCNODE_STEPS,
     *          SLAVEF) .EQ. MYID) THEN
              IPOOL(IIPOOL ) = IN
              IIPOOL = IIPOOL + 1
            ELSE
              PROCDEST = ZMUMPS_275( STEP(IN), PROCNODE_STEPS,
     *                   SLAVEF )
              IF ( .NOT. DEJA_SEND( PROCDEST ) ) THEN
 110            CALL ZMUMPS_78( NRHS, IN, 0, 0,
     *          LIELL, LIELL,
     *          IWCB( PTRICB(STEP(INODE))+2 ),
     *          W( PTRACB(STEP(INODE))),
     *          PROCDEST, NOEUD, COMM, IERR )
                IF ( IERR .EQ. -1 ) THEN
                  CALL ZMUMPS_41(
     *            .FALSE., FLAG,
     *            BUFR, LBUFR, LBUFR_BYTES,
     *            MYID, SLAVEF, COMM,
     *            N, IWCB, LIWW, POSIWCB,
     *            W, LWC, POSWCB,
     *            IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     *            IPOOL, LPOOL, STEP,
     *            FRERE, FILS, PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     *            PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     *            RHS, LRHS, NRHS, MTYPE, 
     *            RHSCOMP, LRHSCOMP, POSINRHSCOMP
     *            )
                  IF ( INFO( 1 ) .LT. 0 ) GOTO 270
                  GOTO 110
                ELSE IF ( IERR .eq. -2 ) THEN
                  INFO(1) = -17
                  INFO(2) = LIELL * NRHS * KEEP(35) +
     *                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                ELSE IF ( IERR .eq. -3 ) THEN
                  INFO(1) = -20
                  INFO(2) = LIELL * NRHS * KEEP(35) +
     *                    ( LIELL + 2 ) * KEEP(34)
                  GOTO 260
                END IF
                DEJA_SEND( PROCDEST ) = .TRUE.
              END IF
            END IF
          IN = FRERE( STEP( IN ) )
          IF ( IN .GT. 0 ) GOTO 300
          DO I=1,(IIPOOL-POOL_FIRST_POS)/2
           TMP=IPOOL(POOL_FIRST_POS+I-1)
           IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
           IPOOL(IIPOOL-I)=TMP
          ENDDO
          IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
          CALL ZMUMPS_151(NRHS, N, KEEP(28),
     *          IWCB, LIWW, W, LWC,
     *          POSWCB, POSIWCB, PTRICB, PTRACB)
        END IF
      ELSE IF (MSGTAG.EQ.TERREUR) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GO TO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     *      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
      ELSE
          INFO(1) = -100
          INFO(2) = MSGTAG
          GOTO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL ZMUMPS_44( MYID, SLAVEF, COMM )
 270  CONTINUE
      RETURN
      END SUBROUTINE ZMUMPS_42
      SUBROUTINE ZMUMPS_286( NRHS, DESCA_PAR, DESCB_PAR,
     *  CNTXT_PAR,LOCAL_M,LOCAL_N,MBLOCK,NBLOCK,
     *  IPIV,LPIV,MASTER_ROOT,MYID,COMM,
     *  RHS_SEQ,SIZE_ROOT,A,INFO,MTYPE,LDLT )
      IMPLICIT NONE
      INTEGER NRHS, MTYPE
      INTEGER DESCA_PAR( 9 ), DESCB_PAR( 9 )
      INTEGER LOCAL_M, LOCAL_N, MBLOCK, NBLOCK
      INTEGER CNTXT_PAR, MASTER_ROOT, SIZE_ROOT
      INTEGER MYID, COMM
      INTEGER LPIV, IPIV( LPIV )
      INTEGER INFO(40), LDLT
      COMPLEX*16 RHS_SEQ( SIZE_ROOT *NRHS)
      COMPLEX*16 A( LOCAL_M, LOCAL_N )
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER IERR, NPROW, NPCOL, MYROW, MYCOL
      INTEGER LOCAL_N_RHS
      COMPLEX*16, ALLOCATABLE, DIMENSION( :,: ) ::RHS_PAR
      EXTERNAL NUMROC
      INTEGER  NUMROC
      INTEGER allocok
      CALL BLACS_GRIDINFO( CNTXT_PAR, NPROW, NPCOL, MYROW, MYCOL )
      LOCAL_N_RHS = NUMROC(NRHS, NBLOCK, MYCOL, 0, NPCOL)
      LOCAL_N_RHS = MAX(1,LOCAL_N_RHS)
      ALLOCATE(RHS_PAR(LOCAL_M, LOCAL_N_RHS),stat=allocok)
      IF (allocok > 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root.'
        WRITE(*,*) ' Reduce number of right hand sides.'
        CALL ZMUMPS_ABORT()
      ENDIF
      CALL ZMUMPS_290( MYID, SIZE_ROOT, NRHS, RHS_SEQ,
     *      LOCAL_M, LOCAL_N_RHS,
     *      MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *      NPROW, NPCOL, COMM )
      IF ( LDLT .eq. 0 .OR. LDLT .eq. 2 ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          CALL PZGETRS('N',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR,1,1,DESCB_PAR,IERR)
        ELSE
          CALL PZGETRS('T',SIZE_ROOT,NRHS,A,1,1,DESCA_PAR,IPIV,
     *      RHS_PAR, 1, 1, DESCB_PAR,IERR)
        END IF
      ELSE
        CALL PZPOTRS( 'L', SIZE_ROOT, NRHS, A, 1, 1, DESCA_PAR,
     *    RHS_PAR, 1, 1, DESCB_PAR, IERR )
      END IF
      IF ( IERR .LT. 0 ) THEN
        WRITE(*,*) ' Problem during solve of the root'
        CALL ZMUMPS_ABORT()
      END IF
      CALL ZMUMPS_156( MYID, SIZE_ROOT, NRHS,
     *    RHS_SEQ, LOCAL_M, LOCAL_N_RHS,
     *    MBLOCK, NBLOCK, RHS_PAR, MASTER_ROOT,
     *    NPROW, NPCOL, COMM )
      DEALLOCATE(RHS_PAR)
      RETURN
      END SUBROUTINE ZMUMPS_286
