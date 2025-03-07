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
      MODULE SMUMPS_LOAD
      implicit none
!      PRIVATE
      PUBLIC :: SMUMPS_188, SMUMPS_185,
     *       SMUMPS_189,
     *       SMUMPS_184, SMUMPS_190,
     *       SMUMPS_183, SMUMPS_187, 
     *       SMUMPS_186,
     *       SMUMPS_409,
     *       SMUMPS_384, SMUMPS_461,
     *       SMUMPS_467, SMUMPS_471,
     *       SMUMPS_472,
     *       SMUMPS_500,
     *       SMUMPS_501,
     *       SMUMPS_520,
     *       SMUMPS_513,
     *       SMUMPS_514, SMUMPS_512
     *       ,SMUMPS_533,
     $       CLEAN_POOL_MEM_INFO,COMPUTE_MAX_MEM,
     $       CHECK_MEM_CONST_FOR_POOL,SMUMPS_554,
     $       SMUMPS_553,
     $       SMUMPS_555
#if  defined (DELTA_MD_STK)
     *       ,SMUMPS_LOAD_UPDATE_MD_INFO
#endif
      DOUBLE PRECISION, DIMENSION(:),
     *       ALLOCATABLE, SAVE, PRIVATE :: LOAD_FLOPS
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: BUF_LOAD_RECV
      INTEGER, SAVE, PRIVATE :: LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES
      INTEGER, SAVE, PRIVATE :: K50, K69, K35,MAX_SURF_MASTER
      LOGICAL, SAVE, PRIVATE :: BDC_MEM, BDC_POOL, BDC_SBTR, 
     $     BDC_POOL_MNG,
     $     BDC_M2_MEM,BDC_M2_FLOPS,BDC_MD,REMOVE_NODE_FLAG,
     $     REMOVE_NODE_FLAG_MEM
      DOUBLE PRECISION, SAVE, PRIVATE :: REMOVE_NODE_COST,
     $     REMOVE_NODE_COST_MEM
      INTEGER, SAVE, PRIVATE :: SBTR_WHICH_M
      DOUBLE PRECISION, DIMENSION(:),
     *       ALLOCATABLE, TARGET, SAVE, PRIVATE :: WLOAD
#if defined(OLD_LOAD_MECHANISM)
#if defined (CHECK_COHERENCE) 
      INTEGER, SAVE, PRIVATE :: NB_LEVEL2
      LOGICAL, PRIVATE :: AMI_CHOSEN,IS_DISPLAYED
#endif
#endif
#if ! defined(OLD_LOAD_MECHANISM)
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: FUTURE_NIV2
      DOUBLE PRECISION, SAVE, PRIVATE :: DELTA_LOAD, DELTA_MEM
#else
      DOUBLE PRECISION, SAVE, PRIVATE :: LAST_LOAD_SENT, 
     $           DM_LAST_MEM_SENT
#endif
      INTEGER, SAVE, PRIVATE :: CHECK_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE, TARGET, PRIVATE :: 
     $          IDWLOAD
      DOUBLE PRECISION, SAVE, PRIVATE :: COST_SUBTREE
      DOUBLE PRECISION, SAVE, PRIVATE :: ALPHA
      DOUBLE PRECISION, SAVE, PRIVATE :: BETA
      INTEGER, SAVE, PRIVATE :: MYID, NPROCS, COMM_LD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE, 
     $           PRIVATE :: POOL_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, PRIVATE, 
     $           SAVE :: SBTR_MEM
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, 
     $           PRIVATE, SAVE :: SBTR_CUR
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, 
     $           PRIVATE, SAVE :: NB_SON
      DOUBLE PRECISION, 
     $           PRIVATE, SAVE :: SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     $           PRIVATE, SAVE :: PEAK_SBTR_CUR_LOCAL
      DOUBLE PRECISION, 
     $           PRIVATE, SAVE :: MAX_PEAK_STK
      DOUBLE PRECISION, SAVE, 
     $           PRIVATE :: POOL_LAST_COST_SENT
      DOUBLE PRECISION, SAVE, 
     $           PRIVATE :: MIN_DIFF
      INTEGER, SAVE :: POS_ID,POS_MEM
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: CB_COST_ID
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE
     $           :: CB_COST_MEM
      PUBLIC :: CB_COST_ID, CB_COST_MEM,POS_MEM,POS_ID
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE,
     $        PRIVATE::MD_MEM,LU_USAGE,
     $                       TAB_MAXS
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::MEM_SUBTREE
      INTEGER  :: NB_SUBTREES,NIV1_FLAG
      INTEGER, PRIVATE  :: INDICE_SBTR,INDICE_SBTR_ARRAY
      INTEGER,SAVE :: INSIDE_SUBTREE
      PUBLIC :: NB_SUBTREES,MEM_SUBTREE,INSIDE_SUBTREE,NIV1_FLAG
      DOUBLE PRECISION, SAVE, PRIVATE :: DM_SUMLU,
     *                   DM_THRES_MEM
      DOUBLE PRECISION, DIMENSION(:),
     *   ALLOCATABLE, SAVE , PRIVATE:: DM_MEM
      INTEGER, SAVE, PRIVATE :: POOL_SIZE,ID_MAX_M2
      DOUBLE PRECISION, SAVE, PRIVATE :: MAX_M2,TMP_M2
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, PRIVATE:: POOL_NIV2
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     $      PRIVATE :: POOL_NIV2_COST, NIV2
      DOUBLE PRECISION, SAVE, PRIVATE  ::      CHK_LD
      INTEGER, DIMENSION(:),POINTER, SAVE, PRIVATE  :: 
     $         KEEP_LOAD, PROCNODE_LOAD, STEP_TO_NIV2_LOAD
      INTEGER*8, DIMENSION(:), POINTER, SAVE, PRIVATE:: KEEP8_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE :: 
     *         FILS_LOAD, STEP_LOAD,
     *         FRERE_LOAD, ND_LOAD,
     $         NE_LOAD,DAD_LOAD
      INTEGER, DIMENSION(:,:),POINTER, SAVE, PRIVATE :: CAND_LOAD
      INTEGER, DIMENSION(:),POINTER, SAVE, 
     $         PRIVATE :: MY_FIRST_LEAF,MY_NB_LEAF, MY_ROOT_SBTR
      INTEGER, DIMENSION(:),ALLOCATABLE,SAVE, 
     $         PRIVATE ::SBTR_FIRST_POS_IN_POOL
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE,SAVE, 
     $         PRIVATE ::SBTR_PEAK_ARRAY,
     $     SBTR_CUR_ARRAY
      INTEGER, SAVE, 
     $         PRIVATE :: N_LOAD
      DOUBLE PRECISION,DIMENSION(:),POINTER, SAVE :: COST_TRAV
      INTEGER, DIMENSION(:),POINTER, SAVE :: DEPTH_FIRST_LOAD
      PUBLIC :: DEPTH_FIRST_LOAD,COST_TRAV, FILS_LOAD,STEP_LOAD,
     $     FRERE_LOAD, ND_LOAD,NE_LOAD,DAD_LOAD
      INTEGER, SAVE     :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     $     SECOND_CURRENT_BEST
      PUBLIC :: ROOT_CURRENT_SUBTREE,CURRENT_BEST,
     $     SECOND_CURRENT_BEST
      CONTAINS
      SUBROUTINE SMUMPS_188( COST_SUBTREE_ARG, K64, K66,
     $     MAXS )
      IMPLICIT NONE
      DOUBLE PRECISION COST_SUBTREE_ARG
      INTEGER K64, K66,MAXS
      DOUBLE PRECISION T64, T66
      T64 = MAX ( DBLE(K64), DBLE(1) )
      T64 = MIN ( T64, DBLE(1000)  )
      T66 = MAX (DBLE(K66), DBLE(100))
      MIN_DIFF     =  ( T64 / DBLE(1000)  )* 
     *                  T66 * DBLE(1000000)
      DM_THRES_MEM = (MAXS/1000)
      COST_SUBTREE = COST_SUBTREE_ARG
      END SUBROUTINE SMUMPS_188
      SUBROUTINE SMUMPS_472(
     *  NCBSON_MAX, SLAVEF,
     *  KEEP,KEEP8,ICNTL,
     *  CAND_OF_NODE,
     *  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,INODE)
       IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER*8 KEEP8(150)
      INTEGER, intent(in) :: ICNTL(40)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER I
      INTEGER LP,MP
      LP=ICNTL(4)
      MP=ICNTL(2)
      IF ( KEEP(48) == 0 .OR. KEEP(48) .EQ. 3 ) THEN
         CALL SMUMPS_499(
     *        SLAVEF,
     *        KEEP,KEEP8,
     *        CAND_OF_NODE,
     *        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      ELSE IF ( KEEP(48) == 4 ) THEN
         CALL SMUMPS_504(
     *        SLAVEF,
     *        KEEP,KEEP8,
     *        CAND_OF_NODE,
     *        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID)
         DO I=1,NSLAVES_NODE
            IF(TAB_POS(I+1)-TAB_POS(I).LE.0)THEN
               WRITE(*,*)'probleme de partition dans 
     *SMUMPS_545'
               CALL SMUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE IF ( KEEP(48) == 5 ) THEN
         CALL SMUMPS_518(
     *        NCBSON_MAX,
     *        SLAVEF,
     *        KEEP,KEEP8,
     *        CAND_OF_NODE,
     *        MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *        TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     *        MP,LP)
         DO I=1,NSLAVES_NODE
            IF(TAB_POS(I+1)-TAB_POS(I).LE.0)THEN
               WRITE(*,*)'problem with partition in
     *SMUMPS_518'
               CALL SMUMPS_ABORT()
            ENDIF
         ENDDO
      ELSE
        WRITE(*,*) "Strategy 6 not implemented"
        CALL SMUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_472
      SUBROUTINE SMUMPS_499(
     *  SLAVEF,
     *  KEEP,KEEP8,
     *  CAND_OF_NODE,
     *  MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     *  TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER*8 KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB
      INTEGER, intent(in) :: CAND_OF_NODE(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER ITEMP, NMB_OF_CAND, NSLAVES_LESS
      DOUBLE PRECISION MSG_SIZE
      LOGICAL force_cand
      INTEGER  SMUMPS_12
      EXTERNAL SMUMPS_12
      IF ( KEEP(48) == 0 .AND. KEEP(50) .NE. 0) THEN
      write(*,*) "Internal error 2 in SMUMPS_499."
      CALL SMUMPS_ABORT()
      END IF
      IF ( KEEP(48) == 3 .AND. KEEP(50) .EQ. 0) THEN
      write(*,*) "Internal error 3 in SMUMPS_499."
      CALL SMUMPS_ABORT()
      END IF
      MSG_SIZE = DBLE( NFRONT - NCB ) * DBLE(NCB)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
      END IF
      IF (force_cand) THEN
        ITEMP=SMUMPS_409
     *       (MEM_DISTRIB(0:SLAVEF-1),
     *        CAND_OF_NODE,
     *
     *        KEEP(69), SLAVEF, MSG_SIZE,
     *        NMB_OF_CAND )
      ELSE
        ITEMP=SMUMPS_186(KEEP(69),MEM_DISTRIB,MSG_SIZE)
        NMB_OF_CAND = SLAVEF - 1
      END IF
      NSLAVES_LESS = MAX(ITEMP,1)
      NSLAVES_NODE = SMUMPS_12(KEEP(10), KEEP(48),
     *          KEEP(50),SLAVEF,
     *          NCB, NFRONT, NSLAVES_less, NMB_OF_CAND)
      CALL SMUMPS_441(
     *            KEEP,KEEP8, SLAVEF,
     *            TAB_POS,
     *            NSLAVES_NODE, NFRONT, NCB
     *             )
      IF (force_cand) THEN
        CALL SMUMPS_384(MEM_DISTRIB(0),
     *       CAND_OF_NODE, SLAVEF, NSLAVES_NODE,
     *       SLAVES_LIST)
      ELSE
        CALL SMUMPS_189(MEM_DISTRIB(0),
     *       MSG_SIZE, SLAVES_LIST, NSLAVES_NODE)
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_499
      SUBROUTINE SMUMPS_185( id, MEMORY_MD_ARG )
      USE SMUMPS_BUFFER
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(SMUMPS_STRUC), TARGET :: id
      INTEGER, intent(in) :: MEMORY_MD_ARG
      INTEGER K34_LOC,K35_LOC
      INTEGER allocok, IERR, I, BUF_LOAD_SIZE
      INTEGER MAX_SBTR
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION MEMORY_SENT
      INTEGER,DIMENSION(:),POINTER:: KEEP
      PARAMETER( ZERO=0.0d0 )
      INTEGER MEMORY_MD,la,WHAT
      STEP_TO_NIV2_LOAD=>id%ISTEP_TO_INIV2
      CAND_LOAD=>id%CANDIDATES
      ND_LOAD=>id%ND_STEPS
      KEEP_LOAD=>id%KEEP
      KEEP     =>id%KEEP
      KEEP8_LOAD=>id%KEEP8
      FILS_LOAD=>id%FILS
      FRERE_LOAD=>id%FRERE_STEPS
      DAD_LOAD=>id%DAD_STEPS
      PROCNODE_LOAD=>id%PROCNODE_STEPS
      STEP_LOAD=>id%STEP
      NE_LOAD=>id%NE_STEPS
      N_LOAD=id%N
      ROOT_CURRENT_SUBTREE=-9999
      MEMORY_MD=MEMORY_MD_ARG
      LA=id%MAXS
      MAX_SURF_MASTER=id%MAX_SURF_MASTER+
     & INT(DBLE(id%KEEP(12))*DBLE(id%MAX_SURF_MASTER)/100.0d0)
      COMM_LD = id%COMM_LOAD
      K69  = KEEP(69)
      IF ( KEEP(47) .le. 0 .OR. KEEP(47) .gt. 4 ) THEN
        write(*,*) "Internal error 1 in SMUMPS_185"
        CALL SMUMPS_ABORT()
      END IF
      CHK_LD=DBLE(0)
      BDC_MEM      = ( KEEP(47) >= 2 )
      BDC_POOL     = ( KEEP(47) >= 3 )
      BDC_SBTR     = ( KEEP(47) >= 4 )
      BDC_M2_MEM   = ( ( KEEP(80) == 2 .OR. KEEP(80) == 3 )
     *             .AND. KEEP(47) == 4 )
      BDC_M2_FLOPS   = ( KEEP(80) == 1 
     *             .AND. KEEP(47) .GE. 1 )
      BDC_MD       = (KEEP(86)==1)
      SBTR_WHICH_M       = KEEP(90)
      REMOVE_NODE_FLAG=.FALSE.
      REMOVE_NODE_FLAG_MEM=.FALSE.
      REMOVE_NODE_COST_MEM=DBLE(0)
      REMOVE_NODE_COST=DBLE(0)
      IF (KEEP(80) .LT. 0 .OR. KEEP(80)>3) THEN
        WRITE(*,*) "Unimplemented KEEP(80) Strategy"
        CALL SMUMPS_ABORT()
      ENDIF
      IF ((KEEP(80) == 2 .OR. KEEP(80)==3 ).AND. KEEP(47).NE.4)
     *  THEN
        WRITE(*,*) "Internal error 3 in SMUMPS_185"
         CALL SMUMPS_ABORT()
      END IF
      IF (KEEP(81) == 1 .AND. KEEP(47) < 2) THEN
        WRITE(*,*) "Internal error 2 in SMUMPS_185"
        CALL SMUMPS_ABORT()
      ENDIF
      BDC_POOL_MNG = ((KEEP(81) == 1).AND.(KEEP(47) >= 2))
      IF(KEEP(76).EQ.4)THEN
         DEPTH_FIRST_LOAD=>id%DEPTH_FIRST
      ENDIF
      IF(KEEP(76).EQ.5)THEN
         COST_TRAV=>id%COST_TRAV
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
         ALLOCATE(NIV2(id%NSLAVES), NB_SON(KEEP(28)),
     $            POOL_NIV2(100),POOL_NIV2_COST(100),
     $            stat=allocok)
         NB_SON=id%NE_STEPS
         NIV2=0
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES + KEEP(28) + 200
           RETURN
         ENDIF
      ENDIF
      K50      = id%SYM
      CALL MPI_COMM_RANK( COMM_LD, MYID, IERR )
      NPROCS = id%NSLAVES
      DM_SUMLU=ZERO
      POOL_SIZE=0
      IF(BDC_MD)THEN
         IF ( ALLOCATED(MD_MEM) ) DEALLOCATE(MD_MEM)
         ALLOCATE( MD_MEM( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in SMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         IF ( ALLOCATED(TAB_MAXS) ) DEALLOCATE(TAB_MAXS)
         ALLOCATE( TAB_MAXS( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in SMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         TAB_MAXS=0
         IF ( ALLOCATED(LU_USAGE) ) DEALLOCATE(LU_USAGE)
         ALLOCATE( LU_USAGE( 0: NPROCS - 1 ), stat=allocok )
         IF ( allocok .gt. 0 ) THEN
            WRITE(*,*) 'PB allocation in SMUMPS_185'
            id%INFO(1) = -13
            id%INFO(2) = NPROCS
            RETURN
         END IF
         LU_USAGE=DBLE(0)
         MD_MEM=DBLE(0)
      ENDIF
      IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
         ALLOCATE(CB_COST_MEM(2*2000*id%NSLAVES), 
     $            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_MEM=DBLE(0)
         ALLOCATE(CB_COST_ID(2000*3), 
     $            stat=allocok)
         IF (allocok > 0) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = id%NSLAVES 
           RETURN
         ENDIF
         CB_COST_ID=0
         POS_MEM=1
         POS_ID=1
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      ALLOCATE(FUTURE_NIV2(NPROCS), stat=allocok)
      IF (allocok > 0 ) THEN
         WRITE(*,*) 'PB allocation in SMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      ENDIF
      DO I = 1, NPROCS
        FUTURE_NIV2(I) = id%FUTURE_NIV2(I)
        IF(BDC_MD)THEN
           IF(FUTURE_NIV2(I).EQ.0)THEN
              MD_MEM(I-1)=999999999
           ENDIF
        ENDIF
      ENDDO
      DELTA_MEM=ZERO
      DELTA_LOAD=ZERO
#endif
      CHECK_MEM=0
#if defined(OLD_LOAD_MECHANISM)
#if defined (CHECK_COHERENCE) 
      NB_LEVEL2=0
      AMI_CHOSEN=.FALSE.
      IS_DISPLAYED=.FALSE.
#endif
#endif
      IF(BDC_SBTR.OR.BDC_POOL_MNG)THEN
         NB_SUBTREES=id%NBSA_LOCAL
         IF (ALLOCATED(MEM_SUBTREE)) DEALLOCATE(MEM_SUBTREE)
         ALLOCATE(MEM_SUBTREE(id%NBSA_LOCAL),stat=allocok)
         DO I=1,id%NBSA_LOCAL
            MEM_SUBTREE(I)=id%MEM_SUBTREE(I)
         ENDDO
         MY_FIRST_LEAF=>id%MY_FIRST_LEAF
         MY_NB_LEAF=>id%MY_NB_LEAF
         MY_ROOT_SBTR=>id%MY_ROOT_SBTR
         IF (ALLOCATED(SBTR_FIRST_POS_IN_POOL))
     $        DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         ALLOCATE(SBTR_FIRST_POS_IN_POOL(id%NBSA_LOCAL),stat=allocok)
         INSIDE_SUBTREE=0
         PEAK_SBTR_CUR_LOCAL = DBLE(0)
         SBTR_CUR_LOCAL      = DBLE(0)
         IF (ALLOCATED(SBTR_PEAK_ARRAY)) DEALLOCATE(SBTR_PEAK_ARRAY)
         ALLOCATE(SBTR_PEAK_ARRAY(id%NBSA_LOCAL),stat=allocok)
         SBTR_PEAK_ARRAY=DBLE(0)
         IF (ALLOCATED(SBTR_CUR_ARRAY)) DEALLOCATE(SBTR_CUR_ARRAY)
         ALLOCATE(SBTR_CUR_ARRAY(id%NBSA_LOCAL),stat=allocok)
         SBTR_CUR_ARRAY=DBLE(0)
         INDICE_SBTR_ARRAY=1
         NIV1_FLAG=0
         INDICE_SBTR=1
      ENDIF
      IF ( ALLOCATED(LOAD_FLOPS) ) DEALLOCATE( LOAD_FLOPS )
      ALLOCATE( LOAD_FLOPS( 0: NPROCS - 1 ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in SMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( ALLOCATED(WLOAD) ) DEALLOCATE( WLOAD )
      ALLOCATE( WLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in SMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( ALLOCATED(IDWLOAD) ) DEALLOCATE( IDWLOAD )
      ALLOCATE( IDWLOAD( NPROCS ), stat=allocok )
      IF ( allocok .gt. 0 ) THEN
         WRITE(*,*) 'PB allocation in SMUMPS_185'
         id%INFO(1) = -13
         id%INFO(2) = NPROCS
         RETURN
      END IF
      IF ( BDC_MEM ) THEN
        IF ( ALLOCATED(DM_MEM) ) DEALLOCATE( DM_MEM )
        ALLOCATE( DM_MEM( 0:NPROCS-1 ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
      END IF
      IF ( BDC_POOL ) THEN
        IF ( ALLOCATED(POOL_MEM) ) DEALLOCATE(POOL_MEM)
        ALLOCATE( POOL_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        POOL_MEM = DBLE(0) ! Array instruction
        POOL_LAST_COST_SENT = DBLE(0)
      END IF
      IF ( BDC_SBTR ) THEN
        IF ( ALLOCATED(SBTR_MEM) ) DEALLOCATE(SBTR_MEM)
        ALLOCATE( SBTR_MEM(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        IF ( ALLOCATED(SBTR_CUR) ) DEALLOCATE(SBTR_CUR)
        ALLOCATE( SBTR_CUR(0: NPROCS -1), stat=allocok)
        IF ( allocok .gt. 0 ) THEN
           WRITE(*,*) 'PB allocation in SMUMPS_185'
           id%INFO(1) = -13
           id%INFO(2) = NPROCS
           RETURN
        END IF
        SBTR_CUR = DBLE(0) ! Array instruction
        SBTR_MEM = DBLE(0) ! Array instruction
      END IF
      CALL SMUMPS_546(K34_LOC,K35_LOC)
      K35  = K35_LOC
      BUF_LOAD_SIZE = K34_LOC * 2 * ( NPROCS - 1 ) +
     *                NPROCS * ( K35_LOC + K34_LOC )
      IF (BDC_MEM) THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35_LOC
      END IF
      IF (BDC_SBTR)THEN
        BUF_LOAD_SIZE = BUF_LOAD_SIZE + NPROCS * K35_LOC
      ENDIF
      LBUF_LOAD_RECV = (BUF_LOAD_SIZE+K34_LOC)/K34_LOC
      LBUF_LOAD_RECV_BYTES = LBUF_LOAD_RECV * K34_LOC
      IF ( ALLOCATED(BUF_LOAD_RECV) ) DEALLOCATE(BUF_LOAD_RECV)
      ALLOCATE( BUF_LOAD_RECV( LBUF_LOAD_RECV), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in SMUMPS_185'
        id%INFO(1) = -13
        id%INFO(2) = LBUF_LOAD_RECV
        RETURN
      ENDIF 
      BUF_LOAD_SIZE = BUF_LOAD_SIZE * 20
      CALL SMUMPS_54( BUF_LOAD_SIZE, IERR )
      IF ( IERR .LT. 0 ) THEN
         id%INFO(1) = -13
         id%INFO(2) = BUF_LOAD_SIZE
         RETURN
      END IF
      DO I = 0, NPROCS - 1
         LOAD_FLOPS( I ) = ZERO
      END DO
#if defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MYID ) = COST_SUBTREE
      LAST_LOAD_SENT = ZERO
#endif
      IF ( BDC_MEM ) THEN
        DO I = 0, NPROCS - 1
          DM_MEM( I )=ZERO
        END DO
#if defined(OLD_LOAD_MECHANISM)
        DM_LAST_MEM_SENT=ZERO
#endif
      ENDIF
      CALL SMUMPS_425(KEEP(69))
      IF(BDC_MD)THEN
         MAX_SBTR=DBLE(0)
         IF(BDC_SBTR)THEN
            DO I=1,id%NBSA_LOCAL
               MAX_SBTR=MAX(id%MEM_SUBTREE(I),MAX_SBTR)
            ENDDO
         ENDIF
         MD_MEM(MYID)=MEMORY_MD
         WHAT=8
         CALL SMUMPS_460( WHAT,
     *        COMM_LD, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *        FUTURE_NIV2,
#endif
     *        DBLE(MEMORY_MD),DBLE(0) ,MYID, iERR  )
         WHAT=9
         MEMORY_SENT = DBLE(LA-MAX_SURF_MASTER-MAX_SBTR)
     *      - MAX( DBLE(LA) * DBLE(3) / DBLE(100),
     *      DBLE(2) * DBLE(MAX(KEEP(5),KEEP(6))) * DBLE(KEEP(127)))
         IF (KEEP(12) > 25) THEN
           MEMORY_SENT = MEMORY_SENT -
     *                   DBLE(KEEP(12))*0.2d0*DBLE(LA)/100.0d0
         ENDIF
         TAB_MAXS(MYID)=INT(MEMORY_SENT)
         CALL SMUMPS_460( WHAT,
     *        COMM_LD, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *        FUTURE_NIV2,
#endif
     *        MEMORY_SENT,
     $        DBLE(0),MYID, iERR  )
      ENDIF
      END SUBROUTINE SMUMPS_185
      SUBROUTINE SMUMPS_190( CHECK_FLOPS,PROCESS_BANDE,
     $     INC_LOAD, KEEP,KEEP8 )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      DOUBLE PRECISION INC_LOAD
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      LOGICAL PROCESS_BANDE
      INTEGER CHECK_FLOPS
      INTEGER IERR
      DOUBLE PRECISION ZERO, SEND_MEM, SEND_LOAD,SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC MAX, ABS
      IF (INC_LOAD == 0) THEN
         IF(REMOVE_NODE_FLAG)THEN
            REMOVE_NODE_FLAG=.FALSE.
         ENDIF
         RETURN
      ENDIF
      IF((CHECK_FLOPS.NE.0).AND.
     $     (CHECK_FLOPS.NE.1).AND.(CHECK_FLOPS.NE.2))THEN
         WRITE(*,*)MYID,': Bad value for CHECK_FLOPS'
         CALL SMUMPS_ABORT()
      ENDIF
      IF(CHECK_FLOPS.EQ.1)THEN
         CHK_LD=CHK_LD+INC_LOAD
      ELSE 
         IF(CHECK_FLOPS.EQ.2)THEN
            RETURN
         ENDIF
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE ) THEN
         RETURN                 ! Info sent by master of node already
      ENDIF
#endif
      LOAD_FLOPS( MYID ) = MAX( LOAD_FLOPS( MYID ) + INC_LOAD, ZERO)
      IF(BDC_M2_FLOPS.AND.REMOVE_NODE_FLAG)THEN
         IF(INC_LOAD.NE.REMOVE_NODE_COST)THEN
            IF(INC_LOAD.GT.REMOVE_NODE_COST)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_LOAD = DELTA_LOAD +
     $              (INC_LOAD-REMOVE_NODE_COST)
               GOTO 888
#else
               GOTO 888
#endif               
            ELSE
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_LOAD = DELTA_LOAD -
     $              (REMOVE_NODE_COST-INC_LOAD)
               GOTO 888
#else
               GOTO 888
#endif               
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      DELTA_LOAD = DELTA_LOAD + INC_LOAD
 888  CONTINUE
      IF ( DELTA_LOAD > MIN_DIFF .OR. DELTA_LOAD < -MIN_DIFF) THEN
         SEND_LOAD = DELTA_LOAD
         IF (BDC_MEM) THEN
           SEND_MEM = DELTA_MEM
         ELSE
           SEND_MEM = ZERO
         END IF
#else
 888   CONTINUE
       IF ( ABS( LOAD_FLOPS ( MYID ) -
     *      LAST_LOAD_SENT ).GT.MIN_DIFF)THEN
          IERR = 0
          SEND_LOAD  = LOAD_FLOPS( MYID )
          IF ( BDC_MEM ) THEN
             SEND_MEM = DM_MEM(MYID)
          ELSE
             SEND_MEM = ZERO
          END IF
#endif
         IF(BDC_SBTR)THEN
           SBTR_TMP=SBTR_CUR(MYID)
         ELSE
           SBTR_TMP=0
         ENDIF
 111     CONTINUE
         CALL SMUMPS_77( BDC_SBTR,BDC_MEM,
     *        BDC_MD,COMM_LD, NPROCS,
     *        SEND_LOAD,
     *        SEND_MEM,SBTR_TMP,
     *        DM_SUMLU,
#if ! defined(OLD_LOAD_MECHANISM)
     *        FUTURE_NIV2,
#endif
     *        MYID, IERR )
           IF ( IERR == -1 )THEN
             CALL SMUMPS_467(COMM_LD, KEEP)
             GOTO 111
           ELSE IF ( IERR .NE.0 ) THEN
             WRITE(*,*) "Internal Error in SMUMPS_190",IERR
             CALL SMUMPS_ABORT()
           ENDIF
         IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
           DELTA_LOAD = ZERO
           IF (BDC_MEM) DELTA_MEM  = ZERO
#else          
           LAST_LOAD_SENT = LOAD_FLOPS( MYID )
           IF ( BDC_MEM ) DM_LAST_MEM_SENT = DM_MEM( MYID )
#endif
         END IF
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG)THEN
         REMOVE_NODE_FLAG=.FALSE.
      ENDIF
      END SUBROUTINE SMUMPS_190
      SUBROUTINE SMUMPS_471( SSARBR,
     *           PROCESS_BANDE_ARG, MEM_VALUE, NEW_LU, INC_MEM_ARG,
     *           KEEP,KEEP8,LRLU)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MEM_VALUE, INC_MEM_ARG, NEW_LU,LRLU
      LOGICAL, INTENT(IN) :: PROCESS_BANDE_ARG, SSARBR
      INTEGER IERR, KEEP(500)
      INTEGER*8 KEEP8(150)
      DOUBLE PRECISION ZERO, SEND_MEM, TMP_MEM, SBTR_TMP
      PARAMETER( ZERO=0.0d0 )
      INTRINSIC MAX, ABS
      INTEGER INC_MEM
      LOGICAL PROCESS_BANDE
      INTEGER I
      PROCESS_BANDE=PROCESS_BANDE_ARG
      INC_MEM = INC_MEM_ARG
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( PROCESS_BANDE .AND. NEW_LU .NE. 0) THEN
        WRITE(*,*) " Internal Error in SMUMPS_471."
        WRITE(*,*) " NEW_LU must be zero if called from PROCESS_BANDE"
        CALL SMUMPS_ABORT()
      ENDIF
#endif
      DM_SUMLU = DM_SUMLU + DBLE(NEW_LU)
      CHECK_MEM = CHECK_MEM + INC_MEM
#if ! defined (SIMULATE_OOC_FAC)
      IF (KEEP(201).NE.0) THEN
#endif
        GOTO 777
#if ! defined (SIMULATE_OOC_FAC)
      ENDIF
#endif
      IF ( MEM_VALUE .NE. CHECK_MEM ) THEN
         WRITE(*,*)MYID,
     *   ':Problem with increments in SMUMPS_471',
     *   CHECK_MEM, MEM_VALUE, INC_MEM
         CALL SMUMPS_ABORT()
      ENDIF
 777  CONTINUE
#if ! defined(OLD_LOAD_MECHANISM)
      IF (PROCESS_BANDE) THEN
         RETURN
      ENDIF
#endif
      IF(BDC_POOL_MNG) THEN
         IF(SBTR_WHICH_M.EQ.0)THEN
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+INC_MEM-NEW_LU
         ELSE
            IF (SSARBR) SBTR_CUR_LOCAL = SBTR_CUR_LOCAL+INC_MEM
         ENDIF
      ENDIF
      IF ( .NOT. BDC_MEM ) THEN
         RETURN
      ENDIF
#if defined(OLD_LOAD_MECHANISM)
      DM_MEM( MYID ) = DBLE(CHECK_MEM) - DM_SUMLU
      TMP_MEM = DM_MEM(MYID)
#endif
      IF (BDC_SBTR .AND. SSARBR) THEN
         IF(SBTR_WHICH_M.EQ.0)THEN
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+INC_MEM-NEW_LU
         ELSE
            SBTR_CUR(MYID) = SBTR_CUR(MYID)+INC_MEM
         ENDIF
        SBTR_TMP = SBTR_CUR(MYID)
      ELSE
        SBTR_TMP=0
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF ( NEW_LU > 0 ) THEN
        INC_MEM = INC_MEM - NEW_LU
      ENDIF
      DM_MEM( MYID ) = DM_MEM(MYID) + INC_MEM
      MAX_PEAK_STK=MAX(MAX_PEAK_STK,DM_MEM(MYID))
      IF(BDC_M2_MEM.AND.REMOVE_NODE_FLAG_MEM)THEN
         IF(INC_MEM.NE.REMOVE_NODE_COST_MEM)THEN
            IF(INC_MEM.GT.REMOVE_NODE_COST_MEM)THEN
               DELTA_MEM = DELTA_MEM +
     $              (INC_MEM-REMOVE_NODE_COST_MEM)
               GOTO 888               
            ELSE
               DELTA_MEM = DELTA_MEM -
     $              (REMOVE_NODE_COST_MEM-INC_MEM)
               GOTO 888
            ENDIF
         ENDIF
         GOTO 333
      ENDIF
      DELTA_MEM = DELTA_MEM + DBLE(INC_MEM)
 888  CONTINUE
      IF ((KEEP(48).NE.5).OR.
     $     ((KEEP(48).EQ.5).AND.(ABS(DELTA_MEM)
     $      .GE.0.1d0*DBLE(LRLU))))THEN
         IF ( ABS(DELTA_MEM) > DM_THRES_MEM ) THEN
            SEND_MEM = DELTA_MEM
#else
      IF(BDC_M2_MEM.AND.REMOVE_NODE_FLAG_MEM)THEN
         IF(INC_MEM.NE.REMOVE_NODE_COST_MEM)THEN
            IF(INC_MEM.EQ.REMOVE_NODE_COST_MEM)THEN               
               GOTO 333
            ELSEIF(INC_MEM.LT.REMOVE_NODE_COST_MEM)THEN
               GOTO 333
            ENDIF
         ENDIF
      ENDIF
      IF ((KEEP(48).NE.5).OR.
     $     ((KEEP(48).EQ.5).AND.
     $      (ABS(TMP_MEM-DM_LAST_MEM_SENT).GE.
     $     0.1d0*DBLE(LRLU))))THEN
         IF ( ABS( TMP_MEM-DM_LAST_MEM_SENT) >
     *     DM_THRES_MEM )  THEN
            IERR = 0
            SEND_MEM = TMP_MEM
#endif
 111        CONTINUE
            CALL SMUMPS_77( 
     *           BDC_SBTR,
     *           BDC_MEM,BDC_MD, COMM_LD,
     *           NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *           DELTA_LOAD,
#else
     *           LOAD_FLOPS( MYID ), 
#endif
     *           SEND_MEM,SBTR_TMP,
     *           DM_SUMLU,
#if ! defined(OLD_LOAD_MECHANISM)
     *           FUTURE_NIV2,
#endif
     $           MYID,IERR )
            IF ( IERR == -1 )THEN
               CALL SMUMPS_467(COMM_LD, KEEP)
               GOTO 111
           ELSE IF ( IERR .NE. 0 ) THEN
             WRITE(*,*) "Internal Error in SMUMPS_471",IERR
             CALL SMUMPS_ABORT()
           ENDIF
           IF ( IERR .EQ. 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
             DELTA_LOAD = ZERO
             DELTA_MEM  = ZERO
#else
             LAST_LOAD_SENT   = LOAD_FLOPS  ( MYID )
             DM_LAST_MEM_SENT = TMP_MEM
#endif
           END IF
         ENDIF
      ENDIF
 333  CONTINUE
      IF(REMOVE_NODE_FLAG_MEM)THEN
         REMOVE_NODE_FLAG_MEM=.FALSE.
      ENDIF
      END SUBROUTINE SMUMPS_471
      INTEGER FUNCTION SMUMPS_186( K69, MEM_DISTRIB,MSG_SIZE )
      IMPLICIT NONE
      INTEGER I, NLESS, K69 
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION LREF
      DOUBLE PRECISION MSG_SIZE
      NLESS = 0
      DO I=1,NPROCS
         IDWLOAD(I) = I - 1
      ENDDO
      WLOAD(1:NPROCS) = LOAD_FLOPS(0:NPROCS-1)
      IF(BDC_M2_FLOPS)THEN
         DO I=1,NPROCS
            WLOAD(I)=WLOAD(I)+NIV2(I)
         ENDDO
      ENDIF
      IF(K69 .gt. 1) THEN
         CALL SMUMPS_426(MEM_DISTRIB,MSG_SIZE,IDWLOAD,NPROCS)
      ENDIF
      LREF = LOAD_FLOPS(MYID)
      DO I=1, NPROCS
         IF (WLOAD(I).LT.LREF) NLESS=NLESS+1
      ENDDO
      SMUMPS_186 = NLESS
      RETURN
      END FUNCTION SMUMPS_186
      SUBROUTINE SMUMPS_189(MEM_DISTRIB,MSG_SIZE,DEST,
     *     NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES
      INTEGER DEST(NSLAVES)
      INTEGER, DIMENSION(0:NPROCS - 1) :: MEM_DISTRIB
      INTEGER I,J,NBDEST
      DOUBLE PRECISION MSG_SIZE
      IF ( NSLAVES.eq.NPROCS-1 ) THEN
        J = MYID+1
        DO I=1,NSLAVES
           J=J+1
           IF (J.GT.NPROCS) J=1
           DEST(I) = J - 1
        ENDDO
      ELSE
        DO I=1,NPROCS
           IDWLOAD(I) = I - 1
        ENDDO
        CALL SMUMPS_558(NPROCS, WLOAD, IDWLOAD)
         NBDEST = 0
         DO I=1, NSLAVES
            J = IDWLOAD(I)
            IF (J.NE.MYID) THEN
               NBDEST = NBDEST+1
               DEST(NBDEST) = J
            ENDIF
         ENDDO
         IF (NBDEST.NE.NSLAVES) THEN
            DEST(NSLAVES) = IDWLOAD(NSLAVES+1)
         ENDIF
         IF(BDC_MD)THEN
            J=NSLAVES+1
            do i=NSLAVES+1,NPROCS
               IF(IDWLOAD(I).NE.MYID)THEN
                  dest(J)= IDWLOAD(i)
                  J=J+1
               ENDIF
            end do
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_189
      SUBROUTINE SMUMPS_184( IPROC, LOAD_IPROC )
      IMPLICIT NONE
      INTEGER, intent(in) :: IPROC
      DOUBLE PRECISION, intent(out) :: LOAD_IPROC
      LOAD_IPROC = LOAD_FLOPS( IPROC )
      RETURN
      END SUBROUTINE SMUMPS_184
      SUBROUTINE SMUMPS_183( INFO1, IERR )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER INFO1, IERR
      DEALLOCATE( LOAD_FLOPS )
      DEALLOCATE( WLOAD )
      DEALLOCATE( IDWLOAD )
#if ! defined(OLD_LOAD_MECHANISM)
      DEALLOCATE(FUTURE_NIV2)
#endif
      IF(BDC_MD)THEN
         DEALLOCATE(MD_MEM)
         DEALLOCATE(LU_USAGE)
         DEALLOCATE(TAB_MAXS)
      ENDIF
      IF ( BDC_MEM ) DEALLOCATE( DM_MEM )
      IF ( BDC_POOL) DEALLOCATE( POOL_MEM )
      IF ( BDC_SBTR) THEN
         DEALLOCATE( SBTR_MEM )
         DEALLOCATE( SBTR_CUR )
         DEALLOCATE(SBTR_FIRST_POS_IN_POOL)
         NULLIFY(MY_FIRST_LEAF)
         NULLIFY(MY_NB_LEAF)
         NULLIFY(MY_ROOT_SBTR)
      ENDIF
      IF(KEEP_LOAD(76).EQ.4)THEN
         NULLIFY(DEPTH_FIRST_LOAD)
      ENDIF
      IF(KEEP_LOAD(76).EQ.5)THEN
         NULLIFY(COST_TRAV)
      ENDIF
      IF (BDC_M2_MEM.OR.BDC_M2_FLOPS) THEN
        DEALLOCATE(NB_SON,POOL_NIV2,POOL_NIV2_COST, NIV2)
      END IF
      IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
         DEALLOCATE(CB_COST_MEM)
         DEALLOCATE(CB_COST_ID)
      ENDIF
      NULLIFY(ND_LOAD)
      NULLIFY(KEEP_LOAD)
      NULLIFY(KEEP8_LOAD)
      NULLIFY(FILS_LOAD)
      NULLIFY(FRERE_LOAD)
      NULLIFY(PROCNODE_LOAD)
      NULLIFY(STEP_LOAD)
      NULLIFY(NE_LOAD)
      NULLIFY(CAND_LOAD)
      NULLIFY(STEP_TO_NIV2_LOAD)
      NULLIFY(DAD_LOAD)
      IF (BDC_SBTR.OR.BDC_POOL_MNG) THEN
         DEALLOCATE(MEM_SUBTREE)
         DEALLOCATE(SBTR_PEAK_ARRAY)
         DEALLOCATE(SBTR_CUR_ARRAY)
      ENDIF
      CALL SMUMPS_58( IERR )
      CALL SMUMPS_150( MYID, COMM_LD,
     *     BUF_LOAD_RECV, LBUF_LOAD_RECV,
     *     LBUF_LOAD_RECV_BYTES )
      DEALLOCATE(BUF_LOAD_RECV)
      END SUBROUTINE SMUMPS_183
#if defined (LAMPORT_)
      RECURSIVE SUBROUTINE SMUMPS_467(COMM, KEEP)
#else
      SUBROUTINE SMUMPS_467(COMM, KEEP)
#endif      
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, MSGTAG, MSGLEN, MSGSOU,COMM
      INTEGER KEEP(500)
      INTEGER STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
 10   CONTINUE
      CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     *     FLAG, STATUS, IERR )
      IF (FLAG) THEN
        KEEP(65)=KEEP(65)+1
        MSGTAG = STATUS( MPI_TAG )
        MSGSOU = STATUS( MPI_SOURCE )
        IF ( MSGTAG .NE. UPDATE_LOAD) THEN
          write(*,*) "Internal error 1 in SMUMPS_467",
     *    MSGTAG
          CALL SMUMPS_ABORT()
        ENDIF
        CALL MPI_GET_COUNT(STATUS, MPI_PACKED, MSGLEN, IERR)
        IF ( MSGLEN > LBUF_LOAD_RECV_BYTES ) THEN
          write(*,*) "Internal error 2 in SMUMPS_467",
     *    MSGLEN, LBUF_LOAD_RECV_BYTES
          CALL SMUMPS_ABORT()
        ENDIF
        CALL MPI_RECV( BUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES,
     *    MPI_PACKED, MSGSOU, MSGTAG, COMM_LD, STATUS, IERR)
        CALL SMUMPS_187( MSGSOU, BUF_LOAD_RECV,
     *  LBUF_LOAD_RECV, LBUF_LOAD_RECV_BYTES )
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_467
      RECURSIVE SUBROUTINE SMUMPS_187
     *   ( MSGSOU, BUFR, LBUFR, LBUFR_BYTES )
      IMPLICIT NONE
      INTEGER MSGSOU, LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INCLUDE 'mpif.h'
      INTEGER POSITION, IERR, WHAT, NSLAVES, I
      DOUBLE PRECISION LOAD_RECEIVED
      INTEGER INODE_RECEIVED,K,J,SON,NSLAVES_TEMP,POS_TEMP,
     $     NCB_RECEIVED,SURF
      INTEGER, POINTER, DIMENSION (:) :: LIST_SLAVES
      DOUBLE PRECISION, POINTER, DIMENSION (:) :: LOAD_INCR
      LOGICAL FLAG_LEADER
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     WHAT, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
      IF ( WHAT == 0 ) THEN
#if ! defined(OLD_LOAD_MECHANISM)
#else
#endif
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
      LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
#else
      LOAD_FLOPS( MSGSOU ) = LOAD_RECEIVED
#endif
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     *       COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
          DM_MEM(MSGSOU)  = DM_MEM(MSGSOU) + LOAD_RECEIVED
#else
          DM_MEM(MSGSOU)  = LOAD_RECEIVED
#endif
          MAX_PEAK_STK=MAX(MAX_PEAK_STK,DM_MEM(MSGSOU))
        END IF
        IF(BDC_SBTR)THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *       LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     *       COMM_LD, IERR )
          SBTR_CUR(MSGSOU)=LOAD_RECEIVED
        ENDIF
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          LOAD_RECEIVED, 1, MPI_DOUBLE_PRECISION,
     *          COMM_LD, IERR )
           LU_USAGE(MSGSOU)=LOAD_RECEIVED
        ENDIF
      ELSEIF (( WHAT == 1).OR.(WHAT.EQ.19)) THEN
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     NSLAVES, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     INODE_RECEIVED, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LIST_SLAVES, NSLAVES, MPI_INTEGER,
     *     COMM_LD, IERR)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_INCR, NSLAVES, MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR)
        DO I = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
          IF ( LIST_SLAVES(I) /= MYID ) THEN
#endif
            LOAD_FLOPS(LIST_SLAVES(I)) =
     *      LOAD_FLOPS(LIST_SLAVES(I)) + 
     *      LOAD_INCR(I)
#if defined(OLD_LOAD_MECHANISM)
          END IF
#endif
        END DO
        IF ( BDC_MEM ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_INCR, NSLAVES, MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR)
          DO I = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
            IF ( LIST_SLAVES(I) /= MYID ) THEN
#endif
              DM_MEM(LIST_SLAVES(I)) = DM_MEM(LIST_SLAVES(I)) + 
     *        LOAD_INCR(I)
              MAX_PEAK_STK=MAX(MAX_PEAK_STK,DM_MEM(LIST_SLAVES(I)))
#if defined(OLD_LOAD_MECHANISM)
            END IF
#endif
          END DO
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          LOAD_INCR, NSLAVES, MPI_DOUBLE_PRECISION,
     *          COMM_LD, IERR)
           CALL CLEAN_POOL_MEM_INFO(INODE_RECEIVED)
           CB_COST_ID(POS_ID)=INODE_RECEIVED
           CB_COST_ID(POS_ID+1)=NSLAVES
           CB_COST_ID(POS_ID+2)=POS_MEM
           POS_ID=POS_ID+3
           DO I=1,NSLAVES
              WRITE(*,*)MYID,':',LIST_SLAVES(I),'->',LOAD_INCR(I)
              CB_COST_MEM(POS_MEM)=LIST_SLAVES(I)
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=LOAD_INCR(I)
              POS_MEM=POS_MEM+1
           ENDDO
        ENDIF
#if defined(OLD_LOAD_MECHANISM)
#if defined (CHECK_COHERENCE) 
 344    CONTINUE
#endif
#endif
        NULLIFY( LIST_SLAVES )
        NULLIFY( LOAD_INCR )
      ELSE IF (WHAT == 2 ) THEN
        IF ( .not. BDC_POOL ) THEN
          WRITE(*,*) "Internal error 2 in SMUMPS_187"
          CALL SMUMPS_ABORT()
        END IF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
        POOL_MEM(MSGSOU)=LOAD_RECEIVED
      ELSE IF ( WHAT == 3 ) THEN
        IF ( .NOT. BDC_SBTR) THEN
          WRITE(*,*) "Internal error 3 in SMUMPS_187"
          CALL SMUMPS_ABORT()
        ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
        SBTR_MEM(MSGSOU)=SBTR_MEM(MSGSOU)+LOAD_RECEIVED
#if ! defined(OLD_LOAD_MECHANISM)
      ELSE IF (WHAT == 4) THEN
        FUTURE_NIV2(MSGSOU+1)=0
        IF(BDC_MD)THEN
           CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *          SURF, 1, MPI_INTEGER,
     *          COMM_LD, IERR )
          MD_MEM(MSGSOU)=999999999
          TAB_MAXS(MSGSOU)=TAB_MAXS(MSGSOU)+SURF
        ENDIF
#endif
        IF(BDC_M2_MEM.OR.BDC_M2_FLOPS)THEN
        ENDIF
      ELSE IF (WHAT == 5) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 7 in SMUMPS_187"
            CALL SMUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     INODE_RECEIVED, 1,
     *     MPI_INTEGER,
     *     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            CALL PROCESS_NIV2_MEM_MSG(INODE_RECEIVED)
         ELSEIF(BDC_M2_FLOPS) THEN
            CALL PROCESS_NIV2_FLOPS_MSG(INODE_RECEIVED)
         ENDIF
         IF((KEEP_LOAD(81).EQ.2).OR.(KEEP_LOAD(81).EQ.3))THEN
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *           INODE_RECEIVED, 1,
     *           MPI_INTEGER,
     *           COMM_LD, IERR )   
               CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *              NCB_RECEIVED, 1,
     *              MPI_INTEGER,
     *              COMM_LD, IERR )   
            IF(SMUMPS_330(STEP_LOAD(INODE_RECEIVED),
     $           PROCNODE_LOAD,NPROCS).EQ.1)THEN
               CB_COST_ID(POS_ID)=INODE_RECEIVED
               CB_COST_ID(POS_ID+1)=1
               CB_COST_ID(POS_ID+2)=POS_MEM
               POS_ID=POS_ID+3
               CB_COST_MEM(POS_MEM)=MSGSOU
               POS_MEM=POS_MEM+1
               CB_COST_MEM(POS_MEM)=NCB_RECEIVED*NCB_RECEIVED
               POS_MEM=POS_MEM+1
            ENDIF
         ENDIF
      ELSE IF ( WHAT == 6 ) THEN
         IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
            WRITE(*,*) "Internal error 8 in SMUMPS_187"
            CALL SMUMPS_ABORT()
         ENDIF
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED
            IF(NIV2(MSGSOU+1).LT.0)THEN
               IF(ABS(NIV2(MSGSOU+1)).LE.
     $              SQRT(EPSILON(LOAD_RECEIVED)))THEN
                  NIV2(MSGSOU+1)=DBLE(0)
               ELSE
                  WRITE(*,*)'problem with NIV2_FLOPS message',
     $                 NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
                  CALL SMUMPS_ABORT()
               ENDIF
            ENDIF
         ENDIF
      ELSEIF(WHAT == 17)THEN
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
         IF(BDC_M2_MEM) THEN
            NIV2(MSGSOU+1) = LOAD_RECEIVED
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *           LOAD_RECEIVED, 1,
     *           MPI_DOUBLE_PRECISION,
     *           COMM_LD, IERR )
            IF(BDC_MD)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DM_MEM(MYID)=DM_MEM(MYID)+LOAD_RECEIVED
#else
               DM_MEM(MYID)=LOAD_RECEIVED
#endif
            ELSEIF(BDC_POOL)THEN
               POOL_MEM(MSGSOU)=LOAD_RECEIVED
            ENDIF
         ELSEIF(BDC_M2_FLOPS) THEN
            NIV2(MSGSOU+1) = NIV2(MSGSOU+1) + LOAD_RECEIVED            
            IF(NIV2(MSGSOU+1).LT.0)THEN
               WRITE(*,*)'problem with NIV2_FLOPS message',
     $              NIV2(MSGSOU+1),MSGSOU,LOAD_RECEIVED
               CALL SMUMPS_ABORT
            ENDIF
            CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *           LOAD_RECEIVED, 1,
     *           MPI_DOUBLE_PRECISION,
     *           COMM_LD, IERR )
#if ! defined(OLD_LOAD_MECHANISM)
            LOAD_FLOPS( MSGSOU ) = LOAD_FLOPS(MSGSOU) + LOAD_RECEIVED
#else
            LOAD_FLOPS( MSGSOU ) = LOAD_RECEIVED
#endif
         ENDIF
      ELSEIF ( WHAT == 7 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 4
     *in SMUMPS_187'
            CALL SMUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     NSLAVES, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     INODE_RECEIVED, 1, MPI_INTEGER,
     *     COMM_LD, IERR )
        LIST_SLAVES => IDWLOAD
        LOAD_INCR => WLOAD
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LIST_SLAVES, NSLAVES, MPI_INTEGER,
     *     COMM_LD, IERR)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_INCR, NSLAVES, MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR)
        DO I = 1, NSLAVES
#if defined(OLD_LOAD_MECHANISM)
          IF ( LIST_SLAVES(I) /= MYID ) THEN
#endif
            MD_MEM(LIST_SLAVES(I)) =
     *      MD_MEM(LIST_SLAVES(I)) + 
     *      LOAD_INCR(I)
#if ! defined(OLD_LOAD_MECHANISM)
            IF(FUTURE_NIV2(LIST_SLAVES(I)+1).EQ.0)THEN
               MD_MEM(LIST_SLAVES(I))=999999999
            ENDIF
#endif
#if defined(OLD_LOAD_MECHANISM)
          END IF
#endif
        END DO
      ELSEIF ( WHAT == 8 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 5
     *in SMUMPS_187'
            CALL SMUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
        MD_MEM(MSGSOU)=MD_MEM(MSGSOU)+LOAD_RECEIVED
#if ! defined(OLD_LOAD_MECHANISM)
        IF(FUTURE_NIV2(MSGSOU+1).EQ.0)THEN
           MD_MEM(MSGSOU)=999999999
        ENDIF
#endif
      ELSEIF ( WHAT == 9 ) THEN
         IF(.NOT.BDC_MD)THEN
            WRITE(*,*)MYID,': Internal error 6
     *in SMUMPS_187'
            CALL SMUMPS_ABORT()
         ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *     LOAD_RECEIVED, 1,
     *     MPI_DOUBLE_PRECISION,
     *     COMM_LD, IERR )
        TAB_MAXS(MSGSOU)=INT(LOAD_RECEIVED)
      ELSE
          WRITE(*,*) "Internal error 1 in SMUMPS_187"
          CALL SMUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE SMUMPS_187
      integer function SMUMPS_409
     *                 (MEM_DISTRIB,CAND,
     *                  K69,
     *                  SLAVEF,MSG_SIZE,
     *                  NMB_OF_CAND )
      implicit none      
      integer, intent(in) :: K69, SLAVEF
      INTEGER, intent(in) :: CAND(SLAVEF+1)
      INTEGER, DIMENSION(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      INTEGER, intent(out) :: NMB_OF_CAND
      integer i,nless
      DOUBLE PRECISION lref
      DOUBLE PRECISION MSG_SIZE
      nless = 0 
      NMB_OF_CAND=CAND(SLAVEF+1)
      do i=1,NMB_OF_CAND
         WLOAD(I)=LOAD_FLOPS(CAND(I))
         IF(BDC_M2_FLOPS)THEN
            WLOAD(I)=WLOAD(I)+NIV2(CAND(I)+1)
         ENDIF
      end do
      IF(K69 .gt. 1) THEN
         CALL SMUMPS_426(MEM_DISTRIB,MSG_SIZE,
     *        cand,NMB_OF_CAND)
      ENDIF
      lref = LOAD_FLOPS(myid)
      do i=1, NMB_OF_CAND
         if (wload(i).lt.lref) nless=nless+1
      end do 
      SMUMPS_409 = nless
      return
      end function SMUMPS_409
      subroutine SMUMPS_384
     *           (MEM_DISTRIB,CAND, 
     *
     *            SLAVEF,
     *            nslaves_inode, dest)
      implicit none
      integer, intent(in) :: nslaves_inode, SLAVEF
      integer, intent(in) :: CAND(SLAVEF+1)
      integer, dimension(0:NPROCS - 1), intent(in) :: MEM_DISTRIB
      integer, intent(out) :: dest(CAND(SLAVEF+1))
      integer i,j,NMB_OF_CAND,deficit,nextproc
      DOUBLE PRECISION lref
      external SMUMPS_558
      NMB_OF_CAND = CAND(SLAVEF+1)
      if(nslaves_inode.ge.NPROCS .or.
     *   nslaves_inode.gt.NMB_OF_CAND) then
         write(*,*)'Internal error in SMUMPS_384',
     *   nslaves_inode, nprocs, NMB_OF_CAND
         CALL SMUMPS_ABORT()
      end if
      if (nslaves_inode.eq.NPROCS-1) then
         j=myid+1
         do i=1,nslaves_inode
            if(j.ge.nprocs) j=0
            dest(i)=j
            j=j+1
         end do
      else
        do i=1,NMB_OF_CAND
               IDWLOAD(i)=i
        end do
        call SMUMPS_558(NMB_OF_CAND,
     *       wload(1:NMB_OF_CAND),IDWLOAD(1:NMB_OF_CAND) )
        do i=1,nslaves_inode
           dest(i)= cand(IDWLOAD(i))
        end do
        IF(BDC_MD)THEN
           do i=nslaves_inode+1,NMB_OF_CAND
              dest(i)= cand(IDWLOAD(i))
           end do
        ENDIF
      end if   
      return
      end subroutine SMUMPS_384
      SUBROUTINE SMUMPS_425(K69)
      IMPLICIT NONE
      INTEGER K69
      IF (K69 .LE. 4) THEN
         ALPHA = 0.0d0
         BETA = 0.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 5) THEN
         ALPHA = 0.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 6) THEN
         ALPHA = 0.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 7) THEN
         ALPHA = 0.5d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 8) THEN
         ALPHA = 1.0d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 9) THEN
         ALPHA = 1.0d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 10) THEN
         ALPHA = 1.0d0
         BETA = 150000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 11) THEN
         ALPHA = 1.5d0
         BETA = 50000.0d0
         RETURN
      ENDIF
      IF (K69 .EQ. 12) THEN
         ALPHA = 1.5d0
         BETA = 100000.0d0
         RETURN
      ENDIF
      ALPHA = 1.5d0
      BETA = 150000.0d0
      RETURN
      END SUBROUTINE SMUMPS_425
      SUBROUTINE SMUMPS_426(MEM_DISTRIB,MSG_SIZE,ARRAY_ADM,LEN)
      IMPLICIT NONE
      INTEGER I,J,LEN,COUNTMORE,COUNTLOCLESS
      INTEGER, DIMENSION(0:NPROCS-1) :: MEM_DISTRIB
      DOUBLE PRECISION MSG_SIZE,FORBIGMSG
      INTEGER ARRAY_ADM(LEN)
      DOUBLE PRECISION MY_LOAD
      FORBIGMSG = 1.0d0
      IF (K69 .lt.2) THEN
         RETURN
      ENDIF
      IF(BDC_M2_FLOPS)THEN
         MY_LOAD=LOAD_FLOPS(MYID)+NIV2(MYID+1)
      ELSE
         MY_LOAD=LOAD_FLOPS(MYID)
      ENDIF
      IF((MSG_SIZE * K35 ) .gt. 3200000.0) THEN
         FORBIGMSG = 2.0d0
      ENDIF
      IF (K69 .le. 4) THEN
         DO I = 1,LEN
            IF ((MEM_DISTRIB(ARRAY_ADM(I)) .EQ. 1) .AND.
     *      WLOAD(I) .LT. MY_LOAD ) THEN
               WLOAD(I) = WLOAD(I)/MY_LOAD
            ELSE
              IF ( MEM_DISTRIB(ARRAY_ADM(I)) .NE. 1 ) THEN
                WLOAD(I) = WLOAD(I) *
     *              MEM_DISTRIB(ARRAY_ADM(I))
     *              * FORBIGMSG
     *              + DBLE(2)
              ENDIF
            ENDIF
         ENDDO
         RETURN
      ENDIF
      DO I = 1,LEN
         IF ((MEM_DISTRIB(ARRAY_ADM(I)) .EQ. 1) .AND.
     *        WLOAD(I) .LT. MY_LOAD ) THEN
            WLOAD(I) = WLOAD(I) /  MY_LOAD
         ELSE
            IF(MEM_DISTRIB(ARRAY_ADM(I)) .NE. 1) THEN     
               WLOAD(I) = (WLOAD(I) +
     *              ALPHA * MSG_SIZE * K35  +
     *              BETA) * FORBIGMSG
            ENDIF
         ENDIF
      ENDDO  
      RETURN
      END SUBROUTINE SMUMPS_426
      SUBROUTINE SMUMPS_461(MYID, SLAVEF, COMM,
     *     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES, NSLAVES,INODE)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: MYID, SLAVEF, COMM, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, INTENT (IN) :: LIST_SLAVES( NSLAVES )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER NCB, NFRONT, NBROWS_SLAVE
      INTEGER I, IERR,WHAT,INODE
      DOUBLE PRECISION MEM_INCREMENT ( NSLAVES )
      DOUBLE PRECISION FLOPS_INCREMENT( NSLAVES )
      DOUBLE PRECISION CB_BAND( NSLAVES )
      IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
         WHAT=1
      ELSE
         WHAT=19
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      FUTURE_NIV2(MYID+1) = FUTURE_NIV2(MYID+1) - 1
      IF ( FUTURE_NIV2(MYID+1) < 0 ) THEN
        WRITE(*,*) "Internal error in SMUMPS_461"
        CALL SMUMPS_ABORT()
      ENDIF
      IF ( FUTURE_NIV2(MYID + 1) == 0 ) THEN
 112    CONTINUE
        CALL SMUMPS_502(COMM,MYID,SLAVEF,
     *       MAX_SURF_MASTER,IERR)
        IF (IERR == -1 ) THEN
          CALL SMUMPS_467(COMM_LD, KEEP)
          GOTO 112
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in SMUMPS_461",
     *    IERR
          CALL SMUMPS_ABORT()
        ENDIF
      TAB_MAXS(MYID) = TAB_MAXS(MYID) + MAX_SURF_MASTER
      ENDIF
#endif
      IF ( NSLAVES /= TAB_POS(SLAVEF + 2) ) THEN
        write(*,*) "Error 1 in SMUMPS_461",
     +             NSLAVES, TAB_POS(SLAVEF+2)
        stop
      ENDIF
      NCB = TAB_POS(NSLAVES+1) - 1
      NFRONT = NCB + NASS
      DO I = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(I+1) - TAB_POS(I)
         IF ( KEEP(50) == 0 ) THEN
            FLOPS_INCREMENT( I ) = (DBLE(NBROWS_SLAVE)*DBLE( NASS ))+
     *           DBLE(NBROWS_SLAVE) * DBLE(NASS) *
     $           DBLE(2*NFRONT-NASS-1)
         ELSE
            FLOPS_INCREMENT( I ) = DBLE(NBROWS_SLAVE) * DBLE(NASS ) *
     *           DBLE( 2 * ( NASS + TAB_POS(I+1) - 1 ) 
     *           - NBROWS_SLAVE - NASS + 1 )
         ENDIF
         IF ( BDC_MEM ) THEN
            IF ( KEEP(50) == 0 ) THEN
               MEM_INCREMENT( I ) = DBLE(NBROWS_SLAVE) *
     $              DBLE(NFRONT)
            ELSE
               MEM_INCREMENT( I ) = DBLE(NBROWS_SLAVE) *
     *              DBLE( NASS + TAB_POS(I+1) - 1 )
            END IF
         ENDIF
         IF((KEEP(81).NE.2).AND.(KEEP(81).NE.3))THEN
            CB_BAND(I)=-999999
         ELSE
            IF ( KEEP(50) == 0 ) THEN
               CB_BAND( I ) = DBLE(NBROWS_SLAVE) *
     $              DBLE(NFRONT-NASS)
            ELSE
               CB_BAND( I ) = DBLE(NBROWS_SLAVE) *
     *              DBLE(TAB_POS(I+1)-1)
            END IF
         ENDIF
      END DO
      IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
         CB_COST_ID(POS_ID)=INODE
         CB_COST_ID(POS_ID+1)=NSLAVES
         CB_COST_ID(POS_ID+2)=POS_MEM
         POS_ID=POS_ID+3
         DO I=1,NSLAVES
            CB_COST_MEM(POS_MEM)=LIST_SLAVES(I)
            POS_MEM=POS_MEM+1
            CB_COST_MEM(POS_MEM)=CB_BAND(I)
            POS_MEM=POS_MEM+1
         ENDDO
      ENDIF
 111  CONTINUE
      CALL SMUMPS_524(BDC_MEM, COMM, MYID, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *     FUTURE_NIV2,
#endif
     *     NSLAVES, LIST_SLAVES,INODE,
     *     MEM_INCREMENT,
     *     FLOPS_INCREMENT,CB_BAND, WHAT,IERR)
        IF ( IERR == -1 ) THEN
          CALL SMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in SMUMPS_461",
     *    IERR
          CALL SMUMPS_ABORT()
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
#endif
        DO I = 1, NSLAVES
          LOAD_FLOPS(LIST_SLAVES(I)) = LOAD_FLOPS(LIST_SLAVES(I))
     *       +  FLOPS_INCREMENT(I)
          IF ( BDC_MEM ) THEN
            DM_MEM(LIST_SLAVES(I)) = DM_MEM(LIST_SLAVES(I))
     *       +  MEM_INCREMENT(I)
          END IF
        ENDDO
#if ! defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      RETURN
      END SUBROUTINE SMUMPS_461
      SUBROUTINE SMUMPS_500(
     $     POOL, LPOOL,
     *     PROCNODE, KEEP,KEEP8, SLAVEF, COMM, MYID, STEP, N,
     *     ND, FILS )
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL, SLAVEF, COMM, MYID 
      INTEGER N, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER POOL( LPOOL ), PROCNODE( KEEP(28) ), STEP( N )
      INTEGER ND( KEEP(28) ), FILS( N )
      INTEGER I, INODE, NELIM, NFR, LEVEL, IERR, WHAT
      DOUBLE PRECISION COST
      INTEGER NBINSUBTREE,NBTOP,INSUBTREE
      LOGICAL SMUMPS_167
      INTEGER SMUMPS_330
      EXTERNAL SMUMPS_167, SMUMPS_330
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF(BDC_MD)THEN
         RETURN
      ENDIF
      IF((KEEP(76).EQ.0).OR.(KEEP(76).EQ.2))THEN
         IF(NBTOP.NE.0)THEN
            DO I = LPOOL-NBTOP-2, MIN(LPOOL-3,LPOOL-NBTOP-2+3)
               INODE = POOL( I )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=DBLE(0) 
            GOTO 30
         ELSE
            DO I = NBINSUBTREE, MAX(1,NBINSUBTREE-3), -1
               INODE = POOL( I )
               IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                  GOTO 20
               END IF
            END DO
            COST=DBLE(0) 
            GOTO 30
         ENDIF
      ELSE
         IF(KEEP(76).EQ.1)THEN
            IF(INSUBTREE.EQ.1)THEN
               DO I = NBINSUBTREE, MAX(1,NBINSUBTREE-3), -1
                  INODE = POOL( I )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=DBLE(0) 
               GOTO 30
            ELSE
               DO I = LPOOL-NBTOP-2, MIN(LPOOL-3,LPOOL-NBTOP-2+3)
                  INODE = POOL( I )
                  IF (INODE .LE. N .AND. INODE .GE. 1 ) THEN
                     GOTO 20
                  END IF
               END DO
               COST=DBLE(0) 
               GOTO 30
            ENDIF
         ELSE
            WRITE(*,*)
     *      'Internal error: Unknown pool management strategy'
            CALL SMUMPS_ABORT()
         ENDIF
      ENDIF
 20   CONTINUE
        I = INODE
        NELIM = 0
 10     CONTINUE
        IF ( I > 0 ) THEN
          NELIM = NELIM + 1
          I = FILS(I)
          GOTO 10
        ENDIF
        NFR = ND( STEP(INODE) )
        LEVEL = SMUMPS_330( STEP(INODE), PROCNODE, SLAVEF )
        IF (LEVEL .EQ. 1) THEN
          COST = DBLE( NFR ) * DBLE( NFR )
        ELSE
          IF ( KEEP(50) == 0 ) THEN
            COST = DBLE( NFR ) * DBLE( NELIM )
          ELSE
            COST = DBLE( NELIM ) * DBLE( NELIM )
          ENDIF
        ENDIF
 30   CONTINUE
      IF ( ABS(POOL_LAST_COST_SENT-COST).GT.DM_THRES_MEM ) THEN
        WHAT = 2
 111    CONTINUE
        CALL SMUMPS_460( WHAT,
     *         COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *               FUTURE_NIV2,
#endif
     *         COST, DBLE(0),MYID, iERR  )
        POOL_LAST_COST_SENT = COST
        POOL_MEM(MYID)=COST
        IF ( IERR == -1 )THEN
          CALL SMUMPS_467(COMM_LD, KEEP)
          GOTO 111
        ELSE IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) "Internal Error in SMUMPS_500",
     *    IERR
          CALL SMUMPS_ABORT()
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_500
      SUBROUTINE SMUMPS_501(
     $     OK,INODE,POOL,LPOOL,MYID,slavef,COMM,KEEP,KEEP8)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER LPOOL,MYID,SLAVEF,COMM,INODE
      INTEGER POOL(LPOOL),KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER WHAT,IERR
      LOGICAL OK
      DOUBLE PRECISION ZERO,COST
      LOGICAL FLAG
      PARAMETER(ZERO = 0.0d0)
      EXTERNAL SMUMPS_283,SMUMPS_170
      LOGICAL SMUMPS_283,SMUMPS_170
      IF((INODE.LE.0).OR.(INODE.GT.N_LOAD)) THEN
         RETURN
      ENDIF
      IF(.NOT.SMUMPS_170(STEP_LOAD(INODE),
     $     PROCNODE_LOAD,NPROCS))THEN
         RETURN         
      ENDIF
      IF(SMUMPS_283(STEP_LOAD(INODE),PROCNODE_LOAD,NPROCS))THEN
         IF(NE_LOAD(STEP_LOAD(INODE)).EQ.0)THEN
            RETURN
         ENDIF
      ENDIF
      FLAG=.FALSE.
      IF(INDICE_SBTR.LE.NB_SUBTREES)THEN
         IF(INODE.EQ.MY_FIRST_LEAF(INDICE_SBTR))THEN
            FLAG=.TRUE.
         ENDIF
      ENDIF
      IF(FLAG)THEN
         SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)=MEM_SUBTREE(INDICE_SBTR)
         SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)=SBTR_CUR(MYID)
         INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY+1
         WHAT = 3
         IF(MEM_SUBTREE(INDICE_SBTR).GE.DM_THRES_MEM)THEN
 111        CONTINUE
            CALL SMUMPS_460(
     *           WHAT, COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *           FUTURE_NIV2,
#endif
     *           DBLE(MEM_SUBTREE(INDICE_SBTR)), DBLE(0),MYID, iERR  )
            IF ( IERR == -1 )THEN
               CALL SMUMPS_467(COMM_LD, KEEP)
               GOTO 111
            ELSE IF ( IERR .NE. 0 ) THEN
               WRITE(*,*)
     *         "Internal Error 1 in SMUMPS_501",
     *         IERR
               CALL SMUMPS_ABORT()
            ENDIF
         ENDIF
         SBTR_MEM(MYID)=SBTR_MEM(MYID)+
     $        DBLE(MEM_SUBTREE(INDICE_SBTR))
         INDICE_SBTR=INDICE_SBTR+1
         IF(INSIDE_SUBTREE.EQ.0)THEN
            INSIDE_SUBTREE=1
         ENDIF
      ELSE 
         IF(INODE.EQ.MY_ROOT_SBTR(INDICE_SBTR-1))THEN
            WHAT = 3
            COST=-SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY-1)
            IF(ABS(COST).GE.DM_THRES_MEM)THEN
 112           CONTINUE
               CALL SMUMPS_460(
     *              WHAT, COMM, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *              FUTURE_NIV2,
#endif
     *              COST, DBLE(0) ,MYID,iERR  )
               IF ( IERR == -1 )THEN
                  CALL SMUMPS_467(COMM_LD, KEEP)
                  GOTO 112
               ELSE IF ( IERR .NE. 0 ) THEN
                  WRITE(*,*)
     *        "Internal Error 3 in SMUMPS_501",
     *        IERR
                  CALL SMUMPS_ABORT()
               ENDIF
            ENDIF
            INDICE_SBTR_ARRAY=INDICE_SBTR_ARRAY-1            
            SBTR_MEM(MYID)=SBTR_MEM(MYID)-
     $           SBTR_PEAK_ARRAY(INDICE_SBTR_ARRAY)
            SBTR_CUR(MYID)=SBTR_CUR_ARRAY(INDICE_SBTR_ARRAY)
            IF(INDICE_SBTR_ARRAY.EQ.1)THEN
               SBTR_CUR(MYID)=0
               INSIDE_SUBTREE=0
            ENDIF
         ENDIF
         ENDIF 
 666  CONTINUE
      END SUBROUTINE SMUMPS_501
      SUBROUTINE SMUMPS_504
     &      (SLAVEF,KEEP,KEEP8,PROCS,MEM_DISTRIB,NCB,NFRONT,
     $       NSLAVES_NODE,TAB_POS,
     $       SLAVES_LIST,SIZE_SLAVES_LIST,MYID)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER*8 KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1)
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47,K10, K48, K50
      INTEGER J, NSLAVES, NSLAVES_REAL, ISLAVE
      INTEGER KMIN, KMAX, BLSIZE, NPIV,
     *        KMINSURF, NSLAVES_max
      DOUBLE PRECISION WK_MASTER, WK_SLAVE
      DOUBLE PRECISION  SMUMPS_45
      EXTERNAL SMUMPS_45
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,AFFECTED
      INTEGER ADDITIONNAL_ROWS,BLOCK_SIZE,I,X,TOTAL_MEM,REF,POS
      LOGICAL force_cand
      DOUBLE PRECISION TEMP(SLAVEF),SUB(SLAVEF),PEAK,MAX_TEMP
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF),SUB_ID(SLAVEF),M2(SLAVEF)
      EXTERNAL MPI_WTIME
      DOUBLE PRECISION MPI_WTIME
      IF (KEEP(10) .GT. 0) THEN
      write(*,*)MYID,
     * ": Internal Error 1 in SMUMPS_504"
      CALL SMUMPS_ABORT()
      ENDIF
      K10=ABS(KEEP(10))
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
      END IF
      IF(K48.NE.4)THEN
         WRITE(*,*)'SMUMPS_COMPUTE_PARTI_ACTV_MEM_K10
     $      should be called with KEEP(48) different from 4'
         CALL SMUMPS_ABORT()
      ENDIF
         KMIN=1
         KMAX=K10/NFRONT
         IF(FORCE_CAND)THEN
            DO I=1,PROCS(SLAVEF+1)
               WLOAD(I)=DM_MEM(PROCS(I))
               IDWLOAD(I)=PROCS(I)
            ENDDO
            NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = DM_MEM(0:NUMBER_OF_PROCS-1)
            DO I=1,NUMBER_OF_PROCS
               IDWLOAD(I) = I - 1
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         NB_ROWS=0
         CALL SMUMPS_558(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         TOTAL_MEM=NCB*NFRONT
         SOMME=0
         J=1
         PEAK=DBLE(0)
         DO I=1,NUMBER_OF_PROCS
            IF((IDWLOAD(I).NE.MYID))THEN
               PEAK=MAX(PEAK,WLOAD(I))
               TEMP_ID(J)=IDWLOAD(I)
               TEMP(J)=WLOAD(I)
                IF(BDC_SBTR)THEN
                   TEMP(J)=TEMP(J)+SBTR_MEM(IDWLOAD(I))-
     $                  SBTR_CUR(IDWLOAD(I))
                ENDIF
                IF(BDC_POOL)THEN
                   TEMP(J)=TEMP(J)+POOL_MEM(TEMP_ID(J))
                ENDIF
                IF(BDC_M2_MEM)THEN
                   TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
                ENDIF
                J=J+1
            ENDIF
         ENDDO
         NUMBER_OF_PROCS=J-1
         CALL SMUMPS_558(NUMBER_OF_PROCS, TEMP, TEMP_ID)
         IF(K50.EQ.0)THEN
           PEAK=MAX(PEAK,
     *       DM_MEM(MYID)+DBLE(NFRONT)*DBLE(NFRONT-NCB))
         ELSE
           PEAK=MAX(PEAK,
     *       DM_MEM(MYID)+DBLE(NFRONT-NCB)*DBLE(NFRONT-NCB))
         ENDIF
         PEAK=MAX(PEAK,TEMP(OTHERS))
         SOMME=0
         DO I=1,NUMBER_OF_PROCS
           SOMME=SOMME+TEMP(OTHERS)-TEMP(I)
         ENDDO               
         IF(SOMME.LE.DBLE(TOTAL_MEM)) THEN
            GOTO 096
         ENDIF
 096     CONTINUE
         SOMME=0
         DO I=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(I)
         ENDDO
         IF(DBLE(TOTAL_MEM).GE.SOMME) THEN
 887        CONTINUE     
            AFFECTED=0
            CHOSEN=0
            ACC=0
            DO I=1,OTHERS
               IF(K50.EQ.0)THEN
                  IF((TEMP(OTHERS)-TEMP(I)).GT.K10)THEN
                     AFFECTED=K10
                  ELSE
                     AFFECTED=TEMP(OTHERS)-TEMP(I)
                  ENDIF
                  x=INT(AFFECTED/DBLE(NFRONT))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
               ENDIF
               IF(K50.NE.0)THEN
                  IF((TEMP(OTHERS)-TEMP(I)).GT.K10)THEN
                     AFFECTED=K10
                  ELSE
                     AFFECTED=TEMP(OTHERS)-TEMP(I)
                  ENDIF
                  x=INT((-(NFRONT-NCB+ACC)
     $                 +SQRT(DBLE(((NFRONT-NCB+ACC)*
     $                 (NFRONT-NCB+ACC))+DBLE(4)*
     *                 (AFFECTED))))/
     $                 DBLE(2))
                  IF((ACC+X).GT.NCB) X=NCB-ACC
                  IF(X.LE.0) THEN
                     WRITE(*,*)"Internal Error 2 in
     $                    SMUMPS_504"
                     CALL SMUMPS_ABORT()
                  ENDIF
               ENDIF
               NB_ROWS(I)=X
               CHOSEN=CHOSEN+1
               ACC=ACC+X
               IF(NCB-ACC.LT.KMIN) GOTO 111
               IF(NCB.EQ.ACC) GOTO 111
               ENDDO
 111           CONTINUE
               IF((ACC.GT.NCB))THEN
                  X=0
                  DO I=1,OTHERS
                     X=X+NB_ROWS(I)
                  ENDDO
                  WRITE(*,*)'NCB=',NCB,',SOMME=',X
                  WRITE(*,*)MYID,
     $               ": Internal Error 3 in SMUMPS_504"
                  CALL SMUMPS_ABORT()
               ENDIF
               IF((NCB.NE.ACC))THEN
                  IF(K50.NE.0)THEN
                     IF(CHOSEN.NE.0)THEN
                        ADDITIONNAL_ROWS=NCB-ACC
                        NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+ADDITIONNAL_ROWS
                     ELSE
                        AFFECTED=TOTAL_MEM/DBLE(NUMBER_OF_PROCS)
                        CHOSEN=0
                        ACC=0
                        DO I=1,OTHERS
                           x=INT((-(NFRONT-NCB+ACC)
     $                          +SQRT(DBLE(((NFRONT-NCB+ACC)*
     $                          (NFRONT-NCB+ACC))+DBLE(4)*
     *                          (AFFECTED))))/
     $                          DBLE(2))
                           IF((ACC+X).GT.NCB) X=NCB-ACC
                           NB_ROWS(I)=X
                           CHOSEN=CHOSEN+1
                           ACC=ACC+X
                           IF(NCB-ACC.LT.KMIN) GOTO 002
                           IF(NCB.EQ.ACC) GOTO 002
                        ENDDO
 002                    CONTINUE
                        IF(ACC.LT.NCB)THEN
                           NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+(NCB-ACC)
                        ENDIF
                     ENDIF
                     GOTO 333
                  ENDIF
                  ADDITIONNAL_ROWS=NCB-ACC
                  DO I=CHOSEN,1,-1
                     IF(INT(DBLE(ADDITIONNAL_ROWS)/
     $                    DBLE(I)).NE.0)THEN
                        GOTO 222 
                     ENDIF
                  ENDDO
 222              CONTINUE
                  X=INT(DBLE(ADDITIONNAL_ROWS)/DBLE(I))
                  DO J=1,I
                     NB_ROWS(J)=NB_ROWS(J)+X
                     ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                  ENDDO 
                  IF(ADDITIONNAL_ROWS.NE.0) THEN
                     NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                  ENDIF    
               ENDIF
 333           CONTINUE
               IF(NB_ROWS(CHOSEN).EQ.0) CHOSEN=CHOSEN-1
               GOTO 889
            ELSE
               DO I=OTHERS,1,-1
                  SOMME=0
                  DO J=1,I
                     SOMME=SOMME+TEMP(J)
                  ENDDO
                  SOMME=(DBLE(I)*TEMP(I))-SOMME
                  IF(TOTAL_MEM.GE.INT(SOMME)) GOTO 444
               ENDDO
 444           CONTINUE
               REF=I
               DO J=1,I
                  IF(TEMP(J).EQ.TEMP(I)) THEN
                     SMALL_SET=J
                     GOTO 123
                  ENDIF
               ENDDO
 123           CONTINUE
               IF(I.EQ.1)THEN 
                  NB_ROWS(I)=NCB
                  CHOSEN=1
                  GOTO 666
               ENDIF
 323           CONTINUE
               AFFECTED=0
               CHOSEN=0
               ACC=0
               DO I=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(I)).GT.K10)THEN
                        AFFECTED=K10
                     ELSE
                        AFFECTED=TEMP(SMALL_SET)-TEMP(I)
                     ENDIF
                     x=INT(AFFECTED/DBLE(NFRONT))
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                     IF((TEMP(SMALL_SET)-TEMP(I)).GT.K10)THEN
                        AFFECTED=K10
                     ELSE
                        AFFECTED=TEMP(SMALL_SET)-TEMP(I)
                     ENDIF
                      x=INT((-(NFRONT-NCB+ACC)
     $                  +SQRT(DBLE(((NFRONT-NCB+ACC)*
     $                  (NFRONT-NCB+ACC))+DBLE(4)*
     *                  (AFFECTED))))/
     $                  DBLE(2))
                     IF(X.LT.0)THEN
                        WRITE(*,*)MYID,
     *             ': Internal error 4 in SMUMPS_504'
                        CALL SMUMPS_ABORT()
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  NB_ROWS(I)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     *            ': Internal error 5 in SMUMPS_504'
                    CALL SMUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=0
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     *           ':Internal error 6 in SMUMPS_504'
                  CALL SMUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LT.OTHERS)THEN
                       SMALL_SET=REF+1
                       REF=SMALL_SET
                       GOTO 323
                     ELSE
                       NB_ROWS(CHOSEN)=NB_ROWS(CHOSEN)+NCB-ACC
                       GOTO 666
                     ENDIF
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
#if ! defined (OLD_PART)
                 I=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                .AND.(I.LE.NUMBER_OF_PROCS))
                    J=1
#if ! defined (PART1_)
                    X=INT(ADDITIONNAL_ROWS/(I-1))
                    IF((X.EQ.0).AND.(ADDITIONNAL_ROWS.NE.0))THEN
                       DO WHILE ((J.LT.I).AND.(ADDITIONNAL_ROWS.GT.0))
                         NB_ROWS(J)=NB_ROWS(J)+1
                         ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                         J=J+1
                       ENDDO
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          WRITE(*,*)MYID,
     *             ':Internal error 7 in SMUMPS_504'
                         CALL SMUMPS_ABORT()
                      ENDIF
                      GOTO 047
                    ENDIF
                    IF((TEMP(1)+DBLE((NB_ROWS(1)+X)*NFRONT)).LE.
     $                   TEMP(I))THEN
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                      .AND.(J.LT.I))
                          AFFECTED=X
                          IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
                    ELSE
#endif
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                      .AND.(J.LE.I))
                          AFFECTED=INT((TEMP(I)-(TEMP(J)+
     $                         DBLE(NB_ROWS(J)*
     *                         NFRONT)))/DBLE(NFRONT))
                          IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(J)
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
#if ! defined (PART1_)
                    ENDIF
#endif
                    I=I+1
                 ENDDO
 047             CONTINUE
                 IF((ADDITIONNAL_ROWS.EQ.0).AND.
     $                (I.LT.NUMBER_OF_PROCS))THEN
                    CHOSEN=I-1
                 ELSE
                    CHOSEN=I-2
                 ENDIF
#if ! defined (PART1_)
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     $                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO I=1,CHOSEN
                       NB_ROWS(I)=NB_ROWS(I)+1
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
#endif
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     $                (ADDITIONNAL_ROWS.NE.0))THEN
                    I=CHOSEN+1
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                   .AND.(I.LE.NUMBER_OF_PROCS))
                       J=1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                      .AND.(J.LE.I))
                          AFFECTED=INT((TEMP(I)-(TEMP(J)+
     $                         DBLE(NB_ROWS(J)*
     *                         NFRONT)))/DBLE(NFRONT))
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          J=J+1
                       ENDDO
                       I=I+1
                    ENDDO
                    CHOSEN=I-2
                 ENDIF
 049             CONTINUE
#else
                 DO I=CHOSEN,1,-1
                    IF(INT(DBLE(ADDITIONNAL_ROWS)/
     $                   DBLE(I)).NE.0)THEN
                       GOTO 555
                    ENDIF
                 ENDDO
 555             CONTINUE
                 X=INT(DBLE(ADDITIONNAL_ROWS)/DBLE(I))
                 DO J=1,I
                    IF(NB_ROWS(J)+X.GT.K10/NCB)THEN
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                      ((K10/NCB)-NB_ROWS(J))
                       NB_ROWS(J)=(K10/NFRONT)
                    ELSE IF ((TEMP(J)+DBLE((NB_ROWS(J)+x))*
     $                      DBLE(NFRONT)).GT.
     *                 PEAK)THEN
                       AFFECTED=INT((PEAK-(TEMP(J)+DBLE(NB_ROWS(J))*
     *                 DBLE(NFRONT))/DBLE(NFRONT))
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                       NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                    ELSE
                       NB_ROWS(J)=NB_ROWS(J)+X
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    ENDIF
                    IF(((TEMP(J)+DBLE(NB_ROWS(J))*DBLE(NFRONT))
     *               .GT. PEAK)
     *               .AND.(SMALL_SET.LT.OTHERS))THEN
                       WRITE(*,*)MYID,
     *           ':Internal error 8 in SMUMPS_504'
                       SMALL_SET=SMALL_SET+1
                       CALL SMUMPS_ABORT()
                    ENDIF
                 ENDDO
                 SOMME=0
                 DO J=1,CHOSEN
                   SOMME=SOMME+NB_ROWS(J)
                 ENDDO
                 IF(ADDITIONNAL_ROWS.NE.0) THEN
                       DO J=1,CHOSEN
                          IF(NB_ROWS(J).LT.0)THEN
                            WRITE(*,*)MYID,
     *           ':Internal error 9 in SMUMPS_504'
                            CALL SMUMPS_ABORT()
                          ENDIF
                          IF ((TEMP(J)+DBLE(NB_ROWS(J))
     *                 *DBLE(NFRONT)).GT.
     *                 PEAK)THEN
                             WRITE(*,*)MYID,
     *         ':Internal error 10 in SMUMPS_504'
                          CALL SMUMPS_ABORT()
                       ENDIF
                          IF ((TEMP(J)+DBLE(NB_ROWS(J)+
     *                 ADDITIONNAL_ROWS)*DBLE(NFRONT)).GE.
     *                 PEAK)THEN
                            AFFECTED=INT((PEAK-(TEMP(J)+
     *                      DBLE(NB_ROWS(J))*
     *                      DBLE(NFRONT))/DBLE(NFRONT))
                            ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                            NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                            IF((TEMP(J)+DBLE(NFRONT)*
     $                           DBLE(NB_ROWS(J))).GT.
     *                         PEAK)THEN
                               WRITE(*,*)MYID,
     *         ':Internal error 11 in SMUMPS_504'
                              CALL SMUMPS_ABORT()
                            ENDIF  
                          ELSE
                            NB_ROWS(J)=NB_ROWS(J)+ADDITIONNAL_ROWS
                            ADDITIONNAL_ROWS=0
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0) GOTO 666
                       ENDDO
                       IF(SMALL_SET.EQ.(NUMBER_OF_PROCS)) THEN
                         NB_ROWS=0
                         GOTO 887
                       ENDIF
                       I=CHOSEN+1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                        .AND.(I.LE.NUMBER_OF_PROCS))
                          IF(NB_ROWS(I)+ADDITIONNAL_ROWS.LT.K10/NFRONT)
     $                         THEN
                             NB_ROWS(I)=NB_ROWS(I)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-(K10/
     $                            NFRONT
     $                            -NB_ROWS(I))
                             NB_ROWS(I)=K10/NFRONT
                          ENDIF
                          I=I+1
                       ENDDO
                       CHOSEN=I-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                         DO I=CHOSEN,1,-1
                          IF(INT(DBLE(ADDITIONNAL_ROWS)/DBLE(I))
     *                          .NE.0)THEN
                             GOTO 372 
                          ENDIF
                         ENDDO
 372                     CONTINUE
                         X=INT(DBLE(ADDITIONNAL_ROWS)/DBLE(I))
                         DO J=1,I
                           NB_ROWS(J)=NB_ROWS(J)+X
                           ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                         ENDDO 
                         IF(ADDITIONNAL_ROWS.NE.0) THEN
                           NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                         ENDIF    
                    ENDIF
                    ENDIF
#endif
              ENDIF
 666          CONTINUE
              SOMME=0
              X=0
              POS=0
              DO I=1,CHOSEN
                 IF(K50.NE.0) THEN
                    IF((TEMP(I)+DBLE(NB_ROWS(I))
     *                   *DBLE(X+NB_ROWS(I)+NFRONT-NCB))
     *                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 IF(K50.EQ.0) THEN
                    IF((TEMP(I)+DBLE(NB_ROWS(I))*DBLE(NFRONT))
     *                   .GT.PEAK)THEN
                       SMALL_SET=SMALL_SET+1
                    ENDIF
                 ENDIF
                 X=X+NB_ROWS(I)
                 SOMME=SOMME+ DBLE(NB_ROWS(I))
              ENDDO
           ENDIF
 889       CONTINUE
           J=CHOSEN
           X=0
           DO I=J,1,-1
             IF(NB_ROWS(I).EQ.0)THEN
                IF(X.EQ.1)THEN
                  WRITE(*,*)MYID,
     *         ':Internal error 12 in SMUMPS_504'
                  CALL SMUMPS_ABORT()
                ENDIF
                CHOSEN=CHOSEN-1
             ELSE 
                  IF(NB_ROWS(I).GT.0)THEN
                    X=1
                  ELSE
                    WRITE(*,*)
     *            'Internal error 13 in SMUMPS_504'
                    CALL SMUMPS_ABORT()
                  ENDIF
             ENDIF
          ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           DO I=1,CHOSEN
              SLAVES_LIST(I)=TEMP_ID(I)
              TAB_POS(I)=POS
              POS=POS+NB_ROWS(I) 
              IF(NB_ROWS(I).LE.0)THEN
                WRITE(*,*)
     *          'Internal error 14 in SMUMPS_504'
                 CALL SMUMPS_ABORT()
              ENDIF
           ENDDO
           DO I=CHOSEN+1,NUMBER_OF_PROCS
              SLAVES_LIST(I)=TEMP_ID(I)
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)
     *        'Internal error 15 in SMUMPS_504'
             CALL SMUMPS_ABORT()
           ENDIF
      END SUBROUTINE SMUMPS_504
      SUBROUTINE SMUMPS_518
     &      (NCBSON_MAX,SLAVEF,KEEP,KEEP8,
     *       PROCS,MEM_DISTRIB,NCB,NFRONT,
     $       NSLAVES_NODE,TAB_POS,
     $       SLAVES_LIST,SIZE_SLAVES_LIST,MYID,INODE,MP,LP)
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER*8 KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: NCBSON_MAX
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(in) :: MP,LP
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE
      INTEGER NUMBER_OF_PROCS,K47,K10, K48, K50,K83,K69
      INTEGER J, NSLAVES, NSLAVES_REAL, ISLAVE
      INTEGER KMIN, KMAX, BLSIZE, NPIV,
     *        KMINSURF, NSLAVES_max
      DOUBLE PRECISION WK_MASTER, WK_SLAVE
      DOUBLE PRECISION  SMUMPS_45
      EXTERNAL SMUMPS_45
      INTEGER OTHERS,CHOSEN,SMALL_SET,ACC
      DOUBLE PRECISION SOMME,AFFECTED,DELTA,A,B,C,MASTER_WORK
      INTEGER ADDITIONNAL_ROWS,BLOCK_SIZE,I,X,REF,POS,NELIM,NB_DIFF
      LOGICAL force_cand,SMP
      INTEGER BANDE_K10,NB_SAT,NB_ZERO
      DOUBLE PRECISION TEMP(SLAVEF),SUB(SLAVEF),PEAK,TOTAL_COST,
     $     MAX_MEM_ALLOW,MAX_FLOPS_ALLOW
      INTEGER TEMP_ID(SLAVEF),NB_ROWS(SLAVEF),M2(SLAVEF),
     $     NB_ROWS_SMP(SLAVEF)!,MIN_GRAN
      INTEGER NSLAVES_REF,NCB_FILS
      EXTERNAL MPI_WTIME,SMUMPS_442
      INTEGER SMUMPS_442
      INTEGER POS_MIN_LOAD,SIZE_MY_SMP,WHAT,ADDITIONNAL_ROWS_SPECIAL
      LOGICAL HAVE_TYPE1_SON
      DOUBLE PRECISION MIN_LOAD,MAX_LOAD,TEMP_MAX_LOAD
      DOUBLE PRECISION MPI_WTIME
      DOUBLE PRECISION BUF_SIZE,NELIM_MEM_SIZE
      DOUBLE PRECISION MEM_SIZE_STRONG(SLAVEF),MEM_SIZE_WEAK(SLAVEF)
      INTEGER SMUMPS_330
      EXTERNAL SMUMPS_330
      K10=ABS(KEEP(10))
      TEMP_MAX_LOAD=DBLE(0)
      K50=KEEP(50)
      K48=KEEP(48)
      K47=KEEP(47)
      K83=KEEP(83)
      K69=0
      NCB_FILS=NCBSON_MAX
      IF(NCB_FILS*MIN(NCB,NCB_FILS).GT.K10)THEN
         HAVE_TYPE1_SON=.TRUE.
      ELSE
         HAVE_TYPE1_SON=.FALSE.
      ENDIF
      SMP=(K69.NE.0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
      END IF
      NELIM=NFRONT-NCB
         KMAX=K10/NCB
         IF(FORCE_CAND)THEN
           DO I=1,PROCS(SLAVEF+1)
              WLOAD(I)=LOAD_FLOPS(PROCS(I))
              IDWLOAD(I)=PROCS(I)
              IF (WLOAD(I) < -0.5d0 ) THEN
                 IF((MP.GT.0).AND.(LP.GE.2))THEN
                    WRITE(MP,*)MYID,': Warning: negative load ',
     $                   WLOAD(I)
                 ENDIF
              ENDIF 
              WLOAD(I)=MAX(WLOAD(I),0.0d0)
           ENDDO
           NUMBER_OF_PROCS=PROCS(SLAVEF+1)
            OTHERS=NUMBER_OF_PROCS
         ELSE
            NUMBER_OF_PROCS=SLAVEF
            WLOAD(1:SLAVEF) = LOAD_FLOPS(0:NUMBER_OF_PROCS-1)
            DO I=1,NUMBER_OF_PROCS
               IDWLOAD(I) = I - 1
               IF (WLOAD(I) < -0.5d0 ) THEN
                  IF((MP.GT.0).AND.(LP.GE.2))THEN
                     WRITE(MP,*)MYID,': Negative load ',
     $                    WLOAD(I)
                  ENDIF
               ENDIF 
               WLOAD(I)=MAX(WLOAD(I),0.0d0)
            ENDDO
            OTHERS=NUMBER_OF_PROCS-1
         ENDIF
         KMAX=INT(NCB/OTHERS)
         KMIN=SMUMPS_442(NCB*KMAX,K50,KMAX,NCB)
         NB_ROWS=0
         CALL SMUMPS_558(NUMBER_OF_PROCS, WLOAD, IDWLOAD)
         IF(K50.EQ.0)THEN
            TOTAL_COST=DBLE( NELIM ) * DBLE ( NCB ) +
     *           DBLE(NCB) * DBLE(NELIM)*DBLE(2*NFRONT-NELIM-1)
         ELSE
            TOTAL_COST=DBLE(NELIM) * DBLE ( NCB ) *
     *           DBLE(NFRONT+1)
         ENDIF
         CALL SMUMPS_511(NFRONT,NELIM,NELIM,K50,
     $        2,MASTER_WORK)
         SOMME=0
         J=1
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.GT.K83))THEN
            MASTER_WORK=DBLE(KEEP(88))*MASTER_WORK/DBLE(100)
         ENDIF
         IF(FORCE_CAND.AND.(NUMBER_OF_PROCS.LE.K83))THEN
            MASTER_WORK=DBLE(KEEP(87))*MASTER_WORK/DBLE(100)
         ENDIF
         IF(MASTER_WORK.LT.DBLE(1))THEN
            MASTER_WORK=DBLE(1)
         ENDIF
         NSLAVES_REF=INT(TOTAL_COST/MASTER_WORK)+1
         IF(FORCE_CAND)THEN
            NSLAVES_REF=MIN(NSLAVES_REF,NUMBER_OF_PROCS)
         ELSE
            NSLAVES_REF=MIN(NSLAVES_REF,NUMBER_OF_PROCS-1)
         ENDIF
        DO I=1,NUMBER_OF_PROCS
           IF((IDWLOAD(I).NE.MYID))THEN
              TEMP_ID(J)=IDWLOAD(I)
              TEMP(J)=WLOAD(I)
              IF(BDC_M2_FLOPS)THEN
                 TEMP(J)=TEMP(J)+NIV2(TEMP_ID(J)+1)
              ENDIF
              J=J+1
           ENDIF
        ENDDO               
        NUMBER_OF_PROCS=J-1
        CALL SMUMPS_558(NUMBER_OF_PROCS, TEMP, TEMP_ID)
        SOMME=0
        AFFECTED=0
        DO I=1,OTHERS
            SOMME=SOMME+TEMP(OTHERS)-TEMP(I)
            AFFECTED=AFFECTED+TEMP(I)
        ENDDO
         AFFECTED=(AFFECTED/DBLE(OTHERS))+
     $        (TOTAL_COST/DBLE(OTHERS))
         SIZE_MY_SMP=OTHERS
         MIN_LOAD=TEMP(1)
         POS_MIN_LOAD=1
         IF(.NOT.SMP) MAX_LOAD=TEMP(OTHERS)
         IF(SMP)THEN
            J=1
            DO I=1,OTHERS
               IF(MEM_DISTRIB(TEMP_ID(I)).EQ.1)THEN
                  IF(TEMP(I).LE.AFFECTED)THEN
                     WLOAD(J)=TEMP(I)
                     IDWLOAD(J)=TEMP_ID(I)
                     J=J+1
                  ELSE
                  ENDIF
               ENDIF
            ENDDO
            MAX_LOAD=WLOAD(J-1)
            SIZE_MY_SMP=J-1
            DO I=1,OTHERS
               IF((MEM_DISTRIB(TEMP_ID(I)).NE.1).OR.
     $              ((MEM_DISTRIB(TEMP_ID(I)).EQ.1).AND.
     $              (TEMP(I).GE.AFFECTED)))THEN
                  WLOAD(J)=TEMP(I)
                  IDWLOAD(J)=TEMP_ID(I)
                  J=J+1
               ENDIF
            ENDDO
            TEMP=WLOAD
            TEMP_ID=IDWLOAD
         ENDIF
        IF(BDC_MD)THEN
           BUF_SIZE=K10
           IF (KEEP(201).NE.0) THEN
              A=INT((DBLE(KEEP(100))/DBLE(2))/DBLE(NELIM))
              IF(K50.EQ.0)THEN
                 BUF_SIZE=MIN(BUF_SIZE,DBLE(INT(A*DBLE(NCB))))               
              ELSE
                 BUF_SIZE=MIN(BUF_SIZE,DBLE(INT(A*A)))
              ENDIF
           ENDIF
           BUF_SIZE=K10
           DO I=1,NUMBER_OF_PROCS
              A=DBLE(INT(MD_MEM(TEMP_ID(I))/
     $             DBLE(NELIM)))
              A=A*DBLE(NFRONT)
              IF(K50.EQ.0)THEN
                 B=DBLE((INT(DBLE(NCB)/DBLE(NUMBER_OF_PROCS))+1)*
     $                DBLE(NFRONT))
              ELSE
                 WHAT = 5 
                 x=OTHERS
                 CALL SMUMPS_503(WHAT, KEEP,KEEP8, NCB,
     $                NFRONT, MIN(NCB,OTHERS), J, X)
                 B=DBLE(X)+(DBLE(J)*DBLE(NELIM))
              ENDIF
              NELIM_MEM_SIZE=A+B
              MEM_SIZE_WEAK(I)=NELIM_MEM_SIZE
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(I)=
     $                 TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                 LU_USAGE(TEMP_ID(I))-NIV2(TEMP_ID(I)+1)
               ELSE
                  MEM_SIZE_STRONG(I)=
     $                 TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                 LU_USAGE(TEMP_ID(I))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(I)=
     $                    TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                    LU_USAGE(TEMP_ID(I))-NIV2(TEMP_ID(I)+1)-
     $                    (SBTR_MEM(TEMP_ID(I))-SBTR_CUR(TEMP_ID(I)))
                  ELSE
                     MEM_SIZE_STRONG(I)=
     $                    TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                    LU_USAGE(TEMP_ID(I))-
     $                    (SBTR_MEM(TEMP_ID(I))-SBTR_CUR(TEMP_ID(I)))
                  ENDIF
               ENDIF
            ENDIF
            IF(MIN(MEM_SIZE_STRONG(I),MEM_SIZE_WEAK(I)).LT.0)THEN
                IF(MEM_SIZE_STRONG(I).LT.0)THEN
                   MEM_SIZE_STRONG(I)=0
                ELSE
                   MEM_SIZE_WEAK(I)=0                   
                ENDIF
             ENDIF
          ENDDO
       ELSE
          BUF_SIZE=K10
          DO I=1,NUMBER_OF_PROCS
            IF((SBTR_WHICH_M.EQ.0).OR.(.NOT.BDC_SBTR))THEN
               IF(BDC_M2_MEM)THEN
                  MEM_SIZE_STRONG(I)=
     $                 TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                 LU_USAGE(TEMP_ID(I))-NIV2(TEMP_ID(I)+1)
               ELSE
                  MEM_SIZE_STRONG(I)=
     $                 TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                 LU_USAGE(TEMP_ID(I))
               ENDIF
            ELSE
               IF(BDC_SBTR)THEN
                  IF(BDC_M2_MEM)THEN
                     MEM_SIZE_STRONG(I)=
     $                    TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                    LU_USAGE(TEMP_ID(I))-NIV2(TEMP_ID(I)+1)-
     $                    (SBTR_MEM(TEMP_ID(I))-SBTR_CUR(TEMP_ID(I)))
                  ELSE
                     MEM_SIZE_STRONG(I)=
     $                    TAB_MAXS(TEMP_ID(I))-DM_MEM(TEMP_ID(I))-
     $                    LU_USAGE(TEMP_ID(I))-
     $                    (SBTR_MEM(TEMP_ID(I))-SBTR_CUR(TEMP_ID(I)))                     
                  ENDIF
               ENDIF
            ENDIF
            MEM_SIZE_STRONG(I)=MAX(DBLE(0),MEM_SIZE_STRONG(I))
            MEM_SIZE_WEAK(I)=HUGE(MEM_SIZE_WEAK(I))
          ENDDO
       ENDIF
       IF((((NUMBER_OF_PROCS.LE.K83).AND.force_cand).AND.
     $      (TOTAL_COST.GE.SOMME)).OR.
     $      (.NOT.force_cand).OR.
     $      (((NUMBER_OF_PROCS+1).GT.K83).AND.force_cand))THEN
               REF=NSLAVES_REF
               SMALL_SET=NSLAVES_REF
               IF(.NOT.SMP)THEN
                  DO I=NSLAVES_REF,1,-1
                     SOMME=0
                     DO J=1,I
                        SOMME=SOMME+TEMP(J)
                     ENDDO
                     SOMME=(DBLE(I)*TEMP(I))-SOMME
                     IF(TOTAL_COST.GE.SOMME) GOTO 444
                  ENDDO
 444              CONTINUE
                  REF=I
                  SMALL_SET=REF
                  MAX_LOAD=TEMP(SMALL_SET)
               ELSE
                  X=MIN(SIZE_MY_SMP,NSLAVES_REF)
 450              CONTINUE
                  SOMME=0
                  DO J=1,X
                     SOMME=SOMME+(TEMP(X)-TEMP(J))
                  ENDDO
                  IF(SOMME.GT.TOTAL_COST)THEN
                     X=X-1
                     GOTO 450
                  ELSE
                     IF(X.LT.SIZE_MY_SMP) THEN
                        REF=X
                        SMALL_SET=REF
                        MAX_LOAD=TEMP(SMALL_SET)
                     ELSE
                        X=MIN(SIZE_MY_SMP,NSLAVES_REF)
                        J=X+1
                        MAX_LOAD=TEMP(X)
                        AFFECTED=MAX_LOAD
                        DO I=X+1,OTHERS
                           IF(TEMP(I).GT.MAX_LOAD)THEN
                              SOMME=SOMME+((I-1)*(TEMP(I)-MAX_LOAD))
                              AFFECTED=MAX_LOAD
                              MAX_LOAD=TEMP(I)
                           ELSE
                              SOMME=SOMME+(MAX_LOAD-TEMP(I))
                           ENDIF
                           IF(I.EQ.NSLAVES_REF)THEN
                              SMALL_SET=NSLAVES_REF
                              REF=SMALL_SET
                              GOTO 323
                           ENDIF
                           IF(SOMME.GT.TOTAL_COST)THEN
                              REF=I-1
                              SMALL_SET=I-1
                              MAX_LOAD=AFFECTED
                              GOTO 323
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
 323           CONTINUE
               MAX_LOAD=DBLE(0)
               DO I=1,SMALL_SET
                  MAX_LOAD=MAX(MAX_LOAD,TEMP(I))
               ENDDO
               TEMP_MAX_LOAD=MAX_LOAD
               NB_ROWS=0
               AFFECTED=0
               CHOSEN=0
               ACC=0
               NB_SAT=0
               NB_ZERO=0
               DO I=1,SMALL_SET
                  IF(K50.EQ.0)THEN
                     X=INT(BUF_SIZE/(NCB+1))-1
                     BANDE_K10=X*NFRONT
                  ELSE
                     A=1
                     B=DBLE(ACC+2)
                     C=DBLE(-BUF_SIZE+ACC+NELIM)
                     DELTA=(B*B)-(4*A*C)
                     x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                     IF(X.GT.NCB-ACC) X=NCB-ACC
                     BANDE_K10=x*(NELIM+ACC+x)
                  ENDIF
                  IF(HAVE_TYPE1_SON)THEN
                     IF(K50.EQ.0)THEN
                        X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                        BANDE_K10=X*NFRONT
                     ELSE
                        A=1
                        B=DBLE(ACC+2+NELIM)
                        C=DBLE(-BUF_SIZE+ACC+NELIM)
                        DELTA=(B*B)-(4*A*C)
                        x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                        IF(X.GT.NCB-ACC) X=NCB-ACC
                        BANDE_K10=x*(NELIM+ACC+x)
                     ENDIF
                  ENDIF
                  MAX_MEM_ALLOW=BANDE_K10
                  IF(BDC_MD)THEN
                     MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                    DBLE(BANDE_K10))
                     MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                  ENDIF
                  IF(K50.EQ.0)THEN
                     KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                     x=INT((MAX_LOAD-TEMP(I))/
     $                    (NELIM*(2*NFRONT-NELIM)))
                     IF(x.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           x=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF                        
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF
                  IF(K50.NE.0)THEN
                        A=1
                        B=ACC+NELIM
                        C=DBLE(-MAX_MEM_ALLOW)
                        DELTA=((B*B)-(4*A*C))
                        KMAX=INT((-B+SQRT(DELTA))/(2*A))
                     A=DBLE(NELIM)
                     B=DBLE(NELIM)*(DBLE(NELIM)+DBLE(2*ACC+1))
                     C=-(MAX_LOAD-TEMP(I))
                     DELTA=(B*B-(4*A*C))
                     x=INT((-B+SQRT(DELTA))/(2*A))
                     IF(X.LT.0) THEN
                        WRITE(*,*)MYID,
     *    ': Internal error 1 in SMUMPS_518'
                        CALL SMUMPS_ABORT()
                     ENDIF
                     IF(x.GE.KMAX)THEN
                        IF(KMAX.GE.KMIN)THEN
                           x=KMAX
                           NB_SAT=NB_SAT+1
                        ELSE
                           X=0
                        ENDIF
                     ELSE
                        IF(X.LT.KMIN)THEN
                           X=0
                        ENDIF
                     ENDIF
                     IF((ACC+X).GT.NCB) X=NCB-ACC
                  ENDIF                  
                  NB_ROWS(I)=X
                  ACC=ACC+X
                  CHOSEN=CHOSEN+1
                  IF(SMP)THEN
                     IF(MIN_LOAD.GT.TEMP(I))THEN
                        MIN_LOAD=TEMP(I)
                        POS_MIN_LOAD=I
                     ENDIF
                  ENDIF
                  AFFECTED=MAX_LOAD
                  IF(K50.EQ.0)THEN
                     MAX_LOAD=MAX(MAX_LOAD,
     $                    (TEMP(I)+(DBLE(NELIM) *
     $                    DBLE(NB_ROWS(I)))+
     $                    (DBLE(NB_ROWS(I))*DBLE(NELIM)*
     $                    DBLE(2*NFRONT-NELIM-1))))
                  ELSE
                     MAX_LOAD=MAX(MAX_LOAD,
     $               TEMP(I)+(DBLE(NELIM) * DBLE(NB_ROWS(I)))*
     $                    DBLE(2*(NELIM+ACC)-NB_ROWS(I)
     $                    -NELIM+1))
                  ENDIF
                  IF(AFFECTED.LT.MAX_LOAD)THEN
                  ENDIF
                  IF(NCB-ACC.LT.KMIN) GOTO 888
                  IF(NCB.EQ.ACC) GOTO 888
                  IF(ACC.GT.NCB) THEN
                    WRITE(*,*)MYID,
     *      ': Internal error 2 in SMUMPS_518'
                    CALL SMUMPS_ABORT()
                  ENDIF
               ENDDO
 888           CONTINUE
               SOMME=0
               X=NFRONT-NCB
               IF((ACC.GT.NCB))THEN
                  WRITE(*,*)MYID,
     *          ': Internal error 3 in SMUMPS_518'
                  CALL SMUMPS_ABORT()
               ENDIF
               IF((ACC.LT.NCB))THEN
                  IF(K50.NE.0)THEN
                     IF(SMALL_SET.LE.OTHERS)THEN
                       IF((NB_SAT.EQ.SMALL_SET).AND.(SMALL_SET.LT.
     $                      NSLAVES_REF))THEN
                          SMALL_SET=REF+1
                          REF=REF+1
                          NB_ROWS=0
                          GOTO 323
                       ENDIF
                       ADDITIONNAL_ROWS_SPECIAL=NCB-ACC
                       DO I=1,SMALL_SET
                          MAX_LOAD=TEMP_MAX_LOAD
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=DBLE(NELIM)*
     $                         DBLE(ADDITIONNAL_ROWS)*
     $                         (2*NFRONT-ADDITIONNAL_ROWS-NELIM+1)
                          SOMME=SOMME/(SMALL_SET-NB_SAT)
                          NB_ROWS=0
                          NB_ZERO=0
                          ACC=0
                          CHOSEN=0
                          NB_SAT=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO J=1,SMALL_SET
                             A=1
                             B=DBLE(ACC+2)
                             C=DBLE(-BUF_SIZE+ACC+NELIM)
                             DELTA=(B*B)-(4*A*C)
                             x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K10=x*(NELIM+ACC+x)
                             IF(HAVE_TYPE1_SON)THEN
                                A=1
                                B=DBLE(ACC+2+NELIM)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                        MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             A=1
                             B=ACC+NELIM
                             C=DBLE(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(4*A*C))
                             KMAX=INT((-B+SQRT(DELTA))/(2*A))
                             A=DBLE(NELIM)
                             B=(DBLE(NELIM)*DBLE(NELIM+2*ACC+1))
                             C=-(MAX_LOAD-TEMP(J)+SOMME)
                             DELTA=(B*B-(4*A*C))
                             x=INT((-B+SQRT(DELTA))/(2*A))
                             X=X+1
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     *    ': Internal error 4 in SMUMPS_518'
                                CALL SMUMPS_ABORT()
                             ENDIF
                             IF(x.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   x=KMAX
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.MIN(KMIN,KMAX))THEN
                                   NB_ZERO=NB_ZERO+1
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(J)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(J))THEN
                                   MIN_LOAD=TEMP(J)
                                   POS_MIN_LOAD=I
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             AFFECTED=MAX_LOAD
                             TEMP_MAX_LOAD=MAX(TEMP_MAX_LOAD,
     $                            TEMP(J)+(DBLE(NELIM) *
     $                            DBLE(NB_ROWS(J)))*
     $                            DBLE(2*(NELIM+
     $                            ACC)-NB_ROWS(J)
     $                            -NELIM+1))
                             IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                   IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                      SMALL_SET=REF+1
                                      REF=REF+1
                                      NB_ROWS=0
                                      GOTO 323
                                   ENDIF
                                ENDIF
                             ENDIF
                             IF(NCB.EQ.ACC) GOTO 666
                          ENDDO
                          IF(NB_SAT.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF(NB_ZERO.EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                          IF((NB_SAT+NB_ZERO).EQ.SMALL_SET)THEN
                             IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                SMALL_SET=REF+1
                                REF=REF+1
                                NB_ROWS=0
                                GOTO 323
                             ELSE
                                GOTO 434
                             ENDIF
                          ENDIF
                       ENDDO
 434                   CONTINUE
                       ADDITIONNAL_ROWS=NCB-ACC
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          IF(ADDITIONNAL_ROWS.LT.KMIN)THEN
                             I=CHOSEN
                             J=ACC
 436                         CONTINUE
                             IF(NB_ROWS(I).NE.0)THEN
                                J=J-NB_ROWS(I)
                                A=1
                                B=DBLE(J+2)
                                C=DBLE(-BUF_SIZE+J+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-J) X=NCB-J
                                BANDE_K10=x*(NELIM+J+x)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=1
                                   B=DBLE(J+2+NELIM)
                                   C=DBLE(-BUF_SIZE+J+NELIM)
                                   DELTA=(B*B)-(4*A*C)
                                   x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                   IF(X.GT.NCB-J) X=NCB-J
                                   BANDE_K10=x*(NELIM+J+x)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K10
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=MIN(
     $                         MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                                  DBLE(BANDE_K10))
                                   MAX_MEM_ALLOW=MAX(DBLE(0),
     $                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=1
                                B=J+NELIM
                                C=DBLE(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(4*A*C))
                                KMAX=INT((-B+SQRT(DELTA))/(2*A))
                                IF(NB_ROWS(I).NE.KMAX)THEN
                                   IF(NCB-J.LE.KMAX)THEN
                                      NB_ROWS(I)=+NCB-J
                                      ADDITIONNAL_ROWS=0
                                   ENDIF
                                ENDIF
                                TEMP_MAX_LOAD=MAX(TEMP_MAX_LOAD,
     $                               TEMP(I)+
     $                               (DBLE(NELIM) * DBLE(NB_ROWS(I)))*
     $                               DBLE(2*(NELIM+
     $                               ACC)-NB_ROWS(I)
     $                               -NELIM+1))
                                IF(REF.LE.NUMBER_OF_PROCS-1)THEN
                                   IF(TEMP_MAX_LOAD.GT.TEMP(REF+1))THEN
                                      IF(SMALL_SET.LT.NSLAVES_REF)THEN
                                         SMALL_SET=REF+1
                                         REF=REF+1
                                         NB_ROWS=0
                                         GOTO 323
                                      ENDIF
                                   ENDIF
                                ENDIF
                             ELSE
                                I=I-1
                                IF(I.NE.0)GOTO 436
                             ENDIF
                             IF(ADDITIONNAL_ROWS.NE.0)THEN
                                I=CHOSEN
                                IF(I.NE.SMALL_SET)THEN
                                   I=I+1
                                   IF(NB_ROWS(I).NE.0)THEN
                                      WRITE(*,*)MYID,
     *    ': Internal error 5 in SMUMPS_518'
                                      CALL SMUMPS_ABORT()
                                   ENDIF
                                ENDIF
                                NB_ROWS(I)=NB_ROWS(I)+ADDITIONNAL_ROWS
                                ADDITIONNAL_ROWS=0
                             ENDIF
                             CHOSEN=I
                          ENDIF
                       ENDIF
                       I=CHOSEN+1
                       DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                      .AND.(I.LE.NUMBER_OF_PROCS))
                          IF((TEMP(I).LE.MAX_LOAD))THEN
                             A=1
                             B=DBLE(ACC+2)
                             C=DBLE(-BUF_SIZE+ACC+NELIM)
                             DELTA=(B*B)-(4*A*C)
                             x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K10=x*(NELIM+ACC+x)
                             IF(HAVE_TYPE1_SON)THEN
                                A=1
                                B=DBLE(ACC+2+NELIM)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             A=1
                             B=ACC+NELIM
                             C=DBLE(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(4*A*C))
                             KMAX=INT((-B+SQRT(DELTA))/(2*A))
                             A=DBLE(NELIM)
                             B=DBLE(NELIM)*DBLE(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(I))
                             DELTA=(B*B-(4*A*C))
                             x=INT((-B+SQRT(DELTA))/(2*A))
                             IF(x.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   x=KMAX
                                ELSE
                                   x=0
                                ENDIF
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=0
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(I)=X
                             ACC=ACC+X
                             ADDITIONNAL_ROWS=NCB-ACC
                          ELSE IF((TEMP(I).GT.MAX_LOAD))THEN
                             MAX_LOAD=TEMP(I)
                             NB_SAT=0
                             ACC=0
                             NB_ROWS=0
                             DO J=1,I
                                A=1
                                B=DBLE(ACC+2)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=1
                                   B=DBLE(ACC+2+NELIM)
                                   C=DBLE(-BUF_SIZE+ACC+NELIM)
                                   DELTA=(B*B)-(4*A*C)
                                   x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K10=x*(NELIM+ACC+x)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K10
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                                  DBLE(BANDE_K10))
                                   MAX_MEM_ALLOW=MAX(DBLE(0),
     $                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=1
                                B=ACC+NELIM
                                C=DBLE(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(4*A*C))
                                KMAX=INT((-B+SQRT(DELTA))/(2*A))
                                A=DBLE(NELIM)
                                B=DBLE(NELIM)*DBLE(NELIM+2*ACC+1)
                                C=-(MAX_LOAD-TEMP(J))
                                DELTA=(B*B-(4*A*C))
                                x=INT((-B+SQRT(DELTA))/(2*A))
                                IF(X.LT.0) THEN
                        WRITE(*,*)MYID,
     *    ': Internal error 6 in SMUMPS_518'
                                   CALL SMUMPS_ABORT()
                                ENDIF
                                IF(x.GE.KMAX)THEN
                                   IF(KMAX.GE.KMIN)THEN
                                      x=KMAX
                                      NB_SAT=NB_SAT+1
                                   ELSE
                                      X=0
                                   ENDIF
                                ELSE
                                   IF(X.LT.MIN(KMIN,KMAX))THEN
                                      X=0
                                   ENDIF
                                ENDIF
                                IF((ACC+X).GT.NCB) X=NCB-ACC
                                NB_ROWS(J)=X
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(J))THEN
                                      MIN_LOAD=TEMP(J)
                                      POS_MIN_LOAD=I
                                   ENDIF
                                ENDIF
                                ACC=ACC+X
                                MAX_LOAD=MAX(MAX_LOAD,
     $                               TEMP(J)+
     $                               (DBLE(NELIM)*DBLE(NB_ROWS(J)))*
     $                               DBLE(2*(NELIM+
     $                               ACC)-NB_ROWS(J)
     $                               -NELIM+1))
                                IF(NCB.EQ.ACC) GOTO 741
                                IF(NCB-ACC.LT.KMIN) GOTO 210
                             ENDDO
 210                         CONTINUE
                          ENDIF
 741                      CONTINUE
                          I=I+1
                          ADDITIONNAL_ROWS=NCB-ACC
                       ENDDO
                       CHOSEN=I-1
                       IF(ADDITIONNAL_ROWS.NE.0)THEN
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=DBLE(NELIM)*DBLE(ADDITIONNAL_ROWS)*
     $                         DBLE(2*NFRONT-ADDITIONNAL_ROWS-NELIM+1)
                          SOMME=SOMME/NUMBER_OF_PROCS
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO I=1,OTHERS
                             A=1
                             B=DBLE(ACC+2)
                             C=DBLE(-BUF_SIZE+ACC+NELIM)
                             DELTA=(B*B)-(4*A*C)
                             x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K10=x*(NELIM+ACC+x)
                             IF(HAVE_TYPE1_SON)THEN
                                A=1
                                B=DBLE(ACC+2+NELIM)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             A=1
                             B=ACC+NELIM
                             C=DBLE(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(4*A*C))
                             KMAX=INT((-B+SQRT(DELTA))/(2*A))
                             A=DBLE(NELIM)
                             B=DBLE(NELIM)*DBLE(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(I)+SOMME)
                             DELTA=(B*B-(4*A*C))
                             x=INT((-B+SQRT(DELTA))/(2*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     *    ': Internal error 7 in SMUMPS_518'
                                CALL SMUMPS_ABORT()
                             ENDIF
                             IF(x.GE.KMAX)THEN
                                IF(KMAX.GE.KMIN)THEN
                                   x=KMAX
                                ELSE
                                   X=0
                                ENDIF
                             ELSE
                                IF(X.LT.MIN(KMIN,KMAX))THEN
                                   X=MIN(KMAX,KMIN)
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(I)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(I))THEN
                                   MIN_LOAD=TEMP(I)
                                   POS_MIN_LOAD=I
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 488
                          ENDDO
 488                      CONTINUE
                          ADDITIONNAL_ROWS=NCB-ACC
                          SOMME=DBLE(NELIM)*
     $                         DBLE(ADDITIONNAL_ROWS)*
     $                         DBLE(2*NFRONT-ADDITIONNAL_ROWS-NELIM+1)
                          SOMME=SOMME/NUMBER_OF_PROCS
                          NB_ROWS=0
                          ACC=0
                          CHOSEN=0
                          IF(SMP)THEN
                             MIN_LOAD=TEMP(1)
                             POS_MIN_LOAD=1
                          ENDIF
                          DO I=1,OTHERS
                             A=1
                             B=DBLE(ACC+2)
                             C=DBLE(-BUF_SIZE+ACC+NELIM)
                             DELTA=(B*B)-(4*A*C)
                             x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                             IF(X.GT.NCB-ACC) X=NCB-ACC
                             BANDE_K10=x*(NELIM+ACC+x)
                             IF(HAVE_TYPE1_SON)THEN
                                A=1
                                B=DBLE(ACC+2+NELIM)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                               MEM_SIZE_STRONG(I))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             A=1
                             B=ACC+NELIM
                             C=DBLE(-MAX_MEM_ALLOW)
                             DELTA=((B*B)-(4*A*C))
                             KMAX=INT((-B+SQRT(DELTA))/(2*A))
                             A=DBLE(NELIM)
                             B=DBLE(NELIM)*DBLE(NELIM+2*ACC+1)
                             C=-(MAX_LOAD-TEMP(I)+SOMME)
                             DELTA=(B*B-(4*A*C))
                             x=INT((-B+SQRT(DELTA))/(2*A))
                             IF(X.LT.0) THEN
                                WRITE(*,*)MYID,
     *    ': Internal error 8 in SMUMPS_518'
                                CALL SMUMPS_ABORT()
                             ENDIF
                             IF(x.GE.KMAX)THEN
                                x=KMAX
                             ELSE
                                IF(X.LT.KMIN)THEN
                                   X=KMIN
                                ENDIF
                             ENDIF
                             IF((ACC+X).GT.NCB) X=NCB-ACC
                             NB_ROWS(I)=X
                             IF(SMP)THEN
                                IF(MIN_LOAD.GT.TEMP(I))THEN
                                   MIN_LOAD=TEMP(I)
                                   POS_MIN_LOAD=I
                                ENDIF
                             ENDIF
                             CHOSEN=CHOSEN+1
                             ACC=ACC+X
                             IF(NCB.EQ.ACC) GOTO 666
                             IF(NCB-ACC.LT.KMIN) GOTO 477
                          ENDDO
 477                      CONTINUE
                          IF(ACC.NE.NCB)THEN
                             NB_SAT=0
                             ACC=0
                             CHOSEN=0
                             IF(SMP)THEN
                                MIN_LOAD=TEMP(1)
                                POS_MIN_LOAD=1
                             ENDIF
                             DO I=1,OTHERS
                                A=1
                                B=DBLE(ACC+2)
                                C=DBLE(-BUF_SIZE+ACC+NELIM)
                                DELTA=(B*B)-(4*A*C)
                                x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                IF(X.GT.NCB-ACC) X=NCB-ACC
                                BANDE_K10=x*(NELIM+ACC+x)
                                IF(HAVE_TYPE1_SON)THEN
                                   A=1
                                   B=DBLE(ACC+2+NELIM)
                                   C=DBLE(-BUF_SIZE+ACC+NELIM)
                                   DELTA=(B*B)-(4*A*C)
                                   x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                   IF(X.GT.NCB-ACC) X=NCB-ACC
                                   BANDE_K10=x*(NELIM+ACC+x)
                                ENDIF
                                MAX_MEM_ALLOW=BANDE_K10
                                IF(BDC_MD)THEN
                                   MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                                  MEM_SIZE_STRONG(I))
                                   MAX_MEM_ALLOW=MAX(DBLE(0),
     $                                  MAX_MEM_ALLOW)
                                ENDIF
                                A=1
                                B=ACC+NELIM
                                C=DBLE(-MAX_MEM_ALLOW)
                                DELTA=((B*B)-(4*A*C))
                                KMAX=INT((-B+SQRT(DELTA))/(2*A))
                                X=KMAX-NB_ROWS(I)
                                IF((ACC+NB_ROWS(I)+X).GT.NCB)
     $                               X=NCB-(ACC+NB_ROWS(I))
                                NB_ROWS(I)=NB_ROWS(I)+X
                                IF((NB_ROWS(I)*(NB_ROWS(I)+ACC)).EQ.
     $                               BANDE_K10)THEN
                                   NB_SAT=NB_SAT+1
                                ENDIF
                                ACC=ACC+NB_ROWS(I)
                                IF(SMP)THEN
                                   IF(MIN_LOAD.GT.TEMP(I))THEN
                                      MIN_LOAD=TEMP(I)
                                      POS_MIN_LOAD=I
                                   ENDIF
                                ENDIF
                                CHOSEN=CHOSEN+1
                                IF(NCB.EQ.ACC) GOTO 666
                                IF(NCB-ACC.LT.KMIN) GOTO 834
                             ENDDO
 834                         CONTINUE
                          ENDIF
                          IF(ACC.NE.NCB)THEN
                            ADDITIONNAL_ROWS=NCB-ACC
                            SOMME=DBLE(NELIM)*
     $                           DBLE(ADDITIONNAL_ROWS)*
     $                           DBLE(2*NFRONT-ADDITIONNAL_ROWS-
     $                           NELIM+1)
                            SOMME=SOMME/(NUMBER_OF_PROCS-NB_SAT)
                            ACC=0
                            DO I=1,CHOSEN
                               A=1
                               B=DBLE(ACC+2)
                               C=DBLE(-BUF_SIZE+ACC+NELIM)
                               DELTA=(B*B)-(4*A*C)
                               x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                               IF(X.GT.NCB-ACC) X=NCB-ACC
                               BANDE_K10=x*(NELIM+ACC+x)
                               IF(HAVE_TYPE1_SON)THEN
                                  A=1
                                  B=DBLE(ACC+2+NELIM)
                                  C=DBLE(-BUF_SIZE+ACC+NELIM)
                                  DELTA=(B*B)-(4*A*C)
                                  x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                                  IF(X.GT.NCB-ACC) X=NCB-ACC
                                  BANDE_K10=x*(NELIM+ACC+x)
                               ENDIF
                               IF((NB_ROWS(I)*(NB_ROWS(I)+ACC)).EQ.
     $                              BANDE_K10)THEN
                                  GOTO 102
                               ENDIF
                               A=DBLE(NELIM)
                               B=DBLE(NELIM)*
     $                              DBLE(NELIM+2*(ACC+NB_ROWS(I))+1)
                               C=-(SOMME)
                               DELTA=(B*B-(4*A*C))
                               x=INT((-B+SQRT(DELTA))/(2*A))
                               A=1
                               B=ACC+NELIM
                               C=DBLE(-BANDE_K10)
                               DELTA=((B*B)-(4*A*C))
                               KMAX=INT((-B+SQRT(DELTA))/(2*A))
                               IF(X.LT.0) THEN
                                  WRITE(*,*)MYID,
     *    ': Internal error 9 in SMUMPS_518'
                                  CALL SMUMPS_ABORT()
                               ENDIF
                               IF((ACC+X+NB_ROWS(I)).GT.NCB)THEN
                                  IF((NCB-ACC).GT.KMAX)THEN
                                     NB_ROWS(I)=KMAX
                                  ELSE
                                     NB_ROWS(I)=NCB-ACC
                                  ENDIF
                               ELSE
                                  IF((NB_ROWS(I)+X).GT.KMAX)THEN
                                     NB_ROWS(I)=KMAX
                                  ELSE
                                     NB_ROWS(I)=NB_ROWS(I)+X
                                  ENDIF
                               ENDIF
 102                           CONTINUE
                               ACC=ACC+NB_ROWS(I)
                               IF(NCB.EQ.ACC) THEN
                                  CHOSEN=I
                                  GOTO 666
                               ENDIF
                               IF(NCB-ACC.LT.KMIN) THEN
                                  CHOSEN=I
                                  GOTO 007
                               ENDIF
                            ENDDO
 007                        CONTINUE
                            DO I=1,CHOSEN
                               NB_ROWS(I)=NB_ROWS(I)+1
                               ACC=ACC+1
                               IF(ACC.EQ.NCB)GOTO 666
                            ENDDO
                            IF(ACC.LT.NCB)THEN
                               IF(SMP)THEN
                                  NB_ROWS(1)=NB_ROWS(1)+NCB-ACC
                               ELSE
                                  NB_ROWS(POS_MIN_LOAD)=
     $                                 NB_ROWS(POS_MIN_LOAD)+NCB-ACC
                               ENDIF
                            ENDIF
                         ENDIF
                         GOTO 666
                     ENDIF
                  ENDIF
                  GOTO 666
                 ENDIF
                 ADDITIONNAL_ROWS=NCB-ACC
                 I=CHOSEN+1
                 IF(NB_SAT.EQ.SMALL_SET) GOTO 777
                 DO I=1,SMALL_SET
                    IDWLOAD(I)=I
                    AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                    BANDE_K10=AFFECTED*NFRONT
                    IF(HAVE_TYPE1_SON)THEN
                       AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                       BANDE_K10=AFFECTED*NFRONT
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K10
                    IF(BDC_MD)THEN
                       MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                      DBLE(BANDE_K10))
                       MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                    ENDIF
                    WLOAD(I)=MAX_MEM_ALLOW
                 ENDDO
                 CALL SMUMPS_558(SMALL_SET, WLOAD, IDWLOAD)
                 NB_ZERO=0
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     $                (SMALL_SET.LT.NSLAVES_REF))THEN
                    SMALL_SET=REF+1
                    REF=REF+1
                    NB_ROWS=0
                    GOTO 323
                 ENDIF
                 IF((NB_SAT.EQ.SMALL_SET).AND.
     $                (SMALL_SET.LE.NUMBER_OF_PROCS))GOTO 777
                 AFFECTED=INT(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                 AFFECTED=MAX(INT(AFFECTED),1)
                 DO I=1,SMALL_SET
                    KMAX=INT(WLOAD(I)/DBLE(NFRONT))
                    IF(NB_ROWS(IDWLOAD(I)).EQ.KMAX)THEN 
                       GOTO 912
                    ENDIF
                    IF((NB_ROWS(IDWLOAD(I))+MIN(INT(AFFECTED),
     $                   ADDITIONNAL_ROWS)).GT.KMAX)THEN
                       IF(NB_ROWS(IDWLOAD(I)).GT.KMAX)THEN
                       ENDIF
                       ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                      (KMAX-NB_ROWS(IDWLOAD(I)))
                       NB_ROWS(IDWLOAD(I))=KMAX
                       NB_SAT=NB_SAT+1
                       IF(NB_SAT.EQ.SMALL_SET)THEN
                          IF(SMALL_SET.NE.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ELSE
                             MAX_LOAD=MAX(MAX_LOAD,
     $                            (TEMP(IDWLOAD(I))+(DBLE(NELIM) *
     $                            DBLE(NB_ROWS(IDWLOAD(I))))+
     $                            (DBLE(NB_ROWS(IDWLOAD(I)))*
     $                            DBLE(NELIM))*
     $                            DBLE(2*NFRONT-NELIM-1)))
                             GOTO 777
                          ENDIF
                       ENDIF
                       AFFECTED=INT(ADDITIONNAL_ROWS/(SMALL_SET-NB_SAT))
                       AFFECTED=MAX(INT(AFFECTED),1)
                    ELSE
                       IF((NB_ROWS(IDWLOAD(I))+MIN(INT(AFFECTED),
     $                      ADDITIONNAL_ROWS)).GE.KMIN)THEN
                          X=MIN(INT(AFFECTED),ADDITIONNAL_ROWS)
                          NB_ROWS(IDWLOAD(I))=NB_ROWS(IDWLOAD(I))+
     $                         X
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                       ELSE
                          x=INT((MAX_LOAD-TEMP(IDWLOAD(I)))/
     $                         (NELIM*(2*NFRONT-NELIM)))
                          IF(X+AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ELSE
                             X=AFFECTED+X
                          ENDIF
                          IF(X.GE.KMIN)THEN
                             NB_ROWS(IDWLOAD(I))=NB_ROWS(IDWLOAD(I))+
     $                            X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                            X
                          ELSE
                             NB_ZERO=NB_ZERO+1
                          ENDIF
                       ENDIF
                    ENDIF
 912                CONTINUE
                    MAX_LOAD=MAX(MAX_LOAD,
     $                   (TEMP(IDWLOAD(I))+(DBLE(NELIM)*
     $                   DBLE(NB_ROWS(IDWLOAD(I))))+
     $                   (DBLE(NB_ROWS(IDWLOAD(I)))*DBLE(NELIM))*
     $                   DBLE(2*NFRONT-NELIM-1)))
                    IF(SMALL_SET.LT.NUMBER_OF_PROCS)THEN
                       IF(MAX_LOAD.GT.TEMP(SMALL_SET+1))THEN
                          IF(SMALL_SET.LT.NSLAVES_REF)THEN
                             SMALL_SET=REF+1
                             REF=REF+1
                             NB_ROWS=0
                             GOTO 323
                          ENDIF
                       ENDIF
                    ENDIF
                    IF(SMALL_SET.EQ.NB_SAT)GOTO 777
                    IF(ADDITIONNAL_ROWS.EQ.0)THEN
                       CHOSEN=SMALL_SET
                       GOTO 049
                    ENDIF
                 ENDDO
 777             CONTINUE
                 IF((NB_ZERO.NE.0).AND.(ADDITIONNAL_ROWS.GE.KMIN))THEN
                    J=NB_ZERO
 732                CONTINUE
                    X=INT(ADDITIONNAL_ROWS/(J))
                    IF(X.LT.KMIN)THEN
                       J=J-1
                       GOTO 732
                    ENDIF
                    IF(X*J.LT.ADDITIONNAL_ROWS)THEN
                       X=X+1
                    ENDIF
                    DO I=1,SMALL_SET
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                         DBLE(BANDE_K10))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                       IF(NB_ROWS(I).EQ.0)THEN
                          IF(X.GT.ADDITIONNAL_ROWS)THEN
                             X=ADDITIONNAL_ROWS
                          ENDIF
                          IF(X.GT.KMAX)THEN
                             X=KMAX
                          ENDIF
                          IF(X.GT.KMIN)THEN
                             NB_ROWS(I)=X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             MAX_LOAD=MAX(MAX_LOAD,
     $                            (TEMP(I)+(DBLE(NELIM) *
     $                            DBLE(NB_ROWS(I)))+
     $                            (DBLE(NB_ROWS(I))*DBLE(NELIM))*
     $                            DBLE(2*NFRONT-NELIM-1)))
                          ENDIF
                       ENDIF
                    ENDDO
                 ENDIF
                 I=CHOSEN+1
                 DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                .AND.(I.LE.NUMBER_OF_PROCS))
                    IF((TEMP(I).LE.MAX_LOAD))THEN
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                         DBLE(BANDE_K10))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                       AFFECTED=INT((MAX_LOAD-TEMP(I))/
     $                      (NELIM*(2*NFRONT-NELIM)))
                       IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                          AFFECTED=ADDITIONNAL_ROWS
                       ENDIF
                       IF(NB_ROWS(I).LT.KMAX)THEN
                          IF((AFFECTED+NB_ROWS(I)).GT.KMAX)THEN
                             AFFECTED=KMAX-NB_ROWS(I)
                             NB_SAT=NB_SAT+1
                          ELSE
                             IF((AFFECTED+NB_ROWS(I)).LT.
     $                            KMIN)THEN
                                AFFECTED=0
                             ENDIF
                          ENDIF
                          NB_ROWS(I)=NB_ROWS(I)+AFFECTED
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                       ENDIF
                    ELSE IF((TEMP(I).GT.MAX_LOAD))THEN
                       IF(NB_SAT.EQ.I-1) GOTO 218
                       X=(ADDITIONNAL_ROWS/(I-1-NB_SAT))
                       ACC=1
                       DO J=1,I-1
                          AFFECTED=(REAL(NELIM * (NB_ROWS(J)+x))
     $                         +REAL((NB_ROWS(J)+x)*NELIM)*
     $                         REAL(2*NFRONT-NELIM-1))
                          IF((TEMP(J)+AFFECTED).GT.MAX_LOAD)THEN
                             ACC=0
                          ENDIF
                       ENDDO
                       IF(ACC.EQ.1)THEN
                          MAX_LOAD=TEMP(I)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                         .AND.(J.LT.I))
                             AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                             BANDE_K10=AFFECTED*NFRONT
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=INT((BUF_SIZE-NFRONT)/
     $                               (NFRONT+1))
                                BANDE_K10=AFFECTED*NFRONT
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     $                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO                          
                       ELSE
                          MAX_LOAD=TEMP(I)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                         .AND.(J.LT.I))
                             AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                             BANDE_K10=AFFECTED*NFRONT
                             IF(HAVE_TYPE1_SON)THEN
                                AFFECTED=INT((BUF_SIZE-NFRONT)/
     $                               (NFRONT+1))
                                BANDE_K10=AFFECTED*NFRONT
                             ENDIF
                             AFFECTED=((DBLE(NELIM)* DBLE(NB_ROWS(J)))
     $                            +(DBLE(NB_ROWS(J))*DBLE(NELIM))*
     $                            DBLE(2*NFRONT-NELIM-1))
                             x=INT((MAX_LOAD-(TEMP(J)+AFFECTED))/
     $                            (DBLE(NELIM)*DBLE(2*NFRONT-NELIM)))
                             IF(X.LT.0)THEN
                                WRITE(*,*)MYID,
     *    ': Internal error 10 in SMUMPS_518'
                                CALL SMUMPS_ABORT()
                             ENDIF
                             AFFECTED=X
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                             IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                                AFFECTED=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((AFFECTED+NB_ROWS(J)).GT.KMAX)THEN
                                   AFFECTED=KMAX-NB_ROWS(J)
                                   NB_SAT=NB_SAT+1
                                ELSE
                                   IF((AFFECTED+NB_ROWS(J)).LT.
     $                                  KMIN)THEN
                                      AFFECTED=0
                                   ENDIF
                                ENDIF
                                NB_ROWS(J)=NB_ROWS(J)+AFFECTED
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                               AFFECTED
                             ENDIF
                             J=J+1
                          ENDDO
                       ENDIF
                    ENDIF
 218                CONTINUE
                    I=I+1
                 ENDDO
                 CHOSEN=I-1
                 IF((CHOSEN.EQ.NUMBER_OF_PROCS-1).AND.
     $                 (ADDITIONNAL_ROWS.NE.0))THEN
                    DO I=1,CHOSEN
                       IF(NB_ROWS(I)+1.GE.KMIN)THEN
                          NB_ROWS(I)=NB_ROWS(I)+1
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-1
                       ENDIF
                       MAX_LOAD=MAX(MAX_LOAD,
     $                      (TEMP(I)+(DBLE(NELIM) *
     $                      DBLE(NB_ROWS(I)))+
     $                      (DBLE(NB_ROWS(I))*DBLE(NELIM))*
     $                      DBLE(2*NFRONT-NELIM-1)))
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 048
                    ENDDO
 048                CONTINUE
                 ENDIF
                 IF((ADDITIONNAL_ROWS.NE.0))THEN
                    IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                       I=CHOSEN+1
                    ELSE
                       IF(CHOSEN.NE.NUMBER_OF_PROCS)THEN
                          WRITE(*,*)MYID,
     *    ': Internal error 11 in SMUMPS_518'
                          CALL SMUMPS_ABORT()
                       ENDIF
                       I=CHOSEN
                    ENDIF
                    DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                   .AND.(I.LE.NUMBER_OF_PROCS))
                       IF(TEMP(I).LE.MAX_LOAD)THEN
                          AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=AFFECTED*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=AFFECTED*NFRONT
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K10
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                            DBLE(BANDE_K10))
                             MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                          AFFECTED=((DBLE(NELIM) * DBLE(NB_ROWS(I)))
     $                         +(DBLE(NB_ROWS(I))*DBLE(NELIM))*
     $                         DBLE(2*NFRONT-NELIM-1))
                          x=INT((MAX_LOAD-(TEMP(I)+AFFECTED))/
     $                         (DBLE(NELIM)*DBLE(2*NFRONT-NELIM)))
                          AFFECTED=X
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     *    ': Internal error 12 in SMUMPS_518'
                             CALL SMUMPS_ABORT()
                          ENDIF
                          IF(AFFECTED.GT.ADDITIONNAL_ROWS)THEN
                             AFFECTED=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(I).LT.KMAX)THEN
                             IF((AFFECTED+NB_ROWS(I)).GT.KMAX)THEN
                                AFFECTED=KMAX-NB_ROWS(I)
                             ELSE
                                IF((AFFECTED+NB_ROWS(I)).LT.
     $                               KMIN)THEN
                                   AFFECTED=0
                                ENDIF
                             ENDIF
                             NB_ROWS(I)=NB_ROWS(I)+AFFECTED
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-AFFECTED
                          ENDIF
                          IF(I.NE.NUMBER_OF_PROCS) GOTO 624
                       ELSE IF((TEMP(I).GT.MAX_LOAD))THEN
                          X=INT(ADDITIONNAL_ROWS/I-1)
                          X=MAX(X,1)
                          IF((MAX_LOAD+((DBLE(NELIM)*
     $                         DBLE(x))+(DBLE(
     $                         x)*DBLE(NELIM))*DBLE(
     $                         (2*NFRONT-NELIM-1)))).LE.TEMP(I))THEN
                             AFFECTED=X
                             POS=1
                          ELSE
                             POS=0
                          ENDIF
                          MAX_LOAD=TEMP(I)
                          J=1
                          DO WHILE ((ADDITIONNAL_ROWS.NE.0)
     *                         .AND.(J.LT.I))
                             X=INT(BUF_SIZE/(NCB+1))-1
                             BANDE_K10=X*NFRONT
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(HAVE_TYPE1_SON)THEN
                                X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                                BANDE_K10=X*NFRONT
                             ENDIF
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                             IF(POS.EQ.0)THEN
                                AFFECTED=((DBLE(NELIM) *
     $                               DBLE(NB_ROWS(J)))
     $                               +(DBLE(NB_ROWS(J))*DBLE(NELIM))*
     $                               DBLE(2*NFRONT-NELIM-1))
                                x=INT((TEMP(I)-(TEMP(J)+AFFECTED))/
     $                               (NELIM*(2*NFRONT-NELIM)))
                             ELSE
                                X=INT(AFFECTED)
                             ENDIF
                             IF(x.GT.ADDITIONNAL_ROWS)THEN
                                x=ADDITIONNAL_ROWS
                             ENDIF
                             IF(NB_ROWS(J).LT.KMAX)THEN
                                IF((x+NB_ROWS(J)).GT.KMAX)THEN
                                   x=KMAX-NB_ROWS(J)
                                ELSE
                                   IF((NB_ROWS(J)+x).LT.
     $                                  KMIN)THEN
                                     x=0
                                  ENDIF
                               ENDIF
                               NB_ROWS(J)=NB_ROWS(J)+x
                               ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-x
                            ENDIF
                            J=J+1
                         ENDDO
                       ENDIF
 624                   CONTINUE
                       I=I+1
                    ENDDO
                    CHOSEN=I-1
                    IF(ADDITIONNAL_ROWS.NE.0)THEN
                       ACC=0
                       DO I=1,CHOSEN
                          X=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=X*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=X*NFRONT
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K10
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                            DBLE(BANDE_K10))
                             MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                          ENDIF
                             KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                          AFFECTED=((DBLE(NELIM) * DBLE(NB_ROWS(I)))
     $                    +(DBLE(NB_ROWS(I))*DBLE(NELIM))*
     $                    DBLE(2*NFRONT-NELIM-1))
                          x=INT((MAX_LOAD-
     $                         (TEMP(I)+AFFECTED))/
     $                         (NELIM*(2*NFRONT-NELIM)))
                          IF(X.LT.0)THEN
                             WRITE(*,*)MYID,
     *    ': Internal error 13 in SMUMPS_518'
                             CALL SMUMPS_ABORT()
                          ENDIF
                          IF(x.GT.ADDITIONNAL_ROWS)THEN
                             x=ADDITIONNAL_ROWS
                          ENDIF
                          IF(NB_ROWS(I).LT.KMAX)THEN
                             IF((x+NB_ROWS(I)).GE.KMAX)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                               (KMAX-NB_ROWS(I))
                                NB_ROWS(I)=KMAX
                             ELSE
                                IF((X+NB_ROWS(I)).GE.
     $                               KMIN)THEN
                                   NB_ROWS(I)=NB_ROWS(I)+x
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                                   ACC=ACC+1
                                ELSE
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(CHOSEN.LT.NUMBER_OF_PROCS)THEN
                          CHOSEN=CHOSEN+1
                       ENDIF
                       IF(ACC.EQ.0)THEN
                          ACC=1
                       ENDIF
                       X=INT(ADDITIONNAL_ROWS/ACC)
                       X=MAX(X,1)
                       ACC=0
                       DO I=1,CHOSEN
                          J=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=J*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             J=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=J*NFRONT
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K10
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                            DBLE(BANDE_K10))
                             MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                          AFFECTED=((DBLE(NELIM) * DBLE(NB_ROWS(I)))
     $                         +(DBLE(NB_ROWS(I))*DBLE(NELIM))*
     $                         DBLE(2*NFRONT-NELIM-1))
                          J=INT((MAX_LOAD-
     $                         (TEMP(I)+AFFECTED))/
     $                         (NELIM*(2*NFRONT-NELIM)))
                          IF(NB_ROWS(I).LT.KMAX)THEN
                             IF((MIN(x,J)+NB_ROWS(I)).GE.KMAX)THEN
                                IF((KMAX-NB_ROWS(I)).GT.
     $                               ADDITIONNAL_ROWS)THEN
                                   NB_ROWS(I)=NB_ROWS(I)+
     $                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ELSE
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                                  (KMAX-NB_ROWS(I))
                                   NB_ROWS(I)=KMAX
                                ENDIF
                             ELSE
                                IF((MIN(x,J)+NB_ROWS(I)).GE.
     $                            KMIN)THEN
                                   NB_ROWS(I)=NB_ROWS(I)+MIN(x,J)
                                   ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                                  MIN(x,J)
                                   ACC=ACC+1
                                ENDIF
                             ENDIF
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       IF(ACC.GT.0)THEN
                          DO I=1,CHOSEN
                             X=INT(BUF_SIZE/(NCB+1))-1
                             BANDE_K10=X*NFRONT
                             IF(HAVE_TYPE1_SON)THEN
                                X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                                BANDE_K10=X*NFRONT
                             ENDIF
                             MAX_MEM_ALLOW=BANDE_K10
                             IF(BDC_MD)THEN
                                MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)),
     $                               DBLE(BANDE_K10))
                                MAX_MEM_ALLOW=MAX(DBLE(0),
     $                               MAX_MEM_ALLOW)
                             ENDIF
                             KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                             IF(KMAX-NB_ROWS(I).LT.
     $                            ADDITIONNAL_ROWS)THEN
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                               (KMAX-NB_ROWS(I))
                                NB_ROWS(I)=KMAX
                             ELSE
                                IF(NB_ROWS(I).EQ.0)THEN
                                   IF(MIN(KMIN,KMAX).LT.
     $                                  ADDITIONNAL_ROWS)THEN
                                      NB_ROWS(I)=MIN(KMIN,KMAX)
                                      ADDITIONNAL_ROWS=
     $                                     ADDITIONNAL_ROWS-
     $                                     MIN(KMIN,KMAX)
                                   ENDIF
                                ELSE
                                   NB_ROWS(I)=NB_ROWS(I)+
     $                                  ADDITIONNAL_ROWS
                                   ADDITIONNAL_ROWS=0
                                ENDIF
                             ENDIF
                             IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                          ENDDO
                       ENDIF
                       DO I=1,CHOSEN
                          IDWLOAD(I)=I
                          AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=AFFECTED*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=AFFECTED*NFRONT
                          ENDIF                          
                          WLOAD(I)=(BANDE_K10-(NB_ROWS(I)*NFRONT))
                       ENDDO
                       CALL SMUMPS_558(NUMBER_OF_PROCS, WLOAD,
     $                      IDWLOAD)
                       NB_SAT=0
                       DO I=1,CHOSEN
                          X=INT(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                          X=MAX(X,1)
                          AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=AFFECTED*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=AFFECTED*NFRONT
                          ENDIF
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                            MEM_SIZE_STRONG(I))
                             MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=INT(MAX_MEM_ALLOW/(NFRONT))
                          IF(NB_ROWS(IDWLOAD(I)).LT.KMAX)THEN
                             IF((NB_ROWS(IDWLOAD(I))+X).LT.KMAX)THEN
                                NB_ROWS(IDWLOAD(I))=
     $                               NB_ROWS(IDWLOAD(I))+X
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                             ELSE
                                ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                               (KMAX-NB_ROWS(IDWLOAD(I)))
                                NB_ROWS(IDWLOAD(I))=KMAX
                             ENDIF
                          ENDIF
                          IF(NB_ROWS(IDWLOAD(I)).EQ.KMAX)THEN
                             NB_SAT=NB_SAT+1
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0) GOTO 049
                       ENDDO
                       DO I=1,CHOSEN
                          X=INT(BUF_SIZE/(NCB+1))-1
                          BANDE_K10=X*NFRONT
                          IF(HAVE_TYPE1_SON)THEN
                             X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                             BANDE_K10=X*NFRONT
                          ENDIF
                          MAX_MEM_ALLOW=BANDE_K10
                          IF(BDC_MD)THEN
                             MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                            MEM_SIZE_STRONG(I))
                             MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                          ENDIF
                          KMAX=INT(MAX_MEM_ALLOW/(NFRONT))
                          IF(KMAX-NB_ROWS(I).LT.ADDITIONNAL_ROWS)THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                            (KMAX-NB_ROWS(I))
                             NB_ROWS(I)=KMAX
                          ELSE
                             NB_ROWS(I)=NB_ROWS(I)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       X=INT(ADDITIONNAL_ROWS/CHOSEN)
                       X=MAX(X,1)
                       DO I=1,CHOSEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(I)=NB_ROWS(I)+X
                          IF(ADDITIONNAL_ROWS.EQ.0)GOTO 049
                       ENDDO
                       NB_ROWS(1)=NB_ROWS(1)+ADDITIONNAL_ROWS
                    ENDIF
                 ENDIF
 049             CONTINUE
              ENDIF
 666          CONTINUE
              SOMME=0
              X=0
              POS=0
              DO I=1,CHOSEN
                 X=X+NB_ROWS(I)
                 SOMME=SOMME+ DBLE(NB_ROWS(I))
              ENDDO
              GOTO 890
           ELSE IF((KEEP(83).GE.NUMBER_OF_PROCS).AND.force_cand)THEN
              MAX_LOAD=DBLE(0)
              DO I=1,OTHERS
                 MAX_LOAD=MAX(MAX_LOAD,TEMP(I))
              ENDDO
              ACC=0
              CHOSEN=0
              X=1
              DO I=1,OTHERS
              ENDDO
              DO I=2,OTHERS 
                 IF(TEMP(I).EQ.TEMP(1))THEN
                    X=X+1
                 ELSE 
                    GOTO 329
                 ENDIF
              ENDDO
 329          CONTINUE
              AFFECTED=TOTAL_COST/DBLE(X)
              TEMP_MAX_LOAD=DBLE(0)
              DO I=1,OTHERS
                 IF(K50.EQ.0)THEN
                    X=INT(BUF_SIZE/(NCB+1))-1
                    BANDE_K10=X*NFRONT
                 ELSE
                    A=1
                    B=DBLE(ACC+2)
                    C=DBLE(-BUF_SIZE+ACC+NELIM)
                    DELTA=(B*B)-(4*A*C)
                    x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K10=x*(NELIM+ACC+x)
                 ENDIF
                 IF(HAVE_TYPE1_SON)THEN
                    IF(K50.EQ.0)THEN
                       X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                       BANDE_K10=X*NFRONT
                    ELSE
                       A=1
                       B=DBLE(ACC+2+NELIM)
                       C=DBLE(-BUF_SIZE+ACC+NELIM)
                       DELTA=(B*B)-(4*A*C)
                       x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K10=x*(NELIM+ACC+x)
                    ENDIF
                 ENDIF
                 MAX_MEM_ALLOW=BANDE_K10
                 IF(BDC_MD)THEN
                    MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                    MIN(MEM_SIZE_WEAK(I),MEM_SIZE_STRONG(I)))
                    MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                 ENDIF
                 IF(K50.EQ.0)THEN
                       KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                    IF(AFFECTED+TEMP(I).GT.MAX_LOAD)THEN
                       SOMME=MAX_LOAD-TEMP(I)
                    ELSE
                       SOMME=AFFECTED
                    ENDIF
                    x=INT(SOMME/
     $                   (NELIM*(2*NFRONT-NELIM)))
                    IF(x.GT.KMAX)THEN
                       x=KMAX
                    ELSE
                       IF(X.LT.KMIN)THEN
                          X=MIN(KMIN,KMAX)
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 IF(K50.NE.0)THEN
                       A=1
                       B=ACC+NELIM
                       C=DBLE(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(4*A*C))
                       KMAX=INT((-B+SQRT(DELTA))/(2*A))
                    A=DBLE(NELIM)
                    B=DBLE(NELIM)*DBLE(NELIM+2*ACC+1)
                    IF(AFFECTED+TEMP(I).GT.MAX_LOAD)THEN
                       C=-(MAX_LOAD-TEMP(I))
                    ELSE
                       C=-AFFECTED
                    ENDIF
                    DELTA=(B*B-(4*A*C))
                    x=INT((-B+SQRT(DELTA))/(2*A))
                    IF(X.LT.0) THEN
                       WRITE(*,*)MYID,
     *    ': Internal error 14 in SMUMPS_518'
                       CALL SMUMPS_ABORT()
                    ENDIF
                    IF(x.GE.KMAX)THEN
                       IF(KMAX.GT.KMIN)THEN
                          x=KMAX
                       ELSE
                          x=0
                       ENDIF
                    ELSE
                       IF(X.LE.MIN(KMIN,KMAX))THEN
                          IF(KMAX.LT.KMIN)THEN
                             X=0
                          ELSE
                             X=MIN(KMIN,KMAX)
                          ENDIF
                       ENDIF
                    ENDIF
                    IF((ACC+X).GT.NCB) X=NCB-ACC
                 ENDIF
                 TEMP_MAX_LOAD=MAX(TEMP_MAX_LOAD,TEMP(I))
                 NB_ROWS(I)=X
                 CHOSEN=CHOSEN+1
                 ACC=ACC+X
                 IF(ACC.EQ.NCB) GOTO 541
              ENDDO
 541          CONTINUE
              IF(ACC.LT.NCB)THEN
                 IF(K50.EQ.0)THEN
                    ADDITIONNAL_ROWS=NCB-ACC
                    DO J=1,CHOSEN
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(
     $                    MIN(MEM_SIZE_WEAK(J),MEM_SIZE_STRONG(J)),
     $                         DBLE(BANDE_K10))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                       IF((NB_ROWS(J)).LT.KMAX)THEN
                          IF(ADDITIONNAL_ROWS.GT.(KMAX-NB_ROWS(J)))THEN
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                            (KMAX-NB_ROWS(J))
                             NB_ROWS(J)=KMAX
                          ELSE
                             NB_ROWS(J)=NB_ROWS(J)+ADDITIONNAL_ROWS
                             ADDITIONNAL_ROWS=0
                          ENDIF
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    X=INT(ADDITIONNAL_ROWS/CHOSEN)
                    X=MAX(X,1)
                    DO J=1,CHOSEN
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                         MEM_SIZE_STRONG(J))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=INT(MAX_MEM_ALLOW/NFRONT)
                       IF((NB_ROWS(J)+X).GT.KMAX)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                         (KMAX-NB_ROWS(J))
                          NB_ROWS(J)=KMAX
                       ELSE
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          NB_ROWS(J)=NB_ROWS(J)+X
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO 
                    DO I=1,CHOSEN
                       X=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=X*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          X=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=X*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                         MEM_SIZE_STRONG(I))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       KMAX=INT(MAX_MEM_ALLOW/(NFRONT))
                       IF(KMAX-NB_ROWS(I).LT.ADDITIONNAL_ROWS)THEN
                          ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                         (KMAX-NB_ROWS(I))
                          NB_ROWS(I)=KMAX
                       ELSE
                          NB_ROWS(I)=NB_ROWS(I)+ADDITIONNAL_ROWS
                          ADDITIONNAL_ROWS=0
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0)GOTO 889
                    ENDDO
                    DO I=1,NUMBER_OF_PROCS
                       IDWLOAD(I)=I
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF                          
                       WLOAD(I)=(BANDE_K10-(NB_ROWS(I)*NFRONT))
                    ENDDO
                    CALL SMUMPS_558(NUMBER_OF_PROCS, WLOAD,
     $                   IDWLOAD)
                    NB_SAT=0
                    DO I=1,CHOSEN
                       X=INT(ADDITIONNAL_ROWS/(CHOSEN-NB_SAT))
                       X=MAX(X,1)
                       AFFECTED=INT(BUF_SIZE/(NCB+1))-1
                       BANDE_K10=AFFECTED*NFRONT
                       IF(HAVE_TYPE1_SON)THEN
                          AFFECTED=INT((BUF_SIZE-NFRONT)/(NFRONT+1))
                          BANDE_K10=AFFECTED*NFRONT
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       KMAX=INT(MAX_MEM_ALLOW/(NFRONT))
                       IF(NB_ROWS(IDWLOAD(I)).LT.KMAX)THEN
                          IF((NB_ROWS(IDWLOAD(I))+X).LT.KMAX)THEN
                             NB_ROWS(IDWLOAD(I))=
     $                            NB_ROWS(IDWLOAD(I))+X
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                          ELSE
                             ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-
     $                            (KMAX-NB_ROWS(IDWLOAD(I)))
                             NB_ROWS(IDWLOAD(I))=KMAX
                          ENDIF
                       ENDIF
                       IF(NB_ROWS(IDWLOAD(I)).EQ.KMAX)THEN
                          NB_SAT=NB_SAT+1
                       ENDIF
                       IF(ADDITIONNAL_ROWS.EQ.0) GOTO 889
                    ENDDO
                    GOTO 994
                 ELSE
                    ACC=0
                    CHOSEN=0
                    DO I=1,OTHERS
                       A=1
                       B=DBLE(ACC+2)
                       C=DBLE(-BUF_SIZE+ACC+NELIM)
                       DELTA=(B*B)-(4*A*C)
                       x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K10=x*(NELIM+ACC+x)
                       IF(HAVE_TYPE1_SON)THEN
                          A=1
                          B=DBLE(ACC+2+NELIM)
                          C=DBLE(-BUF_SIZE+ACC+NELIM)
                          DELTA=(B*B)-(4*A*C)
                          x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                          IF(X.GT.NCB-ACC) X=NCB-ACC
                          BANDE_K10=x*(NELIM+ACC+x)
                       ENDIF
                       MAX_MEM_ALLOW=BANDE_K10
                       IF(BDC_MD)THEN
                          MAX_MEM_ALLOW=MIN(DBLE(BANDE_K10),
     $                         MEM_SIZE_STRONG(I))
                          MAX_MEM_ALLOW=MAX(DBLE(0),MAX_MEM_ALLOW)
                       ENDIF
                       A=1
                       B=ACC+NELIM
                       C=DBLE(-MAX_MEM_ALLOW)
                       DELTA=((B*B)-(4*A*C))
                       KMAX=INT((-B+SQRT(DELTA))/(2*A))
                       X=KMAX-NB_ROWS(I)
                       IF((ACC+NB_ROWS(I)+X).GT.NCB)
     $                            X=NCB-(ACC+NB_ROWS(I))
                       NB_ROWS(I)=NB_ROWS(I)+X
                       ACC=ACC+NB_ROWS(I)
                       CHOSEN=CHOSEN+1
                       IF(NCB.EQ.ACC) GOTO 889
                    ENDDO
                    ADDITIONNAL_ROWS=NCB-ACC
                 ENDIF
 993             CONTINUE
                 ACC=0
                 CHOSEN=0
                 DO I=1,OTHERS
                    A=1
                    B=DBLE(ACC+2)
                    C=DBLE(-BUF_SIZE+ACC+NELIM)
                    DELTA=(B*B)-(4*A*C)
                    x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                    IF(X.GT.NCB-ACC) X=NCB-ACC
                    BANDE_K10=x*(NELIM+ACC+x)
                    IF(HAVE_TYPE1_SON)THEN
                       A=1
                       B=DBLE(ACC+2+NELIM)
                       C=DBLE(-BUF_SIZE+ACC+NELIM)
                       DELTA=(B*B)-(4*A*C)
                       x=INT((-B+SQRT(DELTA))/DBLE(2*A))
                       IF(X.GT.NCB-ACC) X=NCB-ACC
                       BANDE_K10=x*(NELIM+ACC+x)
                    ENDIF
                    MAX_MEM_ALLOW=BANDE_K10
                    A=1
                    B=ACC+NELIM
                    C=DBLE(-MAX_MEM_ALLOW)
                    DELTA=((B*B)-(4*A*C))
                    KMAX=INT((-B+SQRT(DELTA))/(2*A))
                    X=KMAX-NB_ROWS(I)
                    IF((ACC+NB_ROWS(I)+X).GT.NCB)
     $                   X=NCB-(ACC+NB_ROWS(I))
                    NB_ROWS(I)=NB_ROWS(I)+X
                    ACC=ACC+NB_ROWS(I)
                    CHOSEN=CHOSEN+1
                    IF(NCB.EQ.ACC) GOTO 889
                 ENDDO
                 ADDITIONNAL_ROWS=NCB-ACC
 994             CONTINUE
                 X=INT(DBLE(ADDITIONNAL_ROWS)/DBLE(OTHERS))
                 IF((X*OTHERS).LT.ADDITIONNAL_ROWS)THEN
                    X=X+1
                 ENDIF
                 DO I=1,OTHERS
                    NB_ROWS(I)=NB_ROWS(I)+X
                    ADDITIONNAL_ROWS=ADDITIONNAL_ROWS-X
                    IF(ADDITIONNAL_ROWS.LT.X)X=ADDITIONNAL_ROWS
                 ENDDO
                 CHOSEN=OTHERS
              ENDIF
           ENDIF
 889       CONTINUE
           MAX_LOAD=TEMP_MAX_LOAD
 890       CONTINUE
           J=CHOSEN
           X=0
              DO I=J,1,-1
                 IF(NB_ROWS(I).EQ.0)THEN
                    CHOSEN=CHOSEN-1
                    ELSE 
                       IF(NB_ROWS(I).GT.0)THEN
                          X=1
                       ELSE
                          WRITE(*,*)MYID,
     *    ': Internal error 15 in SMUMPS_518'
                          CALL SMUMPS_ABORT()
                       ENDIF
                    ENDIF
                 ENDDO
           NSLAVES_NODE=CHOSEN
           TAB_POS(NSLAVES_NODE+1)= NCB+1
           TAB_POS(SLAVEF+2) = CHOSEN
           POS=1
           x=1
           DO I=1,J
              IF(NB_ROWS(I).NE.0)THEN
                 SLAVES_LIST(X)=TEMP_ID(I)
                 TAB_POS(X)=POS
                 POS=POS+NB_ROWS(I) 
                 IF(NB_ROWS(I).LE.0)THEN
                    WRITE(*,*)MYID,
     *    ': Internal error 16 in SMUMPS_518'
                    CALL SMUMPS_ABORT()
                 ENDIF
                 X=X+1
               ENDIF
           ENDDO
           DO I=CHOSEN+1,NUMBER_OF_PROCS
              SLAVES_LIST(I)=TEMP_ID(I)
           ENDDO
           IF(POS.NE.(NCB+1))THEN
              WRITE(*,*)MYID,
     *    ': Internal error 17 in SMUMPS_518',
     $             POS,NCB+1
             CALL SMUMPS_ABORT()
           ENDIF
      END SUBROUTINE SMUMPS_518
      SUBROUTINE SMUMPS_520
     &      (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     *       STEP,POOL,LPOOL,PROCNODE,N)
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER STEP(KEEP(28)), POOL(LPOOL), PROCNODE(KEEP(28))
      LOGICAL UPPER
      INTEGER MEM_COST,J
      INTEGER NBINSUBTREE,I,NBTOP
      EXTERNAL SMUMPS_167, SMUMPS_283, SMUMPS_508,
     * SMUMPS_170
      LOGICAL SMUMPS_167, SMUMPS_283, SMUMPS_508,
     * SMUMPS_170
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF(KEEP(47).LT.2)THEN
         WRITE(*,*)'SMUMPS_520 must
     $        be called with K47>=2'
         CALL SMUMPS_ABORT()
      ENDIF        
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
      MEM_COST=SMUMPS_543(INODE)
         IF((DM_MEM(MYID)+DBLE(MEM_COST)+ PEAK_SBTR_CUR_LOCAL-
     *        SBTR_CUR_LOCAL)
     *        .GT.MAX_PEAK_STK)THEN
            DO I=NBTOP-1,1,-1
               INODE = POOL( LPOOL - 2 - I)
               MEM_COST=SMUMPS_543(INODE)
               IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
                  DO J=I+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
               IF((DM_MEM(MYID)+DBLE(MEM_COST)+PEAK_SBTR_CUR_LOCAL-
     *              SBTR_CUR_LOCAL).LE.
     *              MAX_PEAK_STK) THEN
                  DO J=I+1,NBTOP,-1
                     POOL(J-1)=POOL(J)
                  ENDDO
                  UPPER=.TRUE.
                  RETURN
               ENDIF
            ENDDO
            IF(NBINSUBTREE.NE.0)THEN
               INODE = POOL( NBINSUBTREE )
               IF(.NOT.SMUMPS_170(STEP(INODE),
     *              PROCNODE,SLAVEF))THEN
                  WRITE(*,*)
     *        'Internal error 1 in SMUMPS_520'
                  CALL SMUMPS_ABORT()
               ENDIF
               UPPER=.FALSE.
               RETURN
            ENDIF
            INODE=POOL(LPOOL-2-NBTOP)
            UPPER=.TRUE.
            RETURN
         ENDIF
      ENDIF
      UPPER=.TRUE.
      END SUBROUTINE SMUMPS_520
      SUBROUTINE SMUMPS_513(WHAT)
      IMPLICIT NONE
      LOGICAL WHAT
      IF(.NOT.BDC_POOL_MNG)THEN
         WRITE(*,*)'SMUMPS_513
     $        should be called when K81>0 and K47>2'
      ENDIF
      IF(WHAT)THEN
         PEAK_SBTR_CUR_LOCAL=PEAK_SBTR_CUR_LOCAL+
     $        DBLE(MEM_SUBTREE(INDICE_SBTR))
         IF(.NOT.BDC_SBTR) INDICE_SBTR=INDICE_SBTR+1
      ELSE
         PEAK_SBTR_CUR_LOCAL=0
         SBTR_CUR_LOCAL=0
      ENDIF
      END SUBROUTINE SMUMPS_513
      INTEGER FUNCTION SMUMPS_543( INODE )
      IMPLICIT NONE
      INTEGER INODE,LEVEL,I,NELIM,NFR
      INTEGER COST
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      I = INODE
      NELIM = 0
 10   CONTINUE
      IF ( I > 0 ) THEN
        NELIM = NELIM + 1
        I = FILS_LOAD(I)
        GOTO 10
      ENDIF
      NFR = ND_LOAD( STEP_LOAD(INODE) )
      LEVEL = SMUMPS_330( STEP_LOAD(INODE), PROCNODE_LOAD, NPROCS )
      IF (LEVEL .EQ. 1) THEN
        COST =  NFR * NFR 
      ELSE
        IF ( K50 == 0 ) THEN
           COST =  NFR * NELIM 
        ELSE
           COST = NELIM * NELIM 
        ENDIF
      ENDIF
      SMUMPS_543=COST
      RETURN
      END FUNCTION SMUMPS_543
      RECURSIVE SUBROUTINE SMUMPS_515(FLAG,COST,COMM)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER COMM,WHAT,IERR
      LOGICAL FLAG
      DOUBLE PRECISION COST
      DOUBLE PRECISION TO_BE_SENT
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      IF(FLAG)THEN
         WHAT=17
         IF(BDC_M2_FLOPS)THEN
#if ! defined(OLD_LOAD_MECHANISM)
            TO_BE_SENT=DELTA_LOAD-COST
            DELTA_LOAD=DBLE(0)
#else
            TO_BE_SENT=LAST_LOAD_SENT-COST
            LAST_LOAD_SENT=LAST_LOAD_SENT-COST
#endif
         ELSE IF(BDC_M2_MEM)THEN
            IF(BDC_POOL.AND.(.NOT.BDC_MD))THEN
               TO_BE_SENT=MAX(TMP_M2,POOL_LAST_COST_SENT)
               POOL_LAST_COST_SENT=TO_BE_SENT
            ELSE IF(BDC_MD)THEN
#if ! defined(OLD_LOAD_MECHANISM)
               DELTA_MEM=DELTA_MEM+TMP_M2
               TO_BE_SENT=DELTA_MEM
#else
               TO_BE_SENT=DM_LAST_MEM_SENT+TMP_M2
               DM_LAST_MEM_SENT=DM_LAST_MEM_SENT+TMP_M2
#endif
            ELSE
               TO_BE_SENT=DBLE(0)
            ENDIF
         ENDIF
      ELSE
         WHAT=6
         TO_BE_SENT=DBLE(0)
      ENDIF
 111  CONTINUE
      CALL SMUMPS_460( WHAT,
     *         COMM, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *               FUTURE_NIV2,
#endif
     *         COST,TO_BE_SENT,MYID, iERR  )
      IF ( IERR == -1 )THEN
         CALL SMUMPS_467(COMM_LD, KEEP_LOAD)
         GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error in SMUMPS_500",
     *   IERR
         CALL SMUMPS_ABORT()
      ENDIF
      END SUBROUTINE SMUMPS_515
      SUBROUTINE SMUMPS_512(INODE,STEP,NSTEPS,PROCNODE,FRERE,
     *     NE,COMM,SLAVEF,MYID,KEEP,KEEP8,N)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER INODE,NSTEPS,MYID,SLAVEF,COMM,N
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER FRERE(NSTEPS),NE(NSTEPS),STEP(N),PROCNODE(NSTEPS)
      EXTERNAL SMUMPS_170,SMUMPS_275
      LOGICAL SMUMPS_170
      INTEGER I,NCB,NELIM
      INTEGER SMUMPS_275
      INTEGER FATHER_NODE,FATHER,WHAT,IERR
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      DOUBLE PRECISION ID
      IF((.NOT.BDC_M2_MEM).AND.(.NOT.BDC_M2_FLOPS))THEN
         WRITE(*,*)MYID,': Problem in SMUMPS_512'
         CALL SMUMPS_ABORT()
      ENDIF
      IF((INODE.LT.0).OR.(INODE.GT.N)) THEN
         RETURN
      ENDIF
      I=INODE
      NELIM = 0
 10   CONTINUE
      IF ( I > 0 ) THEN
         NELIM = NELIM + 1
         I = FILS_LOAD(I)
         GOTO 10
      ENDIF
      NCB=ND_LOAD(STEP_LOAD(INODE))-NELIM
      WHAT=5
      FATHER_NODE=DAD_LOAD(STEP_LOAD(INODE))
      IF (FATHER_NODE.EQ.0) THEN
         RETURN
      ENDIF
      IF((FRERE(STEP(FATHER_NODE)).EQ.0).AND.
     $     (FATHER_NODE.EQ.KEEP(38)))THEN
         RETURN 
      ENDIF
      IF(SMUMPS_170(STEP(FATHER_NODE),
     *   PROCNODE,SLAVEF)) THEN
         RETURN
      ENDIF
      FATHER=SMUMPS_275(STEP(FATHER_NODE),PROCNODE,SLAVEF)
      IF(FATHER.EQ.MYID)THEN
        IF(BDC_M2_MEM)THEN
           CALL PROCESS_NIV2_MEM_MSG(FATHER_NODE)
        ELSEIF(BDC_M2_FLOPS)THEN
           CALL PROCESS_NIV2_FLOPS_MSG(FATHER_NODE)
        ENDIF
        IF((KEEP(81).EQ.2).OR.(KEEP(81).EQ.3))THEN
           IF(SMUMPS_330(STEP_LOAD(INODE),
     $          PROCNODE_LOAD,NPROCS).EQ.1)THEN
              CB_COST_ID(POS_ID)=INODE
              CB_COST_ID(POS_ID+1)=1
              CB_COST_ID(POS_ID+2)=POS_MEM
              POS_ID=POS_ID+3
              CB_COST_MEM(POS_MEM)=MYID
              POS_MEM=POS_MEM+1
              CB_COST_MEM(POS_MEM)=NCB*NCB
              POS_MEM=POS_MEM+1
           ENDIF
        ENDIF
        GOTO 666
      ENDIF
 111  CONTINUE      
      CALL SMUMPS_519(WHAT, COMM, NPROCS,
     *     FATHER_NODE,INODE,NCB, KEEP(81),MYID,
     $     FATHER, IERR)
      IF (IERR == -1 ) THEN
        CALL SMUMPS_467(COMM, KEEP)
        GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
        WRITE(*,*) "Internal Error in SMUMPS_461",
     *  IERR
        CALL SMUMPS_ABORT()
      ENDIF 
 666  CONTINUE
      END SUBROUTINE SMUMPS_512
      SUBROUTINE SMUMPS_514(INODE,NUM_CALL)
      IMPLICIT NONE
      INTEGER I,J,MAXI,IND_MAXI
      INTEGER INODE,NUM_CALL
      INTEGER MAX_ALL_M2
      IF(BDC_M2_MEM)THEN
         IF(((NUM_CALL.EQ.1).AND.(BDC_MD)).OR.
     $       ((NUM_CALL.EQ.2).AND.(.NOT.BDC_MD)))THEN
            RETURN
         ENDIF
      ENDIF
      IF(POOL_SIZE.EQ.0)THEN
         IF(FRERE_LOAD(STEP_LOAD(INODE)).EQ.0) THEN
            RETURN
         ENDIF
         NB_SON(STEP_LOAD(INODE))=-1
         RETURN
      ENDIF
      DO I=POOL_SIZE,1,-1
         IF(POOL_NIV2(I).EQ.INODE) GOTO 666
      ENDDO
 666  CONTINUE
      IF(I.EQ.0)THEN
         NB_SON(STEP_LOAD(INODE))=-1
         RETURN         
      ENDIF
      IF(BDC_M2_MEM)THEN
         IF(POOL_NIV2_COST(I).EQ.MAX_M2)THEN
            TMP_M2=MAX_M2
            MAXI=0
            IND_MAXI=-9999
            DO J=POOL_SIZE,1,-1
               IF(J.NE.I) THEN
                  IF(POOL_NIV2_COST(J).GT.MAXI)THEN
                     MAXI=POOL_NIV2_COST(J)
                     IND_MAXI=J
                  ENDIF
               ENDIF
            ENDDO
            MAX_M2=MAXI
            J=IND_MAXI
            REMOVE_NODE_FLAG_MEM=.TRUE.
            REMOVE_NODE_COST_MEM=TMP_M2
            CALL SMUMPS_515(REMOVE_NODE_FLAG,MAX_M2,COMM_LD)
            NIV2(MYID+1)=MAX_M2
      ENDIF
      ELSEIF(BDC_M2_FLOPS)THEN
         IF(INODE.EQ.KEEP_LOAD(20))THEN
            RETURN
         ENDIF
         REMOVE_NODE_COST=POOL_NIV2_COST(I)
         REMOVE_NODE_FLAG=.TRUE.
         CALL SMUMPS_515(REMOVE_NODE_FLAG,
     $        -POOL_NIV2_COST(I),COMM_LD)
         NIV2(MYID+1)=NIV2(MYID+1)-POOL_NIV2_COST(I)
      ENDIF
      DO J=I+1,POOL_SIZE
         POOL_NIV2(J-1)=POOL_NIV2(J)
         POOL_NIV2_COST(J-1)=POOL_NIV2_COST(J)
      ENDDO
      POOL_SIZE=POOL_SIZE-1
      END SUBROUTINE SMUMPS_514
      RECURSIVE SUBROUTINE PROCESS_NIV2_MEM_MSG(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      INTEGER NBTOP, I
      IF(NB_SON(STEP_LOAD(INODE)).EQ.-1)THEN
         RETURN
      ELSE
         IF(NB_SON(STEP_LOAD(INODE)).LT.0)THEN
            WRITE(*,*)
     *        'Internal error 1 in PROCESS_NIV2_MEM_MSG'
            CALL SMUMPS_ABORT()
         ENDIF
      ENDIF
      NB_SON(STEP_LOAD(INODE))=
     *     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     *        SMUMPS_543(INODE)
         POOL_SIZE=POOL_SIZE+1
         IF(POOL_NIV2_COST(POOL_SIZE).GT.MAX_M2)THEN
            MAX_M2=POOL_NIV2_COST(POOL_SIZE)
            ID_MAX_M2=POOL_NIV2(POOL_SIZE)
            CALL SMUMPS_515(REMOVE_NODE_FLAG_MEM,MAX_M2,COMM_LD)
            NIV2(1+MYID)=MAX_M2
         ENDIF
      ENDIF
      END SUBROUTINE PROCESS_NIV2_MEM_MSG      
      RECURSIVE SUBROUTINE PROCESS_NIV2_FLOPS_MSG(INODE)
      IMPLICIT NONE
      INTEGER INODE
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      NB_SON(STEP_LOAD(INODE))=
     *     NB_SON(STEP_LOAD(INODE))-1
      IF(NB_SON(STEP_LOAD(INODE)).EQ.0)THEN
         IF(INODE.EQ.KEEP_LOAD(20))THEN
            RETURN
         ENDIF
         POOL_NIV2(POOL_SIZE+1)=INODE
         POOL_NIV2_COST(POOL_SIZE+1)=
     *        SMUMPS_542(INODE)
         POOL_SIZE=POOL_SIZE+1
         MAX_M2=POOL_NIV2_COST(POOL_SIZE)
         ID_MAX_M2=POOL_NIV2(POOL_SIZE)
         CALL SMUMPS_515(REMOVE_NODE_FLAG,
     $           POOL_NIV2_COST(POOL_SIZE),
     $        COMM_LD)
         NIV2(MYID+1)=POOL_NIV2_COST(POOL_SIZE)+NIV2(MYID+1)
      ENDIF
      END SUBROUTINE PROCESS_NIV2_FLOPS_MSG
      DOUBLE PRECISION FUNCTION SMUMPS_542(INODE)
      INTEGER INODE
      INTEGER NFRONT,NELIM,I,LEVEL
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      DOUBLE PRECISION COST
      I = INODE
      NELIM = 0
 10   CONTINUE
      IF ( I > 0 ) THEN
        NELIM = NELIM + 1
        I = FILS_LOAD(I)
        GOTO 10
      ENDIF
      NFRONT = ND_LOAD( STEP_LOAD(INODE) )
      LEVEL = SMUMPS_330( STEP_LOAD(INODE), PROCNODE_LOAD, NPROCS )
      COST=DBLE(0)
      CALL SMUMPS_511(NFRONT,NELIM,NELIM,
     *                          KEEP_LOAD(50),LEVEL,COST)
      SMUMPS_542=COST
      RETURN
      END FUNCTION SMUMPS_542
      INTEGER FUNCTION SMUMPS_541( INODE )
      IMPLICIT NONE
      INTEGER INODE,NELIM,NFR,SON,IN,I
      INTEGER COST_CB
      EXTERNAL SMUMPS_330
      INTEGER SMUMPS_330
      COST_CB=0
      I = INODE
 10   CONTINUE
      IF ( I > 0 ) THEN
        I = FILS_LOAD(I)
        GOTO 10
      ENDIF
      SON=-I
      DO I=1, NE_LOAD(STEP_LOAD(INODE))
         NFR = ND_LOAD( STEP_LOAD(SON) )
         IN=SON
         NELIM = 0
 20      CONTINUE
         IF ( IN > 0 ) THEN
            NELIM = NELIM + 1
            IN = FILS_LOAD(IN)
            GOTO 20
         ENDIF
         COST_CB=COST_CB+((NFR-NELIM)*(NFR-NELIM))
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      SMUMPS_541=COST_CB
      RETURN
      END FUNCTION SMUMPS_541
      SUBROUTINE SMUMPS_533(SLAVEF,NMB_OF_CAND,
     *     TAB_POS, NASS, KEEP,KEEP8, LIST_SLAVES,
     $     NSLAVES,INODE)
      USE SMUMPS_BUFFER
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, NASS, NSLAVES
      INTEGER, INTENT (IN) :: TAB_POS(SLAVEF+2)
      INTEGER, intent(in) :: NMB_OF_CAND
      INTEGER, INTENT (IN) :: LIST_SLAVES( NMB_OF_CAND )
      INTEGER KEEP(500),INODE,NELIM
      INTEGER*8 KEEP8(150)
      INTEGER allocok
      DOUBLE PRECISION MEM_COST,FCT_COST
      DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE ::EMPTY_ARRAY,
     $     DELTA_MD,EMPTY_ARRAY2
      INTEGER, DIMENSION(:),ALLOCATABLE ::PROCS
      INTEGER NBROWS_SLAVE,I,WHAT,IERR,NFRONT,NPROCS_LOC
      LOGICAL force_cand
      MEM_COST=DBLE(0)
      FCT_COST=DBLE(0)
      IF ( KEEP(24) == 0 .OR. KEEP(24) == 1 ) THEN
        force_cand = .FALSE.
        NPROCS_LOC=SLAVEF-1
      ELSE
        force_cand = (mod(KEEP(24),2).eq.0)
        NPROCS_LOC=NMB_OF_CAND
      END IF
      IF(force_cand)THEN
         CALL SMUMPS_540(INODE,FCT_COST,
     $        MEM_COST,NPROCS_LOC,NASS)
      ELSE
         CALL SMUMPS_540(INODE,FCT_COST,
     $        MEM_COST,SLAVEF-1,NASS)
      ENDIF
      DO I=1,SLAVEF
         IDWLOAD(I)=I-1
      ENDDO
      ALLOCATE(EMPTY_ARRAY(NPROCS_LOC),DELTA_MD(NPROCS_LOC),
     $     EMPTY_ARRAY2(NPROCS_LOC),
     $    stat=allocok)
#if ! defined (DELTA_MD_STK)
      DO I = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(I+1) - TAB_POS(I)
         DELTA_MD( I ) = FCT_COST - DBLE(NBROWS_SLAVE)*
     *        DBLE(NASS)
      END DO
#else
      NFRONT=ND_LOAD(STEP_LOAD(INODE))
      DO I = 1, NSLAVES
         NBROWS_SLAVE = TAB_POS(I+1) - TAB_POS(I)
         IF ( KEEP(50) == 0 ) THEN
            DELTA_MD( I ) = FCT_COST - DBLE(NBROWS_SLAVE)*
     *           DBLE(NFRONT)
         ELSE
            DELTA_MD ( I ) = FCT_COST - DBLE(NBROWS_SLAVE) *
     *           DBLE( NASS + TAB_POS(I+1) - 1 )
         END IF
      END DO
#endif
      IF(FORCE_CAND)THEN
         DO I=NSLAVES+1,NPROCS_LOC
            DELTA_MD( I ) = FCT_COST
         ENDDO
      ELSE
         DO I=NSLAVES+1,SLAVEF-1
            DELTA_MD( I ) = FCT_COST
         ENDDO
      ENDIF
      WHAT=7
 111  CONTINUE
      CALL SMUMPS_524(.FALSE., COMM_LD, MYID, SLAVEF,
#if ! defined(OLD_LOAD_MECHANISM)
     *     FUTURE_NIV2,
#endif
     *     NPROCS_LOC, LIST_SLAVES,0,
     *     EMPTY_ARRAY,
     *     DELTA_MD,EMPTY_ARRAY2,WHAT, IERR)
      IF ( IERR == -1 ) THEN
          CALL SMUMPS_467(COMM_LD, KEEP)
          GOTO 111
      ELSE IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) "Internal Error in SMUMPS_533",
     *   IERR
         CALL SMUMPS_ABORT()
      ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
      IF (FUTURE_NIV2(MYID+1) .NE. 0) THEN
#endif
        DO I = 1, NSLAVES
           MD_MEM(LIST_SLAVES(I))=MD_MEM(LIST_SLAVES(I))+
     $          DELTA_MD( I )
#if ! defined(OLD_LOAD_MECHANISM)
           IF(FUTURE_NIV2(LIST_SLAVES(I)+1).EQ.0)THEN
              MD_MEM(LIST_SLAVES(I))=999999999
           ENDIF
#endif
        ENDDO
#if ! defined(OLD_LOAD_MECHANISM)
      ENDIF
#endif
      DEALLOCATE(EMPTY_ARRAY)
      DEALLOCATE(DELTA_MD)
      END SUBROUTINE SMUMPS_533
      SUBROUTINE SMUMPS_540(INODE,FCT_COST,
     $     MEM_COST,NSLAVES,NELIM)
      IMPLICIT NONE
      INTEGER INODE,NSLAVES,NFR,NELIM,IN
      DOUBLE PRECISION MEM_COST,FCT_COST
      NFR=ND_LOAD(STEP_LOAD(INODE))
      IN = INODE
      FCT_COST=DBLE((INT(DBLE(NFR-NELIM)/DBLE(NSLAVES))+1)*
     $     DBLE(NELIM))
      MEM_COST=DBLE((INT(DBLE(NFR-NELIM)/DBLE(NSLAVES))+1)*
     $     DBLE(NFR))
      END SUBROUTINE SMUMPS_540
      SUBROUTINE CLEAN_POOL_MEM_INFO(INODE)
      IMPLICIT NONE
      INTEGER INODE
      INTEGER I,J,SON,NSLAVES_TEMP,POS_TEMP,K
      INTEGER SMUMPS_275
      EXTERNAL SMUMPS_275
      IF((INODE.LT.0).OR.(INODE.GT.N_LOAD))THEN
         RETURN
      ENDIF
      IF(POS_ID.GT.1)THEN
         I=INODE
 10      CONTINUE
         IF ( I > 0 ) THEN
            I = FILS_LOAD(I)
            GOTO 10
         ENDIF
         SON=-I
         IF(POS_ID.LT.NE_LOAD(STEP_LOAD(INODE))*3)THEN
            I=1
         ENDIF
         DO I=1, NE_LOAD(STEP_LOAD(INODE))
            J=1
            DO WHILE (J.LT.POS_ID)
               IF(CB_COST_ID(J).EQ.SON)GOTO 295 
               J=J+3
            ENDDO
 295        CONTINUE
            IF(J.GE.POS_ID)THEN
               IF(SMUMPS_275(
     $            STEP_LOAD(INODE),PROCNODE_LOAD,NPROCS).EQ.MYID)THEN
                  IF(INODE.EQ.KEEP_LOAD(38))THEN
                     GOTO 666
#if ! defined(OLD_LOAD_MECHANISM)
                  ELSE
                     IF(FUTURE_NIV2(MYID+1).NE.0)THEN
                        WRITE(*,*)MYID,': I did not find ',SON
                        CALL SMUMPS_ABORT()
                     ENDIF
                     GOTO 666
#endif
                  ENDIF
               ELSE
                  GOTO 666
               ENDIF
            ENDIF
            NSLAVES_TEMP=CB_COST_ID(J+1)       
            POS_TEMP=CB_COST_ID(J+2)
            DO K=J,POS_ID-1
               CB_COST_ID(K)=CB_COST_ID(K+3)
            ENDDO
            K=POS_TEMP
            DO WHILE (K.LE.POS_MEM-1)
               CB_COST_MEM(K)=CB_COST_MEM(K+2*NSLAVES_TEMP)
               K=K+1
            ENDDO
            POS_MEM=POS_MEM-2*NSLAVES_TEMP
            POS_ID=POS_ID-3
            IF((POS_MEM.LT.1).OR.(POS_ID.LT.1))THEN
               WRITE(*,*)MYID,': negative pos_mem or pos_id'
               CALL SMUMPS_ABORT()
            ENDIF
 666        CONTINUE
            SON=FRERE_LOAD(STEP_LOAD(SON))
         ENDDO
      ENDIF
      END SUBROUTINE CLEAN_POOL_MEM_INFO
      SUBROUTINE CHECK_MEM_CONST_FOR_POOL(FLAG)
      IMPLICIT NONE
      LOGICAL FLAG
      INTEGER I
      DOUBLE PRECISION MEM
      FLAG=.FALSE.
      DO I=0,NPROCS-1
         MEM=DM_MEM(I)+LU_USAGE(I)
         IF(BDC_SBTR)THEN
            MEM=MEM+SBTR_MEM(I)-SBTR_CUR(I)
         ENDIF
         IF((MEM/TAB_MAXS(I)).GT.0.8d0)THEN
            FLAG=.TRUE.
            GOTO 666
         ENDIF
      ENDDO
 666  CONTINUE
      END SUBROUTINE CHECK_MEM_CONST_FOR_POOL
      SUBROUTINE SMUMPS_554(NBINSUBTREE,INSUBTREE,NBTOP,
     $           MIN_COST,SBTR)
      IMPLICIT NONE
      INTEGER NBINSUBTREE,INSUBTREE,NBTOP
      DOUBLE PRECISION MIN_COST
      LOGICAL SBTR
      INTEGER I
      DOUBLE PRECISION TMP_COST,TMP_MIN
      TMP_MIN=HUGE(TMP_MIN)
      DO I=0,NPROCS-1
         IF(I.NE.MYID)THEN
            IF(BDC_SBTR)THEN
               TMP_MIN=MIN(TMP_MIN,TAB_MAXS(I)-(DM_MEM(I)+
     $              LU_USAGE(I))-(SBTR_MEM(I)-SBTR_CUR(I)))
            ELSE
               TMP_MIN=MIN(TMP_MIN,TAB_MAXS(I)-
     $              (DM_MEM(I)+LU_USAGE(I)))
            ENDIF
         ENDIF
      ENDDO
      IF(NBINSUBTREE.GT.0)THEN
         IF(INSUBTREE.EQ.1)THEN
            TMP_COST=TAB_MAXS(MYID)-(DM_MEM(MYID)+LU_USAGE(MYID))
     $           -(SBTR_MEM(MYID)-SBTR_CUR(MYID))
         ELSE
            SBTR=.FALSE.
            GOTO 777
         ENDIF
      ENDIF
      TMP_MIN=MIN(TMP_COST,TMP_MIN)
      IF(TMP_MIN.GT.MIN_COST) SBTR=.TRUE.
 777  CONTINUE
      END SUBROUTINE SMUMPS_554
      SUBROUTINE COMPUTE_MAX_MEM(INODE,MAX_MEM,PROC)
      IMPLICIT NONE
      INTEGER INODE,PROC
      INTEGER I,POS,NSLAVES,SLAVE,NCAND,J,NELIM,NCB,NFRONT,SON,K
      INTEGER allocok
      EXTERNAL SMUMPS_330
      INTEGER  SMUMPS_330
      DOUBLE PRECISION MAX_MEM,TMP
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MEM_ON_PROCS,
     $     RECV_BUF
      LOGICAL, DIMENSION(:), ALLOCATABLE :: CONCERNED
      DOUBLE PRECISION MAX_SENT_MSG
      IF((FRERE_LOAD(STEP_LOAD(INODE)).EQ.0)
     $           .AND.(INODE.EQ.KEEP_LOAD(38)))THEN
         RETURN
      ENDIF
#if  defined(NOT_ATM_POOL_SPECIAL)
      IF((INODE.LT.0).OR.(INODE.GT.N_LOAD))THEN
         MAX_MEM=HUGE(MAX_MEM)
         DO I=0,NPROCS-1
            TMP=TAB_MAXS(I)-(DM_MEM(I)+LU_USAGE(I))
            IF(BDC_SBTR)THEN
               TMP=TMP-(SBTR_MEM(I)-SBTR_CUR(I))
            ENDIF
            MAX_MEM=MIN(MAX_MEM,TMP)
         ENDDO
         RETURN
      ENDIF
#endif
      ALLOCATE( MEM_ON_PROCS(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in COMPUTE_MAX_MEM'
        CALL SMUMPS_ABORT()
      ENDIF 
      ALLOCATE( CONCERNED(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in COMPUTE_MAX_MEM'
        CALL SMUMPS_ABORT()
      ENDIF 
      ALLOCATE( RECV_BUF(0:NPROCS-1), stat=allocok)
      IF ( allocok > 0 ) THEN
        WRITE(*,*) 'PB allocation in COMPUTE_MAX_MEM'
        CALL SMUMPS_ABORT()
      ENDIF 
      RECV_BUF=DBLE(0)
      MAX_SENT_MSG=DBLE(0)
      I = INODE
      NELIM = 0
 10   CONTINUE
      IF ( I > 0 ) THEN
        NELIM = NELIM + 1
        I = FILS_LOAD(I)
        GOTO 10
      ENDIF
      SON=-I
      NFRONT=ND_LOAD(STEP_LOAD(INODE))
      NCB=NFRONT-NELIM
      IF(SMUMPS_330(STEP_LOAD(INODE)
     $     ,PROCNODE_LOAD,NPROCS).EQ.2)THEN
         NCAND=CAND_LOAD(NPROCS+1, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
      ENDIF
      DO I=0,NPROCS-1
         IF(I.EQ.MYID)THEN
            MEM_ON_PROCS(I)=TAB_MAXS(I)-(DM_MEM(I)+LU_USAGE(I)+
     $           SMUMPS_543(INODE))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(I)=MEM_ON_PROCS(I)-(SBTR_MEM(I)-SBTR_CUR(I))
            ENDIF
            CONCERNED(I)=.TRUE.
         ELSE
            MEM_ON_PROCS(I)=TAB_MAXS(I)-(DM_MEM(I)+LU_USAGE(I))
            IF(BDC_SBTR)THEN
               MEM_ON_PROCS(I)=MEM_ON_PROCS(I)-(SBTR_MEM(I)-SBTR_CUR(I))
            ENDIF
            IF(BDC_M2_MEM)THEN
               MEM_ON_PROCS(I)=MEM_ON_PROCS(I)-NIV2(I+1)
            ENDIF
         ENDIF
         IF(SMUMPS_330(STEP_LOAD(INODE)
     $        ,PROCNODE_LOAD,NPROCS).EQ.2)THEN
            IF(BDC_MD.AND.(KEEP_LOAD(48).EQ.5))THEN
               DO J=1,NCAND
                  IF(CAND_LOAD(J, STEP_TO_NIV2_LOAD(STEP_LOAD(INODE)))
     $                 .EQ.I)THEN
                     MEM_ON_PROCS(I)=MEM_ON_PROCS(I)-
     $                 ((DBLE(NFRONT)*DBLE(NCB))/DBLE(NCAND))
                     CONCERNED(I)=.TRUE.
                     GOTO 666
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
 666     CONTINUE
      ENDDO
      DO K=1, NE_LOAD(STEP_LOAD(INODE))
         I=1
         DO WHILE (I.LE.POS_ID)
            IF(CB_COST_ID(I).EQ.SON)GOTO 295
            I=I+3
         ENDDO
 295     CONTINUE
         IF(I.GE.POS_ID)THEN
#if ! defined(OLD_LOAD_MECHANISM)
            IF(FUTURE_NIV2(MYID+1).NE.0)THEN
               WRITE(*,*)MYID,': ',SON,'has not been found
     $              in COMPUTE_MAX_MEM'
               CALL SMUMPS_ABORT()
            ENDIF
#endif
            GOTO 777
         ENDIF
         NSLAVES=CB_COST_ID(I+1)
         POS=CB_COST_ID(I+2)
         DO I=1,NSLAVES
            SLAVE=CB_COST_MEM(POS)
            IF(.NOT.CONCERNED(SLAVE))THEN
               MEM_ON_PROCS(SLAVE)=MEM_ON_PROCS(SLAVE)+
     $              CB_COST_MEM(POS+1)
            ENDIF
            DO J=0,NPROCS-1
               IF(CONCERNED(J))THEN
                  IF(SLAVE.NE.J)THEN
                     RECV_BUF(J)=MAX(RECV_BUF(J),
     $                    CB_COST_MEM(POS+1))
                  ENDIF
               ENDIF
            ENDDO
            POS=POS+2
         ENDDO
 777     CONTINUE
         SON=FRERE_LOAD(STEP_LOAD(SON))
      ENDDO
      MAX_MEM=HUGE(MAX_MEM)
      DO I=0,NPROCS-1
         IF(MAX_MEM.GT.MEM_ON_PROCS(I))THEN
            PROC=I
         ENDIF
         MAX_MEM=MIN(MEM_ON_PROCS(I),MAX_MEM)
      ENDDO
      DEALLOCATE(MEM_ON_PROCS)
      DEALLOCATE(CONCERNED)
      DEALLOCATE(RECV_BUF)
      END SUBROUTINE COMPUTE_MAX_MEM
      SUBROUTINE SMUMPS_553(MIN_PROC,POOL,
     $                      LPOOL,INODE)
      IMPLICIT NONE
      INTEGER INODE,LPOOL,MIN_PROC
      INTEGER POOL(LPOOL)
      EXTERNAL SMUMPS_275
      INTEGER SMUMPS_275
      INTEGER I,NBTOP,INSUBTREE,NBINSUBTREE,NODE,FATHER,SON,J
      INTEGER FIRST_LEAF,SBTR_NB_LEAF,POS,K,allocok,L
      INTEGER, ALLOCATABLE, DIMENSION (:) ::  TMP_SBTR
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF((KEEP_LOAD(47).EQ.4).AND.
     $     ((NBINSUBTREE.NE.0)))THEN
         DO J=INDICE_SBTR,NB_SUBTREES
            NODE=MY_ROOT_SBTR(J) 
            FATHER=DAD_LOAD(STEP_LOAD(NODE))
            I=FATHER
 110        CONTINUE
            IF ( I > 0 ) THEN
               I = FILS_LOAD(I)
               GOTO 110
            ENDIF
            SON=-I
            I=SON
 120        CONTINUE
            IF ( I > 0 ) THEN
               IF(SMUMPS_275(STEP_LOAD(I),PROCNODE_LOAD,NPROCS).EQ.
     $              MIN_PROC)THEN
                  SBTR_NB_LEAF=MY_NB_LEAF(J)
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  IF(POOL(POS+SBTR_NB_LEAF).NE.MY_FIRST_LEAF(J))THEN
                     WRITE(*,*)MYID,': The first leaf is not ok'
                     CALL SMUMPS_ABORT()
                  ENDIF
                  ALLOCATE (TMP_SBTR(SBTR_NB_LEAF), stat=allocok)
                  IF (allocok > 0 ) THEN
                     WRITE(*,*)MYID,': Not enough space
     $                    for allocation'
                     CALL SMUMPS_ABORT()
                  ENDIF
                  POS=SBTR_FIRST_POS_IN_POOL(J)
                  DO K=1,SBTR_NB_LEAF
                     TMP_SBTR(K)=POOL(POS+K-1)
                  ENDDO
                  DO K=POS+1,NBINSUBTREE-SBTR_NB_LEAF
                     POOL(K)=POOL(K+SBTR_NB_LEAF)
                  ENDDO
                  POS=1
                  DO K=NBINSUBTREE-SBTR_NB_LEAF+1,NBINSUBTREE
                     POOL(K)=TMP_SBTR(POS)
                     POS=POS+1
                  ENDDO
                  DO K=INDICE_SBTR,J
                     SBTR_FIRST_POS_IN_POOL(K)=SBTR_FIRST_POS_IN_POOL(K)
     $                    -SBTR_FIRST_POS_IN_POOL(J)
                  ENDDO
                  SBTR_FIRST_POS_IN_POOL(J)=NBINSUBTREE-SBTR_NB_LEAF
                  POS=MY_FIRST_LEAF(J)
                  L=MY_NB_LEAF(J)
                  DO K=INDICE_SBTR,J
                     MY_FIRST_LEAF(J)=MY_FIRST_LEAF(J+1)
                     MY_NB_LEAF(J)=MY_NB_LEAF(J+1)
                  ENDDO
                  MY_FIRST_LEAF(INDICE_SBTR)=POS
                  MY_NB_LEAF(INDICE_SBTR)=L
                  INODE=POOL(NBINSUBTREE)
                  DEALLOCATE(TMP_SBTR)
                  RETURN
               ENDIF
               I = FRERE_LOAD(STEP_LOAD(I))
               GOTO 120
            ENDIF           
         ENDDO
      ENDIF
      DO J=NBTOP,1,-1
#if defined(NOT_ATM_POOL_SPECIAL)
         IF ( POOL(LPOOL-2-J) < 0 ) THEN
            NODE=-POOL(LPOOL-2-J)
         ELSE IF ( POOL(LPOOL-2-J) > N_LOAD ) THEN
            NODE = POOL(LPOOL-2-J) - N_LOAD
         ELSE
            NODE = POOL(LPOOL-2-J)
         ENDIF
#else
         NODE=POOL(LPOOL-2-J)
#endif
         FATHER=DAD_LOAD(STEP_LOAD(NODE))
         I=FATHER
 11      CONTINUE
         IF ( I > 0 ) THEN
            I = FILS_LOAD(I)
            GOTO 11
         ENDIF
         SON=-I
         I=SON
 12      CONTINUE
         IF ( I > 0 ) THEN
            IF(SMUMPS_275(STEP_LOAD(I),PROCNODE_LOAD,NPROCS).EQ.
     $           MIN_PROC)THEN
               INODE=NODE
               RETURN
            ENDIF
            I = FRERE_LOAD(STEP_LOAD(I))
            GOTO 12
         ENDIF
      ENDDO
      END SUBROUTINE SMUMPS_553
      SUBROUTINE SMUMPS_555(POOL, LPOOL,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER LPOOL,POOL(LPOOL),KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER I,POS
      EXTERNAL SMUMPS_283
      LOGICAL SMUMPS_283
      IF(.NOT.BDC_SBTR) RETURN
      POS=0
      DO I=NB_SUBTREES,1,-1
         DO WHILE(SMUMPS_283(STEP_LOAD(POOL(POS+1)),
     $        PROCNODE_LOAD,NPROCS))
            POS=POS+1
         ENDDO
         SBTR_FIRST_POS_IN_POOL(I)=POS+1
         POS=POS+MY_NB_LEAF(I)
      ENDDO
      END SUBROUTINE SMUMPS_555
      END MODULE SMUMPS_LOAD 
