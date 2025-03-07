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
      SUBROUTINE CMUMPS_267( BUFR, LBUFR, LBUFR_BYTES,
     *                           N, ID, MYID, STEP,
     *                           PTLUST_S, 
     *                           K10, K28, SLAVEF, IW, LIW,
     *                           A, LA, COMM )
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER K10, K28, LA, LIW, COMM
      INTEGER ID, MYID, SLAVEF, N 
      INTEGER BUFR( LBUFR ) 
      INTEGER STEP(N), PTLUST_S( K28 ), IW( LIW )
      COMPLEX A( LA )
      INCLUDE 'mpif.h'
      INTEGER POSITION, IOLDPS, IERR
      INTEGER NELIM, LCONT, NROW, NPIV, IAPOS, NSLAVES, ISLAVE
      INTEGER INODE, NCOL_L, NROW_L, NCB, BLSIZE, LASTBL
      INTEGER SIZELUMASTER, SIZELUSLAVE, NFRONT
      INCLUDE 'mumps_headers.h'
      ISLAVE = ID - MYID
      IF ( ISLAVE .LT. 1 ) ISLAVE = ISLAVE + SLAVEF 
      POSITION = 0
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 INODE, 1, MPI_INTEGER,
     *                 COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 NCOL_L, 1, MPI_INTEGER,
     *                 COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 NROW_L, 1, MPI_INTEGER,
     *                 COMM, IERR )
      IOLDPS = PTLUST_S(STEP( INODE ))
      LCONT  = IW( IOLDPS +XSIZE)
      NELIM  = IW( IOLDPS + 1 +XSIZE)
      NROW   = IW( IOLDPS + 2 +XSIZE)
      NPIV   = IW( IOLDPS + 3 +XSIZE)
      IAPOS  = IW( IOLDPS + 4 +XSIZE)
      NSLAVES= IW( IOLDPS + 5 +XSIZE)
      NFRONT = LCONT + NPIV
      NCB    = NFRONT - NROW
      IF ( NSLAVES .NE. MIN(NCB/K10,SLAVEF-1)) THEN
        WRITE(*,*) 'Process_factor: error1. Wrong number of slaves'
        WRITE(*,*) 'INODE,NFRONT,NROW,NSLAVES=',
     * INODE,NFRONT,NCB,NSLAVES
        stop
      END IF
      BLSIZE = NCB / NSLAVES
      LASTBL = MOD(NCB,NSLAVES) + BLSIZE
      SIZELUMASTER = ( LCONT + NROW ) * NPIV
      SIZELUSLAVE  = NPIV * BLSIZE
         IF ( NROW_L .NE. LASTBL .AND. NROW_L .NE. BLSIZE ) THEN
           WRITE(*,*) 'Process_factor: error2. wrong block size'
           CALL CMUMPS_ABORT()
         END IF
         IF ( NSLAVES .eq. 0 ) THEN
           WRITE(*,*) 'Process_factor: error3'
           CALL CMUMPS_ABORT()
         END IF
         IF ( NPIV .NE. NCOL_L ) THEN
           WRITE(*,*) 'Process_factor: error4'
           CALL CMUMPS_ABORT()
         END IF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 A(IAPOS + SIZELUMASTER + ( ISLAVE - 1 )
     *                 * SIZELUSLAVE),
     *                 NCOL_L * NROW_L,
     *                 MPI_COMPLEX,
     *                 COMM, IERR )
      RETURN
      END SUBROUTINE CMUMPS_267
      SUBROUTINE CMUMPS_152(SSARBR, MYID, N, IPOSBLOCK,
     *       RPOSBLOCK,
     *       IW, LIW, NIRBDU,
     *       LRLU, LRLUS, IPTRLU,
     *       IWPOSCB, LA, KEEP, KEEP8
     $     )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER IPOSBLOCK, RPOSBLOCK,
     *         LIW, NIRBDU, LRLU, LRLUS,
     *        IPTRLU, IWPOSCB, N, LA
      INTEGER IW( LIW ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER MYID
      LOGICAL SSARBR
      INTEGER SIZFI_BLOCK, SIZFR_BLOCK, SIZFI, SIZFR
      INTEGER IPOSSHIFT
      INTEGER SIZFR_BLOCK_EFF, SIZEHOLE
      INCLUDE 'mumps_headers.h'
      IPOSSHIFT = IPOSBLOCK + XSIZE
      SIZFI_BLOCK=IW(IPOSBLOCK+XXI)
      SIZFR_BLOCK=IW(IPOSBLOCK+XXR)
      IF (KEEP(216).eq.3) THEN
        SIZFR_BLOCK_EFF=SIZFR_BLOCK
      ELSE
        CALL CMUMPS_628( IW(IPOSBLOCK),
     *                     LIW-IPOSBLOCK+1,
     *                     SIZEHOLE)
        SIZFR_BLOCK_EFF=SIZFR_BLOCK-SIZEHOLE
      ENDIF
      IF ( IPOSBLOCK .eq. IWPOSCB + 1 ) THEN
         IPTRLU  = IPTRLU  + SIZFR_BLOCK
         IWPOSCB = IWPOSCB + SIZFI_BLOCK
         LRLU    = LRLU  + SIZFR_BLOCK
         LRLUS   = LRLUS + SIZFR_BLOCK_EFF
      CALL CMUMPS_471(SSARBR,.FALSE.,
     *         LA-LRLUS,0,-SIZFR_BLOCK_EFF,KEEP,KEEP8,LRLU)
 90      IF ( IWPOSCB .eq. NIRBDU ) GO TO 100
         IPOSSHIFT = IWPOSCB + XSIZE
         SIZFI = IW( IWPOSCB+1+XXI )
         SIZFR = IW( IWPOSCB+1+XXR )
         IF ( IW( IWPOSCB+1+XXS ) .EQ. S_FREE ) THEN
              IPTRLU  = IPTRLU + SIZFR
              LRLU    = LRLU + SIZFR
              IWPOSCB = IWPOSCB + SIZFI
              GO TO 90
         ENDIF
 100     CONTINUE
         IW( IWPOSCB+1+XXP)=TOP_OF_STACK
      ELSE
         IW( IPOSBLOCK +XXS)=S_FREE
         LRLUS = LRLUS + SIZFR_BLOCK_EFF
      CALL CMUMPS_471(SSARBR,.FALSE.,
     *            LA-LRLUS,0,-SIZFR_BLOCK_EFF,KEEP,KEEP8,LRLU)
      END IF
      RETURN
      END SUBROUTINE CMUMPS_152
      SUBROUTINE CMUMPS_144( COMM_LOAD, ASS_IRECV, 
     *           N, INODE, FPERE, IW, LIW, A, LA,
     *           UU, NOFFW,
     *           NPVW,
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP, PIMASTER,
     *             PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     *              DKEEP,PIVNUL_LIST,LPN_LIST)
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, LA, NOFFW, NPVW
      INTEGER IW( LIW )
      COMPLEX A( LA )
      REAL UU, SEUIL
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER NBFIN, SLAVEF, 
     *        IFLAG, IERROR, LEAF, LPOOL
      INTEGER POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, 
     *        NIRBDU, COMP 
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     *        ITLOC(N), FILS(N),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR),
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER  INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     *        PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        STEP(N),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     *        PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      LOGICAL AVOID_DELAYED
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEGKJI, NBTLKJ,IBEG_BLOCK
      INTEGER NASS, NEL1, IEND, IOLDPS, POSELT, dummy, allocok
      LOGICAL LASTBL 
      REAL UUTEMP
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INCLUDE 'mumps_headers.h'
      EXTERNAL CMUMPS_224, CMUMPS_233, 
     *         CMUMPS_225, CMUMPS_232,
     *         CMUMPS_294,
     *         CMUMPS_44
      LOGICAL STATICMODE
      REAL SEUIL_LOC
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
      IBEG_BLOCK=1
      dummy  = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = IABS(IW(IOLDPS+2+XSIZE))
      IF (NASS .GT. KEEP(3)) THEN
        NBTLKJ = MIN( KEEP(6), NASS )
      ELSE
        NBTLKJ = MIN( KEEP(5),NASS )
      ENDIF
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID,' : FACTO_NIV2 :failed to allocate ',NASS,
     * ' integers'
        IFLAG  = -13
        IERROR =NASS
        GO TO 490
      END IF
 50   CONTINUE
      IBEGKJI = IBEG_BLOCK
      CALL CMUMPS_224(NFRONT,NASS,IBEGKJI, NASS, IPIV,
     *                N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     *                 DKEEP(1),PIVNUL_LIST(1),LPN_LIST)
      IF (IFLAG.LT.0) GOTO 490
      IF (INOPV.EQ.1) THEN
         IF(STATICMODE) THEN
            INOPV = -1
            GOTO 50
         ENDIF
      ENDIF
      IF (INOPV.GE.1) THEN
          LASTBL = (INOPV.EQ.1)
          IEND = IW(IOLDPS+1+XSIZE)
          CALL CMUMPS_294( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, NFRONT,
     *             IBEGKJI, IEND, IPIV, NASS,LASTBL, dummy, 
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
         CALL CMUMPS_233(IBEG_BLOCK,NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *            IOLDPS,POSELT,NBTLKJ,KEEP(4))
         GOTO 50
      ENDIF
      NPVW = NPVW + 1
      IF (NASS.LE.1) THEN
        IFINB = -1
      ELSE
         CALL CMUMPS_225(IBEG_BLOCK,
     *             NFRONT, NASS, N,INODE,IW,LIW,A,LA,
     *             IOLDPS,POSELT,IFINB,
     *             NBTLKJ,KEEP(4))
      ENDIF
      IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
       IF (IFINB.EQ.0) GOTO 50
       IF ((IFINB.EQ.1).OR.(IFINB.EQ.-1)) THEN
          LASTBL = (IFINB.EQ.-1) 
          IEND = IW(IOLDPS+1+XSIZE)
          CALL CMUMPS_294(COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE, IW, LIW, 
     *             IOLDPS, POSELT, A, LA, NFRONT, 
     *             IBEGKJI, IEND, IPIV, NASS, LASTBL, dummy, 
     *
     *             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     *             IFLAG, IERROR, IPOOL,LPOOL,
     *             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *             LRLUS, NIRBDU, COMP,
     *             PTRIST, PTRAST, PTLUST_S, PTRFAC,
     *             STEP, PIMASTER, PAMASTER,
     *             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     *             OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *             LPTRAR, NELT, FRTPTR, FRTELT, 
     *             ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 500
       ENDIF
       IF (IFINB.EQ.(-1)) GOTO 500
       NPIV   = IW(IOLDPS+1+XSIZE)
       NEL1   = NASS - NPIV
      CALL CMUMPS_232(A,LA,
     *           NFRONT,NPIV,NASS,POSELT,NBTLKJ)
      GO TO 50
 490  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
 500  CONTINUE
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE CMUMPS_144
      SUBROUTINE CMUMPS_176( COMM_LOAD, ASS_IRECV, 
     *    ROOT, FRERE, IROOT, 
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
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: ROOT
      INTEGER IROOT
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     * PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND(KEEP(28)), FRERE(KEEP(28))
      INTEGER INTARR( MAX(1,KEEP(14)) )
      COMPLEX DBLARR( MAX(1,KEEP(13)) )
      INTEGER NELIM, NB_CONTRI_GLOBAL, NUMORG, 
     *        NFRONT, IROW, JCOL, PDEST, HF, IOLDPS,
     *        IN, DEB_ROW, ILOC_ROW, IFSON, ILOC_COL,
     *        IPOS_SON, NELIM_SON, NSLAVES_SON, HS,
     *        IROW_SON, ICOL_SON, ISLAVE, IERR, 
     *        NELIM_SENT
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      NB_CONTRI_GLOBAL = KEEP(41)
      NUMORG    = root%ROOT_SIZE
      NELIM     = KEEP(42)
      NFRONT    = NUMORG + KEEP(42)
      DO IROW = 0, root%NPROW - 1
        DO JCOL = 0, root%NPCOL - 1
            PDEST = IROW * root%NPCOL + JCOL
          IF ( PDEST .NE. MYID ) THEN
           CALL CMUMPS_73(NFRONT, 
     *     NB_CONTRI_GLOBAL, PDEST, COMM, IERR)
              if (ierr.lt.0) then
                write(6,*) ' error detected by ',
     *          'CMUMPS_73'
                stop
               endif
           ENDIF
        END DO
      END DO
      CALL  CMUMPS_270( NFRONT,
     *    NB_CONTRI_GLOBAL, ROOT,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     *    IFLAG, IERROR, COMM, COMM_LOAD,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND )
       IF (IFLAG < 0 ) RETURN
      HF = 6 + XSIZE
      IOLDPS = PTLUST_S(STEP(IROOT))
      IN = IROOT
      DEB_ROW = IOLDPS + HF
      ILOC_ROW    = DEB_ROW
      DO WHILE (IN.GT.0) 
       IW(ILOC_ROW)           = IN
       IW(ILOC_ROW+NFRONT)    = IN
       ILOC_ROW = ILOC_ROW + 1
       IN = FILS(IN)
      END DO
      IFSON = -IN
      ILOC_ROW    = IOLDPS + HF + NUMORG
      ILOC_COL    = ILOC_ROW + NFRONT
      IF ( NELIM.GT.0 ) THEN
        IN = IFSON
        DO WHILE (IN.GT.0)
          IPOS_SON  = PIMASTER(STEP(IN))
          IF (IPOS_SON .EQ. 0) GOTO 100
          NELIM_SON   = IW(IPOS_SON+1+XSIZE)
              if (nelim_son.eq.0) then
                write(6,*) ' error 1 in process_last_rtnelind'
                stop
              endif
          NSLAVES_SON = IW(IPOS_SON+5+XSIZE)
          HS          = 6 + NSLAVES_SON + XSIZE
          IROW_SON    = IPOS_SON + HS
          ICOL_SON    = IROW_SON + NELIM_SON
          IW(ILOC_ROW:ILOC_ROW+NELIM_SON-1) = 
     *            IW(IROW_SON:IROW_SON+NELIM_SON-1) 
          IW(ILOC_COL:ILOC_COL+NELIM_SON-1) = 
     *            IW(ICOL_SON:ICOL_SON+NELIM_SON-1) 
          NELIM_SENT = ILOC_ROW - IOLDPS - HF + 1
          DO ISLAVE = 0,NSLAVES_SON
            IF (ISLAVE.EQ.0) THEN
             PDEST= CMUMPS_275(STEP(IN),PROCNODE_STEPS,SLAVEF)
            ELSE
             PDEST = IW(IPOS_SON + 5 + ISLAVE+XSIZE)
            ENDIF
            IF (PDEST.NE.MYID) THEN
             CALL CMUMPS_74(IN, NELIM_SENT,
     *        PDEST, COMM, IERR )
               if (ierr.lt.0) then
                write(6,*) ' error detected by ',
     *          'CMUMPS_73'
                stop
               endif
            ELSE
             CALL CMUMPS_271( COMM_LOAD, ASS_IRECV,
     *       IN, NELIM_SENT, ROOT,
     *
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
     *       OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *       INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *       LPTRAR, NELT, FRTPTR, FRTELT, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
             IF ( ISLAVE .NE. 0 ) THEN
               CALL CMUMPS_626( N, IN, PTRIST, PTRAST,
     *         IW, LIW, A, LA, LRLU, LRLUS, IWPOSCB,
     *         NIRBDU, IPTRLU, STEP, MYID, KEEP
     *        )
             ENDIF
             IPOS_SON  = PIMASTER(STEP(IN))
            ENDIF
          END DO
          CALL  CMUMPS_152( .FALSE.,MYID,N, IPOS_SON, 
     *       PTRAST(STEP(IN)),
     *       IW, LIW, NIRBDU,
     *       LRLU, LRLUS, IPTRLU,
     *       IWPOSCB, LA, KEEP,KEEP8
     $         )
          ILOC_ROW = ILOC_ROW + NELIM_SON
          ILOC_COL = ILOC_COL + NELIM_SON
 100      CONTINUE
          IN = FRERE(STEP(IN))
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_176
      SUBROUTINE CMUMPS_276( ICNTL, INFO, COMM, ID )
      INTEGER ICNTL(40), INFO(40), COMM, ID
      INCLUDE 'mpif.h'
      INTEGER IN( 2 ), OUT( 2 ), LP, IERR
      LP      = ICNTL( 1 )
      IN( 1 ) = INFO ( 1 )
      IN( 2 ) = ID
      CALL MPI_ALLREDUCE( IN, OUT, 1, MPI_2INTEGER, MPI_MINLOC,
     &                    COMM, IERR)
      IF ( OUT( 1 ) .LT. 0 .and. INFO(1) .GE. 0 ) THEN
        INFO( 1 ) = -001
        INFO( 2 ) = OUT( 2 )
      END IF
      RETURN
      END SUBROUTINE CMUMPS_276
      SUBROUTINE CMUMPS_268(MYID,BUFR, LBUFR, 
     *     LBUFR_BYTES,
     *     PROCNODE_STEPS, SLAVEF,
     *     IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     *     N, IW, LIW, A, LA, NIRBDU,
     *     PTRIST, PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S,
     *     COMP,
     *     IFLAG, IERROR, COMM, COMM_LOAD, NBPROCFILS,
     *     IPOOL, LPOOL, LEAF, KEEP,KEEP8, ND, FILS, FRERE, ITLOC, 
     *
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), STEP(N),
     *        
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) ), ITLOC( N )
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER ND(KEEP(28)), FILS( N ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER POSITION, IFATH, ISON, NROW, NCOL, NELIM, 
     *        NSLAVES
      INTEGER NOREAL, NOINT, INIV2
      DOUBLE PRECISION FLOP1
      LOGICAL PERETYPE2
      INCLUDE 'mumps_headers.h'
      INTEGER  CMUMPS_330
      EXTERNAL CMUMPS_330
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          IFATH, 1, MPI_INTEGER
     *        , COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        ISON , 1, MPI_INTEGER, 
     *        COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        NSLAVES, 1,
     *        MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          NROW , 1, MPI_INTEGER
     *        , COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *          NCOL , 1, MPI_INTEGER
     *        , COMM, IERR)
      NOINT = 6 + NROW + NCOL + NSLAVES + XSIZE
      IF ( NSLAVES .NE. 0 .and. KEEP(50).ne.0 ) THEN
        NOREAL = NROW * NROW
      ELSE
        NOREAL = NROW * NCOL
      ENDIF
      CALL CMUMPS_22(.FALSE.,.FALSE.,
     *   MYID,N,KEEP,KEEP8,IW,LIW,A,LA,NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER, ITLOC,
     *   NOINT, NOREAL, ISON, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR )
      IF ( IFLAG .LT. 0 ) THEN
        RETURN
      ENDIF
      PIMASTER(STEP( ISON )) = IWPOSCB + 1
      PAMASTER(STEP( ISON )) = IPTRLU  + 1
      IW( IWPOSCB + 1 + XSIZE ) = NCOL
      NELIM = NROW
      IW( IWPOSCB + 2 + XSIZE ) = NELIM
      IW( IWPOSCB + 3 + XSIZE ) = NROW
      IF ( NSLAVES .NE. 0 .and. KEEP(50).ne.0 ) THEN
        IW( IWPOSCB + 4 + XSIZE ) = NROW - NCOL
        IF ( NROW - NCOL .GE. 0 ) THEN
          WRITE(*,*) 'Error in PROCESS_MAITRE2:',NROW,NCOL
          CALL CMUMPS_ABORT()
        END IF
      ELSE
        IW( IWPOSCB + 4 + XSIZE ) = 0
      END IF
      IW( IWPOSCB + 5 + XSIZE ) = 1
      IW( IWPOSCB + 6 + XSIZE ) = NSLAVES
      IF (NSLAVES.GT.0) THEN
       CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                 IW( IWPOSCB + 7 + XSIZE ),
     *                 NSLAVES, MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 7 + XSIZE + NSLAVES),
     *        NROW, MPI_INTEGER, COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 7 + XSIZE + NROW + NSLAVES),
     *        NCOL, MPI_INTEGER, COMM, IERR)
       IF ( ( KEEP(48).NE. 0 ).AND.(NSLAVES .GT. 0 )) THEN
          INIV2 = ISTEP_TO_INIV2 ( STEP(ISON) )
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        TAB_POS_IN_PERE(1,INIV2),
     *        NSLAVES+1, MPI_INTEGER, COMM, IERR)
          TAB_POS_IN_PERE(SLAVEF+2,INIV2) = NSLAVES
       ENDIF
      IF (NOREAL.GT.0) THEN
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION, A(IPTRLU + 1),
     *        NOREAL, MPI_COMPLEX, COMM, IERR)
      ENDIF
      PERETYPE2 = ( CMUMPS_330(STEP(IFATH),PROCNODE_STEPS,
     *              SLAVEF) .EQ. 2 )
      NSTK_S( STEP(IFATH ))       = NSTK_S( STEP(IFATH) ) - 1
      IF ( NSTK_S( STEP(IFATH)) .EQ. 0 ) THEN
          CALL CMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     *         SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     *         STEP, IFATH )
          IF (KEEP(47) .GE. 3) THEN
             CALL CMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *            MYID, STEP, N, ND, FILS )
          ENDIF
          CALL CMUMPS_137( IFATH, N, PROCNODE_STEPS,
     *                            SLAVEF, ND,
     *                            FILS,FRERE, STEP, PIMASTER,
     *                            KEEP(28), KEEP(50),
     *                            FLOP1,IW, LIW )
          IF (IFATH.NE.KEEP(20))
     *    CALL CMUMPS_190(1, .FALSE., FLOP1, KEEP,KEEP8)
      END IF
      RETURN
      END SUBROUTINE CMUMPS_268
      SUBROUTINE CMUMPS_242(DATA, LDATA, MPITYPE, ROOT, COMMW, TAG,
     *SLAVEF)
      USE CMUMPS_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER LDATA, ROOT, COMMW, TAG, MPITYPE, SLAVEF
      INTEGER DEST
      INTEGER DATA(LDATA)
      DO 10 DEST = 0, SLAVEF - 1
        IF (DEST .NE. ROOT) THEN
          IF ( LDATA .EQ. 1 .and. MPITYPE .EQ. MPI_INTEGER ) THEN
            CALL CMUMPS_62( DATA(1), DEST, TAG, 
     *                                COMMW, IERR )
          ELSE
            WRITE(*,*) 'Error : bad argument to CMUMPS_242'
            CALL CMUMPS_ABORT()
          END IF
        ENDIF
   10 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_242
      SUBROUTINE CMUMPS_44( MYID, SLAVEF, COMM )
      INTEGER MYID, SLAVEF, COMM
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY (1)
      CALL CMUMPS_242( DUMMY, 1, MPI_INTEGER, MYID,
     *                 COMM, TERREUR, SLAVEF )
      RETURN
      END SUBROUTINE CMUMPS_44
      SUBROUTINE CMUMPS_20( NPROCS, CNTL, ICNTL,
     *                    KEEP,KEEP8,
     *                    INFO, INFOG, RINFO, RINFOG, SYM, PAR,
     *                    DKEEP)
      REAL    DKEEP(30)
      REAL    CNTL(15), RINFO(20), RINFOG(20)
      INTEGER ICNTL(40), KEEP(500), SYM, PAR, NPROCS
      INTEGER INFO(40), INFOG(40)
      INTEGER*8 KEEP8(150)
      KEEP(1:500) = 0
      KEEP8(1:150)= 0
      INFO(1:40)  = 0
      INFOG(1:40) = 0
      ICNTL(1:40) = 0
      RINFO(1:20) = 0.0E0
      RINFOG(1:20)= 0.0E0
      CNTL(1:15)   = 0.0E0
      DKEEP(1:30) = 0.0E0
      KEEP( 50 ) = SYM
      IF (SYM.EQ.1) THEN
          SYM      = 2
          KEEP(50) = 2
      ENDIF
      IF ( KEEP(50).NE.1 .and. KEEP(50).NE.2 ) KEEP( 50 ) = 0
      IF ( KEEP(50) .NE. 1 ) THEN
        CNTL(1)   = 0.01E0
      ELSE
        CNTL(1)   = 0.0E0
      END IF
      CNTL(2) = SQRT(EPSILON(0.0E0))
      CNTL(3) = -1.0E0
      CNTL(4) = -1.0E0
      CNTL(5) = 0.0E0
      CNTL(6) = -1.0E0
      KEEP(46) = PAR
      IF ( KEEP(46) .NE. 0 .AND.
     *     KEEP(46) .NE. 1 ) THEN
           KEEP(46) = 1
      END IF
      ICNTL(1)  = 6
      ICNTL(2)  = 0
      ICNTL(3)  = 6
      ICNTL(4)  = 2
      ICNTL(5)  = 0
      IF (SYM.NE.1) THEN
       ICNTL(6)  = 7
      ELSE
       ICNTL(6)  = 0
      ENDIF
      ICNTL(7) = 7
      IF (SYM.NE.1) THEN
        ICNTL(8)  = 7
      ELSE
        ICNTL(8)  = 0
      ENDIF
      ICNTL(9)  = 1
      ICNTL(10)  = 0
      ICNTL(11)  = 0
      IF(SYM .EQ. 2) THEN
         ICNTL(12)  = 0
      ELSE
         ICNTL(12)  = 1
      ENDIF
      ICNTL(13) = 0
      IF (SYM.eq.1.AND.NPROCS.EQ.1) THEN 
        ICNTL(14) = 5  ! it should work with 0
      ELSE IF (NPROCS .GT. 4) THEN
        ICNTL(14) = 30
      ELSE
        ICNTL(14) = 20
      END IF
      ICNTL(15) = 0
      ICNTL(16) = 0
      ICNTL(17) = 0
      ICNTL(18) = 0
      ICNTL(19) = 0
      ICNTL(20) = 0
      ICNTL(21) = 0
      ICNTL(22) = 0
      ICNTL(23) = 0
      ICNTL(24) = 0
      ICNTL(39) = 1
      ICNTL(40)  = 0 ! Second level parallelism on
      KEEP(12) = 0
      KEEP(11) = 2147483646
      KEEP(24) = 18
      KEEP(68) = 0
      KEEP(36) = 1
      KEEP(1) = 8
      KEEP(7)  = 150
      KEEP(8)  = 120
      KEEP(57) = 500
      KEEP(58) = 250
      KEEP(59) = 120
#if defined(t3e) || defined(sgi)
      KEEP(9)  = 400
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(85) =  160
        KEEP(62) =  50
        IF (NPROCS.GE.128) KEEP(62)=200
        IF (NPROCS.GE.128) KEEP(9)=500
        IF (NPROCS.GE.256) KEEP(9)=600
      ELSE
        KEEP(4)  = 24
        KEEP(3) =  96
        KEEP(5) =  16
        KEEP(6) =  48
        KEEP(85) =  80
        KEEP(62) = 100
        IF (NPROCS.GE.128) KEEP(62)=150
        IF (NPROCS.GE.64) KEEP(9)=500
        IF (NPROCS.GE.128) KEEP(9)=600
        IF (NPROCS.GE.256) KEEP(9)=700
      END IF
#elif defined(SP_)
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 800
        KEEP(85) = 400
      ELSE
        KEEP(4)  = 24
        KEEP(3) =  96
        KEEP(5) =  16
        KEEP(6) =  48
        KEEP(9)  = 400
        KEEP(85) = 100
      END IF
      KEEP(62) = 150
#else
      IF ( SYM .eq. 0 ) THEN
        KEEP(4)  = 32
        KEEP(3)  = 96
        KEEP(5)  = 16
        KEEP(6)  = 32
        KEEP(9)  = 700
        KEEP(85) =  300
        KEEP(62) =  50
        IF (NPROCS.GE.128) KEEP(62)=200
        IF (NPROCS.GE.128) KEEP(9)=800
        IF (NPROCS.GE.256) KEEP(9)=900
      ELSE
        KEEP(4)  = 24 
        KEEP(3)  = 96  
        KEEP(5)  = 16
        KEEP(6)  = 48
        KEEP(9)  = 400
        KEEP(85) = 100
        KEEP(62) = 100
        IF (NPROCS.GE.128) KEEP(62)=150
        IF (NPROCS.GE.64) KEEP(9)=800
        IF (NPROCS.GE.128) KEEP(9)=900
      END IF
#endif
      KEEP(63) = 60
      KEEP(48) = 5
      KEEP(17) = 0
      CALL CMUMPS_464( KEEP(34), KEEP(35))
      KEEP(37) = MAX(800, 2*NPROCS)
      IF ( NPROCS > 256 ) THEN
        KEEP(39) = 10000
      ELSEIF ( NPROCS > 128 ) THEN
        KEEP(39) = 20000
      ELSEIF ( NPROCS > 64 ) THEN
        KEEP(39) = 40000
      ELSEIF ( NPROCS > 16 ) THEN
        KEEP(39) = 80000
      ELSE
        KEEP(39) = 160000
      END IF
      KEEP(40) = -1 - 456789
      KEEP(45) = 0
      KEEP(47) = 2
#if defined(SP_)
      KEEP( 51 )  = 70
#elif defined(valuesV416sgi)
      KEEP( 51 )  = 48
#elif defined(t3e) || defined(sgi)
      KEEP( 51 )  = 24
#else
      KEEP( 51 )  = 48
#endif
      KEEP(64) = 10
      KEEP(69) = 4
      KEEP(75) = 1
      KEEP(76) = 2
      KEEP(77) = 30
      IF (NPROCS.GT.4) THEN
          KEEP(78)=MAX(
     &       int(log(REAL(NPROCS))/log(REAL(2))) - 2 
     &       , 0         )
      ENDIF
      KEEP(210) = 2 ! extra splitting is off (default for 4.6.3)
      KEEP(79) = -10
      KEEP(80) = 1
      KEEP(81) = 0
      KEEP(82) = 5
      KEEP(83) = MIN(8,NPROCS/4)
      KEEP(83) = MAX(MIN(4,NPROCS),MAX(KEEP(83),1))
      KEEP(84) = -8
      KEEP(86)=1
      KEEP(87)=0
      KEEP(88)=0
      KEEP(90)=1
      KEEP(91)=MIN(8, NPROCS)
      KEEP(91) = MAX(MIN(4,NPROCS),MIN(KEEP(83),KEEP(91)))
      KEEP(207)=1
      IF(NPROCS.LT.48)THEN
         KEEP(102)=150
      ELSEIF(NPROCS.LT.128)THEN
         KEEP(102)=150
      ELSEIF(NPROCS.LT.256)THEN
         KEEP(102)=200
      ELSEIF(NPROCS.LT.512)THEN
         KEEP(102)=300
      ELSEIF(NPROCS.GT.512)THEN
         KEEP(102)=400
      ENDIF
      KEEP(99)=0
      KEEP(100)=0
      KEEP(204)=0
      KEEP(205)=0
      KEEP(209)=-1 ! Size /strategy for MAXS for solve
      KEEP(104) = 16
      KEEP(107)=0
      KEEP(213)=201
      KEEP(217)=0
      KEEP(215)=0
      KEEP(216)=1
      KEEP(218)=50
      KEEP(219)=1
      RETURN
      END SUBROUTINE CMUMPS_20
      SUBROUTINE CMUMPS_195(N, NZ, IRN, ICN, LIW, IKEEP, PTRAR,
     *            IORD, NFSIZ, FILS, FRERE, LISTVAR_SCHUR, SIZE_SCHUR,
     *            ICNTL, INFO, KEEP,KEEP8, NSLAVES, SYM, PIV, id)
      USE CMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N,NZ,LIW,IORD,SIZE_SCHUR, NSLAVES, SYM
      INTEGER PTRAR(N,4), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER IKEEP(N,3)
      INTEGER  LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER INFO(40), ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      TYPE (CMUMPS_STRUC) :: id      
      INTEGER IRN(NZ), ICN(NZ)  
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER IERR
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW, IN, IFSON
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG, ITEMP, symmetry
      INTEGER MedDens, NBQD, AvgDens
      LOGICAL PROK, COMPRESS_SCHUR
      INTEGER NBBUCK
      INTEGER, DIMENSION(:), ALLOCATABLE :: HEAD
      INTEGER NUMFLAG
      INTEGER OPT_METIS_SIZE
      INTEGER, DIMENSION(:), ALLOCATABLE :: OPTIONS_METIS
      INTEGER THRESH, IVersion
      LOGICAL AGG6
      INTEGER MINSYM
      PARAMETER (MINSYM=50)
      INTEGER K79REF
      PARAMETER(K79REF=12000000)
      INTEGER PIV(N)
      INTEGER MTRANS, COMPRESS,NCMP,IERROR,J,JPERM,NCST,NLOCKED,P11,
     *     K1,K2,P1,P2
      INTEGER TOTEL
      LOGICAL IDENT,V1,V2, REBUILTG, SPLITROOT
      EXTERNAL CMUMPS_197, CMUMPS_198, 
     *         CMUMPS_199, CMUMPS_351,
     *         CMUMPS_557, CMUMPS_201
#if defined(OLDDFS)
      EXTERNAL CMUMPS_200
#endif
      EXTERNAL CMUMPS_623
      EXTERNAL CMUMPS_547, CMUMPS_550,
     *     CMUMPS_565, CMUMPS_566,
     *     CMUMPS_549,CMUMPS_556
#if(SP_DEBUG >= 3)
      INTERFACE
         SUBROUTINE CMUMPS_547(
     *        N,NZ, IRN, ICN, PERM,
     *        NCMP, IW, LW, IPE, LEN, IQ, 
     *        FLAG, ICMP, IWFR,
     *        IERROR, KEEP,KEEP8, ICNTL)
         INTEGER N,NZ,NCMP,LW,IWFR,IERROR
         INTEGER ICNTL(40),KEEP(500)
         INTEGER*8 KEEP8(150)
         INTEGER IRN(NZ),ICN(NZ),IW(LW),PERM(N),IPE(N+1)
         INTEGER LEN(N),IQ(N),FLAG(N),ICMP(N) 
         END SUBROUTINE CMUMPS_547
         SUBROUTINE CMUMPS_550(N,NCMP,N11,N22,PIV,
     *        PERM,IPERM)
         IMPLICIT NONE
         INTEGER N11,N22,N,NCMP
         INTEGER PIV(N),PERM(N),IPERM(N)
         END SUBROUTINE CMUMPS_550
         SUBROUTINE CMUMPS_565(N,NCMP,N11,N22,PIV,
     *        PE,NV,PE_EXP,NV_EXP)
         INTEGER N11,N22,N,NCMP
         INTEGER PIV(N),NV(N),PE(N),NV_EXP(N),PE_EXP(N)
         END SUBROUTINE CMUMPS_565
         SUBROUTINE CMUMPS_566(N,LIW,IPE,IW,N11,N22,PIV,
     *        PE,NV,WORK)
         INTEGER N11,N22,N,LIW
         INTEGER PIV(N),NV(N),PE(N),IPE(N+1),IW(LIW),WORK(N)
         END SUBROUTINE CMUMPS_566
         SUBROUTINE CMUMPS_549(N,PE,PERM,NFILS,WORK)
         INTEGER N
         INTEGER PE(N),PERM(N),NFILS(N),WORK(N)
         END SUBROUTINE CMUMPS_549
         SUBROUTINE CMUMPS_556(
     *        N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     *        NCST,KEEP,KEEP8,id)
         USE CMUMPS_STRUC_DEF
         TYPE (CMUMPS_STRUC) :: id
         INTEGER N,NCST
         INTEGER PIV(N),FRERE(N),FILS(N),NFSIZ(N),IKEEP(N,3)
         INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
         END SUBROUTINE CMUMPS_556
      END INTERFACE
#endif
      ALLOCATE( IW ( LIW ), stat = IERR )
      IF ( IERR .GT. 0 ) THEN
         INFO( 1 ) = -7
         INFO( 2 ) = LIW
         RETURN
      ENDIF
      LLIW = LIW - 2*N - 1
      L1 = LLIW + 1
      L2 = L1 + N
      MPRINT= ICNTL(3)
      PROK  = (MPRINT.GT.0)
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      LDIAG = ICNTL(4)
      REBUILTG = .FALSE.
      COMPRESS_SCHUR = .FALSE.
      IF (KEEP(1).LT.0) KEEP(1) = 0
      NEMIN = KEEP(1)
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       WRITE (MP,99999) N, NZ, LIW, INFO(1)
       K = MIN0(10,NZ)
       IF (LDIAG.EQ.4) K = NZ
       IF (K.GT.0) WRITE (MP,99998) (IRN(I),ICN(I),I=1,K)
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       ENDIF
      ENDIF
         NCMP    = N   
         IF (KEEP(60).NE.0) THEN
           IF ((SIZE_SCHUR.LE.0 ).OR.
     &          (SIZE_SCHUR.GE.N) ) GOTO 90
         ENDIF
#if defined(metis)
          IF  ( ( KEEP(60).NE.0).AND.(SIZE_SCHUR.GT.0)
     &         .AND.
     &          ((IORD.EQ.7).OR.(IORD.EQ.5))
     &        )THEN
          COMPRESS_SCHUR=.TRUE.
          NCMP          = N-SIZE_SCHUR
          CALL CMUMPS_623(N,NCMP,NZ,IRN, ICN, IW(1), LLIW, 
     *         IW(L2), PTRAR(1,2),
     *         PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *         INFO(1), INFO(2), ICNTL, symmetry, 
     *         SYM, MedDens, NBQD, AvgDens, 
     *         LISTVAR_SCHUR, SIZE_SCHUR, 
     *         FRERE,FILS)
           IORD = 5
           KEEP(95) = 1
           NBQD     = 0    ! Not clear how to use it
          ELSE
#endif
          CALL CMUMPS_351(N,NZ,IRN, ICN, IW(1), LLIW, 
     *         IW(L2), PTRAR(1,2),
     *         PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *         INFO(1), INFO(2), ICNTL, symmetry, 
     *         SYM, MedDens, NBQD, AvgDens)
#if defined(metis)
          ENDIF
#endif
          INFO(8) = symmetry
          IF(NBQD .GT. 0) THEN
             IF( KEEP(50) .EQ. 2 .AND. ICNTL(12) .LE. 1 ) THEN
                IF(KEEP(95) .NE. 1) THEN
                   IF ( PROK ) 
     *                  WRITE( MPRINT,*) 
     *                  'Compressed/constrained ordering set OFF'
                   KEEP(95) = 1   
                ENDIF
             ENDIF
          ENDIF
      IF ( (KEEP(60).NE.0) .AND. (IORD.GT.1) .AND.
     &     .NOT. COMPRESS_SCHUR ) THEN
            IORD = 0    ! AMD based ordering is used
      ENDIF 
      IF (IORD.EQ.7) 
     *     CALL CMUMPS_423( N, SYM, NSLAVES, IORD, 
     *                       symmetry,MedDens, NBQD, AvgDens)
#if ! defined(metis)
      IF ( IORD .EQ. 5 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: METIS not available. Ordering set to default.'
        CALL CMUMPS_423( N, SYM, NSLAVES, IORD,
     *                       symmetry,MedDens, NBQD, AvgDens)
      END IF
#endif
#if ! defined(pord)
      IF ( IORD .EQ. 4 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: PORD not available. Ordering set to default.'
        CALL CMUMPS_423( N, SYM, NSLAVES, IORD,
     *                       symmetry,MedDens, NBQD, AvgDens)
      END IF
#endif
#if ! defined(scotch)
      IF ( IORD .EQ. 3 ) THEN
        IF (PROK) WRITE(MP,*)
     *  'WARNING: SCOTCH not available. Ordering set to default.'
        CALL CMUMPS_423( N, SYM, NSLAVES, IORD,
     *                       symmetry,MedDens, NBQD, AvgDens)
      END IF
#endif
      IF(SYM .EQ. 2) THEN
         IF(KEEP(95) .EQ. 3 .AND. IORD .NE. 2) THEN
            IF (PROK) WRITE(MP,*)
     *       'WARNING: ANAL_F constrained ordering not available ', 
     *       'with selected ordering'
            KEEP(95) = 2
         ENDIF
         IF(KEEP(95) .EQ. 2 .AND. IORD .EQ. 0) THEN 
            IF (PROK) WRITE(MP,*)
     *           'WARNING: ANAL_F AMD not available with ', 
     *           'compressed ordering -> move to QAMD'
            IORD = 6
         ENDIF
      ELSE
         KEEP(95) = 1
      ENDIF
      MTRANS = KEEP(23)
      COMPRESS = KEEP(95) - 1
      IF(COMPRESS .GT. 0 .AND. KEEP(52) .EQ. -2) THEN
         IF(id%CNTL(4) .GE. 0.0E0) THEN
            IF (KEEP(1).LE.8) THEN
             KEEP(1) = 16   
            ELSE
             KEEP(1) = 2*KEEP(1)
            ENDIF
            IF (PROK) 
     *       WRITE(MP,*) 'Setting static pivoting ON, COMPRESS =',
     *           COMPRESS
         ENDIF
      ENDIF
      IF(MTRANS .GT. 0 .AND. KEEP(50) .EQ. 2) THEN
         KEEP(23) = 0
      ENDIF
#if(SP_DEBUG >= 2)
      WRITE(*,*) 'ICNTL(12),MTRANS,COMPRESS',KEEP(95),MTRANS,COMPRESS
#endif
      IF ((IORD.NE.1) .AND. (IORD.NE.5)) THEN
          IF(COMPRESS .EQ. 2) THEN
#if(SP_DEBUG >= 2)
             WRITE(*,*) 'before 2x2 splitting,KEEP(93),KEEP(94)',
     *            KEEP(93),KEEP(94)
#endif
             CALL  CMUMPS_556(
     *            N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     *            NCST,KEEP,KEEP8,id)
#if(SP_DEBUG >= 2)
             WRITE(*,*) 'after 2x2 splitting,KEEP(93),NCST,KEEP(94)',
     *            KEEP(93),NCST,KEEP(94)
#endif
          ENDIF
 100      CONTINUE
          IF(COMPRESS .GE. 1) THEN
             CALL CMUMPS_547(
     *            N,NZ, IRN, ICN, PIV,
     *            NCMP, IW(1), LLIW, IW(L2),PTRAR(1,2),PTRAR, 
     *            IW(L1), FILS, IWFR,
     *            IERROR, KEEP,KEEP8, ICNTL)
             symmetry = 1.0E0
          ELSE IF (REBUILTG) THEN
             CALL CMUMPS_351
     *            (N,NZ,IRN, ICN, IW(1), LLIW, IW(L2), PTRAR(1,2), 
     *            PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114), 
     *            INFO(1), INFO(2), ICNTL, symmetry, SYM,
     *            MedDens, NBQD, AvgDens)
             INFO(8) = symmetry
             NCMP = N
          ENDIF
          IF ( (symmetry.LT.MINSYM).AND.(SYM.EQ.0) ) THEN
             IF(KEEP(23) .EQ. 7 ) THEN
                KEEP(23) = -5
                DEALLOCATE (IW)
                RETURN
             ELSE IF(KEEP(23) .EQ. -9876543) THEN
                IDENT = .TRUE.
                KEEP(23) = 5
                IF (PROK) WRITE(MPRINT,'(A)')
     *               ' ... Apply column permutation (already computed)'
                DO J=1,N
                   JPERM = PIV(J)
                   FILS(JPERM) = J
                   IF (JPERM.NE.J) IDENT = .FALSE.
                ENDDO
                IF (.NOT.IDENT) THEN
                   DO K=1,NZ
                      J = ICN(K)
                      IF ((J.LE.0).OR.(J.GT.N)) CYCLE
                      ICN(K) = FILS(J)
                   ENDDO
                   IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *                  WRITE(MP,'(/A)')
     *                  ' WARNING input matrix data modified'
                   REBUILTG = .TRUE.
                   GOTO 100
                ELSE
                   KEEP(23) = 0
                ENDIF
             ENDIF
          ELSE IF (KEEP(23) .EQ. 7 .OR. KEEP(23) .EQ. -9876543 ) THEN
             IF (PROK) WRITE(MPRINT,'(A)')
     *               ' ... No column permutation'
             KEEP(23) = 0
          ENDIF
          IF (PROK) THEN
             IF (IORD.EQ.2) THEN
                WRITE(MPRINT,'(A)') ' Ordering based on AMF '
#if defined(scotch)
             ELSE IF (IORD.EQ.3) THEN
                WRITE(MPRINT,'(A)') ' Ordering based on SCOTCH '
#endif
#if defined(pord)
             ELSE IF (IORD.EQ.4) THEN
                WRITE(MPRINT,'(A)') ' Ordering based on PORD '
#endif
             ELSE IF (IORD.EQ.6) THEN
                WRITE(MPRINT,'(A)') ' Ordering based on QAMD '
             ELSE
                WRITE(MPRINT,'(A)') ' Ordering based on AMD '
             ENDIF
          ENDIF
       IF ( KEEP(60) .NE. 0 ) THEN
         CALL CMUMPS_162(N, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1), 
     *   IW(L1), IKEEP, 
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3),
     *   LISTVAR_SCHUR, SIZE_SCHUR)
         IF (KEEP(60)==1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSE
           KEEP(38) = LISTVAR_SCHUR(1)
         ENDIF
       ELSE
      IF ( .FALSE. ) THEN
#if defined(pord)
      ELSEIF (IORD .EQ. 4) THEN
         IF(COMPRESS .EQ. 1) THEN
            DO I=L1,L1-1+KEEP(93)/2
               IW(I) = 2
            ENDDO
            DO I=L1+KEEP(93)/2,L1+NCMP-1
               IW(I) = 1
            ENDDO
            CALL CMUMPS_PORDF_WND(NCMP, IWFR-1, IW(L2), IW, 
     *           IW(L1), NCMPA, N)
         ELSE
            CALL CMUMPS_PORDF(NCMP, IWFR-1, IW(L2), IW(1), 
     *           IW(L1), NCMPA)
         ENDIF
         IF ( NCMPA .NE. 0 ) THEN
            write(6,*) ' Out PORD, NCMPA=', NCMPA
            INFO( 1 ) = -9999
            INFO( 2 ) = 4
            RETURN
         ENDIF
#endif    
#if defined(scotch)
      ELSEIF (IORD .EQ. 3) THEN
         CALL CMUMPS_SCOTCH(NCMP, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1), 
     *        IW(L1), IKEEP, 
     *        IKEEP(1,2), NCMPA)
        IF ( NCMPA .NE. 0 ) THEN
          write(6,*) ' Out SCTOCH, NCMPA=', NCMPA
          INFO( 1 ) = -9999
          INFO( 2 ) = 3
          RETURN
        ENDIF
#endif
      ELSEIF (IORD .EQ. 2) THEN
         NBBUCK = 2*N
         ALLOCATE( HEAD ( 0: NBBUCK + 1), stat = IERR )
         IF ( IERR .GT. 0 ) THEN
            INFO( 1 ) = -7
            INFO( 2 ) = NBBUCK+2
            RETURN
         ENDIF
         IF(COMPRESS .GE. 1) THEN
            DO I=L1,L1-1+KEEP(93)/2
               IW(I) = 2
            ENDDO
            DO I=L1+KEEP(93)/2,L1+NCMP-1
               IW(I) = 1
            ENDDO
         ELSE
            IW(L1) = -1
         ENDIF
         IF(COMPRESS .LE. 1) THEN
            CALL CMUMPS_337(NCMP, NBBUCK, LLIW, IW(L2),
     *           IWFR, PTRAR(1,2),
     *           IW(1), IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS,
     *           IKEEP(1,3), PTRAR, PTRAR(1,3), HEAD)
         ELSE
            IF(PROK) WRITE(MPRINT,'(A)') 
     *           ' Constrained Ordering based on AMF'
            CALL CMUMPS_560(NCMP, NBBUCK, LLIW, IW(L2),
     *           IWFR, PTRAR(1,2), 
     *           IW(1), IW(L1), IKEEP, IKEEP(1,2), NCMPA, FILS, 
     *           IKEEP(1,3), PTRAR, PTRAR(1,3), HEAD,
     *           NFSIZ, FRERE)
         ENDIF
         DEALLOCATE(HEAD)
      ELSEIF (IORD .EQ. 6) THEN
        ALLOCATE( HEAD ( N ), stat = IERR )
        IF ( IERR .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = N
          RETURN
        ENDIF
         THRESH = 1
         IVersion = 2
         IF(COMPRESS .EQ. 1) THEN
            DO I=L1,L1-1+KEEP(93)/2
               IW(I) = 2
            ENDDO
            DO I=L1+KEEP(93)/2,L1+NCMP-1
               IW(I) = 1
            ENDDO
            TOTEL = KEEP(93)+KEEP(94)
         ELSE
            IW(L1) = -1
            TOTEL = N
         ENDIF
         CALL CMUMPS_421(TOTEL,IVersion, THRESH, HEAD,
     *        NCMP, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1),
     *        IW(L1), IKEEP,
     *        IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3))
        DEALLOCATE(HEAD)
      ELSE
         CALL CMUMPS_197(NCMP, LLIW, IW(L2), IWFR, PTRAR(1,2), IW(1), 
     *        IW(L1), IKEEP, 
     *        IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3))
      ENDIF
      ENDIF
      IF(COMPRESS .GE. 1) THEN
         CALL CMUMPS_548(NCMP,IW(L2),IW(L1),FILS)
#if(SP_DEBUG >= 2)
         WRITE(*,*) 'BEFORE CMUMPS_549'
#endif
         CALL CMUMPS_549(NCMP,IW(L2),IKEEP(1,1),
     *        FRERE,PTRAR(1,1))
#if(SP_DEBUG >= 2)
         WRITE(*,*) 'AFTER CMUMPS_549'
#endif
         DO I=1,NCMP
            IKEEP(IKEEP(I,1),2)=I
         ENDDO
#if(SP_DEBUG >= 2)
#endif
         CALL CMUMPS_550(N,NCMP,KEEP(94),KEEP(93),
     *        PIV,IKEEP(1,1),IKEEP(1,2))
#if(SP_DEBUG >= 2)
#endif
         COMPRESS = -1
      ENDIF
      ENDIF
      IF ((IORD.EQ.1) .OR. (IORD.EQ.5)) THEN
#if defined(metis)
         IF (IORD.EQ.5) THEN
            IF (PROK) THEN
               WRITE(MPRINT,'(A)') ' Ordering based on METIS '
            ENDIF
 200        CONTINUE
            IF(COMPRESS .EQ. 1) THEN
               CALL CMUMPS_547(
     *              N,NZ, IRN, ICN, PIV,
     *              NCMP, IW(1), LLIW, IW(L2),PTRAR(1,2),PTRAR, 
     *              IW(L1), FILS, IWFR,
     *              IERROR, KEEP,KEEP8, ICNTL)
               symmetry = 1.0E0
            ELSE IF (REBUILTG) THEN
               CALL CMUMPS_351
     *              (N,NZ,IRN, ICN, IW(1), LLIW, IW(L2), PTRAR(1,2), 
     *              PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114), 
     *              INFO(1), INFO(2), ICNTL, symmetry, SYM,
     *              MedDens, NBQD, AvgDens)
               INFO(8) = symmetry
               NCMP = N
            ENDIF
            IF ( (symmetry.LT.MINSYM).AND.(SYM.EQ.0) ) THEN
               IF(KEEP(23) .EQ. 7 ) THEN
                  KEEP(23) = -5
                  DEALLOCATE (IW)
                  RETURN
               ELSE IF(KEEP(23) .EQ. -9876543) THEN
                  IDENT = .TRUE.
                  KEEP(23) = 5
                  IF (PROK) WRITE(MPRINT,'(A)')
     *        ' ... Apply column permutation (already computed)'
                  DO J=1,N
                     JPERM = PIV(J)
                     FILS(JPERM) = J
                     IF (JPERM.NE.J) IDENT = .FALSE.
                  ENDDO
                  IF (.NOT.IDENT) THEN
                     DO K=1,NZ
                        J = ICN(K)
                        IF ((J.LE.0).OR.(J.GT.N)) CYCLE
                        ICN(K) = FILS(J)
                     ENDDO
                     IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *                    WRITE(MP,'(/A)')
     *                    ' WARNING input matrix data modified'
                     REBUILTG = .TRUE.
                     GOTO 200
                  ELSE
                     KEEP(23) = 0
                  ENDIF
               ENDIF
            ELSE IF (KEEP(23) .EQ. 7 .OR. KEEP(23) .EQ. -9876543 ) THEN
               IF (PROK) WRITE(MPRINT,'(A)')
     *              ' ... No column permutation'
               KEEP(23) = 0
            ENDIF
            NUMFLAG = 1
            OPT_METIS_SIZE = 8
            ALLOCATE( OPTIONS_METIS (OPT_METIS_SIZE ), stat = IERR )
            IF ( IERR .GT. 0 ) THEN
               INFO( 1 ) = -7
               INFO( 2 ) = OPT_METIS_SIZE
               RETURN
            ENDIF
            OPTIONS_METIS(1) = 0
            IF (COMPRESS .EQ. 1) THEN
               DO I=1,KEEP(93)/2
                  FILS(I) = 2
               ENDDO
               DO I=KEEP(93)/2+1,NCMP
                  FILS(I) = 1
               ENDDO
               CALL METIS_NODEWND(NCMP, IW(L2), IW(1),FILS, 
     *              NUMFLAG, OPTIONS_METIS,
     *              IKEEP(1,2), IKEEP(1,1) )
            ELSE
               CALL METIS_NODEND(NCMP, IW(L2), IW(1), NUMFLAG, 
     *              OPTIONS_METIS,
     *              IKEEP(1,2), IKEEP(1,1) )
            ENDIF
            DEALLOCATE (OPTIONS_METIS)
            IF ( COMPRESS_SCHUR ) THEN
               CALL CMUMPS_622(
     *         N, NCMP, IKEEP(1,1),IKEEP(1,2),
     *         LISTVAR_SCHUR, SIZE_SCHUR, FILS)
            ENDIF
            IF(COMPRESS .EQ. 1) THEN
#if(SP_DEBUG >= 5)
               WRITE(*,*) 'IKEEP = ',IKEEP(1:NCMP,1)
               WRITE(*,*) 'PIV',PIV
#endif
               CALL CMUMPS_550(N,NCMP,KEEP(94),
     *              KEEP(93),PIV,IKEEP(1,1),IKEEP(1,2))
#if(SP_DEBUG >= 5)
               WRITE(*,*) 'IKEEP = ',IKEEP(1:N,1)
#endif
            ENDIF
         ELSE IF (IORD.NE.1) THEN
            WRITE(*,*) IORD
            WRITE(*,*)  'bad option for ordering'
            GO TO 40
         ENDIF
#endif
         IF (PROK) THEN
            IF (IORD.EQ.1) THEN
               WRITE(MPRINT,'(A)') ' Ordering given is used'
            ENDIF
         ENDIF
         DO K=1,N
            PTRAR(K,1) = 0
         ENDDO
         DO K=1,N
            IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     *           GO TO 40
            IF (PTRAR(IKEEP(K,1),1).EQ.1) THEN
               GOTO 40
            ELSE
               PTRAR(IKEEP(K,1),1) = 1
            ENDIF
         ENDDO
      ENDIF
       IF (IORD.EQ.1 .OR. IORD.EQ.5 .OR. COMPRESS.EQ.-1) THEN
          COMPRESS = 0
        IF (KEEP(106)==1) THEN
        IF (IORD.NE.1) THEN ! GNEW already done
          CALL CMUMPS_351(N,NZ,IRN, ICN, IW(1), LLIW, 
     *         IW(L2), PTRAR(1,2),
     *         PTRAR, IW(L1), IWFR, KEEP(113), KEEP(114),
     *         INFO(1), INFO(2), ICNTL, symmetry, SYM,
     *         MedDens, NBQD, AvgDens)
          INFO(8) = symmetry
        ENDIF
          ALLOCATE( HEAD ( 2*N ), stat = IERR )
          IF ( IERR .GT. 0 ) THEN
             INFO( 1 ) = -7
             INFO( 2 ) = 2*N
             RETURN
          ENDIF
          THRESH = -1
          IF (KEEP(60) == 0) THEN
             ITEMP = 0 
          ELSE
             ITEMP = SIZE_SCHUR
             IF (KEEP(60)==1) THEN
                KEEP(20) = LISTVAR_SCHUR(1)
             ELSE
                KEEP(38) = LISTVAR_SCHUR(1)
             ENDIF
          ENDIF
          AGG6 = ( NSLAVES == 1 )
          CALL CMUMPS_422(THRESH, HEAD,
     *         N, LLIW, IW(L2), IWFR, PTRAR(1,2), IW,
     *         IW(L1), HEAD(N+1),
     *         IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), PTRAR, PTRAR(1,3), 
     *         IKEEP(1,1), LISTVAR_SCHUR, ITEMP, AGG6)
          DEALLOCATE(HEAD)
        ELSE
          CALL CMUMPS_198(N, NZ, IRN, ICN, IKEEP, IW(1),
     *         LLIW, IW(L2),
     *         PTRAR(1,2), IW(L1), IWFR,
     *         INFO(1),INFO(2), KEEP(11), MP)
          IF (KEEP(60) .EQ. 0) THEN
             ITEMP = 0 
             CALL CMUMPS_199(N, IW(L2), IW, LLIW, IWFR, IKEEP,
     *            IKEEP(1,2), IW(L1),
     *            PTRAR, NCMPA, ITEMP)
          ELSE
             CALL CMUMPS_199(N, IW(L2), IW, LLIW, IWFR, IKEEP,
     *            IKEEP(1,2), IW(L1),
     *            PTRAR, NCMPA, SIZE_SCHUR)
             IF (KEEP(60) .EQ. 1) THEN
                KEEP(20) = LISTVAR_SCHUR(1)
             ELSE
                KEEP(38) = LISTVAR_SCHUR(1)
             ENDIF
          ENDIF
        ENDIF
       ENDIF
#if defined(OLDDFS)
      CALL CMUMPS_200
     *    (N, IW(L2), IW(L1), IKEEP(1,1), IKEEP(1,2), IKEEP(1,3),
     *     NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, KEEP(60))
#else
      CALL CMUMPS_557
     *    (N, IW(L2), IW(L1), IKEEP(1,1), IKEEP(1,2), IKEEP(1,3),
     *     NFSIZ, PTRAR, INFO(6), FILS, FRERE, 
     *     PTRAR(1,3), NEMIN, PTRAR(1,4), KEEP(60),
     *     KEEP(20),KEEP(38),PTRAR(1,2),KEEP(104),IW(1),KEEP(50))
#endif
      IF (KEEP(60).NE.0)  THEN
         IF (KEEP(60)==1) THEN
           IN = KEEP(20)
         ELSE
           IN = KEEP(38)
         ENDIF
         DO WHILE (IN.GT.0) 
          IN = FILS (IN)
         END DO
         IFSON = -IN
         IF (KEEP(60)==1) THEN
           IN = KEEP(20)
         ELSE
           IN = KEEP(38)
         ENDIF
         DO I=2,SIZE_SCHUR
          FILS(IN) = LISTVAR_SCHUR (I)
          IN       = FILS(IN)
          FRERE (IN) = N+1
         ENDDO
         FILS(IN) = -IFSON
      ENDIF
      CALL CMUMPS_201(IKEEP(1,2),
     *  PTRAR(1,3), INFO(6),
     *  INFO(5), KEEP(2), KEEP(50), INFO(20), 
     *  KEEP(101),KEEP8(1),KEEP(108))
      IF ( KEEP(53) .NE. 0 ) THEN
        CALL CMUMPS_209( N, FRERE, FILS, NFSIZ, KEEP(20) )
      END IF
      IF (  (KEEP(48) == 4 .AND. KEEP(10).GT.0)
     *     .OR.
     *      (KEEP (48)==5 .AND. KEEP(10) .GT. 0 )
     *     .OR.
     *      (KEEP(24).NE.0.AND.KEEP(10).GT.0) ) THEN 
          CALL CMUMPS_510(KEEP(10), KEEP(2),
     *    KEEP(48), KEEP(50), NSLAVES)
      END IF
      IF (KEEP(210).LT.0.OR.KEEP(210).GT.2) KEEP(210)=0
      IF (KEEP(210).EQ.0.AND.KEEP(201).NE.0) KEEP(210)=1 !On
      IF (KEEP(210).EQ.0.AND.KEEP(201).EQ.0) KEEP(210)=2 !Off
      IF (KEEP(210).EQ.2) KEEP(79)=HUGE(KEEP(79))
      IF (KEEP(210).EQ.1.AND.KEEP(79).LE.0) THEN
        IF ( HUGE(KEEP(79)) / K79REF + 1 .GE. NSLAVES ) THEN
        KEEP(79)=HUGE(KEEP(79))
        ELSE
        KEEP(79)=K79REF * NSLAVES
        ENDIF
      ENDIF
      IF (KEEP(210).EQ.1) THEN
       SPLITROOT = .FALSE. 
       IF ( KEEP(62).GE.1) 
     *  CALL CMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     *                       NSLAVES, KEEP,KEEP8, SPLITROOT,
     *                       MP, LDIAG)
      ENDIF
      SPLITROOT = ((ICNTL(13).GT.0 .AND. NSLAVES.GT.ICNTL(13)) .OR.
     &              ICNTL(13).EQ.-1 )
     &              .AND. (KEEP(60).EQ.0)
      IF (SPLITROOT)
     *   CALL CMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     *                    NSLAVES, KEEP,KEEP8, SPLITROOT,
     *                    MP, LDIAG)
      IF ( KEEP(24) == 0 ) THEN
         CALL CMUMPS_202(N, NZ, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2),
     *      IRN, ICN, PTRAR(1,3), PTRAR(1,4), PTRAR,
     *      ICNTL, INFO, KEEP,KEEP8)
      ELSE
        CALL CMUMPS_348(N, NZ, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2), 
     *      IRN, ICN, PTRAR, ICNTL, INFO, KEEP,KEEP8)
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99986) (PTRAR(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99985) (PTRAR(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99984) (PTRAR(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
   40 INFO(1) = -4
      INFO(2) = K
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
      GOTO 90
   45 INFO(1) = -32    ! Schur size is out of range
      INFO(2) = SIZE_SCHUR
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
   90 CONTINUE 
      DEALLOCATE(IW)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     * '                N         NZ         LIW       INFO(1)'/,
     * 9X, I8, I11, I12, I14)
99998 FORMAT ('Matrix entries:    IRN()   ICN()'/
     * (I12, I7, I12, I7, I12, I7))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99986 FORMAT ('PTRAR(.,1)=', 10I6/(12X, 10I6))
99985 FORMAT ('PTRAR(.,2)=', 10I6/(12X, 10I6))
99984 FORMAT ('PTRAR(.,3)=', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE CMUMPS_195
#if defined(NOT_USED)
      SUBROUTINE CMUMPS_196(N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     *  IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM)
      INTEGER N,NZ,LW,IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER symmetry, SYM
      INTEGER ICNTL(40)
      INTEGER  IRN(NZ), ICN(NZ) 
      INTEGER LEN(N), IQ(N)
      INTEGER IPE(N)
      INTEGER FLAG(N), IW(LW)
      INTEGER MP, MPG
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR
      INTEGER NZOFFA, NDIAGA
      REAL RSYM
      INTRINSIC NINT
      MP = ICNTL(2)
      MPG = ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      DO 10 I=1,N
        IPE(I) = 0
   10 CONTINUE
      DO 50 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
   50 CONTINUE
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (MOD(IFLAG,2).EQ.0) IFLAG  = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (MOD(K,10).GT.3 .OR. MOD(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (MOD(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = MAX(IPE(N)+IQ(N)-1,IQ(N))
      DO 115 I = 1,N
         FLAG(I) = 0
         IPE(I)  = IQ(I)
  115 CONTINUE
      DO 130 K=1,LAST
        IW(K) = 0
  130 CONTINUE
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IPE(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         K1 = IPE(I) 
         IF (K1.EQ.0) GO TO 280
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IF (SYM.EQ.0) THEN
      RSYM =  REAL(NDIAGA+2*NZOFFA - (IWFR-1))/
     *            REAL(NZOFFA+NDIAGA) 
      symmetry = NINT (100.0*RSYM)
      ELSE
       symmetry = 100
      ENDIF
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE CMUMPS_196
#endif
      SUBROUTINE CMUMPS_197(N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER NEXT(N), LEN(N), 
     $        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     $        W(N)
      INTEGER IW(IWLEN), NV(N), PE(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER*8 HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG1)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
#if defined (NOAGG1)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = MIN (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_197
      SUBROUTINE CMUMPS_199(N,IPE,IW, LW, IWFR, IPS, IPV, NV, FLAG,
     *                  NCMPA, SIZE_SCHUR)
      INTEGER N,LW,IWFR,NCMPA,SIZE_SCHUR
      INTEGER FLAG(N)
      INTEGER IPS(N), IPV(N)
      INTEGER IW(LW), NV(N), IPE(N)
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N-SIZE_SCHUR 
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL CMUMPS_194(N, IPE, IW, IP-1, LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = MIN0(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      IF (SIZE_SCHUR == 0) RETURN
      DO ML = N-SIZE_SCHUR+1,N
        ME = IPV(ML)
        IE = ME
        DO KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 160
          LN = IW(JP)
  160     IPE(IE) = -IPV(N-SIZE_SCHUR+1)
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 190
        ENDDO
  190   NV(ME) = 0
        IPE(ME) = -IPV(N-SIZE_SCHUR+1)
      ENDDO
      ME = IPV(N-SIZE_SCHUR+1)
      IPE(ME) = 0
      NV(ME) = SIZE_SCHUR
      RETURN
      END SUBROUTINE CMUMPS_199
      SUBROUTINE CMUMPS_198(N, NZ, IRN, ICN, PERM,
     * IW, LW, IPE, IQ, FLAG,
     * IWFR, IFLAG, IERROR, IOVFLO, MP)
      INTEGER N,NZ,LW,IWFR,IFLAG,IERROR
      INTEGER PERM(N)
      INTEGER IQ(N)
      INTEGER IRN(NZ), ICN(NZ) 
      INTEGER IPE(N), IW(LW), FLAG(N)
      INTEGER MP
      INTEGER IOVFLO
      INTEGER I,J,K,LBIG,L,ID,IN,LEN,JDUMMY,K1,K2
      IERROR = 0
      DO 10 I=1,N
        IQ(I) = 0
   10 CONTINUE
      DO 80 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF (I.EQ.J) GOTO 40
        IF (I.GT.J) GOTO 30
        IF (I.GE.1 .AND. J.LE.N) GO TO 60
        GO TO 50
   30   IF (J.GE.1 .AND. I.LE.N) GO TO 60
        GO TO 50
   40   IW(K) = 0
        IF (I.GE.1 .AND. I.LE.N) GO TO 80
   50   IERROR = IERROR + 1
        IW(K) = 0
        IF (IERROR.LE.1 .AND. MP.GT.0) WRITE (MP,99999) 
        IF (IERROR.LE.10 .AND. MP.GT.0) WRITE (MP,99998) K, I, J
        GO TO 80
   60   IF (PERM(J).GT.PERM(I)) GO TO 70
        IQ(J) = IQ(J) + 1
        GO TO 80
   70   IQ(I) = IQ(I) + 1
   80 CONTINUE
      IF (IERROR.GE.1) THEN
        IF (MOD(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
      ENDIF
      IWFR = 1
      LBIG = 0
      DO 100 I=1,N
        L = IQ(I)
        LBIG = MAX0(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  100 CONTINUE
      DO 140 K=1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 140
        L = K
        IW(K) = 0
        DO 130 ID=1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 110
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 120
  110     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  120     I = -IN
          IF (I.LE.0) GO TO 140
  130   CONTINUE
  140 CONTINUE
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 170 I=1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 160
        DO 150 JDUMMY=1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  150   CONTINUE
  160   IPE(J) = L
        L = L - 1
  170 CONTINUE
      IF (LBIG.GE.IOVFLO) GO TO 190
      DO 180 I=1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  180 CONTINUE
      GO TO 230
  190 IWFR = 1
      DO 220 I=1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 200
        IPE(I) = 0
        GO TO 220
  200   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 210 K=K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 210
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  210   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  220 CONTINUE
  230 RETURN
99999 FORMAT (' *** WARNING MESSAGE FROM CMUMPS_198 ***' )
99998 FORMAT (I6, ' NON-ZERO (IN ROW, I6, 11H AND COLUMN ', I6,
     * ') IGNORED')
      END SUBROUTINE CMUMPS_198
      SUBROUTINE CMUMPS_194(N, IPE, IW, LW, IWFR,NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER   IW(LW)
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END SUBROUTINE CMUMPS_194
#if defined(OLDDFS)
      SUBROUTINE CMUMPS_200(N, IPE, NV, IPS, NE, NA, NFSIZ, 
     *                  NSTEPS,
     *                  FILS, FRERE,NDD,NEMIN, KEEP60)
      INTEGER N,NSTEPS
      INTEGER NDD(N)
      INTEGER FILS(N), FRERE(N)
      INTEGER IPS(N), NE(N), NA(N), NFSIZ(N)
      INTEGER IPE(N), NV(N)
      INTEGER NEMIN, KEEP60
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,INP,IFSON,INC,INO
      INTEGER INOS,IB,IL,INT
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
      DO 20 I=1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF (NV(I).LE.0) GO TO 50
        IF = -IPE(I)
        IF (IF.NE.0) THEN
         IS = -IPS(IF)
         IF (IS.GT.0) IPE(I) = IS
         IPS(IF) = -I
        ELSE
         NR = NR - 1
         NE(NR) = I
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
      IF (NV(INB).EQ.0) THEN
       INS = INB
       GO TO 1070
      ENDIF
      IF (NV(INB).GE.NV(INS)) THEN
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
       FILS(INF) = -INB
       IPS(INF)  = -INB
       IPE(INS)  = IPE(INB)
       IPE(INB)  = INS
       INS       = INB
       GO TO 1070
      ENDIF
      INSW = INFS
 1100 INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
       INSW = INFS
       GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
      INS      =INB
      GO TO 1070
 1151 CONTINUE
      DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I  = 0
      IL = 0
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
        IF (NV(I).GT.0) GO TO 89
      IN = I
 81   IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 81
      IF = -IN
      IN = IF
 82   INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 82
      IFSON = -IN
      FILS(INL) = I
      IN = I
 83   INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 83
      IF (IFSON .EQ. I) GO TO 86
      FILS(INP) = -IFSON
      IN = IFSON
 84   INC =IN
      IN = FRERE(IN)
      IF (IN.NE.I) GO TO 84
      FRERE(INC) = FRERE(I)
      GO TO 120
 86   IF (FRERE(I).LT.0) FILS(INP) = 0
      IF (FRERE(I).GT.0) FILS(INP) = -FRERE(I)
      GO TO 120
   89   IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        NDD(IS) = NV(I)
        NFSIZ(I) = NV(I)
        IF (NA(IS).LT.1) GO TO 110
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.NDD(IS)) ) GOTO 110
        IF (NDD(IS-1)-NE(IS-1).EQ.NDD(IS)) GO TO 100
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(NDD(IS)-NDD(IS-1)+NE(IS-1)).GE.
     &    ((NDD(IS)+NE(IS-1))*
     &    (NDD(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        NDD(IS-1) = NDD(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
      IN=I
 101  INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 101
      IFSON = -IN
      IN = IFSON
 102  INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
      FILS(INL) = INO
      NFSIZ(I) = NDD(IS-1)
      IN = INO
 103  INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 103
      INOS = -IN
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
      FILS(INP) = -IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
      IF (INOS.EQ.0) FRERE(INS) = -I
      IF (INOS.NE.0) FRERE(INS) =  INOS
      IF (INOS.EQ.0) GO TO 109
 107  IN = INOS
      IF (IN.EQ.0) GO TO 109
 108  INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
 109  CONTINUE
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.LT.0) GOTO 150
        IF (IB.EQ.0) GOTO 140
        NA(IL) = 0
  140   I = IB
        GO TO 160
  150   I = -IB
        IL = IL + 1
  160 CONTINUE
      NSTEPS = IS - 1
      DO 170 I=1,N
        K = FILS(I)
        IF (K.GT.0) THEN
          FRERE(K)  = N + 1
          NFSIZ(K)  = 0
        ENDIF
 170  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_200
#else
      SUBROUTINE CMUMPS_557(N, IPE, NV, IPS, NE, NA, NFSIZ,
     *               NODE, NSTEPS,
     *               FILS, FRERE, ND, NEMIN, SUBORD, KEEP60, 
     *               KEEP20, KEEP38, NAMALG,NAMALGMAX,
     *               CUMUL,KEEP50)
      IMPLICIT NONE
      INTEGER N, NSTEPS, KEEP60, KEEP20, KEEP38, KEEP50
      INTEGER ND(N), NFSIZ(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N)
      INTEGER NEMIN,AMALG_COUNT
      INTEGER NAMALG(N),NAMALGMAX, CUMUL(N)
      DOUBLE PRECISION ACCU, FLOPS_FATHER, FLOPS_SON,
     *                  FLOPS_AVANT, FLOPS_APRES
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL,INT
      INTEGER IPERM
#if ! defined(NOAMALGTOFATHER)
      INTEGER DADI
      LOGICAL AMALG_TO_father_OK
#endif
      AMALG_COUNT = 0
      DO 10 I=1,N
        CUMUL(I)= 0
        IPS(I)  = 0
        NE(I)   = 0
        NODE(I) = 1
        SUBORD(I) = 0
        NAMALG(I) = 0
   10 CONTINUE
      FRERE(1:N) = IPE(1:N)
      NR = N + 1
      DO 50 I=1,N
        IF = -FRERE(I)
        IF (NV(I).EQ.0) THEN
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          NODE(IF) = NODE(IF)+1
        ELSE
          IF (IF.NE.0) THEN
            IS = -IPS(IF)
            IF (IS.GT.0) FRERE(I) = IS
            IPS(IF) = -I
          ELSE
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (FRERE(INS).LT.0) THEN
       INS       = -FRERE(INS)
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (FRERE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = FRERE(INS)
#if ! defined(NOAMALGTOFATHER)
      IF (NV(INB).GE.NV(INS)) THEN
#else
      IF (NV(INB).GE.NV(INS)) THEN
#endif
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = FRERE(INF)
      IF (INF.GT.0) GO TO 1090
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
        FILS(INF) = -INB
        IPS(INF)  = -INB
        FRERE(INS)  = FRERE(INB)
        FRERE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = FRERE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        FRERE(INS) = FRERE(INB)
        FRERE(INB) = INS
        FRERE(INSW)= INB
      ENDIF
        INS      = INB
        GO TO 1070
 1151 DO 51 I=1,N
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      IPERM = 1
      DO 160 K=1,N
        IF (I.LE.0) THEN
         IF (NR.GT.N) EXIT
         I = NE(NR)
         NE(NR) = 0
         NR = NR + 1
         IL = N
         NA(N) = 0
        ENDIF
        DO 70 L=1,N
          IF (IPS(I).GE.0) EXIT
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
#if ! defined(NOAMALGTOFATHER)
        DADI = -IPE(I)  ! IPE was saved
        IF ( (DADI.NE.0) .AND.
     *      (
     *       (KEEP60.EQ.0).OR.
     *       ( (KEEP20.NE.DADI).AND.(KEEP38.NE.DADI) )
     *      ) 
     *     ) THEN
           ACCU = 
     &     ( DBLE(20000)*
     &       DBLE(NODE(I))*DBLE(NV(DADI)-NV(I)+NODE(I)) 
     &      )
     &      /
     &          ( DBLE(NV(DADI)+NODE(I))*
     &              DBLE(NV(DADI)+NODE(I))  )
           ACCU = ACCU + CUMUL(I) 
           AMALG_TO_father_OK = ( (NODE(I).LE.NEMIN).OR.
     &         (NODE(DADI).LE.NEMIN) )
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &       ( 
     &        ( DBLE(2*(NODE(I)))*
     &         DBLE((NV(DADI)-NV(I)+NODE(I)))
     &        ) .LT.
     &        (  DBLE(NV(DADI)+NODE(I))*
     &         DBLE(NV(DADI)+NODE(I))*DBLE(NEMIN)/DBLE(100)
     &        ) 
     &       ) )
           AMALG_TO_father_OK = ( AMALG_TO_father_OK .AND.
     &     ( ACCU .LE. DBLE(NEMIN)*DBLE(100) )
     &                           )
#if 0
           AMALG_TO_father_OK = (AMALG_TO_father_OK 
     *          .AND. NAMALG(I).LT.NAMALGMAX
     &          .AND. NAMALG(DADI).LT.NAMALGMAX)
#endif
#if 1
           IF (AMALG_TO_father_OK) THEN
              CALL CMUMPS_511(NV(I),NODE(I),NODE(I),
     *                                  KEEP50,1,FLOPS_SON)
              CALL CMUMPS_511(NV(DADI),NODE(DADI),
     *                             NODE(DADI),
     *                             KEEP50,1,FLOPS_FATHER)
              FLOPS_AVANT = FLOPS_FATHER+FLOPS_SON
     *                      + MAX(DBLE(200.0) * DBLE(NV(I)-NODE(I))
     *                            * DBLE(NV(I)-NODE(I)),
     *                            DBLE(10000.0))
              CALL CMUMPS_511(NV(DADI)+NODE(I),
     *                             NODE(DADI)+NODE(I),
     *                             NODE(DADI)+NODE(I),
     *                             KEEP50,1,FLOPS_APRES)
              IF(FLOPS_APRES .GT. FLOPS_AVANT) THEN
                 AMALG_TO_father_OK = .FALSE.
              ENDIF
           ENDIF
#endif     
#if 1
           IF (AMALG_TO_father_OK) THEN
              NAMALG(DADI) = NAMALG(DADI) +  NAMALG(I) + 1
              CUMUL(DADI)  =  CUMUL(DADI)+NINT(ACCU)
           ELSE
              AMALG_TO_father_OK = (AMALG_TO_father_OK .OR.
     &                       ( NV(I)-NODE(I).EQ.NV(DADI)) )
              IF (AMALG_TO_father_OK) THEN
                NAMALG(DADI) = NAMALG(DADI) +  NAMALG(I)
                CUMUL(DADI)  = CUMUL(DADI)+NINT(ACCU)
              ENDIF
           ENDIF
#endif
           IF (AMALG_TO_father_OK) THEN
             AMALG_COUNT = AMALG_COUNT+1
             IN = DADI
 75          IF (SUBORD(IN).EQ.0) GOTO 76
               IN = SUBORD(IN)
               GOTO 75
 76          CONTINUE
             SUBORD(IN) = I
             NV(I)      = 0
             IFSON = -FILS(DADI)
             IF (IFSON.EQ.I) THEN
              IF (FILS(I).LT.0) THEN
                FILS(DADI) =  FILS(I)
                GOTO 78
              ELSE
                IF (FRERE(I).GT.0) THEN
                  FILS(DADI) = -FRERE(I)  
                ELSE
                  FILS(DADI) = 0
                ENDIF
                GOTO 90
              ENDIF
             ENDIF
             IN = IFSON
  77         INS = IN
             IN = FRERE(IN)
             IF (IN.NE.I) GOTO 77
             IF (FILS(I) .LT.0) THEN
               FRERE(INS) = -FILS(I)
             ELSE
               FRERE(INS) = FRERE(I)  ! skip I which is a leave
               GOTO 90
             ENDIF
  78         CONTINUE
             IN = -FILS(I)
  79         INO = IN
             IN = FRERE(IN)
             IF (IN.GT.0) GOTO 79
             FRERE(INO) = FRERE(I)
  90         CONTINUE
             NODE(DADI) = NODE(DADI)+ NODE(I) 
             NV(DADI)   = NV(DADI) +  NODE(I) 
             NA(IL+1)   = NA(IL+1) + NA(IL)
             GOTO 120
           ENDIF
        ENDIF
#endif
        NE(IS) = NE(IS) + NODE(I) 
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        IPS(I) = IPERM
        IPERM = IPERM + 1
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          IPS(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
  778   IF (NA(IS).LE.0) GO TO 110
#if defined(NOAMALGTOFATHER)
        IF (   (KEEP60.NE.0).AND.
     &         (NE(IS).EQ.ND(IS)) ) GOTO 110
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) THEN
           GO TO 100
        ENDIF
        IF(NAMALG(IS-1) .GE. NAMALGMAX) THEN
           GOTO 110
        ENDIF
        IF ((NE(IS-1).GE.NEMIN).AND.
     &         (NE(IS).GE.NEMIN) ) GO TO 110
        IF (2*NE(IS-1)*(ND(IS)-ND(IS-1)+NE(IS-1)).GE.
     &    ((ND(IS)+NE(IS-1))*
     &    (ND(IS)+NE(IS-1))*NEMIN/100)) GO TO 110
        NAMALG(IS-1) = NAMALG(IS-1)+1
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
        IFSON = -FILS(I)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
        NV(INO) = 0
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        GO TO 888
  889   SUBORD(IN) = INO
      INOS = -FILS(INO)
      IF (IFSON.EQ.INO) THEN 
         FILS(I) = -INOS
         GO TO 107
      ENDIF
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
        IF (INOS.EQ.0) THEN
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT) = -I
        GO TO 120
#endif
  110   IS = IS + 1
  120   IB = FRERE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
      DO I=1, N
        IF (NV(I).EQ.0) THEN
          FRERE(I) = N+1
        ELSE
          NFSIZ(I) = ND(NODE(I))
          IF (SUBORD(I) .NE.0) THEN
           INOS = -FILS(I)  ! Real son
           INO = I
           DO WHILE (SUBORD(INO).NE.0) 
             IS = SUBORD(INO)
             FILS(INO) = IS
             INO = IS
           END DO
           FILS(INO) = -INOS
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_557
#endif
      SUBROUTINE CMUMPS_201(NE, ND, NSTEPS,
     * MAXFR, MAXELIM, K50, SIZEFAC, MAXFAC, SIZEFAC8,MAXNPIV )
      IMPLICIT NONE
      INTEGER NSTEPS,MAXNPIV
      INTEGER  MAXFR, MAXELIM, K50, SIZEFAC, MAXFAC
      INTEGER  NE(NSTEPS), ND(NSTEPS)
      INTEGER*8 SIZEFAC8, TEMP8
      INTEGER ITREE, NFR, NELIM
      MAXFR   = 0
      MAXFAC  = 0
      MAXELIM = 0
      MAXNPIV = 0
      SIZEFAC8 = 0
      DO ITREE=1,NSTEPS
        NELIM = NE(ITREE)
        NFR = ND(ITREE)
        IF (NFR.GT.MAXFR)         MAXFR   = NFR
        IF (NFR-NELIM.GT.MAXELIM) MAXELIM = NFR - NELIM
        IF (NELIM .GT. MAXNPIV) THEN
           IF(NFR .NE. NELIM) MAXNPIV = NELIM
        ENDIF
        IF (K50.EQ.0) THEN
          SIZEFAC8 = SIZEFAC8 + (2*NFR - NELIM)*NELIM 
          MAXFAC = MAX(MAXFAC, (2*NFR - NELIM)*NELIM )
        ELSE
         SIZEFAC8 = SIZEFAC8 + (NELIM*(NELIM+1))/2 + 
     &                       (NFR-NELIM)*NELIM
         MAXFAC = MAX(MAXFAC, NFR * NELIM)
        ENDIF
      END DO
      IF (SIZEFAC8.GT.HUGE(SIZEFAC)) THEN
         TEMP8 = SIZEFAC8 / 1000000
         IF (TEMP8.GT.HUGE(SIZEFAC)) THEN
           SIZEFAC=-1
         ELSE
           SIZEFAC = -TEMP8
         ENDIF
      ELSE
         SIZEFAC=SIZEFAC8
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_201
      SUBROUTINE CMUMPS_348( N, NZ, PERM, FILS, FRERE, ND,
     * NSTK, NA, IRN, ICN,
     * IW4, ICNTL, INFO, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER ICNTL(40),NBLEAF
      INTEGER FILS(N), FRERE(N), PERM(N), ND(N)
      INTEGER NSTK(N), NA(N) 
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER IW4(N)
      INTEGER NBROOT,ILEAF,I,IN,ISON
      INTEGER IOLD,K,JOLD,INEW,JNEW,ISHIFT
      NA   = 0
      NSTK = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
      IF (FRERE(I).EQ. N+1) CYCLE
      IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
      IN = I
 12   IN = FILS(IN)
      IF (IN.GT.0) GO TO 12
      IF (IN.EQ.0) THEN 
        NA(ILEAF) = I
        ILEAF     = ILEAF + 1
        CYCLE
      ENDIF
      ISON = -IN
 13   NSTK(I) = NSTK(I) + 1
      ISON = FRERE(ISON)
      IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IF (N.GT.1) THEN
       IF (NBLEAF.GT.N-2) THEN
        IF (NBLEAF.EQ.N-1) THEN
         NA(N-1) = -NA(N-1)-1
         NA(N)   = NBROOT
        ELSE
         NA(N) = -NA(N)-1
        ENDIF
       ELSE
        NA(N-1) = NBLEAF
        NA(N)   = NBROOT
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_348
      SUBROUTINE CMUMPS_202( N, NZ, PERM, FILS, FRERE, ND,
     * NSTK, NA, IRN, ICN,
     * PTRAIW, PTRARW, IW4, ICNTL, INFO, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER ICNTL(40),NBLEAF
      INTEGER FILS(N), FRERE(N), PERM(N), ND(N)
      INTEGER NSTK(N), NA(N) 
      INTEGER PTRAIW(N), PTRARW(N), IW4(N,2)
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER NBROOT,ILEAF,I,IN,ISON
      INTEGER IOLD,K,JOLD,INEW,JNEW,ISHIFT
      NA   = 0
      NSTK = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
      IF (FRERE(I).EQ. N+1) CYCLE
      IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
      IN = I
 12   IN = FILS(IN)
      IF (IN.GT.0) GO TO 12
      IF (IN.EQ.0) THEN 
        NA(ILEAF) = I
        ILEAF     = ILEAF + 1
        CYCLE
      ENDIF
      ISON = -IN
 13   NSTK(I) = NSTK(I) + 1
      ISON = FRERE(ISON)
      IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IW4(1:N,1) = NA(1:N)
      IF (N.GT.1) THEN
       IF (NBLEAF.GT.N-2) THEN
        IF (NBLEAF.EQ.N-1) THEN
         NA(N-1) = -NA(N-1)-1
         NA(N)   = NBROOT
        ELSE
         NA(N) = -NA(N)-1
        ENDIF
       ELSE
        NA(N-1) = NBLEAF
        NA(N)   = NBROOT
       ENDIF
      ENDIF
        IW4(1:N,2) = NSTK(1:N)
        CALL CMUMPS_181( N, NBROOT, NBLEAF, IW4(1,1), N,
     &        PERM, FILS, FRERE, IW4(1,2) )
      IF (KEEP(55).NE.0) RETURN
      DO 50 IOLD=1,N
        IW4(IOLD,1) = 0
        IW4(IOLD,2) = 0
   50 CONTINUE
      DO 70 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) GOTO 70
        IF (IOLD.NE.JOLD) THEN
         INEW = PERM(IOLD)
         JNEW = PERM(JOLD)
         IF ( KEEP( 50 ) .EQ. 0 ) THEN
           IF (INEW.LT.JNEW) THEN
             IW4(IOLD,2) = IW4(IOLD,2) + 1
           ELSE
             IW4(JOLD,1) = IW4(JOLD,1) + 1
           ENDIF
         ELSE
           IF ( INEW .LT. JNEW ) THEN
             IW4( IOLD, 1 ) = IW4( IOLD, 1 ) + 1
           ELSE 
             IW4( JOLD, 1 ) = IW4( JOLD, 1 ) + 1
           END IF
         ENDIF
        ENDIF
   70 CONTINUE
      PTRARW(1) = 1
      PTRAIW(1) = 1
      IF (N.GT.1) THEN
      DO 90 I=1,N-1
       ISHIFT = IW4(I,1) + IW4(I,2)
       PTRARW(I+1) = PTRARW(I) + ISHIFT + 1
       PTRAIW(I+1) = PTRAIW(I) + ISHIFT + 3
  90  CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_202
      SUBROUTINE CMUMPS_353 
     * ( N, PERM, FILS, FRERE, ND,
     * NSTK, NA, IW4, ICNTL, INFO, KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER N
      INTEGER INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER ICNTL(40),NBLEAF
      INTEGER FILS(N), FRERE(N), PERM(N), ND(N)
      INTEGER NSTK(N), NA(N) 
      INTEGER IW4(N,2)
      INTEGER NBROOT,ILEAF,I,IN,ISON
      INTEGER IOLD,K,JOLD,INEW,JNEW,ISHIFT
      NA   = 0
      NSTK = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
      IF (FRERE(I).EQ. N+1) CYCLE
      IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
      IN = I
 12   IN = FILS(IN)
      IF (IN.GT.0) GO TO 12
      IF (IN.EQ.0) THEN 
        NA(ILEAF) = I
        ILEAF     = ILEAF + 1
        CYCLE
      ENDIF
      ISON = -IN
 13   NSTK(I) = NSTK(I) + 1
      ISON = FRERE(ISON)
      IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IW4(1:N,1) = NA(1:N)
      IF (N.GT.1) THEN
       IF (NBLEAF.GT.N-2) THEN
        IF (NBLEAF.EQ.N-1) THEN
         NA(N-1) = -NA(N-1)-1
         NA(N)   = NBROOT
        ELSE
         NA(N) = -NA(N)-1
        ENDIF
       ELSE
        NA(N-1) = NBLEAF
        NA(N)   = NBROOT
       ENDIF
      ENDIF
        IW4(1:N,2) = NSTK(1:N)
        CALL CMUMPS_181( N, NBROOT, NBLEAF, IW4(1,1), N,
     &        PERM, FILS, FRERE, IW4(1,2) )
      RETURN
      END SUBROUTINE CMUMPS_353
      SUBROUTINE CMUMPS_203( N, NZ, MTRANS, PERM,
     &     id, ICNTL, INFO)
      USE CMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (CMUMPS_STRUC) :: id
      INTEGER N, NZ, LIWG
      INTEGER PERM(N)
      INTEGER MTRANS
      INTEGER ICNTL(40), INFO(40)
      INTEGER  allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IW
      REAL, ALLOCATABLE, DIMENSION(:) :: S2
      TARGET :: S2
      INTEGER LS2,LSC
      INTEGER ICNTL64(10), INFO64(10)
      INTEGER ICNTL_SYM_MWM(10),INFO_SYM_MWM(10)
      REAL CNTL64(10)
      INTEGER LDW, LDWMIN
      INTEGER MPRINT,LP, MP, IPIW, LIW, LIWMIN
      INTEGER JPERM
      INTEGER NUMNZ, I, J, JPOS, K, NZREAL
      INTEGER PLENR, IP, IRNW,RSPOS,CSPOS
      LOGICAL PROK, IDENT, DUPPLI
      INTEGER NZTOT, K50, KER_SIZE, NZER_DIAG, MTRANSLOC,RZ_DIAG
      LOGICAL SCALINGLOC
      INTEGER,POINTER,DIMENSION(:) :: ZERODIAG
      INTEGER,POINTER,DIMENSION(:) :: STR_KER
      INTEGER,POINTER,DIMENSION(:) :: MARKED
      INTEGER,POINTER,DIMENSION(:) :: FLAG
      INTEGER,POINTER,DIMENSION(:) :: PIV_OUT
      REAL THEMIN, THEMAX, COLNORM,MAXDBL
      REAL ZERO,TWO,ONE
      PARAMETER(ZERO = 0.0E0,TWO = 2.0E0,ONE = 1.0E0)
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK = (MPRINT.GT.0)
      IF (PROK) WRITE(MPRINT,101)
 101  FORMAT(/'****** Preprocessing of original matrix '/)
      K50 = id%KEEP(50)
      SCALINGLOC = .FALSE.
      IF(id%KEEP(52) .EQ. -2) THEN
         IF(.not.ASSOCIATED(id%A)) THEN
            INFO(1) = -22
            INFO(2) = 4
            GOTO 500
         ELSE
            SCALINGLOC = .TRUE.
         ENDIF
      ELSE IF(id%KEEP(52) .EQ. 7) THEN
         SCALINGLOC = .TRUE.
         IF(K50 .NE. 2) THEN
            IF( MTRANS .NE. 5 .AND. MTRANS .NE. 6 
     *           .AND. MTRANS .NE. 7) THEN
               SCALINGLOC = .FALSE.
               IF (PROK) 
     *              WRITE(MPRINT,*) 'Analysis: auto scaling set OFF'
            ENDIF
         ENDIF
         IF(.not.ASSOCIATED(id%A)) THEN
            SCALINGLOC = .FALSE.
            IF (PROK) 
     *           WRITE(MPRINT,*) 'Analysis: auto scaling set OFF'
         ENDIF
      ENDIF
      IF(SCALINGLOC) THEN
         IF (PROK) WRITE(MPRINT,*) 
     *        'Scaling will be computed during analysis'
      ENDIF
      MTRANSLOC = MTRANS
      IF (MTRANS.LT.0 .OR. MTRANS.GT.7) GO TO 500
      IF (K50 .EQ. 0) THEN
         IF(.NOT. SCALINGLOC .AND. MTRANS .EQ. 7) THEN 
            GO TO 500
         ENDIF
         IF(SCALINGLOC) THEN
            MTRANSLOC = 5
         ENDIF
      ELSE
         IF (MTRANS .EQ. 7) MTRANSLOC = 5
      ENDIF
      IF(SCALINGLOC .AND. MTRANSLOC .NE. 5 .AND.
     *     MTRANSLOC .NE. 6 ) THEN
         IF (PROK) WRITE(MPRINT,*)
     *        'WARNING scaling required: set MTRANS option to 5'
         MTRANSLOC = 5
      ENDIF
      IF (N.EQ.1) THEN
        MTRANS=0
        GO TO 500
      ENDIF
      IF(K50 .EQ. 2) THEN
         NZTOT = 2*NZ+N
      ELSE
         NZTOT = NZ
      ENDIF
      ZERODIAG => id%IS1(N+1:2*N)
      STR_KER => id%IS1(2*N+1:3*N)
      CALL CMUMPS_448(ICNTL64,CNTL64)
      ICNTL64(1) = ICNTL(1)
      ICNTL64(2) = ICNTL(2)
      ICNTL64(3) = ICNTL(2)
      ICNTL64(4) = -1
      IF (ICNTL(4).EQ.3) ICNTL64(4) = 0
      IF (ICNTL(4).EQ.4) ICNTL64(4) = 1
      ICNTL64(5) = -1
      IF (PROK) THEN
         WRITE(MPRINT,'(A,I3)')
     *     'Compute maximum matching (Maximum Transversal):',
     *        MTRANSLOC
         IF (MTRANSLOC.EQ.1)
     *   WRITE(MPRINT,'(A,I3)')' ... JOB =',MTRANSLOC
         IF (MTRANSLOC.EQ.2)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': BOTTLENECK THESIS'
         IF (MTRANSLOC.EQ.3)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': BOTTLENECK SIMAX'
         IF (MTRANSLOC.EQ.4)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,': MAXIMIZE SUM DIAGIONAL'
         IF (MTRANSLOC.EQ.5 .OR. MTRANSLOC.EQ.6)
     *   WRITE(MPRINT,'(A,I3,A)')
     *     ' ... JOB =',MTRANSLOC,
     *     ': MAXIMIZE PRODUCT DIAGONAL AND SCALE'
      ENDIF
      id%INFOG(23) = MTRANSLOC
      CNTL64(2) = HUGE(CNTL64(2))
      IRNW = 1
      IP = IRNW + NZTOT
      PLENR = IP + N + 1
      IPIW = PLENR
      IF (MTRANSLOC.EQ.1) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.2) LIWMIN = 4*N
      IF (MTRANSLOC.EQ.3) LIWMIN = 10*N + NZTOT
      IF (MTRANSLOC.EQ.4) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.5) LIWMIN = 5*N
      IF (MTRANSLOC.EQ.6) LIWMIN = 5*N + NZTOT
      LIW = LIWMIN
      LIWG  = LIW + (NZTOT + N + 1)
      ALLOCATE(IW(LIWG), stat=allocok)
      IF (allocok .GT. 0 ) GOTO 410
      IF (MTRANSLOC.EQ.1) THEN
       LDWMIN = N+3
      ENDIF
      IF (MTRANSLOC.EQ.2) LDWMIN = MAX(N+NZTOT,N+3)
      IF (MTRANSLOC.EQ.3) LDWMIN = MAX(NZTOT+1,N+3)
      IF (MTRANSLOC.EQ.4) LDWMIN = 2*N + MAX(NZTOT,N+3)
      IF (MTRANSLOC.EQ.5) LDWMIN = 3*N + NZTOT
      IF (MTRANSLOC.EQ.6) LDWMIN = 4*N + NZTOT
      LDW   = LDWMIN
      IF (allocok .GT. 0 ) GOTO 420
      ALLOCATE(S2(LDW), stat=allocok)
      IF(MTRANSLOC .NE. 1) LDW = LDW-NZTOT
      RSPOS = NZTOT
      CSPOS = RSPOS+N
      IF (allocok .GT. 0 ) GOTO 430
      NZREAL = 0
      DO 5 J=1,N
        IW(PLENR+J-1) = 0
  5   CONTINUE
      IF(K50 .EQ. 0) THEN
         DO 10 K=1,NZ
            I = id%IRN(K)
            J = id%JCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1) ) THEN
               IW(PLENR+J-1) = IW(PLENR+J-1) + 1
               NZREAL = NZREAL + 1
            ENDIF
 10      CONTINUE
      ELSE
         ZERODIAG = 0
         NZER_DIAG = N
         RZ_DIAG = 0
         DO K=1,NZ
            I = id%IRN(K)
            J = id%JCN(K)
            IF ( (J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1) ) THEN
               IW(PLENR+J-1) = IW(PLENR+J-1) + 1
               NZREAL = NZREAL + 1
               IF(I .NE. J) THEN
                  IW(PLENR+I-1) = IW(PLENR+I-1) + 1
                  NZREAL = NZREAL + 1
               ELSE
                  IF(ZERODIAG(I) .EQ. 0) THEN
                     ZERODIAG(I) = K
                     IF(ASSOCIATED(id%A)) THEN
                        IF(ABS(id%A(K)) .EQ. REAL(0.0E0)) THEN
                           RZ_DIAG = RZ_DIAG + 1
                        ENDIF
                     ENDIF
                     NZER_DIAG = NZER_DIAG - 1                     
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
#if(SP_DEBUG >= 1)
         WRITE(*,*) 'NZER_DIAG,RZ_DIAG,TOTZER',NZER_DIAG,RZ_DIAG,
     *        NZER_DIAG+RZ_DIAG
#endif
         IF(MTRANSLOC .GE. 4) THEN
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  IW(PLENR+I-1) = IW(PLENR+I-1) + 1
                  NZREAL = NZREAL + 1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IW(IP)   = 1
      DO 20 J=1,N
        IW(IP+J)   = IW(IP+J-1)+IW(PLENR+J-1)
  20  CONTINUE
      DO 25 J=1, N
        IW(PLENR+J-1 ) = IW(IP+J-1 )
#if (SP_DEBUG >= 2)
        IF(IW(PLENR+J-1) .LE. 0) THEN
           WRITE(*,*) 'J,IW(PLENR+J-1)',J,IW(PLENR+J-1)
           CALL CMUMPS_ABORT()
        ENDIF
#endif
  25  CONTINUE
      IF(K50 .EQ. 0) THEN
         IF (MTRANSLOC.EQ.1) THEN
            DO 30 K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
              ENDIF
 30         CONTINUE
         ELSE
            IF ( .not.ASSOCIATED(id%A)) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            DO 35 K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = ABS(id%A(K))
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
               ENDIF
 35         CONTINUE
         ENDIF
      ELSE
         IF (MTRANSLOC.EQ.1) THEN
            DO K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
                  IW(IRNW+JPOS-1) = I
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
                  IF(I.NE.J) THEN
                     JPOS            = IW(PLENR+I-1)
                     IW(IRNW+JPOS-1) = J
                     IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            IF ( .not.ASSOCIATED(id%A)) THEN
               INFO(1) = -22
               INFO(2) = 4
               GOTO 500
            ENDIF
            K = 1
            THEMIN = ZERO
            DO
               IF(THEMIN .NE. ZERO) EXIT
               THEMIN = ABS(id%A(K))
               K = K+1
            ENDDO
            THEMAX = THEMIN
            DO K=1,NZ
               I = id%IRN(K)
               J = id%JCN(K)
               IF ( (J.LE.N).AND.(J.GE.1) .AND.
     *              (I.LE.N).AND.(I.GE.1)) THEN
                  JPOS            = IW(PLENR+J-1)
#if (SP_DEBUG >= 2)
                  IF(JPOS .LE. 0) THEN
                     WRITE(*,*) 'I,J,JPOS,K',I,J,JPOS,K
                     CALL CMUMPS_ABORT()
                  ENDIF
#endif
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = ABS(id%A(K))
                  IW(PLENR+J-1)   = IW(PLENR+J-1) + 1
#if (SP_DEBUG >= 2)
                  IF(IW(PLENR+J-1) .LE. 0) THEN
                     WRITE(*,*) 'C0 J,IW(PLENR+J-1)',J,IW(PLENR+J-1)
                     CALL CMUMPS_ABORT()
                  ENDIF
#endif
                  IF(ABS(id%A(K)) .GT. THEMAX) THEN
                     THEMAX = ABS(id%A(K))
                  ELSE IF(ABS(id%A(K)) .LT. THEMIN 
     *                    .AND. ABS(id%A(K)).GT. ZERO) THEN
                     THEMIN = ABS(id%A(K))
                  ENDIF
                  IF(I.NE.J) THEN
                     JPOS            = IW(PLENR+I-1)
#if (SP_DEBUG >= 2)
                     IF(JPOS .LE. 0) THEN
                        WRITE(*,*) 'C1 I,J,JPOS,K,PLENR',
     *                       I,J,JPOS,K,PLENR
                        CALL CMUMPS_ABORT()
                     ENDIF
#endif
                     IW(IRNW+JPOS-1) = J
                     S2(JPOS)         = ABS(id%A(K))
                     IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
#if (SP_DEBUG >= 2)
                     IF(IW(PLENR+J-1) .LE. 0) THEN
                        WRITE(*,*) 'C2 J,IW(PLENR+J-1)',J,IW(PLENR+J-1)
                        CALL CMUMPS_ABORT()
                     ENDIF
#endif
                  ENDIF
               ENDIF
#if (SP_DEBUG >= 3)
               DO JPOS = 1,N
                  IF(IW(PLENR+JPOS-1) .LT. IW(IP+JPOS-1)
     *                 .OR. IW(PLENR+JPOS-1) .GT. IW(IP+JPOS)) THEN
                     WRITE(*,*) 
     *            'JPOS,IW(PLENR+JPOS-1),IW(IP+JPOS-1),IW(IP+JPOS)',
     *                 JPOS,IW(PLENR+JPOS-1),IW(IP+JPOS-1),IW(IP+JPOS)
                     CALL CMUMPS_ABORT()
                  ENDIF
               ENDDO
               IF(IW(PLENR+3289-1) .LE. 0) THEN
                  WRITE(*,*) 'bbb I,J,JPOS,K',I,J,JPOS,K
                  WRITE(*,*) 'IW(PLENR+J-1)',IW(PLENR+3289-1)
                  CALL CMUMPS_ABORT()
               ENDIF
#endif
            ENDDO
            DO I =1, N
               IF(ZERODIAG(I) .EQ. 0) THEN
                  JPOS            = IW(PLENR+I-1)
                  IW(IRNW+JPOS-1) = I
                  S2(JPOS)         = ZERO
                  IW(PLENR+I-1)   = IW(PLENR+I-1) + 1
               ENDIF
            ENDDO
            CNTL64(2) = (LOG(THEMAX/THEMIN))*(CMPLX(N))
     *           - LOG(THEMIN) + ONE
         ENDIF
      ENDIF
      DUPPLI = .FALSE.
      I = NZREAL
      FLAG => id%IS1(3*N+1:4*N)
      IF(MTRANSLOC.NE.1) THEN
         CALL CMUMPS_563(N,NZREAL,IW(IP),IW(IRNW),S2,
     *        PERM,FLAG)
      ELSE
         CALL CMUMPS_562(N,NZREAL,IW(IP),IW(IRNW),
     *        PERM,FLAG)
      ENDIF
      IF(NZREAL .NE. I) DUPPLI = .TRUE.
      LS2 = NZTOT
      IF ( MTRANSLOC .EQ. 1 ) THEN
         LS2 = 1
         LDW = 1
      ENDIF
      CALL CMUMPS_559(MTRANSLOC ,N, N, NZREAL, 
     *     IW(IP), IW(IRNW), S2(1), LS2,
     *     NUMNZ, PERM, LIW, IW(IPIW), LDW, S2(LS2+1),
     *     ICNTL64, CNTL64, INFO64)
      IF (INFO64(1).LT.0) THEN
         IF (LP.GT.0 .AND. ICNTL(4).GE.1)
     *        WRITE(LP,'(A,I5)')
     *   ' INTERNAL ERROR in MAXTRANS INFO(1)=',INFO64(1)
         INFO(1) = -9964
         INFO(2) = INFO64(1)
         GO TO 500
      ENDIF
      IF (INFO64(1).GT.0) THEN
         IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *        WRITE(MP,'(A,I5)')
     *        ' WARNING in MAXTRANS INFO(1)=',INFO64(1)
      ENDIF
      KER_SIZE = 0
      IF(K50 .EQ. 2) THEN
         DO I=1,N
            IF(ZERODIAG(I) .EQ. 0) THEN
               IF(PERM(I) .EQ. I) THEN
                  KER_SIZE = KER_SIZE + 1
                  PERM(I) = -I
                  STR_KER(KER_SIZE) = I
               ENDIF
            ENDIF
         ENDDO
      ENDIF
      IF (NUMNZ.LT.N) GO TO 400
      IF(KER_SIZE .GT. 0) GOTO 400
      IF(K50 .EQ. 0) THEN
         IDENT = .TRUE.
         IF(MTRANS .EQ. 0) GOTO 102
         DO 80 J=1,N
            JPERM = PERM(J)
            IW(PLENR+JPERM-1) = J
            IF (JPERM.NE.J) IDENT = .FALSE.
 80      CONTINUE
         IF(IDENT) THEN 
            MTRANS = 0
         ELSE
            IF(MTRANS .EQ. 7) THEN
               MTRANS = -9876543
               GOTO 102
            ENDIF
         ENDIF
         IF (.NOT.IDENT) THEN
            IF (PROK) WRITE(MPRINT,'(A)')
     *           ' ... Apply column permutation'
            DO 100 K=1,NZ
               J = id%JCN(K)
               IF ((J.LE.0).OR.(J.GT.N)) GO TO 100
               id%JCN(K) = IW(PLENR+J-1)
 100        CONTINUE
            IF (MP.GT.0 .AND. ICNTL(4).GE.2)
     *           WRITE(MP,'(/A)')
     *           ' WARNING input matrix data modified'
         ENDIF
 102     CONTINUE
         IF (SCALINGLOC) THEN
            IF ( associated(id%COLSCA))
     *           DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     *           DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in CMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in CMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            id%KEEP(52) = -2
            id%KEEP(74) = 1
            MAXDBL = LOG(HUGE(MAXDBL))
            DO J=1,N
               IF(S2(RSPOS+J) .GT. MAXDBL) THEN
                  S2(RSPOS+J) = ZERO
               ENDIF
               IF(S2(CSPOS+J) .GT. MAXDBL) THEN
                  S2(CSPOS+J)= ZERO
               ENDIF
            ENDDO
            DO 105 J=1,N
               id%ROWSCA(J) = CMPLX(EXP(S2(RSPOS+J)))
               IF(id%ROWSCA(J) .EQ. ZERO) THEN
                  id%ROWSCA(J) = ONE
               ENDIF
               id%COLSCA(IW(PLENR+J-1))= CMPLX(EXP(S2(CSPOS+J)))
               IF(id%COLSCA(IW(PLENR+J-1)) .EQ. ZERO) THEN
                  id%COLSCA(IW(PLENR+J-1)) = ONE
               ENDIF
 105        CONTINUE
         ENDIF
      ELSE
         IDENT = .FALSE.         
#if (SP_DEBUG >= 3)
         WRITE(*,*) PERM
         WRITE(*,*) S2(RSPOS+1:RSPOS+N)
         WRITE(*,*) S2(CSPOS+1:CSPOS+N)
#endif
         IF(SCALINGLOC) THEN
            IF ( associated(id%COLSCA)) DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA)) DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in CMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of COLSCA'
                  GOTO 500
               ENDIF
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=allocok)
            IF (allocok .GT.0) THEN
               id%INFO(1)=-5
               id%INFO(2)=N
               IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
                  WRITE (LP,'(/A)') '** Error in CMUMPS_203'
                  WRITE (LP,'(A)') 
     *                 '** Failure during allocation of ROWSCA'
                  GOTO 500
               ENDIF
            ENDIF
            id%KEEP(52) = -2
            id%KEEP(74) = 1
#if 0
            J = 5000
            I = 1
            DO 
               IF(PERM(J) .EQ. 5000) EXIT
               J = PERM(J)
               I = I+1
            ENDDO
            WRITE(*,*) 'cycle 5000 lenght = ',I
#endif
            MAXDBL = LOG(HUGE(MAXDBL))
            DO J=1,N
               IF(S2(RSPOS+J)+S2(CSPOS+J) .GT. MAXDBL) THEN
                  S2(RSPOS+J) = ZERO
                  S2(CSPOS+J)= ZERO
               ENDIF
            ENDDO
            DO J=1,N
               IF(PERM(J) .GT. 0) THEN
                  id%ROWSCA(J) = 
     *                 CMPLX(EXP((S2(RSPOS+J)+S2(CSPOS+J))/TWO))
                  IF(id%ROWSCA(J) .EQ. ZERO) THEN
                     id%ROWSCA(J) = ONE
                  ENDIF
                  id%COLSCA(J)= id%ROWSCA(J)
#if (SP_DEBUG >= 3)
                  WRITE(*,*) 'J,ROWSCA(J)',J,id%ROWSCA(J)
#endif
               ENDIF
            ENDDO
            DO JPOS=1,KER_SIZE
               I = STR_KER(JPOS)
               COLNORM = ZERO
               DO J = IW(IP+I-1),IW(IP+J) - 1
                  COLNORM = MAX(COLNORM,S2(J))
               ENDDO
               COLNORM = EXP(COLNORM) 
               id%ROWSCA(I) = ONE / COLNORM
               id%COLSCA(I) = id%ROWSCA(I)
#if (SP_DEBUG >= 3)
               WRITE(*,*) 'I,ROWSCA(I)',I,id%ROWSCA(I)
#endif
            ENDDO
#if (SP_DEBUG >= 2)
            WRITE(*,*) 'check scaling factors'
            S2 = ZERO
            DO K=1,NZ
               I=id%IRN(K)
               J=id%JCN(K)
               S2(I) = MAX(ABS(S2(I)),
     *              ABS(id%ROWSCA(I)*id%A(K)*id%COLSCA(J)))
               S2(J) = MAX(ABS(S2(J)),
     *              ABS(id%ROWSCA(I)*id%A(K)*id%COLSCA(J)))
            ENDDO
            COLNORM = SQRT(EPSILON(COLNORM))
            DO K = 1,N
               IF(PERM(K) .GT. 0 .AND. .NOT. DUPPLI) THEN
                  IF( ABS(S2(K)) .LT. (ONE-COLNORM)
     *                 .OR. ABS(S2(K)) .GT. (ONE+COLNORM) ) THEN
                     WRITE(*,*) 'pb 1 K,S2(K),PERM(K),PERM(PERM(K))',
     *                   K,S2(K),PERM(K),PERM(PERM(K))
                     WRITE(*,*) 'id%ROWSCA(K),id%ROWSCA(PERM(K))',
     *                           id%ROWSCA(K),id%ROWSCA(PERM(K)),
     *                           'id%COLSCA(K),id%COLSCA(PERM(K))',
     *                            id%COLSCA(K),id%COLSCA(PERM(K)) 
                     CALL CMUMPS_ABORT()
                  ENDIF
               ELSE
                  IF( ABS(S2(K)) .GT. (ONE+COLNORM) ) THEN
                     WRITE(*,*) 'pb 2 K,S2(K),PERM(K)',K,S2(K),PERM(K)
                     CALL CMUMPS_ABORT()
                  ENDIF
               ENDIF
            ENDDO
            IF(DUPPLI) THEN
               S2 = ZERO
               DO K=1,NZ
                  I=id%IRN(K)
                  J=id%JCN(K)
                  IF(I .EQ. PERM(J)) THEN
                     S2(I) = S2(I) + 
     *                    ABS(id%ROWSCA(I)*id%A(K)*id%COLSCA(J))
                     IF(I.NE. J) THEN
                        S2(J) = S2(J) + 
     *                       ABS(id%ROWSCA(I)*id%A(K)*id%COLSCA(J))
                     ENDIF
                  ENDIF
               ENDDO
               COLNORM = SQRT(EPSILON(COLNORM))
               DO K = 1,N
                  IF(PERM(K) .GT. 0) THEN
                     IF( ABS(S2(K)) .LT. (ONE-COLNORM)
     *                    .OR. ABS(S2(K)) .GT. (ONE+COLNORM) ) THEN
                        WRITE(*,*) 'pb 3 K,S2(K),PERM(K),PERM(PERM(K))',
     *                       K,S2(K),PERM(K),PERM(PERM(K))
                        WRITE(*,*) 'id%ROWSCA(K),id%ROWSCA(PERM(K))',
     *                       id%ROWSCA(K),id%ROWSCA(PERM(K)),
     *                       'id%COLSCA(K),id%COLSCA(PERM(K))',
     *                       id%COLSCA(K),id%COLSCA(PERM(K)) 
                        CALL CMUMPS_ABORT()
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            WRITE(*,*) 'scaling factors OK'
#endif
         ENDIF
         IF(MTRANS .EQ. 7 .OR. id%KEEP(95) .EQ. 0) THEN
            IF( (NZER_DIAG+RZ_DIAG) .LT. (N/10) 
     *           .AND. id%KEEP(95) .EQ. 0) THEN
               MTRANS = 0
               id%KEEP(95) = 1
               GOTO 390
            ELSE
               IF(id%KEEP(95) .EQ. 0) THEN
                 IF(SCALINGLOC) THEN
                  id%KEEP(95) = 3
                 ELSE
                  id%KEEP(95) = 2   
                 ENDIF
               ENDIF
               IF(MTRANS .EQ. 7) MTRANS = 5
            ENDIF
         ENDIF
         IF(MTRANS .EQ. 0) GOTO 390
         ICNTL_SYM_MWM = 0
         INFO_SYM_MWM = 0
         IF(MTRANS .EQ. 5 .OR. MTRANS .EQ. 6 .OR.
     *        MTRANS .EQ. 7) THEN
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ELSE IF(MTRANS .EQ. 4) THEN
            ICNTL_SYM_MWM(1) = 2
            ICNTL_SYM_MWM(2) = 1
         ELSE
            ICNTL_SYM_MWM(1) = 0
            ICNTL_SYM_MWM(2) = 1
         ENDIF
         MARKED => id%IS1(2*N+1:3*N)
         FLAG => id%IS1(3*N+1:4*N)
         PIV_OUT => id%IS1(4*N+1:5*N)
         IF(MTRANSLOC .LT. 4) THEN
            LSC = 1
         ELSE
            LSC = 2*N
         ENDIF
         CALL CMUMPS_551(
     *        N, NZREAL, IW(IP), IW(IRNW), S2(1),LSC, PERM, ZERODIAG,
     *        ICNTL_SYM_MWM, S2(LSC+1),MARKED,FLAG,
     *        PIV_OUT, INFO_SYM_MWM)
         IF(INFO_SYM_MWM(1) .NE. 0) THEN
            WRITE(*,*) '** Error in CMUMPS_203'
            RETURN
         ENDIF
         IF(INFO_SYM_MWM(3) .EQ. N) THEN
            IDENT = .TRUE.
         ELSEIF( (N-INFO_SYM_MWM(4)-INFO_SYM_MWM(3)) .GT. N/10
     *           ) THEN
            IDENT = .TRUE.
            id%KEEP(95) = 1
         ELSE
            DO I=1,N
               PERM(I) = PIV_OUT(I)
            ENDDO
         ENDIF
         id%KEEP(93) = INFO_SYM_MWM(4)
         id%KEEP(94) = INFO_SYM_MWM(3)
         IF (IDENT) MTRANS=0
      ENDIF
 390  IF(MTRANS .EQ. 0) THEN
         id%KEEP(95) = 1 
         IF (PROK) THEN
           WRITE (MPRINT,'(A)')
     *  ' ... Column permutation not used'
         ENDIF
      ENDIF
      GO TO 500
 400  IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     *   WRITE (LP,'(/A)') '** Error: Matrix is structurally singular'
      INFO(1) = -6
      INFO(2) = NUMNZ
      GOTO 500
 410  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in CMUMPS_203'
       WRITE (LP,'(A,I9)')
     & '** Failure during allocation of INTEGER array of size ',
     & LIWG
      ENDIF
      INFO(1) = -5
      INFO(2) = LIWG
      GOTO 500
 420  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in CMUMPS_203'
       WRITE (LP,'(A)') '** Failure during allocation of S1'
      ENDIF
      INFO(1) = -5
      INFO(2) = NZTOT
      GOTO 500
 430  IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
       WRITE (LP,'(/A)') '** Error in CMUMPS_203'
       WRITE (LP,'(A)') '** Failure during allocation of S2'
      ENDIF
      INFO(1) = -5
      INFO(2) = LDW
 500  CONTINUE
      IF (ALLOCATED(IW)) DEALLOCATE(IW)
      IF (ALLOCATED(S2)) DEALLOCATE(S2)
      RETURN
      END SUBROUTINE CMUMPS_203
      SUBROUTINE CMUMPS_21( N, MYID,PROCNODE, ND,
     *                         SSARBR, NBSA, FRERE, FILS,
     *                         SLAVEF,
     *                         KEEP,KEEP8, ICNTL, INFO )
      IMPLICIT NONE
      INTEGER KEEP(500), ICNTL(40), INFO(40)
      INTEGER*8 KEEP8(150)
      INTEGER N, NBSA, SLAVEF, MYID
      INTEGER PROCNODE( N ), ND( N )
      INTEGER SSARBR( NBSA ), FRERE( N ), FILS( N )
      LOGICAL  CMUMPS_167, CMUMPS_283
      INTEGER  CMUMPS_330, CMUMPS_275, CMUMPS_497
      EXTERNAL CMUMPS_330, CMUMPS_275,
     *         CMUMPS_167, CMUMPS_283, CMUMPS_497
      INTEGER I, INODE, INODEROOT, IN, NPIV, NFRONT, NFR, NB_NIV2
      INTEGER IROOTTREE, SIZEROOT, PROCVAL, NSLAVES, KMAX
      INTEGER NNODE, ICT, TPN, MP, allocok
      LOGICAL ONE_LEVEL2
      INTEGER, ALLOCATABLE :: NUMNODE( : ), COST( : ),
     *                        MEMUSED( : ), NODEMAP(:)
      MP = ICNTL(3)
      SIZEROOT  = -1
      IROOTTREE = -1
      NB_NIV2  = 0
      IF ( SLAVEF .eq. 1 ) THEN
       DO INODE = 1, N
          PROCNODE(INODE) = PROCNODE( INODE ) - SLAVEF
       END DO
      ELSE
       DO I = 1, NBSA
         INODEROOT = SSARBR( I )
         PROCNODE( INODEROOT ) = PROCNODE( INODEROOT )
     *                         - SLAVEF
         INODE = INODEROOT
 30      CONTINUE
         IN = INODE
         DO WHILE ( IN .NE. 0 )
           INODE = IN
           DO WHILE ( IN .GT. 0 )
             IN = FILS( IN )
           END DO
           IF ( IN .LT. 0 ) IN = - IN
         END DO
  10     CONTINUE
         IF ( INODE .EQ. INODEROOT ) GOTO 100
         PROCNODE( INODE ) = PROCNODE( INODE ) - 2 * SLAVEF
         IN = FRERE( INODE )
         INODE = ABS( IN )
         IF ( IN .LT. 0 ) THEN
             GO TO 10
         ELSE
             GO TO 30
         END IF
 100     CONTINUE
        END DO
      END IF
      ONE_LEVEL2 = .FALSE.
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 .AND.
     *     .NOT. CMUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
          IN = INODE
          NPIV = 0
          DO WHILE ( IN.GT.0 )
            IN = FILS( IN )
            NPIV = NPIV + 1
          END DO
          NFRONT = ND( INODE )
          IF ( (CMUMPS_283( INODE,PROCNODE,SLAVEF).AND.
     &     (IN.NE.0.OR.NFRONT.NE.NPIV )) .OR. SLAVEF.EQ.1 ) GOTO 991
          IF ( NFRONT-NPIV > KEEP(9) .AND.
     *         ( NPIV > KEEP(4) .OR. KEEP(45) .EQ. 0 ) 
     *         .AND. ICNTL(40) .EQ. 0 ) THEN
            PROCNODE( INODE ) = PROCNODE( INODE ) + SLAVEF
            ONE_LEVEL2 = .TRUE.
            NB_NIV2 = NB_NIV2 + 1
          END IF
          IF ( NFRONT .EQ. NPIV ) THEN
            IF ( NFRONT .GT. SIZEROOT ) THEN
              IROOTTREE = INODE
              SIZEROOT  = NFRONT
            END IF
          END IF
        END IF
 991    CONTINUE
        IF ( ( SLAVEF .eq. 1 .OR.
     *         CMUMPS_283( INODE, PROCNODE, SLAVEF ) )
     *         .AND. KEEP(53) .NE. 0 .AND.
     *         FRERE( INODE ) .LT. N + 1 ) THEN
          IN = INODE
          NPIV = 0
          DO WHILE ( IN.GT.0 )
            IN = FILS( IN )
            NPIV = NPIV + 1
          END DO
          NFRONT = ND( INODE )
          IF ( NFRONT .EQ. NPIV ) THEN
            IF ( NFRONT .GT. SIZEROOT ) THEN
              IROOTTREE = INODE
              SIZEROOT  = NFRONT
            END IF
          END IF
        END IF
      END DO
      KEEP ( 56 ) = NB_NIV2
      IF ( SIZEROOT .GT. KEEP( 37 )
     *       .and. (ICNTL(13) .LE. 0 )
     *       .and. KEEP(60) .EQ. 0
#ifndef null_space_old
     *       .and. KEEP(53) .eq. 0
#endif
     *  ) THEN
         IF ( MP .GT. 0 )
     *   WRITE( MP, * ) 'A root of estimated size ',SIZEROOT,
     *   ' has been selected for Scalapack'
        KEEP(38) = IROOTTREE
        PROCNODE( IROOTTREE ) = 1 + 2 * SLAVEF
      ELSE IF (KEEP(60) .EQ. 2 .OR. KEEP(60).EQ.3) THEN
        PROCNODE( KEEP(38) ) = 1 + 2 * SLAVEF
        IF ( MP .GT. 0 )
     *  WRITE( MP, * ) ' Largest root node of size ',SIZEROOT
      ELSE
        KEEP(38) = 0
        IF ( MP .GT. 0 )
     *  WRITE( MP, * ) ' Largest root node of size ',SIZEROOT
      END IF
#if defined(try_nullspace)
      IF (  KEEP( 38 ) .eq. 0  .and. KEEP(53).NE.0 ) THEN
       KEEP(20) = IROOTTREE
      ELSE IF ( KEEP(60) .eq. 0 ) THEN
       KEEP(20) = 0
      ENDIF
#endif
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          IN = FILS( INODE )
          DO WHILE ( IN > 0 )
            PROCNODE( IN ) = PROCNODE( INODE )
            IN = FILS( IN )
          END DO
        END IF
      END DO
      IF ( KEEP( 45 ) .NE. 0 ) RETURN
      NNODE = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          TPN = CMUMPS_330( INODE, PROCNODE, SLAVEF )
          IF ( TPN .NE. 3 .AND.
     *         .not. CMUMPS_283( INODE, PROCNODE, SLAVEF ) .AND.
     *         .not. CMUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
               NNODE = NNODE + 1
          END IF
        END IF
      END DO
      IF ( NNODE .GT. 0 ) THEN
        ALLOCATE( NUMNODE( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( NODEMAP( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( COST   ( NNODE  ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = NNODE
          RETURN
        endif
        ALLOCATE( MEMUSED( SLAVEF ), stat = allocok )
        if (allocok .gt. 0) THEN
          INFO(1) = -7
          INFO(2) = SLAVEF
          RETURN
        endif
        MEMUSED( 1: SLAVEF ) = 0
        ICT = 0
        DO INODE = 1, N
        IF ( FRERE( INODE ) .LT. N + 1 ) THEN
          TPN = CMUMPS_330( INODE, PROCNODE, SLAVEF )
          IF ( TPN .NE. 3 .AND.
     *       .not. CMUMPS_283( INODE, PROCNODE, SLAVEF ) .AND.
     *       .not. CMUMPS_167( INODE, PROCNODE, SLAVEF ) ) THEN
             ICT            = ICT + 1
             NUMNODE( ICT ) = INODE
             NFR  = ND( INODE )
             NPIV = 0
             IN   = INODE
             DO WHILE ( IN .GT. 0 )
               NPIV = NPIV + 1
               IN = FILS( IN )
             END DO
             IF ( TPN .eq. 1 ) THEN
               IF ( KEEP( 50 ) .eq. 0 ) THEN
                 COST( ICT ) = NPIV * ( 2 * NFR - NPIV )
               ELSE
                 COST( ICT ) = NPIV * NFR
               END IF
             ELSE
               KMAX  = CMUMPS_497( KEEP(10), NFR - NPIV)
               NSLAVES     = MAX( MIN((NFR-NPIV)/KMAX, SLAVEF-1),
     *                            1 )
               IF ( KEEP(50) .eq. 0 ) THEN
                 COST( ICT ) = NPIV * NFR - NPIV *
     *           ( ( NFR - NPIV ) / NSLAVES )
               ELSE
                 COST( ICT ) = NPIV * NPIV - NPIV *
     *           ( ( NFR - NPIV ) / NSLAVES )
               END IF
             END IF
          ELSE IF ( CMUMPS_283( INODE, PROCNODE, SLAVEF ) .OR.
     *              CMUMPS_167( INODE, PROCNODE, SLAVEF )  ) THEN
               NFR  = ND( INODE )
               NPIV = 0
               IN   = INODE
               DO WHILE ( IN .GT. 0 )
                 NPIV = NPIV + 1
                 IN = FILS( IN )
               END DO
               IF ( KEEP( 50 ) .eq. 0 ) THEN
                 MEMUSED(1+CMUMPS_275( INODE, PROCNODE, SLAVEF))
     *           = MEMUSED(1+CMUMPS_275( INODE, PROCNODE, SLAVEF))
     *           + NPIV * ( 2 * NFR - NPIV )
               ELSE
                 MEMUSED(1+CMUMPS_275( INODE, PROCNODE, SLAVEF))
     *           = MEMUSED(1+CMUMPS_275( INODE, PROCNODE, SLAVEF))
     *           + NPIV * NFR
               END IF
          END IF
        END IF
        END DO
        CALL CMUMPS_212( N, NNODE, NUMNODE, COST, MEMUSED,
     *                               PROCNODE, SLAVEF, NODEMAP )
        DO I = 1, NNODE
          INODE   = NUMNODE( I )
          PROCVAL = PROCNODE( INODE ) + ( NODEMAP( I ) - 1 )
     *            - CMUMPS_275( INODE, PROCNODE, SLAVEF )
          IN = INODE
          DO WHILE ( IN > 0 )
            PROCNODE( IN ) = PROCVAL
            IN = FILS( IN )
          END DO
        END DO
        DEALLOCATE( NUMNODE )
        DEALLOCATE( NODEMAP )
        DEALLOCATE( COST    )
        DEALLOCATE( MEMUSED )
      END IF
      RETURN
      END SUBROUTINE CMUMPS_21
      SUBROUTINE CMUMPS_212( N, NNODE, NUMNODE, COST,
     *       MEMUSED, PROCNODE, SLAVEF, NODEMAP )
      IMPLICIT NONE
      INTEGER N, NNODE, SLAVEF
      INTEGER PROCNODE( N )
      INTEGER NUMNODE( NNODE ), COST( NNODE ), NODEMAP( NNODE )
      INTEGER MEMUSED( SLAVEF )
      INTEGER IMAXLOAD, IMINLOAD, IPROC, INODE
      CALL CMUMPS_149
     *       ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      CALL CMUMPS_260( NNODE, NUMNODE, COST )
      DO INODE = 1, NNODE
        IF ( COST( INODE ) > 0 ) THEN
          IPROC = IMINLOAD
        ELSE
          IPROC = IMAXLOAD
        END IF
        MEMUSED( IPROC ) = MEMUSED( IPROC ) + COST( INODE )
        NODEMAP( INODE ) = IPROC
        CALL CMUMPS_149
     *       ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      END DO
      RETURN
      END SUBROUTINE CMUMPS_212
      SUBROUTINE CMUMPS_149
     *           ( MEMUSED, SLAVEF, IMINLOAD, IMAXLOAD )
      INTEGER SLAVEF, IMINLOAD, IMAXLOAD
      INTEGER MEMUSED( SLAVEF )
      INTEGER MINLOAD, MAXLOAD, IPROC
      MINLOAD = MEMUSED( 1 )
      MAXLOAD = MEMUSED( 1 )
      IMINLOAD = 1
      IMAXLOAD = 1
      DO IPROC = 2, SLAVEF
        IF ( MEMUSED( IPROC ) .GT. MAXLOAD ) THEN
          MAXLOAD  = MEMUSED( IPROC )
          IMAXLOAD = IPROC
        END IF
        IF ( MEMUSED( IPROC ) .LT. MINLOAD ) THEN
          MINLOAD  = MEMUSED( IPROC )
          IMINLOAD = IPROC
        END IF
      END DO  
      RETURN
      END SUBROUTINE CMUMPS_149
      SUBROUTINE CMUMPS_260( NNODE, NUMNODE, COST )
      IMPLICIT NONE
      INTEGER NNODE
      INTEGER NUMNODE( NNODE ), COST( NNODE )
      LOGICAL FIN
      INTEGER I, LOC_COST, LOC_NUM
      FIN = .FALSE.
      DO WHILE ( .NOT. FIN )
        FIN = .TRUE.
        DO I = 1, NNODE - 1
          IF ( ABS( COST( I ) ) .LT. ABS( COST( I + 1 ) ) ) THEN
            LOC_COST      = COST( I )
            COST( I )     = COST( I + 1 )
            COST( I + 1 ) = LOC_COST
            LOC_NUM          = NUMNODE( I )
            NUMNODE( I )     = NUMNODE( I + 1 )
            NUMNODE( I + 1 ) = LOC_NUM
            FIN = .FALSE.
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_260
      SUBROUTINE CMUMPS_100
     *( MYID, COMM, KEEP,KEEP8, INFO, INFOG, RINFO, RINFOG, ICNTL )
      IMPLICIT NONE
      INTEGER COMM, MYID, KEEP(500), INFO(40), ICNTL(40), INFOG(40)
      INTEGER*8 KEEP8(150)
      REAL RINFO(20), RINFOG(20)
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR_MPI, MPG
      PARAMETER( MASTER = 0 )
      IF ( KEEP(46) .eq. 0 .and. MYID .eq. MASTER ) THEN
        INFO(3)  = 0
        INFO(4)  = 0
        INFO(5)  = 0
        RINFO(1) = 0.0E0
      END IF
      MPG = ICNTL(3)
      CALL CMUMPS_621( INFO(3), INFOG(3), KEEP8(5), COMM )
      CALL MPI_REDUCE( INFO(3), INFOG(3), 1, MPI_INTEGER,
     *  MPI_SUM, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO(4), INFOG(4), 1, MPI_INTEGER,
     *  MPI_SUM, MASTER, COMM, IERR_MPI )
      IF (MYID.eq.MASTER) KEEP(126) = INFOG(4)
      CALL MPI_REDUCE( INFO(5), INFOG(5), 1, MPI_INTEGER,
     *  MPI_MAX, MASTER, COMM, IERR_MPI )
      IF (MYID.EQ.MASTER) KEEP(127) = INFOG(5)
      CALL MPI_REDUCE( RINFO(1), RINFOG(1), 1, MPI_REAL,
     *  MPI_SUM, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq.MASTER.and.MPG.GT.0) THEN
       WRITE(MPG, 99992) INFO(1), INFO(2), 
     *  KEEP8(1), KEEP8(5), INFOG(4),
     *  INFOG(5), KEEP(28), INFOG(7), KEEP(23), ICNTL(7), 
     *  KEEP(12), KEEP(56), KEEP(61), RINFOG(1)
       IF (KEEP(95).GT.1)             ! meaningful only with sym=2
     *                    WRITE(MPG, 99993) KEEP(95)  ! ICNTL(12)
       IF (KEEP(54).GT.0) WRITE(MPG, 99994) KEEP(54)  ! Distributed matrix
       IF (KEEP(60).GT.0) WRITE(MPG, 99995) KEEP(60)  ! Schur option
      ENDIF
      RETURN
99992 FORMAT(/'Leaving analysis phase with  ...'/
     1       'INFOG(1)                                       =',I12/
     2       'INFOG(2)                                       =',I12/
     3       ' -- (20) Number of entries in factors (estim.) =',I12/
     3       ' --  (3) Storage of factors  (REAL, estimated) =',I12/
     3       ' --  (4) Storage of factors  (INT , estimated) =',I12/
     5       ' --  (5) Maximum frontal size      (estimated) =',I12/
     6       ' --  (6) Number of nodes in the tree           =',I12/
     6       ' --  (7) Ordering option effectively used      =',I12/
     7       'ICNTL(6) Maximum transversal option            =',I12/
     8       'ICNTL(7) Pivot order option                    =',I12/
     8       'Percentage of memory relaxation (effective)    =',I12/
     8       'Number of level 2 nodes                        =',I12/
     8       'Number of splitted nodes                       =',I12/
     9   'RINFO(1) Operations during elimination (estim) =  ',1PD10.3)
99993 FORMAT('Ordering compressed/constrained (ICNTL(12))    =',I12)
99994 FORMAT('Distributed matrix entry format (ICNTL(18))    =',I12)
99995 FORMAT('Effective Schur option (ICNTL(19))             =',I12)
      END SUBROUTINE CMUMPS_100
      SUBROUTINE CMUMPS_23(N,IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER PE(N), LEN(N), 
     $        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     $        W(N)
      INTEGER IW(IWLEN), NV(N), NEXT(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X, 
     $        NPRINC
      INTEGER*8 HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      NPRINC = 0
      DO I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
      ENDDO
      DO I=1, N
        IF (LEN (I).GE.0) THEN
           DEGREE (I) = LEN (I)
           NPRINC = NPRINC + 1
        ELSE
           J        = -LEN (I)
           DEGREE (I) = - 1
           IF ( PE(I) .NE. 0 ) THEN
             LEN (I) = LEN(J)
           ELSE
             LEN (I) = 0
           ENDIF
           PE (I)   = -J
           NV (J)   = NV (J) + NV (I)
           NV (I)   = 0
           ELEN (I) = 0
        ENDIF
      ENDDO
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE IF ( DEG.EQ. 0) THEN
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG2)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
#if defined (NOAGG2)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = MIN (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_23
      SUBROUTINE CMUMPS_162(N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     $                   LISTVAR_SCHUR, SIZE_SCHUR)
      INTEGER SIZE_SCHUR
      INTEGER LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER N, IWLEN, PFREE, NCMPA
      INTEGER LEN(N), 
     $        ELEN(N), LAST(N), DEGREE(N), HEAD(N), 
     $        W(N), NEXT(N)
      INTEGER IW(IWLEN), NV(N), PE(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     $        NBFLAG, NREAL, LASTD, NELME
      INTEGER*8 HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      NBFLAG = 0
      LASTD  = 0
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE(I) = LEN(I)
   10 CONTINUE
      NBFLAG = SIZE_SCHUR
      DO K=1,SIZE_SCHUR
       I = LISTVAR_SCHUR(K)
       DEGREE(I) = N+1
       IF ((LEN(I) .EQ.0).OR.(LEN(I).EQ.-N-1)) THEN
         PE (I)     = 0
         LEN(I)     = 0
       ENDIF
       DEG = N
       IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
       ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
       ENDIF
      ENDDO
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N+1)  GOTO 20
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = N-NEL
   30 IF (NEL .LT. NREAL) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          write (*,*) ' Error 1 in HALO_AMD '
          NCMPA = -N
          GOTO 500
        ENDIF
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG3)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N+1) DEG = N+1
#if defined (NOAGG3)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N+1) THEN
                 DEG        = MIN (DEG, NLEFT)
                 DEGREE (I) = MIN (DEGREE (I), DEG)
            ENDIF
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.LE.N) THEN
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I) .LE. N) ) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
            DEGREE (I) = DEG
              ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, N
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           IF (ME.NE.LISTVAR_SCHUR(1)) THEN
             write(6,*) ' error 1 in CMUMPS_162 '
             write(6,*) ' wrong principal var for Schur !!'
             NCMPA = -N - 2
             CALL CMUMPS_ABORT()
           ENDIF
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         write(*,*) ' Error 2 in CMUMPS_162 NEL, N=', NEL,N
         NCMPA = -N - 1
         CALL CMUMPS_ABORT()
        ENDIF
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_162
      SUBROUTINE CMUMPS_97
     *           ( N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES, 
     *             KEEP, KEEP8, SPLITROOT, MP, LDIAG )
      IMPLICIT NONE
      INTEGER N, NSTEPS, NSLAVES, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N ),
     *        IPOOL( NSTEPS + 1 )
      LOGICAL SPLITROOT
      INTEGER MP, LDIAG  ! to control printing icntl(3,4)
      INTEGER INODE, DEPTH, I, IBEG, IEND, IIPOOL, NROOT
      INTEGER MAX_DEPTH, ISON, TOT_CUT, MAX_CUT, STRAT, K79
      INTEGER NFRONT, K82
      K79  = KEEP(79)
      K82  = ABS(KEEP(82))  
      STRAT=KEEP(62)
      IF (KEEP(210).EQ.1) THEN
        MAX_DEPTH = 2*NSLAVES*K82
        STRAT     = STRAT/4
      ELSE
        IF (( NSLAVES .eq. 1 ).AND. (.NOT. SPLITROOT) ) RETURN
        MAX_DEPTH = INT( LOG( REAL( NSLAVES - 1 ) ) 
     &                 / LOG(2.0E0) )
      ENDIF
      NROOT = 0
      DO INODE = 1, N
        IF ( FRERE(INODE) .eq. 0 ) THEN
          NROOT = NROOT + 1
          IPOOL( NROOT ) = INODE
        END IF
      END DO
      IBEG = 1
      IEND = NROOT
      IIPOOL   = NROOT + 1
      IF (SPLITROOT) MAX_DEPTH=1
      DO DEPTH = 1, MAX_DEPTH
        DO I = IBEG, IEND
          INODE = IPOOL( I )
          ISON = INODE
          DO WHILE ( ISON .GT. 0 )
            ISON = FILS( ISON )
          END DO
          ISON = - ISON
          DO WHILE ( ISON .GT. 0 )
            IPOOL( IIPOOL ) = ISON
            IIPOOL = IIPOOL + 1
            ISON = FRERE( ISON )
          END DO
        END DO
        IPOOL( IBEG ) = -IPOOL( IBEG )
        IBEG = IEND + 1
        IEND = IIPOOL - 1
      END DO
      IPOOL( IBEG ) = -IPOOL( IBEG )
      TOT_CUT = 0
      IF (SPLITROOT) THEN
        MAX_CUT = NROOT*MAX(K82,2)
        INODE = ABS(IPOOL(1))
        NFRONT = NFSIZ( INODE )
        K79     = MAX(REAL(NFRONT*NFRONT)/
     &            REAL((K82+1)*(K82+1)), REAL(1))
      ELSE
         MAX_CUT = 2 * NSLAVES
         IF (KEEP(210).EQ.1) THEN
            MAX_CUT = 4 * (MAX_CUT + 4)
         ENDIF
      ENDIF
      DEPTH   = -1
      DO I = 1, IIPOOL - 1
        INODE = IPOOL( I )
        IF ( INODE .LT. 0 ) THEN
          INODE = -INODE
          DEPTH = DEPTH + 1
        END IF
        CALL CMUMPS_313
     *           ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES,
     *             KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *             K79, SPLITROOT, MP, LDIAG )
        IF ( TOT_CUT > MAX_CUT )  EXIT
      END DO
      KEEP(61) = TOT_CUT
      RETURN
      END SUBROUTINE CMUMPS_97
      RECURSIVE SUBROUTINE CMUMPS_313
     * ( INODE, N, FRERE, FILS, NFSIZ, NSTEPS, NSLAVES, KEEP,KEEP8,
     *   TOT_CUT, STRAT, DEPTH, K79, SPLITROOT, MP, LDIAG )
      IMPLICIT NONE
      INTEGER INODE, N, NSTEPS, NSLAVES, KEEP(500), STRAT, 
     &        DEPTH, TOT_CUT, K79, MP, LDIAG
      INTEGER*8 KEEP8(150)
      INTEGER FRERE( N ), FILS( N ), NFSIZ( N )
      LOGICAL SPLITROOT
      INTEGER I, IN, NPIV, NFRONT, NSLAVES_ESTIM
      REAL WK_SLAVE, WK_MASTER
      INTEGER INODE_SON, INODE_FATH, IN_SON, IN_FATH, IN_GRANDFATH
      INTEGER NPIV_SON, NPIV_FATH, MAX_CUT
      INTEGER NCB, KAVERAGE, NSLAVESMIN, NSLAVESMAX
      INTEGER  CMUMPS_50,
     *         CMUMPS_52
      EXTERNAL  CMUMPS_50,
     *         CMUMPS_52
      IF  ( (KEEP(210).EQ.1.AND.KEEP(60).EQ.0) .OR.
     &       (SPLITROOT) ) THEN
        IF ( FRERE ( INODE ) .eq. 0 ) THEN 
          NFRONT = NFSIZ( INODE )
          NPIV = NFRONT
          NCB = 0
          IF (NFRONT*NFRONT.GT.K79) THEN 
           GOTO 333
          ELSE IF ( (SPLITROOT).AND.
     &      (MP.GT.0 .AND. LDIAG.GE.2)) THEN
             WRITE(MP,'(A,I11)') 
     &        ' Order of root node after splitting : ',
     &        NFRONT
          ENDIF
        ENDIF
      ENDIF
      IF ( FRERE ( INODE ) .eq. 0 ) RETURN
      NFRONT = NFSIZ( INODE )
      IN = INODE
      NPIV = 0
      DO WHILE( IN > 0 )
        IN = FILS( IN )
        NPIV = NPIV + 1
      END DO
      NCB = NFRONT - NPIV
      IF ( (NFRONT - (NPIV/2)) .LE. KEEP(9)) RETURN
      IF ((KEEP(50) == 0.and.NFRONT * NPIV > K79 ) .OR.
     &(KEEP(50) .NE.0.and.NPIV * NPIV > K79 )) GOTO 333
      IF (KEEP(210).EQ.1) THEN
        NSLAVESMIN    = 1   ! not used
        NSLAVESMAX    = 64  ! not used
        NSLAVES_ESTIM = 32+NSLAVES
      ELSE
        NSLAVESMIN = CMUMPS_50 
     *         ( NSLAVES, KEEP(48), KEEP(10), KEEP(50),
     *         NFRONT, NCB)
        NSLAVESMAX = CMUMPS_52 
     *        ( NSLAVES, KEEP(48), KEEP(10), KEEP(50),
     *          NFRONT, NCB)
        NSLAVES_ESTIM = MAX (1, 
     *   NINT( REAL(NSLAVESMAX-NSLAVESMIN)/REAL(3) )
     *                    )
        NSLAVES_ESTIM = MIN (NSLAVES_ESTIM, NSLAVES-1)
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
       WK_MASTER = 0.6667E0 * 
     *                REAL(NPIV)*REAL(NPIV)*REAL(NPIV) +
     *                REAL(NPIV)*REAL(NPIV)*REAL(NCB)
       WK_SLAVE  = REAL( NPIV ) * REAL( NCB ) *
     *         ( 2.0E0 * REAL(NFRONT) - REAL(NPIV) )
     *         / REAL(NSLAVES_ESTIM)
      ELSE
       WK_MASTER = REAL(NPIV)*REAL(NPIV)*REAL(NPIV) / REAL(3)
       WK_SLAVE  = 
     *           (REAL(NPIV)*REAL(NCB)*REAL(NFRONT)) 
     *           / REAL(NSLAVES_ESTIM)
      ENDIF
      IF (KEEP(210).EQ.1) THEN
        IF ( ( 100 + STRAT )
     &        * WK_SLAVE / 100 .GE. WK_MASTER ) RETURN
      ELSE
        IF ( ( 100 + STRAT * MAX( DEPTH-1, 1 ) )
     &        * WK_SLAVE / 100 .GE. WK_MASTER ) RETURN
      ENDIF
 333  CONTINUE
      IF (NPIV .LE. 1 ) RETURN
       NSTEPS  = NSTEPS + 1
       TOT_CUT = TOT_CUT + 1
       NPIV_SON  = MAX(NPIV/2,1)
       NPIV_FATH = NPIV - NPIV_SON
       INODE_SON = INODE
       IN_SON = INODE
       DO I = 1, NPIV_SON - 1
         IN_SON = FILS( IN_SON )
       END DO
       INODE_FATH = FILS( IN_SON )
       IF ( INODE_FATH .LT. 0 ) THEN
       write(*,*) 'Error: INODE_FATH < 0 ', INODE_FATH
       END IF
       IN_FATH = INODE_FATH
       DO WHILE ( FILS( IN_FATH ) > 0 )
         IN_FATH = FILS( IN_FATH )
       END DO
       FRERE( INODE_FATH ) = FRERE( INODE_SON )
       FRERE( INODE_SON  ) = - INODE_FATH
       FILS ( IN_SON     ) = FILS( IN_FATH )
       FILS ( IN_FATH    ) = - INODE_SON
       IN = FRERE( INODE_FATH )
       DO WHILE ( IN > 0 )
           IN = FRERE( IN )
       END DO
       IF ( IN .eq. 0 )  GO TO 10
       IN = -IN
       DO WHILE ( FILS( IN ) > 0 )
           IN = FILS( IN )
       END DO
       IN_GRANDFATH = IN
       IF ( FILS( IN_GRANDFATH ) .eq. - INODE_SON ) THEN
           FILS( IN_GRANDFATH ) = -INODE_FATH
       ELSE
           IN = IN_GRANDFATH
           IN = - FILS ( IN )
           DO WHILE ( FRERE( IN ) > 0 )
             IF ( FRERE( IN ) .eq. INODE_SON ) THEN
               FRERE( IN ) = INODE_FATH
               GOTO 10
             END IF
             IN = FRERE( IN )
           END DO
           WRITE(*,*) 'ERROR 2 in SPLIT NODE',
     *          IN_GRANDFATH, IN, FRERE(IN)
       END IF
 10    CONTINUE
       NFSIZ(INODE_SON) = NFRONT
       NFSIZ(INODE_FATH) = NFRONT - NPIV_SON
       KEEP(2) = MAX( KEEP(2), NFRONT - NPIV_SON )
      CALL CMUMPS_313
     *  ( INODE_FATH, N, FRERE, FILS, NFSIZ, NSTEPS,
     *   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *   K79, SPLITROOT, MP, LDIAG )
      IF (.NOT. SPLITROOT) THEN
        CALL CMUMPS_313
     *   ( INODE_SON, N, FRERE, FILS, NFSIZ, NSTEPS,
     *   NSLAVES, KEEP,KEEP8, TOT_CUT, STRAT, DEPTH, 
     *   K79, SPLITROOT, MP, LDIAG )
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_313
      SUBROUTINE CMUMPS_351
     * (N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     * IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM, MedDens, NBQD, AvgDens)
      INTEGER N,NZ,LW,IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER symmetry, SYM
      INTEGER MedDens, NBQD, AvgDens
      INTEGER ICNTL(40)
      INTEGER  IRN(NZ), ICN(NZ) 
      INTEGER LEN(N)
      INTEGER IPE(N+1)
      INTEGER FLAG(N), IW(LW)
      INTEGER IQ(N)
      INTEGER MP, MPG
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR, THRESH
      INTEGER NZOFFA, NDIAGA
      REAL RSYM
      INTRINSIC NINT
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      DO 10 I=1,N
        IPE(I) = 0
   10 CONTINUE
      DO 50 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
   50 CONTINUE
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (MOD(IFLAG,2) .EQ. 0) IFLAG = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (MOD(K,10).GT.3 .OR. MOD(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (MOD(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I) 
  110   CONTINUE
      ENDIF
      LAST = MAX(IPE(N)+IQ(N)-1,IQ(N))
      FLAG(1:N) = 0
      IPE(1:N)  = IQ(1:N)
      IW(1:LAST) = 0
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I) 
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J) 
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         IF (IQ(I).EQ.0) THEN
             LEN(I) = 0
            IPE(I) = IWFR
            GOTO 280
         ENDIF
         K1 = IPE(I) 
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L 
  280  CONTINUE
      ENDIF
      IPE(N+1) = IPE(N) + LEN(N)
      IWFR = IPE(N+1)
      IF (SYM.EQ.0) THEN
      RSYM =  REAL(NDIAGA+2*NZOFFA - (IWFR-1))/
     *            REAL(NZOFFA+NDIAGA) 
      symmetry = NINT (100.0*RSYM)
         IF (MPG .GT. 0)
     &  write(MPG,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,I5)') 
     &  ' ... Structural symmetry (in percent)=', symmetry
      ELSE
       symmetry = 100
      ENDIF
      AvgDens = NINT(REAL(IWFR-1)/REAL(N))
      THRESH  = AvgDens*50 - AvgDens/10 + 1
      NBQD    = 0
      IF (N.GT.2) THEN
        IQ(1:N) = 0
        DO I= 1, N
          K = MAX(LEN(I),1)
          IQ(K) = IQ(K) + 1
          IF (K.GT.THRESH) NBQD = NBQD+1
        ENDDO
        K = 0
        MedDens = 0
        DO WHILE (K .LT. (N/2))
         MedDens = MedDens + 1
         K       = K+IQ(MedDens)
        ENDDO
      ELSE
        MedDens = AvgDens
      ENDIF
         IF (MPG .GT. 0)
     &  write(MPG,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
        IF (MP.GT.0 .AND. MPG.NE.MP)
     &  write(MP,'(A,3I5)') 
     &  ' Density: NBdense, Average, Median   =',
     &  NBQD, AvgDens, MedDens
      RETURN
99999 FORMAT (/'*** Warning message from analysis routine ***')
      END SUBROUTINE CMUMPS_351
      SUBROUTINE CMUMPS_337(N, NBBUCK, 
     $                   IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD)
      IMPLICIT NONE
      INTEGER N, IWLEN, PFREE, LEN(N), 
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), NEXT(N),
     $        W(N) 
      INTEGER PE(N), IW(IWLEN), NV(N)
      INTEGER NBBUCK
      INTEGER HEAD(0:NBBUCK+1), WF(N)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     $        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER*8 HASH, HMOD
       REAL    RMF, RMF1 
       REAL    dummy
       INTEGER idummy
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD, HUGE
      INTEGER TOTEL
      LOGICAL COMPRESS
      idummy = HUGE(idummy) - 1
      dummy = REAL(idummy)
      N2 = -NBBUCK-1
      PAS = MAX((N/8), 1)
      WFLG = 2
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, NBBUCK-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
        LAST(I) = 0
        W(I) = 1
        ELEN (I) = 0
   10 CONTINUE
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF(COMPRESS) THEN
         TOTEL = 0
         DO I=1,N
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-N-1) THEN
                  LEN (I)    = 0
                  PE (I)     = 0
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
            ELSE
               TOTEL = TOTEL + NV(I)
               DEGREE(I) = 0
               DO J= PE(I) , PE(I)+LEN(I)-1 
                  DEGREE(I) = DEGREE(I) + NV(IW(J))
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-N-1) THEN
                  LEN (I)    = 0
                  PE (I)     = 0
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
            ELSE
               DEGREE (I) = LEN (I)
            ENDIF
         ENDDO
         TOTEL = N - NBFLAG
      ENDIF
#if(SP_DEBUG >= 2)
         write(*,*) 'deginit fin totel =',TOTEL
#endif
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
        IF (DEG .GT. 0) THEN
          WF(I) = DEG
           IF (DEG.GT.N) THEN
            DEG = MIN(((DEG-N)/PAS) + N , NBBUCK)
           ENDIF
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = TOTEL-NEL
#if(SP_DEBUG >= 1)
      IF(NREAL .GT. TOTEL) THEN
         CALL CMUMPS_ABORT()
      ENDIF
#endif
   30 IF (NEL .LT. TOTEL) THEN
        DO 40 DEG = MINDEG, NBBUCK
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          NCMPA = -N
          CALL CMUMPS_ABORT()
        ENDIF
       IF (DEG.GT.N) THEN
         J = NEXT(ME)
         K = WF(ME)
   55    CONTINUE
         IF (J.GT.0) THEN
          IF (WF(J).LT.K) THEN
           ME = J
           K  = WF(ME)
          ENDIF
          J= NEXT(J)
          GOTO 55
         ENDIF
         ILAST = LAST(ME)
         INEXT = NEXT(ME)
         IF (INEXT .NE. 0) LAST (INEXT) = ILAST
         IF (ILAST .NE. 0) THEN
           NEXT (ILAST) = INEXT
         ELSE
           HEAD (DEG) = INEXT
         ENDIF
        ELSE
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).NE.N2) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                IF (WF(I).GT.N) THEN
                 DEG = MIN(((WF(I)-N)/PAS) + N , NBBUCK)
                ELSE
                 DEG = WF(I)
                ENDIF
                HEAD (DEG) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.N) THEN
                   DEG = MIN(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N2) DEG = N2
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N2) THEN
                 IF ( DEGREE (I).LT.DEG ) THEN
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
            WF(I)      = WF4 + 2*NVI*WF3
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.NE.N2) THEN
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                WF(I)  = MAX(WF(I),WF(J))
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).NE.N2) THEN
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
              DEG = DEGREE(I)
              RMF1  = REAL(DEG)*REAL( (DEG-1) + 2*DEGME )
     &              - REAL(WF(I))
              DEGREE(I) = NLEFT - NVI
              DEG       = DEGREE(I) 
              RMF = REAL(DEG)*REAL(DEG-1) 
     &         -  REAL(DEGME-NVI)*REAL(DEGME-NVI-1)
              RMF = MIN(RMF, RMF1)
            ELSE 
              DEG = DEGREE(I)
              DEGREE(I) = DEGREE (I) + DEGME - NVI
              RMF  = REAL(DEG)*REAL( (DEG-1) + 2*DEGME ) 
     &              - REAL(WF(I))
            ENDIF
            RMF =  RMF / REAL(NVI+1)
            IF (RMF.LT.dummy) THEN
             WF(I) = INT ( ANINT( RMF ))
            ELSEIF (RMF / REAL(N) .LT. dummy) THEN 
             WF(I) = INT ( ANINT( RMF/REAL(N) ))
            ELSE
             WF(I) = idummy
            ENDIF
            WF(I) = MAX(1,WF(I))
            DEG = WF(I)
            IF (DEG.GT.N) THEN
              DEG = MIN(((DEG-N)/PAS) + N , NBBUCK)
            ENDIF
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
              ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N2) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         NCMPA = -N - 1
         GOTO 500
        ENDIF
      ENDIF
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF(.NOT. COMPRESS) THEN
         DO 300 I = 1, N
            K = ABS (ELEN (I))
            LAST (K) = I
            ELEN (I) = K
 300     CONTINUE
      ENDIF
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_337
      SUBROUTINE CMUMPS_421 
     $                (TOTEL, IVersion, THRESH, NDENSE, 
     $                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     $                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W) 
      INTEGER TOTEL
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N), 
     $     ELEN(N), NCMPA, DEGREE(N), 
     *     LAST(TOTEL), HEAD(TOTEL), NEXT(N), 
     $     W(N)
      INTEGER NDENSE(N)
      INTEGER IVersion, THRESH
      INTEGER THRESM, MINDEN, MAXDEN, NDME
      INTEGER NBD,NBED, NBDM, LASTD, NELME, DEG1
      LOGICAL IDENSE
      REAL RELDEN
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER*8 HASH, HMOD
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD
      LOGICAL COMPRESS
      IF (THRESH.GT.0) THEN 
         THRESM  = MIN(N,THRESH)
         DO I=1,N
             THRESM = MAX(THRESM, LEN(I))
          ENDDO
           RELDEN = REAL(PFREE-1)/REAL(N)
         THRESM =  INT(RELDEN)*10 + (THRESM-INT(RELDEN))/10 + 1
      ELSE
         THRESM = TOTEL
      ENDIF
      IF (THRESM.GE.0) THEN
       IF ((THRESM.GT.TOTEL).OR.(THRESM.LT.2)) THEN 
          THRESM = TOTEL
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        W (I) = 1
        ELEN (I) = 0
   10 CONTINUE
      HEAD(N:TOTEL) = 0
      LAST(N:TOTEL) = 0
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF(COMPRESS) THEN
         DO I=1,N
            DEGREE(I) = 0
            DO J= PE(I) , PE(I)+LEN(I)-1 
               DEGREE(I) = DEGREE(I) + NV(IW(J))
            ENDDO
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            DEGREE (I) = LEN (I)
         ENDDO
      ENDIF
#if(SP_DEBUG >= 2)
         write(*,*) 'deginit fin totel =',TOTEL
#endif
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+NV(I).GE.THRESM) ) THEN
               NBD = NBD+1
               IF (DEG+NV(I).NE.TOTEL-NEL) THEN
                  DEGREE(I) = DEGREE(I)+TOTEL+1
                  DEG = TOTEL
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  NEXT (I) = INEXT
                  HEAD (DEG) = I 
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) LASTD=I
               ELSE
                  NBED = NBED+1
                  DEGREE(I) = TOTEL+1
                  DEG = TOTEL
                  IF (LASTD.EQ.0) THEN
                     LASTD     = I 
                     HEAD(DEG) = I
                     NEXT(I)   = 0 
                     LAST(I)   = 0
                  ELSE
                     NEXT(LASTD) = I
                     LAST(I)     = LASTD
                     LASTD       = I
                     NEXT(I)     = 0
                  ENDIF
               ENDIF
            ELSE
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ELSE
            NEL = NEL + NV(I)
            ELEN (I) = -NEL
            PE (I) = 0
            W (I) = 0
         ENDIF
 20   CONTINUE
#if (SP_DEBUG>=2)
      DO I=1,N
         IF(NEXT(I) .GT. N .AND. ELEN(I) .GT. 0) THEN
            WRITE(*,*) 'ERROR I, NEXT(I)',I,NEXT(I)
            CALL CMUMPS_ABORT()
         ENDIF
      ENDDO
#endif
          IF (NBD.EQ.0) THRESM = TOTEL
          NLEFT = TOTEL - NEL
 30       IF (NEL .LT. TOTEL) THEN
        DO 40 DEG = MINDEG, TOTEL
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (DEG.LT.TOTEL)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          NBDM = MAX(NBDM,NBD)
          IF (DEGREE(ME).GT.TOTEL+1) THEN
            MINDEN = NBD
            MAXDEN = 0
            IF (WFLG+NBD+1 .LE. WFLG) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            MINDEN = MIN (MINDEN, NDME)
            MAXDEN = MAX (MAXDEN, NDME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            IF (IVersion.EQ.1) THEN
              DEG = MAX (DEGREE(ME)-(TOTEL+1), 1)
            ELSE
              DEG = NDENSE(ME)
            ENDIF
            DEGREE(ME) = DEG
            MINDEG = MIN(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
#if(SP_DEBUG >= 2)
            IF(JNEXT .GT. N .OR. JNEXT .LT. 0) THEN
               WRITE(*,*) 'pb JNEXT',JNEXT
               CALL CMUMPS_ABORT()
            ENDIF
#endif
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(TOTEL+1) ) GOTO 51
            ENDIF
            HEAD (TOTEL) = ME
            IF (IVersion .EQ.1 ) THEN
             THRESM = TOTEL
            ELSE
             THRESM=MAX(THRESM*2,MINDEN+(MAXDEN-MINDEN)/2)
             THRESM = MIN(THRESM,NBD)
             IF (THRESM.GE.NBD) THRESM=TOTEL
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.TOTEL+1) THEN
           IF (NBD.NE.NBED) THEN
            write(6,*) ' Internal ERROR quasi dense rows remains'
            CALL CMUMPS_ABORT()
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.TOTEL+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.TOTEL) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL CMUMPS_ABORT()
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.TOTEL) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
#if(SP_DEBUG >= 2)
                 IF(INEXT .GT. N .OR. INEXT .LT. 0) THEN
                    WRITE(*,*) 'pb INEXT',INEXT
                    CALL CMUMPS_ABORT()
                 ENDIF
#endif
                 NEXT (ILAST) = INEXT
              ELSE
                 HEAD (DEGREE (I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.TOTEL) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
#if(SP_DEBUG >= 2)
                   IF(INEXT .GT. N .OR. INEXT .LT. 0) THEN
                      WRITE(*,*) 'pb INEXT 2',INEXT
                      CALL CMUMPS_ABORT()
                   ENDIF
#endif
                   NEXT (ILAST) = INEXT
                ELSE
                   HEAD (DEGREE (I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#if defined (NOAGG5)
            ELSE IF (DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
            ELSE IF ((DEXT .EQ. 0) .AND.
     &                (NDENSE(ME).EQ.NBD)) THEN
                PE (E) = -ME
                W (E)  = 0
            ELSE IF (DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + E
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.TOTEL) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
#if defined (NOAGG5)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF ((DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD)) THEN
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = MIN (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
#if(SP_DEBUG >= 2)
               IF(-J .GT. N .OR. -J .LT. 0) THEN
                  WRITE(*,*) 'pb next J',J
                  CALL CMUMPS_ABORT()
               ENDIF
#endif
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
#if(SP_DEBUG >= 2)
               IF(LAST(J) .GT. N .OR. LAST(J) .LT. 0) THEN
                  WRITE(*,*) 'pb LAST(J)',LAST(J)
                  CALL CMUMPS_ABORT()
               ENDIF
#endif
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.TOTEL) ) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
#if(SP_DEBUG >= 2)
                IF(J .GT. N .OR. J .LT. 0) THEN
                   WRITE(*,*) 'pb3 J',J
                   CALL CMUMPS_ABORT()
                ENDIF
#endif
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.TOTEL) THEN
            DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
       IF ( (Iversion .NE. 1).AND. (THRESM.GE.0)) THEN
          IF (DEG+NVI .GE. THRESM) THEN
             IF (THRESM.EQ.TOTEL) THEN
                IF ((ELEN(I).LE.2) .AND. ((DEG+NVI).EQ.NLEFT) ) THEN
                   DEGREE(I) = TOTEL+1
                   IDENSE = .TRUE.
                ENDIF
             ELSE
                IDENSE = .TRUE.
                IF ((ELEN(I).LE.2).AND.((DEG+NVI).EQ.NLEFT) ) THEN
                   DEGREE(I) = TOTEL+1
                ELSE
                   DEGREE(I) = TOTEL+1+DEGREE(I)
                ENDIF
             ENDIF
          ENDIF
          IF (IDENSE) THEN
             P1 = PE(I)
             P2 = P1 + ELEN(I) - 1
             IF (P2.GE.P1) THEN
                DO 264 PJ=P1,P2
                   E= IW(PJ)
                   NDENSE (E) = NDENSE(E) + NVI
 264            CONTINUE
             ENDIF
             NBD = NBD+NVI
             DEG = TOTEL
             IF (DEGREE(I).EQ.TOTEL+1) THEN
                NBED = NBED +NVI
                IF (LASTD.EQ.0) THEN
                   LASTD     = I
                   HEAD(DEG) = I
                   NEXT(I)   = 0
                   LAST(I)   = 0
                ELSE
#if(SP_DEBUG >= 2)
                   IF(I .GT. N .OR. I .LT. 0) THEN
                      WRITE(*,*) 'pb4 I',I
                      CALL CMUMPS_ABORT()
                   ENDIF
#endif
                   NEXT(LASTD) = I
                   LAST(I)     = LASTD
                   LASTD       = I
                   NEXT(I)     = 0
                ENDIF
             ELSE
                INEXT = HEAD(DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
#if(SP_DEBUG >= 2)
                IF(INEXT .GT. N .OR. INEXT .LT. 0) THEN
                   WRITE(*,*) 'pb5 INEXT',INEXT
                   CALL CMUMPS_ABORT()
                ENDIF
#endif
                NEXT (I) = INEXT
                HEAD (DEG) = I
                LAST(I)    = 0
                IF (LASTD.EQ.0) LASTD=I
             ENDIF
          ENDIF
       ENDIF
       IF (.NOT.IDENSE) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
#if(SP_DEBUG >= 2)
          IF(INEXT .GT. N .OR. INEXT .LT. 0) THEN
             WRITE(*,*) 'pb7 INEXT,N,DEG',INEXT,N,DEG
             CALL CMUMPS_ABORT()
          ENDIF
#endif
          NEXT (I) = INEXT
          LAST (I) = 0
          HEAD (DEG) = I
       ENDIF
       MINDEG = MIN (MINDEG, DEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF(.NOT. COMPRESS) THEN
         DO 300 I = 1, N
            K = ABS (ELEN (I))
            LAST (K) = I
            ELEN (I) = K
 300     CONTINUE
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_421
      SUBROUTINE CMUMPS_560(N, NBBUCK, 
     $     IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $     LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD,
     *     CONSTRAINT,THESON )
      IMPLICIT NONE
      INTEGER N, IWLEN, PFREE, LEN(N), 
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), NEXT(N),
     $        W(N) 
      INTEGER PE(N), IW(IWLEN), NV(N)
      INTEGER NBBUCK
      INTEGER HEAD(0:NBBUCK+1), WF(N)
      INTEGER CONSTRAINT(N),THESON(N)
      INTEGER PREV,TOTO
      LOGICAL TROUVE
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     $        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER*8 HASH, HMOD
       REAL    RMF, RMF1 
       REAL    dummy
       INTEGER idummy
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
      INTRINSIC MAX, MIN, MOD, HUGE
      INTEGER TOTEL
      idummy = HUGE(idummy) - 1
      dummy = REAL(idummy)
      N2 = -NBBUCK-1
      PAS = MAX((N/8), 1)
      WFLG = 2
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, NBBUCK-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
         THESON(I) = 0
         LAST (I) = 0
         W (I) = 1
         ELEN (I) = 0
   10 CONTINUE
      TOTEL = 0
      DO I=1,N
         IF (LEN(I).LT.0) THEN 
            DEGREE (I) = N2
            NBFLAG     = NBFLAG +1
            IF (LEN(I).EQ.-N-1) THEN
               LEN (I)    = 0
               PE (I)     = 0
            ELSE
               LEN (I)    = - LEN(I)
            ENDIF
         ELSE
            TOTEL = TOTEL + NV(I)
            DEGREE(I) = 0
            DO J= PE(I) , PE(I)+LEN(I)-1 
               DEGREE(I) = DEGREE(I) + NV(IW(J))
            ENDDO
         ENDIF
      ENDDO
#if(SP_DEBUG >= 2)
         write(*,*) 'deginit fin totel =',TOTEL
#endif
      NREAL = N - NBFLAG
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
        IF (DEG .GT. 0) THEN
           WF(I) = DEG
           IF (DEG.GT.N) THEN
            DEG = MIN(((DEG-N)/PAS) + N , NBBUCK)
           ENDIF
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
           NEL = NEL + NV(I)
           ELEN (I) = -NEL
           PE (I) = 0
           W (I) = 0
        ENDIF
   20 CONTINUE
      NLEFT = TOTEL-NEL
#if(SP_DEBUG >= 2)
      WRITE(*,*) 'before 30, NEL = ',NEL
#endif
#if(SP_DEBUG >= 1)
      IF(NREAL .GT. TOTEL) THEN
         CALL CMUMPS_ABORT()
      ENDIF
#endif
 30   IF (NEL .LT. TOTEL) THEN
         DO 40 DEG = MINDEG, NBBUCK
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (ME.LE.0) THEN
            NCMPA = -N
            CALL CMUMPS_ABORT()
         ENDIF
         IF (DEG.GT.N) THEN
            J = NEXT(ME)
            K = WF(ME)
            IF(CONSTRAINT(ME) .LT. 0) THEN
               K = -1
            ENDIF
 55         CONTINUE
            IF (J.GT.0) THEN
               IF(CONSTRAINT(J) .GE. 0) THEN
                  IF (WF(J).LT.K .OR. K .LT. 0) THEN
                     ME = J
                     K  = WF(ME)
                  ENDIF
               ENDIF
               J= NEXT(J)
               GOTO 55
            ENDIF
            ILAST = LAST(ME)
            INEXT = NEXT(ME)
            IF (INEXT .NE. 0) LAST (INEXT) = ILAST
            IF (ILAST .NE. 0) THEN
               NEXT (ILAST) = INEXT
            ELSE
               HEAD (DEG) = INEXT
            ENDIF
         ELSE
            IF(CONSTRAINT(ME) .GE. 0) GOTO 59
 56         CONTINUE
            IF(NEXT(ME) .NE. 0) THEN
               ME = NEXT(ME)
               IF(CONSTRAINT(ME) .GE. 0) THEN
                  GOTO 59
               ELSE
                  GOTO 56
               ENDIF
            ELSE
 57            DEG = DEG+1
               ME = HEAD(DEG) 
               IF(ME .GT. 0) THEN
                  IF(CONSTRAINT(ME) .GE. 0) THEN
                     GOTO 59
                  ELSE
                     GOTO 56
                  ENDIF
               ELSE
                  GOTO 57
               ENDIF
            ENDIF   
 59         PREV = LAST (ME)
            INEXT = NEXT (ME)
            IF(PREV .NE. 0) THEN
               NEXT(PREV) = INEXT
            ELSE
               HEAD (DEG) = INEXT
            ENDIF
            IF (INEXT .NE. 0) LAST (INEXT) = PREV
         ENDIF
         TOTO = ME
 5910    IF(TOTO .NE. 0) THEN
            J = CONSTRAINT(TOTO)
            IF(J .GT. 0) THEN
               CONSTRAINT(J) = 0
            ENDIF
            TOTO = THESON(TOTO)
            GOTO 5910
         ENDIF
            ELENME = ELEN (ME)
            ELEN (ME) = - (NEL + 1)
            NVPIV = NV (ME)
#if(SP_DEBUG >= 4)
            TOTO = ME
            I = 0
            DO 
               IF(TOTO .LE. 0) EXIT
               TOTO = THESON(TOTO)
               I = I+1
            ENDDO
            IF(I .NE. NVPIV) THEN
               WRITE(*,*) 'ME,NEL,I,NVPIV',ME,NEL,I,NVPIV
               CALL CMUMPS_ABORT()
            ENDIF
#endif
            NEL = NEL + NVPIV
            NV (ME) = -NVPIV
            DEGME = 0
            IF (ELENME .EQ. 0) THEN
               PME1 = PE (ME)
               PME2 = PME1 - 1
               DO 60 P = PME1, PME1 + LEN (ME) - 1
                  I = IW (P)
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
                     DEGME = DEGME + NVI
                     NV (I) = -NVI
                     PME2 = PME2 + 1
                     IW (PME2) = I
                     IF (DEGREE(I).NE.N2) THEN
                        ILAST = LAST (I)
                        INEXT = NEXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           NEXT (ILAST) = INEXT
                        ELSE
                           IF (WF(I).GT.N) THEN
                              DEG = MIN(((WF(I)-N)/PAS) + N , NBBUCK)
                           ELSE
                              DEG = WF(I)
                           ENDIF
                           HEAD (DEG) = INEXT
                        ENDIF
                     ENDIF
                  ENDIF
 60            CONTINUE
               NEWMEM = 0
            ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.N) THEN
                   DEG = MIN(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
#if (SP_DEBUG >= 4 )
             IF(ME .GT. N .OR. ME .LE. 0) THEN
                WRITE(*,*) 'pb1111111111111111111111'
                CALL CMUMPS_ABORT()
             ENDIF
#endif
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
#if (SP_DEBUG >= 4 )
             IF(ME .GT. N .OR. ME .LE. 0) THEN
                WRITE(*,*) 'pb22222222222222222222222'
                CALL CMUMPS_ABORT()
             ENDIF
#endif
              PE (E) = -ME
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (DEGREE(I).EQ.N2) DEG = N2
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
#if (SP_DEBUG >= 4 )
             IF(ME .GT. N .OR. ME .LE. 0) THEN
                WRITE(*,*) 'pb3333333333333333333333'
                CALL CMUMPS_ABORT()
             ENDIF
#endif
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            IF (DEGREE(I).NE.N2) THEN
                 IF ( DEGREE (I).LT.DEG ) THEN
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
            WF(I)      = WF4 + 2*NVI*WF3
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            IF (DEG.NE.N2) THEN
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                 IF(CONSTRAINT(J) .LT. 0
     *                .AND. CONSTRAINT(I) .LT. 0) THEN
                    GOTO 240
                 ENDIF
                 IF(CONSTRAINT(I) .GE. 0) THEN
                    IF(CONSTRAINT(J) .LT. 0) THEN
                       TOTO = I
 221                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. J) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 221
                       ENDIF
                    ELSE
                       GOTO 225
                    ENDIF
                 ELSE
                    IF(CONSTRAINT(J) .GE. 0) THEN
                       TOTO = J
 222                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. I) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 222
                       ENDIF
                    ENDIF
                 ENDIF
                 GOTO 240
 225             CONTINUE
                 IF (LEN (J) .NE. LN) GO TO 240
                 IF (ELEN (J) .NE. ELN) GO TO 240
                 DO 230 P = PE (J) + 1, PE (J) + LN - 1
                    IF (W (IW (P)) .NE. WFLG) GO TO 240
 230             CONTINUE
                 TOTO = I
 231             IF(THESON(TOTO) .NE. 0) THEN
                    TOTO = THESON(TOTO)
                    GOTO 231
                 ENDIF
                 THESON(TOTO) = J
                 IF(CONSTRAINT(I) .LT. 0) THEN
                    CONSTRAINT(I) = 0
                 ENDIF
#if (SP_DEBUG >= 4 )
                 IF(I .GT. N .OR. I .LE. 0) THEN
                    WRITE(*,*) 'pb444444444444444444'
                    CALL CMUMPS_ABORT()
                 ENDIF
#endif
                 PE (J) = -I
                 WF(I)  = MAX(WF(I),WF(J))
                 NV (I) = NV (I) + NV (J)
                 NV (J) = 0
                 ELEN (J) = 0
                 J = NEXT (J)
                 NEXT (JLAST) = J
                 GO TO 220
 240             CONTINUE
                 JLAST = J
                 J = NEXT (J)
                 GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
           ENDIF
          ENDIF
 250   CONTINUE
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
           I = IW (PME)
           NVI = -NV (I)
           IF (NVI .GT. 0) THEN
              NV (I) = NVI
              IF (DEGREE(I).NE.N2) THEN
                 DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
                 IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
                  DEG = DEGREE(I)
                  RMF1  = REAL(DEG)*REAL( (DEG-1) + 2*DEGME )
     &                 - REAL(WF(I))
                  DEGREE(I) = NLEFT - NVI
                  DEG       = DEGREE(I) 
                  RMF = REAL(DEG)*REAL(DEG-1) 
     &                 -  REAL(DEGME-NVI)*REAL(DEGME-NVI-1)
                  RMF = MIN(RMF, RMF1)
               ELSE 
                  DEG = DEGREE(I)
                  DEGREE(I) = DEGREE (I) + DEGME - NVI
                  RMF  = REAL(DEG)*REAL( (DEG-1) + 2*DEGME ) 
     &                 - REAL(WF(I))
               ENDIF
               RMF =  RMF / REAL(NVI+1)
               IF (RMF.LT.dummy) THEN
                  WF(I) = INT ( ANINT( RMF ))
               ELSEIF (RMF / REAL(N) .LT. dummy) THEN 
                  WF(I) = INT ( ANINT( RMF/REAL(N) ))
               ELSE
                  WF(I) = idummy
               ENDIF
               WF(I) = MAX(1,WF(I))
               DEG = WF(I)
               IF (DEG.GT.N) THEN
                  DEG = MIN(((DEG-N)/PAS) + N , NBBUCK)
               ENDIF
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               LAST (I) = 0
               HEAD (DEG) = I
               MINDEG = MIN (MINDEG, DEG)
            ENDIF
            IW (P) = I
            P = P + 1
         ENDIF
 260  CONTINUE
      NV (ME) = NVPIV + DEGME
      LEN (ME) = P - PME1
      IF (LEN (ME) .EQ. 0) THEN
         PE (ME) = 0
         W (ME) = 0
      ENDIF
      IF (NEWMEM .NE. 0) THEN
         PFREE = P
         MEM = MEM - NEWMEM + LEN (ME)
      ENDIF
      GO TO 30
      ENDIF
#if(SP_DEBUG >= 2)
      WRITE(*,*) '30 END, NEL =',NEL
#endif
      IF (NEL.LT.N) THEN 
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
#if (SP_DEBUG >= 4 )
             IF(ME .GT. N .OR. ME .LE. 0) THEN
                WRITE(*,*) 'pb55555555555555555555'
                CALL CMUMPS_ABORT()
             ENDIF
#endif
             PE(X) = -ME
            ELSEIF (DEGREE(X).EQ.N2) THEN
             NEL   = NEL + NV(X)
#if (SP_DEBUG >= 4 )
             IF(ME .GT. N .OR. ME .LE. 0) THEN
                WRITE(*,*) 'pb66666666666666666666666'
                CALL CMUMPS_ABORT()
             ENDIF
#endif
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
           ENDDO
           ELEN(ME) = NELME
           NV(ME)   = N-NREAL
           PE(ME)   = 0
        IF (NEL.NE.N) THEN
         NCMPA = -N - 1
         GOTO 500
        ENDIF
      ENDIF
      DO 290 I = 1, N
         IF (ELEN (I) .EQ. 0) THEN
            J = -PE (I)
 270        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               J = -PE (J)
               GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
 280        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               JNEXT = -PE (J)
               PE (J) = -E
               IF (ELEN (J) .EQ. 0) THEN
                  ELEN (J) = K
                  K = K + 1
               ENDIF
               J = JNEXT
               GO TO 280
            ENDIF
            ELEN (E) = -K
         ENDIF
 290  CONTINUE
#if (SP_DEBUG >= 10)
      write(*,*) 'PE',PE
      write(*,*) 'NV',NV
#endif
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_560
      SUBROUTINE CMUMPS_422 
     $                ( THRESH, NDENSE, 
     $                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     $                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     $                 PERM, LISTVAR_SCHUR, SIZE_SCHUR, AGG6 ) 
      IMPLICIT NONE
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N), 
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N), 
     $        W(N), SIZE_SCHUR
      LOGICAL AGG6
      INTEGER NDENSE(N), LISTVAR_SCHUR(MAX(1,SIZE_SCHUR))
      INTEGER PERM(N)
      INTEGER THRESH
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME, DEG1
      LOGICAL IDENSE
      INTEGER FDEG, ThresMin, ThresPrev, IBEGSchur, NbSchur, 
     $        ThresMinINIT
      LOGICAL SchurON
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER*8 HASH, HMOD 
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
#if (SP_DEBUG >= 1)
      REAL FILL
#endif
      INTRINSIC MAX, MIN, MOD
#if (SP_DEBUG >= 1)
      FILL = 0.0E0
#endif
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0
           NV(1) = 1
           RETURN
        ENDIF
        SIZE_SCHUR = MIN(N,SIZE_SCHUR)
        SIZE_SCHUR = MAX(0,SIZE_SCHUR)
        SchurON   = (SIZE_SCHUR > 0)
        IBEGSchur = N-SIZE_SCHUR+1
        IF (THRESH.GT.N) THRESH = N
        IF (THRESH.LT.0) THRESH = 0
        IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0
                ENDIF
             ENDIF
           ENDDO
        ENDIF
        IF (SchurON) THEN
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
        ELSE
             THRESM    = MAX(INT(31*N/32),THRESH)
             THRESM    = MAX(THRESM,1)
             ThresMin  = MAX( 3*THRESM / 4, 1)
             ThresPrev = THRESM
        ENDIF
        ThresMinINIT = ThresMin/4
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
          THRESM = N
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
            NBD = NBD+1
            IF (FDEG.NE.N+1) THEN
             DEGREE(I) = DEGREE(I)+N+2
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+1
             DEGREE(I) = N+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
   30 IF (NEL .LT. N) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.N)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          MINDEG = 1
          NBDM = MAX(NBDM,NBD)
          IF (DEGREE(ME).GT.N+1) THEN
            IF (WFLG+NBD+1 .LE. WFLG) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + LEN(ME) -1
            LN       = P1
            ELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0) THEN
                X = E
  53            X = -PE(X)
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(LN) = IW(ELN)
               IW(ELN) = E
               LN  = LN+1
               ELN = ELN + 1
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(LN)=E
               LN = LN+1
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = LN-P1
            ELEN(ME) = ELN- P1
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = MIN(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(N+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            IF (THRESM.LT.N) THEN
             ThresMin  = MAX(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = MIN(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = MAX(
     $         THRESM + INT(SQRT(REAL(ThresMin)))+ ThresMinINIT ,
     $         ThresPrev)
             THRESM    = MIN(THRESM,N) 
             ThresMin  = MIN(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.N+1) THEN
           IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in CMUMPS_422 quasi dense rows remains'
            CALL CMUMPS_ABORT()
           ENDIF
           NbSchur = 0   ! Only for checking
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -LISTVAR_SCHUR(1)
            ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
             NbSchur = NbSchur+ 1
            ENDIF
   59      CONTINUE
           IF (NbSchur.NE.SIZE_SCHUR) then
             write(6,*) ' Internal error 2 in QAMD :',
     &         ' Schur size expected:',SIZE_SCHUR, 'Real:', NbSchur
             STOP
           ENDIF
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL CMUMPS_ABORT()
           ENDIF
           IF (ME.NE. LISTVAR_SCHUR(1)) THEN
             DO I=1, SIZE_SCHUR
               PE(LISTVAR_SCHUR(I)) = -LISTVAR_SCHUR(1)
             ENDDO
             PE(LISTVAR_SCHUR(1)) = 0
             NV( LISTVAR_SCHUR(1))= NV(ME)
             NV(ME)               = 0
             ELEN( LISTVAR_SCHUR(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.N) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.N) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.N) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (.NOT. AGG6 .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (AGG6 .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
                PE (E) = -ME
                W (E)  = 0
             ELSE IF (AGG6 .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + E
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.N) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     *     .OR.
     *         (AGG6.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     *       )
     *    THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = MIN (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
#if(SP_DEBUG >= 1)
            IF(HASH .LT. 0) THEN
               WRITE(*,*) 'in symqamd hash =',HASH
               CALL CMUMPS_ABORT()
            ENDIF
#endif
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.N) ) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                IF (PERM(J).GT.PERM(X)) THEN
                  PE (J) = -X
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                  PE (X) = -J
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        IF ( (THRESM .GT. 0).AND.(THRESM.LT.N) ) THEN 
          THRESM = MAX(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.N) THEN
            DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
               IDENSE = .TRUE.
               DEGREE(I) = DEGREE(I)+N+2
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + ELEN(I) - 1
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               FDEG = N
               DEG = N
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
            MINDEG = MIN (MINDEG, FDEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
        PERM(PERMeqN) = N
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE CMUMPS_422 
      SUBROUTINE CMUMPS_423(N, SYM, NPROCS, IORD,
     *                       symmetry,MedDens, NBQD, AvgDens)
      IMPLICIT NONE
      INTEGER N, NPROCS, IORD, SYM
      INTEGER symmetry,MedDens, NBQD, AvgDens
      INTEGER MAXQD
      INTEGER SMALLSYM, SMALLUNS
      PARAMETER (SMALLUNS=5000, SMALLSYM=10000)
      MAXQD = 2
      IF (SYM.NE.0) THEN
         IF ( N.LE.SMALLSYM ) THEN 
            IF (NBQD.GE.MAXQD) THEN
               IORD = 6         ! QAMD
            ELSE
               IORD = 2         ! AMF
            ENDIF
         ELSE
            IF (NBQD.GE.MedDens*NPROCS) THEN
                  IORD = 6      ! QAMD
                  RETURN
            ENDIF
#if  defined(metis)
            IORD = 5
#else
#if defined(pord)
               IORD = 4
#else
               IORD = 6
#endif
#endif
         ENDIF
      ELSE
       IF ( N.LE.SMALLUNS ) THEN
         IF (NBQD.GE.MAXQD) THEN
          IORD = 6  ! QAMD
         ELSE
          IORD = 2   ! AMF
         ENDIF
       ELSE
        IF (NBQD.GE.MedDens*NPROCS) THEN
           IORD = 6      ! QAMD
            RETURN
         ENDIF
#if  defined(metis)
        IORD = 5
#else
#if defined(pord)
         IORD = 4
#else
         IORD = 6
#endif
#endif
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_423
      SUBROUTINE CMUMPS_510
     *     (KEEP10, KEEP2, KEEP48 ,KEEP50, NSLAVES)
      IMPLICIT NONE
      INTEGER NSLAVES, KEEP10, KEEP2, KEEP48, KEEP50
      KEEP10 = MAX(KEEP10*KEEP2,1)
#if defined(t3e) 
      KEEP10 = MIN(1500000, KEEP10)
#elif defined(SP_)
      KEEP10 = MIN(3000000, KEEP10)
#else
      KEEP10 = MIN(2000000, KEEP10)
#endif
#if defined(t3e) 
      IF (NSLAVES .GT. 64) THEN
         KEEP10 = 
     &        MIN((8*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ELSE
         KEEP10 = 
     &        MIN((4*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ENDIF 
#else
      IF (NSLAVES.GT.64) THEN
         KEEP10 = 
     &        MIN((6*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ELSE
         KEEP10 = 
     &        MIN((4*(KEEP2*KEEP2)/NSLAVES)+1, KEEP10)
      ENDIF
#endif
         IF (KEEP50 .EQ. 0 ) THEN
            KEEP10 = MAX(KEEP10,((7*KEEP2*KEEP2)/
     *           4/MAX(NSLAVES-1,1)) + KEEP2)
         ELSE
            KEEP10 = MAX(KEEP10,((7*KEEP2*KEEP2)/
     *           4/(MAX(NSLAVES-1,1))) +KEEP2)
         ENDIF
      IF (KEEP50 .EQ. 0 ) THEN
#if defined(t3e)
         KEEP10 = MAX(KEEP10,200000)
#else 
         KEEP10 = MAX(KEEP10,300000)
#endif
      ELSE
#if defined(t3e)
         KEEP10 = MAX(KEEP10,40000)
#else 
         KEEP10 = MAX(KEEP10,80000)
#endif
      ENDIF
      KEEP10 = -KEEP10 
      RETURN
      END SUBROUTINE CMUMPS_510
      SUBROUTINE CMUMPS_559(JOB,M,N,NE,
     *     IP,IRN,A,LA,NUM,PERM,LIW,IW,LDW,DW,
     *     ICNTL,CNTL,INFO)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL=10, NCNTL=10, NINFO=10)
      INTEGER JOB,M,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),PERM(M),IW(LIW)
      INTEGER ICNTL(NICNTL),INFO(NINFO)
      INTEGER LA
      REAL A(LA)
      REAL DW(LDW),CNTL(NCNTL)
      INTEGER I,J,K,WARN1,WARN2,WARN4
      REAL FACT,ZERO,ONE,RINF,RINF2,RINF3
      PARAMETER (ZERO=0.0E+00,ONE=1.0E+0)
      EXTERNAL CMUMPS_457,CMUMPS_444,CMUMPS_451,CMUMPS_452,CMUMPS_454
      INTRINSIC ABS,LOG
      RINF = CNTL(2)
      RINF2 = HUGE(RINF2)/REAL(2*N)
      RINF3 = 0.0E0
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
      IF (JOB.LT.1 .OR. JOB.GT.6) THEN
         INFO(1) = -1
         INFO(2) = JOB
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
         GO TO 99
      ENDIF
      IF (M.LT.1 .OR. M.LT.N) THEN
         INFO(1) = -2
         INFO(2) = M
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'M',M
         GO TO 99
      ENDIF
      IF (N.LT.1) THEN
         INFO(1) = -2
         INFO(2) = N
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
         GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
         INFO(1) = -3
         INFO(2) = NE
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
         GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = 4*N +   M
      IF (JOB.EQ.2) K = 2*N + 2*M
      IF (JOB.EQ.3) K = 8*N + 2*M + NE
      IF (JOB.EQ.4) K = 3*N + 2*M
      IF (JOB.EQ.5) K = 3*N + 2*M
      IF (JOB.EQ.6) K = 3*N + 2*M + NE
      IF (LIW.LT.K) THEN
         INFO(1) = -4
         INFO(2) = K
         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
         GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
         IF (JOB.EQ.2) K =       M
         IF (JOB.EQ.3) K = 1
         IF (JOB.EQ.4) K =     2*M
         IF (JOB.EQ.5) K = N + 2*M
         IF (JOB.EQ.6) K = N + 3*M
         IF (LDW.LT.K) THEN
            INFO(1) = -5
            INFO(2) = K
            IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
            GO TO 99
         ENDIF
      ENDIF
      IF (ICNTL(5).EQ.0) THEN
         DO 3 I = 1,M
            IW(I) = 0
 3       CONTINUE
         DO 6 J = 1,N
            DO 4 K = IP(J),IP(J+1)-1
               I = IRN(K)
               IF (I.LT.1 .OR. I.GT.M) THEN
                  INFO(1) = -6
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
                  GO TO 99
               ENDIF
               IF (IW(I).EQ.J) THEN
                  INFO(1) = -7
                  INFO(2) = J
                  IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I 
                  GO TO 99
               ELSE
                  IW(I) = J
               ENDIF
 4          CONTINUE
 6       CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9020) JOB,M,N,NE
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,MIN(10,N+1))
               WRITE(ICNTL(3),9022) (IRN(J),J=1,MIN(10,NE))
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,MIN(10,NE))
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
               WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
               IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
            ENDIF
            WRITE(ICNTL(3),9024) (ICNTL(J),J=1,NICNTL)
            WRITE(ICNTL(3),9025) (CNTL(J),J=1,NCNTL)
         ENDIF
      ENDIF
      DO 8 I=1,NINFO
         INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
         DO 10 J = 1,N
            IW(J) = IP(J+1) - IP(J)
 10      CONTINUE
         CALL CMUMPS_457(M,N,IRN,NE,IP,IW(1),PERM,NUM,
     &        IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1))
         GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
         DW(1) = MAX(ZERO,CNTL(1))
         CALL CMUMPS_444(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IW(1),IW(N+1),IW(2*N+1),IW(2*N+M+1),DW,RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
         DO 20 K = 1,NE
            IW(K) = IRN(K)
 20      CONTINUE
         CALL CMUMPS_451(N,NE,IP,IW,A)
         FACT = MAX(ZERO,CNTL(1))
         CALL CMUMPS_452(M,N,NE,IP,IW(1),A,PERM,NUM,IW(NE+1),
     &        IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &        IW(NE+5*N+1),IW(NE+5*N+M+1),FACT,RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
         DO 50 J = 1,N
            FACT = ZERO
            DO 30 K = IP(J),IP(J+1)-1
               IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
 30         CONTINUE
            IF(FACT .GT. RINF3) RINF3 = FACT
            DO 40 K = IP(J),IP(J+1)-1
               A(K) = FACT - ABS(A(K))
 40         CONTINUE
 50      CONTINUE
         DW(1) = MAX(ZERO,CNTL(1))
         DW(2) = RINF3
         IW(1) = JOB
         CALL CMUMPS_454(M,N,NE,IP,IRN,A,PERM,NUM,
     &        IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &        DW(1),DW(M+1),RINF2)
         GO TO 90
      ENDIF
      IF (JOB.EQ.5 .or. JOB.EQ.6) THEN
         RINF3=ONE
         IF (JOB.EQ.5) THEN
            DO 75 J = 1,N
               FACT = ZERO
               DO 60 K = IP(J),IP(J+1)-1
                  IF (A(K).GT.FACT) FACT = A(K)
 60            CONTINUE
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = LOG(FACT)
                  IF(FACT .GT. RINF3) RINF3=FACT
                  DO 70 K = IP(J),IP(J+1)-1
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - LOG(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 70               CONTINUE
               ELSE
                  DO 71 K = IP(J),IP(J+1)-1
                     A(K) = ONE
 71               CONTINUE
               ENDIF
 75         CONTINUE
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 175 K = 1,NE
               IW(3*N+2*M+K) = IRN(K)
 175        CONTINUE
            DO 61 I = 1,M
               DW(2*M+N+I) = ZERO
 61         CONTINUE
            DO 63 J = 1,N
               DO 62 K = IP(J),IP(J+1)-1
                  I = IRN(K)
                  IF (A(K).GT.DW(2*M+N+I)) THEN
                     DW(2*M+N+I) = A(K)
                  ENDIF
 62            CONTINUE
 63         CONTINUE
            DO 64 I = 1,M
               IF (DW(2*M+N+I).NE.ZERO) THEN
                  DW(2*M+N+I) = 1/DW(2*M+N+I)
               ENDIF
 64         CONTINUE
            DO 66 J = 1,N
               DO 65 K = IP(J),IP(J+1)-1
                  I = IRN(K)
                  A(K) = DW(2*M+N+I) * A(K)
 65            CONTINUE
 66         CONTINUE
            CALL CMUMPS_451(N,NE,IP,IW(3*N+2*M+1),A)
            DO 176 J = 1,N
               IF (IP(J).NE.IP(J+1)) THEN
                  FACT = A(IP(J))
               ELSE
                  FACT = ZERO
               ENDIF
               DW(2*M+J) = FACT
               IF (FACT.NE.ZERO) THEN
                  FACT = LOG(FACT)
                  DO 170 K = IP(J),IP(J+1)-1
                     IF (A(K).NE.ZERO) THEN
                        A(K) = FACT - LOG(A(K))
                        IF(A(K) .GT. RINF3) RINF3=A(K)
                     ELSE
                        A(K) = FACT + RINF
                     ENDIF
 170              CONTINUE
               ELSE
                  DO 171 K = IP(J),IP(J+1)-1
                     A(K) = ONE
 171              CONTINUE
               ENDIF
 176        CONTINUE
         ENDIF
         DW(1) = MAX(ZERO,CNTL(1))
         RINF3 = RINF3+ONE
         DW(2) = RINF3
         IW(1) = JOB
         IF (JOB.EQ.5) THEN
            CALL CMUMPS_454(M,N,NE,IP,IRN,A,PERM,NUM,
     &           IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &           DW(1),DW(M+1),RINF2)
         ENDIF
         IF (JOB.EQ.6) THEN
            CALL CMUMPS_454(M,N,NE,IP,IW(3*N+2*M+1),A,PERM,NUM,
     &           IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &           DW(1),DW(M+1),RINF2)
         ENDIF
         IF (JOB.EQ.6) THEN
            DO 79 I = 1,M
               IF (DW(2*M+N+I).NE.0) THEN
                  DW(I) = DW(I) + LOG(DW(2*M+N+I))
               ENDIF
 79         CONTINUE
         ENDIF
         IF (NUM.EQ.N) THEN
            DO 80 J = 1,N
               IF (DW(2*M+J).NE.ZERO) THEN
                  DW(M+J) = DW(M+J) - LOG(DW(2*M+J))
               ELSE
                  DW(M+J) = ZERO
               ENDIF
 80         CONTINUE
         ENDIF
         FACT = 0.5*LOG(RINF2)
         DO 86 I = 1,M
            IF (DW(I).LT.FACT) GO TO 86
            WARN2 = 2
            GO TO 90
 86      CONTINUE 
         DO 87 J = 1,N
            IF (DW(M+J).LT.FACT) GO TO 87
            WARN2 = 2
            GO TO 90
 87      CONTINUE 
      ENDIF
 90   IF (NUM.LT.N) WARN1 = 1
      IF (JOB.EQ.4 .OR. JOB.EQ.5 .OR. JOB.EQ.6) THEN 
         IF (CNTL(1).LT.ZERO) WARN4 = 4
      ENDIF
      IF (INFO(1).EQ.0) THEN
         INFO(1) = WARN1 + WARN2 + WARN4
         IF (INFO(1).GT.0 .AND. ICNTL(2).GE.0) THEN
            WRITE(ICNTL(2),9010) INFO(1)
            IF (WARN1.EQ.1) WRITE(ICNTL(2),9011)
            IF (WARN2.EQ.2) WRITE(ICNTL(2),9012)
            IF (WARN4.EQ.4) WRITE(ICNTL(2),9014)
         ENDIF
      ENDIF
      IF (ICNTL(3).GE.0) THEN
         IF (ICNTL(4).EQ.0 .OR. ICNTL(4).EQ.1) THEN
            WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
            WRITE(ICNTL(3),9031) NUM
            IF (ICNTL(4).EQ.0) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,MIN(10,M))
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,MIN(10,M))
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,MIN(10,N))
               ENDIF
            ELSEIF (ICNTL(4).EQ.1) THEN
               WRITE(ICNTL(3),9032) (PERM(J),J=1,M)
               IF (JOB.EQ.5 .OR. JOB.EQ.6) THEN
                  WRITE(ICNTL(3),9033) (DW(J),J=1,M)
                  WRITE(ICNTL(3),9034) (DW(M+J),J=1,N)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 99   RETURN
 9001 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2,
     &     ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2/
     &     '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2/
     &     '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2/
     &     '        Column ',I8,
     &     ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in CMUMPS_443. INFO(1) = ',I2/
     &     '        CNTL(2) is negative or zero.')
 9010 FORMAT (' ****** Warning from CMUMPS_443. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')
 9020 FORMAT (' ****** Input parameters for CMUMPS_443:'/
     &     ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)   = ',8I8/(15X,8I8))
 9022 FORMAT (' IRN(1:NE)   = ',8I8/(15X,8I8))
 9023 FORMAT (' A(1:NE)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10) = ',8I8/(15X,2I8))
 9025 FORMAT (' CNTL(1:10)  = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for CMUMPS_443:'/
     &     ' INFO(1:2)   = ',2I8)
 9031 FORMAT (' NUM         = ',I8)
 9032 FORMAT (' PERM(1:M)   = ',8I8/(15X,8I8))
 9033 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END SUBROUTINE CMUMPS_559
      SUBROUTINE CMUMPS_563(N,NZ,IP,IRN,A,FLAG,POSI)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER IP(N+1),IRN(NZ)
      REAL A(NZ)
      INTEGER WR_POS,BEG_COL,ROW,COL,K,SV_POS
      INTEGER FLAG(N), POSI(N)
      INTEGER ierr
      FLAG = 0
      WR_POS = 1
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               A(WR_POS) = A(K)
               FLAG(ROW) = COL
               POSI(ROW) = WR_POS
               WR_POS = WR_POS+1
            ELSE
               SV_POS = POSI(ROW)
               A(SV_POS) = A(SV_POS) + A(K)
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1
      RETURN
      END SUBROUTINE CMUMPS_563
      SUBROUTINE CMUMPS_562(N,NZ,IP,IRN,FLAG,POSI)
      IMPLICIT NONE
      INTEGER N,NZ
      INTEGER IP(N+1),IRN(NZ)
      INTEGER WR_POS,BEG_COL,ROW,COL,K,SV_POS
      INTEGER FLAG(N), POSI(N)
      INTEGER ierr
      FLAG = 0
      WR_POS = 1
      DO COL=1,N
         BEG_COL = WR_POS
         DO K=IP(COL),IP(COL+1)-1
            ROW = IRN(K)
            IF(FLAG(ROW) .NE. COL) THEN
               IRN(WR_POS) = ROW
               FLAG(ROW) = COL
               POSI(ROW) = WR_POS
               WR_POS = WR_POS+1
            ENDIF
         ENDDO
         IP(COL) = BEG_COL
      ENDDO
      IP(N+1) = WR_POS
      NZ = WR_POS-1
      RETURN
      END SUBROUTINE CMUMPS_562
