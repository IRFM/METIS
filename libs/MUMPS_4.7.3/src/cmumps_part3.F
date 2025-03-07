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
      RECURSIVE SUBROUTINE CMUMPS_210( COMM_LOAD, ASS_IRECV,
     *  BUFR, LBUFR, LBUFR_BYTES,
     *
     *  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     *  NFRONT_PERE, NASS_PERE, NFS4FATHER,LMAP, TROW,
     *  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *  LRLUS, N, IW,
     *  LIW, A, LA, NIRBDU,
     *  PTRIST, PTLUST_S, PTRFAC,
     *  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *  ICNTL, KEEP,KEEP8,
     *  root, OPASSW, OPELIW,
     *  ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     *  LPTRAR, NELT, FRTPTR, FRTELT, 
     *
     *  ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *  )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC ) :: ROOT
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF, NBFIN, POSFAC
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER COMP
      INTEGER NSTK( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE( * )
      INTEGER LMAP 
      INTEGER TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER INDICE_PERE
      INTEGER PDEST, PDEST_MASTER
      INTEGER NFRONT, SIZFR
      LOGICAL FLAG
      INTEGER SHIFTCB_SON, LDA_SON
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND, IROW_SON,
     *        NPIV, NROWS_TO_STACK, II, COLLIST, POSROW
      INTEGER NBCOLS_EFF
      INTEGER PDEST_MASTER_ISON, IPOS_IN_SLAVE
      INTEGER ROW_LENGTH
      LOGICAL DESCLU, SLAVE_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED 
      INTEGER MSGSOU, MSGTAG
      INTEGER LP
      INTEGER ITMP
      LOGICAL SAME_PROC, COMPRESSCB
      INCLUDE 'mumps_headers.h'
      INTEGER CMUMPS_275, CMUMPS_330
      EXTERNAL CMUMPS_275, CMUMPS_330
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAX_ROW_LENGTH,
     *                                      SIZE_SLAVE_REALS
      INTEGER MAX_ROW_LENGTH_I, SIZE_SLAVE_REALS_I, ITYPE
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) myid,
     *  ' : PB allocation SLAVES_PERE in CMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      IF (NSLAVES_PERE.GT.0) 
     &SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = CMUMPS_275( STEP(INODE_PERE),
     *                 PROCNODE_STEPS, SLAVEF )
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) write(LP,*) myid,
     *  ' : PB allocation NBROW in CMUMPS_210'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) THEN
        write(LP,*) myid, ' : PB allocation LMAP in CMUMPS_210'
        ENDIF
        IFLAG  =-13
        IERROR = LMAP
        GOTO 700
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      PDEST_MASTER_ISON = CMUMPS_275(STEP(ISON),
     *                    PROCNODE_STEPS,SLAVEF)
      SLAVE_ISON = PDEST_MASTER_ISON .NE. MYID
      IF (SLAVE_ISON) THEN
        DO WHILE ( PTRIST(STEP( ISON )) .EQ. 0 )
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL CMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    PDEST_MASTER_ISON, MAITRE_DESC_BANDE,
     *    STATUS,
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
          IF ( IFLAG .LT. 0 ) GOTO 600
        END DO
        DO WHILE (
     *     ( IW( PTRIST(STEP(ISON)) + 1 + XSIZE ) .NE.
     *       IW( PTRIST(STEP(ISON)) + 3 + XSIZE ) ) .OR.
     *     ( KEEP(50) .NE. 0 .AND.
     *       IW( PTRIST(STEP(ISON)) + 6 + XSIZE ) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1 + XSIZE ) .NE.
     *           IW( PTRIST(STEP(ISON)) + 3 + XSIZE ) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL CMUMPS_329( COMM_LOAD,
     *    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MSGTAG,
     *    STATUS, 
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *    IFLAG, IERROR, COMM,
     *    NBPROCFILS,
     *    IPOOL, LPOOL, LEAF,
     *    NBFIN, MYID, SLAVEF,
     *
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) GOTO 600
        END DO
      ENDIF
      IF ( NSLAVES_PERE .EQ. 0 ) THEN
        NBROW( 0 ) = LMAP
      ELSE
        DO I = 0, NSLAVES_PERE
          NBROW( I ) = 0
        END DO
        DO I = 1, LMAP_LOC
          INDICE_PERE = MAP( I )
          CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *           NASS_PERE,
     *           NFRONT_PERE - NASS_PERE,
     *           NSLAVES_PERE,
     *           INDICE_PERE,
     *           NOSLA,
     *           IPOS_IN_SLAVE )
          NBROW( NOSLA ) = NBROW( NOSLA ) + 1
        END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
          IF (LP.GT.0) THEN
          write(LP,*) MYID,': PB allocation PERM in CMUMPS_210'
          ENDIF
          IFLAG  =-13
          IERROR = LMAP_LOC
          GOTO 700
      endif
      IF (keep(50) .NE. 0 ) THEN
         ALLOCATE(SIZE_SLAVE_REALS(0:NSLAVES_PERE),
     *            MAX_ROW_LENGTH(0:NSLAVES_PERE),stat=allocok)
         if (allocok .GT. 0) THEN
             write(*,*) MYID,
     *       ': PB allocation SIZE_SLV_PERE in CMUMPS_210'
             IFLAG  =-13
             IERROR = 2*NSLAVES_PERE + 2
             GOTO 700
         endif
         SIZE_SLAVE_REALS(0:NSLAVES_PERE)=0
         MAX_ROW_LENGTH  (0:NSLAVES_PERE)=0
      ENDIF
         ISTCHK     = PTRIST(STEP(ISON))
         NBCOLS     = IW(ISTCHK+XSIZE)
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *           NASS_PERE,
     *           NFRONT_PERE - NASS_PERE,
     *           NSLAVES_PERE,
     *           INDICE_PERE,
     *           NOSLA,
     *           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
          IF (KEEP(50) .NE. 0) THEN
            ROW_LENGTH = NBCOLS + I - LMAP_LOC
            SIZE_SLAVE_REALS(NOSLA) = SIZE_SLAVE_REALS(NOSLA) +
     *                                ROW_LENGTH
            MAX_ROW_LENGTH(NOSLA) =
     *                    MAX(MAX_ROW_LENGTH(NOSLA),ROW_LENGTH)
          ENDIF
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = SLAVES_PERE(0)
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .EQ. MYID ) THEN 
            NBPROCFILS(STEP(INODE_PERE)) =
     *                 NBPROCFILS(STEP(INODE_PERE)) - 1
            IF ( PDEST .EQ. PDEST_MASTER ) THEN
              NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
            ENDIF
            ISTCHK     = PTRIST(STEP(ISON))
            NBCOLS     = IW(ISTCHK+XSIZE)
            NROW       = IW(ISTCHK+2+XSIZE)
            NPIV       = IW(ISTCHK+3+XSIZE)
            NSLSON     = IW(ISTCHK+5+XSIZE)
            NFRONT     = NPIV + NBCOLS
            COMPRESSCB = (IW(ISTCHK+XXS).EQ.S_CB1COMP)
            SIZFR      = IW(ISTCHK+XXR)
            IF (IW(ISTCHK+XXS).EQ.S_NOLCBCONTIG) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = NPIV*NROW
            ELSE IF (IW(ISTCHK+XXS).EQ.S_NOLCLEANED) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = 0
            ELSE
               LDA_SON     = NFRONT
               SHIFTCB_SON = NPIV
            ENDIF
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
            ELSE
              NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
            ENDIF
            IF (PDEST .NE. PDEST_MASTER) THEN
               IF ( KEEP(55) .eq. 0 ) THEN
                 CALL CMUMPS_539
     *           (N, INODE_PERE, IW, LIW,
     *           A, LA, NROWS_TO_STACK, NBCOLS,
     *           OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     *           KEEP,KEEP8, MYID )
               ELSE
                 CALL CMUMPS_123(NELT, FRTPTR, FRTELT,
     *           N, INODE_PERE, IW, LIW,
     *           A, LA, NROWS_TO_STACK, NBCOLS, 
     *           OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *           FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     *           KEEP, KEEP8, MYID )
               ENDIF
            ENDIF
            DO II = 1,NROWS_TO_STACK
              IROW_SON = PERM(NBROW(I)+II-1)
              INDICE_PERE=MAP(IROW_SON)
             CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
              INDICE_PERE = IPOS_IN_SLAVE
              IF ( COMPRESSCB ) THEN
                IF (NBCOLS - NROW .EQ. 0 ) THEN
                  ITMP = IROW_SON 
                  POSROW = PTRAST(STEP(ISON))+ ITMP * ( ITMP - 1 ) / 2
                ELSE
                  ITMP = IROW_SON + NBCOLS - NROW
                  POSROW = PTRAST(STEP(ISON))+ ITMP * ( ITMP - 1 ) / 2
     *                   - (NBCOLS - NROW) * (NBCOLS - NROW+1)/2
                ENDIF
              ELSE
              POSROW = PTRAST(STEP(ISON)) + SHIFTCB_SON
     *               +(IROW_SON-1)*LDA_SON
              ENDIF
                IF (PDEST == PDEST_MASTER) THEN
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 CALL CMUMPS_39(N, INODE_PERE, IW, LIW, 
     *            A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE,
     *            A(POSROW), PTLUST_S, PTRAST,
     *            STEP, PIMASTER, OPASSW,
     *            IWPOSCB, MYID, KEEP,KEEP8)
                ELSE
                 ISTCHK  = PTRIST(STEP(ISON))
                 COLLIST = ISTCHK + 6 + XSIZE 
     &                     + IW( ISTCHK + 5 +XSIZE) + NROW + NPIV
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 CALL CMUMPS_40(N, INODE_PERE,
     *            IW, LIW,
     *            A, LA, 1, NBCOLS_EFF, INDICE_PERE,
     *            IW( COLLIST ), A(POSROW),
     *            OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *            FILS, ICNTL, KEEP,KEEP8,
     *            MYID )
                ENDIF
            ENDDO
            IF (PDEST.EQ.PDEST_MASTER) THEN 
             IF (KEEP(219).NE.0) THEN
               IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
                  IF (COMPRESSCB) THEN
                    WRITE(*,*) "Error 1 in PARPIV/CMUMPS_210"
                    CALL CMUMPS_ABORT()
                  ELSE
                  POSROW = PTRAST(STEP(ISON))+SHIFTCB_SON+
     *                 (NBROW(1)-1)*LDA_SON
                  ENDIF
                  CALL CMUMPS_617(NFS4FATHER,IERR)
                  IF (IERR .NE.0) THEN
                    IF (LP .GT. 0) THEN
                      WRITE(LP, *) "MAX_ARRAY allocation failed"
                    ENDIF
                    IFLAG=-13
                    IERROR=NFS4FATHER
                    GOTO 700
                  ENDIF
                  ITMP=-9999
                  IF ( LMAP_LOC-NBROW(1)+1 .NE. 0 ) THEN
                  CALL CMUMPS_618(
     *                 A(POSROW),SIZFR-SHIFTCB_SON-(NBROW(1)-1)*LDA_SON,
     *                 LDA_SON, LMAP_LOC-NBROW(1)+1,
     *                 BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,ITMP)
                  ENDIF
                  CALL CMUMPS_619(N, INODE_PERE, IW, LIW, 
     *                 A, LA, ISON, NFS4FATHER,
     *                 BUF_MAX_ARRAY, PTLUST_S, PTRAST,
     *                 STEP, PIMASTER,
     *                 OPASSW,IWPOSCB,MYID, KEEP,KEEP8)
               ENDIF
             ENDIF !PARPIV
             IF ( NBPROCFILS(STEP(ISON)) .EQ. 0) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL CMUMPS_530(N, ISON, INODE_PERE,
     *             IWPOSCB, PIMASTER, PTLUST_S, IW, LIW, STEP,
     *             KEEP,KEEP8)
               ENDIF
               IF (SAME_PROC) THEN
                 ISTCHK_LOC = PTRIST(STEP(ISON))
                 PTRIST(STEP( ISON) ) = -99999999
               ELSE
                 PIMASTER(STEP( ISON )) = -99999999
               ENDIF
               CALL CMUMPS_152(.FALSE., MYID, N,
     *            ISTCHK_LOC,
     *            PAMASTER(STEP(ISON)),
     *            IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU, IWPOSCB,
     *            LA, KEEP,KEEP8
     *            )
             ENDIF
             IF ( NBPROCFILS(STEP(INODE_PERE)) .EQ. 0 ) THEN
               CALL CMUMPS_507( N, IPOOL, LPOOL,
     *           PROCNODE_STEPS,
     *           SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     *           KEEP(47), STEP, INODE_PERE+N )
               IF (KEEP(47) .GE. 3) THEN
                 CALL CMUMPS_500(
     $          IPOOL, LPOOL, 
     *          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *          MYID, STEP, N, ND, FILS )
               ENDIF
             END IF
            ELSE
             CALL CMUMPS_531
     *       (N, INODE_PERE, IW, LIW,
     *       NBROW(I), STEP, PTRIST, ITLOC, KEEP,KEEP8)
            END IF
        END IF
      END DO
      DO I = NSLAVES_PERE, 0, -1
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
 95         DESCLU = .FALSE.
            IF ( PTRIST(STEP(ISON)) .lt.0 .or.
     *         IW ( PTRIST(STEP(ISON) )+XSIZE ) .GT. N ) THEN
              WRITE(*,*) MYID,': Internal error in Maplig'
              WRITE(*,*) MYID,': PTRIST(STEP(ISON))=',
     *                           PTRIST(STEP(ISON))
              WRITE(*,*) MYID,': I, NBROW(I)=',I, NBROW(I)
              WRITE(*,*) MYID,': NSLAVES_PERE=',NSLAVES_PERE
              WRITE(*,*) MYID,': ISON, INODE_PERE=',ISON,INODE_PERE
              WRITE(*,*) MYID,': Son header=',
     *                IW(PTRIST(STEP(ISON)): PTRIST(STEP(ISON))+5+XSIZE)
              CALL CMUMPS_ABORT()
            END IF
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
            ELSE
              NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
            ENDIF
            IF (KEEP(50).NE.0) THEN
              SIZE_SLAVE_REALS_I = SIZE_SLAVE_REALS(I)
              MAX_ROW_LENGTH_I   = MAX_ROW_LENGTH(I)
            ELSE
              SIZE_SLAVE_REALS_I = -99999
              MAX_ROW_LENGTH_I   = -99999
            ENDIF
            COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS).EQ.S_CB1COMP)
            CALL CMUMPS_67( DESCLU, INODE_PERE,
     *      NFRONT_PERE, NASS_PERE, NFS4FATHER,
     *           NSLAVES_PERE, ISON,
     *      NROWS_TO_SEND, LMAP_LOC, MAP,
     *      PERM(MIN(LMAP_LOC,NBROW(I))),
     *      IW( PTRIST(STEP(ISON))),
     *      A(PTRAST(STEP(ISON))), I, PDEST, PDEST_MASTER, 
     *      COMM, IERR, 
     *
     *      KEEP,KEEP8, STEP, N, SLAVEF,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE, MAX_ROW_LENGTH_I,
     *      SIZE_SLAVE_REALS_I, COMPRESSCB )
            IF ( IERR .EQ. -2 ) THEN
              IFLAG  = -17
              IF (LP .GT. 0) THEN
                WRITE(LP,*)
     *          "FAILURE: SEND BUFFER TOO SMALL IN CMUMPS_210"
              ENDIF
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     *        NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+XSIZE) * KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP .GT. 0) THEN
                WRITE(LP,*)
     *          "FAILURE: RECV BUFFER TOO SMALL IN CMUMPS_210"
              ENDIF
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     *        NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+XSIZE) * KEEP( 35 )
              GOTO 700
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
              IERROR = NFS4FATHER
              IF (LP .GT. 0) THEN
                WRITE(LP, *)
     * "FAILURE: MAX_ARRAY allocation failed IN CMUMPS_210"
              ENDIF
              GO TO 700
             END IF
            END IF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL CMUMPS_329( COMM_LOAD,
     *          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *          MPI_ANY_SOURCE, MPI_ANY_TAG,
     *          STATUS,
     *          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *          IWPOS, IWPOSCB, IPTRLU,
     *          LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, 
     *          PTRIST, PTLUST_S, PTRFAC,
     *          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *          IFLAG, IERROR, COMM,
     *          NBPROCFILS,
     *          IPOOL, LPOOL, LEAF,
     *          NBFIN, MYID, SLAVEF,
     *
     *          root, OPASSW, OPELIW, ITLOC, FILS, 
     *          PTRARW, PTRAIW,
     *          INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     *          NELT, FRTPTR, FRTELT, 
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ITYPE = CMUMPS_330(STEP(ISON), PROCNODE_STEPS,SLAVEF)
      IF (KEEP(214) .EQ. 2) THEN
        CALL CMUMPS_314( N, ISON,
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, ITYPE
     *     )
      ENDIF
      CALL CMUMPS_626( N, ISON, PTRIST, PTRAST, IW, LIW,
     $             A, LA, LRLU, LRLUS, IWPOSCB, NIRBDU, IPTRLU,
     $             STEP, MYID, KEEP
     $)
 600  CONTINUE
      DEALLOCATE(SLAVES_PERE)
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      DEALLOCATE(PERM)
      IF (KEEP(50).NE.0) THEN
        DEALLOCATE(SIZE_SLAVE_REALS,MAX_ROW_LENGTH)
      ENDIF
      RETURN
 700  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_210
      SUBROUTINE CMUMPS_211( COMM_LOAD, ASS_IRECV, 
     *  BUFR, LBUFR, LBUFR_BYTES,
     *
     *  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     *  NFRONT_PERE, NASS_PERE, NFS4FATHER, LMAP, TROW,
     *  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *  LRLUS, N, IW,
     *  LIW, A, LA, NIRBDU,
     *  PTRIST, PTLUST_S, PTRFAC,
     *  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *  ICNTL, KEEP,KEEP8, root,
     *  OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *  ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     *
     *  ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *  )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER SLAVEF, NBFIN, POSFAC
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      COMPLEX A( LA )
      INTEGER COMP
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE(NSLAVES_PERE)
      INTEGER NELIM, LMAP, TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER LPTRAR, NELT, ROW_LENGTH
      INTEGER IW( LIW )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK( KEEP(28) ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *PTRAST(KEEP(28)), STEP(N), PIMASTER(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER LP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER APOS
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER INDICE_PERE, LREQI, LREQA
      INTEGER PDEST, PDEST_MASTER, NFRONT
      LOGICAL FLAG, SAME_PROC, FREE, DESCLU
      INTEGER SIZFI, SIZFR, ASIZE
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND,
     *        NPIV, NROWS_TO_STACK, II, HS, IROW_SON,
     *        POSROW, IPOS_IN_SLAVE
      INTEGER NBCOLS_EFF
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL COMPRESSCB
      INCLUDE 'mumps_headers.h'
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM,
     *SIZE_SLAVE_REALS, MAX_ROW_LENGTH
      INTEGER SIZE_SLAVE_REALS_I, MAX_ROW_LENGTH_I
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      if (NSLAVES_PERE.le.0) then
       write(6,*) ' error 2 in maplig_fils_niv1 ', NSLAVES_PERE
       stop
      endif
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      IF (allocok .GT. 0) THEN
        IF (LP > 0)
     *  write(LP,*) myid,
     *  ' : PB allocation NBROW in CMUMPS_211'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat =allocok)
      IF ( allocok .GT. 0 ) THEN
        IF (LP > 0) write(LP,*) myid,
     *  ' : PB allocation SLAVES_PERE in CMUMPS_211'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = CMUMPS_275( STEP(INODE_PERE),
     *      PROCNODE_STEPS, SLAVEF )
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) myid,
     *   ' : PB allocation LMAP in CMUMPS_211'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
      MAP( 1 : LMAP_LOC ) = TROW( 1 : LMAP_LOC )
      DO I = 0, NSLAVES_PERE
        NBROW( I ) = 0
      END DO
      IF (NSLAVES_PERE == 0) THEN
        NBROW(0) = LMAP_LOC
      ELSE
       DO I = 1, LMAP_LOC
        INDICE_PERE = MAP( I )
        CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
        NBROW( NOSLA ) = NBROW( NOSLA ) + 1
       END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     *  ': PB allocation PERM in CMUMPS_211'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
      IF (keep(50) .NE. 0 ) THEN
           ALLOCATE(SIZE_SLAVE_REALS(0:NSLAVES_PERE), 
     *              MAX_ROW_LENGTH(0:NSLAVES_PERE),stat=allocok)
           if (allocok .GT. 0) THEN
             IF (LP > 0) write(LP,*) MYID,
     *       ': PB allocation SIZE_SLV_PERE in CMUMPS_211'
             IFLAG  =-13
             IERROR = 2*NSLAVES_PERE + 2
             GOTO 700
           endif
           SIZE_SLAVE_REALS(0:NSLAVES_PERE)=0
           MAX_ROW_LENGTH(0:NSLAVES_PERE)=0
      ENDIF
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+XSIZE)
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *           NASS_PERE,
     *           NFRONT_PERE - NASS_PERE,
     *           NSLAVES_PERE,
     *           INDICE_PERE,
     *           NOSLA,
     *           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
          IF (KEEP(50) .NE. 0) THEN
            ROW_LENGTH = NBCOLS + I - LMAP_LOC
            SIZE_SLAVE_REALS(NOSLA) = SIZE_SLAVE_REALS(NOSLA) +
     *                                ROW_LENGTH
            MAX_ROW_LENGTH(NOSLA) =
     *                    MAX(MAX_ROW_LENGTH(NOSLA),ROW_LENGTH)
          ENDIF
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = MYID
      IF ( SLAVES_PERE(0) .NE. MYID ) THEN
        WRITE(*,*) 'Error 1 in MAPLIG_FILS_NIV1:',MYID, SLAVES_PERE
        stop
      END IF
      PDEST        = PDEST_MASTER
        I = 0
        NBPROCFILS(STEP(INODE_PERE))=NBPROCFILS(STEP(INODE_PERE))-1
        NBPROCFILS(STEP(ISON))  = NBPROCFILS(STEP(ISON)) - 1
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+XSIZE)
        NELIM      = IW(ISTCHK+1+XSIZE)
        NROW       = IW(ISTCHK+2+XSIZE)
        NPIV       = IW(ISTCHK+3+XSIZE)
        IF (NPIV.LT.0) THEN
         write(6,*) ' Error 2 in CMUMPS_211 ', NPIV
         CALL CMUMPS_ABORT()
        ENDIF
        NSLSON     = IW(ISTCHK+5+XSIZE)
        NFRONT     = NPIV + NBCOLS
        COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS) .eq. S_CB1COMP)
        IF (I == NSLAVES_PERE) THEN
          NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
        ELSE
          NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
        ENDIF
        DO II = 1,NROWS_TO_STACK
          IROW_SON=PERM(NBROW(I)+II-1)
          INDICE_PERE = MAP(IROW_SON)
          CALL CMUMPS_47(
     *         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *         NASS_PERE,
     *         NFRONT_PERE - NASS_PERE,
     *         NSLAVES_PERE,
     *         INDICE_PERE,
     *         NOSLA,
     *         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF (COMPRESSCB) THEN
            IF (NELIM.EQ.0) THEN
            POSROW = PAMASTER(STEP(ISON)) +
     *         IROW_SON*(IROW_SON-1)/2
            ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     *         (NELIM+IROW_SON)*(NELIM+IROW_SON-1)/2
            ENDIF
          ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     *             (NELIM+IROW_SON-1)*NBCOLS
          ENDIF
          IF (KEEP(50).NE.0) THEN
            NBCOLS_EFF = NELIM + IROW_SON
          ELSE
            NBCOLS_EFF = NBCOLS
          ENDIF
          CALL CMUMPS_39(N, INODE_PERE, IW, LIW, 
     *    A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE,
     *    A(POSROW), PTLUST_S, PTRAST,
     *    STEP, PIMASTER, OPASSW, IWPOSCB, 
     *    MYID, KEEP,KEEP8)
        ENDDO
        IF (KEEP(219).NE.0) THEN
         IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
           IF (COMPRESSCB) THEN
              POSROW = PAMASTER(STEP(ISON))+
     *                 (NELIM+NBROW(1))*(NELIM+NBROW(1)-1)/2
              ASIZE  = (LMAP_LOC+NELIM)*(NELIM+LMAP_LOC+1)/2
     *               - (NELIM+NBROW(1)) * (NELIM+NBROW(1)-1)/2
           ELSE
              POSROW = PAMASTER(STEP(ISON)) +
     *                 (NELIM+NBROW(1)-1)*NBCOLS
              ASIZE  = (LMAP_LOC-NBROW(1)+1) * NBCOLS
           ENDIF
           CALL CMUMPS_617(NFS4FATHER,IERR)
           IF (IERR .NE.0) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     *    ": PB allocation MAX_ARRAY during CMUMPS_211"
              IFLAG=-13
              IERROR=NFS4FATHER
              GOTO 700
           ENDIF
           IF  ( LMAP_LOC-NBROW(1)+1 .NE. 0 ) THEN
           CALL CMUMPS_618(
     *          A(POSROW),ASIZE,NBCOLS,LMAP_LOC-NBROW(1)+1,
     *          BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,
     *          NELIM+NBROW(1))
           ENDIF
           CALL CMUMPS_619(N, INODE_PERE, IW, LIW, 
     *          A, LA, ISON, NFS4FATHER,
     *          BUF_MAX_ARRAY, PTLUST_S, PTRAST,
     *          STEP, PIMASTER, OPASSW,
     *          IWPOSCB,MYID, KEEP,KEEP8)
         ENDIF
        ENDIF !PARPIV
          IF ( NBPROCFILS(STEP(ISON)) .EQ. 0) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL CMUMPS_530(N, ISON, INODE_PERE,
     *            IWPOSCB, PIMASTER, PTLUST_S, IW, LIW, STEP,
     *            KEEP,KEEP8)
               ENDIF
          ENDIF
          IF ( NBPROCFILS(STEP(INODE_PERE)) .EQ. 0 ) THEN
            CALL CMUMPS_507( N, IPOOL, LPOOL,
     *        PROCNODE_STEPS,
     *        SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     *        KEEP(47), STEP, INODE_PERE+N )
            IF (KEEP(47) .GE. 3) THEN
              CALL CMUMPS_500(
     $       IPOOL, LPOOL, 
     *       PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *       MYID, STEP, N, ND, FILS )
            ENDIF
          END IF
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
 95        NFRONT = IW(PIMASTER(STEP(ISON))+XSIZE)
           NELIM  = IW(PIMASTER(STEP(ISON))+1+XSIZE)
           APOS = PAMASTER(STEP(ISON))
           DESCLU = .TRUE.
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
            IF (KEEP(50).NE.0) THEN
              SIZE_SLAVE_REALS_I = SIZE_SLAVE_REALS(I)
              MAX_ROW_LENGTH_I   = MAX_ROW_LENGTH(I)
            ELSE
              SIZE_SLAVE_REALS_I = -99999
              MAX_ROW_LENGTH_I   = -99999
            ENDIF
           CALL CMUMPS_67( DESCLU, INODE_PERE,
     *      NFRONT_PERE, NASS_PERE, NFS4FATHER, 
     *           NSLAVES_PERE,
     *      ISON, NROWS_TO_SEND, LMAP_LOC,
     *      MAP, PERM(MIN(LMAP_LOC,NBROW(I))),
     *      IW(PIMASTER(STEP(ISON))),
     *      A(APOS), I, PDEST, PDEST_MASTER, COMM, IERR,
     *
     *      KEEP,KEEP8, STEP, N, SLAVEF,
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *      MAX_ROW_LENGTH_I, SIZE_SLAVE_REALS_I, COMPRESSCB)
            IF ( IERR .EQ. -2 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     *": FAILURE, SEND BUFFER TOO SMALL DURING CMUMPS_211"
              IFLAG  = -17
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     *        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     *": FAILURE, RECV BUFFER TOO SMALL DURING CMUMPS_211"
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     *        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
               IERROR = BUF_LMAX_ARRAY
              IF (LP > 0) WRITE(LP,*) MYID,
     *": FAILURE, MAX_ARRAY ALLOC FAILED DURING CMUMPS_211"
               GO TO 700
             ENDIF
            ENDIF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .FALSE.
              MESSAGE_RECEIVED = .FALSE.
              CALL CMUMPS_329( COMM_LOAD,
     *          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *          MPI_ANY_SOURCE, MPI_ANY_TAG,
     *          STATUS, 
     *          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *          IWPOS, IWPOSCB, IPTRLU,
     *          LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *          PTLUST_S, PTRFAC,
     *          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     *          IFLAG, IERROR, COMM,
     *          NBPROCFILS,
     *          IPOOL, LPOOL, LEAF,
     *          NBFIN, MYID, SLAVEF,
     *          root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *          INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *          LPTRAR, NELT, FRTPTR, FRTELT, 
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP( ISON )) = -77777777
            IF ( IW(ISTCHK+XSIZE) .GE. 0 ) THEN
              WRITE(*,*) 'error 3 in CMUMPS_211'
              CALL CMUMPS_ABORT()
            ENDIF
      CALL CMUMPS_152(.FALSE., MYID, N, ISTCHK,
     *     PAMASTER(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP,KEEP8
     *     )
 600  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      DEALLOCATE(PERM)
      DEALLOCATE(SLAVES_PERE)
      IF (KEEP(50).NE.0) THEN
        DEALLOCATE(SIZE_SLAVE_REALS,MAX_ROW_LENGTH)
      ENDIF
      RETURN
 700  CONTINUE
      CALL CMUMPS_44(MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_211
      SUBROUTINE CMUMPS_209( N, FRERE, FILS, NFSIZ, THEROOT )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: N
      INTEGER, intent( in    )  :: NFSIZ( N )
      INTEGER, intent( inout )  :: FRERE( N ), FILS( N )
      INTEGER, intent( out   )  :: THEROOT
      INTEGER INODE, IROOT, IFILS, IN, IROOTLAST, SIZE
      IROOT = -9999
      SIZE  = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .EQ. 0 )  THEN
          IF ( NFSIZ( INODE ) .GT. SIZE ) THEN
            SIZE  = NFSIZ( INODE )
            IROOT = INODE
          END IF
        ENDIF
      END DO
      IN = IROOT
      DO WHILE ( FILS( IN ) .GT. 0 )
        IN = FILS( IN )
      END DO
      IROOTLAST = IN
      IFILS     = - FILS ( IN )
      DO INODE = 1, N
        IF ( FRERE( INODE ) .eq. 0 .and. INODE .ne. IROOT ) THEN
          IF ( IFILS .eq. 0 ) THEN
            FILS( IROOTLAST ) = - INODE
            FRERE( INODE )    = -IROOT
            IFILS             = INODE
          ELSE
            FRERE( INODE ) = -FILS( IROOTLAST )
            FILS( IROOTLAST ) = - INODE
          END IF
        END IF
      END DO
      THEROOT = IROOT
      RETURN
      END SUBROUTINE CMUMPS_209
      SUBROUTINE CMUMPS_93(MYID,N,IOLDPS,TYPE,IW, LIW, A, LA,
     $POSFAC, LRLU, LRLUS, IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8,
     $SSARBR,INODE,IERR)
      USE CMUMPS_LOAD
      USE CMUMPS_OOC
      IMPLICIT NONE
      INTEGER MYID
      INTEGER IOLDPS, TYPE, LIW, LA, N, KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER POSFAC, LRLU, LRLUS, IWPOS, LDLT
      INTEGER PTRAST(KEEP(28)), PTRFAC(KEEP(28)), STEP( N )
      LOGICAL SSARBR
      INTEGER IOLDSHIFT, IPSSHIFT
      INCLUDE 'mumps_headers.h'
      INTEGER LCONT, NELIM, NROW, NPIV, IAPOS, INTSIZ
      INTEGER NFRONT, NSLAVES
      INTEGER IPS, IPSIZE, SIZELU, SIZECB, I
      LOGICAL MOVEPTRAST
      INTEGER INODE
      INTEGER IERR
      IERR=0
      LDLT = KEEP(50)
      IOLDSHIFT = IOLDPS + XSIZE
      IF ( IW( IOLDSHIFT ) < 0 ) THEN
        write(*,*) ' ERROR 1 compressLU:Should not point to a band.'
        stop
      ELSE IF ( IW( IOLDSHIFT + 2 ) < 0 ) THEN
        write(*,*) ' ERROR 2 compressLU:Stack not performed yet',
     *  IW(IOLDSHIFT + 2)
        stop
      ENDIF
      LCONT  = IW( IOLDSHIFT )
      NELIM  = IW( IOLDSHIFT + 1 )
      NROW   = IW( IOLDSHIFT + 2 )
      NPIV   = IW( IOLDSHIFT + 3 )
      IAPOS  = PTRFAC(IW( IOLDSHIFT + 4 ))
      NSLAVES= IW( IOLDSHIFT + 5 )
      NFRONT = LCONT + NPIV
      INTSIZ = IW(IOLDPS+XXI)
      IF ( (NSLAVES > 0  .AND. TYPE .NE. 2) .OR. 
     *   (NSLAVES .eq. 0 .AND. TYPE .EQ. 2 ) ) THEN
          WRITE(*,*) ' ERROR 3 compressLU: problem with level of inode'
          stop
      END IF
      IF (LDLT.EQ.0) THEN
        SIZELU = ( LCONT + NROW ) * NPIV
      ELSE
        SIZELU =   NROW * NPIV
      ENDIF
      IF ( TYPE .EQ. 2 ) THEN
        IF (LDLT.EQ.0) THEN
         SIZECB = NELIM * LCONT
        ELSE
         IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
            SIZECB = (NELIM+1) * ( NELIM + NPIV )
          ELSE
            SIZECB = NELIM * ( NELIM + NPIV )
          ENDIF
        ENDIF
      ELSE
        IF (LDLT.EQ.0) THEN
         SIZECB = LCONT * LCONT
        ELSE
         SIZECB = NROW * LCONT
        ENDIF
      END IF
      IW(IOLDPS+XXR)=IW(IOLDPS+XXR)-SIZECB
      IF ((SIZECB.EQ.0).AND.(KEEP(201).EQ.0)) THEN
         GOTO 500
      ENDIF
      IF ( IOLDPS + INTSIZ .NE. IWPOS ) THEN
         IF (KEEP(201).NE.0) THEN
            KEEP(31)=KEEP(31)+SIZELU
            CALL CMUMPS_576(INODE,PTRFAC,KEEP,KEEP8,
     *           A,LA,SIZELU, IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
         ENDIF
         IPS = IOLDPS + INTSIZ
         MOVEPTRAST = .FALSE.
         DO WHILE ( IPS .NE. IWPOS )
           IPSIZE = IW(IPS+XXI)    
           IPSSHIFT = IPS + XSIZE
           IF ( IW( IPSSHIFT + 2 ) < 0 ) THEN
             NFRONT = IW( IPSSHIFT )
#if ! defined(SIMULATE_OOC_FAC)
             IF(KEEP(201).EQ.0)THEN
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                      PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
#endif
                PTRFAC(IW(IPSSHIFT+4))=PTRFAC(IW(IPSSHIFT+4)) -
     $               SIZECB - SIZELU
#if ! defined(SIMULATE_OOC_FAC)
             ENDIF
#endif
             MOVEPTRAST = .TRUE.
#if ! defined(SIMULATE_OOC_FAC)
             IF(KEEP(201).EQ.0)THEN
              PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
             ELSE
#endif
              PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
     $               - SIZELU
#if ! defined(SIMULATE_OOC_FAC)
             ENDIF
#endif
           ELSE IF ( IW( IPSSHIFT ) < 0 ) THEN
#if ! defined(SIMULATE_OOC_FAC)
             IF(KEEP(201).EQ.0)THEN
                PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))-SIZECB
             ELSE
#endif
                PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))
     &                                   -SIZECB-SIZELU
#if ! defined(SIMULATE_OOC_FAC)
             ENDIF
#endif
           ELSE
             NFRONT = IW( IPSSHIFT ) + IW( IPSSHIFT + 3 )
#if ! defined(SIMULATE_OOC_FAC)
             IF(KEEP(201).EQ.0)THEN
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                    PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
#endif
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &               PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
     $               - SIZELU
#if ! defined(SIMULATE_OOC_FAC)
             ENDIF
#endif
           END IF
           IPS = IPS + IPSIZE
         END DO
         IF ((SIZECB .NE. 0).OR.(KEEP(201).NE.0)) THEN
#if ! defined (SIMULATE_OOC_FAC)
            IF (KEEP(201).NE.0) THEN
#endif
               DO I=IAPOS, POSFAC - SIZECB - SIZELU - 1
                  A( I ) = A( I + SIZECB + SIZELU)
               END DO
#if ! defined (SIMULATE_OOC_FAC)
            ELSE
               DO I=IAPOS + SIZELU, POSFAC - SIZECB - 1
                  A( I ) = A( I + SIZECB )
               END DO
            ENDIF
#endif     
         END IF
      ELSE
         IF (KEEP(201).NE.0) THEN
            KEEP(31)=KEEP(31)+SIZELU
            CALL CMUMPS_576(INODE,PTRFAC,KEEP,KEEP8,A,LA,SIZELU,
     $           IERR)
            IF(IERR.LT.0)THEN
               RETURN
            ENDIF
         ENDIF
      ENDIF
      POSFAC = POSFAC - SIZECB
      LRLU   = LRLU   + SIZECB
      LRLUS  = LRLUS  + SIZECB
#if ! defined (SIMULATE_OOC_FAC)
      IF (KEEP(201).NE.0) THEN
#endif
       POSFAC = POSFAC + SIZECB - (SIZECB+SIZELU)
       LRLU   = LRLU   - SIZECB + (SIZECB+SIZELU)
       LRLUS  = LRLUS  - SIZECB + (SIZECB+SIZELU)
#if ! defined (SIMULATE_OOC_FAC)
      ENDIF
#endif      
 500  CONTINUE
      CALL CMUMPS_471(SSARBR,.FALSE.,
     *         LA-LRLUS,SIZELU,-SIZECB,KEEP,KEEP8,LRLU)
      RETURN
      END SUBROUTINE CMUMPS_93
      SUBROUTINE CMUMPS_314( N, ISON, 
     *    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA, 
     *    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     *    NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, TYPE_SON
     *     )
      USE CMUMPS_OOC
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N, ISON, LIW, LA, LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC,
     *        COMP, NIRBDU, IPTRLU, IFLAG, IERROR, SLAVEF, MYID, COMM,
     *        TYPE_SON
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)), PAMASTER(KEEP(28)), IW(LIW)
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)), ITLOC(N)
      DOUBLE PRECISION OPELIW
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      COMPLEX A( LA )
      INTEGER  NFRONT, NCOL_L, NROW_L, LREQI, LREQA,NSLAVES_L,
     *         POSA, POSI, POSALOC, ISTCHK, SIZFR, SIZFI, I,
     *         NSLSON, HS, IROW_L, ICOL_L, OLDPOS, NROW, NCOL,
     *         LDA_BAND, NASS
      LOGICAL SAME_PROC, FREE
      INCLUDE 'mumps_headers.h'
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      FLOP1 = ZERO
      NCOL_L = IW( PTRIST(STEP( ISON )) + 3 + XSIZE )
      NROW_L = IW( PTRIST(STEP( ISON )) + 2 + XSIZE )
      NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 + XSIZE )
      LDA_BAND = NCOL_L + IW( PTRIST(STEP( ISON )) + XSIZE )
      IF  ( KEEP(50) .eq. 0 ) THEN
        NFRONT = LDA_BAND
      ELSE
        NFRONT = IW( PTRIST(STEP( ISON )) + 7 + XSIZE )
      END IF
      IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
         GOTO 80
      ENDIF
      LREQI   = 4 + NCOL_L + NROW_L + XSIZE
      LREQA   = NCOL_L * NROW_L
      IF ( LRLU .LT. LREQA .OR.
     *  IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG  = -9
          IERROR = lREQA - LRLUS
          GO TO 700
        END IF
        CALL CMUMPS_94( N,KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS,IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
        COMP = COMP+1
        IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB compress Stack_band:LRLU,LRLUS=',
     *            LRLU, LRLUS
                  IFLAG = -9
                  IERROR = LREQA - LRLU
                  GOTO 700
        END IF
        IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
          IFLAG  = -8
          IERROR = IWPOS + LREQI - 1 - IWPOSCB
          GOTO 700
        END IF
      END IF
      POSA = POSFAC
      POSFAC = POSFAC + LREQA
      LRLU = LRLU - LREQA
      LRLUS = LRLUS - LREQA
      KEEP(67) = MIN(LRLUS, KEEP(67))
#if ! defined (SIMULATE_OOC_FAC)
      IF (KEEP(201).NE.0) THEN
#endif
      POSA = POSFAC
      POSFAC = POSFAC - LREQA
      LRLU = LRLU + LREQA
      LRLUS = LRLUS + LREQA
#if ! defined (SIMULATE_OOC_FAC)
      ENDIF
#endif      
      CALL CMUMPS_471(.FALSE.,.FALSE.,
     *             LA-LRLUS,LREQA,LREQA,KEEP,KEEP8,LRLU)
#if ! defined (SIMULATE_OOC_FAC)
      IF (KEEP(201).NE.0) THEN
#endif
       POSA = POSFAC
       POSFAC = POSFAC + LREQA
       LRLU = LRLU - LREQA
       LRLUS = LRLUS - LREQA
#if ! defined (SIMULATE_OOC_FAC)
      ENDIF
#endif      
      POSI = IWPOS
      IWPOS = IWPOS + LREQI
      PTLUST_S(STEP( ISON )) = POSI
      IW(POSI+XXI)=LREQI
      IW(POSI+XXR)=LREQA
      IW(POSI+XXS)=-9999
      POSI=POSI+XSIZE
      IW( POSI     ) = - NCOL_L
      IW( POSI + 1 ) =   NROW_L
      IW( POSI + 2 ) =   NFRONT - NCOL_L
      IW( POSI + 3 ) =   STEP(ISON)
      PTRFAC(STEP(ISON)) = POSA
      IROW_L = PTRIST(STEP( ISON )) + 6 + NSLAVES_L + XSIZE
      ICOL_L = PTRIST(STEP( ISON )) + 6 + NROW_L + NSLAVES_L + XSIZE
      IW( POSI + 4 : POSI + 3 + NROW_L ) =
     *                       IW( IROW_L: IROW_L + NROW_L - 1 )
      IW( POSI + NROW_L + 4: POSI + NROW_L + NCOL_L + 3 ) =
     *                       IW( ICOL_L: ICOL_L + NCOL_L - 1 )
      POSALOC = POSA
      DO I = 1, NROW_L
        OLDPOS =  PTRAST( STEP(ISON)) + ( I - 1 ) * LDA_BAND
        A( POSALOC : POSALOC + NCOL_L - 1 ) = 
     *                   A( OLDPOS: OLDPOS + NCOL_L - 1 )
        POSALOC = POSALOC + NCOL_L
      END DO
      IF (KEEP(201).NE.0) THEN
       KEEP(31)=KEEP(31)+LREQA
       CALL CMUMPS_576(ISON,PTRFAC,KEEP,KEEP8,A,LA,LREQA,IFLAG)
       IF(IFLAG.LT.0)THEN
         WRITE(*,*)MYID,': Internal error in CMUMPS_576'
         IERROR=0
         GOTO 700
       ENDIF
      ENDIF
#if ! defined (SIMULATE_OOC_FAC)
      IF (KEEP(201).NE.0) THEN
#endif
       POSA = POSFAC
       POSFAC = POSFAC - LREQA
       LRLU = LRLU + LREQA
       LRLUS = LRLUS + LREQA
#if ! defined (SIMULATE_OOC_FAC)
      ENDIF
#endif      
  80  CONTINUE
      IF (TYPE_SON == 1) THEN
         GOTO 90
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = DBLE( NCOL_L * NROW_L) +
     *     DBLE(NROW_L*NCOL_L)*DBLE(2*NFRONT-NCOL_L-1)
      ELSE
         FLOP1 = DBLE( NCOL_L ) * DBLE( NROW_L )
     *            * ( 2 * LDA_BAND - NROW_L - NCOL_L + 1)
      END IF
      OPELIW = OPELIW + FLOP1
      FLOP1_EFFECTIVE = FLOP1
      NASS = IW( PTRIST(STEP( ISON )) + 4 + XSIZE )
      IF ( NCOL_L .NE. NASS ) THEN
        IF ( KEEP(50).eq.0 ) THEN
           FLOP1 = DBLE( NASS * NROW_L) +
     *     DBLE(NROW_L*NASS)*DBLE(2*NFRONT-NASS-1)
        ELSE
           FLOP1 = DBLE( NASS ) * DBLE( NROW_L ) *
     *     DBLE( 2 * LDA_BAND - NROW_L - NASS + 1)
        END IF
      END IF
      CALL CMUMPS_190(1,.FALSE.,FLOP1_EFFECTIVE-FLOP1,
     *                        KEEP,KEEP8)
      CALL CMUMPS_190(2,.FALSE.,-FLOP1,KEEP,KEEP8)
 90   CONTINUE
      RETURN
 700  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_314
      SUBROUTINE CMUMPS_626( N, ISON, 
     *    PTRIST, PTRAST, IW, LIW, A, LA, 
     *    LRLU, LRLUS, IWPOSCB,
     *    NIRBDU, IPTRLU, STEP, MYID, KEEP
     $     )
      IMPLICIT NONE
      include 'mumps_headers.h'
      INTEGER ISON, MYID, N, LRLU, LRLUS, IWPOSCB, NIRBDU, IPTRLU
      INTEGER KEEP(500), STEP(N)
      INTEGER*8 KEEP8(150)
      INTEGER PTRAST(KEEP(28)), PTRIST(KEEP(28))
      INTEGER LA, LIW
      INTEGER IW(LIW)
      COMPLEX A(LA)
      INTEGER ISTCHK
      ISTCHK = PTRIST(STEP(ISON))
      CALL CMUMPS_152(.FALSE.,MYID, N, ISTCHK,
     *     PTRAST(STEP(ISON)),
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP,KEEP8
     *     )
      PTRIST(STEP( ISON )) = -9999888
      PTRAST(STEP( ISON )) = -9999888
      RETURN
      END SUBROUTINE CMUMPS_626
      INTEGER FUNCTION CMUMPS_330(INODE, PROCNODE, SLAVEF)
      IMPLICIT NONE
      INTEGER INODE, SLAVEF 
      INTEGER PROCNODE(*)
      INTEGER MYPROC, TPN
      MYPROC = PROCNODE( INODE )
      IF (MYPROC <= SLAVEF ) THEN
        CMUMPS_330 = 1
      ELSE
        TPN = (MYPROC-1+2*SLAVEF)/SLAVEF - 1
        IF ( TPN .LT. 1 ) TPN = 1
        CMUMPS_330 = TPN
      END IF
      RETURN 
      END
      INTEGER FUNCTION CMUMPS_275(INODE, PROCNODE, SLAVEF)
      IMPLICIT NONE
      INTEGER INODE, SLAVEF 
      INTEGER PROCNODE(*)
      INTEGER MYPROC
      IF (SLAVEF == 1) THEN
        CMUMPS_275 = 0
      ELSE
        MYPROC = PROCNODE( INODE )
        CMUMPS_275=MOD(2*SLAVEF+MYPROC-1,SLAVEF)
      END IF
      RETURN
      END
      LOGICAL FUNCTION CMUMPS_283( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+2*SLAVEF)/SLAVEF - 1
      CMUMPS_283 = ( TPN .eq. 0 )
      RETURN
      END
      LOGICAL FUNCTION CMUMPS_167( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+SLAVEF+SLAVEF)/SLAVEF - 1
      CMUMPS_167 = ( TPN .eq. -1 )
      RETURN 
      END
      LOGICAL FUNCTION CMUMPS_170
     *        ( INODE, PROCNODE, SLAVEF )
      IMPLICIT NONE
      INTEGER INODE, SLAVEF
      INTEGER  PROCNODE(*)
      INTEGER TPN, MYPROC
      MYPROC = PROCNODE( INODE )
      TPN = (MYPROC-1+SLAVEF+SLAVEF)/SLAVEF - 1
      CMUMPS_170 =
     *           ( TPN .eq. -1 .OR. TPN .eq. 0 )
      RETURN
      END
      LOGICAL FUNCTION CMUMPS_358( MYID, SLAVEF, INODE,
     *                 NMB_PAR2, ISTEP_TO_INIV2 , K71, STEP, N, 
     *                 CANDIDATES, KEEP24 )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, INODE, NMB_PAR2, KEEP24, I
      INTEGER K71, N
      INTEGER ISTEP_TO_INIV2 ( K71 ), STEP ( N )
      INTEGER CANDIDATES(SLAVEF+1, MAX(NMB_PAR2,1))
      INTEGER NCAND, POSINODE
      CMUMPS_358 = .FALSE.
      IF (KEEP24 .eq. 0) RETURN
      POSINODE = ISTEP_TO_INIV2 ( STEP (INODE) )
      NCAND = CANDIDATES( SLAVEF+1, POSINODE )
      DO I = 1, NCAND
        IF (MYID .EQ. CANDIDATES( I, POSINODE ))
     *     CMUMPS_358 = .TRUE.
      END DO
      RETURN
      END
      SUBROUTINE CMUMPS_291(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()
      RETURN
      END SUBROUTINE CMUMPS_291
      SUBROUTINE CMUMPS_292(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()-T
      RETURN
      END SUBROUTINE CMUMPS_292
      SUBROUTINE CMUMPS_312( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      REAL VAL( N )
      INTEGER I, ISWAP
      REAL SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_312
      SUBROUTINE CMUMPS_558( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_558
      SUBROUTINE CMUMPS_158
     *( N, MYID, STEP, PROCNODE_STEPS, K28, SLAVEF, NA, LNA,  MYLEAF )
      IMPLICIT NONE
      INTEGER, intent(in) :: N, SLAVEF, MYID, K28, LNA
      INTEGER, intent(in) :: PROCNODE_STEPS( K28 ), STEP(N), NA( LNA )
      INTEGER, intent(out) :: MYLEAF
      INTEGER INODE, NBLEAF, I
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      MYLEAF = 0
      NBLEAF = NA(1)
      DO I = 1, NBLEAF
        INODE = NA(2+I)
        IF ( CMUMPS_275(STEP(INODE),PROCNODE_STEPS,SLAVEF)
     *       .eq. MYID ) THEN
          MYLEAF = MYLEAF + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_158
#if defined (PESSL)
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     $                     LLD, INFO )
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
#if defined(DESC8)
      PARAMETER          ( DLEN_ = 8, DTYPE_ = 1,
     $                     CTXT_ = 7, M_ = 1, N_ = 2, MB_ = 3, NB_ = 4,
     $                     RSRC_ = 5, CSRC_ = 6, LLD_ = 8 )
#else
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
#endif
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           BLACS_GRIDINFO, PXERBLA
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      INTRINSIC          MAX, MIN
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.MAX( 1, NUMROC( M, MB, MYROW, IRSRC,
     $                                NPROW ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 )
     $   CALL PXERBLA( ICTXT, 'DESCINIT', -INFO )
#ifndef DESC8
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
#endif
      DESC( M_ )  = MAX( 0, M )
      DESC( N_ )  = MAX( 0, N )
      DESC( MB_ ) = MAX( 1, MB )
      DESC( NB_ ) = MAX( 1, NB )
      DESC( RSRC_ ) = MAX( 0, MIN( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = MAX( 0, MIN( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = MAX( LLD, MAX( 1, NUMROC( DESC( M_ ), DESC( MB_ ),
     $                              MYROW, DESC( RSRC_ ), NPROW ) ) )
      RETURN
      END SUBROUTINE DESCINIT
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
      INTEGER            ICTXT, INFO
      CHARACTER*(*)      SRNAME
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           BLACS_GRIDINFO
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number', I4, ' had an illegal value' )
      END SUBROUTINE PXERBLA
#endif
      SUBROUTINE CMUMPS_256( N, NZ, IRN, ICN, ASPK, X, Y,
     *                         LDLT, MTYPE, MAXTRANS, PERM )
      INTEGER N, NZ, LDLT, MTYPE, MAXTRANS
      INTEGER IRN( NZ ), ICN( NZ ) 
      INTEGER PERM( N )
      COMPLEX ASPK( NZ ), X( N ), Y( N )
      INTEGER K, I, J
      COMPLEX PX( N )
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      Y = CMPLX(ZERO)
      IF ( MAXTRANS .eq. 1 .and. MTYPE .eq. 1) THEN
        DO I = 1, N
          PX(I) = X( PERM( I ) )
        END DO
      ELSE
        PX = X
      END IF
      IF ( LDLT .eq. 0 ) THEN
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(J) = Y(J) + ASPK(K) * PX(I)
        ENDDO
       ENDIF
      ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     *        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
          IF (J.NE.I) THEN
            Y(J) = Y(J) + ASPK(K) * PX(I)
          ENDIF
        ENDDO
      END IF
      IF ( MAXTRANS .EQ. 1 .AND. MTYPE .eq. 0 ) THEN
      PX = Y
      DO I = 1, N
        Y( PERM( I ) ) = PX( I )
      END DO
      END IF
      RETURN
      END SUBROUTINE CMUMPS_256
      SUBROUTINE CMUMPS_214( KEEP,KEEP8,
     *           MYID, N, NELT, LNA, NZ, NA_ELT, NSLAVES,
     *           MEMORY_MBYTES, EFF, OOC, PERLU_ON, MEMORY )
      IMPLICIT NONE
      LOGICAL,   INTENT(IN)  :: EFF, OOC, PERLU_ON
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER MYID, N, NELT, NSLAVES, LNA, NZ, NA_ELT
      INTEGER*8, INTENT(OUT) :: MEMORY
      INTEGER,   INTENT(OUT) :: MEMORY_MBYTES
      LOGICAL I_AM_SLAVE, I_AM_MASTER
      INTEGER NB_INT, NB_REAL, PERLU, NBRECORDS
      REAL DMEMORY_BYTES, TEMP, NB_BYTES
      INTEGER  CMUMPS_LBUF_INT, CMUMPS_LBUFR_BYTES, CMUMPS_LBUF
      INTEGER TEMPI, TEMPR, NBUFS
      INTEGER MAXS_MIN, MIN_PERLU
      PERLU     = KEEP(12)
      IF (.NOT.PERLU_ON) PERLU = 0
      I_AM_MASTER = ( MYID .eq. 0 )
      I_AM_SLAVE  = ( KEEP(46).eq. 1 .or. MYID .ne. 0 )
      TEMP    = 0
      NB_REAL = 0
      NB_BYTES = 0.0
      NB_INT  = 0
      NB_INT = NB_INT + 5*KEEP(28)
      NB_INT = NB_INT +KEEP(28)+KEEP(56)*(NSLAVES+2)
      NB_INT = NB_INT + 3*N
      IF (KEEP(23).ne.0 .and. I_AM_MASTER) NB_INT=NB_INT + N
      IF (KEEP(55).eq.0) THEN
        NB_INT = NB_INT + 2 * N
      ELSE
        NB_INT = NB_INT + 2 * ( NELT+1 )
      ENDIF
      IF (KEEP(55) .ne. 0 ) THEN
        NB_INT = NB_INT + N+1 + NELT
      END IF
      NB_INT = NB_INT + LNA
      IF ( OOC ) THEN
        MAXS_MIN = KEEP(200)
        NB_REAL = NB_REAL + KEEP(100)
        NB_INT=NB_INT+4*KEEP(28)+KEEP(28)
      ELSE
        MAXS_MIN = KEEP(16)
      ENDIF
      IF ( .NOT. EFF ) THEN
        IF (KEEP(50).EQ.1) THEN
         NB_REAL = NB_REAL + MAXS_MIN +PERLU*(MAXS_MIN / 100 + 1 )
        ELSE
         NB_REAL = NB_REAL + MAXS_MIN +PERLU*(MAXS_MIN / 100 + 1 )
        ENDIF
      ELSE
        NB_REAL = NB_REAL + KEEP(67)
      ENDIF
      NB_REAL = NB_REAL + KEEP(13)
      IF ( .not. ( I_AM_SLAVE .and. I_AM_MASTER .and. KEEP(52) .eq. 0
     *         .and. KEEP(55) .ne. 0 ) ) THEN
        NB_INT  = NB_INT  + KEEP(14)
      END IF
      IF ( I_AM_SLAVE .and. KEEP(38) .ne. 0 ) THEN
        NB_INT = NB_INT + 2 * N
      END IF
      TEMPI= 0
      TEMPR = 0
      NBRECORDS = KEEP(39)
      IF (KEEP(55).eq.0) THEN
        NBRECORDS = MIN(KEEP(39), NZ)
      ELSE
        NBRECORDS = MIN(KEEP(39), NA_ELT)
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( I_AM_MASTER ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            NBUFS = NSLAVES 
          ELSE
            NBUFS = NSLAVES - 1
            IF (KEEP(55) .eq. 0 )
     *      TEMPI = TEMPI + 2 * N
          END IF
          TEMPI = TEMPI + 2 * NBRECORDS * NBUFS
          TEMPR = TEMPR + NBRECORDS * NBUFS
        ELSE
          IF ( KEEP(55) .eq. 0 )THEN
            TEMPI = TEMPI + 2 * NBRECORDS
            TEMPR = TEMPR + NBRECORDS
          END IF
        END IF
      ELSE
        IF ( I_AM_SLAVE ) THEN
          TEMPI = TEMPI + ( 1 + 4 * NSLAVES ) * NBRECORDS
          TEMPR = TEMPR + ( 1 + 2 * NSLAVES ) * NBRECORDS
        END IF
      END IF
      TEMP = MAX( NB_BYTES + REAL(NB_INT +TEMPI) * KEEP(34)
     *           + REAL(NB_REAL+TEMPR) * KEEP(35) 
     *            , TEMP )
      IF ( I_AM_SLAVE ) THEN
        CMUMPS_LBUFR_BYTES = KEEP( 44 ) * KEEP( 35 )
        CMUMPS_LBUFR_BYTES = MAX( CMUMPS_LBUFR_BYTES,
     *                      100000 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        CMUMPS_LBUFR_BYTES = CMUMPS_LBUFR_BYTES
     *        + INT( REAL(MAX(PERLU,MIN_PERLU))*
     *        REAL(CMUMPS_LBUFR_BYTES)/100E0)
        IF (KEEP(50).NE.1) THEN
           CMUMPS_LBUFR_BYTES = CMUMPS_LBUFR_BYTES
     *            + INT( REAL(MAX(PERLU,0))*
     *            REAL(CMUMPS_LBUFR_BYTES)/100E0)
        ENDIF
        IF (KEEP(48).EQ.5) THEN
          CMUMPS_LBUFR_BYTES=CMUMPS_LBUFR_BYTES+5*KEEP(34)
        ENDIF
        CMUMPS_LBUFR_BYTES = CMUMPS_LBUFR_BYTES +
     *        (6 + NSLAVES - 2 + 2 * KEEP(127))
     *        * KEEP(34)
        IF(KEEP(219).NE.0.AND.
     *  KEEP(50) .EQ. 2)
     *  CMUMPS_LBUFR_BYTES=CMUMPS_LBUFR_BYTES+
     *              KEEP(108)*KEEP(35)+KEEP(34)
        NB_BYTES = NB_BYTES + CMUMPS_LBUFR_BYTES
        IF (NSLAVES.EQ.2) THEN
         CMUMPS_LBUF = INT( 1.01E0 * KEEP( 43 ) *KEEP( 35 )  )
        ELSE
         CMUMPS_LBUF = INT( REAL(KEEP(213)) / 100.0E0
     *                     * KEEP( 43 ) * KEEP( 35 )  )
        ENDIF
        CMUMPS_LBUF = MAX( CMUMPS_LBUF, 100000 )
        CMUMPS_LBUF = CMUMPS_LBUF
     *                 + INT( 2.0E0 * REAL(MAX(PERLU,0))*
     *                   REAL(CMUMPS_LBUF)/100E0)
        NB_BYTES = NB_BYTES + CMUMPS_LBUF
        CMUMPS_LBUF_INT = 
     *         NSLAVES * NSLAVES * 4
     *               * KEEP(34)
        NB_BYTES = NB_BYTES + CMUMPS_LBUF_INT
        IF ( EFF ) THEN
          NB_INT = NB_INT + KEEP(15)
        ELSE
          NB_INT = NB_INT + 
     *           KEEP(15) + 2 * MAX(PERLU,10) *
     *           ( KEEP(15) / 100 + 1 )
        ENDIF
        NB_INT = NB_INT + KEEP(28)
        NB_INT = NB_INT + KEEP(28)
        NB_INT = NB_INT + N + 7 * KEEP(28) + 1
      END IF
      DMEMORY_BYTES = NB_BYTES + REAL(NB_INT) * KEEP(34) +
     *               REAL(NB_REAL) * KEEP(35)
      DMEMORY_BYTES = MAX( DMEMORY_BYTES, TEMP )
      MEMORY = DMEMORY_BYTES
      MEMORY_MBYTES = INT( DMEMORY_BYTES / 1 000 000 )  + 1
      RETURN
      END SUBROUTINE CMUMPS_214
      SUBROUTINE CMUMPS_243(MYID, COMM, INFO, INFOG, IRANK)
      IMPLICIT NONE
      INTEGER MYID, COMM, IRANK, INFO(1), INFOG(2)
      INCLUDE 'mpif.h'
      INTEGER IERR_MPI, MASTER
      INTEGER TEMP1(2), TEMP2(2)
      PARAMETER( MASTER = 0 )
      CALL MPI_REDUCE( INFO(1), INFOG(1), 1, MPI_INTEGER,
     *                 MPI_MAX, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO(1), INFOG(2), 1, MPI_INTEGER,
     *                 MPI_SUM, MASTER, COMM, IERR_MPI )
      TEMP1(1) = INFO(1)
      TEMP1(2) = MYID
      CALL MPI_REDUCE( TEMP1, TEMP2, 1, MPI_2INTEGER,
     *                 MPI_MAXLOC, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq. MASTER ) THEN
        IF ( INFOG(1) .ne. TEMP2(1) ) THEN
          write(*,*) 'Error in CMUMPS_243'
          stop
        END IF
        IRANK    = TEMP2(2)
      ELSE
        IRANK    = -1
      END IF
      RETURN
      END SUBROUTINE CMUMPS_243
      SUBROUTINE CMUMPS_362(N, LEAF, NBROOT, NROOT_LOC,
     *           MYID_NODES,
     *           SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     *           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, NROOT_LOC, NBROOT, MYID_NODES,
     *        SLAVEF, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), NA(LNA),
     *        IPOOL(LPOOL)
      INTEGER NBLEAF, INODE, I
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      NBLEAF = NA(1)
      NBROOT = NA(2)
      LEAF = 1
      DO I = 1, NBLEAF
        INODE = NA(I+2)
        IF (CMUMPS_275(STEP(INODE),PROCNODE_STEPS,SLAVEF)
     *   .EQ.MYID_NODES) THEN
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
          ENDIF
      ENDDO
      NROOT_LOC = 0
      DO I = 1, NBROOT
        INODE = NA(I+2+NBLEAF)
        IF (CMUMPS_275(STEP(INODE),PROCNODE_STEPS,
     *    SLAVEF).EQ.MYID_NODES) THEN
            NROOT_LOC = NROOT_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_362
      LOGICAL FUNCTION CMUMPS_438(TAB1,TAB2,LEN1,LEN2)
      IMPLICIT NONE
      INTEGER LEN1 , LEN2 ,I
      INTEGER TAB1(LEN1)
      INTEGER TAB2(LEN2)
      CMUMPS_438=.FALSE.
      IF(LEN1 .NE. LEN2) THEN
         RETURN
      ENDIF
      DO I=1 , LEN1
         IF(TAB1(I) .NE. TAB2(I)) THEN
            RETURN
         ENDIF
      ENDDO
      CMUMPS_438=.TRUE.
      RETURN
      END FUNCTION CMUMPS_438
      SUBROUTINE CMUMPS_463( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_463
      SUBROUTINE CMUMPS_466( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE CMUMPS_466
      SUBROUTINE CMUMPS_ABORT()
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
      IERRCODE = -99
      CALL MPI_ABORT(MPI_COMM_WORLD, IERRCODE, IERR)
      RETURN
      END SUBROUTINE CMUMPS_ABORT
      SUBROUTINE CMUMPS_618(
     *     A,ASIZE,NCOL,NROW,
     *     M_ARRAY,NMAX,COMPRESSCB,LROW1)
      IMPLICIT NONE
      INTEGER NROW,NCOL,NMAX,ASIZE,LROW1
      LOGICAL COMPRESSCB
      COMPLEX A(ASIZE)
      REAL M_ARRAY(NMAX)
      INTEGER APOS,I,J, LROW
      REAL ZERO,TMP
      PARAMETER (ZERO=0.0E0)
      M_ARRAY(1:NMAX) = ZERO
      APOS = 0
      IF (COMPRESSCB) THEN
        LROW=LROW1
      ELSE
        LROW=NCOL
      ENDIF
      DO I=1,NROW
         DO J=1,NMAX
            TMP = ABS(A(APOS+J))
            IF(TMP.GT.M_ARRAY(J)) M_ARRAY(J) = TMP
         ENDDO
         APOS = APOS + LROW
         IF (COMPRESSCB) LROW=LROW+1
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_618
      SUBROUTINE CMUMPS_621( ILOC4, ISUM4, ISUM8, COMM)
      IMPLICIT NONE
      INTEGER, intent(in)    :: COMM
      INTEGER, intent(in)    :: ILOC4
      INTEGER, intent(out)   :: ISUM4
      INTEGER*8, intent(out) :: ISUM8
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER*8 IDIFF8
      DOUBLE PRECISION DLOC
      DOUBLE PRECISION DSUM
      CALL MPI_ALLREDUCE(ILOC4,ISUM4,1,MPI_INTEGER,
     *                   MPI_SUM, COMM, IERR)
      DLOC = DBLE(ILOC4)
      CALL MPI_ALLREDUCE(DLOC,DSUM,1,MPI_DOUBLE_PRECISION,
     *                   MPI_SUM, COMM, IERR)
      ISUM8  = INT(DSUM+.5E0,8)
      IDIFF8 = ISUM8 - ISUM4
      IF ( IDIFF8 .NE. 0 ) THEN
        ISUM4 = -( ISUM8 / 1000000 )
      ENDIF
      RETURN
      END SUBROUTINE CMUMPS_621
      SUBROUTINE  CMUMPS_633(KEEP12,ICNTL14,
     *     KEEP50,KEEP54,ICNTL6,ICNTL8)
      IMPLICIT NONE
      INTEGER, intent(out)::KEEP12
      INTEGER, intent(in)::ICNTL14,KEEP50,KEEP54,ICNTL6,ICNTL8
      KEEP12 = ICNTL14 
      IF(ICNTL6.EQ.0 .AND. ICNTL8.EQ.0) RETURN
      IF ( (KEEP54.NE.0).AND. (KEEP50.NE.1)
     *     .AND. (KEEP12 .GT. 0) ) KEEP12= KEEP12+5
      RETURN
      END SUBROUTINE CMUMPS_633
      SUBROUTINE CMUMPS_22( SSARBR, PROCESS_BANDE,
     *   MYID,N, KEEP,KEEP8,
     *   IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP,PIMASTER,PAMASTER, ITLOC,
     *   LREQ, LREQCB, NODE_ARG, STATE_ARG, SET_HEADER,
     *   COMP, LRLUS, IFLAG, IERROR )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU, KEEP(500),
     &        IPTRLU,IWPOS,IWPOSCB
      INTEGER NODE_ARG, STATE_ARG
      INTEGER*8 KEEP8(150)
      INTEGER IW(LIW),PTRIST(KEEP(28)),PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     * PAMASTER(KEEP(28)), ITLOC(N)
      INTEGER MYID, IXXP
      COMPLEX A(LA)
      LOGICAL PROCESS_BANDE, SSARBR, SET_HEADER
      INTEGER COMP, LREQ, LREQCB, LRLUS, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER INODE_LOC,NPIV,NASS,NROW,NCB,MEM_GAIN
      INTEGER ISIZEHOLE, RSIZEHOLE
      LOGICAL DONE
      IF (IWPOSCB.EQ.NIRBDU) THEN
        IF (LREQ.NE.XSIZE.OR.LREQCB.NE.0
     *      .OR. .NOT. SET_HEADER) THEN
          WRITE(*,*) "Internal error in CMUMPS_22",
     *      SET_HEADER, LREQ, LREQCB
          CALL CMUMPS_ABORT()
        ENDIF
        IF (IWPOSCB-IWPOS+1 .LT. XSIZE) THEN
          WRITE(*,*) "Internal error 2 in CMUMPS_22"
          CALL CMUMPS_ABORT()
        ENDIF
        IWPOSCB=IWPOSCB-XSIZE
        IW(IWPOSCB+1+XXI)=XSIZE
        IW(IWPOSCB+1+XXR)=0
        IW(IWPOSCB+1+XXN)=-919191
        IW(IWPOSCB+1+XXS)=S_NOTFREE
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
        RETURN
      ENDIF
      IF (KEEP(214).EQ.1.AND.
     *    KEEP(216).EQ.1.AND.
     *    IWPOSCB.NE.NIRBDU) THEN
       IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG.OR.
     *     IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
        NCB  = IW( IWPOSCB+1 + XSIZE )
        NROW = IW( IWPOSCB+1 + XSIZE + 2)
        NPIV = IW( IWPOSCB+1 + XSIZE + 3)
        INODE_LOC= IW( IWPOSCB+1 + XXN)
        CALL CMUMPS_632(IWPOSCB+1,IW,NIRBDU,
     *                          ISIZEHOLE,RSIZEHOLE)
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG) THEN
          CALL CMUMPS_627(A,LA,IPTRLU+1,
     *                           NROW,NCB,NPIV+NCB,0,
     *                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED
          MEM_GAIN            = NROW*NPIV
        ENDIF
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
          NASS = IW( IWPOSCB+1 + XSIZE + 4)
          CALL CMUMPS_627(A,LA,IPTRLU+1,
     *                           NROW,NCB,NPIV+NCB,NASS-NPIV,
     *                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED38
          MEM_GAIN = NROW*(NPIV+NCB-(NASS-NPIV))
        ENDIF
        IF (ISIZEHOLE.NE.0) THEN
          CALL CMUMPS_630( IW,NIRBDU,IWPOSCB+1,
     *                       IWPOSCB+IW(IWPOSCB+1+XXI),
     *                       ISIZEHOLE )
          IWPOSCB=IWPOSCB+ISIZEHOLE
          IW(IWPOSCB+1+XXP+IW(IWPOSCB+1+XXI))=IWPOSCB+1
          PTRIST(STEP(INODE_LOC))=PTRIST(STEP(INODE_LOC))+
     *    ISIZEHOLE
        ENDIF
        IW(IWPOSCB+1 + XXR) = IW(IWPOSCB+1 + XXR)-MEM_GAIN
        IPTRLU              = IPTRLU+MEM_GAIN+RSIZEHOLE
        LRLU                = LRLU+MEM_GAIN+RSIZEHOLE
        PTRAST(STEP(INODE_LOC))=
     *  PTRAST(STEP(INODE_LOC))+MEM_GAIN+RSIZEHOLE
       ENDIF
      ENDIF
      DONE =.FALSE.
      IF ((IPTRLU.LT.LREQCB).OR.(LRLU.LT.LREQCB)) THEN
        IF (LRLUS.LT.LREQCB) THEN
          GOTO 620
        ELSE
          CALL CMUMPS_94(N,KEEP(28),IW,LIW,A,LA,NIRBDU,
     *                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     *                    PTRIST,PTRAST,
     *                    STEP, PIMASTER,PAMASTER,ITLOC,KEEP(216),LRLUS)
          IF ( LRLU .NE. LRLUS ) THEN
            WRITE(*,*) 'PB compress... alloc_cb',
     *      'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 620
          END IF
          DONE = .TRUE.
          COMP = COMP + 1
        ENDIF
      ENDIF
      IF (IWPOSCB-IWPOS+1 .LT. LREQ) THEN
       IF (DONE) GOTO 600
                 CALL CMUMPS_94(N,KEEP(28),IW,LIW,A,LA,NIRBDU,
     *                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     *                    PTRIST,PTRAST,
     *                    STEP, PIMASTER,PAMASTER,ITLOC,KEEP(216),LRLUS)
                 IF ( LRLU .NE. LRLUS ) THEN
                   WRITE(*,*) 'PB compress... alloc_cb',
     *             'LRLU,LRLUS=',LRLU,LRLUS
                   GOTO 620
                 END IF
          COMP = COMP + 1
          IF (IWPOSCB-IWPOS+1 .LT. LREQ) GOTO 600
      ENDIF
      IXXP=IWPOSCB+XXP+1
      IF (IXXP.GT.NIRBDU) THEN
        WRITE(*,*) "Internal error 3 in CMUMPS_22",IXXP
      ENDIF
      IF (IW(IXXP).GT.0) THEN
        WRITE(*,*) "Internal error 2 in CMUMPS_22",IW(IXXP),IXXP
      ENDIF
      IWPOSCB = IWPOSCB - LREQ
      IF (SET_HEADER) THEN
        IW(IXXP)= IWPOSCB + 1
        IW(IWPOSCB+1+XXI)=LREQ
        IW(IWPOSCB+1+XXR)=LREQCB
        IW(IWPOSCB+1+XXS)=STATE_ARG
        IW(IWPOSCB+1+XXN)=NODE_ARG
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
      ENDIF
      IPTRLU = IPTRLU - LREQCB
      LRLU   = LRLU - LREQCB
      LRLUS  = LRLUS - LREQCB
      KEEP(67) = MIN(LRLUS, KEEP(67))
#if ! defined(OLD_LOAD_MECHANISM)
      CALL CMUMPS_471(SSARBR,PROCESS_BANDE,
     *              LA-LRLUS,0,LREQCB,KEEP,KEEP8,LRLU)
#else
#if defined (CHECK_COHERENCE)
      CALL CMUMPS_471(SSARBR,PROCESS_BANDE,
     *              LA-LRLUS,0,LREQCB,KEEP,KEEP8,LRLU)
#else
      CALL CMUMPS_471(SSARBR,.FALSE.,
     *              LA-LRLUS,0,LREQCB,KEEP,KEEP8,LRLU)
#endif
#endif
      RETURN
 600  IFLAG  = -8
      IERROR = LREQ
      RETURN
 620  IFLAG  = -9
      IERROR = LREQCB - LRLUS
      RETURN
      END SUBROUTINE CMUMPS_22
      SUBROUTINE CMUMPS_244(N, NZ, NSTEPS,
     * A, LA, IW, LIW, SYM_PERM, NA, LNA,
     * NE_STEPS, NFSIZ, FILS,
     * STEP, FRERE, DAD, CAND, 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     * PTRAR, LDPTRAR,
     * PTRIST, PTLUST_S, PTRFAC, IW1, IW2, POOL, LPOOL,  
     * CNTL, ICNTL, INFO, RINFO, KEEP,KEEP8,PROCNODE_STEPS,
     * SLAVEF,
     * COMM_NODES, MYID, MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,
     * root, NELT, FRTPTR, FRTELT, COMM_LOAD,
     * ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     * MEM_DISTRIB,
     * DKEEP,PIVNUL_LIST,LPN_LIST)
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER N,NZ,NSTEPS,LA,LIW,LPOOL,SLAVEF,COMM_NODES
      INTEGER MYID, MYID_NODES,LNA
      COMPLEX A(LA)
      REAL RINFO(20)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER BUFR( LBUFR )
      INTEGER NELT, LDPTRAR
      INTEGER FRTPTR(*), FRTELT(*)
      REAL CNTL(15)
      INTEGER   ICNTL(40)
      INTEGER   INFO(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER   IW(LIW), SYM_PERM(N), NA(LNA),
     *          NE_STEPS(KEEP(28)), FILS(N),
     *          FRERE(KEEP(28)), NFSIZ(KEEP(28)), 
     *          DAD(KEEP(28))
      INTEGER   CAND(SLAVEF+1, MAX(1,KEEP(56)))
      INTEGER   STEP(N)
      INTEGER   PTRAR(LDPTRAR,2)
      INTEGER   PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER   IW1(5*KEEP(28)), IW2(N), POOL(LPOOL)
      INTEGER   PROCNODE_STEPS(KEEP(28))
      INTEGER   COMM_LOAD, ASS_IRECV
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INTEGER   INTARR(MAX(1,KEEP(14)))
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      REAL SEUIL, SEUIL_LDLT_NIV2
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      REAL UULOC 
      INTEGER LP, MP
      INTEGER NSTK,PTRAST, NBPROCFILS
      INTEGER PIMASTER, PAMASTER
      INTEGER NIRBDU
      LOGICAL PROK
      REAL ZERO, ONE
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
      KEEP(41)=0
      KEEP(42)=0
      KEEP(31)= 0
      INFO(11)= 0
      INFO(12)= 0
      INFO(13)= 0
      INFO(14)  = 0
      RINFO(2) = ZERO
      RINFO(3) = ZERO
      NSTEPS   = 0
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK   = (MP.GT.0)
      UULOC = CNTL(1)
      IF (UULOC.GT.ONE)   UULOC=ONE
      IF (UULOC.LT.ZERO)  UULOC=ZERO
      IF (KEEP(50).NE.0.AND.UULOC.GT.0.5E0) THEN
        UULOC = 0.5E0
      ENDIF
      PTRAST     = 1
      PIMASTER   = PTRAST + KEEP(28)
      PAMASTER   = PIMASTER + KEEP(28) 
      NSTK       = PAMASTER + KEEP(28)
      NBPROCFILS = NSTK + KEEP(28)
      NIRBDU = LIW
      CALL CMUMPS_250(N, NZ, NSTEPS, A, LA, IW,
     * LIW, SYM_PERM, NA, LNA, NE_STEPS, IW1(NSTK),
     * IW1( NBPROCFILS ),
     * INFO(11),NFSIZ,FILS,STEP, FRERE, DAD, CAND, 
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * PTRAR(1,2),PTRAR(1,1),IW1(PTRAST),PTRIST,
     * PTLUST_S, PTRFAC,
     * IW1(PIMASTER), IW1(PAMASTER),
     * IW2,NIRBDU,
     * POOL(1), LPOOL ,
     * UULOC,ICNTL, INFO, RINFO, KEEP,KEEP8,PROCNODE_STEPS,
     * SLAVEF,
     * COMM_NODES,MYID,MYID_NODES,
     * BUFR,LBUFR,LBUFR_BYTES,
     * INTARR,DBLARR,root,NELT, FRTPTR, 
     * FRTELT, LDPTRAR, COMM_LOAD, ASS_IRECV, SEUIL,
     * SEUIL_LDLT_NIV2, MEM_DISTRIB,
     *     DKEEP(1),PIVNUL_LIST(1),LPN_LIST
     *     )
      IF (KEEP(201) .NE. 0) THEN
        INFO(9)=KEEP(31)
      ENDIF
      IF (PROK) THEN
        WRITE (MP,99980) INFO(1), INFO(2),
     *       KEEP(28), INFO(9), INFO(10), INFO(11), INFO(12),
     *       INFO(13), INFO(14), RINFO(2), RINFO(3)
      ENDIF
      RETURN
99980 FORMAT (/' LEAVING FACTORIZATION PHASE WITH ...'/
     1      ' INFO (1)                                      =',I12/
     2      '  --- (2)                                      =',I12/
     5      '           NUMBER OF NODES IN THE TREE         =',I12/
     6      ' INFO (9)  REAL SPACE FOR FACTORS              =',I12/
     7      '  --- (10) INTEGER SPACE FOR FACTORS           =',I12/
     8      '  --- (11) MAXIMUM SIZE OF FRONTAL MATRICES    =',I12/
     9      '  --- (12) NUMBER OF OFF DIAGONAL PIVOTS       =',I12/
     1      '  --- (13) NUMBER OF DELAYED PIVOTS            =',I12/
     2      '  --- (14) NUMBER OF MEMORY COMPRESSES         =',I12/
     3  ' RINFO(2)  OPERATIONS DURING NODE ASSEMBLY     =',1PD10.3/
     4  ' -----(3)  OPERATIONS DURING NODE ELIMINATION  =',1PD10.3)
      END SUBROUTINE CMUMPS_244
      SUBROUTINE CMUMPS_269( MYID,KEEP,KEEP8,
     *           BUFR, LBUFR, LBUFR_BYTES,
     *           IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     *           N, IW, LIW, A, LA, NIRBDU,
     *           PTRIST, PTRAST, STEP, PIMASTER, PAMASTER,
     *           NSTK_S, COMP,
     *           FPERE, FLAG, IFLAG, IERROR, COMM, ITLOC )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER KEEP(500), BUFR( LBUFR )
      INTEGER*8 KEEP8(150)
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER PTRIST( KEEP(28) ), PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP, FPERE
      LOGICAL FLAG
      INTEGER NSTK_S( KEEP(28) ), ITLOC( N )
      INTEGER IFLAG, IERROR, COMM
      INTEGER POSITION, FINODE, FLCONT, LREQ, LREQCB
      INCLUDE 'mumps_headers.h'
      LOGICAL COMPRESSCB
      FLAG = .FALSE.
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FINODE, 1, MPI_INTEGER,
     *                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FPERE, 1, MPI_INTEGER, 
     *                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *                FLCONT, 1, MPI_INTEGER,
     *                COMM, IERR)
      COMPRESSCB = (FLCONT.LT.0) 
      IF (COMPRESSCB) THEN
        FLCONT   = -FLCONT
        LREQCB  = (FLCONT * (FLCONT+1)) /2
      ELSE
        LREQCB  = FLCONT * FLCONT
      ENDIF
      LREQ    = 2 * FLCONT + 6 + XSIZE
      IF (IPTRLU.LT.0) WRITE(*,*) 'before alloc_cb:IPTRLU = ',IPTRLU
      CALL CMUMPS_22( .FALSE.,.FALSE.,
     *   MYID,N, KEEP,KEEP8, IW, LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER,PAMASTER, ITLOC,
     *   LREQ, LREQCB, FINODE, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR )
      IF (IPTRLU.LT.0) WRITE(*,*) 'after alloc_cb:IPTRLU = ',IPTRLU
      IF ( IFLAG .LT. 0 ) RETURN
      PIMASTER(STEP( FINODE )) = IWPOSCB + 1
      PAMASTER(STEP( FINODE )) = IPTRLU  + 1
      IF (COMPRESSCB)  IW(IWPOSCB + 1 + XXS ) = S_CB1COMP
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     *        IW(IWPOSCB + 1+XSIZE), LREQ-XSIZE,
     *        MPI_INTEGER, COMM, IERR)
      IF ( LREQCB .ne. 0 ) THEN
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION, A(IPTRLU + 1),
     *        LREQCB, MPI_COMPLEX, COMM, IERR)
      END IF
      NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
      IF ( NSTK_S(STEP(FPERE)).EQ.0 ) THEN
          FLAG = . TRUE.
      END IF
      RETURN
      END SUBROUTINE CMUMPS_269
      SUBROUTINE CMUMPS_270( TOT_ROOT_SIZE,
     *    TOT_CONT_TO_RECV, ROOT,
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
      USE CMUMPS_LOAD
      USE CMUMPS_OOC        
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'cmumps_root.h'
      TYPE (CMUMPS_ROOT_STRUC) :: ROOT
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      INTEGER TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      COMPLEX A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) ), ND( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC(N), FILS(N), PTRARW(N), PTRAIW(N)
      INTEGER INTARR(MAX(1,KEEP(14)))
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER NEW_LOCAL_M, NEW_LOCAL_N
      INTEGER OLD_LOCAL_M, OLD_LOCAL_N
      INTEGER LREQI, LREQA, IROOT
      INTEGER POSHEAD, IPOS_SON,IERR
      LOGICAL MASTER_OF_ROOT
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mumps_headers.h'
      INTEGER NUMROC, CMUMPS_275
      EXTERNAL NUMROC, CMUMPS_275
      IROOT = KEEP( 38 )
      root%TOT_ROOT_SIZE = TOT_ROOT_SIZE
      MASTER_OF_ROOT = ( MYID .EQ. CMUMPS_275( STEP(IROOT),
     *                   PROCNODE_STEPS, SLAVEF ) )
      NEW_LOCAL_M  = NUMROC( TOT_ROOT_SIZE, root%MBLOCK,
     *               root%MYROW, 0, root%NPROW )
      NEW_LOCAL_M  = MAX( 1, NEW_LOCAL_M )
      NEW_LOCAL_N  = NUMROC( TOT_ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
      IF (KEEP(60) .NE. 0) THEN
        IF (root%yes) THEN
        IF ( NEW_LOCAL_M .NE. root%SCHUR_MLOC .OR.
     *       NEW_LOCAL_N .NE. root%SCHUR_NLOC ) THEN
          WRITE(*,*) "Internal error 1 in CMUMPS_270"
          CALL CMUMPS_ABORT()
        ENDIF
        ENDIF
        PTLUST_S(STEP(IROOT)) = -4444
        PTRFAC(STEP(IROOT)) = -4445
        PTRIST(STEP(IROOT)) = 0
        IF ( MASTER_OF_ROOT ) THEN
          LREQI=6+2*TOT_ROOT_SIZE+XSIZE
          LREQA=0
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
           CALL CMUMPS_94( N, KEEP(28), IW, LIW, A, LA,
     *           NIRBDU, LRLU, IPTRLU,
     *           IWPOS, IWPOSCB, PTRIST, PTRAST,
     *           STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS )
           COMP = COMP + 1
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB1 compress root2slave:LRLU,LRLUS=',
     *            LRLU, LRLUS
                  IFLAG = -9
                  IERROR = LREQA - LRLUS
                  GOTO 700
           END IF
          ENDIF
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
            IFLAG = -8
            IERROR = IWPOS + LREQI - 1 - IWPOSCB
            GOTO 700
          ENDIF
          PTLUST_S(STEP(IROOT))= IWPOS
          IWPOS = IWPOS + LREQI
          POSHEAD = PTLUST_S( STEP(IROOT))
          IW( POSHEAD + XXI )=LREQI
          IW( POSHEAD + XXR )=LREQA
          IW( POSHEAD + XXS )=-9999
          IW( POSHEAD +XSIZE) = 0
          IW( POSHEAD + 1 +XSIZE) = -1
          IW( POSHEAD + 2 +XSIZE) = -1
          IW( POSHEAD + 4 +XSIZE) = STEP(IROOT)
          IW( POSHEAD + 5 +XSIZE) = 0
          IW( POSHEAD + 3 +XSIZE) = TOT_ROOT_SIZE
        ENDIF
        GOTO 100
      ENDIF
      IF ( MASTER_OF_ROOT ) THEN
        LREQI = 6 + 2 * TOT_ROOT_SIZE+XSIZE
      ELSE
        LREQI = 6+XSIZE
      END IF
      LREQA = NEW_LOCAL_M * NEW_LOCAL_N
      IF ( LRLU . LT. LREQA .OR.
     *     IWPOS + LREQI - 1. GT. IWPOSCB )THEN
           IF ( LRLUS .LT. LREQA ) THEN
             IFLAG  = -9
             IERROR = LREQA - LRLUS
             GOTO 700
           END IF
           CALL CMUMPS_94( N, KEEP(28), IW, LIW, A, LA,
     *           NIRBDU, LRLU, IPTRLU,
     *           IWPOS, IWPOSCB, PTRIST, PTRAST,
     *           STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS )
           COMP = COMP + 1
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB2 compress root2slave:LRLU,LRLUS=',
     *            LRLU, LRLUS
                  IFLAG = -9
                  IERROR = LREQA - LRLUS
                  GOTO 700
           END IF
           IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
              IFLAG  = -8
              IERROR = IWPOS + LREQI - 1 - IWPOSCB
              GOTO 700
           END IF
      END IF
      PTLUST_S(STEP( IROOT )) = IWPOS
      IWPOS           = IWPOS + LREQI
      PTRAST (STEP(IROOT)) = POSFAC
      PTRFAC (STEP(IROOT)) = POSFAC
      POSFAC           = POSFAC + LREQA
      LRLU   = LRLU  - LREQA
      LRLUS  = LRLUS - LREQA
      KEEP(67) = MIN(KEEP(67), LRLUS)
      CALL CMUMPS_471(.FALSE.,.FALSE.,
     *                           LA-LRLUS,0,LREQA,KEEP,KEEP8,LRLU)
      POSHEAD = PTLUST_S( STEP(IROOT))
      IW( POSHEAD + XXI )     = LREQI
      IW( POSHEAD + XXR ) = LREQA
      IW( POSHEAD + XXS ) = S_NOTFREE
      IW( POSHEAD + XSIZE ) = 0
      IW( POSHEAD + 1 + XSIZE ) = NEW_LOCAL_N
      IW( POSHEAD + 2 + XSIZE ) = NEW_LOCAL_M
      IW( POSHEAD + 4 + XSIZE ) = STEP(IROOT)
      IW( POSHEAD + 5 + XSIZE ) = 0
      IF ( MASTER_OF_ROOT ) THEN
        IW( POSHEAD + 3 + XSIZE ) = TOT_ROOT_SIZE
      ELSE
        IW( POSHEAD + 3 + XSIZE ) = 0
      ENDIF
      IF ( KEEP( 50 ) .eq. 0 .OR. KEEP(50) .eq. 2 ) THEN
      OPELIW = OPELIW + ( DBLE(2*TOT_ROOT_SIZE) *
     *         DBLE(TOT_ROOT_SIZE) * DBLE(TOT_ROOT_SIZE ) / DBLE(3)
     *       - 0.5E0 * DBLE( TOT_ROOT_SIZE ) * DBLE( TOT_ROOT_SIZE )
     *       - DBLE( TOT_ROOT_SIZE ) / DBLE( 6 ) )
     *       / DBLE( root%NPROW * root%NPCOL )
      ELSE
      OPELIW = OPELIW + ( DBLE(TOT_ROOT_SIZE) *
     *         DBLE( TOT_ROOT_SIZE) * 
     *         DBLE( TOT_ROOT_SIZE + 1 ) )
     *          / DBLE( 3 * root%NPROW * root%NPCOL )
      END IF
      IF ( PTRIST(STEP( IROOT )) .LT. 0 ) THEN
        A(PTRAST(STEP(IROOT)):PTRAST(STEP(IROOT))+LREQA-1)=CMPLX(ZERO)
        PTRIST(STEP( IROOT )) = 0
        PAMASTER(STEP( IROOT )) = 0
      ELSE IF ( PTRIST( STEP(IROOT)) .eq. 0 ) THEN
        A(PTRAST(STEP(IROOT)):PTRAST(STEP(IROOT))+LREQA-1)=CMPLX(ZERO)
        CALL CMUMPS_35( N, root, IROOT,
     *       A(PTRAST(STEP(IROOT))), NEW_LOCAL_M, NEW_LOCAL_N,
     *       FILS, PTRAIW, PTRARW,
     *       INTARR, DBLARR, KEEP,KEEP8,
     *       MYID )
      ELSE
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + XSIZE )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + XSIZE)
        IF ( TOT_ROOT_SIZE .eq. root%ROOT_SIZE ) THEN
          IF ( LREQA .NE. OLD_LOCAL_M * OLD_LOCAL_N ) THEN
            write(*,*) 'error 1 in PROCESS_ROOT2SLAVE',
     *     OLD_LOCAL_M, OLD_LOCAL_N
            stop
          END IF
          CALL CCOPY(LREQA, A(PAMASTER(STEP(IROOT))), 1,
     *                      A(PTRAST(STEP(IROOT))), 1)
        ELSE
          CALL CMUMPS_96( A( PTRAST(STEP(IROOT))), 
     *        NEW_LOCAL_M,
     *        NEW_LOCAL_N, A( PAMASTER( STEP(IROOT)) ), OLD_LOCAL_M,
     *        OLD_LOCAL_N )
        END IF
        IF ( PTRIST( STEP( IROOT ) ) .NE. 0 ) THEN
           IPOS_SON= PTRIST( STEP(IROOT))
           CALL CMUMPS_152(.FALSE., MYID, N, IPOS_SON,
     *          PAMASTER(STEP(IROOT)),
     *          IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *          IWPOSCB, LA, KEEP,KEEP8
     $         )
           PTRIST(STEP( IROOT ))   = 0
           PAMASTER(STEP( IROOT )) = 0
        END IF
      END IF
 100  CONTINUE
      NBPROCFILS(STEP(IROOT))=NBPROCFILS(STEP(IROOT)) + TOT_CONT_TO_RECV
      IF ( NBPROCFILS(STEP(IROOT)) .eq. 0 ) THEN
         IF (KEEP(201) .NE. 0) THEN
           CALL CMUMPS_580(IERR)
         ENDIF
        CALL CMUMPS_507( N, IPOOL, LPOOL, PROCNODE_STEPS,
     *       SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     *       STEP, IROOT + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL CMUMPS_500(
     $          IPOOL, LPOOL, 
     *          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
 700  CONTINUE
      CALL CMUMPS_44( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE CMUMPS_270
      SUBROUTINE CMUMPS_96
     *( NEW, M_NEW, N_NEW,OLD, M_OLD, N_OLD )
      INTEGER M_NEW, N_NEW, M_OLD, N_OLD
      COMPLEX NEW( M_NEW, N_NEW ), OLD( M_OLD, N_OLD )
      INTEGER J
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      DO J = 1, N_OLD
        NEW( 1: M_OLD, J ) = OLD( 1: M_OLD, J )
        NEW( M_OLD + 1: M_NEW, J ) = CMPLX(ZERO)
      END DO
      NEW( 1: M_NEW,N_OLD + 1: N_NEW ) = CMPLX(ZERO)
      RETURN
      END SUBROUTINE CMUMPS_96
      INTEGER FUNCTION CMUMPS_505(KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      CMUMPS_505 = KEEP(28) + 1 + 3
      RETURN
      END FUNCTION CMUMPS_505
      SUBROUTINE CMUMPS_506(IPOOL, LPOOL, LEAF)
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER LPOOL, LEAF
      INTEGER IPOOL(LPOOL)
      IPOOL(LPOOL-2) = 0
      IPOOL(LPOOL-1) = 0
      IPOOL(LPOOL)   = LEAF-1
      RETURN
      END SUBROUTINE CMUMPS_506
      SUBROUTINE CMUMPS_507
     *           (N, POOL, LPOOL, PROCNODE, SLAVEF,
     *           K28, K76, K80, K47, STEP, INODE)
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N, INODE, LPOOL, K28, SLAVEF, K76, K80, K47
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(K28)
      EXTERNAL CMUMPS_170
      LOGICAL CMUMPS_170, ATM_CURRENT_NODE
      INTEGER NBINSUBTREE, NBTOP, INODE_EFF,POS_TO_INSERT
      INTEGER IPOS1, IPOS2, ISWAP
      INTEGER NODE,J,I
      ATM_CURRENT_NODE = ( K76 == 2 .OR. K76 ==3 .OR.
     $     K76==4 .OR. K76==5)
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      IF (INODE > N ) THEN
        INODE_EFF = INODE - N
      ELSE IF (INODE < 0) THEN
        INODE_EFF = - INODE
      ELSE
        INODE_EFF = INODE
      ENDIF
      IF(((INODE.GT.0).AND.(INODE.LE.N)).AND.(.NOT.
     $     CMUMPS_170(STEP(INODE_EFF),
     *     PROCNODE, SLAVEF))) THEN
         IF ((K80 == 1 .AND. K47 .GE. 1) .OR.
     $     (( K80 == 2 .OR. K80==3 ) .AND.
     *          ( K47 == 4 ))) THEN
            CALL CMUMPS_514(INODE,1)
         ENDIF
      ENDIF
      IF ( CMUMPS_170(STEP(INODE_EFF),
     *                             PROCNODE, SLAVEF) ) THEN
        POOL(NBINSUBTREE + 1 ) = INODE
        NBINSUBTREE = NBINSUBTREE + 1
      ELSE
         POS_TO_INSERT=NBTOP+1
         IF((K76.EQ.4).OR.(K76.EQ.5))THEN
#if defined(NOT_ATM_POOL_SPECIAL)
            J=NBTOP
#else
            IF((INODE.GT.N).OR.(INODE.LE.0))THEN
               DO J=NBTOP,1,-1
                  IF((POOL(LPOOL-2-J).GT.0)
     $                 .AND.(POOL(LPOOL-2-J).LE.N))THEN
                     GOTO 333
                  ENDIF
                  IF ( POOL(LPOOL-2-J) < 0 ) THEN
                     NODE=-POOL(LPOOL-2-J)
                  ELSE IF ( POOL(LPOOL-2-J) > N ) THEN
                     NODE = POOL(LPOOL-2-J) - N
                  ELSE
                     NODE = POOL(LPOOL-2-J)
                  ENDIF
                  IF(K76.EQ.4)THEN
                     IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     $                    DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  IF(K76.EQ.5)THEN
                     IF(COST_TRAV(STEP(NODE)).LE.
     $                    COST_TRAV(STEP(INODE_EFF)))THEN
                        GOTO 333
                     ENDIF
                  ENDIF
                  POS_TO_INSERT=POS_TO_INSERT-1
               ENDDO
               IF(J.EQ.0) J=1
 333           CONTINUE
               DO I=NBTOP,POS_TO_INSERT,-1
                  POOL(LPOOL-2-I-1)=POOL(LPOOL-2-I)
               ENDDO
               POOL(LPOOL-2-POS_TO_INSERT)=INODE
               NBTOP = NBTOP + 1
               GOTO 20
            ENDIF
            DO J=NBTOP,1,-1
               IF((POOL(LPOOL-2-J).GT.0).AND.(POOL(LPOOL-2-J).LE.N))THEN
                  GOTO 888
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
 888        CONTINUE
#endif
            DO I=J,1,-1
#if defined(NOT_ATM_POOL_SPECIAL)
               IF ( POOL(LPOOL-2-I) < 0 ) THEN
                  NODE=-POOL(LPOOL-2-I)
               ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                  NODE = POOL(LPOOL-2-I) - N
               ELSE
                  NODE = POOL(LPOOL-2-I)
               ENDIF
#else
               NODE=POOL(LPOOL-2-I)
#endif
               IF(K76.EQ.4)THEN
                  IF(DEPTH_FIRST_LOAD(STEP(NODE)).GE.
     $                 DEPTH_FIRST_LOAD(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               IF(K76.EQ.5)THEN
                  IF(COST_TRAV(STEP(NODE)).LE.
     $                 COST_TRAV(STEP(INODE_EFF)))THEN
                     GOTO 999
                  ENDIF
               ENDIF
               POS_TO_INSERT=POS_TO_INSERT-1
            ENDDO
            IF(I.EQ.0) I=1
 999        CONTINUE
            DO J=NBTOP,POS_TO_INSERT,-1
               POOL(LPOOL-2-J-1)=POOL(LPOOL-2-J)
            ENDDO
            POOL(LPOOL-2-POS_TO_INSERT)=INODE
            NBTOP = NBTOP + 1
            GOTO 20
         ENDIF
         POOL( LPOOL - 2 - ( NBTOP + 1 ) ) = INODE
         NBTOP = NBTOP + 1
        IPOS1 = LPOOL - 2 - NBTOP
        IPOS2 = LPOOL - 2 - NBTOP + 1
 10     CONTINUE
        IF ( IPOS2 == LPOOL - 2 ) GOTO 20
        IF ( POOL(IPOS1) < 0 ) GOTO 20
        IF ( POOL(IPOS2) < 0 ) GOTO 30
        IF ( ATM_CURRENT_NODE ) THEN
          IF ( POOL(IPOS1) > N ) GOTO 20
          IF ( POOL(IPOS2) > N ) GOTO 30
        END IF
        GOTO 20
 30     CONTINUE
        ISWAP = POOL(IPOS1)
        POOL(IPOS1) = POOL(IPOS2)
        POOL(IPOS2) = ISWAP
        IPOS1 = IPOS1 + 1
        IPOS2 = IPOS2 + 1
        GOTO 10
 20     CONTINUE
      ENDIF
      POOL(LPOOL) = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      RETURN
      END SUBROUTINE CMUMPS_507
      LOGICAL FUNCTION CMUMPS_508(POOL, LPOOL)
      IMPLICIT NONE
      INTEGER LPOOL
      INTEGER POOL(LPOOL)
      INTEGER NBINSUBTREE, NBTOP
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      CMUMPS_508 = (NBINSUBTREE + NBTOP == 0)
      RETURN
      END FUNCTION CMUMPS_508
      SUBROUTINE CMUMPS_509( N, POOL, LPOOL, PROCNODE, SLAVEF,
     *           STEP, INODE, KEEP,KEEP8, MYID, ND )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE, LPOOL, SLAVEF, N
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), POOL(LPOOL), PROCNODE(KEEP(28)),
     *        ND(KEEP(28))
      EXTERNAL CMUMPS_167, CMUMPS_283, CMUMPS_508
      LOGICAL CMUMPS_167, CMUMPS_283, CMUMPS_508
      EXTERNAL CMUMPS_275
      INTEGER CMUMPS_275
      INTEGER NBINSUBTREE, NBTOP, INSUBTREE, INODE_EFF, MYID
      LOGICAL LEFT, ATOMIC_SUBTREE,UPPER,FLAG_MEM,SBTR_FLAG,PROC_FLAG
      INTEGER NODE_TO_EXTRACT,POS_TO_EXTRACT,I,J,K,SON,
     $     NSLAVES_TEMP,POS_TEMP,MIN_PROC
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      IF ( KEEP(76) > 6 .OR. KEEP(76) < 0 ) THEN
         WRITE(*,*) "Error 2 in CMUMPS_509: unknown strategy"
         CALL CMUMPS_ABORT()
      ENDIF
      ATOMIC_SUBTREE =  ( KEEP(76) == 1 .OR. KEEP(76) == 3)
      IF ( CMUMPS_508(POOL, LPOOL) ) THEN
         WRITE(*,*) "Error 1 in CMUMPS_509"
         CALL CMUMPS_ABORT()
      ENDIF
      IF ( .NOT. ATOMIC_SUBTREE ) THEN
         LEFT = (NBTOP == 0)
         IF(.NOT.LEFT)THEN
            IF((KEEP(76).EQ.4).OR.(KEEP(76).EQ.5))THEN
               IF(NBINSUBTREE.EQ.0)THEN
                  LEFT=.FALSE.
               ELSE
                  IF ( POOL(NBINSUBTREE) < 0 ) THEN
                     I = -POOL(NBINSUBTREE)
                  ELSE IF ( POOL(NBINSUBTREE) > N ) THEN
                     I = POOL(NBINSUBTREE) - N
                  ELSE
                     I = POOL(NBINSUBTREE)
                  ENDIF
                  IF ( POOL(LPOOL-2-NBTOP) < 0 ) THEN
                     J = -POOL(LPOOL-2-NBTOP)
                  ELSE IF ( POOL(LPOOL-2-NBTOP) > N ) THEN
                     J = POOL(LPOOL-2-NBTOP) - N
                  ELSE
                     J = POOL(LPOOL-2-NBTOP)
                  ENDIF
                  IF(KEEP(76).EQ.4)THEN
                     IF(DEPTH_FIRST_LOAD(STEP(J)).GE.
     $                    DEPTH_FIRST_LOAD(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
                  IF(KEEP(76).EQ.5)THEN
                     IF(COST_TRAV(STEP(J)).LE.
     $                    COST_TRAV(STEP(I)))THEN
                        LEFT=.TRUE.
                     ELSE
                        LEFT=.FALSE.
                     ENDIF
                  ENDIF
               ENDIF           
            ENDIF
         ENDIF
      ELSE
         IF ( INSUBTREE == 1 ) THEN
            IF (NBINSUBTREE == 0) THEN
               WRITE(*,*) "Error 3 in CMUMPS_509"
               CALL CMUMPS_ABORT()
            ENDIF
            LEFT = .TRUE.
         ELSE
            LEFT = ( NBTOP == 0)
         ENDIF
      ENDIF
 222  CONTINUE
      IF ( LEFT ) THEN
         INODE = POOL( NBINSUBTREE )
         IF(KEEP(81).EQ.2)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
               CALL CMUMPS_561(INODE,POOL,LPOOL,N,
     $              STEP,KEEP,KEEP8,PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     $              PROC_FLAG,MIN_PROC)
               IF(.NOT.SBTR_FLAG)THEN
                  WRITE(*,*)MYID,': ca a change pour moi'
                  LEFT=.FALSE.
                  GOTO 222
               ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
            ENDIF
#endif
         ELSEIF(KEEP(81).EQ.3)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
            IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
               NODE_TO_EXTRACT=INODE
               FLAG_MEM=.FALSE.
               CALL CHECK_MEM_CONST_FOR_POOL(FLAG_MEM)
               IF(FLAG_MEM)THEN
                  CALL CMUMPS_561(INODE,POOL,LPOOL,N,
     $                 STEP,KEEP,KEEP8,
     *                 PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     $                 PROC_FLAG,MIN_PROC)
                  IF(.NOT.SBTR_FLAG)THEN
                     LEFT=.FALSE.
                     GOTO 222
                  ENDIF
               ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
            ENDIF
#endif
         ENDIF
         NBINSUBTREE = NBINSUBTREE - 1
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
         IF ( CMUMPS_167( STEP(INODE_EFF), PROCNODE, SLAVEF) ) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     *           (INSUBTREE.EQ.0))THEN
               CALL CMUMPS_513(.TRUE.)
            ENDIF
            INSUBTREE = 1
         ELSE IF ( CMUMPS_283( STEP(INODE_EFF), PROCNODE,
     *           SLAVEF)) THEN
            IF((KEEP(47).GE.2.AND.KEEP(81).EQ.1).AND.
     *           (INSUBTREE.EQ.1))THEN
               CALL CMUMPS_513(.FALSE.)
            ENDIF
            INSUBTREE = 0
         END IF
      ELSE
         IF (NBTOP < 1 ) THEN
            WRITE(*,*) "Error 5 in CMUMPS_509", NBTOP
            CALL CMUMPS_ABORT()
         ENDIF
         INODE = POOL( LPOOL - 2 - NBTOP )
         IF(KEEP(81).EQ.1)THEN
            CALL CMUMPS_520
     &           (INODE,UPPER,SLAVEF,KEEP,KEEP8,
     *            STEP,POOL,LPOOL,PROCNODE,N)
            IF(UPPER)THEN
               GOTO 666
            ELSE
               NBINSUBTREE=NBINSUBTREE-1
               IF ( CMUMPS_167( STEP(INODE), PROCNODE,
     *              SLAVEF) ) THEN
                  INSUBTREE = 1
               ELSE IF ( CMUMPS_283( STEP(INODE), PROCNODE,
     *                 SLAVEF)) THEN
                  INSUBTREE = 0
               ENDIF
               GOTO 777
            ENDIF
         ENDIF
         IF(KEEP(81).EQ.2)THEN
            CALL CMUMPS_561(INODE,POOL,LPOOL,N,STEP,
     *           KEEP,KEEP8,
     $           PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
            IF(SBTR_FLAG)THEN
               LEFT=.TRUE.
               GOTO 222
            ENDIF
         ELSE
#if defined(POOL_EXTRACT_MNG)
            IF(KEEP(76).EQ.4)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif                 
                  POS_TO_EXTRACT=-1
                  NODE_TO_EXTRACT=-1
                  DO I=NBTOP,1,-1
                     IF(NODE_TO_EXTRACT.LT.0)THEN
                        POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                        INODE_EFF = POOL(LPOOL-2-I)
                        IF ( POOL(LPOOL-2-I) < 0 ) THEN
                           NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                        ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                        ELSE
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                        ENDIF
#else
                        NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                     ELSE
                        IF(DEPTH_FIRST_LOAD(STEP(POOL(LPOOL-2-I))).LT.
     $                       DEPTH_FIRST_LOAD(STEP(NODE_TO_EXTRACT)))
     $                       THEN
                           POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                           INODE_EFF = POOL(LPOOL-2-I)
                           IF ( POOL(LPOOL-2-I) < 0 ) THEN
                              NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                           ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                           ELSE
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                           ENDIF
#else
                           NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                        ENDIF
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
                  INODE = NODE_TO_EXTRACT
#else
                  INODE = INODE_EFF
#endif
                  DO I=POS_TO_EXTRACT,NBTOP
                     IF(I.NE.NBTOP)THEN
                        POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
            IF(KEEP(76).EQ.5)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
                  POS_TO_EXTRACT=-1
                  NODE_TO_EXTRACT=-1
                  DO I=NBTOP,1,-1
                     IF(NODE_TO_EXTRACT.LT.0)THEN
                        POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                        INODE_EFF = POOL(LPOOL-2-I)
                        IF ( POOL(LPOOL-2-I) < 0 ) THEN
                           NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                        ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                        ELSE
                           NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                        ENDIF
#else
                        NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                     ELSE
                        IF(COST_TRAV(STEP(POOL(LPOOL-2-I))).GT.
     $                       COST_TRAV(STEP(NODE_TO_EXTRACT)))THEN
                           POS_TO_EXTRACT=I
#if defined(NOT_ATM_POOL_SPECIAL)
                           INODE_EFF = POOL(LPOOL-2-I)
                           IF ( POOL(LPOOL-2-I) < 0 ) THEN
                              NODE_TO_EXTRACT=-POOL(LPOOL-2-I)
                           ELSE IF ( POOL(LPOOL-2-I) > N ) THEN
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I) - N
                           ELSE
                              NODE_TO_EXTRACT = POOL(LPOOL-2-I)
                           ENDIF
#else
                           NODE_TO_EXTRACT=POOL(LPOOL-2-I)
#endif
                        ENDIF
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)              
                  INODE = NODE_TO_EXTRACT
#else
                  INODE = INODE_EFF
#endif
                  DO I=POS_TO_EXTRACT,NBTOP
                     IF(I.NE.NBTOP)THEN
                        POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
                     ENDIF
                  ENDDO
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
#endif
            IF(KEEP(81).EQ.3)THEN
#if ! defined(NOT_ATM_POOL_SPECIAL)
               IF((INODE.GE.0).AND.(INODE.LE.N))THEN
#endif
                  NODE_TO_EXTRACT=INODE
                  FLAG_MEM=.FALSE.
                  CALL CHECK_MEM_CONST_FOR_POOL(FLAG_MEM)
                  IF(FLAG_MEM)THEN
                     CALL CMUMPS_561(INODE,POOL,LPOOL,N,
     $                    STEP,KEEP,KEEP8,
     *                    PROCNODE,SLAVEF,MYID,SBTR_FLAG,
     $                    PROC_FLAG,MIN_PROC)
                     IF(SBTR_FLAG)THEN
                        LEFT=.TRUE.
                        GOTO 222
                     ENDIF
                  ELSE
                     CALL CLEAN_POOL_MEM_INFO(INODE)
                  ENDIF
#if ! defined(NOT_ATM_POOL_SPECIAL)
               ENDIF
#endif
            ENDIF
         ENDIF
 666     CONTINUE
         NBTOP = NBTOP - 1
         IF((INODE.GT.0).AND.(INODE.LE.N))THEN
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     *           ( KEEP(47) == 4 ))) THEN
               CALL CMUMPS_514(INODE,2)
            ENDIF
         ENDIF
         IF ( INODE < 0 ) THEN
            INODE_EFF = -INODE
         ELSE IF ( INODE > N ) THEN
            INODE_EFF = INODE - N
         ELSE
            INODE_EFF = INODE
         ENDIF
      END IF
 777  CONTINUE
      POOL(LPOOL)     = NBINSUBTREE 
      POOL(LPOOL - 1) = NBTOP
      POOL(LPOOL - 2) = INSUBTREE
      RETURN
      END SUBROUTINE CMUMPS_509
      SUBROUTINE CMUMPS_552(INODE,POOL,LPOOL,N,STEP,
     *     KEEP,KEEP8,
     $     PROCNODE,SLAVEF,MYID,SBTR,FLAG_SAME_PROC,MIN_PROC)
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,MYID,SLAVEF,PROC,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),STEP(N),PROCNODE(KEEP(28))
      INTEGER*8 KEEP8(150)
      INTEGER CMUMPS_275
      EXTERNAL CMUMPS_275
      LOGICAL SBTR,FLAG_SAME_PROC
      INTEGER POS_TO_EXTRACT,NODE_TO_EXTRACT,NBTOP,I,INSUBTREE,
     $     NBINSUBTREE
      DOUBLE PRECISION MIN_COST, TMP_COST
      NBINSUBTREE = POOL(LPOOL)
      NBTOP       = POOL(LPOOL - 1)
      INSUBTREE   = POOL(LPOOL - 2)
      MIN_COST=HUGE(MIN_COST) 
      TMP_COST=HUGE(TMP_COST)
#if ! defined(NOT_ATM_POOL_SPECIAL)
      IF((INODE.GT.0).AND.(INODE.LE.N))THEN
#endif
         POS_TO_EXTRACT=-1
         NODE_TO_EXTRACT=-1
         DO I=NBTOP,1,-1
            IF(NODE_TO_EXTRACT.LT.0)THEN
               POS_TO_EXTRACT=I
               NODE_TO_EXTRACT=POOL(LPOOL-2-I)
               CALL COMPUTE_MAX_MEM(NODE_TO_EXTRACT,TMP_COST,PROC)
               MIN_COST=TMP_COST
               MIN_PROC=PROC
            ELSE
               CALL COMPUTE_MAX_MEM(POOL(LPOOL-2-I),TMP_COST,PROC)
               IF((PROC.NE.MIN_PROC).OR.(TMP_COST.NE.MIN_COST))THEN
                  FLAG_SAME_PROC=.TRUE.
               ENDIF
               IF(TMP_COST.GT.MIN_COST)THEN
                  POS_TO_EXTRACT=I
                  NODE_TO_EXTRACT=POOL(LPOOL-2-I)
                  MIN_COST=TMP_COST
                  MIN_PROC=PROC
               ENDIF
            ENDIF
         ENDDO
         IF((KEEP(47).EQ.4).AND.(NBINSUBTREE.NE.0))THEN
            CALL CMUMPS_554(NBINSUBTREE,INSUBTREE,NBTOP,
     $           MIN_COST,SBTR)
            IF(SBTR)THEN
               WRITE(*,*)MYID,': selecting from subtree'
               RETURN
            ENDIF
         ENDIF
         IF((.NOT.SBTR).AND.(.NOT.FLAG_SAME_PROC))THEN
            WRITE(*,*)MYID,': I must search for a task
     $           to save My friend'
            RETURN
         ENDIF
         INODE = NODE_TO_EXTRACT
         DO I=POS_TO_EXTRACT,NBTOP
            IF(I.NE.NBTOP)THEN
               POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
            ENDIF
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
         CALL CLEAN_POOL_MEM_INFO(INODE)
#if ! defined(NOT_ATM_POOL_SPECIAL)
      ENDIF
#endif
      END SUBROUTINE CMUMPS_552
      SUBROUTINE CMUMPS_561(INODE,POOL,LPOOL,N,STEP,
     *     KEEP,KEEP8,
     $     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER INODE,LPOOL,N,SLAVEF,MYID,MIN_PROC
      INTEGER POOL(LPOOL),KEEP(500),PROCNODE(KEEP(28)),STEP(N)
      INTEGER*8 KEEP8(150)
      LOGICAL SBTR_FLAG,PROC_FLAG
      EXTERNAL CMUMPS_167
      LOGICAL CMUMPS_167
      INTEGER NODE_TO_EXTRACT,I,POS_TO_EXTRACT,NBTOP
      NBTOP= POOL(LPOOL - 1)
      SBTR_FLAG=.FALSE.
      PROC_FLAG=.FALSE.
      CALL CMUMPS_552(INODE,POOL,LPOOL,N,STEP,KEEP,KEEP8,
     $     PROCNODE,SLAVEF,MYID,SBTR_FLAG,PROC_FLAG,MIN_PROC)
      IF(SBTR_FLAG)THEN
         RETURN
      ENDIF
      IF(.NOT.PROC_FLAG)THEN
         NODE_TO_EXTRACT=INODE
         IF((INODE.GE.0).AND.(INODE.LE.N))THEN
            CALL CMUMPS_553(MIN_PROC,POOL,
     $           LPOOL,INODE)
            IF(CMUMPS_167(STEP(INODE),PROCNODE,
     *           SLAVEF))THEN
               WRITE(*,*)MYID,': Extracting from a subtree
     $              for helping',MIN_PROC
               SBTR_FLAG=.TRUE.
               RETURN
            ELSE
               IF(NODE_TO_EXTRACT.NE.INODE)THEN
                  WRITE(*,*)MYID,': Extracting from top
     $                 inode=',inode,'for helping',MIN_PROC
               ENDIF
               CALL CLEAN_POOL_MEM_INFO(INODE)
            ENDIF
         ENDIF
         DO I=1,NBTOP
            IF (POOL(LPOOL-2-I).EQ.INODE)THEN
               GOTO 452
            ENDIF
         ENDDO
 452     CONTINUE
         POS_TO_EXTRACT=I
         DO I=POS_TO_EXTRACT,NBTOP-1
            POOL(LPOOL-2-I)=POOL(LPOOL-2-I-1)
         ENDDO
         POOL(LPOOL-2-NBTOP)=INODE
      ENDIF
      END SUBROUTINE CMUMPS_561
      SUBROUTINE CMUMPS_128(N, NELT, ELTPTR, ELTVAR, LIW,
     *            IKEEP, PTRAR,
     *            IORD, NFSIZ, FILS, FRERE, 
     *            LISTVAR_SCHUR, SIZE_SCHUR,
     *            ICNTL, INFO, KEEP,KEEP8,
     *            ELTNOD, NSLAVES, 
     *            XNODEL, NODEL)
      IMPLICIT NONE
      INTEGER N,NELT,LIW,IORD, SIZE_SCHUR, NSLAVES
      INTEGER PTRAR(N,3), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER ELTPTR(NELT+1) 
      INTEGER XNODEL(N+1), NODEL(ELTPTR(NELT+1)-1)
      INTEGER ELTVAR(ELTPTR(NELT+1)-1)
      INTEGER IKEEP(N,3)
      INTEGER LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER INFO(40), ICNTL(40), KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER ELTNOD(NELT)
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW,IFSON,IN
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG
      INTEGER NZ, allocok, ITEMP
      LOGICAL PROK, NOSUPERVAR
      INTEGER K79REF
      PARAMETER(K79REF=12000000)
      LOGICAL SPLITROOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW2
          INTEGER OPT_METIS_SIZE, NUMFLAG
          PARAMETER(OPT_METIS_SIZE = 8, NUMFLAG = 1)
          INTEGER OPTIONS_METIS(OPT_METIS_SIZE)
      INTEGER IDUM 
      INTEGER IDUMARR(1)
      INTEGER IDUMARR4(1)
      EXTERNAL CMUMPS_197, CMUMPS_130, CMUMPS_131,
     *         CMUMPS_129, CMUMPS_132, 
     *         CMUMPS_133, CMUMPS_134,
     *         CMUMPS_199,
     *         CMUMPS_557, CMUMPS_201
#if defined(OLDDFS)
      EXTERNAL CMUMPS_200
#endif
        ALLOCATE( IW ( LIW ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO( 1 ) = -7
          INFO( 2 ) = LIW
          RETURN
        ENDIF
      MPRINT= ICNTL(3)
      PROK  = (MPRINT.GT.0)
      LP    = ICNTL(1)
      MP    = ICNTL(3)
      LDIAG = ICNTL(4)
      IF (KEEP(60).NE.0) THEN
       NOSUPERVAR=.TRUE.
       IF (IORD.GT.1) IORD = 0
      ELSE
       NOSUPERVAR=.FALSE.
      ENDIF
      IF (IORD == 7) THEN
         IF ( N < 10000 ) THEN
           IORD = 0
         ELSE
#if defined(metis)
           IORD = 5
#else
           IORD = 0
#endif
         ENDIF
      END IF
#if ! defined(metis)
      IF (IORD == 5) IORD = 0
#endif
      IF (KEEP(1).LT.1) KEEP(1) = 1
      NEMIN = KEEP(1)
      IF (LDIAG.LE.2 .OR. MP.LE.0) GO TO 10
      WRITE (MP,99999) N, NELT, LIW, INFO(1)
      K = MIN0(10,NELT+1)
      IF (LDIAG.EQ.4) K = NELT+1
      IF (K.GT.0) WRITE (MP,99998) (ELTPTR(I),I=1,K)
      K = MIN0(10,ELTPTR(NELT+1)-1)
      IF (LDIAG.EQ.4) K = ELTPTR(NELT+1)-1
      IF (K.GT.0) WRITE (MP,99995) (ELTVAR(I),I=1,K)
      K = MIN0(10,N)
      IF (LDIAG.EQ.4) K = N
      IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
      ENDIF
   10 L1 = 1
      L2 = L1 + N
      IF (LIW .LT. 3*N) THEN
          INFO(1)= -2002
          INFO(2) = LIW
      ENDIF
#if defined(metis)
      IF ( IORD == 5 ) THEN
        IF (LIW .LT. N+N+1) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        ENDIF
      ELSE
#endif
      IF (NOSUPERVAR) THEN
        IF ( LIW .LT. 2*N ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ELSE
        IF ( LIW .LT.  4*N+4 ) THEN
          INFO(1)= -2002
          INFO(2) = LIW
          RETURN
        END IF
      ENDIF
#if defined(metis)
      ENDIF
#endif
      IDUM=0
      CALL CMUMPS_258(NELT, N, ELTPTR(NELT+1)-1, ELTPTR, ELTVAR,
     *           XNODEL, NODEL, IW(L1), IDUM, ICNTL)
      IF (IORD.NE.1 .AND. IORD .NE. 5) THEN
        IORD = 0
        IF (NOSUPERVAR) THEN
          CALL CMUMPS_129(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              PTRAR(1,2), IW(L1))
        ELSE
         CALL CMUMPS_130(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              PTRAR(1,2), 4*N+4, IW(L1))
        ENDIF
        LLIW = MAX(NZ,N)
        ALLOCATE( IW2(LLIW), stat = allocok )
        IF (allocok.GT.0) THEN
          INFO(1) = -7
          INFO(2) = LLIW
          RETURN
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL CMUMPS_132(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              IW2, LLIW, PTRAR, PTRAR(1,2),
     *              IW(L1), IWFR)
        ELSE
         CALL CMUMPS_131(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              IW2, LLIW, PTRAR, PTRAR(1,2),
     *              IW(L1), IWFR)
        ENDIF
        IF (NOSUPERVAR) THEN
         CALL CMUMPS_162(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2,
     *   IW(L1), IKEEP,
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3),
     *   LISTVAR_SCHUR, SIZE_SCHUR)
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
         ELSE
           WRITE(*,*) "Internal error in CMUMPS_128",KEEP(60)
           CALL CMUMPS_ABORT()
         ENDIF
        ELSE
         CALL CMUMPS_23(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW2, 
     *   IW(L1), IKEEP, 
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3))
        ENDIF
      ELSE
#if defined(metis)
        IF (IORD.EQ.5) THEN
         IF (PROK) THEN
          WRITE(MPRINT,'(A)') ' Ordering based on METIS '
         ENDIF
         CALL CMUMPS_129(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              PTRAR(1,2), IW(L1))
         LLIW = MAX(NZ,N)
         ALLOCATE( IW2(LLIW), stat = allocok )
         IF (allocok.GT.0) THEN
           INFO(1) = -7
           INFO(2) = LLIW
           RETURN
         ENDIF
         CALL CMUMPS_538(N, NZ, NELT, ELTPTR(NELT+1)-1, 
     *              ELTPTR, ELTVAR, XNODEL, NODEL,
     *              IW2, LLIW, IW(L2), PTRAR(1,2),
     *              IW(L1), IWFR)
          OPTIONS_METIS(1) = 0
          CALL METIS_NODEND(N, IW(L2), IW2(1), NUMFLAG, OPTIONS_METIS,
     *       IKEEP(1,2), IKEEP(1,1) )
           DEALLOCATE(IW2)
        ELSE IF (IORD.NE.1) THEN
          WRITE(*,*) IORD
          WRITE(*,*)  'bad option for ordering'
          CALL CMUMPS_ABORT()
        ENDIF
#endif
       DO K=1,N
         IW(L1+K) = 0
       ENDDO
       DO K=1,N
         IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N)) 
     *    GO TO 40
         IF (IW(L1+IKEEP(K,1)).EQ.1) THEN
          GOTO 40
         ELSE
          IW(L1+IKEEP(K,1)) = 1
         ENDIF
       ENDDO
       CALL CMUMPS_133(N, NZ, NELT, ELTPTR(NELT+1)-1,
     *             ELTPTR, ELTVAR, XNODEL, NODEL, 
     *             IKEEP, PTRAR(1,2), IW(L1))
       LLIW = NZ+N
       ALLOCATE( IW2(LLIW), stat = allocok )
       IF (allocok.GT.0) THEN
         INFO(1) = -7
         INFO(2) = LLIW
         RETURN
       ENDIF
       CALL CMUMPS_134(N, NZ, NELT, ELTPTR(NELT+1)-1,
     *             ELTPTR, ELTVAR, XNODEL, NODEL, 
     *             IKEEP, IW2, LLIW, PTRAR, PTRAR(1,2),
     *             IW(L1), IWFR)
       IF (KEEP(60) == 0) THEN
         ITEMP = 0 
       ELSE
         ITEMP = SIZE_SCHUR
         IF (KEEP(60) == 1) THEN
           KEEP(20) = LISTVAR_SCHUR(1)
         ELSEIF (KEEP(60) == 2 .OR. KEEP(60) == 3 ) THEN
           KEEP(38) = LISTVAR_SCHUR(1)
         ELSE
           WRITE(*,*) "Internal error in CMUMPS_128",KEEP(60)
           CALL CMUMPS_ABORT()
         ENDIF
       ENDIF
       CALL CMUMPS_199(N, PTRAR, IW2, LLIW, IWFR, IKEEP, 
     *    IKEEP(1,2), IW(L1),
     *    IW(L2), NCMPA, ITEMP)
      ENDIF
#if defined(OLDDFS)
      CALL CMUMPS_200(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2),
     *     IKEEP(1,3),
     *     NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, KEEP(60))
#else
      CALL CMUMPS_557(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2),
     *     IKEEP(1,3),
     *     NFSIZ, PTRAR(1,2), 
     *     INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN, 
     *     IW(L2),  KEEP(60), KEEP(20), KEEP(38),
     *     IW2,KEEP(104),IW(L2+N),KEEP(50))
#endif
      DEALLOCATE(IW2)
      IF (KEEP(60).NE.0) THEN
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
     *  INFO(5), KEEP(2),KEEP(50),INFO(20),
     *  KEEP(101), KEEP8(1),KEEP(108))
      IF ( KEEP(53) .NE. 0 ) THEN
        CALL CMUMPS_209( N, FRERE, FILS, NFSIZ, KEEP(20) )
      END IF
      IF ( KEEP(48) == 4 .OR.
     *   ( (KEEP(24).NE.0).AND.(KEEP(10).GT.0) ) ) THEN
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
      SPLITROOT = ((ICNTL(13).GT.0) .AND. (NSLAVES.GE.ICNTL(13)))
      IF (SPLITROOT)
     *   CALL CMUMPS_97(N, FRERE, FILS, NFSIZ,INFO(6),
     *                    NSLAVES, KEEP,KEEP8, SPLITROOT,
     *                    MP, LDIAG)
      IF (KEEP(24)==0) THEN
       CALL CMUMPS_353(N, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2),
     *      PTRAR, ICNTL, INFO, KEEP,KEEP8)
      ELSE
      CALL CMUMPS_348(N, IDUM, IKEEP(1,1), FILS, FRERE,
     *      NFSIZ, IKEEP(1,3), IKEEP(1,2),
     *      IDUMARR4, IDUMARR4, PTRAR, ICNTL, INFO, KEEP,KEEP8)
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
   40 INFO(1) = -4
      INFO(2) = K
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
   90 CONTINUE
      DEALLOCATE(IW)
      RETURN
99999 FORMAT (/'Entering analysis phase with ...'/
     * '                N         NELT       LIW       INFO(1)'/,
     * 9X, I8, I11, I12, I14)
99998 FORMAT ('Element pointers:  ELTPTR()   '/(9X, 7I10))
99995 FORMAT ('Element variables: ELTVAR()   '/(9X, 7I10))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from Analysis   *  INFO(1)=', I3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99982 FORMAT ('Error in permutation array KEEP   INFO(2)=', I3)
      END SUBROUTINE CMUMPS_128
      SUBROUTINE CMUMPS_258( NELT, N, NELNOD, XELNOD, ELNOD,
     $                        XNODEL, NODEL, FLAG, IERROR, ICNTL ) 
      IMPLICIT NONE
      INTEGER NELT, N, NELNOD, IERROR, ICNTL(40)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     $        FLAG(N)
      INTEGER I, J, K, MP, NBERR
      MP = ICNTL(2)
      FLAG(1:N) = 0
      XNODEL(1:N) = 0
      IERROR = 0
      DO I = 1, NELT
        DO K = XELNOD(I), XELNOD(I+1)-1
          J = ELNOD(K)
          IF ( J.LT.1 .OR. J.GT.N ) THEN
            IERROR = IERROR + 1
          ELSE
            IF ( FLAG(J).NE.I ) THEN
              XNODEL(J) = XNODEL(J) + 1
              FLAG(J) = I
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF ( IERROR.GT.0 .AND. MP.GT.0 .AND. ICNTL(4).GE.2 ) THEN
        NBERR = 0
        WRITE(MP,99999)
        DO I = 1, NELT
          DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF ( J.LT.1 .OR. J.GT.N ) THEN
              NBERR = NBERR + 1
              IF (NBERR.LE.10) THEN
                WRITE(MP,'(A,I8,A,I8,A)')
     *          'Element ',I,' variable ',J,' ignored.' 
              ELSE
                GO TO 100
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
  100 CONTINUE
      K = 1
      DO I = 1, N
         K = K + XNODEL(I)
         XNODEL(I) = K
      ENDDO
      XNODEL(N+1) = XNODEL(N)
      FLAG(1:N) = 0
      DO I = 1, NELT
         DO K = XELNOD(I), XELNOD(I+1)-1
            J = ELNOD(K)
            IF (FLAG(J).NE.I) THEN   
              XNODEL(J) = XNODEL(J) - 1
              NODEL(XNODEL(J)) = I
              FLAG(J) = I
            ENDIF
         ENDDO
      ENDDO
      RETURN
99999 FORMAT (/'*** Warning message from subroutine CMUMPS_258 ***')
      END SUBROUTINE CMUMPS_258
      SUBROUTINE CMUMPS_129(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  LEN, FLAG)
      IMPLICIT NONE
      INTEGER N, NELT, NELNOD, NZ
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER  XNODEL(N+1), NODEL(NELNOD),
     $        FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN 
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                LEN(I) = LEN(I) + 1
                LEN(J) = LEN(J) + 1
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0
      DO I = 1,N
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_129
      SUBROUTINE CMUMPS_538(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER LEN(N)
      INTEGER IPE(N+1)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     *          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IWFR = IWFR + LEN(I)
          IPE(I) = IWFR
      ENDDO 
      IPE(N+1)=IPE(N)
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_538
      SUBROUTINE CMUMPS_132(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER  XNODEL(N+1), NODEL(NELNOD), 
     *          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IWFR = IWFR + LEN(I)
        IF (LEN(I).GT.0) THEN
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.LT.J) .AND. (FLAG(J).NE.I)) THEN
                IPE(I) = IPE(I) - 1
                IW(IPE(I)) = J
                IPE(J) = IPE(J) - 1
                IW(IPE(J)) = I
                FLAG(J) = I
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_132
      SUBROUTINE CMUMPS_133(N, NZ, NELT, NELNOD,
     * XELNOD, ELNOD, XNODEL, NODEL, 
     * PERM, LEN, FLAG)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), FLAG(N)
      INTEGER I,J,K1,K2,K3
      FLAG(1:N) = 0
      LEN(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  LEN(I) = LEN(I) + 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NZ = 0  
      DO I = 1,N 
        NZ = NZ + LEN(I)
      ENDDO   
      RETURN  
      END SUBROUTINE CMUMPS_133
      SUBROUTINE CMUMPS_134(N, NZ, NELT, NELNOD,
     * XELNOD, ELNOD, XNODEL, NODEL, 
     * PERM, IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER  XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER PERM(N)
      INTEGER IPE(N), LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD), IW(LW), 
     *          FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 0
      DO I = 1,N
        IWFR = IWFR + LEN(I) + 1
        IPE(I) = IWFR 
      ENDDO
      IWFR = IWFR + 1
      FLAG(1:N) = 0
      DO I = 1,N
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                IF (PERM(J).GT.PERM(I)) THEN
                  IW(IPE(I)) = J
                  IPE(I) = IPE(I) - 1
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DO I = 1,N
        J = IPE(I)
        IW(J) = LEN(I)
        IF (LEN(I).EQ.0) IPE(I) = 0
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_134
      SUBROUTINE CMUMPS_25( MYID, SLAVEF, N, 
     *           PROCNODE, STEP, PTRAIW, PTRARW, 
     *           NELT, FRTPTR, FRTELT, 
     *           KEEP,KEEP8, ICNTL, SYM )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, N, NELT, SYM
      INTEGER KEEP( 500 ), ICNTL( 40 )
      INTEGER*8 KEEP8(150)
      INTEGER PTRAIW( NELT+1 ), PTRARW( NELT+1 )
      INTEGER STEP( N )
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER PROCNODE( KEEP(28) )
      INTEGER MASTER
      PARAMETER (MASTER=0)
      INTEGER CMUMPS_330, CMUMPS_275
      EXTERNAL CMUMPS_330, CMUMPS_275
      INTEGER ELT, I, K, IPTRI, IPTRR, NVAR
      INTEGER TYPE_PARALL, ITYPE, IRANK
      TYPE_PARALL = KEEP(46)
      PTRAIW( 1:NELT ) = 0
      DO I = 1, N
        IF (STEP(I).LT.0) CYCLE
        ITYPE = CMUMPS_330( ABS(STEP(I)), PROCNODE, SLAVEF )
        IRANK = CMUMPS_275( ABS(STEP(I)), PROCNODE, SLAVEF )
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF ( (ITYPE .EQ. 2) .OR.
     *       (ITYPE .EQ. 1 .AND. IRANK .EQ. MYID) ) THEN
          DO K = FRTPTR(I),FRTPTR(I+1)-1
            ELT = FRTELT(K)
            PTRAIW( ELT ) = PTRARW(ELT+1) - PTRARW(ELT)
          ENDDO
        ELSE 
        END IF
      END DO
      IPTRI = 1
      DO ELT = 1,NELT
        NVAR = PTRAIW( ELT )
        PTRAIW( ELT ) = IPTRI
        IPTRI = IPTRI + NVAR
      ENDDO
      PTRAIW( NELT+1 ) = IPTRI
      KEEP( 14 ) = IPTRI - 1
      IF ( .TRUE. ) THEN  
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRAIW( ELT+1 ) - PTRAIW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ELSE
        IF (SYM .EQ. 0) THEN
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + NVAR*NVAR 
          ENDDO
          PTRARW( NELT+1 ) = IPTRR
        ELSE
          IPTRR = 1
          DO ELT = 1,NELT
            NVAR = PTRARW( ELT+1 ) - PTRARW( ELT )
            PTRARW( ELT ) = IPTRR
            IPTRR = IPTRR + (NVAR*(NVAR+1))/2  
          ENDDO 
          PTRARW( NELT+1 ) = IPTRR
        ENDIF
      ENDIF 
      KEEP( 13 ) = IPTRR - 1
      RETURN
      END SUBROUTINE CMUMPS_25
      SUBROUTINE CMUMPS_120( N, NELT, ELTPROC, SLAVEF, PROCNODE )
      IMPLICIT NONE
      INTEGER N, NELT, SLAVEF
      INTEGER PROCNODE( N ), ELTPROC( NELT )
      INTEGER ELT, I, ITYPE, CMUMPS_330, CMUMPS_275
      EXTERNAL CMUMPS_330, CMUMPS_275
      DO ELT = 1, NELT
          I = ELTPROC(ELT)
          IF ( I .NE. 0) THEN
           ITYPE = CMUMPS_330(I,PROCNODE,SLAVEF)
           IF (ITYPE.EQ.1) THEN
             ELTPROC(ELT) = CMUMPS_275(I,PROCNODE,SLAVEF)
           ELSE IF (ITYPE.EQ.2) THEN
             ELTPROC(ELT) = -1
           ELSE
            ELTPROC(ELT) = -2
           ENDIF
          ELSE
           ELTPROC(ELT) = -3
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_120
      SUBROUTINE CMUMPS_153(N, NELT, NELNOD, FRERE, FILS, NA, NE,
     *           XNODEL, NODEL, FRTPTR, FRTELT, ELTNOD) 
      IMPLICIT NONE
      INTEGER  N, NELT, NELNOD
      INTEGER  FRERE(N), FILS(N), NA(N), NE(N)
      INTEGER  FRTPTR(N+1), FRTELT(NELT), ELTNOD(NELT)
      INTEGER  XNODEL(N+1), NODEL(NELNOD) 
      INTEGER TNSTK( N ), IPOOL( N )
      INTEGER I, K, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, III, IN
      TNSTK = NE
      LEAF = 1
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF = 1
        IPOOL(1) = 1
        LEAF = LEAF + 1
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
        DO 20 I=1,NBLEAF-1
           INODE = NA(I)
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
 20     CONTINUE
        INODE = -NA(N)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
        IF (NBLEAF-1.GT.0) THEN
         DO 30 I=1,NBLEAF-1
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 30      CONTINUE
        ENDIF
        INODE = -NA(N-1)-1
        IPOOL(LEAF) = INODE
        LEAF        = LEAF + 1
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
        DO 40 I = 1,NBLEAF
          INODE = NA(I)
          IPOOL(LEAF) = INODE
          LEAF        = LEAF + 1
 40     CONTINUE
      ENDIF
      ELTNOD(1:NELT) = 0
      III = 1
 90   CONTINUE
        IF (III.NE.LEAF) THEN
           INODE=IPOOL(III)
           III = III + 1
        ELSE 
           WRITE(6,*) ' ERROR 1 in file CMUMPS_153 '
           CALL CMUMPS_ABORT()
        ENDIF
 95     CONTINUE 
        IN = INODE
 100    CONTINUE
        DO K = XNODEL(IN),XNODEL(IN+1)-1
          I = NODEL(K)
          IF (ELTNOD(I).EQ.0) ELTNOD(I) = INODE
        ENDDO
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
        IN = INODE
 110    IN = FRERE(IN)
        IF (IN.GT.0) GO TO 110
        IF (IN.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
        ELSE
         IFATH = -IN
        ENDIF
        TNSTK(IFATH) = TNSTK(IFATH) - 1
        IF ( TNSTK(IFATH) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
        ELSE
            GOTO 90
        ENDIF
  115 CONTINUE
      FRTPTR(1:N) = 0
      DO I = 1,NELT
        IF (ELTNOD(I) .NE. 0) THEN
         FRTPTR(ELTNOD(I)) = FRTPTR(ELTNOD(I)) + 1
        ENDIF
      ENDDO
      K = 1
      DO I = 1,N
        K = K + FRTPTR(I)
        FRTPTR(I) = K
      ENDDO
      FRTPTR(N+1) = FRTPTR(N)
      DO K = 1,NELT
        INODE = ELTNOD(K)
        IF (INODE .NE. 0) THEN
         FRTPTR(INODE) = FRTPTR(INODE) - 1
         FRTELT(FRTPTR(INODE)) = K
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_153
      SUBROUTINE CMUMPS_130(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  LEN, LW, IW)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     $        IW(LW)
      INTEGER I,J,K1,K2,K3,LP,NSUP,SUPVAR
      INTEGER INFO44(6)
      EXTERNAL CMUMPS_315
      LP = 6
      CALL CMUMPS_315(N,NELT,XELNOD(NELT+1)-1,ELNOD,XELNOD,
     *           NSUP,IW(3*N+3+1),3*N+3,IW,LP,INFO44)
      IF (INFO44(1) .LT. 0) THEN
        IF (LP.GE.0) WRITE(LP,*) 
     *     'Error return from CMUMPS_315. INFO(1) = ',INFO44(1) 
      ENDIF
      IW(1:NSUP) = 0
      LEN(1:N) = 0
      DO I = 1,N
        SUPVAR = IW(3*N+3+1+I)
        IF (SUPVAR .EQ. 0) CYCLE
        IF (IW(SUPVAR).NE.0) THEN
          LEN(I) = -IW(SUPVAR)
        ELSE
          IW(SUPVAR) = I
        ENDIF
      ENDDO
      IW(N+1:2*N) = 0
      NZ = 0
      DO SUPVAR = 1,NSUP
        I = IW(SUPVAR)
        DO K1 = XNODEL(I),XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2),XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J).GE.0) THEN
                IF ((I.NE.J) .AND. (IW(N+J).NE.I)) THEN
                  IW(N+J) = I
                  LEN(I) = LEN(I) + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        NZ = NZ + LEN(I)
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_130
      SUBROUTINE CMUMPS_131(N, NZ, NELT, NELNOD,
     *  XELNOD, ELNOD, XNODEL, NODEL, 
     *  IW, LW, IPE, LEN, FLAG, IWFR)
      IMPLICIT NONE
      INTEGER N,NZ,NELT,NELNOD,LW,IWFR
      INTEGER XELNOD(NELT+1), ELNOD(NELNOD)
      INTEGER LEN(N)
      INTEGER IPE(N)
      INTEGER XNODEL(N+1), NODEL(NELNOD),
     *          IW(LW), FLAG(N)
      INTEGER I,J,K1,K2,K3
      IWFR = 1
      DO I = 1,N
        IF (LEN(I).GT.0) THEN
          IWFR = IWFR + LEN(I)
          IPE(I) = IWFR
        ELSE
          IPE(I) = 0
        ENDIF
      ENDDO 
      FLAG(1:N) = 0
      DO I = 1,N
        IF (LEN(I).LE.0) CYCLE
        DO K1 = XNODEL(I), XNODEL(I+1)-1
          K2 = NODEL(K1)
          DO K3 = XELNOD(K2), XELNOD(K2+1)-1
            J = ELNOD(K3)
            IF ((J.GE.1) .AND. (J.LE.N)) THEN
              IF (LEN(J) .GT. 0) THEN
                IF ((I.NE.J) .AND. (FLAG(J).NE.I)) THEN
                  IPE(I) = IPE(I) - 1
                  IW(IPE(I)) = J
                  FLAG(J) = I
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CMUMPS_131
      SUBROUTINE CMUMPS_315(N,NELT,NZ,ELTVAR,ELTPTR,NSUP,SVAR,
     +                 LIW,IW,LP,INFO)
      INTEGER LIW,LP,N,NELT,NSUP,NZ
      INTEGER INFO(6)
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER IW(LIW),SVAR(0:N)
      INTEGER FLAG,NEW,VARS
      EXTERNAL CMUMPS_316
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = 0
      INFO(4) = 0
      IF (N.LT.1) GO TO 10
      IF (NELT.LT.1) GO TO 20
      IF (NZ.LT.ELTPTR(NELT+1)-1) GO TO 30
      IF (LIW.LT.6) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      END IF
      NEW = 1
      VARS = NEW + LIW/3
      FLAG = VARS + LIW/3
      CALL CMUMPS_316(N,NELT,ELTPTR,NZ,ELTVAR,SVAR,NSUP,LIW/3-1,
     +           IW(NEW),IW(VARS),IW(FLAG),INFO)
      IF (INFO(1).EQ.-4) THEN
         INFO(4) = 3*N + 3
         GO TO 40
      ELSE
         INFO(4) = 3*NSUP + 3
      END IF
      GO TO 50
   10 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   20 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   30 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
      GO TO 50
   40 INFO(1) = -4
      IF (LP.GT.0) THEN
         WRITE (LP,FMT=9000) INFO(1)
         WRITE (LP,FMT=9010) INFO(4)
      END IF
   50 RETURN
 9000 FORMAT (/3X,'Error message from CMUMPS_315: INFO(1) = ',I2)
 9010 FORMAT (3X,'LIW is insufficient. Upper bound on required work',
     +       'space is ',I8)
      END
      SUBROUTINE CMUMPS_316( N, NELT, ELTPTR, NZ, ELTVAR,
     *           SVAR, NSUP, MAXSUP, NEW, VARS, FLAG, INFO )
      INTEGER MAXSUP,N,NELT,NSUP,NZ
      INTEGER ELTPTR(NELT+1),ELTVAR(NZ)
      INTEGER INFO(6)
      INTEGER FLAG(0:MAXSUP), NEW(0:MAXSUP),SVAR(0:N),
     *          VARS(0:MAXSUP)
      INTEGER I,IS,J,JS,K,K1,K2
      DO 10 I = 0,N
         SVAR(I) = 0
   10 CONTINUE
      VARS(0) = N + 1
      NEW(0) = -1
      FLAG(0) = 0
      NSUP = 0
      DO 40 J = 1,NELT
         K1 = ELTPTR(J)
         K2 = ELTPTR(J+1) - 1
         DO 20 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
               INFO(2) = INFO(2) + 1
               GO TO 20
            END IF
            IS = SVAR(I)
            IF (IS.LT.0) THEN
               ELTVAR(K) = 0
               INFO(3) = INFO(3) + 1
               GO TO 20
            END IF
            SVAR(I) = SVAR(I) - N - 2
            VARS(IS) = VARS(IS) - 1
   20    CONTINUE
         DO 30 K = K1,K2
            I = ELTVAR(K)
            IF (I.LT.1 .OR. I.GT.N) GO TO 30
            IS = SVAR(I) + N + 2
            IF (FLAG(IS).LT.J) THEN
               FLAG(IS) = J
               IF (VARS(IS).GT.0) THEN
                  NSUP = NSUP + 1
                  IF (NSUP.GT.MAXSUP) THEN
                     INFO(1) = -4
                     RETURN
                  END IF
                  VARS(NSUP) = 1
                  FLAG(NSUP) = J
                  NEW(IS) = NSUP
                  SVAR(I) = NSUP
               ELSE
                  VARS(IS) = 1
                  NEW(IS) = IS
                  SVAR(I) = IS
               END IF
            ELSE
               JS = NEW(IS)
               VARS(JS) = VARS(JS) + 1
               SVAR(I) = JS
            END IF
   30    CONTINUE
   40 CONTINUE
      RETURN
      END
      SUBROUTINE CMUMPS_125(
     *           NUMELT, LIST_ELT, 
     *           MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF, PERM, 
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, NELT, 
     *           IW, LIW, 
     *           INTARR, LINTARR, ITLOC, FILS, FRERE_STEPS, 
     *           KEEP28, SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     *        LINTARR, NUMELT, NFRONT_EFF, KEEP28
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), PIMASTER(KEEP28),
     *  PTRAIW(NELT+1), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP28),
     *        NBPROCFILS(KEEP28), PERM(N)
      INTEGER INTARR(LINTARR)
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     *        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     *        LAST_J_ASS, JMIN, MIN_PERM
      INTEGER ELTI
      INCLUDE 'mumps_headers.h'
      LOGICAL LEVEL1
      INTEGER allocok
      INTEGER , ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      IF (NUMSTK.GT.0) THEN
        ALLOCATE(PTTRI(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF
        ALLOCATE(PTLAST(NUMSTK), stat=allocok)
        IF (allocok .GT. 0) THEN
         IFLAG = -13
         GOTO 800
        ENDIF        
      ENDIF
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1
      ICT11  = IOLDP2 + NFRONT
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN)        = INEW1
       IW(INEW)         = IN
       IW(INEW+NFRONT)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
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
      MIN_PERM = N + 1
      JMIN     = -1
      DO IELL = 1, NUMSTK 
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
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
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT) THEN 
       DO IELL = 1,NUMELT
        ELTI = LIST_ELT(IELL)
         J1= PTRAIW(ELTI)
         J2= PTRAIW(ELTI+1)-1
         DO JJ=J1,J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
       ENDDO
      END IF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
       IW( IOLDP2+NASS1+1 : IOLDP2+NFRONT_EFF ) = 
     *  IW( ICT11+NASS1+1 : ICT11+NFRONT_EFF ) 
      ELSE
          CALL CMUMPS_308( N, PERM, 
     *            IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
          CALL CMUMPS_309( N, NASS1, PERM, ITLOC,
     *               IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     *               IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     *               IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
      END IF
      IF ( NFRONT_EFF .NE. NEWEL1_SAVE ) THEN
        IP1 = IOLDPS + NASS1 + HF
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
          NELIM = IW(J2 + 1 +XSIZE)
          NPIVS = IW(J2 + 3 +XSIZE)
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = IW(J2 + 5 +XSIZE)
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
      DO IELL=1,NUMELT
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          INTARR(JJ) = ITLOC(J)
       END DO
      ENDDO
        K1 = IOLDPS + HF + NUMORG
        K2 = K1 + NFRONT_EFF - 1 + NASS
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
  800 CONTINUE
      IF (ALLOCATED(PTTRI)) DEALLOCATE(PTTRI)
      IF (ALLOCATED(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE CMUMPS_125
      SUBROUTINE CMUMPS_124(
     *           NUMELT, LIST_ELT,
     *           MYID, INODE, N, IOLDPS,
     *           HF, NFRONT, NFRONT_EFF,
     *           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *           IFSON, STEP, PIMASTER, PTRAIW, NELT,
     *           IW, LIW, 
     *           INTARR, LINTARR, ITLOC, FILS, FRERE_STEPS, 
     *           KEEP28,
     *           SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      IMPLICIT NONE
      INTEGER NELT, INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     *        NUMSTK, NUMORG, IWPOSCB, IFSON, MYID, IFLAG,
     *        LINTARR, NUMELT, NFRONT_EFF, KEEP28
      INTEGER LIST_ELT(*)
      INTEGER STEP(N), 
     * PIMASTER(KEEP28),
     *  PTRAIW(NELT+1), IW(LIW), 
     *        ITLOC(N), FILS(N), FRERE_STEPS(KEEP28),
     *        NBPROCFILS(KEEP28)
      INTEGER INTARR(LINTARR)
      LOGICAL SON_LEVEL2, NIV1
      INTEGER NEWEL, INEW, IOLDP2, INEW1,
     *        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     *        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     *        NROWS, HS, IP1, IP2, K1, K2, 
     *        I, K, ELTI
      LOGICAL LEVEL1
      INCLUDE 'mumps_headers.h'
      SON_LEVEL2 = .FALSE.
      IOLDP2 = IOLDPS + HF - 1    ! XSIZE already in HF
      NEWEL = IOLDP2 + NASS1
      NFRONT_EFF = NASS1
      IN = INODE
      INEW = IOLDPS + HF 
      INEW1 = 1
      DO WHILE (IN.GT.0)
       ITLOC(IN) = INEW1
       IW(INEW)  = IN
       INEW1     = INEW1 + 1
       INEW      = INEW + 1
       IN = FILS(IN)
      END DO
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
          HS = NSLSON + 6 +XSIZE
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
          ISON = FRERE_STEPS(STEP(ISON))
  100   CONTINUE
      ENDIF
      DO IELL=1,NUMELT 
       ELTI = LIST_ELT(IELL)
       J1= PTRAIW(ELTI)
       J2= PTRAIW(ELTI+1)-1
       DO JJ=J1,J2
          J = INTARR(JJ)
          IF (ITLOC(J) .EQ. 0) THEN
              NEWEL = NEWEL + 1
              NFRONT_EFF = NFRONT_EFF + 1
              IW(NEWEL)  = J
              INTARR(JJ) = NFRONT_EFF
              ITLOC(J)   = NFRONT_EFF
          ELSE
              INTARR(JJ) = ITLOC(J)
          ENDIF
       END DO
      ENDDO
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
      END SUBROUTINE CMUMPS_124
      SUBROUTINE CMUMPS_36( COMM_LOAD, ASS_IRECV,
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, 
     *    NIRBDU, FILS, FRERE_STEPS, MAXFRW, root,
     *    OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC, PTRAST, 
     *    STEP, PIMASTER, PAMASTER,PTRARW, 
     *    PTRAIW, ITLOC, NSTEPS, SON_LEVEL2,
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     *    ICNTL, KEEP,KEEP8,INTARR,DBLARR, 
     *
     *    NSTK_S,NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     *    PERM, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *    )
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER NELT,N,LIW,LA,NIRBDU,NSTEPS
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP, IERR_MPI
      INTEGER IDUMMY(1)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)), PERM(N), 
     *        FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     *        PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     * PAMASTER(KEEP(28))
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      LOGICAL SON_LEVEL2
      COMPLEX A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(MAX(1,KEEP(14)))
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NBPROCFILS(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER ETATASS
      INCLUDE 'mumps_headers.h'
      LOGICAL COMPRESSCB
      INTEGER LCB
      INTEGER CMUMPS_275, CMUMPS_330
      EXTERNAL CMUMPS_275, CMUMPS_330
      INTEGER LP, HS, HF
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER SIZFI, SIZFR, NCB
      INTEGER LAPOS2,J1,J2
      INTEGER NCOL, NROW, NCOLS, NROWS, LDA_SON
      INTEGER NELIM,JJ,JJ1,JJ2,J3,
     *        IORG
      INTEGER IACHK,JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4,
     *        NUMELT, ELBEG
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12,
     *        AII, ICT21, J
      INTEGER NSLAVES, NSLSON, NPIVS, NPIV_ANA, NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, K, I
      LOGICAL FLAG, BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER NCBSON
      LOGICAL FREE, SAME_PROC
      INTRINSIC REAL
      REAL ZERO
      DATA ZERO /0.0E0/
      LOGICAL  CMUMPS_167, SSARBR
      EXTERNAL CMUMPS_167
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL CMUMPS_170
      LOGICAL CMUMPS_170
      NFS4FATHER = -1
      ETATASS    = 0  ! full assembly, only option for element format
      COMPRESSCB=.FALSE.
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = CMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF (LEVEL.NE.1) THEN 
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL CMUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES + XSIZE
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .ne. 0 ) THEN
        ELBEG  = FRT_PTR(INODE)
      ELSE
        ELBEG  = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NPIV_ANA=NUMORG
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + XSIZE)
         ISON = FRERE_STEPS(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS
      NASS1 = NASS + NUMORG
      LREQ = 2 * NFRONT + HF   ! HF includes XSIZE
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL CMUMPS_94(N, KEEP(28),
     *        IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
      IF (KEEP(50).EQ.0) THEN
        CALL  CMUMPS_124(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF, NFRONT, NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      ELSE
        CALL CMUMPS_125( 
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
        IF (IFLAG.LT.0) GOTO 300
      END IF
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
           IF(CMUMPS_170(STEP(INODE),PROCNODE_STEPS,
     $          SLAVEF))THEN
              NPIV=NASS1-(NFRONT_EFF-ND(STEP(INODE)))
              CALL CMUMPS_511(ND(STEP(INODE)),NPIV,NPIV,
     *                                 KEEP(50),1,FLOP1)             
              NPIV=NPIV_ANA
              CALL CMUMPS_511(ND(STEP(INODE)),NPIV,NPIV,
     *                                 KEEP(50),1,FLOP1_EFF)
              CALL CMUMPS_190(0,.FALSE.,FLOP1-FLOP1_EFF,
     $             KEEP,KEEP8)
           ENDIF
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        LREQ = 2 * NFRONT + HF   ! HF includes XSIZE
        ELSE
         Write(*,*) ' ERROR 1 during ass_niv1_ELT'
         GOTO 270
        ENDIF
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = MAX0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      LAELL = NFRONT * NFRONT
      IF (LRLU .LT. LAELL) THEN
        IF (LRLUS .LT. LAELL) THEN
          GOTO 280
        ELSE
          CALL CMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *         NIRBDU, LRLU, IPTRLU,
     *         IWPOS, IWPOSCB, PTRIST, PTRAST,
     *         STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
          COMP = COMP + 1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv1.F'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          END IF
        END IF
      END IF
      LRLU = LRLU - LAELL
      LRLUS = LRLUS - LAELL
      KEEP(67) = MIN(LRLUS, KEEP(67))
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL
      SSARBR=CMUMPS_167(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      CALL CMUMPS_471(SSARBR,.FALSE.,LA-LRLUS,0,LAELL,
     *    KEEP,KEEP8,
     $     LRLU)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL - 1
      A(POSELT:LAPOS2) = CMPLX(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = CMPLX(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, NFRONT -1
          A(APOS:APOS+JJ) = CMPLX(ZERO)
          APOS = APOS + NFRONT
        END DO
      END IF
#endif
      NASS = NASS1
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS+XXI) = LREQ  ! size of header + integer lists
      IW(IOLDPS+XXR) = LAELL ! number of reals entries in front
      IW(IOLDPS+XXS) = -9999 ! Status of node  to be defined
      IW(IOLDPS+XSIZE)   = NFRONT
      IW(IOLDPS+XSIZE+ 1) = 0
      IW(IOLDPS+XSIZE + 2) = -NASS1
      IW(IOLDPS+XSIZE + 3) = -NASS1
      IW(IOLDPS+XSIZE + 4) = STEP(INODE)
      IW(IOLDPS+XSIZE+5)   = NSLAVES
      IF (NUMSTK.NE.0) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          LSTK      = IW(ISTCHK+XSIZE)
          NELIM     = IW(ISTCHK + 1+XSIZE)
          NPIVS     = IW(ISTCHK + 3+XSIZE)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = IW(ISTCHK + 5+XSIZE)
          HS        = 6 + NSLSON + XSIZE
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
             COMPRESSCB =
     *           ( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
             COMPRESSCB =  ( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW( ISTCHK + 2+XSIZE)
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          J1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1) THEN
           J2 = J1 + LSTK - 1
           SIZFR  = LSTK*LSTK
           IF (COMPRESSCB) SIZFR = (LSTK*(LSTK+1))/2
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR = NELIM * LSTK
           ELSE
             SIZFR = NELIM * NELIM
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          OPASSW = OPASSW + DBLE(SIZFR)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT
            IF (J2.GE.J1) THEN
              DO 170 JJ = J1, J2
                APOS = POSEL1 + IW(JJ) * NFRONT
                DO 160 JJ1 = 1, LSTK
                  JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
                  A(JJ2) = A(JJ2) + A(IACHK + JJ1 - 1)
  160           CONTINUE
                IACHK = IACHK + LSTK
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = SIZFR
            ELSE
              LCB = LDA_SON*( J2 - J1 + 1)
            ENDIF
            CALL CMUMPS_624(
     *           A( PTRAST(STEP( INODE )) ), NFRONT, NASS1,
     *           A( IACHK ), LDA_SON, LCB,
     *           IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     *           COMPRESSCB
     *          )
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF (SAME_PROC) THEN
             IF (KEEP(50).NE.0) THEN
              J2 = J1 + LSTK - 1
              DO JJ = J1, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               IW(JJ) = IW(JJ - NROWS)
              END DO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = IW(JJ) + ICT11
                 IW(JJ) = IW(JPOS)
                END DO
              ENDIF
             ENDIF
           ENDIF
           IF ( SAME_PROC ) THEN
               PTRIST(STEP( ISON )) = -99999999
           ELSE
               PIMASTER(STEP( ISON )) = -99999999
           ENDIF
           CALL CMUMPS_152(SSARBR, MYID, N, ISTCHK,
     *     IACHK,
     *     IW, LIW, NIRBDU, LRLU, LRLUS, IPTRLU,
     *     IWPOSCB, LA, KEEP,KEEP8
     *     )
          ELSE
           PDEST = ISTCHK + 6 + XSIZE
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL CMUMPS_49( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
               CALL CMUMPS_210( COMM_LOAD, ASS_IRECV, 
     *                           BUFR, LBUFR, LBUFR_BYTES,
     *                           INODE, ISON, NSLAVES, IDUMMY,
     *                           NFRONT, NASS1,NFS4FATHER,
     *         TROW_SIZE, IW( INDX ),
     *         PROCNODE_STEPS,
     *         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *         LRLUS, N, IW,
     *         LIW, A, LA, NIRBDU,
     *         PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     *         PIMASTER, PAMASTER, NSTK_S, COMP,
     *         IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     *         LEAF, ICNTL, KEEP,KEEP8, root,
     *         OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *         INTARR, DBLARR, ND, FRERE_STEPS,
     *         NELT+1, NELT, FRT_PTR, FRT_ELT,
     *   
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     *         )
               IF ( IFLAG .LT. 0 ) GOTO 500
               EXIT
             ENDIF
           END DO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + XSIZE
            CALL  CMUMPS_71( INODE, NFRONT, 
     *       NASS1, NFS4FATHER,ISON, MYID,
     *       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     *       COMM, IERR, IW(PDEST), NSLSON, 
     *       SLAVEF, 
     *       KEEP,KEEP8, STEP, N, 
     *       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
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
     *         NELT+1, NELT, FRT_PTR, FRT_ELT, 
     *         ISTEP_TO_INIV2, TAB_POS_IN_PERE )
               IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
  210   ISON = FRERE_STEPS(STEP(ISON))
  220 CONTINUE
      END IF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          AINPUT    = AII + II - J1
          ICT12 = POSELT + (I-1) * NFRONT
          DO JJ=J1,J2
           APOS2 = ICT12 + INTARR(JJ) - 1
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
           AINPUT = AINPUT + SIZE_ELTI
          END DO
         ELSE
          ICT12 = POSELT - NFRONT + I - 1 
          ICT21 = POSELT + (I-1)*NFRONT - 1 
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF (I.LT.J) THEN
              APOS2 = ICT12 + J*NFRONT
           ELSE
              APOS2 = ICT21 + J 
           ENDIF
           A(APOS2) = A(APOS2) + DBLARR(AII)
           AII = AII + 1
          END DO
         END IF
        END DO
      END DO
      GOTO 500
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_36'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_36'
      ENDIF
      IFLAG = -9
      IERROR = LAELL - LRLUS
      GOTO 500
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *  ' FAILURE, SEND BUFFER TOO SMALL DURING CMUMPS_36'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *  ' FAILURE, RECV BUFFER TOO SMALL DURING CMUMPS_36'
      ENDIF
      IFLAG = -20
      LREQ = NCBSON + 6+NSLSON+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION DURING CMUMPS_36'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK 
  490 CALL  CMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_36
      SUBROUTINE CMUMPS_37( COMM_LOAD, ASS_IRECV,
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, IFLAG,
     *    IERROR, ND, NIRBDU, FILS, FRERE_STEPS,
     *    CAND, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *    MAXFRW, root,
     *    OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     *    PTRAIW, ITLOC, NSTEPS, 
     *    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     *    ICNTL, KEEP,KEEP8,INTARR,DBLARR, 
     *    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     *    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     *    PERM,
     *    MEM_DISTRIB)
      USE CMUMPS_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER IERR, STATUS( MPI_STATUS_SIZE )
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER NELT, N,LIW,LA,NIRBDU,NSTEPS, NBFIN
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      INTEGER IFLAG,IERROR,INODE,MAXFRW,
     *        LPOOL, LEAF, 
     *        LRLU, IPTRLU,IWPOS, LRLUS,
     *        POSFAC, IWPOSCB, COMP, SLAVEF
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER IW(LIW), ITLOC(N),
     *        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)),
     *        FILS(N), FRERE_STEPS(KEEP(28)),
     *        PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *        PTRAST(KEEP(28)), STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)),
     *        NSTK_S(KEEP(28)), PERM(N)
      INTEGER   CAND(SLAVEF+1,MAX(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      COMPLEX A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(MAX(1,KEEP(14)))
      COMPLEX DBLARR(MAX(1,KEEP(13)))
      INTEGER MYID, COMM
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF, HF_OLD, NSLAVES_OLD,NCBSON
      INTEGER NCBSON_MAX
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      LOGICAL COMPRESSCB
      INTEGER LCB
      INTEGER NFS4FATHER
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ,LAELL
      INTEGER NCB, IERR_MPI
      INTEGER LAPOS2,J1,J2,J3
      INTEGER NELIM,JJ,JJ1,JJ2,NPIVS,NCOLS,NROWS,
     *        IORG
      INTEGER IACHK, LDAFS, LDA_SON
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12, ICT21
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, K, ELBEG, NUMELT, I, J, AII
      LOGICAL FLAG, SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER NSLAVES_less, ITEMP, NMB_OF_CAND
      logical :: force_cand
      INTEGER APOSMAX
      REAL  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST
      INTEGER IZERO 
      INTEGER IDUMMY(1)
      INTEGER PDEST1(1)
      INTEGER ETATASS
      PARAMETER( IZERO = 0 )
      INTEGER CMUMPS_275, CMUMPS_330
      EXTERNAL CMUMPS_275, CMUMPS_330
      INTRINSIC REAL
      REAL DATA_SIZE
      REAL ZERO
      DATA ZERO /0.0E0/
      COMPRESSCB=.FALSE.
      ETATASS = 0  ! full assembly (only option available here)
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .NE. 0 ) THEN
        ELBEG = FRT_PTR(INODE)
      ELSE
        ELBEG = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      NCBSON_MAX  = 0
      DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         IF ( KEEP(48)==5 .AND. CMUMPS_330(STEP(ISON),
     *        PROCNODE_STEPS,SLAVEF) .EQ. 1) THEN
            NCBSON_MAX = MAX(NCBSON_MAX,IW(PIMASTER(STEP(ISON))+XSIZE))
         END IF
         NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + XSIZE)
         ISON = FRERE_STEPS(STEP(ISON))
      END DO
      NFRONT = ND(STEP(INODE)) + NASS
      MAXFRW = MAX0(MAXFRW, NFRONT)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      DATA_SIZE = NASS1 * NCB
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      end if
      IF (force_cand) THEN
         INIV2 = ISTEP_TO_INIV2( STEP( INODE ))
         SIZE_TMP_SLAVES_LIST = CAND( SLAVEF+1, INIV2 )
      ELSE
         INIV2 = 1
         SIZE_TMP_SLAVES_LIST = SLAVEF - 1
      ENDIF
      ALLOCATE(TMP_SLAVES_LIST(SIZE_TMP_SLAVES_LIST),stat=allocok)
      IF (allocok > 0 ) THEN
        GOTO 265
      ENDIF
      CALL CMUMPS_472( NCBSON_MAX, SLAVEF,
     *     KEEP,KEEP8,ICNTL,
     *     CAND(1,INIV2),
     *     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     *     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *     TMP_SLAVES_LIST,
     *     SIZE_TMP_SLAVES_LIST,INODE )
      HF   = NSLAVES + 6 + XSIZE
      LREQ = 2 * NFRONT + HF
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL CMUMPS_94(N, KEEP(28),
     *        IW, LIW, A, LA, NIRBDU,
     *        LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
          COMP = COMP+1
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .FALSE.
      IF (KEEP(50).EQ.0) THEN
        CALL  CMUMPS_124(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF, NFRONT,NFRONT_EFF,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS, KEEP(28),
     *        SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
      ELSE
        CALL CMUMPS_125(
     *        NUMELT, FRT_ELT(ELBEG),
     *        MYID, INODE, N, IOLDPS, HF,
     *        NFRONT, NFRONT_EFF, PERM,
     *        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     *        IFSON, STEP, PIMASTER, PTRAIW, NELT, IW, LIW,
     *        INTARR, KEEP(14), ITLOC, FILS, FRERE_STEPS,
     *        KEEP(28), SON_LEVEL2, NIV1, NBPROCFILS, IFLAG)
        IF (IFLAG.LT.0) GOTO 250
      ENDIF
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            CALL CMUMPS_472( NCBSON_MAX, 
     *      SLAVEF, KEEP,KEEP8,ICNTL,
     *      CAND(1,INIV2),
     *      MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     *      TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *      TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST,INODE )
            HF = NSLAVES + 6 + XSIZE
            IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF)) -
     &                   (NSLAVES_OLD - NSLAVES)
            IF (NSLAVES_OLD .NE. NSLAVES) THEN
              IF (NSLAVES_OLD > NSLAVES) THEN
               IW(IOLDPS+HF: IOLDPS+HF+2*NFRONT_EFF-1) =
     &         IW(IOLDPS+HF_OLD: IOLDPS+HF_OLD+2*NFRONT_EFF-1)
              ELSE
               IF (IWPOS - 1 > IWPOSCB ) GOTO 270
               DO JJ=2*NFRONT_EFF-1, 0, -1
                 IW(IOLDPS+HF+JJ) = IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              END IF
            END IF
            NFRONT = NFRONT_EFF
            LREQ = 2 * NFRONT + HF
        ELSE
          Write(*,*) ' ERROR 2 during ass_niv2'
          GOTO 270
        ENDIF
      ENDIF
      MAXFRW = MAX0(MAXFRW, NFRONT)
      PTLUST_S(STEP(INODE)) = IOLDPS
      IW(IOLDPS + 1+XSIZE) = 0
      IW(IOLDPS + 2+XSIZE) = -NASS1
      IW(IOLDPS + 3+XSIZE) = -NASS1
      IW(IOLDPS + 4+XSIZE) = STEP(INODE)
      IW(IOLDPS+XSIZE)   = NFRONT
      IW(IOLDPS+5+XSIZE) = NSLAVES
      IW(IOLDPS+6+XSIZE:IOLDPS+5+NSLAVES+XSIZE)=
     &                     TMP_SLAVES_LIST(1:NSLAVES)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE)
      IF (KEEP(73) .EQ. 0) THEN
#endif
#endif
        CALL CMUMPS_461(MYID, SLAVEF, COMM_LOAD,
     *     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *     NASS1, KEEP, KEEP8, IW(IOLDPS+6+XSIZE), NSLAVES,INODE)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE) 
      ENDIF
#endif
#endif
      IF(KEEP(86).EQ.1)THEN
         IF(mod(KEEP(24),2).eq.0)THEN
            CALL CMUMPS_533(SLAVEF,CAND(SLAVEF+1,INIV2),
     *           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, NSLAVES,INODE)
         ELSEIF((KEEP(24).EQ.0).OR.(KEEP(24).EQ.1))THEN
            CALL CMUMPS_533(SLAVEF,SLAVEF-1,
     *           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     *           NASS1, KEEP,KEEP8,TMP_SLAVES_LIST, NSLAVES,INODE)
         ENDIF
      ENDIF         
      DEALLOCATE(TMP_SLAVES_LIST)
      IF (KEEP(50).EQ.0) THEN
        LAELL = NASS1 * NFRONT
        LDAFS = NFRONT
      ELSE
        LAELL = NASS1**2
        IF (KEEP(219).NE.0) THEN
          IF(KEEP(50) .EQ. 2) LAELL = LAELL+NASS1
        ENDIF
        LDAFS = NASS1
      ENDIF
      IF (LRLU .LT. LAELL) THEN
        IF (LRLUS .LT. LAELL) THEN
          GOTO 280
        ELSE
          CALL CMUMPS_94(N, KEEP(28), IW, LIW, A, LA,
     *        NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, PTRIST, PTRAST,
     *        STEP, PIMASTER, PAMASTER, ITLOC,KEEP(216),LRLUS)
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress ass..mpi51f_niv2'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL
      LRLUS = LRLUS - LAELL
       KEEP(67) = MIN(LRLUS, KEEP(67))
      POSELT = POSFAC
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      POSFAC = POSFAC + LAELL
      IW(IOLDPS+XXI)   = LREQ  ! size of header + integer lists
      IW(IOLDPS+XXR) = LAELL ! number of reals entries in front
      IW(IOLDPS+XXS) = -9999 ! Status of node  to be defined
      CALL CMUMPS_471(.FALSE.,.FALSE.,LA-LRLUS,0,LAELL,
     * KEEP,KEEP8,
     $LRLU)
      POSEL1 = POSELT - LDAFS
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL - 1
      A(POSELT:LAPOS2) = CMPLX(ZERO)
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL - 1
        A(POSELT:LAPOS2) = CMPLX(ZERO)
      ELSE
        APOS = POSELT
        DO JJ = 0, LDAFS - 1
          A(APOS:APOS+JJ) = CMPLX(ZERO)
          APOS = APOS + LDAFS
        END DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          A(APOS:APOS+LDAFS-1)=CMPLX(ZERO)
        ENDIF
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + XSIZE + 1)
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK + XSIZE)
          NPIVS   = IW(ISTCHK + XSIZE + 3)
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + XSIZE + 5)
          HS      = 6 + XSIZE + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LE.IWPOS)
          IF ( SAME_PROC ) THEN
            COMPRESSCB=( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
            COMPRESSCB=( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + XSIZE + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          OPASSW = OPASSW + DBLE(NELIM*LSTK)
          J1 = ISTCHK + HS + NROWS + NPIVS
          J2 = J1 + NELIM - 1
          IACHK = PAMASTER(STEP(ISON))
          IF (KEEP(50).eq.0) THEN
           DO 170 JJ = J1, J2
            APOS = POSEL1 + IW(JJ) * LDAFS
            DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + IW(J1 + JJ1 - 1) - 1
              A(JJ2) = A(JJ2) + A(IACHK + JJ1 - 1)
  160       CONTINUE
            IACHK = IACHK + LSTK
  170      CONTINUE
          ELSE
            IF (NSLSON.EQ.0) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = (NELIM*(NELIM+1))/2
            ELSE
              LCB = LDA_SON*( NELIM )
            ENDIF
            CALL CMUMPS_624(
     *           A( POSELT ), LDAFS, NASS1,
     *           A( IACHK ), LDA_SON, LCB,
     *           IW( J1 ), NELIM, NELIM, ETATASS,
     *          COMPRESSCB)
          ENDIF
  210     ISON = FRERE_STEPS(STEP(ISON))
  220   CONTINUE
      ENDIF
      APOSMAX = POSELT + NASS1*NASS1
      IF (KEEP(219).NE.0) THEN
        IF (KEEP(50).EQ.2) THEN
          A( APOSMAX: APOSMAX+NASS1-1)=ZERO
        ENDIF
      ENDIF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          IF (I.LE.NASS1) THEN
           AINPUT    = AII + II - J1
           ICT12 = POSELT + (I-1) * LDAFS
           DO JJ=J1,J2
            APOS2 = ICT12 + INTARR(JJ) - 1
            A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            AINPUT = AINPUT + SIZE_ELTI
           END DO
          ENDIF
         ELSE
          ICT12 = POSELT - LDAFS + I - 1
          ICT21 = POSELT + (I-1)*LDAFS - 1
          IF ( I .GT. NASS1 ) THEN
           IF (KEEP(219).NE.0) THEN
            IF (KEEP(50).EQ.2) THEN
              AINPUT=AII
              DO JJ=II,J2
               J=INTARR(JJ)
               IF (J.LE.NASS1) THEN
                A(APOSMAX+J-1)=CMPLX(
     *              MAX(REAL(A(APOSMAX+J-1)),
     *                  ABS(DBLARR(AINPUT)))
     *                 )
               ENDIF
               AINPUT=AINPUT+1
              ENDDO
            ELSE
              AII = AII + J2 - II + 1
              CYCLE
            ENDIF
           ELSE
             AII = AII + J2 - II + 1
             CYCLE
           ENDIF
          ENDIF ! I .GT. NASS1
          IF (KEEP(219).NE.0) THEN
            MAXARR = ZERO
          ENDIF
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF ( J .LE. NASS1) THEN
              IF (I.LT.J) THEN
                  APOS2 = ICT12 + J*LDAFS
              ELSE
                  APOS2 = ICT21 + J
              ENDIF
              A(APOS2) = A(APOS2) + DBLARR(AII)
           ELSE IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
              MAXARR = MAX(MAXARR,ABS(DBLARR(AII)))
           ENDIF
           AII = AII + 1
          END DO
          IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
            A(APOSMAX+I-1) = CMPLX(
     *           MAX( MAXARR, REAL(A(APOSMAX+I-1)))
     *         )
          ENDIF
         END IF
        END DO
      END DO
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6 + XSIZE
      DO ISLAVE = 1, NSLAVES
              CALL CMUMPS_49( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB,
     &                NSLAVES, 
     &                NBLIG, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
        IERR = -1
        DO WHILE (IERR .EQ.-1)
         IF ( KEEP(50) .eq. 0 ) THEN
           NBCOL =  NFRONT
           CALL CMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      IZERO, IDUMMY,
     *      IW(PDEST), NFRONT, COMM, IERR)
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL CMUMPS_68( INODE,
     *      NBPROCFILS(STEP(INODE)),
     *      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     *      NSLAVES-ISLAVE, IW( PTLUST_S(STEP(INODE))+6+XSIZE+ISLAVE),
     *      IW(PDEST), NFRONT, COMM, IERR)
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *     MPI_ANY_SOURCE, MPI_ANY_TAG,
     *     STATUS, BUFR, LBUFR,
     *     LBUFR_BYTES,
     *     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *     LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *     PTLUST_S, PTRFAC,
     *     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *     IERROR, COMM,
     *     NBPROCFILS,
     *     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *     root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *     INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *     NELT+1, NELT, FRT_PTR, FRT_ELT,
     *     ISTEP_TO_INIV2, TAB_POS_IN_PERE)
           IF ( IFLAG .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST_S(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        END DO
        IF (IERR .EQ. -2) GOTO 300
        IF (IERR .EQ. -3) GOTO 305
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      END DO
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1 + XSIZE)
        LSTK    = IW(ISTCHK + XSIZE)
        NPIVS   = IW(ISTCHK + 3 + XSIZE)
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5 + XSIZE)
        HS      = 6 + NSLSON + XSIZE
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LE.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2 + XSIZE )
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6 + XSIZE
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
           NFS4FATHER = NCBSON
           DO I=0,NCBSON-1
              IF(IW(PTRCOL+I) .GT. NASS1) THEN
                 NFS4FATHER = I
                 EXIT
              ENDIF
           ENDDO
           NFS4FATHER=NFS4FATHER + NELIM
        ELSE
          NFS4FATHER = 0
        ENDIF
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = CMUMPS_275(STEP(ISON),
     *                 PROCNODE_STEPS, SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL CMUMPS_211( COMM_LOAD, ASS_IRECV, 
     *      BUFR, LBUFR, LBUFR_BYTES,
     *      INODE, ISON, NSLAVES, 
     *      IW( PTLUST_S(STEP(INODE)) + 6 +XSIZE),
     *      NFRONT, NASS1, NFS4FATHER,NCBSON, IW( PTRCOL ),
     *      PROCNODE_STEPS,
     *      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *      LRLUS, N, IW, LIW, A, LA, NIRBDU,
     *      PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     *      PIMASTER, PAMASTER, NSTK_S, COMP,
     *      IFLAG, IERROR, MYID, COMM, NBPROCFILS, 
     *      IPOOL, LPOOL, LEAF, ICNTL, KEEP,KEEP8, root,
     *      OPASSW, OPELIW,
     *      ITLOC, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     *      ND, FRERE_STEPS, NELT+1, NELT, 
     *      FRT_PTR, FRT_ELT, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
           IF ( IFLAG .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            CALL  CMUMPS_71( 
     *           INODE, NFRONT,NASS1,NFS4FATHER,
     *           ISON, MYID,
     *      NSLAVES, IW( PTLUST_S(STEP(INODE)) + 6 +XSIZE),
     *      IW(PTRCOL), NCBSON,
     *      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     *      KEEP,KEEP8, STEP, N, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *        PTLUST_S, PTRFAC,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE)
              IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ELSE
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
               CALL CMUMPS_49( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL CMUMPS_210( COMM_LOAD, ASS_IRECV, 
     *        BUFR, LBUFR, LBUFR_BYTES,
     *        INODE, ISON, NSLAVES, 
     *        IW( PTLUST_S(STEP(INODE))+6+XSIZE),
     *        NFRONT, NASS1,NFS4FATHER,
     *        TROW_SIZE, IW( INDX ),
     *        PROCNODE_STEPS,
     *        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     *        LRLUS, N, IW,
     *        LIW, A, LA, NIRBDU,
     *        PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     *        PIMASTER, PAMASTER, NSTK_S, COMP,
     *        IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     *        ICNTL, KEEP,KEEP8, root,
     *        OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     * 
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *        )
              IF ( IFLAG .LT. 0 ) GOTO 500
              EXIT
            ENDIF
          END DO
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + XSIZE
            CALL  CMUMPS_71( 
     *           INODE, NFRONT, NASS1, NFS4FATHER,
     *           ISON, MYID,
     *      NSLAVES, IW(PTLUST_S(STEP(INODE))+6+XSIZE),
     *      IW(PTRCOL), NCBSON,
     *      COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     *      KEEP,KEEP8, STEP, N, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_329( COMM_LOAD, ASS_IRECV,
     *        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *        MPI_ANY_SOURCE, MPI_ANY_TAG,
     *        STATUS, BUFR, LBUFR,
     *        LBUFR_BYTES,
     *        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *        LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *        PTLUST_S, PTRFAC,
     *        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *        IERROR, COMM,
     *        NBPROCFILS,
     *        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *        root,OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE_STEPS,
     *        NELT+1, NELT, FRT_PTR, FRT_ELT,
     *        ISTEP_TO_INIV2, TAB_POS_IN_PERE )
             IF ( IFLAG .LT. 0 ) GOTO 500
            ENDIF
          END DO
          IF (IERR .EQ. -2) GOTO 290
          IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ENDIF
       ISON = FRERE_STEPS(STEP(ISON))
      END DO
      GOTO 500
  250 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG   = -13
      IERROR  = NUMSTK + 1
      GOTO 490
  265 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DYNAMIC ALLOCATION during assembly'
      ENDIF
      IFLAG  = -13
      IERROR = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      IFLAG = -8
      IERROR = LREQ
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *  ' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_37'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *  ' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_37'
      ENDIF
      IFLAG = -9
      IERROR = LAELL - LRLUS
      GOTO 490
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE, SEND BUFFER TOO SMALL (1) DURING CMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NCBSON + 6+NSLSON+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE, RECV BUFFER TOO SMALL (1) DURING CMUMPS_37'
      ENDIF
      IFLAG = -20
      LREQ = NCBSON + 6+NSLSON+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE, SENDBUFFER TOO SMALL (2) DURING CMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     *' FAILURE, RECVBUFFER TOO SMALL (2) DURING CMUMPS_37'
      ENDIF
      IFLAG = -17
      LREQ = NBLIG + NBCOL + 4+XSIZE
      IERROR =  LREQ  * KEEP( 34 ) 
      GOTO 490
  490 CALL CMUMPS_44( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_37
      SUBROUTINE CMUMPS_123(
     *    NELT, FRT_PTR, FRT_ELT,
     *    N, INODE, IW, LIW, A, LA, 
     *    NBROWS, NBCOLS,
     *    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     *    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     *    ICNTL, KEEP, KEEP8, MYID)
      IMPLICIT NONE
      INTEGER NELT, N,LIW,LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER*8 KEEP8(150)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER IW(LIW), ITLOC(N), STEP(N),
     *        PTRIST(KEEP(28)), PTRAST(KEEP(28)),
     *        FILS(N), PTRARW(NELT+1), 
     *        PTRAIW(NELT+1)
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      COMPLEX A(LA),
     *        DBLARR(MAX(1,KEEP(13)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER IOLDPS, POSELT, NBCOLF, NBROWF, NSLAVES, HF,
     *        K1,K2,K,I,J,POSEL1,APOS,JPOS,NASS,JJ,
     *        IN,AINPUT,JK,J1,J2,IJROW,ICT12, ILOC, 
     *        ELBEG, NUMELT, ELTI, SIZE_ELTI, IPOS, 
     *        IPOS1, IPOS2, AII, APOS2, II, IELL
      REAL  ZERO
      PARAMETER (ZERO=0.0E0)
      INCLUDE 'mumps_headers.h'
      INTRINSIC REAL
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+XSIZE)
      NBROWF  = IW(IOLDPS+2+XSIZE)
      NASS    = IW(IOLDPS+1+XSIZE)
       IF ( NBROWS .GT. NBROWF ) THEN
          WRITE(*,*) ' ERR: ERROR : NBROWS > NBROWF'
          WRITE(*,*) ' ERR: INODE =', INODE
          WRITE(*,*) ' ERR: NBROW=',NBROWS,'NBROWF=',NBROWF
          CALL CMUMPS_ABORT()
       END IF
      NSLAVES = IW(IOLDPS+5+XSIZE)
      HF      = 6 + NSLAVES+XSIZE
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1+XSIZE) = NASS
          A(POSELT:POSELT+NBROWF*NBCOLF-1) = CMPLX(ZERO)
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          END DO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
           JPOS     = JPOS + 1
          END DO
          ELBEG  = FRT_PTR(INODE)
          NUMELT = FRT_PTR(INODE+1) - ELBEG
          DO IELL=ELBEG,ELBEG+NUMELT-1
           ELTI = FRT_ELT(IELL)
           J1= PTRAIW(ELTI)
           J2= PTRAIW(ELTI+1)-1
           AII = PTRARW(ELTI)
           SIZE_ELTI = J2 - J1 + 1
           DO II=J1,J2
            I = ITLOC(INTARR(II))
            IF (KEEP(50).EQ.0) THEN
             IF (I.LE.0) CYCLE
             AINPUT    = AII + II - J1
             IPOS = MOD(I,NBCOLF)
             ICT12 = POSELT + (IPOS-1) * NBCOLF
             DO JJ = J1, J2
              JPOS = ITLOC(INTARR(JJ))
              IF (JPOS.LE.0) THEN 
                   JPOS = -JPOS
              ELSE
                   JPOS = JPOS/NBCOLF
              END IF
              APOS2    = ICT12 + JPOS - 1
              A(APOS2) = A(APOS2) +  DBLARR(AINPUT)
              AINPUT   = AINPUT + SIZE_ELTI
             END DO
            ELSE
              IF ( I .EQ. 0 ) THEN 
               AII = AII + J2 - II + 1
               CYCLE
              ENDIF
              IF ( I .LE. 0 ) THEN 
               IPOS1 = -I
               IPOS2 = 0
              ELSE 
               IPOS1 = I/NBCOLF
               IPOS2 = MOD(I,NBCOLF)
              END IF
              ICT12 =  POSELT + (IPOS2-1)*NBCOLF
              DO JJ=II,J2
               AII = AII + 1
               J = ITLOC(INTARR(JJ))
               IF ( J .EQ. 0 ) CYCLE
               IF ( IPOS2.EQ.0 .AND. J.LE.0) CYCLE
               IF ( J .LE. 0 ) THEN
                JPOS = -J
               ELSE
                JPOS = J/NBCOLF
               END IF
               IF ( (IPOS1.GE.JPOS) .AND. (IPOS2.GT.0) ) THEN
                 APOS2 = ICT12  + JPOS - 1
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
               IF ( (IPOS1.LT.JPOS) .AND. (J.GT.0) ) THEN
                 IPOS = MOD(J,NBCOLF)
                 JPOS = IPOS1
                 APOS2 = POSELT + (IPOS-1)*NBCOLF + JPOS - 1
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
              END DO
            END IF
           END DO
          END DO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END IF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          END DO
      END IF
 500  CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_123
      SUBROUTINE CMUMPS_126(
     *            N, NELT, NA_ELT,
     *            COMM, MYID, SLAVEF,
     *            IELPTR_LOC, RELPTR_LOC,
     *            ELTVAR_LOC, ELTVAL_LOC,
     *            KEEP,KEEP8, MAXELT_SIZE,
     *            FRTPTR, FRTELT, A, LA, FILS,
     *            id, root )
      USE CMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, NELT, NA_ELT
      INTEGER COMM, MYID, SLAVEF, MAXELT_SIZE, MSGLEN, LA
      INTEGER FRTPTR( N+1 )
      INTEGER FRTELT( NELT ), FILS ( N )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IELPTR_LOC( NELT + 1 ), RELPTR_LOC( NELT + 1 )
      INTEGER ELTVAR_LOC( MAX(1,KEEP(14)) )
      COMPLEX ELTVAL_LOC( MAX(1,KEEP(13)) )
      COMPLEX A( LA )
      TYPE(CMUMPS_STRUC)     :: id
      TYPE(CMUMPS_ROOT_STRUC) :: root
      INTEGER NUMROC
      EXTERNAL NUMROC
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER STATUS( MPI_STATUS_SIZE ), IERR_MPI
      INTEGER MSGTAG
      INTEGER allocok
      INTEGER I, DEST, MAXELT_REAL_SIZE, MPG, IEL, SIZEI, SIZER
      INTEGER NBRECORDS, NBUF
      INTEGER RECV_IELTPTR, RECV_RELTPTR
      INTEGER IELTPTR, RELTPTR, INODE
      LOGICAL FINI, PROKG, I_AM_SLAVE
      INTEGER PTR_ROOT, LOCAL_M, LOCAL_N, LP, IBEG, IGLOB, JGLOB
      INTEGER ARROW_ROOT
      INTEGER IELT, J, K, NB_REC, IREC
      INTEGER ILOCROOT, JLOCROOT, IPOSROOT, JPOSROOT, IPTR
      INTEGER JCOL_GRID, IROW_GRID
      INTEGER IVALPTR
      INTEGER NBELROOT
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      COMPLEX  VAL
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, DIMENSION( :, : ), ALLOCATABLE :: BUFI
      COMPLEX, DIMENSION( :, : ), ALLOCATABLE :: BUFR
      COMPLEX, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_R
      INTEGER, DIMENSION( : ), ALLOCATABLE :: TEMP_ELT_I
      INTEGER, DIMENSION( : ), ALLOCATABLE :: ELROOTPOS
      INTEGER, DIMENSION( : ), ALLOCATABLE, TARGET :: RG2LALLOC
      INTEGER, DIMENSION( : ), POINTER     :: RG2L
      MPG = id%ICNTL(3)
      LP  = id%ICNTL(1)
      I_AM_SLAVE = ( KEEP(46) .eq. 1 .or. MYID .ne.MASTER )
      PROKG = ( MPG > 0 .and. MYID .eq. MASTER )
      KEEP(49) = 0
      ARROW_ROOT = 0
      IF ( MYID .eq. MASTER ) THEN
        IF ( KEEP(46) .eq. 0 ) THEN
          NBUF = SLAVEF
        ELSE
          NBUF = SLAVEF - 1
        END IF
        NBRECORDS = MIN(KEEP(39),NA_ELT)
        IF ( KEEP(50) .eq. 0 ) THEN
          MAXELT_REAL_SIZE = MAXELT_SIZE * MAXELT_SIZE
        ELSE
          MAXELT_REAL_SIZE = MAXELT_SIZE * (MAXELT_SIZE+1)/2
        END IF
        IF ( MAXELT_REAL_SIZE .GT. KEEP(39) ) THEN
          NBRECORDS = MAXELT_REAL_SIZE
          IF ( MPG .GT. 0 ) THEN
            WRITE(MPG,*)
     * ' ** Warning : For element distrib NBRECORDS set to ',
     * MAXELT_REAL_SIZE,' because one element is large'
          END IF
        END IF
        ALLOCATE( BUFI( 2*NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = 2*NBRECORDS + 1
          GOTO 100
        END IF
        ALLOCATE( BUFR( NBRECORDS+1, NBUF ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = NBRECORDS + 1
          GOTO 100
        END IF
        IF ( KEEP(52) .ne. 0 ) THEN
          ALLOCATE( TEMP_ELT_R( MAXELT_REAL_SIZE ), stat =allocok )
          IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_REAL_SIZE
            GOTO 100
          END IF
        END IF
        ALLOCATE( TEMP_ELT_I( MAXELT_SIZE ), stat=allocok )
        IF ( allocok .gt. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = MAXELT_SIZE
            GOTO 100
        END IF
        IF ( KEEP(38) .ne. 0 ) THEN
          NBELROOT = FRTPTR(KEEP(38)+1)-FRTPTR(KEEP(38))
          ALLOCATE( ELROOTPOS( MAX(NBELROOT,1) ),
     *              stat = allocok )
          IF ( allocok .gt. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = NBELROOT
              GOTO 100
          END IF
          IF (KEEP(46) .eq. 0 ) THEN
           ALLOCATE( RG2LALLOC( N ), stat = allocok )
           IF ( allocok .gt. 0 ) THEN
               id%INFO(1) = -13
               id%INFO(2) = N
               GOTO 100
           END IF
           INODE = KEEP(38)
           I     = 1
           DO WHILE ( INODE .GT. 0 )
             RG2LALLOC( INODE ) = I
             INODE = FILS( INODE )
             I = I + 1
           END DO
           RG2L => RG2LALLOC
          ELSE 
           RG2L => root%RG2L_ROW
          END IF
        END IF
        DO I = 1, NBUF
          BUFI( 1, I ) = 0
          BUFR( 1, I ) = CMPLX(0)
        END DO
      END IF
 100  CONTINUE
      CALL CMUMPS_276( id%ICNTL, id%INFO, COMM, MYID )
      IF ( id%INFO(1) .LT. 0 ) RETURN
      CALL MPI_BCAST( NBRECORDS, 1, MPI_INTEGER, MASTER,
     *                COMM, IERR_MPI )
      RECV_IELTPTR = 1
      RECV_RELTPTR = 1
      IF ( MYID .eq. MASTER ) THEN
        NBELROOT = 0
        RELTPTR = 1
        RELPTR_LOC(1) = 1
        DO IEL = 1, NELT
          IELTPTR = id%ELTPTR( IEL )
          SIZEI   = id%ELTPTR( IEL + 1 ) - IELTPTR
          IF ( KEEP( 50 ) .eq. 0 ) THEN
            SIZER = SIZEI * SIZEI
          ELSE
            SIZER = SIZEI * ( SIZEI + 1 ) / 2
          END IF
          DEST = id%ELTPROC( IEL )
          IF ( DEST .eq. -2 ) THEN
            NBELROOT = NBELROOT + 1
            FRTELT( FRTPTR(KEEP(38)) + NBELROOT - 1 ) = IEL
            ELROOTPOS( NBELROOT ) = RELTPTR
            GOTO 200
          END IF
          IF ( DEST .ge. 0 .and. KEEP(46) .eq. 0 ) DEST = DEST + 1
          IF ( KEEP(52) .ne. 0 ) THEN
            CALL CMUMPS_288( N, SIZEI, SIZER,
     *               id%ELTVAR( IELTPTR ), id%A_ELT( RELTPTR ),
     *               TEMP_ELT_R(1), MAXELT_REAL_SIZE,
     *               id%ROWSCA, id%COLSCA, KEEP(50) )
          END IF
          IF ( DEST .eq. 0 .or. ( DEST .eq. -1 .and. KEEP(46) .ne. 0 ) )
     *      THEN
            ELTVAR_LOC( RECV_IELTPTR: RECV_IELTPTR + SIZEI - 1 )
     *      = id%ELTVAR( IELTPTR: IELTPTR + SIZEI - 1 )
            RECV_IELTPTR = RECV_IELTPTR + SIZEI
            IF ( KEEP(52) .ne. 0 ) THEN
              ELTVAL_LOC( RECV_RELTPTR: RECV_RELTPTR + SIZER - 1)
     *        = TEMP_ELT_R( 1: SIZER )
              RECV_RELTPTR = RECV_RELTPTR + SIZER
            END IF
          END IF
          IF ( DEST .NE. 0 .AND. DEST. NE. -3 ) THEN
            IF ( KEEP(52) .eq. 0 ) THEN
              CALL CMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           id%A_ELT (RELTPTR),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
            ELSE
              CALL CMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           TEMP_ELT_R( 1 ),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
            END IF
          END IF
 200      CONTINUE
          RELTPTR = RELTPTR + SIZER
          IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
            RELPTR_LOC( IEL + 1 ) = RELTPTR
          ELSE
            RELPTR_LOC( IEL + 1 ) = RECV_RELTPTR
          ENDIF
        END DO
        IF ( KEEP(46) .eq. 0 .OR. KEEP(52) .eq. 0 ) THEN
          KEEP(13) = RELTPTR - 1
        ELSE
          KEEP(13) = RECV_RELTPTR - 1
        ENDIF
        IF ( RELTPTR - 1 .ne. id%NA_ELT ) THEN
          WRITE(*,*) ' ** ERROR ELT DIST: RELPTR - 1 / id%NA_ELT=',
     *               RELTPTR - 1,id%NA_ELT
          CALL CMUMPS_ABORT()
        END IF
        DEST = -2
        IELTPTR = 1
        RELTPTR = 1
        SIZEI   = 1
        SIZER   = 1
        CALL CMUMPS_127(
     *           id%ELTVAR(IELTPTR),
     *           id%A_ELT (RELTPTR),
     *           SIZEI, SIZER,
     *
     *           DEST, NBUF, NBRECORDS,
     *           BUFI, BUFR, COMM )
        IF ( KEEP(52) .NE. 0 ) DEALLOCATE( TEMP_ELT_R )
      ELSE
        FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     *     .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        DO WHILE ( .not. FINI )
          CALL MPI_PROBE( MASTER, MPI_ANY_TAG,
     *                    COMM, STATUS, IERR_MPI )
          MSGTAG = STATUS( MPI_TAG    )
          SELECT CASE ( MSGTAG )
             CASE( ELT_INT )
               CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,
     *                             MSGLEN, IERR_MPI )
               CALL MPI_RECV( ELTVAR_LOC( RECV_IELTPTR ), MSGLEN,
     *            MPI_INTEGER, MASTER, ELT_INT,
     *            COMM, STATUS, IERR_MPI )
               RECV_IELTPTR = RECV_IELTPTR + MSGLEN
             CASE( ELT_REAL )
                CALL MPI_GET_COUNT( STATUS, MPI_COMPLEX,
     *                              MSGLEN, IERR_MPI )
                CALL MPI_RECV( ELTVAL_LOC( RECV_RELTPTR ), MSGLEN,
     *            MPI_COMPLEX, MASTER, ELT_REAL,
     *            COMM, STATUS, IERR_MPI )
                RECV_RELTPTR = RECV_RELTPTR + MSGLEN
          END SELECT
          FINI = ( RECV_IELTPTR .eq. IELPTR_LOC( NELT+1 )
     *       .and. RECV_RELTPTR .eq. RELPTR_LOC( NELT+1 ) )
        END DO
      END IF
      IF ( KEEP(38) .NE. 0 ) THEN
         IF ( I_AM_SLAVE .and. root%yes ) THEN
          IF (KEEP(60)==0) THEN
           LOCAL_M = NUMROC( root%ROOT_SIZE, root%MBLOCK,
     *             root%MYROW, 0, root%NPROW )
           LOCAL_M = MAX( 1, LOCAL_M )
           LOCAL_N = NUMROC( root%ROOT_SIZE, root%NBLOCK,
     *               root%MYCOL, 0, root%NPCOL )
           PTR_ROOT = LA - LOCAL_M * LOCAL_N + 1
           IF ( PTR_ROOT .LE. LA ) THEN
             A( PTR_ROOT:LA ) = CMPLX(ZERO)
           END IF
          ELSE
           DO I = 1, root%SCHUR_NLOC
            root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     *      (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC)=CMPLX(ZERO)
           ENDDO
          ENDIF
         END IF
        IF ( MYID .NE. MASTER ) THEN
          ALLOCATE( BUFI( NBRECORDS * 2 + 1, 1 ), stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS * 2 + 1
            GOTO 250
          END IF
          ALLOCATE( BUFR( NBRECORDS, 1 )        , stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            id%INFO(1) = -13
            id%INFO(2) = NBRECORDS
          END IF
        END IF
 250    CONTINUE
        CALL CMUMPS_276( id%ICNTL, id%INFO, COMM, MYID )
        IF ( id%INFO(1) .LT. 0 ) RETURN
        IF ( MYID .eq. MASTER ) THEN
        DO IPTR = FRTPTR(KEEP(38)), FRTPTR(KEEP(38)+1) - 1
          IELT = FRTELT( IPTR )
          SIZEI = id%ELTPTR( IELT + 1 ) - id%ELTPTR( IELT )
          DO I = 1, SIZEI
            TEMP_ELT_I( I ) = RG2L
     *              ( id%ELTVAR( id%ELTPTR(IELT) + I - 1 ) )
          END DO
          IVALPTR = ELROOTPOS( IPTR - FRTPTR(KEEP(38)) + 1 ) - 1
          K = 1
          DO J = 1, SIZEI
            JGLOB = id%ELTVAR( id%ELTPTR( IELT ) + J - 1 )
            IF ( KEEP(50).eq. 0 ) THEN
              IBEG = 1
            ELSE
              IBEG = J
            END IF
            DO I = IBEG, SIZEI
              IGLOB = id%ELTVAR( id%ELTPTR( IELT ) + I - 1 )
              IF ( KEEP(52) .eq. 0 ) THEN
                VAL = id%A_ELT( IVALPTR + K )
              ELSE
                VAL = id%A_ELT( IVALPTR + K ) *
     *                id%ROWSCA( IGLOB ) * id%COLSCA( JGLOB )
              END IF
              IF ( KEEP(50).eq.0 ) THEN
                IPOSROOT = TEMP_ELT_I( I )
                JPOSROOT = TEMP_ELT_I( J )
              ELSE
                IF ( TEMP_ELT_I(I) .GT. TEMP_ELT_I(J) ) THEN
                  IPOSROOT = TEMP_ELT_I(I)
                  JPOSROOT = TEMP_ELT_I(J)
                ELSE
                  IPOSROOT = TEMP_ELT_I(J)
                  JPOSROOT = TEMP_ELT_I(I)
                END IF
              END IF
              IROW_GRID = MOD( ( IPOSROOT - 1 )/root%MBLOCK,
     *                           root%NPROW )
              JCOL_GRID = MOD( ( JPOSROOT - 1 )/root%NBLOCK,
     *                           root%NPCOL )
              IF ( KEEP(46) .eq. 0 ) THEN
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
              ELSE
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID
              END IF
              IF ( DEST .eq. MASTER ) THEN
                ARROW_ROOT = ARROW_ROOT + 1
                ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
                JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60)==0) THEN
                A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *        + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *        * LOCAL_M + ILOCROOT - 1 ) + VAL
              ELSE
                root%SCHUR_POINTER( (JLOCROOT-1) * root%SCHUR_LLD
     *        + ILOCROOT ) = root%SCHUR_POINTER( (JLOCROOT-1) *
     *        root%SCHUR_LLD + ILOCROOT ) + VAL
              ENDIF
              ELSE
                CALL CMUMPS_34(
     *          IPOSROOT, JPOSROOT, VAL, DEST, BUFI, BUFR, NBRECORDS,
     *          NBUF, LP, COMM, KEEP(46) )
              END IF
              K = K + 1
            END DO
          END DO
        END DO
        CALL CMUMPS_18(
     *          BUFI, BUFR, NBRECORDS,
     *          NBUF, LP, COMM, KEEP(46) )
        ELSE
          FINI = .FALSE.
          DO WHILE ( .not. FINI )
            CALL MPI_RECV( BUFI(1,1), 2*NBRECORDS+1,
     *                MPI_INTEGER, MASTER,
     *                ARROWHEAD,
     *                COMM, STATUS, IERR_MPI )
            NB_REC = BUFI(1,1)
            IF (NB_REC.LE.0) THEN
              FINI = .TRUE.
              NB_REC = -NB_REC
            ENDIF
            IF (NB_REC.EQ.0) EXIT
            CALL MPI_RECV( BUFR(1,1), NBRECORDS, MPI_COMPLEX,
     *                     MASTER, ARROWHEAD,
     *                     COMM, STATUS, IERR_MPI )
            ARROW_ROOT = ARROW_ROOT + NB_REC
            DO IREC = 1, NB_REC
              IPOSROOT = BUFI( IREC * 2, 1 )
              JPOSROOT = BUFI( IREC * 2 + 1, 1 )
              VAL      = BUFR( IREC, 1 )
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     *                 ( root%MBLOCK * root%NPROW ) )
     *               + MOD( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     *                 ( root%NBLOCK * root%NPCOL ) )
     *               + MOD( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60).eq.0) THEN
              A( PTR_ROOT + ( JLOCROOT - 1 ) * LOCAL_M
     *      + ILOCROOT - 1 ) =  A( PTR_ROOT + (JLOCROOT - 1)
     *      * LOCAL_M + ILOCROOT - 1 ) + VAL
              ELSE
              root%SCHUR_POINTER(( JLOCROOT - 1 ) * root%SCHUR_LLD
     *      + ILOCROOT ) = root%SCHUR_POINTER( (JLOCROOT - 1)
     *      * root%SCHUR_LLD + ILOCROOT) + VAL
              ENDIF
            END DO
          END DO
          DEALLOCATE( BUFI )
          DEALLOCATE( BUFR )
        END IF
      END IF
      IF ( MYID .eq. MASTER ) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
        IF (KEEP(38).ne.0) THEN 
          DEALLOCATE(ELROOTPOS)
          IF (KEEP(46) .eq. 0 ) THEN
             DEALLOCATE(RG2LALLOC)
          ENDIF
        ENDIF
        DEALLOCATE( TEMP_ELT_I )
      END IF
      KEEP(49) = ARROW_ROOT
      RETURN
      END SUBROUTINE CMUMPS_126
      SUBROUTINE CMUMPS_127(
     *       ELNODES, ELVAL, SIZEI, SIZER,
     *       DEST, NBUF, NBRECORDS, BUFI, BUFR, COMM )
      IMPLICIT NONE
      INTEGER SIZEI, SIZER, DEST, NBUF, NBRECORDS, COMM
      INTEGER ELNODES( SIZEI ), BUFI( 2*NBRECORDS + 1, NBUF )
      COMPLEX ELVAL( SIZER ), BUFR( NBRECORDS + 1, NBUF )
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER I, IBEG, IEND, IERR_MPI, NBRECR
      INTEGER NBRECI
      IF ( DEST .lt. 0 ) THEN
        IBEG = 1
        IEND = NBUF
      ELSE
        IBEG = DEST
        IEND = DEST
      END IF
      DO I = IBEG, IEND
        NBRECI = BUFI(1,I)
        IF ( NBRECI .ne.0  .and.
     *       ( DEST.eq.-2 .or.
     *         NBRECI + SIZEI .GT. 2*NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFI(2, I), NBRECI, MPI_INTEGER,
     *                    I, ELT_INT, COMM, IERR_MPI )
           BUFI(1,I) = 0
           NBRECI    = 0
        END IF
        NBRECR = INT(BUFR(1,I)+0.5E0)
        IF ( NBRECR .ne.0  .and.
     *       ( DEST.eq.-2 .or.
     *         NBRECR + SIZER .GT. NBRECORDS ) ) THEN
           CALL MPI_SEND( BUFR(2, I), NBRECR, MPI_COMPLEX,
     *                    I, ELT_REAL, COMM, IERR_MPI )
           BUFR(1,I) = 0.0E0
           NBRECR    = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          BUFI( 2 + NBRECI : 2 + NBRECI + SIZEI - 1, I ) =
     *    ELNODES( 1: SIZEI )
          BUFR( 2 + NBRECR : 2 + NBRECR + SIZER - 1, I ) =
     *    ELVAL( 1: SIZER )
          BUFI(1,I) = NBRECI + SIZEI
          BUFR(1,I) = CMPLX( NBRECR + SIZER )
        END IF
      END DO
      RETURN
      END SUBROUTINE CMUMPS_127
      SUBROUTINE CMUMPS_213( ELTPTR, NELT, MAXELT_SIZE )
      INTEGER NELT, MAXELT_SIZE
      INTEGER ELTPTR( NELT + 1 )
      INTEGER I, S
      MAXELT_SIZE = 0
      DO I = 1, NELT
        S = ELTPTR( I + 1 ) - ELTPTR( I )
        MAXELT_SIZE = MAX( S, MAXELT_SIZE )
      END DO
      RETURN
      END SUBROUTINE CMUMPS_213
      SUBROUTINE CMUMPS_288( N, SIZEI, SIZER,
     *               ELTVAR, ELTVAL,
     *               SELTVAL, LSELTVAL,
     *               ROWSCA, COLSCA, K50 )
      INTEGER N, SIZEI, SIZER, LSELTVAL, K50
      INTEGER ELTVAR( SIZEI )
      COMPLEX ELTVAL( SIZER )
      COMPLEX SELTVAL( LSELTVAL )
      COMPLEX ROWSCA( N ), COLSCA( N )
      INTEGER I, J, K
      K = 1
      IF ( K50 .eq. 0 ) THEN
        DO J = 1, SIZEI
          DO I = 1, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     *                   ROWSCA(ELTVAR(I)) *
     *                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      ELSE
        DO J = 1, SIZEI
          DO I = J, SIZEI
            SELTVAL(K) = ELTVAL(K) *
     *                   ROWSCA(ELTVAR(I)) *
     *                   COLSCA(ELTVAR(J))
            K = K + 1
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE CMUMPS_288
      SUBROUTINE CMUMPS_257( N, NELT, ELTPTR, ELTVAR, A_ELT,
     *                         X, Y, K50, MTYPE )
      IMPLICIT NONE
      INTEGER N, NELT, K50, MTYPE
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( * )
      COMPLEX A_ELT( * ), X( N ), Y( N )
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      COMPLEX TEMP
      REAL ZERO 
      PARAMETER( ZERO = 0.0E0 )
      Y = CMPLX(ZERO)
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
     *          Y( ELTVAR( IELPTR + I ) ) +
     *             A_ELT( K ) * TEMP
                K = K + 1
              END DO
            END DO
          ELSE
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP + 
     *          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
            END DO
          END IF
        ELSE
          DO J = 1, SIZEI
            Y( ELTVAR( IELPTR + J ) ) =
     *      Y( ELTVAR( IELPTR + J ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            K = K + 1
            DO I = J+1, SIZEI
              Y( ELTVAR( IELPTR + I ) ) =
     *        Y( ELTVAR( IELPTR + I ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     *        Y( ELTVAR( IELPTR + J ) ) +
     *           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE CMUMPS_257
      SUBROUTINE CMUMPS_F77( JOB, SYM, PAR, COMM_F77, N, ICNTL, CNTL,
     *                      NZ, IRN, IRNhere, JCN, JCNhere, A, Ahere,
     *                      NZ_loc, IRN_loc, IRN_lochere,
     *                      JCN_loc, JCN_lochere,
     *                      A_loc, A_lochere,
     *                      NELT, ELTPTR, ELTPTRhere,  ELTVAR,
     *                      ELTVARhere, A_ELT, A_ELThere,
     *                      PERM_IN, PERM_INhere,
     *                      RHS, RHShere, REDRHS, REDRHShere,
     *                      INFO, RINFO, INFOG, RINFOG,
     *                      DEFICIENCY, SIZE_SCHUR, LISTVAR_SCHUR,
     *                      LISTVAR_SCHURhere, SCHUR, SCHURhere,
     *                      COLSCA, COLSCAhere, ROWSCA, ROWSCAhere,
     *                      INSTANCE_NUMBER, NRHS, LRHS, LREDRHS,
     *
     *                      RHS_SPARSE, RHS_SPARSEhere,
     *                      SOL_LOC, SOL_LOChere,
     *                      IRHS_SPARSE, IRHS_SPARSEhere,
     *                      IRHS_PTR, IRHS_PTRhere,
     *                      ISOL_LOC, ISOL_LOChere,
     *                      NZ_RHS, LSOL_LOC
     *                      , 
     * SCHUR_MLOC,
     * SCHUR_NLOC,
     * SCHUR_LLD,
     * MBLOCK,
     * NBLOCK,
     * NPROW,
     * NPCOL,
     *
     * OOC_TMPDIR,
     * OOC_PREFIX,
     * TMPDIRLEN,
     * PREFIXLEN
     *
     * )
      USE CMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER JOB, SYM, PAR, COMM_F77, N, NZ, NZ_LOC, NELT,
     *        DEFICIENCY, SIZE_SCHUR, INSTANCE_NUMBER, NRHS, LRHS,
     *        NZ_RHS, LSOL_LOC, LREDRHS
      INTEGER ICNTL(40), INFO(40), INFOG(40)
      INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
      INTEGER MBLOCK, NBLOCK, NPROW, NPCOL, TMPDIRLEN, PREFIXLEN
      REAL CNTL(15), RINFO(20), RINFOG(20)
      INTEGER, TARGET :: IRN(*), JCN(*), ELTPTR(*), ELTVAR(*)
      INTEGER, TARGET :: PERM_IN(*), IRN_loc(*), JCN_loc(*)
      INTEGER, TARGET :: LISTVAR_SCHUR(*)
      INTEGER, TARGET :: IRHS_PTR(*), IRHS_SPARSE(*), ISOL_LOC(*)
      COMPLEX, TARGET :: A(*), A_ELT(*), A_LOC(*), RHS(*)
      COMPLEX, TARGET :: REDRHS(*)
      COMPLEX, TARGET :: ROWSCA(*), COLSCA(*)
      COMPLEX, TARGET :: SCHUR(*)
      COMPLEX, TARGET :: RHS_SPARSE(*), SOL_LOC(*)
      INTEGER, INTENT(in) :: OOC_TMPDIR(150)
      INTEGER, INTENT(in) :: OOC_PREFIX(150)
      INTEGER IRNhere, JCNhere, Ahere, ELTPTRhere, ELTVARhere,
     *        A_ELThere, PERM_INhere, RHShere, REDRHShere, IRN_lochere,
     *        JCN_lochere, A_lochere, LISTVAR_SCHURhere,
     *        SCHURhere, COLSCAhere, ROWSCAhere, RHS_SPARSEhere,
     *        SOL_LOChere, IRHS_PTRhere, IRHS_SPARSEhere, ISOL_LOChere,
     *        CMUMPS_OOC_PREFIXhere, CMUMPS_OOC_TMDIRhere
      INCLUDE 'mpif.h'
      TYPE CMUMPS_STRUC_PTR
          TYPE (CMUMPS_STRUC), POINTER :: ptr
      END TYPE CMUMPS_STRUC_PTR
      TYPE (CMUMPS_STRUC), POINTER :: mumps_par
      TYPE (CMUMPS_STRUC_PTR), DIMENSION (:), POINTER, SAVE ::
     *  mumps_par_array
      TYPE (CMUMPS_STRUC_PTR), DIMENSION (:), POINTER ::
     *  mumps_par_array_bis
      INTEGER, SAVE :: CMUMPS_STRUC_ARRAY_SIZE = 0
      INTEGER, SAVE :: N_INSTANCES = 0
      INTEGER A_ELT_SIZE, I, Np, IERR
      INTEGER CMUMPS_STRUC_ARRAY_SIZE_INIT
      PARAMETER (CMUMPS_STRUC_ARRAY_SIZE_INIT=10)
      EXTERNAL CMUMPS_AFFECT_MAPPING,
     *         CMUMPS_AFFECT_PIVNUL_LIST,
     *         CMUMPS_AFFECT_SYM_PERM,
     *         CMUMPS_AFFECT_UNS_PERM
      IF (JOB == -1) THEN
        DO I = 1, CMUMPS_STRUC_ARRAY_SIZE
          IF ( .NOT. ASSOCIATED(mumps_par_array(I)%PTR) ) GOTO 10
        END DO
        ALLOCATE( mumps_par_array_bis(CMUMPS_STRUC_ARRAY_SIZE +
     *  CMUMPS_STRUC_ARRAY_SIZE_INIT), stat=ierr)
        IF (ierr /= 0) THEN
          WRITE(*,*) ' ** Allocation Error 1 in CMUMPS_F77.'
          CALL CMUMPS_ABORT()
        END IF
        DO I = 1, CMUMPS_STRUC_ARRAY_SIZE
          mumps_par_array_bis(I)%PTR=>mumps_par_array(I)%PTR
        ENDDO
        IF (ASSOCIATED(mumps_par_array)) DEALLOCATE(mumps_par_array)
        mumps_par_array=>mumps_par_array_bis
        NULLIFY(mumps_par_array_bis)
        DO I = CMUMPS_STRUC_ARRAY_SIZE+1, CMUMPS_STRUC_ARRAY_SIZE +
     *  CMUMPS_STRUC_ARRAY_SIZE_INIT
          NULLIFY(mumps_par_array(I)%PTR)
        ENDDO
        I = CMUMPS_STRUC_ARRAY_SIZE+1
        CMUMPS_STRUC_ARRAY_SIZE = CMUMPS_STRUC_ARRAY_SIZE +
     *  CMUMPS_STRUC_ARRAY_SIZE_INIT
 10     CONTINUE
        INSTANCE_NUMBER = I
        N_INSTANCES = N_INSTANCES+1
        ALLOCATE( mumps_par_array(INSTANCE_NUMBER)%PTR,stat=IERR )
        IF (IERR /= 0) THEN
          WRITE(*,*) '** Allocation Error 2 in CMUMPS_F77.'
          CALL CMUMPS_ABORT()
        ENDIF
        mumps_par_array(INSTANCE_NUMBER)%PTR%KEEP(40) = 0
        mumps_par_array(INSTANCE_NUMBER)%PTR%INSTANCE_NUMBER =
     *  INSTANCE_NUMBER
      END IF
      IF ( INSTANCE_NUMBER .LE. 0 .OR. INSTANCE_NUMBER .GT.
     *     CMUMPS_STRUC_ARRAY_SIZE ) THEN
        WRITE(*,*) ' ** Instance Error 1 in CMUMPS_F77',
     *             INSTANCE_NUMBER
        CALL CMUMPS_ABORT()
      END IF
      IF ( .NOT. ASSOCIATED ( mumps_par_array(INSTANCE_NUMBER)%PTR ) )
     *  THEN
        WRITE(*,*) ' Instance Error 2 in CMUMPS_F77',
     *             INSTANCE_NUMBER
        CALL CMUMPS_ABORT()
      END IF
      mumps_par => mumps_par_array(INSTANCE_NUMBER)%PTR
      mumps_par%SYM = SYM
      mumps_par%PAR = PAR
      mumps_par%JOB = JOB
      mumps_par%N   = N
      mumps_par%NZ  = NZ
      mumps_par%NZ_loc  = NZ_loc
      mumps_par%SIZE_SCHUR  = SIZE_SCHUR
      mumps_par%NELT= NELT
      mumps_par%ICNTL(1:40)=ICNTL(1:40)
      mumps_par%CNTL(1:15)=CNTL(1:15)
      mumps_par%NRHS  = NRHS
      mumps_par%LRHS  = LRHS
      mumps_par%LREDRHS = LREDRHS
      mumps_par%NZ_RHS   = NZ_RHS
      mumps_par%LSOL_LOC = LSOL_LOC
      mumps_par%SCHUR_MLOC   = SCHUR_MLOC
      mumps_par%SCHUR_NLOC   = SCHUR_NLOC
      mumps_par%SCHUR_LLD    = SCHUR_LLD
      mumps_par%MBLOCK = MBLOCK
      mumps_par%NBLOCK = NBLOCK
      mumps_par%NPROW  = NPROW
      mumps_par%NPCOL  = NPCOL
      IF ( COMM_F77 .NE. -987654 ) THEN
        mumps_par%COMM = COMM_F77
      ELSE
        mumps_par%COMM = MPI_COMM_WORLD
      ENDIF
      CALL MPI_BCAST(NRHS,1,MPI_INTEGER,0,mumps_par%COMM,IERR)
      IF ( IRNhere /= 0 ) mumps_par%IRN => IRN(1:NZ)
      IF ( JCNhere /= 0 ) mumps_par%JCN => JCN(1:NZ)
      IF ( Ahere /= 0 )   mumps_par%A   => A(1:NZ)
      IF ( IRN_lochere /= 0 ) mumps_par%IRN_loc => IRN_loc(1:NZ_loc)
      IF ( JCN_lochere /= 0 ) mumps_par%JCN_loc => JCN_loc(1:NZ_loc)
      IF ( A_lochere /= 0 )   mumps_par%A_loc   => A_loc(1:NZ_loc)
      IF ( ELTPTRhere /= 0 ) mumps_par%ELTPTR => ELTPTR(1:NELT+1)
      IF ( ELTVARhere /= 0 ) mumps_par%ELTVAR =>
     *   ELTVAR(1:ELTPTR(NELT+1)-1)
      IF ( A_ELThere /= 0 ) THEN
        A_ELT_SIZE = 0
        DO I = 1, NELT
          Np = ELTPTR(I+1) -ELTPTR(I)
          IF (SYM == 0) THEN
            A_ELT_SIZE = A_ELT_SIZE + Np * Np
          ELSE
            A_ELT_SIZE = A_ELT_SIZE + Np * ( Np + 1 ) / 2
          END IF
        END DO
        mumps_par%A_ELT => A_ELT(1:A_ELT_SIZE)
      END IF
      IF ( PERM_INhere /= 0) mumps_par%PERM_IN => PERM_IN(1:N)
      IF ( LISTVAR_SCHURhere /= 0)
     *   mumps_par%LISTVAR_SCHUR =>LISTVAR_SCHUR(1:SIZE_SCHUR)
      IF ( SCHURhere /= 0 ) THEN
        mumps_par%SCHUR_CINTERFACE=>SCHUR(1:1)
      ENDIF
      IF (NRHS .NE. 1) THEN
        IF ( RHShere /= 0 ) mumps_par%RHS => RHS(1:NRHS*LRHS)
        IF (REDRHShere /= 0)mumps_par%REDRHS=>REDRHS(1:NRHS*LREDRHS)
      ELSE
        IF ( RHShere /= 0 ) mumps_par%RHS => RHS(1:N)
        IF (REDRHShere /= 0)mumps_par%REDRHS=>REDRHS(1:SIZE_SCHUR)
      ENDIF
      IF ( COLSCAhere /= 0) mumps_par%COLSCA => COLSCA(1:N)
      IF ( ROWSCAhere /= 0) mumps_par%ROWSCA => ROWSCA(1:N)
      IF ( RHS_SPARSEhere /=0 ) mumps_par%RHS_SPARSE=>
     *                          RHS_SPARSE(1:NZ_RHS)
      IF ( IRHS_SPARSEhere /=0 ) mumps_par%IRHS_SPARSE=>
     *                          IRHS_SPARSE(1:NZ_RHS)
      IF ( SOL_LOChere /=0 ) mumps_par%SOL_LOC=>
     *                          SOL_LOC(1:LSOL_LOC*NRHS)
      IF ( ISOL_LOChere /=0 ) mumps_par%ISOL_LOC=>
     *                          ISOL_LOC(1:LSOL_LOC)
      IF ( IRHS_PTRhere /=0 ) mumps_par%IRHS_PTR=>
     *                          IRHS_PTR(1:NRHS+1)
      DO I=1,TMPDIRLEN
        mumps_par%OOC_TMPDIR(I:I)=CHAR(OOC_TMPDIR(I))
      ENDDO
      DO I=TMPDIRLEN+1,150
        mumps_par%OOC_TMPDIR(I:I)=' '
      ENDDO
      DO I=1,PREFIXLEN
        mumps_par%OOC_PREFIX(I:I)=CHAR(OOC_PREFIX(I))
      ENDDO
      DO I=PREFIXLEN+1,150
        mumps_par%OOC_PREFIX(I:I)=' '
      ENDDO
      CALL CMUMPS( mumps_par )
      INFO(1:40)=mumps_par%INFO(1:40)
      INFOG(1:40)=mumps_par%INFOG(1:40)
      RINFO(1:20)=mumps_par%RINFO(1:20)
      RINFOG(1:20)=mumps_par%RINFOG(1:20)
      ICNTL(1:40) = mumps_par%ICNTL(1:40)
      CNTL(1:15) = mumps_par%CNTL(1:15)
      SYM = mumps_par%SYM
      PAR = mumps_par%PAR
      JOB = mumps_par%JOB
      N   = mumps_par%N
      NZ  = mumps_par%NZ
      NRHS = mumps_par%NRHS
      LRHS = mumps_par%LRHS
      LREDRHS = mumps_par%LREDRHS
      NZ_loc  = mumps_par%NZ_loc
      NZ_RHS  = mumps_par%NZ_RHS
      LSOL_LOC= mumps_par%LSOL_LOC
      SIZE_SCHUR  = mumps_par%SIZE_SCHUR
      NELT= mumps_par%NELT
      DEFICIENCY = mumps_par%Deficiency
      SCHUR_MLOC   = mumps_par%SCHUR_MLOC
      SCHUR_NLOC   = mumps_par%SCHUR_NLOC
      SCHUR_LLD    = mumps_par%SCHUR_LLD
      MBLOCK       = mumps_par%MBLOCK
      NBLOCK       = mumps_par%NBLOCK
      NPROW        = mumps_par%NPROW
      NPCOL        = mumps_par%NPCOL
      IF ( ASSOCIATED (mumps_par%MAPPING) ) THEN
         CALL CMUMPS_AFFECT_MAPPING(mumps_par%MAPPING(1))
      ELSE
         CALL CMUMPS_NULLIFY_C_MAPPING()
      ENDIF
      IF ( ASSOCIATED (mumps_par%PIVNUL_LIST) ) THEN
         CALL CMUMPS_AFFECT_PIVNUL_LIST(mumps_par%PIVNUL_LIST(1))
      ELSE
         CALL CMUMPS_NULLIFY_C_PIVNUL_LIST()
      ENDIF
      IF ( ASSOCIATED (mumps_par%SYM_PERM) ) THEN
         CALL CMUMPS_AFFECT_SYM_PERM(mumps_par%SYM_PERM(1))
      ELSE
         CALL CMUMPS_NULLIFY_C_SYM_PERM()
      ENDIF
      IF ( ASSOCIATED (mumps_par%UNS_PERM) ) THEN
         CALL CMUMPS_AFFECT_UNS_PERM(mumps_par%UNS_PERM(1))
      ELSE
         CALL CMUMPS_NULLIFY_C_UNS_PERM()
      ENDIF
      IF ( JOB == -2 ) THEN
         IF (ASSOCIATED(mumps_par_array(INSTANCE_NUMBER)%PTR))THEN
           DEALLOCATE(mumps_par_array(INSTANCE_NUMBER)%PTR)
           NULLIFY   (mumps_par_array(INSTANCE_NUMBER)%PTR)
           N_INSTANCES = N_INSTANCES - 1
           IF ( N_INSTANCES == 0 ) THEN
             DEALLOCATE(mumps_par_array)
             CMUMPS_STRUC_ARRAY_SIZE = 0
           END IF
         ELSE
           WRITE(*,*) "** Warning: instance already freed"
           WRITE(*,*) "            this should normally not happen."
         ENDIF
      END IF
      RETURN
      END SUBROUTINE CMUMPS_F77
