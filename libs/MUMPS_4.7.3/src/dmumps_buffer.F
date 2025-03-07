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
        MODULE DMUMPS_BUFFER
        PRIVATE
        PUBLIC :: DMUMPS_61, DMUMPS_528,
     *   DMUMPS_53 ,       DMUMPS_57 ,
     *   DMUMPS_55, DMUMPS_59,
     *   DMUMPS_54,DMUMPS_58,
     *   DMUMPS_66,     DMUMPS_78,
     *   DMUMPS_62,       DMUMPS_68,
     *   DMUMPS_71, DMUMPS_70,
     *   DMUMPS_67,
     *   DMUMPS_65,
     *   DMUMPS_64,
     *   DMUMPS_72,
     *   DMUMPS_75,
     *   DMUMPS_76,
     *   DMUMPS_73, DMUMPS_74,
     *   DMUMPS_63,DMUMPS_77, 
     *   DMUMPS_60,
     *   DMUMPS_524, DMUMPS_469,
     *   DMUMPS_460, DMUMPS_502,
     *   DMUMPS_519
     *   ,DMUMPS_620
     *   ,DMUMPS_617
        INTEGER NEXT, REQ, CONTENT, OVHSIZE
        PARAMETER( NEXT = 0, REQ = 1, CONTENT = 2, OVHSIZE = 2 )
        INTEGER, SAVE :: SIZEofINT, BUF_MYID
        TYPE DMUMPS_BUFFER_TYPE
          INTEGER LBUF, HEAD, TAIL,LBUF_INT, ILASTMSG
          INTEGER, DIMENSION(:),POINTER :: CONTENT
        END TYPE DMUMPS_BUFFER_TYPE
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_CB
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_SMALL
        TYPE ( DMUMPS_BUFFER_TYPE ), SAVE :: BUF_LOAD
        INTEGER, SAVE :: SIZE_RBUF_BYTES
        INTEGER BUF_LMAX_ARRAY
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BUF_MAX_ARRAY
        PUBLIC :: BUF_LMAX_ARRAY, BUF_MAX_ARRAY
      CONTAINS
        SUBROUTINE DMUMPS_528( MYID )
        IMPLICIT NONE
        INTEGER MYID
        BUF_MYID  = MYID
        RETURN
        END SUBROUTINE DMUMPS_528
        SUBROUTINE DMUMPS_61( IntSize )
        IMPLICIT NONE
        INTEGER IntSize
        SIZEofINT = IntSize
        NULLIFY(BUF_CB  %CONTENT)
        NULLIFY(BUF_SMALL%CONTENT)
        NULLIFY(BUF_LOAD%CONTENT)
        BUF_CB%LBUF     = 0
        BUF_CB%LBUF_INT = 0
        BUF_CB%HEAD     = 1
        BUF_CB%TAIL     = 1
        BUF_CB%ILASTMSG = 1
        BUF_SMALL%LBUF     = 0
        BUF_SMALL%LBUF_INT = 0
        BUF_SMALL%HEAD     = 1
        BUF_SMALL%TAIL     = 1
        BUF_SMALL%ILASTMSG = 1
        BUF_LOAD%LBUF     = 0
        BUF_LOAD%LBUF_INT = 0
        BUF_LOAD%HEAD     = 1
        BUF_LOAD%TAIL     = 1
        BUF_LOAD%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_61
        SUBROUTINE DMUMPS_53( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL DMUMPS_2( BUF_CB, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_53
        SUBROUTINE DMUMPS_55( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        CALL DMUMPS_2( BUF_SMALL, SIZE, IERR )
        RETURN
        END SUBROUTINE DMUMPS_55
        SUBROUTINE DMUMPS_54( SIZE, IERR )
        IMPLICIT NONE
        INTEGER SIZE, IERR
        INTEGER stat
        CALL DMUMPS_2( BUF_LOAD, SIZE, IERR )        
        RETURN
        END SUBROUTINE DMUMPS_54
        SUBROUTINE DMUMPS_58( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_LOAD, IERR )
        RETURN
        END SUBROUTINE DMUMPS_58
        SUBROUTINE DMUMPS_620()
        IMPLICIT NONE
        IF (ALLOCATED( BUF_MAX_ARRAY)) DEALLOCATE( BUF_MAX_ARRAY )
        RETURN
        END SUBROUTINE DMUMPS_620
        SUBROUTINE DMUMPS_617(NFS4FATHER,IERR)
        IMPLICIT NONE
        INTEGER IERR, NFS4FATHER
        IERR = 0
        IF (ALLOCATED( BUF_MAX_ARRAY)) THEN
          IF (BUF_LMAX_ARRAY .GE. NFS4FATHER) RETURN
          DEALLOCATE( BUF_MAX_ARRAY )
        ENDIF
        ALLOCATE(BUF_MAX_ARRAY(NFS4FATHER),stat=IERR)
        BUF_LMAX_ARRAY=NFS4FATHER
        RETURN
        END SUBROUTINE DMUMPS_617
        SUBROUTINE DMUMPS_57( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_CB, IERR )
        RETURN
        END SUBROUTINE DMUMPS_57
        SUBROUTINE DMUMPS_59( IERR )
        IMPLICIT NONE
        INTEGER IERR
        CALL DMUMPS_3( BUF_SMALL, IERR )
        RETURN
        END SUBROUTINE DMUMPS_59
        SUBROUTINE DMUMPS_2( BUF, SIZE, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INTEGER SIZE, IERR
        IERR         = 0
        BUF%LBUF     = SIZE
        BUF%LBUF_INT = ( SIZE + SIZEOFINT - 1 ) / SIZEOFINT
        IF ( ASSOCIATED ( BUF%CONTENT ) ) DEALLOCATE( BUF%CONTENT )
        ALLOCATE( BUF%CONTENT( BUF%LBUF_INT ), stat = ierr )
        IF (IERR .NE. 0) THEN
          NULLIFY( BUF%CONTENT )
          IERR         = -1
          BUF%LBUF     =  0
          BUF%LBUF_INT =  0
        END IF
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_2
        SUBROUTINE DMUMPS_3( BUF, IERR )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER STATUS( MPI_STATUS_SIZE )
        LOGICAL FLAG
        IF ( .NOT. ASSOCIATED ( BUF%CONTENT ) ) THEN
          BUF%HEAD     = 1
          BUF%LBUF     = 0
          BUF%LBUF_INT = 0
          BUF%TAIL     = 1
          BUF%ILASTMSG = 1
          RETURN
        END IF
        DO WHILE ( BUF%HEAD.NE.0 .AND. BUF%HEAD .NE. BUF%TAIL )
          CALL MPI_TEST(BUF%CONTENT( BUF%HEAD + REQ ), FLAG,
     *                  STATUS, IERR)
          IF ( .not. FLAG ) THEN
            WRITE(*,*) '** Warning: trying to cancel a request.'
            WRITE(*,*) '** This might be problematic on SGI'
            CALL MPI_CANCEL( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
            CALL MPI_REQUEST_FREE( BUF%CONTENT( BUF%HEAD + REQ ), IERR )
          END IF
          BUF%HEAD = BUF%CONTENT( BUF%HEAD + NEXT )
        END DO
        DEALLOCATE( BUF%CONTENT )
        NULLIFY( BUF%CONTENT )
        BUF%LBUF     = 0
        BUF%LBUF_INT = 0
        BUF%HEAD     = 1
        BUF%TAIL     = 1
        BUF%ILASTMSG = 1
        RETURN
        END SUBROUTINE DMUMPS_3
        SUBROUTINE DMUMPS_66( INODE, FPERE, NFRONT, LCONT,
     *                               NASS, NPIV,
     *                               IWROW, IWCOL, A, COMPRESSCB,
     *                               SIZECB,
     *       DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER DEST, TAG, COMM, IERR
        INTEGER INODE, FPERE, NFRONT, LCONT, NASS, NPIV, SIZECB
        INTEGER IWROW( LCONT ), IWCOL( LCONT )
        DOUBLE PRECISION A( * )
        LOGICAL COMPRESSCB
        INCLUDE 'mpif.h'
        INTEGER POSITION, IREQ, IPOS, I, J1
        INTEGER SIZE1, SIZE2, SIZE
        INTEGER IZERO, IONE
        INTEGER LCONT_SENT
        INTEGER DEST2(1)
        PARAMETER( IZERO = 0, IONE = 1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 9 + LCONT + LCONT, MPI_INTEGER,
     *                      COMM, SIZE1,  IERR)
        CALL MPI_PACK_SIZE( SIZECB, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2,  IERR)
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             RETURN
          ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF (COMPRESSCB) THEN
          LCONT_SENT=-LCONT
        ELSE
          LCONT_SENT=LCONT
        ENDIF
        CALL MPI_PACK( LCONT_SENT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCONT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NASS-NPIV, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LCONT , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IONE,  1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IZERO, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IWROW, LCONT, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IWCOL, LCONT, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LCONT .NE. 0 ) THEN
          J1 = 1
          IF (COMPRESSCB) THEN
           DO I = 1, LCONT
            CALL MPI_PACK( A( J1 ), I, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
             J1 = J1 + NFRONT
           END DO
          ELSE
           DO I = 1, LCONT
            CALL MPI_PACK( A( J1 ), LCONT, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
             J1 = J1 + NFRONT
           END DO
          ENDIF
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_cb: SIZE, POSITION = ', SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_66
        SUBROUTINE DMUMPS_72( NRHS, INODE, IFATH,
     *             EFF_CB_SIZE, LD_CB, NPIV, CB, SOL,
     *             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE, IFATH, EFF_CB_SIZE, LD_CB, NPIV 
        INTEGER DEST, COMM, IERR
        DOUBLE PRECISION CB( LD_CB*NRHS),
     *                   SOL( MAX(1,NPIV*NRHS) )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, SIZE1, SIZE2, K
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        CALL MPI_PACK_SIZE( 4, MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( NRHS * (EFF_CB_SIZE + NPIV),
     *                      MPI_DOUBLE_PRECISION, COMM,
     *                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IFATH, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( EFF_CB_SIZE  , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NPIV , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        DO K = 1, NRHS
               CALL MPI_PACK( CB ( 1 + LD_CB * (K-1) ),
     *                        EFF_CB_SIZE, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        END DO
        IF ( NPIV .GT. 0 ) THEN
          DO K=1, NRHS
          CALL MPI_PACK( SOL(1+NPIV*(K-1)),
     *                         NPIV, MPI_DOUBLE_PRECISION,
     *                         BUF_CB%CONTENT( IPOS ), SIZE,
     *                         POSITION, COMM, IERR )
          ENDDO
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MASTER2SLAVE, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_send_master2slave: SIZE, POSITION = ',
     *               SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_72
        SUBROUTINE DMUMPS_78( NRHS, NODE1, NODE2, NCB, LDW,
     *             LONG,
     *             IW, W,
     *             DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER LDW, DEST, TAG, COMM, IERR
        INTEGER NRHS, NODE1, NODE2, NCB, LONG
        INTEGER IW( MAX( 1, LONG ) )
        DOUBLE PRECISION W( MAX( 1, LDW * NRHS ) )
        INCLUDE 'mpif.h'
        INTEGER POSITION, IREQ, IPOS
        INTEGER SIZE1, SIZE2, SIZE, K
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        IF ( NODE2 .EQ. 0 ) THEN
         CALL MPI_PACK_SIZE( 2+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        ELSE
         CALL MPI_PACK_SIZE( 4+LONG, MPI_INTEGER, COMM, SIZE1, IERR )
        END IF
        SIZE2 = 0
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK_SIZE( NRHS*LONG, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        END IF
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( NODE1, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( NODE2 .NE. 0 ) THEN
          CALL MPI_PACK( NODE2, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        END IF
        CALL MPI_PACK( LONG,  1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LONG .GT. 0 ) THEN
          CALL MPI_PACK( IW, LONG, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          DO K=1, NRHS
          CALL MPI_PACK( W(1+(K-1)*LDW), LONG, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
          END DO
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_78
        SUBROUTINE DMUMPS_62( I, DEST, TAG, COMM, IERR )
        IMPLICIT NONE
        INTEGER I
        INTEGER DEST, TAG, COMM, IERR
        INCLUDE 'mpif.h'
        INTEGER IPOS, IREQ, MSG_SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        MSG_SIZE = SIZEofINT
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, MSG_SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN 
         write(6,*) ' Internal error in DMUMPS_62',
     *       ' Buf size (bytes)= ',BUF_SMALL%LBUF
         RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS ) = I
        CALL MPI_ISEND( BUF_SMALL%CONTENT(IPOS), MSG_SIZE,
     *                  MPI_PACKED, DEST, TAG, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_62
        SUBROUTINE DMUMPS_469(FLAG)
        LOGICAL FLAG
        LOGICAL FLAG1, FLAG2, FLAG3
        CALL DMUMPS_468( BUF_SMALL, FLAG1 )
        CALL DMUMPS_468( BUF_CB, FLAG2 )
        CALL DMUMPS_468( BUF_LOAD, FLAG3 )
        FLAG = FLAG1 .AND. FLAG2 .AND. FLAG3
        RETURN
        END SUBROUTINE DMUMPS_469
        SUBROUTINE DMUMPS_468( B, FLAG )
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        LOGICAL :: FLAG
        INTEGER SIZE_AVAIL
        CALL DMUMPS_79(B, SIZE_AVAIL)
        FLAG = ( B%HEAD == B%TAIL )
        RETURN
        END SUBROUTINE DMUMPS_468
        SUBROUTINE DMUMPS_79( B, SIZE_AV )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        INTEGER SIZE_AV
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER STATUS( MPI_STATUS_SIZE )
        LOGICAL FLAG
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          SIZE_AV = MAX( B%LBUF_INT - B%TAIL + 1, B%HEAD - 1 )
        ELSE
          SIZE_AV = B%HEAD - B%TAIL
        END IF
        SIZE_AV = SIZE_AV * SIZEofINT
        RETURN
        END SUBROUTINE DMUMPS_79
        SUBROUTINE DMUMPS_4( B, IPOS, IREQ, MSG_SIZE, IERR, 
     &    NDEST , PDEST
     &         )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: B
        INTEGER, INTENT(IN)        :: MSG_SIZE
        INTEGER, INTENT(OUT)       :: IPOS, IREQ, IERR
        INTEGER NDEST
        INTEGER, INTENT(IN)        :: PDEST(MAX(1,NDEST))
        INCLUDE 'mpif.h'
        INTEGER MSG_SIZE_INT
        INTEGER IBUF
        LOGICAL FLAG
        INTEGER STATUS( MPI_STATUS_SIZE )
        IERR = 0
        IF ( B%HEAD .NE. B%TAIL ) THEN
 10       CONTINUE
          CALL MPI_TEST( B%CONTENT( B%HEAD + REQ ), FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            B%HEAD = B%CONTENT( B%HEAD + NEXT )
            IF ( B%HEAD .EQ. 0 ) B%HEAD = B%TAIL
            IF ( B%HEAD .NE. B%TAIL ) GOTO 10
          END IF
        END IF
        IF ( B%HEAD .EQ. B%TAIL ) THEN
          B%HEAD = 1
          B%TAIL = 1
          B%ILASTMSG = 1
        END iF
        MSG_SIZE_INT = ( MSG_SIZE + ( SIZEofINT - 1 ) ) / SIZEofINT
        MSG_SIZE_INT = MSG_SIZE_INT + OVHSIZE
        FLAG = (     ( B%HEAD .LE. B%TAIL )
     *               .AND. (
     *                 ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL )
     *                 .OR. ( MSG_SIZE_INT .LE. B%HEAD - 2 ) ) )
     *         .OR.
     *               ( ( B%HEAD .GT. B%TAIL )
     *               .AND. ( MSG_SIZE_INT .LE. B%HEAD - B%TAIL - 1 ) )
        IF ( .NOT. FLAG
     *    ) THEN
          IERR = -1
          IF ( MSG_SIZE_INT .GT. B%LBUF_INT - 1 ) IERR = -2
          IPOS = -1
          IREQ = -1
          RETURN
        END IF
        IF ( B%HEAD .LE. B%TAIL ) THEN
          IF ( MSG_SIZE_INT .LE. B%LBUF_INT - B%TAIL + 1 ) THEN
            IBUF = B%TAIL
          ELSE IF ( MSG_SIZE_INT .LE. B%HEAD - 1 ) THEN
            IBUF = 1
          END IF
        ELSE
          IBUF = B%TAIL
        END IF
        B%CONTENT( B%ILASTMSG + NEXT ) = IBUF
        B%ILASTMSG = IBUF
        B%TAIL = IBUF + MSG_SIZE_INT
        B%CONTENT( IBUF + NEXT ) = 0
        IPOS = IBUF + CONTENT
        IREQ = IBUF + REQ
        RETURN
        END SUBROUTINE DMUMPS_4
        SUBROUTINE DMUMPS_1( BUF, SIZE )
        IMPLICIT NONE
        TYPE ( DMUMPS_BUFFER_TYPE ) :: BUF
        INTEGER SIZE
        INTEGER SIZE_INT
        SIZE_INT = ( SIZE + SIZEofINT - 1 ) / SIZEofINT
        SIZE_INT = SIZE_INT + OVHSIZE
        BUF%TAIL = BUF%ILASTMSG + SIZE_INT
        RETURN
        END SUBROUTINE DMUMPS_1
      SUBROUTINE DMUMPS_68(
     *             INODE, NBPROCFILS, NLIG, ILIG, NCOL, ICOL,
     *             NASS, NSLAVES, LIST_SLAVES,
     *             DEST, NFRONT, COMM, IERR )
      IMPLICIT NONE
        INTEGER COMM, IERR, NFRONT
        INTEGER INODE
        INTEGER NLIG, NCOL, NASS, NSLAVES
        INTEGER NBPROCFILS, DEST
        INTEGER ILIG( NLIG )
        INTEGER ICOL( NCOL )
        INTEGER LIST_SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 6 + NLIG + NCOL + NSLAVES + 1 ) * SIZEofINT
        IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
         IERR = -2
         RETURN
        END IF
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = INODE
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NBPROCFILS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NLIG
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NCOL
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NASS
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NFRONT
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = 
     *  LIST_SLAVES( 1: NSLAVES )
        POSITION = POSITION + NSLAVES
        BUF_CB%CONTENT( POSITION:POSITION + NLIG - 1 ) = ILIG
        POSITION = POSITION + NLIG
        BUF_CB%CONTENT( POSITION:POSITION + NCOL - 1 ) = ICOL
        POSITION = POSITION + NCOL
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_68 :',
     *               ' wrong estimated size'
          CALL DMUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, MPI_PACKED,
     *                  DEST, MAITRE_DESC_BANDE, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_68
        SUBROUTINE DMUMPS_70( IPERE, ISON, NROW,
     *  IROW, NCOL, ICOL, VAL, LDA, NELIM, TYPE_SON,
     *  NSLAVES, SLAVES, DEST, COMM, IERR, 
     * 
     *  SLAVEF, KEEP,KEEP8, INIV2, TAB_POS_IN_PERE )
        IMPLICIT NONE
        INTEGER LDA, NELIM, TYPE_SON
        INTEGER IPERE, ISON, NROW, NCOL, NSLAVES
        INTEGER IROW( NROW )
        INTEGER ICOL( NCOL )
        INTEGER SLAVES( NSLAVES )
        DOUBLE PRECISION VAL(LDA, *)
        INTEGER IPOS, IREQ, DEST, COMM, IERR
        INTEGER SLAVEF, KEEP(500), INIV2
      INTEGER*8 KEEP8(150)
        INTEGER TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE1, SIZE2, SIZE3, SIZE, POSITION, I
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        IF ( NELIM .NE. NROW ) THEN
          WRITE(*,*) 'Error in TRY_SEND_MAITRE2:',NELIM, NROW
          CALL DMUMPS_ABORT()
        END IF
        CALL MPI_PACK_SIZE( NROW+NCOL+5+NSLAVES, MPI_INTEGER,
     *                      COMM, SIZE1, IERR )
        IF ( ( KEEP(48).NE. 0 ).AND.(TYPE_SON .eq. 2 )) THEN
          CALL MPI_PACK_SIZE( NSLAVES+1, MPI_INTEGER,
     *                        COMM, SIZE3, IERR )
        ELSE
          SIZE3 = 0
        ENDIF
        IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
          CALL MPI_PACK_SIZE( NROW * NROW, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        ELSE
          CALL MPI_PACK_SIZE( NROW * NCOL, MPI_DOUBLE_PRECISION,
     *                        COMM, SIZE2, IERR )
        END IF
        SIZE = SIZE1 + SIZE3 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        IF (NELIM.NE.0) THEN
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -2
             RETURN
          ENDIF
        ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON,  1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF (NSLAVES.GT.0) THEN
         CALL MPI_PACK( SLAVES, NSLAVES, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ENDIF
        CALL MPI_PACK( IROW, NROW, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ICOL, NCOL, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF ( ( KEEP(48).NE. 0 ).AND.(TYPE_SON .eq. 2 ) ) THEN
         CALL MPI_PACK( TAB_POS_IN_PERE(1,INIV2), NSLAVES+1, 
     *                 MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ENDIF
        IF (NROW.GE.1) THEN
          IF ( KEEP(50).ne.0  .AND. TYPE_SON .eq. 2 ) THEN
            DO I=1,NROW
              CALL MPI_PACK( VAL(1,I), NROW, 
     *                 MPI_DOUBLE_PRECISION,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
            ENDDO
          ELSE
            DO I=1,NROW
              CALL MPI_PACK( VAL(1,I), NCOL, 
     *                 MPI_DOUBLE_PRECISION,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
            ENDDO
          ENDIF
        ENDIF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MAITRE2, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          write(*,*) 'Try_send_maitre2, SIZE,POSITION=',SIZE,POSITION
          Stop
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_70
        SUBROUTINE DMUMPS_67(DESC_IN_LU,
     *  IPERE, NFRONT_PERE, NASS_PERE, NFS4FATHER,
     * NSLAVES_PERE,
     *  ISON, NBROW, LMAP, MAPROW, PERM, IW_CBSON, A_CBSON,
     *  ISLAVE, PDEST, PDEST_MASTER, COMM, IERR, 
     *  
     * KEEP,KEEP8, STEP, N, SLAVEF,
     * ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     * MAX_ROW_LENGTH, SIZE_SLAVE_REALS, COMPRESSCB )
        IMPLICIT NONE
        INTEGER IPERE, ISON, NBROW 
        INTEGER PDEST, ISLAVE, COMM, IERR
        INTEGER PDEST_MASTER, NASS_PERE, NSLAVES_PERE,
     *       NFRONT_PERE, LMAP
        INTEGER MAPROW( LMAP ), PERM( MAX(1, NBROW ))
        INTEGER MAX_ROW_LENGTH, SIZE_SLAVE_REALS
        INTEGER IW_CBSON( * )
        DOUBLE PRECISION A_CBSON( * )
        LOGICAL DESC_IN_LU, COMPRESSCB
       INTEGER   KEEP(500), N , SLAVEF
       INTEGER*8 KEEP8(150)
       INTEGER   STEP(N), 
     *          ISTEP_TO_INIV2(KEEP(71)), 
     *          TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
       INTEGER NFS4FATHER,SIZE3,PS1,NCA,OFFSET,ASIZE,LROW1
       LOGICAL COMPUTE_MAX
       DOUBLE PRECISION TMP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER LROW, NELIM, ITMP, SIZFR
      INTEGER NPIV, NFRONT, HS
      INTEGER SIZE, SIZE1, SIZE2, POSITION,I
      INTEGER APOS, IPOS_IN_SLAVE
      INTEGER LDA_SON, SHIFTCB_SON, STATE_SON
      INTEGER INDICE_PERE, NROW, IPOS, IREQ, NOSLA
      INTEGER IONE, J, THIS_ROW_LENGTH
      INTEGER PDEST2(1)
      PARAMETER ( IONE=1 )
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0D0)
      IF (KEEP(219).NE.0) THEN
        CALL DMUMPS_617(NFS4FATHER,IERR)
        IF (IERR .NE. 0) THEN
          IERR         = -4
          RETURN
        ENDIF
      ENDIF
        PDEST2(1) = PDEST
        IERR   = 0
        LROW   = IW_CBSON( 1 + XSIZE)
        NELIM  = IW_CBSON( 2 + XSIZE)
        NPIV   = IW_CBSON( 4 + XSIZE)
        IF ( NPIV .LT. 0 ) THEN
          NPIV = 0
        END IF
        NROW   = IW_CBSON( 3 + XSIZE)
        NFRONT = LROW + NPIV
        HS     = 6 + IW_CBSON( 6 + XSIZE) + XSIZE
        SIZFR  = IW_CBSON( 1 + XXR)
        STATE_SON = IW_CBSON(1+XXS)
        IF (STATE_SON .EQ. S_NOLCBCONTIG) THEN
               LDA_SON     = LROW
               SHIFTCB_SON = NPIV*NROW
        ELSE IF (STATE_SON .EQ. S_NOLCLEANED) THEN
               LDA_SON     = LROW
               SHIFTCB_SON = 0
        ELSE
               LDA_SON      = NFRONT
               SHIFTCB_SON = NPIV
        ENDIF
        IF ( PDEST .EQ. PDEST_MASTER ) THEN
          IF (KEEP(50) == 0 ) THEN
          CALL MPI_PACK_SIZE( NBROW + 4, MPI_INTEGER,
     *                        COMM, SIZE1, IERR )
          ELSE
          CALL MPI_PACK_SIZE( 2*NBROW + 4, MPI_INTEGER,
     *                        COMM, SIZE1, IERR )
          ENDIF
        ELSE
          IF (KEEP(50) == 0 ) THEN
            CALL MPI_PACK_SIZE( NBROW + 4 + LROW, MPI_INTEGER,
     *                        COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 2*NBROW + 4 + MAX_ROW_LENGTH,
     *                       MPI_INTEGER, COMM, SIZE1, IERR )
          ENDIF
        END IF
        IF (KEEP(50)==0) THEN
         IF (NBROW > 0) THEN
         CALL MPI_PACK_SIZE( LROW * NBROW, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
         ELSE
         SIZE2=0
         ENDIF
        ELSE
         IF (SIZE_SLAVE_REALS > 0) THEN
         CALL MPI_PACK_SIZE( SIZE_SLAVE_REALS, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
         ELSE
         SIZE2=0
         ENDIF
        ENDIF
        IF (KEEP(219).NE.0) THEN
          COMPUTE_MAX = KEEP(50).EQ.2 .AND.
     *         PDEST.EQ.PDEST_MASTER
          IF(COMPUTE_MAX) THEN
             CALL MPI_PACK_SIZE(1, MPI_INTEGER,
     *          COMM, PS1, IERR )
             IF(NFS4FATHER .GT. 0) THEN
                CALL MPI_PACK_SIZE( NFS4FATHER, MPI_DOUBLE_PRECISION,
     *             COMM, SIZE3, IERR )
             ELSE
                SIZE3=0
             ENDIF
             SIZE3 = SIZE3+PS1
          ELSE
             SIZE3 = 0
          ENDIF
          SIZE = SIZE1 + SIZE2 + SIZE3
        ELSE
          SIZE = SIZE1 + SIZE2
        ENDIF
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE , PDEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             RETURN
          ENDIF
        POSITION = 0
        CALL MPI_PACK( IPERE, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( NBROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        IF (KEEP(50)==0) THEN
        CALL MPI_PACK( LROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ELSE
        CALL MPI_PACK( MAX_ROW_LENGTH, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
        ENDIF
        IF ( PDEST .NE. PDEST_MASTER ) THEN
          IF (KEEP(50)==0) THEN
          CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ), LROW,
     *                 MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
          ELSE
           IF (MAX_ROW_LENGTH > 0) THEN
           CALL MPI_PACK( IW_CBSON( HS + NROW +  NPIV + 1 ),
     *                 MAX_ROW_LENGTH,
     *                 MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ), SIZE,
     *                 POSITION, COMM, IERR )
           ENDIF
          ENDIF
        END IF
        DO J=1, NBROW
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL DMUMPS_47(
     *          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *
     *          NASS_PERE,
     *          NFRONT_PERE - NASS_PERE,
     *          NSLAVES_PERE,
     *          INDICE_PERE,
     *          NOSLA,
     *          IPOS_IN_SLAVE )
           INDICE_PERE = IPOS_IN_SLAVE
           CALL MPI_PACK( INDICE_PERE, 1, MPI_INTEGER,
     *          BUF_CB%CONTENT( IPOS ), SIZE,
     *          POSITION, COMM, IERR )
        ENDDO
        DO J=1, NBROW
           I = PERM(J)
           INDICE_PERE=MAPROW(I)
           CALL DMUMPS_47(
     *          KEEP,KEEP8, IPERE, STEP, N, SLAVEF,
     *          ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     *          
     *          NASS_PERE,
     *          NFRONT_PERE - NASS_PERE,
     *          NSLAVES_PERE,
     *          INDICE_PERE,
     *          NOSLA,
     *          IPOS_IN_SLAVE )
          IF (KEEP(50).ne.0) THEN
            THIS_ROW_LENGTH = LROW + I - LMAP
            CALL MPI_PACK( THIS_ROW_LENGTH, 1, MPI_INTEGER,
     *                      BUF_CB%CONTENT( IPOS ), SIZE,
     *                      POSITION, COMM, IERR )
         ELSE
            THIS_ROW_LENGTH = LROW
         ENDIF
         IF (DESC_IN_LU) THEN 
            IF ( COMPRESSCB ) THEN
             IF (NELIM.EQ.0) THEN
               APOS = I * (I-1) / 2 + 1
             ELSE
               APOS = (NELIM+I)*(NELIM+I-1)/2 + 1
             ENDIF
            ELSE
             APOS = ( I + NELIM - 1 ) * LROW + 1
            ENDIF
         ELSE
            IF ( COMPRESSCB ) THEN
             IF ( LROW .EQ. NROW )  THEN
               APOS = I * (I-1)/2 + 1
             ELSE
               ITMP = I + LROW - NROW
               APOS = ITMP * (ITMP-1)/2 + 1 -
     *              ( LROW - NROW ) * ( LROW - NROW + 1 )/2
             ENDIF
            ELSE
             APOS = ( I - 1 ) * LDA_SON + SHIFTCB_SON + 1
            ENDIF
         ENDIF
         CALL MPI_PACK( A_CBSON( APOS ), THIS_ROW_LENGTH,
     *        MPI_DOUBLE_PRECISION,
     *        BUF_CB%CONTENT( IPOS ), SIZE,
     *        POSITION, COMM, IERR )
        ENDDO
      IF (KEEP(219).NE.0) THEN
        IF(COMPUTE_MAX) THEN
           CALL MPI_PACK(NFS4FATHER,1,
     *          MPI_INTEGER,
     *          BUF_CB%CONTENT( IPOS ), SIZE,
     *          POSITION, COMM, IERR )
           IF(NFS4FATHER .GT. 0) THEN
              BUF_MAX_ARRAY(1:NFS4FATHER) = ZERO
              IF(MAPROW(NROW) .GT. NASS_PERE) THEN
                 DO PS1=1,NROW
                    IF(MAPROW(PS1).GT.NASS_PERE) EXIT
                 ENDDO
                 IF (DESC_IN_LU) THEN
                   IF (COMPRESSCB) THEN
                    APOS   = (NELIM+PS1) * (NELIM+PS1-1) / 2 + 1
                    NCA    = -44444
                    ASIZE  = NROW * (NROW+1)/2 -
     *                       (NELIM+PS1) * (NELIM+PS1-1) / 2
                    LROW1  = PS1
                   ELSE
                    APOS = ( PS1 + NELIM - 1 ) * LROW + 1
                    NCA = LROW 
                    ASIZE = NCA * (NROW-PS1+1)
                    LROW1 = LROW
                   ENDIF
                 ELSE
                    IF (COMPRESSCB) THEN
                      IF (NPIV.NE.0) THEN
         WRITE(*,*) "Error in PARPIV/DMUMPS_67"
                        CALL DMUMPS_ABORT()
                      ENDIF
                      LROW1=LROW-NROW+1
                      ITMP = PS1 + LROW - NROW
                      APOS = ITMP * (ITMP - 1) / 2 + 1 -
     *                      (LROW-NROW)*(LROW-NROW+1)/2
                      ASIZE = LROW*(LROW+1)/2-ITMP*(ITMP-1)/2
                      NCA   = -555555
                    ELSE
                      APOS = ( PS1 - 1 ) * LDA_SON + 1 + SHIFTCB_SON
                      NCA = LDA_SON
                      ASIZE = SIZFR - (SHIFTCB_SON - (PS1-1) * LDA_SON)
                      LROW1=-666666
                    ENDIF
                 ENDIF
                 IF ( NROW-PS1+1 .NE. 0 ) THEN
                   CALL DMUMPS_618(
     *                A_CBSON(APOS),ASIZE,NCA,NROW-PS1+1,
     *                BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,LROW1)
                 ENDIF
              ENDIF
              CALL MPI_PACK(BUF_MAX_ARRAY, NFS4FATHER,
     *             MPI_DOUBLE_PRECISION,
     *             BUF_CB%CONTENT( IPOS ), SIZE,
     *             POSITION, COMM, IERR )
           ENDIF
        ENDIF
      ENDIF !PARPIV
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  PDEST, CONTSIPERENIV2, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' contniv2: SIZE, POSITION =',SIZE, POSITION
          WRITE(*,*) ' NBROW, LROW = ', NBROW, LROW
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_67
        SUBROUTINE DMUMPS_71( 
     *                INODE, NFRONT, NASS1, NFS4FATHER,
     *                ISON, MYID, NSLAVES, SLAVES_PERE,
     *                TROW, NCBSON,
     *                COMM, IERR,
     *                DEST, NDEST, SLAVEF, 
     * 
     *                KEEP,KEEP8, STEP, N, 
     *                ISTEP_TO_INIV2, TAB_POS_IN_PERE
     *
     *                                  )
        IMPLICIT NONE
      INTEGER INODE, NFRONT, NASS1, NCBSON, NSLAVES, 
     *          NDEST
      INTEGER SLAVEF, MYID, ISON
      INTEGER TROW( NCBSON )
      INTEGER DEST( NDEST )
      INTEGER SLAVES_PERE( NSLAVES )
      INTEGER COMM, IERR
      INTEGER KEEP(500), N
      INTEGER*8 KEEP8(150)
      INTEGER STEP(N), 
     *        ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        INTEGER SIZE_AV, IDEST, NSEND, SIZE, NFS4FATHER
        INTEGER TROW_SIZE, POSITION, INDX, INIV2
        INTEGER IPOS, IREQ
        INTEGER IONE
        PARAMETER ( IONE=1 )
        INTEGER NASS_SON
        NASS_SON = -99998
        IERR = 0
        IF ( NDEST .eq. 1 ) THEN
          IF ( DEST(1).EQ.MYID )  GOTO 500
          SIZE = SIZEofINT * ( 7 + NSLAVES + NCBSON )
          IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST
     &                 )
          IF (IERR .LT. 0 ) THEN
             RETURN
          ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             RETURN
          END IF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NCBSON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     *          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     *          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+NCBSON-1 ) =
     *        TROW( 1: NCBSON )
              POSITION = POSITION + NCBSON
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
                WRITE(*,*) 'Error in DMUMPS_71 :',
     *                     ' wrong estimated size'
                CALL DMUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     *                        MPI_PACKED,
     *                        DEST( NDEST ), MAPLIG, COMM,
     *                        BUF_CB%CONTENT( IREQ ),
     *                        IERR )
        ELSE
          NSEND = 0
          DO IDEST = 1, NDEST
            IF ( DEST( IDEST ) .ne. MYID ) NSEND = NSEND + 1
          END DO
          SIZE = SIZEofINT * 
     *         ( ( OVHSIZE + 7 + NSLAVES )* NSEND + NCBSON )
          IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
           SIZE = SIZE + SIZEofINT * NSEND*( NSLAVES + 1 )
          ENDIF
          CALL DMUMPS_79( BUF_CB, SIZE_AV )
          IF ( SIZE_AV .LT. SIZE ) THEN
            IERR = -1
            RETURN
          END IF
          DO IDEST= 1, NDEST
            CALL DMUMPS_49( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                IDEST, NCBSON, 
     &                NDEST, 
     &                TROW_SIZE, INDX  )
            SIZE = SIZEofINT * ( NSLAVES + TROW_SIZE + 7 )
            IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
             SIZE = SIZE + SIZEofINT * ( NSLAVES + 1 )
            ENDIF
            IF ( MYID .NE. DEST( IDEST ) ) THEN
              CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR,
     &                       IONE, DEST(IDEST)
     &                     )
              IF ( IERR .LT. 0 )  THEN
                 WRITE(*,*) 'Problem in DMUMPS_4: IERR<0'
                 CALL DMUMPS_ABORT()
              END IF
              IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
                 IERR = -3
                 RETURN
              ENDIF
              POSITION = IPOS
              BUF_CB%CONTENT( POSITION ) = INODE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = ISON
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NSLAVES
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFRONT
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NASS1
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = TROW_SIZE
              POSITION = POSITION + 1
              BUF_CB%CONTENT( POSITION ) = NFS4FATHER
              POSITION = POSITION + 1
              IF (( NSLAVES.GT.0 ) .AND. (KEEP(48).NE.0) ) THEN
                INIV2 = ISTEP_TO_INIV2 ( STEP(INODE) )
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES )
     *          =  TAB_POS_IN_PERE(1:NSLAVES+1,INIV2)
                POSITION = POSITION + NSLAVES + 1
              ENDIF
              IF ( NSLAVES .NE. 0 ) THEN
                BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 )
     *          = SLAVES_PERE( 1: NSLAVES )
                POSITION = POSITION + NSLAVES
              END IF
              BUF_CB%CONTENT( POSITION:POSITION+TROW_SIZE-1 ) =
     *        TROW( INDX: INDX + TROW_SIZE - 1 )
              POSITION = POSITION + TROW_SIZE
              POSITION = POSITION - IPOS
              IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
               WRITE(*,*) ' ERROR 1 in TRY_SEND_MAPLIG:',
     *          'Wrong estimated size'
               CALL DMUMPS_ABORT()
              END IF
              CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE,
     *                        MPI_PACKED,
     *                        DEST( IDEST ), MAPLIG, COMM,
     *                        BUF_CB%CONTENT( IREQ ),
     *                        IERR )
            END IF
          END DO
        END IF
 500    CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_71
        SUBROUTINE DMUMPS_65( INODE, NFRONT,
     *             NCOL, NPIV, FPERE, LASTBL, IPIV, VAL,
     *             PDEST, NDEST, KEEP50, NB_BLOC_FAC, COMM, IERR )
        IMPLICIT NONE
        INTEGER INODE, NCOL, NPIV, FPERE, NFRONT, NDEST
        INTEGER IPIV( NPIV )
        DOUBLE PRECISION VAL( NFRONT, * )
        INTEGER PDEST( NDEST ) 
        INTEGER   KEEP50, NB_BLOC_FAC, COMM, IERR
        LOGICAL LASTBL
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE,
     *          IDEST, IPOSMSG, I
        INTEGER NPIVSENT
        INTEGER SSS, SS2
        IERR = 0
        IF ( LASTBL ) THEN
          IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        ELSE
          IF ( KEEP50 .eq. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                        MPI_INTEGER, COMM, SIZE1, IERR )
          ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV + ( NDEST - 1 ) * OVHSIZE,
     *                          MPI_INTEGER, COMM, SIZE1, IERR )
          END IF
        END IF
        SIZE2 = 0
        IF (NPIV.GT.0) 
     *    CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 NDEST , PDEST
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
          SSS = 0 
          IF ( LASTBL ) THEN
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 4 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 6 + NPIV , 
     *                           MPI_INTEGER, COMM, SSS, IERR )
           END IF
          ELSE
           IF ( KEEP50 .eq. 0 ) THEN
            CALL MPI_PACK_SIZE( 3 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           ELSE
            CALL MPI_PACK_SIZE( 4 + NPIV ,
     *                        MPI_INTEGER, COMM, SSS, IERR )
           END IF
          END IF
          IF (NPIV.GT.0) 
     *    CALL MPI_PACK_SIZE( NPIV*NCOL, MPI_DOUBLE_PRECISION,
     *                      COMM, SS2, IERR )
          SSS = SSS + SS2
          IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2 
           RETURN
          ENDIF
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        NPIVSENT = NPIV
        IF (LASTBL) NPIVSENT = -NPIV
        CALL MPI_PACK( NPIVSENT, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( LASTBL .or. KEEP50.ne.0 ) THEN
          CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        IF ( LASTBL .AND. KEEP50 .NE. 0 ) THEN
            CALL MPI_PACK( NDEST, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
            CALL MPI_PACK( NB_BLOC_FAC, 1, MPI_INTEGER,
     *                   BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        CALL MPI_PACK( NCOL, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        IF ( NPIV.GT.0) THEN
          CALL MPI_PACK( IPIV, NPIV, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
          DO I = 1, NPIV
            CALL MPI_PACK( VAL(1,I), NCOL,
     *                        MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
          END DO
        ENDIF
        DO IDEST = 1, NDEST
        IF ( KEEP50.eq.0) THEN
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        ELSE
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO_SYM, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blocfacto : size < position'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_65
        SUBROUTINE DMUMPS_64( INODE,
     *             NPIV, FPERE, IPOSK, JPOSK, UIP21K, NCOLU,
     *             NDEST, PDEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER INODE, NCOLU, IPOSK, JPOSK, NPIV, NDEST, FPERE
        DOUBLE PRECISION UIP21K( NPIV, NCOLU )
        INTEGER PDEST( NDEST ) 
        INTEGER   COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE,
     *          IDEST, IPOSMSG, SSS, SS2
        IERR = 0
        CALL MPI_PACK_SIZE( 6 + ( NDEST - 1 ) * OVHSIZE,
     *                      MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( ABS(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     *                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
         CALL MPI_PACK_SIZE( 6 ,
     *                      MPI_INTEGER, COMM, SSS, IERR )
         CALL MPI_PACK_SIZE( ABS(NPIV)*NCOLU, MPI_DOUBLE_PRECISION,
     *                      COMM, SS2, IERR )
         SSS = SSS+SS2
         IF (SSS.GT.SIZE_RBUF_BYTES) THEN
           IERR = -2
           RETURN
         ENDIF
        END IF
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 NDEST, PDEST
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_CB%ILASTMSG = BUF_CB%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_CB%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_CB%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( IPOSK, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( JPOSK, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NPIV, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( FPERE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( NCOLU, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( UIP21K, ABS(NPIV) * NCOLU,
     *                        MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOSMSG ), SIZE,
     *                        POSITION, COMM, IERR )
        DO IDEST = 1, NDEST
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOSMSG ), POSITION, MPI_PACKED,
     *                  PDEST(IDEST), BLOC_FACTO_SYM_SLAVE, COMM,
     *                  BUF_CB%CONTENT( IREQ + ( IDEST-1 ) * OVHSIZE ),
     *                  IERR )
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending blfac slave : size < position'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_64
        SUBROUTINE DMUMPS_75( N, ISON,
     *             NBCOL_SON, NBROW_SON, INDCOL_SON, INDROW_SON,
     *             LD_SON, VAL_SON, TAG, SUBSET_ROW, SUBSET_COL,
     *             NSUBSET_ROW, NSUBSET_COL,
     *             NPROW, NPCOL, MBLOCK, RG2L_ROW, RG2L_COL,
     *             NBLOCK, PDEST, COMM, IERR , 
     *             TAB, TABSIZE, INVERT, SIZE_MSG )
        IMPLICIT NONE
        INTEGER N, ISON, NBCOL_SON, NBROW_SON, NSUBSET_ROW, NSUBSET_COL
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK, LD_SON
        INTEGER PDEST, TAG, COMM, IERR
        INTEGER INDCOL_SON( NBCOL_SON ), INDROW_SON( NBROW_SON )
        INTEGER SUBSET_ROW( * ), SUBSET_COL( * )
        INTEGER RG2L_ROW( N ), RG2L_COL( N ) 
        INTEGER   TABSIZE, SIZE_MSG
        DOUBLE PRECISION VAL_SON( LD_SON, * ), TAB(*)
        LOGICAL INVERT
        INCLUDE 'mpif.h'
        INTEGER SIZE1, SIZE2, SIZE, POSITION
        INTEGER IREQ, IPOS, ITAB
        INTEGER ISUB, JSUB, I, J 
        INTEGER ILOC_ROOT, JLOC_ROOT
        INTEGER IPOS_ROOT, JPOS_ROOT
        INTEGER IONE
        INTEGER PDEST2(1)
        PARAMETER ( IONE=1 )
        PDEST2(1) = PDEST
        IERR = 0
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          CALL MPI_PACK_SIZE( 3 + NSUBSET_ROW + NSUBSET_COL,
     *                      MPI_INTEGER, COMM, SIZE1, IERR )
          CALL MPI_PACK_SIZE( NSUBSET_ROW * NSUBSET_COL,
     *                      MPI_DOUBLE_PRECISION, COMM, SIZE2, IERR )
          SIZE = SIZE1 + SIZE2
        ELSE
          CALL MPI_PACK_SIZE(3,MPI_INTEGER, COMM, SIZE, IERR )
        END IF
        SIZE_MSG = SIZE
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, PDEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        IF (SIZE.GT.SIZE_RBUF_BYTES ) THEN
             IERR = -3
             RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( ISON, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        IF ( .NOT. INVERT ) THEN
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        ELSE
        CALL MPI_PACK( NSUBSET_COL, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        CALL MPI_PACK( NSUBSET_ROW, 1, MPI_INTEGER,
     *                 BUF_CB%CONTENT( IPOS ),
     *                 SIZE, POSITION, COMM, IERR )
        END IF
        IF ( NSUBSET_ROW * NSUBSET_COL .NE. 0 ) THEN
          IF ( .NOT. INVERT ) THEN
            DO ISUB = 1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              IPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              ILOC_ROOT = MBLOCK
     *                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *                 + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     *                      BUF_CB%CONTENT( IPOS ),
     *                      SIZE, POSITION, COMM, IERR )
            END DO
            DO JSUB = 1, NSUBSET_COL
              J         = SUBSET_COL( JSUB )
              JPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              JLOC_ROOT = NBLOCK
     *                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                  + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     *                       BUF_CB%CONTENT( IPOS ),
     *                       SIZE, POSITION, COMM, IERR )
            END DO
          ELSE
            DO JSUB = 1, NSUBSET_COL
              J         = SUBSET_COL( JSUB )
              IPOS_ROOT = RG2L_COL( INDCOL_SON( J ) )
              ILOC_ROOT = MBLOCK
     *                 * ( ( IPOS_ROOT - 1 ) / ( MBLOCK * NPROW ) )
     *                 + MOD( IPOS_ROOT - 1, MBLOCK ) + 1
              CALL MPI_PACK( ILOC_ROOT, 1, MPI_INTEGER,
     *                       BUF_CB%CONTENT( IPOS ),
     *                       SIZE, POSITION, COMM, IERR )
            END DO
            DO ISUB = 1, NSUBSET_ROW
              I         = SUBSET_ROW( ISUB )
              JPOS_ROOT = RG2L_ROW( INDROW_SON( I ) )
              JLOC_ROOT = NBLOCK
     *                  * ( ( JPOS_ROOT - 1 ) / ( NBLOCK * NPCOL ) )
     *                  + MOD( JPOS_ROOT - 1, NBLOCK ) + 1
              CALL MPI_PACK( JLOC_ROOT, 1, MPI_INTEGER,
     *                      BUF_CB%CONTENT( IPOS ),
     *                      SIZE, POSITION, COMM, IERR )
            END DO
          END IF
          IF ( TABSIZE.GE.NSUBSET_ROW*NSUBSET_COL ) THEN
            IF ( .NOT. INVERT ) THEN
              ITAB = 1
              DO ISUB = 1, NSUBSET_ROW
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL
                  J              = SUBSET_COL( JSUB )
                  TAB( ITAB )    = VAL_SON(J,I)
                  ITAB           = ITAB + 1
                END DO
               END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL*NSUBSET_ROW,
     *         MPI_DOUBLE_PRECISION, 
     *         BUF_CB%CONTENT( IPOS ),
     *         SIZE, POSITION, COMM, IERR )
            ELSE
              ITAB = 1
              DO JSUB = 1, NSUBSET_COL
                J = SUBSET_COL( JSUB )
                DO ISUB = 1, NSUBSET_ROW
                  I         = SUBSET_ROW( ISUB )
                  TAB( ITAB ) = VAL_SON( J, I )
                  ITAB = ITAB + 1
                END DO
              END DO
              CALL MPI_PACK(TAB(1), NSUBSET_COL*NSUBSET_ROW,
     *         MPI_DOUBLE_PRECISION, 
     *         BUF_CB%CONTENT( IPOS ),
     *         SIZE, POSITION, COMM, IERR )
            END IF
          ELSE
            IF ( .NOT. INVERT ) THEN
              DO ISUB = 1, NSUBSET_ROW
                I         = SUBSET_ROW( ISUB )
                DO JSUB = 1, NSUBSET_COL
                  J         = SUBSET_COL( JSUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     *            MPI_DOUBLE_PRECISION,
     *            BUF_CB%CONTENT( IPOS ),
     *            SIZE, POSITION, COMM, IERR )
                END DO
              END DO
            ELSE
              DO JSUB = 1, NSUBSET_COL
                J = SUBSET_COL( JSUB )
                DO ISUB = 1, NSUBSET_ROW
                  I         = SUBSET_ROW( ISUB )
                  CALL MPI_PACK( VAL_SON( J, I ), 1,
     *            MPI_DOUBLE_PRECISION,
     *            BUF_CB%CONTENT( IPOS ),
     *            SIZE, POSITION, COMM, IERR )
                END DO
              END DO
            END IF
          ENDIF
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                PDEST, TAG, COMM, BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error sending contribution to root:Size<positn'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .LT. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_75
        SUBROUTINE DMUMPS_76( ISON, NELIM,
     *             NELIM_ROW, NELIM_COL, NSLAVES, SLAVES,
     *             DEST, COMM, IERR )
        INTEGER ISON, NELIM
        INTEGER NSLAVES, DEST, COMM, IERR
        INTEGER NELIM_ROW( NELIM ), NELIM_COL( NELIM )
        INTEGER SLAVES( NSLAVES )
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, POSITION, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1) = DEST
        IERR = 0
        SIZE = ( 3 + NSLAVES + 2 * NELIM ) * SIZEofINT
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
          IF (SIZE.GT.SIZE_RBUF_BYTES) THEN
             IERR = -3
             RETURN
          ENDIF
        POSITION = IPOS
        BUF_CB%CONTENT( POSITION ) = ISON
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NELIM
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION ) = NSLAVES
        POSITION = POSITION + 1
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_ROW
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NELIM - 1 ) = NELIM_COL
        POSITION = POSITION + NELIM
        BUF_CB%CONTENT( POSITION: POSITION + NSLAVES - 1 ) = SLAVES
        POSITION = POSITION + NSLAVES
        POSITION = POSITION - IPOS
        IF ( POSITION * SIZEofINT .NE. SIZE ) THEN
          WRITE(*,*) 'Error in DMUMPS_BUF_SEND_ROOT_NELIM_INDICES:',
     *               'wrong estimated size'
           CALL DMUMPS_ABORT()
        END IF
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_NELIM_INDICES, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_76
        SUBROUTINE DMUMPS_74( ISON, NELIM_ROOT,
     *             DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER ISON, NELIM_ROOT, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, SIZE
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        DEST2(1)=DEST
        IERR = 0
        SIZE = 2 * SIZEofINT
  10    CONTINUE
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
          WRITE(*,*) 'Increase small buffer size. Trying again.'
          GOTO 10
        END IF
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS )     = ISON
        BUF_SMALL%CONTENT( IPOS + 1 ) = NELIM_ROOT
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_2SON, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_74
        SUBROUTINE DMUMPS_73
     *  ( TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR )
        IMPLICIT NONE
        INTEGER TOT_ROOT_SIZE, TOT_CONT2RECV, DEST, COMM, IERR
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER SIZE, IPOS, IREQ
        INTEGER IONE
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        SIZE = 2 * SIZEofINT
        CALL DMUMPS_4( BUF_SMALL, IPOS, IREQ, SIZE, IERR,
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           WRITE(*,*)  'Increase small buf. size'
           CALL DMUMPS_ABORT()
        END IF
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_SMALL%CONTENT( IPOS     ) = TOT_ROOT_SIZE
        BUF_SMALL%CONTENT( IPOS + 1 ) = TOT_CONT2RECV
        CALL MPI_ISEND( BUF_SMALL%CONTENT( IPOS ), SIZE, 
     *                  MPI_PACKED,
     *                  DEST, ROOT_2SLAVE, COMM,
     *                  BUF_SMALL%CONTENT( IREQ ), IERR )
        RETURN
        END SUBROUTINE DMUMPS_73
        SUBROUTINE DMUMPS_63
     *             ( NRHS, INODE, W, LW, LD_W, DEST,MSGTAG,COMM,IERR )
        IMPLICIT NONE
        INTEGER NRHS, INODE,LW,COMM,IERR,DEST,MSGTAG, LD_W
        DOUBLE PRECISION W(LD_W, *)
        INCLUDE 'mpif.h'
        INTEGER SIZE, SIZE1, SIZE2
        INTEGER POSITION, IREQ, IPOS
        INTEGER IONE, K
        INTEGER DEST2(1)
        PARAMETER ( IONE=1 )
        IERR = 0
        DEST2(1) = DEST
        CALL MPI_PACK_SIZE( 2 , MPI_INTEGER, COMM, SIZE1, IERR )
        CALL MPI_PACK_SIZE( LW*NRHS, MPI_DOUBLE_PRECISION, COMM,
     *                      SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_CB, IPOS, IREQ, SIZE, IERR, 
     &                 IONE, DEST2
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        POSITION = 0
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        CALL MPI_PACK( LW   , 1, MPI_INTEGER,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        DO K=1, NRHS
        CALL MPI_PACK( W(1,K), LW, MPI_DOUBLE_PRECISION,
     *                        BUF_CB%CONTENT( IPOS ), SIZE,
     *                        POSITION, COMM, IERR )
        END DO
        CALL MPI_ISEND( BUF_CB%CONTENT( IPOS ), POSITION, MPI_PACKED,
     *                  DEST, MSGTAG, COMM,
     *                  BUF_CB%CONTENT( IREQ ), IERR )
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) 'Try_update: SIZE, POSITION = ',
     *               SIZE, POSITION
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION ) CALL DMUMPS_1( BUF_CB, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_63
        SUBROUTINE DMUMPS_77
     *             ( BDC_SBTR,BDC_MEM,BDC_MD, COMM, NPROCS, LOAD,
     *               MEM,SBTR_CUR,
     *               LU_USAGE,
#if ! defined(OLD_LOAD_MECHANISM)
     *               FUTURE_NIV2,
#endif
     *               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        DOUBLE PRECISION LU_USAGE
        DOUBLE PRECISION LOAD
        DOUBLE PRECISION MEM,SBTR_CUR
        LOGICAL BDC_MEM,BDC_SBTR,BDC_MD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, WHAT, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#endif
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        NREALS = 1
        IF (BDC_MEM) THEN
          NREALS = 2
        ENDIf
        IF (BDC_SBTR)THEN
          NREALS = 3
        ENDIF
        IF(BDC_MD)THEN
           NREALS=NREALS+1
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        WHAT = 0  ! 0: absolute load, only sender makes the update
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM, 1, MPI_DOUBLE_PRECISION,
     *                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        IF (BDC_SBTR) THEN
          CALL MPI_PACK( SBTR_CUR, 1, MPI_DOUBLE_PRECISION,
     *                   BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                   POSITION, COMM, IERR )
        END IF
        IF(BDC_MD)THEN
           CALL MPI_PACK( LU_USAGE, 1, MPI_DOUBLE_PRECISION,
     *          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .ne. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_77'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_77
        SUBROUTINE DMUMPS_460
     *             ( WHAT, COMM, NPROCS, 
#if ! defined(OLD_LOAD_MECHANISM)
     *               FUTURE_NIV2,
#endif
     *               LOAD,UPD_LOAD,
     *               MYID, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT
        DOUBLE PRECISION LOAD,UPD_LOAD
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE1, SIZE2, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG, NREALS
        INTEGER IZERO
        INTEGER MYID2(1)
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        PARAMETER ( IZERO=0 )
        IERR = 0
        IF (WHAT .NE. 2 .AND. WHAT .NE. 3 .AND.
     $       WHAT.NE.6.AND. WHAT.NE.8 .AND.WHAT.NE.9.AND.
     $       WHAT.NE.17) THEN
          WRITE(*,*)
     *  "Internal error 1 in DMUMPS_460",WHAT
        END IF
        MYID2(1) = MYID
        NDEST = NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#endif
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        CALL MPI_PACK_SIZE( 1 + (NDEST-1) * OVHSIZE, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        IF((WHAT.NE.17).AND.(WHAT.NE.10))THEN
           NREALS = 1
        ELSE
           NREALS = 2
        ENDIF
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     &                      COMM, SIZE2, IERR )
        SIZE = SIZE1 + SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( LOAD, 1, MPI_DOUBLE_PRECISION,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        IF((WHAT.EQ.17).OR.(WHAT.EQ.10))THEN
           CALL MPI_PACK( UPD_LOAD, 1, MPI_DOUBLE_PRECISION,
     *          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
          IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
          IF ( I .ne. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_460'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_460
        SUBROUTINE DMUMPS_519
     *             ( WHAT, COMM, NPROCS,
     *               FATHER_NODE,INODE,NCB,K81,
     *               MYID,REMOTE, IERR)
        IMPLICIT NONE
        INTEGER COMM, NPROCS, MYID, IERR, WHAT,REMOTE
        INTEGER FATHER_NODE,INODE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER POSITION, IREQ, IPOS, SIZE
        INTEGER I, NDEST, IDEST, IPOSMSG
        INTEGER IZERO,NCB,K81
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1) = MYID
        NDEST = 1
        IF ( NDEST .eq. 0 ) THEN
           RETURN
        ENDIF
        IF((K81.EQ.2).OR.(K81.EQ.3))THEN
           CALL MPI_PACK_SIZE( 4 + OVHSIZE, 
     *          MPI_INTEGER, COMM,
     *          SIZE, IERR )
        ELSE
           CALL MPI_PACK_SIZE( 2 + OVHSIZE, 
     *          MPI_INTEGER, COMM,
     *          SIZE, IERR )
        ENDIF
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR, 
     &                  IZERO, MYID2 
     &               )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        CALL MPI_PACK( FATHER_NODE, 1, MPI_INTEGER,
     *                 BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *                 POSITION, COMM, IERR )
        IF((K81.EQ.2).OR.(K81.EQ.3))THEN
           CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *          POSITION, COMM, IERR )
           CALL MPI_PACK( NCB, 1, MPI_INTEGER,
     *          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *          POSITION, COMM, IERR )
        ENDIF
        IDEST = 1
        CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                 POSITION, MPI_PACKED, REMOTE,
     *                 UPDATE_LOAD, COMM,
     *                 BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                 IERR )
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_519'
          WRITE(*,*) ' Size,position=',size,position
          stop
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_519
        SUBROUTINE DMUMPS_502( COMM, MYID, NPROCS,
     *  MAX_SURF_MASTER,IERR)
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        INTEGER IPOS, IREQ, IDEST, IPOSMSG, POSITION, I
        INTEGER COMM, MYID, IERR, NPROCS,MAX_SURF_MASTER
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        INTEGER NDEST, NINTS, SIZE
        INTEGER WHAT
        IERR = 0
        MYID2(1) = MYID
        NDEST = NPROCS - 1
        NINTS = 1 + 1 + ( NDEST-1 ) * OVHSIZE
        CALL MPI_PACK_SIZE( NINTS, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE, IERR )
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     *       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        WHAT = 4
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( MAX_SURF_MASTER, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        IDEST = 0
        DO I = 0, NPROCS - 1
           IF ( I .ne. MYID ) THEN
              IDEST = IDEST + 1
              CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *             POSITION, MPI_PACKED, I,
     *             UPDATE_LOAD, COMM,
     *             BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *             IERR )
           END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_524'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_502
        SUBROUTINE DMUMPS_524( BDC_MEM,
     *      COMM, MYID, NPROCS,
#if ! defined(OLD_LOAD_MECHANISM)
     *      FUTURE_NIV2,
#endif
     *      NSLAVES,
     *      LIST_SLAVES,INODE,
     *      MEM_INCREMENT, FLOPS_INCREMENT,CB_BAND, WHAT,
     *      IERR )
        IMPLICIT NONE
        INCLUDE 'mpif.h'
        INCLUDE 'mumps_tags.h'
        LOGICAL BDC_MEM
        INTEGER COMM, MYID, NPROCS, NSLAVES, IERR
#if ! defined(OLD_LOAD_MECHANISM)
        INTEGER FUTURE_NIV2(NPROCS)
#endif
        INTEGER LIST_SLAVES(NSLAVES),INODE
        DOUBLE PRECISION MEM_INCREMENT(NSLAVES)
        DOUBLE PRECISION FLOPS_INCREMENT(NSLAVES)
        DOUBLE PRECISION CB_BAND(NSLAVES)
        INTEGER NDEST, NINTS, NREALS, SIZE1, SIZE2, SIZE
        INTEGER IPOS, IPOSMSG, IREQ, POSITION
        INTEGER I, IDEST, WHAT
        INTEGER IZERO
        INTEGER MYID2(1)
        PARAMETER ( IZERO=0 )
        MYID2(1)=MYID
        IERR = 0
#if ! defined(OLD_LOAD_MECHANISM)
        NDEST = 0
        DO I = 1, NPROCS
          IF ( I .NE. MYID + 1 .AND. FUTURE_NIV2(I).NE.0) THEN
            NDEST = NDEST + 1
          ENDIF
        ENDDO
#else
        NDEST = NPROCS - 1
#endif
        IF ( NDEST == 0 ) THEN
           RETURN
        ENDIF
        NINTS = 2 +  NSLAVES + ( NDEST - 1 ) * OVHSIZE + 1
        NREALS = NSLAVES
        IF (BDC_MEM) NREALS = NREALS + NSLAVES
        IF(WHAT.EQ.19) THEN 
           NREALS = NREALS + NSLAVES
        ENDIF
        CALL MPI_PACK_SIZE( NINTS, 
     *                       MPI_INTEGER, COMM,
     *                       SIZE1, IERR )
        CALL MPI_PACK_SIZE( NREALS, MPI_DOUBLE_PRECISION,
     *       COMM, SIZE2, IERR )
        SIZE = SIZE1+SIZE2
        CALL DMUMPS_4( BUF_LOAD, IPOS, IREQ, SIZE, IERR,
     *       IZERO, MYID2 )
        IF ( IERR .LT. 0 ) THEN
           RETURN
        ENDIF
        BUF_LOAD%ILASTMSG = BUF_LOAD%ILASTMSG + ( NDEST - 1 ) * OVHSIZE
        IPOS = IPOS - OVHSIZE
        DO IDEST = 1, NDEST - 1
          BUF_LOAD%CONTENT( IPOS + ( IDEST - 1 ) * OVHSIZE ) =
     *    IPOS + IDEST * OVHSIZE
        END DO
        BUF_LOAD%CONTENT( IPOS + ( NDEST - 1 ) * OVHSIZE ) = 0
        IPOSMSG = IPOS + OVHSIZE * NDEST
        POSITION = 0
        CALL MPI_PACK( WHAT, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( NSLAVES, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( INODE, 1, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( LIST_SLAVES, NSLAVES, MPI_INTEGER,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        CALL MPI_PACK( FLOPS_INCREMENT, NSLAVES,
     *      MPI_DOUBLE_PRECISION,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        IF (BDC_MEM) THEN
          CALL MPI_PACK( MEM_INCREMENT, NSLAVES,
     *      MPI_DOUBLE_PRECISION,
     *      BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *      POSITION, COMM, IERR )
        END IF
        IF(WHAT.EQ.19)THEN
           CALL MPI_PACK( CB_BAND, NSLAVES,
     *          MPI_DOUBLE_PRECISION,
     *          BUF_LOAD%CONTENT( IPOSMSG ), SIZE,
     *          POSITION, COMM, IERR )
        ENDIF
        IDEST = 0
        DO I = 0, NPROCS - 1
#if ! defined(OLD_LOAD_MECHANISM)
        IF ( I .NE. MYID .AND. FUTURE_NIV2(I+1) .NE. 0) THEN
#else
        IF ( I .NE. MYID ) THEN
#endif
            IDEST = IDEST + 1
            CALL MPI_ISEND( BUF_LOAD%CONTENT( IPOSMSG ),
     *                      POSITION, MPI_PACKED, I,
     *                      UPDATE_LOAD, COMM,
     *                      BUF_LOAD%CONTENT( IREQ+(IDEST-1)*OVHSIZE ),
     *                      IERR )
          END IF
        END DO
        SIZE = SIZE - ( NDEST - 1 ) * OVHSIZE * SIZEofINT
        IF ( SIZE .LT. POSITION ) THEN
          WRITE(*,*) ' Error in DMUMPS_524'
          WRITE(*,*) ' Size,position=',size,position
          CALL DMUMPS_ABORT()
        END IF
        IF ( SIZE .NE. POSITION )
     *  CALL DMUMPS_1( BUF_LOAD, POSITION )
        RETURN
        END SUBROUTINE DMUMPS_524
        SUBROUTINE DMUMPS_60
     *             ( DMUMPS_LBUFR_BYTES)
        IMPLICIT NONE
        INTEGER DMUMPS_LBUFR_BYTES 
        SIZE_RBUF_BYTES = DMUMPS_LBUFR_BYTES
        RETURN
      END SUBROUTINE DMUMPS_60
      END MODULE DMUMPS_BUFFER
