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
C
C  This file contains dummy MPI/BLACS/ScaLAPACK libraries to allow
C  linking/running MUMPS on a platform where MPI is not installed.
C
C******************************************************************
      SUBROUTINE MPI_GATHER( SENDBUF, COUNT, 
     *         DATATYPE, RECVBUF, RECCOUNT, RECTYPE,
     *         ROOT, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, RECCOUNT, RECTYPE,
     *        ROOT, COMM, IERR
      INTEGER SENDBUF(* ), RECVBUF( * )
      IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
        CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
        CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION ) THEN
        CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
        CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
        CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
        CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
        CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE
        WRITE(*,*) 'ERROR in MPI_GATHER=',DATATYPE
        STOP
      END IF
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GATHERV( SENDBUF, COUNT, 
     *         DATATYPE, RECVBUF, RECCOUNT, DISPLS, RECTYPE,
     *         ROOT, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, RECTYPE,
     *        ROOT, COMM, IERR
      INTEGER RECCOUNT( * )
      INTEGER SENDBUF( * ), RECVBUF( * )
      INTEGER DISPLS( * )
C
C     Note that DISPLS is ignored in this version. One may
C     want to copy in reception buffer with a shift DISPLS(1).
C     This requires passing the offset DISPLS(1) to
C     "MUMPS_COPY_DATATYPE" routines.
C
      IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
        CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
        CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION ) THEN
        CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
        CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
        CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
        CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
        CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE
        WRITE(*,*) 'ERROR in MPI_GATHERV=',DATATYPE
        STOP
      END IF
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_ALLREDUCE( SENDBUF, RECVBUF, COUNT, DATATYPE,
     *                          OPERATION, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, OPERATION, COMM, IERR
      INTEGER SENDBUF( * ), RECVBUF( * )
      IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
        CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
        CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION ) THEN
        CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
        CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
        CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
        CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
        CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE
        WRITE(*,*) 'ERROR in MPI_ALLREDUCE. DATATYPE=',DATATYPE
        STOP
      END IF
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_REDUCE( SENDBUF, RECVBUF, COUNT, DATATYPE, OP,
     *           ROOT, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, OP, ROOT, COMM, IERR
      INTEGER SENDBUF(*), RECVBUF(*)
      IF ( DATATYPE .EQ. MPI_INTEGER ) THEN
        CALL MUMPS_COPY_INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_REAL ) THEN
        CALL MUMPS_COPY_REAL( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION ) THEN
        CALL MUMPS_COPY_DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
        CALL MUMPS_COPY_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_COMPLEX ) THEN
        CALL MUMPS_COPY_DOUBLE_COMPLEX( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2DOUBLE_PRECISION) THEN
        CALL MUMPS_COPY_2DOUBLE_PRECISION( SENDBUF, RECVBUF, COUNT )
      ELSE IF ( DATATYPE .EQ. MPI_2INTEGER) THEN
        CALL MUMPS_COPY_2INTEGER( SENDBUF, RECVBUF, COUNT )
      ELSE
        WRITE(*,*) 'ERROR in MPI_REDUCE. DATATYPE=',DATATYPE
        STOP
      END IF
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_BARRIER( COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GET_PROCESSOR_NAME( NAME, RESULTLEN, IERROR)
      CHARACTER (LEN=*) NAME
      INTEGER RESULTLEN,IERROR
      RESULTLEN = 1
      IERROR = 0
      NAME = 'X'
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_BCAST( BUFFER, COUNT, DATATYPE, ROOT, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, ROOT, COMM, IERR
      INTEGER BUFFER( * )
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_CANCEL( IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_CREATE( COMM, GROUP, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, GROUP, COMM2, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_DUP( COMM, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, COMM2, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_FREE( COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_GROUP( COMM, GROUP, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, GROUP, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_RANK( COMM, RANK, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, RANK, IERR
      RANK = 0
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_SIZE( COMM, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, SIZE, IERR
      SIZE = 1
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_COMM_SPLIT( COMM, COLOR, KEY, COMM2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, COLOR, KEY, COMM2, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
c     SUBROUTINE MPI_ERRHANDLER_SET( COMM, ERRHANDLER, IERR )
c     IMPLICIT NONE
c     INCLUDE 'mpif.h'
c     INTEGER COMM, ERRHANDLER, IERR
c     IERR = 0
c     RETURN
c     END
C***********************************************************************
      SUBROUTINE MPI_FINALIZE( IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GET_COUNT( STATUS, DATATYPE, COUNT, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER DATATYPE, COUNT, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_GET_COUNT should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GROUP_FREE( GROUP, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GROUP_RANGE_EXCL( GROUP, N, RANGES, GROUP2, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, N, GROUP2, IERR
      INTEGER RANGES(*)
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_GROUP_SIZE( GROUP, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER GROUP, SIZE, IERR
      SIZE = 1 ! Or should it be zero ?
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_INIT(IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_INITIALIZED( FLAG, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL FLAG
      INTEGER IERR
      FLAG = .TRUE.
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_IPROBE( SOURCE, TAG, COMM, FLAG, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER SOURCE, TAG, COMM, IERR
      INTEGER STATUS(MPI_STATUS_SIZE)
      LOGICAL FLAG
      FLAG = .FALSE.
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_IRECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM,
     *           IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, IREQ, IERR
      INTEGER BUF(*)
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_ISEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM,
     *           IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR, IREQ
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_ISEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_PACK( INBUF, INCOUNT, DATATYPE, OUTBUF, OUTCOUNT,
     *           POSITION, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INCOUNT, DATATYPE, OUTCOUNT, POSITION, COMM, IERR
      INTEGER INBUF(*), OUTBUF(*)
      WRITE(*,*) 'Error. MPI_PACKED should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_PACK_SIZE( INCOUNT, DATATYPE, COMM, SIZE, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INCOUNT, DATATYPE, COMM, SIZE, IERR
      WRITE(*,*) 'Error. MPI_PACK_SIZE should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_PROBE( SOURCE, TAG, COMM, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER SOURCE, TAG, COMM, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_PROBE should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_RECV( BUF, COUNT, DATATYPE, SOURCE, TAG, COMM,
     *           STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, SOURCE, TAG, COMM, IERR
      INTEGER BUF(*), STATUS(MPI_STATUS_SIZE)
      WRITE(*,*) 'Error. MPI_RECV should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_REQUEST_FREE( IREQ, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_SEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_SEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_SSEND( BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, DATATYPE, DEST, TAG, COMM, IERR
      INTEGER BUF(*)
      WRITE(*,*) 'Error. MPI_SSEND should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_TEST( IREQ, FLAG, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL FLAG
      FLAG = .FALSE.
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_UNPACK( INBUF, INSIZE, POSITION, OUTBUF, OUTCOUNT,
     *           DATATYPE, COMM, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER INSIZE, POSITION, OUTCOUNT, DATATYPE, COMM, IERR
      INTEGER INBUF(*), OUTBUF(*)
      WRITE(*,*) 'Error. MPI_UNPACK should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_WAIT( IREQ, STATUS, IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IREQ, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      WRITE(*,*) 'Error. MPI_WAIT should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_WAITANY( COUNT, ARRAY_OF_REQUESTS, INDEX, STATUS,
     *           IERR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COUNT, INDEX, IERR
      INTEGER STATUS( MPI_STATUS_SIZE )
      INTEGER ARRAY_OF_REQUESTS( * )
      WRITE(*,*) 'Error. MPI_WAITANY should not be called.'
      STOP
      IERR = 0
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION MPI_WTIME( )
C     elapsed time
      DOUBLE PRECISION VAL
C     write(*,*) 'Entering MPI_WTIME'
      CALL ELAPSE( VAL )
      MPI_WTIME = VAL
C     write(*,*) 'Exiting MPI_WTIME'
      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_ATTR_PUT( COMM, KEY, VAL, IERR )
      IMPLICIT NONE
      INTEGER COMM, KEY, VAL, IERR

      RETURN
      END
C***********************************************************************
      SUBROUTINE MPI_ABORT( COMM, IERRCODE, IERR )
      IMPLICIT NONE
      INTEGER COMM, IERRCODE, IERR
      WRITE(*,*) "** Abort called"
      STOP
      END SUBROUTINE MPI_ABORT
C***********************************************************************
C     BLACS and ScaLAPACK
C
      SUBROUTINE BLACS_GRIDINIT( CNTXT, C, NPROW, NPCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL
      CHARACTER C
        WRITE(*,*) 'Error. BLACS_GRIDINIT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE BLACS_GRIDINFO( CNTXT, NPROW, NPCOL, MYROW, MYCOL )
      IMPLICIT NONE
      INTEGER CNTXT, NPROW, NPCOL, MYROW, MYCOL
        WRITE(*,*) 'Error. BLACS_GRIDINFO should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE BLACS_GRIDEXIT( CNTXT )
      IMPLICIT NONE
      INTEGER CNTXT
        WRITE(*,*) 'Error. BLACS_GRIDEXIT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC,
     *           ICTXT, LLD, INFO )
      IMPLICIT NONE
      INTEGER ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER DESC( * )
        WRITE(*,*) 'Error. DESCINIT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS ) 
      INTEGER N, NB, IPROC, ISRCPROC, NPROCS
C     Can be called
      IF ( NPROCS .ne. 1 ) THEN
        WRITE(*,*) 'Error. Last parameter from NUMROC should be 1'
        STOP
      ENDIF
      IF ( IPROC .ne. 0 ) THEN
        WRITE(*,*) 'Error. IPROC should be 0 in NUMROC.'
        STOP
      ENDIF
      NUMROC = N
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCPOTRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX            A( * )
        WRITE(*,*) 'Error. PCGETRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     *                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      DOUBLE COMPLEX     A( * )
        WRITE(*,*) 'Error. PZPOTRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE COMPLEX     A( * )
        WRITE(*,*) 'Error. PZGETRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     *                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE COMPLEX     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSPOTRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL               A( * )
        WRITE(*,*) 'Error. PSGETRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     *                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDPOTRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      DOUBLE PRECISION   A( * )
        WRITE(*,*) 'Error. PDGETRF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     *                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW,
     *                    MYCOL, LRINDX, LCINDX, RSRC, CSRC )
      IMPLICIT NONE
      INTEGER            CSRC, GCINDX, GRINDX, LRINDX, LCINDX, MYCOL,
     *                   MYROW, NPCOL, NPROW, RSRC
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
      INTEGER            INDXGLOB, IPROC, ISRCPROC, NB, NPROCS
        INDXG2P = 0
        WRITE(*,*) 'Error. INFOG2L should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCSCAL(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      COMPLEX            ALPHA
      COMPLEX            X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PCSCAL should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZSCAL(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      DOUBLE COMPLEX     ALPHA
      DOUBLE COMPLEX     X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PZSCAL should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDSCAL(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      DOUBLE PRECISION   ALPHA
      DOUBLE PRECISION   X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PDSCAL should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSSCAL(N, ALPHA, X, IX, JX, DESCX, INCX)
      IMPLICIT NONE
      INTEGER            INCX, N, IX, JX
      REAL               ALPHA
      REAL               X( * )
      INTEGER            DESCX( * )
        WRITE(*,*) 'Error. PSSCAL should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PCDOT( CONTXT, N, X, Y )
      INTEGER CONTXT, N
      COMPLEX X(*), Y(*)
        PCDOT = (0.0e0,0.0e0)
        WRITE(*,*) 'Error. PCDOT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PZDOT( CONTXT, N, X, Y )
      INTEGER CONTXT, N
      DOUBLE COMPLEX X(*), Y(*)
        PZDOT = (0.0d0,0.0d0)
        WRITE(*,*) 'Error. PZDOT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PDDOT( CONTXT, N, X, Y )
      INTEGER CONTXT, N
      DOUBLE PRECISION X(*), Y(*)
        PDDOT = 0.0d0
        WRITE(*,*) 'Error. PDDOT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PSDOT( CONTXT, N, X, Y )
      INTEGER CONTXT, N
      REAL X(*), Y(*)
        PSDOT = 0.0e0
        WRITE(*,*) 'Error. PSDOT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE ZGEBS2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBS2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE CGEBS2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBS2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE SGEBS2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBS2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE DGEBS2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBS2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE ZGEBR2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. ZGEBR2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE CGEBR2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      COMPLEX A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. CGEBR2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE SGEBR2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      REAL A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. SGEBR2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE DGEBR2D( CONTXT, SCOPE, TOP, M, N, A, LDA )
      IMPLICIT NONE
      INTEGER CONTXT, M, N, LDA
      DOUBLE PRECISION A(*)
      CHARACTER SCOPE, TOP
        WRITE(*,*) 'Error. DGEBR2D should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     *                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCGETRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     *                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE COMPLEX     A( * ), B( * )
        WRITE(*,*) 'Error. PZGETRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     *                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSGETRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,
     *                    IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDGETRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     *           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      COMPLEX         A( * ), B( * )
        WRITE(*,*) 'Error. PCPOTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     *           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      DOUBLE COMPLEX  A( * ), B( * )
        WRITE(*,*) 'Error. PZPOTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     *           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      REAL            A( * ), B( * )
        WRITE(*,*) 'Error. PSPOTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
     *           DESCB, INFO )
      IMPLICIT NONE
      CHARACTER       UPLO
      INTEGER         IA, IB, INFO, JA, JB, N, NRHS
      INTEGER         DESCA( * ), DESCB( * )
      DOUBLE          PRECISION A( * ), B( * )
        WRITE(*,*) 'Error. PDPOTRS should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSCNRM2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL NORM2
      COMPLEX X( * )
        WRITE(*,*) 'Error. PCNRM2 should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDZNRM2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2
      DOUBLE COMPLEX  X( * )
        WRITE(*,*) 'Error. PZNRM2 should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSNRM2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      REAL    NORM2, X( * )
        WRITE(*,*) 'Error. PSNRM2 should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDNRM2( N, NORM2, X, IX, JX, DESCX, INCX )
      IMPLICIT NONE
      INTEGER N, IX, JX, INCX
      INTEGER DESCX(*)
      DOUBLE PRECISION NORM2, X( * )
        WRITE(*,*) 'Error. PDNRM2 should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      REAL FUNCTION PCLANGE( NORM, M, N, A, IA,  JA,
     *                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      COMPLEX      A( * ), WORK( * )
      PCLANGE = 0.0e0
        WRITE(*,*) 'Error. PCLANGE should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PZLANGE( NORM, M, N, A, IA,  JA,
     *                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PZLANGE = 0.0d0
        WRITE(*,*) 'Error. PZLANGE should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      REAL FUNCTION PSLANGE( NORM, M, N, A, IA,  JA,
     *                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      REAL         A( * ), WORK( * )
      PSLANGE = 0.0e0
        WRITE(*,*) 'Error. PSLANGE should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      DOUBLE PRECISION FUNCTION PDLANGE( NORM, M, N, A, IA,  JA,
     *                 DESCA, WORK )
      CHARACTER    NORM
      INTEGER      IA, JA, M, N
      INTEGER      DESCA( * )
      DOUBLE       PRECISION A( * ), WORK( * )
      PDLANGE = 0.0d0
        WRITE(*,*) 'Error. PDLANGE should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCGECON( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     *           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      COMPLEX         A( * ), WORK( * )
        WRITE(*,*) 'Error. PCGECON should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZGECON( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     *           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      DOUBLE COMPLEX  A( * ), WORK( * )
        WRITE(*,*) 'Error. PZGECON should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSGECON( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     *           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      REAL            ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      REAL            A( * ), WORK( * )
        WRITE(*,*) 'Error. PSGECON should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDGECON( NORM, N,  A,  IA,  JA,  DESCA,  ANORM,
     *           RCOND,  WORK,  LWORK,  IWORK,  LIWORK, INFO )
      IMPLICIT NONE

      CHARACTER       NORM
      INTEGER         IA, INFO, JA, LIWORK, LWORK, N
      DOUBLE          PRECISION ANORM, RCOND
      INTEGER         DESCA( * ), IWORK( * )
      DOUBLE          PRECISION A( * ), WORK( * )
        WRITE(*,*) 'Error. PDGECON should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCGEQPF( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     *           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
      COMPLEX    A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCGEQPF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZGEQPF( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     *           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER    IA, JA, INFO, LWORK, M, N
      INTEGER    DESCA( * ), IPIV( * )
      DOUBLE COMPLEX A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZGEQPF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSGEQPF( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     *           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      REAL       A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSGEQPF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDGEQPF( M,  N,  A,  IA,  JA, DESCA, IPIV, TAU,
     *           WORK, LWORK, INFO )
      IMPLICIT NONE
      INTEGER         IA, JA, INFO, LWORK, M, N
      INTEGER         DESCA( * ), IPIV( * )
      DOUBLE PRECISION A( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDGEQPF should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCAXPY(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     *           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      COMPLEX A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PCAXPY should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZAXPY(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     *           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE COMPLEX A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PZAXPY should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSAXPY(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     *           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      REAL A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PSAXPY should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDAXPY(N, A, X, IX, JX, DESCX, INCX, Y, IY, JY,
     *           DESCY, INCY)
      IMPLICIT NONE
      INTEGER N, IX, IY, JX, JY, INCX, INCY
      INTEGER DESCX(*), DESCY(*)
      DOUBLE PRECISION A(*),X(*),Y(*)
        WRITE(*,*) 'Error. PDAXPY should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      COMPLEX            ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX            A( * ), B( * )
        WRITE(*,*) 'Error. PCTRSM should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      DOUBLE COMPLEX     ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE COMPLEX     A( * ), B( * )
        WRITE(*,*) 'Error. PZTRSM should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      REAL               ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      REAL               A( * ), B( * )
        WRITE(*,*) 'Error. PSTRSM should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, IA,
     $                   JA, DESCA, B, IB, JB, DESCB )
      IMPLICIT NONE
      CHARACTER          SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, IA, JA, IB, JB
      DOUBLE PRECISION   ALPHA
      INTEGER            DESCA( * ), DESCB( * )
      DOUBLE PRECISION   A( * ), B( * )
        WRITE(*,*) 'Error. PDTRSM should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCUNMQR( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     *                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     *                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      COMPLEX   A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PCUNMQR should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PZUNMQR( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     *                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     *                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      DOUBLE COMPLEX A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PZUNMQR should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PSORMQR( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     *                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     *                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER   IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER   DESCA( * ), DESCC( * )
      REAL      A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PSORMQR should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PDORMQR( SIDE,  TRANS,  M,  N,  K,  A,  IA, JA,
     *                    DESCA, TAU, C, IC,  JC,  DESCC,  WORK,
     *                    LWORK, INFO )
      IMPLICIT NONE
      CHARACTER SIDE, TRANS
      INTEGER         IA, IC, INFO, JA, JC, K, LWORK, M, N
      INTEGER         DESCA( * ), DESCC( * )
      DOUBLE PRECISION  A(  *  ), C( * ), TAU( * ), WORK( * )
        WRITE(*,*) 'Error. PDORMQR should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE CHK1MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     *                    DESCAPOS0, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA, NAPOS0
      INTEGER            DESCA( * )
        WRITE(*,*) 'Error. CHK1MAT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PCHK2MAT( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     *                     DESCAPOS0, MB, MBPOS0, NB, NBPOS0, IB, JB,
     *                     DESCB, DESCBPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, DESCBPOS0, IA, IB, INFO, JA, JB, MA,
     *                   MAPOS0, MB, MBPOS0, NA, NAPOS0, NB, NBPOS0,
     *                   NEXTRA
      INTEGER            DESCA( * ), DESCB( * ), EX( NEXTRA ),
     *                   EXPOS( NEXTRA )
        WRITE(*,*) 'Error. PCHK2MAT should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE PXERBLA( CONTXT, SRNAME, INFO )
      IMPLICIT NONE
      INTEGER CONTXT, INFO
      CHARACTER SRNAME
        WRITE(*,*) 'Error. PXERBLA should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
      SUBROUTINE DESCSET( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     *                    LLD )
      IMPLICIT NONE
      INTEGER            ICSRC, ICTXT, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
        WRITE(*,*) 'Error. DESCSET should not be called.'
        STOP
      RETURN
      END
C***********************************************************************
C***********************************************************************
C
C  Utilities to copy data
C
      SUBROUTINE MUMPS_COPY_INTEGER( S, R, N )
      IMPLICIT NONE
      INTEGER N
      INTEGER S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_2INTEGER( S, R, N )
      IMPLICIT NONE
      INTEGER N
      INTEGER S(N+N),R(N+N)
      INTEGER I
      DO I = 1, N+N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_REAL( S, R, N )
      IMPLICIT NONE
      INTEGER N
      REAL S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_2DOUBLE_PRECISION( S, R, N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION S(N+N),R(N+N)
      INTEGER I
      DO I = 1, N+N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_DOUBLE_PRECISION( S, R, N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_COMPLEX( S, R, N )
      IMPLICIT NONE
      INTEGER N
      COMPLEX S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END
      SUBROUTINE MUMPS_COPY_DOUBLE_COMPLEX( S, R, N )
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I) = S(I)
      END DO
      RETURN
      END

