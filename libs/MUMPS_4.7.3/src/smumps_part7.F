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
#if defined(try_null_space)
      SUBROUTINE SMUMPS_295
     *(NN,A,ROOT,WR03,LWR03,KEEP,KEEP8,INFO,LP,DKEEP,
     * GLOBK109,OPELIW,NNEG)
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER NN,LWR03,LP
      REAL A(NN*NN)
      INTEGER INFO(2),KEEP(500)
      INTEGER*8 KEEP8(150)
      REAL DKEEP(30)
      REAL WR03(LWR03)
      DOUBLE PRECISION OPELIW
      INTEGER NNEG
      INTEGER GLOBK109
      INTEGER LDLT,RRSTRAT,DEFICIENCY
      INTEGER I,J, IND, LDA
      INTEGER IERR
      REAL TEMP,TEMP1
      INTEGER ALLOCOK,MAXDEF,MINDEF
      LOGICAL GAP_FIND_IS_ON
      REAL ALPHA,SUM
      REAL FIXA
      REAL TOL,EPS,GAP_SIZE,AVG
      REAL ZERO
      PARAMETER(ZERO=0.0E0)
      PARAMETER(GAP_SIZE=1.0E2)
      REAL SMUMPS_104
      EXTERNAL SMUMPS_104
      EPS = EPSILON(ZERO)
      LDLT=KEEP(50)
      RRSTRAT=KEEP(19)
      IF (RRSTRAT .NE. 1) THEN
        WRITE(*,*) " *** Option ",RRSTRAT,
     *             " for null space no more available."
        CALL SMUMPS_ABORT()
      ENDIF
      IF (  LWR03 .LT. 3 * NN ) THEN
        INFO(1) = -9
        INFO(2) = 3 * NN - LWR03
        IF ( LP .GT. 0 )
     *  WRITE(LP,*) 'SMUMPS_295: LWR03 too small'
        RETURN
      END IF
      MAXDEF=KEEP(21)
      IF ( MAXDEF .LE. 0 ) THEN
        MAXDEF = NN
      ELSE
        MAXDEF = MAX(MAXDEF - GLOBK109,0)
      ENDIF
      MINDEF = MAX(KEEP(22) - GLOBK109,0)
      MINDEF = MIN(MINDEF,NN)
      MAXDEF = MIN(MAXDEF,NN)
      IF (LDLT.NE.0) THEN
        CALL SMUMPS_297(NN,A)
      END IF
      DEFICIENCY=0
      KEEP(17)=DEFICIENCY
      LDA=NN
        OPELIW = OPELIW +
     1                DBLE(4)*DBLE(NN)*DBLE(NN)*
     2                DBLE(NN)/DBLE(3)
        IF(ASSOCIATED(root%QR_TAU)) DEALLOCATE(root%QR_TAU)
        ALLOCATE(root%QR_TAU(NN),stat=ALLOCOK )
        IF(ALLOCOK.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          IF ( LP .GT. 0 )
     *    WRITE(LP,*) 'SMUMPS_295: Pb allocating QR_TAU'
          RETURN
        END IF
        IF (ASSOCIATED(root%IPIV)) DEALLOCATE(root%IPIV)
        ALLOCATE(root%IPIV(NN), stat=ALLOCOK)
        IF (ALLOCOK.GT.0) THEN
          INFO(1)=-13
          INFO(2)=NN
          IF ( LP .GT. 0 )
     *    WRITE(LP,*) 'SMUMPS_295:Problem allocating IPIV'
          RETURN
        END IF
          root%IPIV = 0
          CALL SGEQPF(NN,NN,A(1),LDA,root%IPIV,
     *      root%QR_TAU(1),WR03,IERR)
          IF(IERR.NE.0) THEN
            WRITE(*,*) ' Problem in SGEQPF : IERR = ', IERR
            CALL SMUMPS_ABORT()
          END IF
          IF(MAXDEF.EQ.0) GO TO 5488
          TOL=1.1*SQRT(REAL(NN))*EPS*SMUMPS_104(NN,A,WR03)
          TEMP=ABS(A(LDA*NN))   
          IF(TEMP.GE.10*TOL) THEN
            GAP_FIND_IS_ON=.FALSE.
          ELSE
            GAP_FIND_IS_ON=.TRUE.
            AVG=0.02*TEMP
            SUM= REAL(ZERO)
            GO TO 132
          END IF
          TEMP1=TEMP  ! Rnn,nn saved in TEMP1 
          DO I=1,MIN(NN-1,23)
            TEMP=ABS(A(LDA*NN-I*(NN+1)))
            IF(ABS(TEMP1).LE.10E-4*TEMP) THEN
              GAP_FIND_IS_ON=.TRUE.
              AVG=0.02*ABS(A(LDA*NN))
              SUM= REAL(ZERO)
              GO TO 132
            END IF
            TEMP1=TEMP  ! save previous diagonal entry for comparison
          END DO
          GAP_FIND_IS_ON=.FALSE.
 132      CONTINUE
          IF(GAP_FIND_IS_ON) THEN
            DO IND=0,NN-1
              TEMP=A(LDA*NN-IND*(NN+1))
              IF(ABS(TEMP).GT.AVG*GAP_SIZE+EPS) GO TO 1487
              SUM=SUM+ABS(TEMP)
              AVG=SUM/(IND+1)
              IF(IND.GE.MAXDEF) GOTO 1487
            END DO
          END IF
 1487     CONTINUE
          DEFICIENCY=IND
          IF(GAP_FIND_IS_ON) GO TO 5488
          DEFICIENCY=0
 1488     CONTINUE
          IND=(NN-DEFICIENCY)*LDA-DEFICIENCY
          ALPHA=A(IND)
          TOL=1.8*SQRT(REAL(NN))*EPS*SMUMPS_104(NN,A,WR03)
          IF( ( ABS(ALPHA).LE.TOL )
     *     .OR. ( MINDEF.EQ.MAXDEF )
     *     ) THEN
            DEFICIENCY=DEFICIENCY+1
            IF(DEFICIENCY.LT.MAXDEF) GO TO 1488
          END IF
 5488     CONTINUE
          IF(DEFICIENCY.GT.0) THEN
            CALL STRTRS('U','N','N',NN-DEFICIENCY,DEFICIENCY,
     *      A,LDA,A(LDA*(NN-DEFICIENCY)+1),LDA,IERR)
            IF ( IERR .NE. 0 ) THEN
              WRITE(*,*) ' Problem in STRTRS : IERR = ', IERR
              CALL SMUMPS_ABORT()
            END IF
          END IF
      KEEP(17)=DEFICIENCY
      RETURN
      END SUBROUTINE SMUMPS_295
      SUBROUTINE SMUMPS_296
     *  (NRHS,NN,A,ROOT,
     *   IBEG_ROOT_DEF, IEND_ROOT_DEF,
     *   RHS,WK,LWK,KEEP,KEEP8,MTYPE,INFO)
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER NN, NRHS, LWK
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF
      REAL A(NN*NN)
      INTEGER INFO(2),KEEP(500)
      INTEGER*8 KEEP8(150)
      REAL WK(LWK),RHS(NN,NRHS)
      INTEGER I,J,K,LDA,LDRHS,IERR,MTYPE_SAVE
      INTEGER LDLT,RRSTRAT,MTYPE,DEFICIENCY
      REAL MINUSONE
      PARAMETER(MINUSONE=-1.0E0)
      LDLT      = KEEP(50)
      RRSTRAT   = KEEP(19)
      DEFICIENCY=KEEP(17)
      MTYPE_SAVE=MTYPE
      IF(MTYPE.EQ.1) THEN
        MTYPE=0
      ELSE
        MTYPE=1
      END IF
      LDA=NN
      LDRHS = NN
      IERR = 0
      IF (RRSTRAT .NE. 1) THEN
        WRITE(*,*) " *** Option ",RRSTRAT,
     *             " for null space no more available."
        CALL SMUMPS_ABORT()
      ENDIF
      IF (KEEP(111).EQ.0) THEN
        IF(MTYPE.EQ.1) THEN
          CALL SORMQR('L','T',NN,NRHS,NN,A,LDA,root%QR_TAU,
     *      RHS,LDRHS,WK(1),LWK,IERR)
          IF(IERR.LT.0) THEN
            WRITE(*,*) 'Error return from SORMQR: IERR=', IERR
            RETURN
          END IF
          CALL STRTRS('U','N','N',NN-DEFICIENCY,NRHS,A,LDA,
     *                 RHS,LDRHS,IERR)
          DO I=1,NRHS
            RHS( NN - DEFICIENCY + 1: NN, I ) = 0.0E0 ! Min norm solution?
          ENDDO
          IF ( IERR .LT. 0 ) THEN
            WRITE(*,*) ' Problem during solve of the root: IERR =',IERR
            RETURN
          END IF
          DO I=1,NRHS
            CALL SMUMPS_333(NN,root%IPIV,RHS(1,I),WK)
          ENDDO
        ELSEIF(MTYPE.EQ.0) THEN
          DO I=1,NRHS
            CALL SMUMPS_334(NN,root%IPIV,RHS(1,I),WK)
          ENDDO
          CALL STRTRS('U','T','N',NN-DEFICIENCY,1,
     *                 A,LDA,RHS,LDA,IERR)
          DO I=1,NRHS
          RHS( NN - DEFICIENCY + 1: NN, I ) = 0.0E0 ! 1 sol among many
          ENDDO
          IF(IERR.LT.0) THEN
              WRITE(*,*) 'Error return from STRTRS: IERR=', IERR
              RETURN
          END IF
          CALL SORMQR( 'L','N',NN,NRHS,NN,A,LDA,root%QR_TAU,
     *      RHS,LDRHS,WK,LWK,IERR)
          IF(IERR.LT.0) THEN
              WRITE(*,*) 'Error return from SORMQR: IERR=', IERR
              RETURN
          END IF
        ENDIF
      ELSE
        DO I=IBEG_ROOT_DEF,IEND_ROOT_DEF
          RHS(NN-KEEP(17)+I,I-IBEG_ROOT_DEF+1)=REAL(MINUSONE)
          DO K=1,NN-KEEP(17)
            RHS(K,I-IBEG_ROOT_DEF+1)=A(K + LDA*(NN-KEEP(17)+I-1))
          ENDDO
        ENDDO
      ENDIF
      MTYPE=MTYPE_SAVE
      RETURN
      END SUBROUTINE SMUMPS_296
      SUBROUTINE SMUMPS_297(N,A)
      INTEGER N
      REAL A( N, N )
      INTEGER I,J
      DO I = 2, N
        DO J = 1, I - 1
          A( I, J ) = A( J, I )
        END DO
      END DO
      RETURN
      END SUBROUTINE SMUMPS_297
#endif
      SUBROUTINE SMUMPS_635(N,KEEP,ICNTL,MPG)
      IMPLICIT NONE
          INTEGER N, KEEP(500), ICNTL(40), MPG
          KEEP(19)=0
#if defined(try_null_space)
          KEEP(19)=ICNTL(16)
          IF (KEEP(19).NE.1) KEEP(19)=0
          IF ( KEEP(53) .LE. 0 .and.
     *         KEEP(19) .NE. 0 ) THEN
             KEEP(19) = 0
             IF ( MPG .GT. 0 ) THEN
      WRITE( MPG,'(A)') '** Warning:  ICNTL(16) null space option'
      WRITE( MPG,'(A)') '** disabled (incompatibility with analysis)'
             END IF
          END IF
          KEEP(21) = MIN(ICNTL(17),N)
          KEEP(22) = MAX(ICNTL(15),0)
          IF ( KEEP(19) .ne. 0 .and. KEEP(60) .ne. 0 ) THEN
             IF ( MPG .GT. 0 ) THEN
      WRITE( MPG,'(A)') '** Warning: ICNTL(16)  null space option'
      WRITE( MPG,'(A)') '** disabled (incompatibility with Schur)'
             END IF
            KEEP(19) = 0
          END IF
#endif
          RETURN
      END SUBROUTINE SMUMPS_635
      SUBROUTINE SMUMPS_634(KEEP,ICNTL,MPG)
      IMPLICIT NONE
      INTEGER ICNTL(40), KEEP(500), MPG
      KEEP(111)=ICNTL(25)
      IF (KEEP(19).EQ.0.AND.KEEP(110).EQ.0) THEN
        IF(KEEP(111).NE.0.AND.MPG.GT.0) THEN
          WRITE( MPG,'(A)')
     *'** Warning: ICNTL(25) option disabled because'
          WRITE( MPG,'(A)')
     *'** null space was not required during factorization'
        ENDIF
        KEEP(111)=0
      ENDIF
      IF (ICNTL(9).NE.1) THEN
        IF (KEEP(111).NE.0.AND. MPG.GT.0) THEN
          WRITE(MPG,'(A)')
     *'** Warning: ICNTL(25) option disabled because'
          WRITE( MPG,'(A)')
     *'** it is not available for the transposed system'
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_634
      SUBROUTINE SMUMPS_637(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) id
      NULLIFY(id%root%QR_TAU)
#if defined(null_space_old)
      NULLIFY(id%root%GROW)
      NULLIFY(id%root%GSIN)
      NULLIFY(id%root%GCOS)
      NULLIFY(id%root%PERM_ROW)
      NULLIFY(id%root%PERM_COL)
      NULLIFY(id%root%ELROW)
      NULLIFY(id%root%EUROW)
      NULLIFY(id%root%PTREL)
      NULLIFY(id%root%PTREU)
      NULLIFY(id%root%ELELG)
      NULLIFY(id%root%EUELG)
      NULLIFY(id%root%DL)
#endif
      RETURN
      END SUBROUTINE SMUMPS_637
      SUBROUTINE SMUMPS_636(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) id
      IF (ASSOCIATED(id%root%QR_TAU))  THEN
        DEALLOCATE(id%root%QR_TAU)
        NULLIFY(id%root%QR_TAU)
      ENDIF
#if defined(null_space_old)
      IF (ASSOCIATED(id%root%GROW))    THEN
        DEALLOCATE(id%root%GROW)
        NULLIFY(id%root%GROW)
      ENDIF
      IF (ASSOCIATED(id%root%GSIN))    THEN
        DEALLOCATE(id%root%GSIN)
        NULLIFY(id%root%GSIN)
      ENDIF
      IF (ASSOCIATED(id%root%GCOS))    THEN
        DEALLOCATE(id%root%GCOS)
        NULLIFY(id%root%GCOS)
      ENDIF
      IF (ASSOCIATED(id%root%PERM_ROW))THEN
        DEALLOCATE(id%root%PERM_ROW)
        NULLIFY(id%root%PERM_ROW)
      ENDIF
      IF (ASSOCIATED(id%root%PERM_COL))THEN
        DEALLOCATE(id%root%PERM_COL)
        NULLIFY(id%root%PERM_COL)
      ENDIF
      IF (ASSOCIATED(id%root%ELROW))   THEN
        DEALLOCATE(id%root%ELROW)
        NULLIFY(id%root%ELROW)
      ENDIF
      IF (ASSOCIATED(id%root%EUROW))   THEN
        DEALLOCATE(id%root%EUROW)
        NULLIFY(id%root%EUROW)
      ENDIF
      IF (ASSOCIATED(id%root%PTREL))   THEN
        DEALLOCATE(id%root%PTREL)
        NULLIFY(id%root%PTREL)
      ENDIF
      IF (ASSOCIATED(id%root%PTREU))   THEN
        DEALLOCATE(id%root%PTREU)
        NULLIFY(id%root%PTREU)
      ENDIF
      IF (ASSOCIATED(id%root%ELELG))   THEN
        DEALLOCATE(id%root%ELELG)
        NULLIFY(id%root%ELELG)
      ENDIF
      IF (ASSOCIATED(id%root%EUELG))   THEN
        DEALLOCATE(id%root%EUELG)
        NULLIFY(id%root%EUELG)
      ENDIF
      IF (ASSOCIATED(id%root%DL))      THEN
        DEALLOCATE(id%root%DL)
        NULLIFY(id%root%DL)
      ENDIF
#endif
      RETURN
      END SUBROUTINE SMUMPS_636
      SUBROUTINE SMUMPS_279( PHASE, MBLOCK, NBLOCK, 
     *           SIZE_ROOT,
     *           LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     *           LIWK_RR, LWK_RR )
      INTEGER, INTENT(IN) :: PHASE, SIZE_ROOT
      INTEGER, INTENT(IN) :: MBLOCK, NBLOCK, LOCAL_M, LOCAL_N
      LOGICAL, INTENT(IN) :: ROOT_OWNER
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER*8, INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(OUT):: LIWK_RR, LWK_RR
      INTEGER NBPOSPONED_ESTIM
      PARAMETER (NBPOSPONED_ESTIM=2000)
      INTEGER QR,PAR_ROOT
      QR      =KEEP(19)
      PAR_ROOT=KEEP(38)
      LIWK_RR = 0
      LWK_RR = 0
      IF (PAR_ROOT.EQ.0) THEN
       IF(ROOT_OWNER) THEN
          IF(QR.EQ.1) THEN
           IF (PHASE.EQ.0) THEN
            LWK_RR=3*(SIZE_ROOT+NBPOSPONED_ESTIM)
           ELSE
            LWK_RR=3*SIZE_ROOT
           ENDIF
          END IF
       END IF
      ENDIF
#if defined(null_space_old)
      IF(PAR_ROOT.NE.0) THEN
        IF(QR.GE.1.AND.QR.LE.4) THEN
          LWK_RR = MAX(5*MBLOCK*NBLOCK+
     *     5*MAX(LOCAL_M,LOCAL_N)+6+SIZE_ROOT,4*SIZE_ROOT+5)
        ELSE
          LWK_RR = MAX(5*MBLOCK*NBLOCK+
     *     5*MAX(LOCAL_M,LOCAL_N)+6+SIZE_ROOT,4*SIZE_ROOT+5)
        END IF
      ELSE
        IF(ROOT_OWNER) THEN
          IF(QR.GE.1.AND.QR.LE.4) THEN
            LWK_RR=10000+5*SIZE_ROOT
          ELSE
            LWK_RR=10000+5*SIZE_ROOT
          END IF
        ELSE
          LWK_RR=SIZE_ROOT+1
        END IF
      END IF
#endif
      RETURN
      END
      SUBROUTINE SMUMPS_333(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      REAL RN01(N),X(N)
      DO 100 I=1,N
      RN01(PERM(I))=X(I)
100   CONTINUE
      DO 200 I=1,N
      X(I)=RN01(I)
200   CONTINUE
      RETURN
      END
      SUBROUTINE SMUMPS_334(N,PERM,X,RN01)
      INTEGER N,PERM(N),I
      REAL RN01(N),X(N)
      DO 100 I=1,N
      RN01(I)=X(PERM(I))
100   CONTINUE
      DO 200 I=1,N
      X(I)=RN01(I)
200   CONTINUE
      RETURN
      END
      REAL FUNCTION SMUMPS_104(N,AA,WR01)
      INTEGER N
      REAL AA(N*N),WR01(N)
      INTEGER I,J,JSTRT
      REAL ZERO
      REAL VALUE
      PARAMETER(ZERO=0.0E0)
      VALUE=ZERO
      CALL SMUMPS_117(N,ZERO,WR01,1)
      JSTRT=1
      DO I=1,N
        DO J=JSTRT,JSTRT+N-1
          WR01(I)=WR01(I)+ABS(AA(J))
        END DO
        JSTRT=JSTRT+N
      END DO
      DO I=1,N
        VALUE=MAX(ABS(VALUE),ABS(WR01(I)))
      END DO
      SMUMPS_104=VALUE
      RETURN
      END
      SUBROUTINE SMUMPS_117(N,ALPHA,DX,INCX)
      INTEGER N,INCX
      REAL ALPHA,DX(*)
      INTEGER I,IX,M,MP1
      INTRINSIC MOD
      IF(N.LE.0) RETURN
      IF(INCX.NE.1) THEN
        IX=1
        IF(INCX.LT.0) IX=(-N+1)*INCX+1
        DO I=1,N
          DX(IX)=ALPHA
          IX=IX+INCX
        END DO
        RETURN
      ELSE
        M=MOD(N,7)
        IF(M.NE.0) THEN
          DO I=1,M
            DX(I)=ALPHA
          END DO
          IF(N.LT.7) RETURN
        END IF
        MP1=M+1
        DO I=MP1,N,7
          DX(I)=ALPHA
          DX(I+1)=ALPHA
          DX(I+2)=ALPHA
          DX(I+3)=ALPHA
          DX(I+4)=ALPHA
          DX(I+5)=ALPHA
          DX(I+6)=ALPHA
        END DO
        RETURN
      END IF
      END
      SUBROUTINE SMUMPS_146( MYID, root, N, IROOT,
     *           COMM, IW, LIW, IFREE,
     *           A, LA, PTRAST, PTLUST_S, PTRFAC,
     *           STEP, INFO, LDLT, QR,
     *           WK, LWK, KEEP,KEEP8)
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER N, IROOT, COMM, LIW, LA, MYID, LWK, LIWK, IFREE
      REAL WK( LWK )
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER PTRAST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *STEP(N), IW( LIW )
      INTEGER INFO( 2 ), LDLT, QR
      REAL A( LA )
      INTEGER IOLDPS, IAPOS
      INTEGER LOCAL_M, LOCAL_N, LPIV, IERR, allocok,i
      INCLUDE 'mumps_headers.h'
        IF ( .NOT. root%yes ) RETURN
        IF ( KEEP(60) .NE. 0 ) THEN
          IF ((LDLT == 1 .OR. LDLT == 2) .AND. KEEP(60) == 3 ) THEN
            CALL SMUMPS_320( WK, root%MBLOCK,
     *      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     *      root%SCHUR_POINTER(1),
     *      root%SCHUR_LLD, root%SCHUR_NLOC,
     *      root%TOT_ROOT_SIZE, MYID, COMM )
          ENDIF
        RETURN
        ENDIF
        IOLDPS  = PTLUST_S(STEP(IROOT))+XSIZE
        IAPOS   = PTRAST(STEP(IROOT))
        LOCAL_M = IW( IOLDPS + 2 )
        LOCAL_N = IW( IOLDPS + 1 )
        IAPOS = PTRFAC(IW ( IOLDPS + 4 ))
        IF ( LDLT.EQ.0 .OR. LDLT.EQ.2 .OR. QR.ne.0 ) THEN
         LPIV = LOCAL_M + root%MBLOCK
        ELSE
         LPIV = 1
        END IF
        IF (ASSOCIATED( root%IPIV )) DEALLOCATE(root%IPIV)
        root%LPIV = LPIV
        ALLOCATE( root%IPIV( LPIV ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = LPIV
          WRITE(*,*) MYID,': problem allocating IPIV(',LPIV,') in root'
          CALL SMUMPS_ABORT()
        END IF
        CALL DESCINIT( root%DESCRIPTOR, root%TOT_ROOT_SIZE,
     *      root%TOT_ROOT_SIZE, root%MBLOCK, root%NBLOCK,
     *      0, 0, root%CNTXT_BLACS, LOCAL_M, IERR )
        IF ( LDLT.EQ.2 ) THEN
            IF(root%MBLOCK.NE.root%NBLOCK) THEN
              WRITE(*,*) ' Error: symmetrization only works for'
              WRITE(*,*) ' square block sizes, MBLOCK/NBLOCK=',
     *        root%MBLOCK, root%NBLOCK
              CALL SMUMPS_ABORT()
            END IF
            IF ( LWK .LT.
     *        MIN(root%MBLOCK * root%NBLOCK,root%TOT_ROOT_SIZE*
     *        root%TOT_ROOT_SIZE ) ) THEN
              WRITE(*,*) 'Not enough workspace for symmetrization.'
              CALL SMUMPS_ABORT()
            END IF
            CALL SMUMPS_320( WK, root%MBLOCK,
     *      root%MYROW, root%MYCOL, root%NPROW, root%NPCOL,
     *      A( IAPOS ), LOCAL_M, LOCAL_N,
     *      root%TOT_ROOT_SIZE, MYID, COMM )
        END IF
        IF (LDLT.EQ.0.OR.LDLT.EQ.2) THEN
          CALL PSGETRF( root%TOT_ROOT_SIZE, root%TOT_ROOT_SIZE,
     *      A( IAPOS ),
     *      1, 1, root%DESCRIPTOR, root%IPIV(1), IERR )
          IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
          END IF
        ELSE
          CALL PSPOTRF('L',root%TOT_ROOT_SIZE,A(IAPOS),
     *      1,1,root%DESCRIPTOR,IERR)
            IF ( IERR .GT. 0 ) THEN
              INFO(1)=-10
              INFO(2)=IERR-1
            END IF
        END IF
        RETURN
      END SUBROUTINE SMUMPS_146
      SUBROUTINE SMUMPS_556(
     *     N,PIV,FRERE,FILS,NFSIZ,IKEEP,
     *     NCST,KEEP,KEEP8,id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (SMUMPS_STRUC) :: id
      INTEGER N,NCST
      INTEGER PIV(N),FRERE(N),FILS(N),NFSIZ(N),IKEEP(N,3)
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER I,P11,P1,P2,K1,K2,NLOCKED
      LOGICAL V1,V2
      NCST = 0
      NLOCKED = 0
      P11 = KEEP(93)
      DO I=KEEP(93)-1,1,-2
         P1 = PIV(I)
         P2 = PIV(I+1)
         K1 = IKEEP(P1,1)
         IF(K1 .GT. 0) THEN
            V1 = (ABS(id%A(K1)*(id%ROWSCA(P1)**2)).GE.1.0E-1)
         ELSE
            V1 = .FALSE.
         ENDIF
         K2 = IKEEP(P2,1)
         IF(K2 .GT. 0) THEN
            V2 = (ABS(id%A(K2)*(id%ROWSCA(P2)**2)).GE.1.0E-1)
         ELSE
            V2 = .FALSE.
         ENDIF
         IF(V1 .AND. V2) THEN
            PIV(P11) = P1
            P11 = P11 - 1
            PIV(P11) = P2
            P11 = P11 - 1
         ELSE IF(V1) THEN
            NCST = NCST+1
            FRERE(NCST) = P1
            NCST = NCST+1
            FRERE(NCST) = P2
         ELSE IF(V2) THEN
            NCST = NCST+1
            FRERE(NCST) = P2                
            NCST = NCST+1
            FRERE(NCST) = P1
         ELSE
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P1
            NLOCKED = NLOCKED + 1
            FILS(NLOCKED) = P2                   
         ENDIF
      ENDDO
      DO I=1,NLOCKED
         PIV(I) = FILS(I)
      ENDDO
      KEEP(94) = KEEP(94) + KEEP(93) - NLOCKED
      KEEP(93) = NLOCKED
      DO I=1,NCST
         NLOCKED = NLOCKED + 1
         PIV(NLOCKED) = FRERE(I)
      ENDDO
      DO I=1,KEEP(93)/2
         NFSIZ(I) = 0
      ENDDO
      DO I=(KEEP(93)/2)+1,(KEEP(93)/2)+NCST,2
         NFSIZ(I) = I+1
         NFSIZ(I+1) = -1
      ENDDO
      DO I=(KEEP(93)/2)+NCST+1,(KEEP(93)/2)+KEEP(94)
         NFSIZ(I) = 0
      ENDDO
      END SUBROUTINE SMUMPS_556
      SUBROUTINE SMUMPS_548(N,PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),NV(N),WORK(N)
      INTEGER I,FATHER,LEN,K,NEWSON,NEWFATHER
      DO I=1,N
         IF(NV(I) .GT. 0) CYCLE
         LEN = 1
         WORK(LEN) = I
         FATHER = -PE(I)
         DO
            IF(NV(FATHER) .GT. 0) THEN
               NEWSON = FATHER
               EXIT
            ENDIF
            LEN = LEN + 1
            WORK(LEN) = FATHER
            NV(FATHER) = 1
            FATHER = -PE(FATHER)
         ENDDO
         NEWFATHER = -PE(FATHER)
         PE(WORK(LEN)) = -NEWFATHER
         PE(NEWSON) = -WORK(1)
      ENDDO      
      END SUBROUTINE SMUMPS_548
      SUBROUTINE SMUMPS_549(N,PE,PERM,NFILS,WORK)
      IMPLICIT NONE
      INTEGER N
      INTEGER PE(N),PERM(N),NFILS(N),WORK(N)
      INTEGER I,FATHER,STKLEN,STKPOS,PERMPOS,CURVAR
      NFILS = 0
      DO I=1,N
         FATHER = -PE(I)
         IF(FATHER .NE. 0) NFILS(FATHER) = NFILS(FATHER) + 1
      ENDDO
      STKLEN = 0
      PERMPOS = 1
      DO I=1,N
         IF(NFILS(I) .EQ. 0) THEN
            STKLEN = STKLEN + 1
            WORK(STKLEN) = I
            PERM(I) = PERMPOS
            PERMPOS = PERMPOS + 1
         ENDIF
      ENDDO
      DO STKPOS = 1,STKLEN
         CURVAR = WORK(STKPOS)
         FATHER = -PE(CURVAR)
         DO
            IF(FATHER .EQ. 0) EXIT
            IF(NFILS(FATHER) .EQ. 1) THEN
               PERM(FATHER) = PERMPOS
               FATHER = -PE(FATHER)
               PERMPOS = PERMPOS + 1
            ELSE
               NFILS(FATHER) = NFILS(FATHER) - 1
               EXIT
            ENDIF
         ENDDO
      ENDDO
      END SUBROUTINE SMUMPS_549
      SUBROUTINE SMUMPS_567(N,NCMP,N11,N22,PIV,
     *     PE,PE_EXP)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),PE(NCMP),PE_EXP(N)
      INTEGER I,J,N2,K,THEROOT,FATHER,K1,K2
      N2 = N22/2
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
#endif
      IF(N22+N11 .EQ. N) THEN
         THEROOT = 0
      ELSE
         THEROOT = PIV(N)
         PE_EXP(THEROOT) = 0
      ENDIF
      K = 1
      DO I=1,N2
         K1 = PIV(K)
         K = K+1
         K2 = PIV(K)
         K = K+1
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -K2
         PE_EXP(K2) = -FATHER
      ENDDO
#if(SP_DEBUG >=1)
      IF(K .NE. N22+1) THEN
         WRITE(*,*) 'K NE N22+1',K,N22+1
         CALL SMUMPS_ABORT()
      ENDIF
#endif
      DO I=N2+1,NCMP
         K1 = PIV(K)
         K = K+1
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -FATHER
      ENDDO
      DO I=N11+N22+1,N-1
         K1 = PIV(I)
         PE(K1) = -THEROOT
      ENDDO
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PE_EXP',PE_EXP
#endif
      END SUBROUTINE SMUMPS_567
      SUBROUTINE SMUMPS_566(N,LIW,COLPTR,IW,N11,N22,PIV,
     *     PE,NV,WORK)
      IMPLICIT NONE
      INTEGER N11,N22,N,LIW
      INTEGER PIV(N),NV(N),PE(N),COLPTR(N+1),IW(LIW),WORK(N)
      INTEGER I,J,K,COL,FLAG,THEROOT
      FLAG = 1
      WORK = 0
      THEROOT = PIV(N)
      DO I = N11+N22+1,N
         J = PIV(I)
         WORK(THEROOT) = FLAG
         DO K=COLPTR(J),COLPTR(J+1)-1
            COL = IW(K)
            DO
               IF(NV(COL) .GT. 0) EXIT
               COL = -PE(COL)
            ENDDO
            DO
               IF(WORK(COL) .EQ. FLAG) EXIT
               NV(COL) = NV(COL) + 1
               WORK(COL) = FLAG
               COL = -PE(COL)
            ENDDO
         ENDDO
         IF(FLAG+1 .GT. FLAG) THEN
            FLAG = FLAG+1
         ELSE
            WORK = 0
            FLAG = 1
         ENDIF
      ENDDO
#if (SP_DEBUG >= 2)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
      WRITE(*,*) 'NV',NV
#endif
      END SUBROUTINE SMUMPS_566
      SUBROUTINE SMUMPS_565(N,NCMP,N11,N22,PIV,
     *     PE,NV,PE_EXP,NV_EXP)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),NV(NCMP),PE(NCMP),NV_EXP(N),PE_EXP(N)
      INTEGER I,J,N2,K,THEROOT,FATHER,K1,K2
      N2 = N22/2
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PE',PE
      WRITE(*,*) 'NV',NV
#endif
      IF(N22+N11 .EQ. N) THEN
         THEROOT = 0
      ELSE
         THEROOT = PIV(N)
         PE_EXP(THEROOT) = 0
         NV_EXP(THEROOT) = N-N22-N11
      ENDIF
      K = 1
      DO I=1,N2
         K1 = PIV(K)
         K = K+1
         K2 = PIV(K)
         K = K+1
         NV_EXP(K1) = 0
         NV_EXP(K2) = NV(I)
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -K2
         PE_EXP(K2) = -FATHER
      ENDDO
#if(SP_DEBUG >=1)
      IF(K .NE. N22+1) THEN
         WRITE(*,*) 'K NE N22+1',K,N22+1
         CALL SMUMPS_ABORT()
      ENDIF
#endif
      DO I=N2+1,NCMP
         K1 = PIV(K)
         K = K+1
         NV_EXP(K1) = NV(I)
         FATHER = -PE(I)
         IF(FATHER .EQ. 0) THEN
            FATHER = THEROOT
         ELSE IF(FATHER .LE. N2) THEN
            FATHER = PIV(2*FATHER)
         ELSE
            FATHER = PIV(N2+FATHER)
         ENDIF
         PE_EXP(K1) = -FATHER
      ENDDO
      DO I=N11+N22+1,N-1
         K1 = PIV(I)
         NV(K1) = 0
         PE(K1) = -THEROOT
      ENDDO
#if (SP_DEBUG >= 5)
      WRITE(*,*) 'PE_EXP',PE_EXP
      WRITE(*,*) 'NV_EXP',NV_EXP
#endif
      END SUBROUTINE SMUMPS_565
      SUBROUTINE SMUMPS_550(N,NCMP,N11,N22,PIV,
     *     PERM,IPERM)
      IMPLICIT NONE
      INTEGER N11,N22,N,NCMP
      INTEGER PIV(N),PERM(N),IPERM(N)
      INTEGER CMP_POS,EXP_POS,I,J,N2,K
      N2 = N22/2
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'PERM',PERM(1:NCMP)
      WRITE(*,*) 'IPERM',IPERM(1:NCMP)
#endif
      EXP_POS = 1
      DO CMP_POS=1,NCMP
         J = IPERM(CMP_POS)
         IF(J .LE. N2) THEN
            K = 2*J-1
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
            K = K+1
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ELSE
            K = N2 + J
            I = PIV(K)
            PERM(I) = EXP_POS
            EXP_POS = EXP_POS+1
         ENDIF
      ENDDO
      DO K=N22+N11+1,N
         I = PIV(K)
         PERM(I) = EXP_POS
         EXP_POS = EXP_POS+1
      ENDDO
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PERM',PERM
#endif
      END SUBROUTINE SMUMPS_550
      SUBROUTINE SMUMPS_547(
     *     N,NZ, IRN, ICN, PIV,
     *     NCMP, IW, LW, IPE, LEN, IQ, 
     *     FLAG, ICMP, IWFR,
     *     IERROR, KEEP,KEEP8, ICNTL)
      IMPLICIT NONE
      INTEGER N,NZ,NCMP,LW,IWFR,IERROR
      INTEGER ICNTL(40),KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER IRN(NZ),ICN(NZ),IW(LW),PIV(N),IPE(N+1)
      INTEGER LEN(N),IQ(N),FLAG(N),ICMP(N) 
      INTEGER MP,N11,N22,NDUP
      INTEGER I,K,J,N1,LAST,K1,K2,L
      MP = ICNTL(2)
      IERROR = 0
      N22 = KEEP(93)
      N11 = KEEP(94)
      NCMP = N22/2 + N11
      DO I=1,NCMP
         IPE(I) = 0
      ENDDO
      K = 1
      DO I=1,N22/2
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
         J = PIV(K)
         ICMP(J) = I
         K = K + 1
      ENDDO
      K = N22/2 + 1
      DO I=N22+1,N22+N11
         J = PIV(I)
         ICMP(J) = K
         K = K + 1
      ENDDO
      DO I=N11+N22+1,N
         J = PIV(I)
         ICMP(J) = 0
      ENDDO
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
         IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *        .OR.(J.LT.1)) THEN
            IERROR = IERROR + 1
         ELSE
            IF (I.NE.J) THEN
               IPE(I) = IPE(I) + 1
               IPE(J) = IPE(J) + 1
            ENDIF
         ENDIF
      ENDDO
      IQ(1) = 1
      N1 = NCMP - 1
      IF (N1.GT.0) THEN
         DO I=1,N1
            IQ(I+1) = IPE(I) + IQ(I)
         ENDDO 
      ENDIF
      LAST = MAX(IPE(NCMP)+IQ(NCMP)-1,IQ(NCMP))
      DO I = 1,NCMP
         FLAG(I) = 0
         IPE(I)  = IQ(I)
      ENDDO
      DO K=1,LAST
        IW(K) = 0
      ENDDO
      IWFR = LAST + 1
      DO K=1,NZ
         I = IRN(K)
         J = ICN(K)
         I = ICMP(I)
         J = ICMP(J)
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
      ENDDO
      NDUP = 0
      DO I=1,NCMP
         K1 = IPE(I) 
         K2 = IQ(I) -1
         IF (K1.GT.K2) THEN
            LEN(I) = 0
            IQ(I)  = 0
         ELSE
            DO K=K1,K2
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
            ENDDO
 250        IQ(I) = IQ(I) - IPE(I)
            IF (NDUP.EQ.0) LEN(I) = IQ(I)
         ENDIF
      ENDDO
#if (SP_DEBUG >= 2)
      WRITE(*,*) 'NDUP =',NDUP
#endif
      IF (NDUP.NE.0) THEN
         IWFR = 1
         DO I=1,NCMP
            K1 = IPE(I) 
            IF (IQ(I).EQ.0) THEN
               LEN(I) = 0
               IPE(I) = IWFR
               CYCLE
            ENDIF
            K2 = K1 + IQ(I) - 1
            L = IWFR
            IPE(I) = IWFR
            DO K=K1,K2
               IF (IW(K).NE.0) THEN
                  IW(IWFR) = IW(K)
                  IWFR     = IWFR + 1
               ENDIF
            ENDDO
            LEN(I) = IWFR - L 
         ENDDO
      ENDIF
      IPE(NCMP+1) = IPE(NCMP) + LEN(NCMP)
      IWFR = IPE(NCMP+1)
#if(SP_DEBUG >= 5)
      WRITE(*,*) 'PIV',PIV
      WRITE(*,*) 'ICMP',ICMP
      WRITE(*,*) 'IPE',IPE(1:NCMP+1)
      WRITE(*,*) 'LEN',LEN(1:NCMP)
      WRITE(*,*) 'IW',IW(1:IWFR-1)
#endif
      RETURN
      END SUBROUTINE SMUMPS_547
      SUBROUTINE SMUMPS_551(
     *     N, NE, IP, IRN, SCALING,LSC,CPERM, DIAG,
     *     ICNTL, WEIGHT,MARKED,FLAG,
     *     PIV_OUT, INFO)
      IMPLICIT NONE
      INTEGER N, NE, ICNTL(10), INFO(10),LSC,LWEIGHT
      INTEGER CPERM(N),PIV_OUT(N),IP(N+1),IRN(NE), DIAG(N)
      REAL SCALING(LSC),WEIGHT(N+2)
      INTEGER MARKED(N),FLAG(N)
      INTEGER NUM1,NUM2,NUMTOT,PATH_LENGTH,NLAST
      INTEGER I,J,K,L,CUR_EL,CUR_EL_PATH,CUR_EL_PATH_NEXT,BEST_BEG
      INTEGER L1,L2,PTR_SET1,PTR_SET2,TUP,T22,INTER,MERGE
      REAL BEST_SCORE,CUR_VAL,TMP,VAL
      REAL INITSCORE, SMUMPS_UPDATESCORE, 
     *     SMUMPS_UPDATE_INVERSE, SMUMPS_METRIC2x2
#if(SP_DEBUG >= 1)
      INTERFACE
         FUNCTION SMUMPS_UPDATESCORE(A,B,T)
         REAL A,B
         INTEGER T
         END FUNCTION
         FUNCTION SMUMPS_UPDATE_INVERSE(A,B,T)
         REAL A,B
         INTEGER T
         END FUNCTION SMUMPS_UPDATE_INVERSE
         FUNCTION SMUMPS_METRIC2x2(CUR_EL,CUR_EL_PATH,
     *        SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
         INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
         INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
         REAL VAL
         LOGICAL FLAGON
         INTEGER T
         END FUNCTION SMUMPS_METRIC2x2
      END INTERFACE
#endif
      LOGICAL VRAI,FAUX,MAX_CARD_DIAG,USE_SCALING
      INTEGER JOB_DEF,DONE
      INTEGER SUM,PROD,STRUCT,MA47,MAGNITUDE
      INTEGER ZERO,ONE
      PARAMETER (JOB_DEF = 5, DONE = -1,
     *     SUM = 1, PROD = 2, STRUCT=1, MA47=2, MAGNITUDE=3,
     *     VRAI = .TRUE., FAUX = .FALSE.)
      PARAMETER(ZERO = 0.0E0, ONE = 1.0E0)
      MAX_CARD_DIAG = .TRUE.
      NUM1 = 0
      NUM2 = 0
      NUMTOT = 0
      NLAST = N
      INFO = 0      
      MARKED = 1
      FLAG = 0
      VAL = ONE
      IF(LSC .GT. 1) THEN
         USE_SCALING = .TRUE.
      ELSE
         USE_SCALING = .FALSE.
      ENDIF
      TUP = ICNTL(2)
      IF(TUP .EQ. SUM) THEN
        INITSCORE = ZERO
      ELSE
        INITSCORE = ONE
      ENDIF
      IF(ICNTL(2) .GT. 2 .OR. ICNTL(2) .LE. 0) THEN
         WRITE(*,*)
     *        'ERROR: WRONG VALUE FOR ICNTL(2) = ',ICNTL(2)
         INFO(1) = -1
         RETURN
      ENDIF
      T22 = ICNTL(1)
      IF(ICNTL(1) .LT. 0 .OR. ICNTL(1) .GT. 2) THEN
         WRITE(*,*)
     *        'ERROR: WRONG VALUE FOR ICNTL(1) = ',ICNTL(1)
         INFO(1) = -1
         RETURN
      ENDIF
      DO CUR_EL=1,N
         IF(MARKED(CUR_EL) .LE. 0) THEN
            CYCLE
         ENDIF
         IF(CPERM(CUR_EL) .LT. 0) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         PATH_LENGTH = 2
         CUR_EL_PATH = CPERM(CUR_EL)
         IF(CUR_EL_PATH .EQ. CUR_EL) THEN
            MARKED(CUR_EL) = -1
            CYCLE
         ENDIF
         MARKED(CUR_EL) = 0
         WEIGHT(1) = INITSCORE
         WEIGHT(2) = INITSCORE
         L1 = IP(CUR_EL+1)-IP(CUR_EL)
         L2 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
         PTR_SET1 = IP(CUR_EL)
         PTR_SET2 = IP(CUR_EL_PATH)
         IF(USE_SCALING) THEN
            VAL = -SCALING(CUR_EL_PATH) - SCALING(CUR_EL+N)
         ENDIF
         CUR_VAL = SMUMPS_METRIC2x2(
     *        CUR_EL,CUR_EL_PATH,
     *        IRN(PTR_SET1),IRN(PTR_SET2),
     *        L1,L2,
     *        VAL,DIAG,N,FLAG,FAUX,T22)
         WEIGHT(PATH_LENGTH+1) = 
     *        SMUMPS_UPDATESCORE(WEIGHT(1),CUR_VAL,TUP)
         DO
            IF(CUR_EL_PATH .EQ. CUR_EL) EXIT
            PATH_LENGTH = PATH_LENGTH+1
            MARKED(CUR_EL_PATH) = 0
            CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
            L1 = IP(CUR_EL_PATH+1)-IP(CUR_EL_PATH)
            L2 = IP(CUR_EL_PATH_NEXT+1)-IP(CUR_EL_PATH_NEXT)
            PTR_SET1 = IP(CUR_EL_PATH)
            PTR_SET2 = IP(CUR_EL_PATH_NEXT)
            IF(USE_SCALING) THEN
               VAL = -SCALING(CUR_EL_PATH_NEXT) 
     *              - SCALING(CUR_EL_PATH+N)
            ENDIF
            CUR_VAL = SMUMPS_METRIC2x2(
     *           CUR_EL_PATH,CUR_EL_PATH_NEXT,
     *           IRN(PTR_SET1),IRN(PTR_SET2),
     *           L1,L2,
     *           VAL,DIAG,N,FLAG,VRAI,T22)
            WEIGHT(PATH_LENGTH+1) = 
     *           SMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH-1),CUR_VAL,TUP)
            CUR_EL_PATH = CUR_EL_PATH_NEXT
         ENDDO
         IF(MOD(PATH_LENGTH,2) .EQ. 1) THEN
            IF(WEIGHT(PATH_LENGTH+1) .GE. WEIGHT(PATH_LENGTH)) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
            ELSE
               CUR_EL_PATH = CUR_EL
            ENDIF
            DO I=1,(PATH_LENGTH-1)/2
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 1
         ELSE
            IF(MAX_CARD_DIAG) THEN
               CUR_EL_PATH = CPERM(CUR_EL)
               IF(DIAG(CUR_EL) .NE. 0) THEN 
                  BEST_BEG = CUR_EL_PATH
                  GOTO 1000
               ENDIF
               DO I=1,(PATH_LENGTH/2)
                  CUR_EL_PATH_NEXT = CPERM(CUR_EL_PATH)
                  IF(DIAG(CUR_EL_PATH) .NE. 0) THEN 
                     BEST_BEG = CUR_EL_PATH_NEXT
                     GOTO 1000
                  ENDIF
               ENDDO
            ENDIF
            BEST_BEG = CUR_EL
            BEST_SCORE = WEIGHT(PATH_LENGTH-1)
            CUR_EL_PATH = CPERM(CUR_EL)
            DO I=1,(PATH_LENGTH/2)-1
               TMP = SMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH),
     *              WEIGHT(2*I-1),TUP)
               TMP = SMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               TMP = SMUMPS_UPDATESCORE(WEIGHT(PATH_LENGTH+1),
     *              WEIGHT(2*I),TUP)
               TMP = SMUMPS_UPDATE_INVERSE(TMP,WEIGHT(2*I+1),TUP)
               IF(TMP .GT. BEST_SCORE) THEN
                  BEST_SCORE = TMP
                  BEST_BEG = CUR_EL_PATH
               ENDIF
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
 1000       CUR_EL_PATH = BEST_BEG
            DO I=1,(PATH_LENGTH/2)-1
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
               NUM2 = NUM2+1
               PIV_OUT(NUM2) = CUR_EL_PATH
               CUR_EL_PATH = CPERM(CUR_EL_PATH)
            ENDDO
            NUMTOT = NUMTOT + PATH_LENGTH - 2
            MARKED(CUR_EL_PATH) = -1
         ENDIF
      ENDDO
      DO I=1,N
         IF(MARKED(I) .LT. 0) THEN
            IF(DIAG(I) .EQ. 0) THEN
               PIV_OUT(NLAST) = I
               NLAST = NLAST - 1
            ELSE
               NUM1 = NUM1 + 1
               PIV_OUT(NUM2+NUM1) = I
               NUMTOT = NUMTOT + 1
            ENDIF
         ENDIF
      ENDDO
      INFO(2) = NUMTOT
      INFO(3) = NUM1
      INFO(4) = NUM2
      RETURN
      END SUBROUTINE SMUMPS_551
      FUNCTION SMUMPS_UPDATESCORE(A,B,T)
      IMPLICIT NONE
      REAL SMUMPS_UPDATESCORE
      REAL A,B
      INTEGER T
      INTEGER SUM,PROD
      PARAMETER(SUM = 1,PROD = 2)
      IF(T .EQ. SUM) THEN
         SMUMPS_UPDATESCORE = A+B
      ELSE
         SMUMPS_UPDATESCORE = A*B
      ENDIF
      END FUNCTION SMUMPS_UPDATESCORE
      FUNCTION SMUMPS_UPDATE_INVERSE(A,B,T)
      IMPLICIT NONE
      REAL SMUMPS_UPDATE_INVERSE
      REAL A,B
      INTEGER T
      INTEGER SUM,PROD
      PARAMETER(SUM = 1,PROD = 2)
      IF(T .EQ. SUM) THEN
         SMUMPS_UPDATE_INVERSE = A-B
      ELSE
         SMUMPS_UPDATE_INVERSE = A/B
      ENDIF
      END FUNCTION SMUMPS_UPDATE_INVERSE
      FUNCTION SMUMPS_METRIC2x2(CUR_EL,CUR_EL_PATH,
     *     SET1,SET2,L1,L2,VAL,DIAG,N,FLAG,FLAGON,T)
      IMPLICIT NONE
      REAL SMUMPS_METRIC2x2
      INTEGER CUR_EL,CUR_EL_PATH,L1,L2,N
      INTEGER SET1(L1),SET2(L2),DIAG(N),FLAG(N)
      REAL VAL
      LOGICAL FLAGON
      INTEGER T
      INTEGER I,INTER,MERGE
      INTEGER STRUCT,MA47,MAGNITUDE
      PARAMETER(STRUCT=0,MA47=1,MAGNITUDE=2)
      IF(T .EQ. STRUCT) THEN
         IF(.NOT. FLAGON) THEN
            DO I=1,L1
               FLAG(SET1(I)) = CUR_EL
            ENDDO            
         ENDIF
         INTER = 0
         DO I=1,L2
            IF(FLAG(SET2(I)) .EQ. CUR_EL) THEN
               INTER = INTER + 1
               FLAG(SET2(I)) = CUR_EL_PATH
            ENDIF
         ENDDO
         MERGE = L1 + L2 - INTER
         SMUMPS_METRIC2x2 = REAL(INTER) / REAL(MERGE)
      ELSE IF (T .EQ. MA47) THEN
         MERGE = 3
         IF(DIAG(CUR_EL) .NE. 0) MERGE = 2
         IF(DIAG(CUR_EL_PATH) .NE. 0) MERGE = MERGE - 2
         IF(MERGE .EQ. 0) THEN
            SMUMPS_METRIC2x2 = L1+L2-2
            SMUMPS_METRIC2x2 = -(SMUMPS_METRIC2x2**2)/2
         ELSE IF(MERGE .EQ. 1) THEN
            SMUMPS_METRIC2x2 = - REAL(L1+L2-4) * REAL(L1-2)
         ELSE IF(MERGE .EQ. 2) THEN
            SMUMPS_METRIC2x2 = - REAL(L1+L2-4) * REAL(L2-2)
         ELSE
            SMUMPS_METRIC2x2 = - REAL(L1-2) * REAL(L2-2)
         ENDIF
      ELSE
         SMUMPS_METRIC2x2 = VAL
      ENDIF
      END FUNCTION 
      SUBROUTINE SMUMPS_622(NA, NCMP,
     *      PERM,IPERM, 
     *      LISTVAR_SCHUR, SIZE_SCHUR, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN):: NA, NCMP
      INTEGER, INTENT(IN):: AOTOA(NCMP), IPERM(NCMP)
      INTEGER, INTENT(OUT):: PERM(NA) ! not used a input
      INTEGER CMP_POS, IO, I, K, IPOS
      DO CMP_POS=1, NCMP
        IO              = IPERM(CMP_POS)
        PERM(AOTOA(IO)) = CMP_POS
      ENDDO
      IPOS = NCMP
      DO K =1,  SIZE_SCHUR
        I       = LISTVAR_SCHUR(K)
        IPOS    = IPOS+1
        PERM(I) = IPOS
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_622
      SUBROUTINE SMUMPS_623
     * (NA,N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     * IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL, 
     * symmetry, SYM, MedDens, NBQD, AvgDens, 
     * LISTVAR_SCHUR, SIZE_SCHUR, ATOAO, AOTOA)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NA,N,NZ,LW
      INTEGER, INTENT(IN)  :: SIZE_SCHUR, LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)  :: IRN(NZ), ICN(NZ) 
      INTEGER, INTENT(IN)  :: ICNTL(40), SYM
      INTEGER, INTENT(OUT) :: IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER, INTENT(OUT) :: AOTOA(N)
      INTEGER, INTENT(OUT) :: ATOAO(NA)  ! only used internallly
      INTEGER, INTENT(OUT) :: LEN(N), IPE(N+1)
      INTEGER, INTENT(OUT) :: symmetry, 
     &                        MedDens, NBQD, AvgDens
      INTEGER, INTENT(OUT)  :: FLAG(N), IW(LW), IQ(N)
      INTEGER MP, MPG, NAO
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR, THRESH, IAO
      INTEGER NZOFFA, NDIAGA
      REAL RSYM
      INTRINSIC NINT
      ATOAO(1:NA) = 0
      DO I = 1, SIZE_SCHUR
        ATOAO(LISTVAR_SCHUR(I)) = -1
      ENDDO
      IAO = 0  
      DO I= 1, NA
        IF (ATOAO(I).LT.0) CYCLE
        IAO = IAO +1   
        ATOAO(I)   = IAO
        AOTOA(IAO) = I
      ENDDO
      MP = ICNTL(2)
      MPG= ICNTL(3)
      NIORM  = 3*N
      NDIAGA = 0
      IERROR = 0
      IPE(1:N+1) = 0
      DO K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          I = ATOAO(I)
          J = ATOAO(J)
          IF ((I.LT.0).OR.(J.LT.0)) CYCLE  ! in schur
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ELSE
           NDIAGA = NDIAGA + 1
          ENDIF
        ENDIF
      ENDDO
      NZOFFA  = NIORM - 3*N
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IF (MOD(IFLAG,2).EQ.0) IFLAG  = IFLAG+1
         IF ((MP.GT.0).AND.(ICNTL(4).GE.2))  THEN 
          WRITE (MP,99999) 
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
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
         IF ((I.GT.NA).OR.(J.GT.NA).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) CYCLE
         I = ATOAO(I)
         J = ATOAO(J)
         IF ((I.LT.0).OR.(J.LT.0)) CYCLE  ! in schur
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1 
          ELSE
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
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
      END SUBROUTINE SMUMPS_623
