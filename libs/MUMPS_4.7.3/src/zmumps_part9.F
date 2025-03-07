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
      SUBROUTINE ZMUMPS_448(ICNTL,CNTL)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL
      PARAMETER (NICNTL=10, NCNTL=10)
      INTEGER ICNTL(NICNTL)
      DOUBLE PRECISION CNTL(NCNTL)
      INTEGER I
      ICNTL(1) =  6
      ICNTL(2) =  6
      ICNTL(3) = -1
      ICNTL(4) = -1
      ICNTL(5) =  0
      DO 10 I = 6,NICNTL
        ICNTL(I) = 0
   10 CONTINUE
      CNTL(1) = 0.0D0
      CNTL(2) = 0.0D0
      DO 20 I = 3,NCNTL
        CNTL(I) = 0.0D0
   20 CONTINUE
      RETURN
      END
#if defined(perfnal)
#if defined(NOT_USED)
      SUBROUTINE ZMUMPS_443(JOB,M,N,NE,IP,IRN,A,NUM,PERM,LIW,IW,LDW,DW,
     &           ICNTL,CNTL,INFO)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL=10, NCNTL=10, NINFO=10)
      INTEGER JOB,M,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),PERM(M),IW(LIW)
      INTEGER ICNTL(NICNTL),INFO(NINFO)
      DOUBLE PRECISION A(NE)
      DOUBLE PRECISION DW(LDW),CNTL(NCNTL)
      INTEGER I,J,K,wARN1,WARN2,WARN4
      DOUBLE PRECISION FACT,ZERO,ONE,RINF
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+0)
      EXTERNAL ZMUMPS_457,ZMUMPS_444,ZMUMPS_451,ZMUMPS_452,ZMUMPS_454
      INTRINSIC ABS,LOG
       RINF = CNTL(2)
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
        IF (JOB.EQ.3) K =           NE
        IF (JOB.EQ.4) K =     2*M + NE
        IF (JOB.EQ.5) K = N + 2*M + NE
        IF (JOB.EQ.6) K = N + 3*M + NE
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
    3   CONTINUE
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
    4     CONTINUE
    6   CONTINUE
      ENDIF
      IF (CNTL(2).LE.ZERO) THEN
        INFO(1) = -8
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9008) INFO(1) 
        GO TO 99
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
   10   CONTINUE
        CALL ZMUMPS_457(M,N,IRN,NE,IP,IW(1),PERM,NUM,
     &       IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
        DW(1) = MAX(ZERO,CNTL(1))
        CALL ZMUMPS_444(M,N,NE,IP,IRN,A,PERM,NUM,
     &       IW(1),IW(N+1),IW(2*N+1),IW(2*N+M+1),DW,RINF)
        GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
        CALL ZMUMPS_451(N,NE,IP,IW,DW)
        FACT = MAX(ZERO,CNTL(1))
        CALL ZMUMPS_452(M,N,NE,IP,IW(1),DW,PERM,NUM,IW(NE+1),
     &       IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &       IW(NE+5*N+1),IW(NE+5*N+M+1),FACT,RINF)
        GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*M+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
        DW(1) = MAX(ZERO,CNTL(1))
        IW(1) = JOB
        CALL ZMUMPS_454(M,N,NE,IP,IRN,DW(2*M+1),PERM,NUM,
     &       IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &       DW(1),DW(M+1),RINF)
        GO TO 90
      ENDIF
      IF (JOB.EQ.5 .or. JOB.EQ.6) THEN
        IF (JOB.EQ.5) THEN
          DO 75 J = 1,N
            FACT = ZERO
            DO 60 K = IP(J),IP(J+1)-1
              DW(2*M+N+K) = ABS(A(K))
              IF (DW(2*M+N+K).GT.FACT) FACT = DW(2*M+N+K)
   60       CONTINUE
            DW(2*M+J) = FACT
            IF (FACT.NE.ZERO) THEN
              FACT = LOG(FACT)
              DO 70 K = IP(J),IP(J+1)-1
                IF (DW(2*M+N+K).NE.ZERO) THEN
                  DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
                ELSE
                  DW(2*M+N+K) = RINF/N
                ENDIF
   70         CONTINUE
            ELSE
              DO 71 K = IP(J),IP(J+1)-1
                DW(2*M+N+K) = ONE
   71         CONTINUE
            ENDIF
   75     CONTINUE
        ENDIF
        IF (JOB.EQ.6) THEN
          DO 175 K = 1,NE
            IW(3*N+2*M+K) = IRN(K)
            DW(2*M+N+K) = ABS(A(K))
  175     CONTINUE
          DO 61 I = 1,M
            DW(2*M+N+NE+I) = ZERO
   61     CONTINUE
          DO 63 J = 1,N
            DO 62 K = IP(J),IP(J+1)-1
              I = IRN(K)
              IF (DW(2*M+N+K).GT.DW(2*M+N+NE+I)) THEN
                DW(2*M+N+NE+I) = DW(2*M+N+K)
              ENDIF
   62       CONTINUE
   63     CONTINUE
          DO 64 I = 1,M
            IF (DW(2*M+N+NE+I).NE.ZERO) THEN
              DW(2*M+N+NE+I) = 1/DW(2*M+N+NE+I)
            ENDIF
   64     CONTINUE
          DO 66 J = 1,N
            DO 65 K = IP(J),IP(J+1)-1
              I = IRN(K)
              DW(2*M+N+K) = DW(2*M+N+NE+I) * DW(2*M+N+K)
   65       CONTINUE
   66     CONTINUE
          CALL ZMUMPS_451(N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1))
          DO 176 J = 1,N
            IF (IP(J).NE.IP(J+1)) THEN
              FACT = DW(2*M+N+IP(J))
            ELSE
              FACT = ZERO
            ENDIF
            DW(2*M+J) = FACT
            IF (FACT.NE.ZERO) THEN
              FACT = LOG(FACT)
              DO 170 K = IP(J),IP(J+1)-1
                IF (DW(2*M+N+K).NE.ZERO) THEN
                  DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
                ELSE
                  DW(2*M+N+K) = RINF/N
                ENDIF
  170         CONTINUE
            ELSE
              DO 171 K = IP(J),IP(J+1)-1
                DW(2*M+N+K) = ONE
  171         CONTINUE
            ENDIF
  176     CONTINUE
        ENDIF
        DW(1) = MAX(ZERO,CNTL(1))
        IW(1) = JOB
        IF (JOB.EQ.5) THEN
          CALL ZMUMPS_454(M,N,NE,IP,IRN,DW(2*M+N+1),PERM,NUM,
     &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &         DW(1),DW(M+1),RINF)
        ENDIF
        IF (JOB.EQ.6) THEN
          CALL ZMUMPS_454(M,N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1),PERM,NUM,
     &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
     &         DW(1),DW(M+1),RINF)
        ENDIF
        IF (JOB.EQ.6) THEN
          DO 79 I = 1,M
            IF (DW(2*M+N+NE+I).NE.0) THEN
              DW(I) = DW(I) + LOG(DW(2*M+N+NE+I))
            ENDIF
   79     CONTINUE
        ENDIF
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*M+J).NE.ZERO) THEN
              DW(M+J) = DW(M+J) - LOG(DW(2*M+J))
            ELSE
              DW(M+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
        FACT = 0.5*LOG(RINF)
        DO 86 I = 1,M
          IF (DW(I).LT.FACT) GO TO 86
          WARN2 = 2
          GO TO 90
   86   CONTINUE 
        DO 87 J = 1,N
          IF (DW(M+J).LT.FACT) GO TO 87
          WARN2 = 2
          GO TO 90
   87   CONTINUE 
      ENDIF
   90 IF (NUM.LT.N) WARN1 = 1
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
   99 RETURN
 9001 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9008 FORMAT (' ****** Error in DMUMPS_443. INFO(1) = ',I2/
     &        '        CNTL(2) is negative or zero.')
 9010 FORMAT (' ****** Warning from DMUMPS_443. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')
 9020 FORMAT (' ****** Input parameters for DMUMPS_443:'/
     &        ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)   = ',8I8/(15X,8I8))
 9022 FORMAT (' IRN(1:NE)   = ',8I8/(15X,8I8))
 9023 FORMAT (' A(1:NE)     = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10) = ',8I8/(15X,2I8))
 9025 FORMAT (' CNTL(1:10)  = ',4(1PD14.4)/(15X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for DMUMPS_443:'/
     &        ' INFO(1:2)   = ',2I8)
 9031 FORMAT (' NUM         = ',I8)
 9032 FORMAT (' PERM(1:M)   = ',8I8/(15X,8I8))
 9033 FORMAT (' DW(1:M)     = ',5(F11.3)/(15X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N) = ',5(F11.3)/(15X,5(F11.3)))
      END
#endif
#endif
      SUBROUTINE ZMUMPS_444
     &(M,N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),PR(N),Q(M),L(M)
      DOUBLE PRECISION A(NE)
      DOUBLE PRECISION D(M), RINF
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV,TBV,RLX
      DOUBLE PRECISION ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC ABS,MIN
      EXTERNAL ZMUMPS_445,ZMUMPS_446,ZMUMPS_447,ZMUMPS_455
      RLX = D(1)
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
   10 CONTINUE
      DO 12 K = 1,M
        IPERM(K) = 0
        D(K) = ZERO
   12 CONTINUE
      DO 30 J = 1,N
        A0 = MINONE
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 20
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 20
            JPERM(J) = I 
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 20
            A0 = AI
            I0 = I
          ENDIF
   20   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 30
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   30 CONTINUE
      IF (M.EQ.N) THEN
        DO 35 I = 1,M
          BV = MIN(BV,D(I))
   35   CONTINUE
      ENDIF
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      TBV = BV * (1-RLX)
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.TBV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.TBV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL ZMUMPS_445(I,M,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            TBV = BV * (1-RLX)
            DO 152 IDUM = 1,M
              CALL ZMUMPS_446(QLEN,M,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).LT.TBV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.TBV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.TBV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.TBV) THEN
                IF (DI.NE.MINONE) THEN
                  CALL ZMUMPS_447(L(I),QLEN,M,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL ZMUMPS_445(I,M,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = MIN(BV,CSP)
        TBV = BV * (1-RLX)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE 
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL ZMUMPS_455(M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE ZMUMPS_445(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI
      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
      SUBROUTINE ZMUMPS_446(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
      SUBROUTINE ZMUMPS_447(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34 
          Q(POS) = QK
          L(QK) = POS 
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END
#if defined(perfnal)
#if defined(NOT_USED)
      SUBROUTINE ZMUMPS_458(M,N,NE,IP,IRN,A,MODE,PERM,IW,INFO)
      IMPLICIT NONE
      INTEGER M,N,NE,MODE,INFO
      DOUBLE PRECISION A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(*),IW(*)
      INTEGER I, K
      EXTERNAL ZMUMPS_456,ZMUMPS_449
      INFO = 0
      IF (M.LE.0 .OR. N.LE.0) THEN
        INFO = -2
        GO TO 299
      ENDIF
      IF (NE.LE.0) THEN
        INFO = -3
        GO TO 299
      ENDIF
      IF (MODE.LT.1 .or. MODE.GT.2) THEN
        INFO = -10
        GO TO 299
      ENDIF
      IF (MODE.EQ.1) THEN
        CALL ZMUMPS_456(M,PERM,IW,INFO)
      ELSE
        CALL ZMUMPS_456(N,PERM,IW,INFO)
      ENDIF
      IF (INFO.NE.0) THEN
        INFO = -11
        GO TO 299
      ENDIF
      IF (MODE.EQ.1) THEN
        K = 0
        DO 100 I = 1,M
          IF (PERM(I).EQ.I) GO TO 100
          K = 1
          GO TO 120
  100   CONTINUE
  120   IF (K.EQ.0) GO TO 299
        DO 150 K = 1,NE
          I = IRN(K)
          IRN(K) = PERM(I) 
  150   CONTINUE
        GO TO 299
      ENDIF
      CALL ZMUMPS_449(M,N,NE,IP,IRN,A,PERM,IW(1),IW(N+1))
  299 RETURN
      END
      SUBROUTINE ZMUMPS_449(M,N,NE,IP,IRN,A,PERM,IW,DEST)
      IMPLICIT NONE
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(N),IW(N),DEST(NE)
      DOUBLE PRECISION AOLD,ANEW
      INTEGER I,J,K,IDUM,IOLD,INEW,I0,I1,JOLD,JNEW,LEN
      K = 0
      DO 200 J = 1,N
        IF (PERM(J).EQ.J) GO TO 200
        K = 1
        GO TO 220
  200 CONTINUE
  220 IF (K.EQ.0) GO TO 299
      K = 1
      DO 240 JNEW = 1,N
        JOLD = PERM(JNEW)
        LEN = IP(JOLD+1) - IP(JOLD)
        IW(JNEW) = LEN
        J = IP(JOLD)
        DO 230 I = 0,LEN-1
          DEST(J+I) = K+I
  230   CONTINUE
        K = K + LEN
  240 CONTINUE
      IP(1) = 1
      DO 250 JNEW = 2,N+1
        IP(JNEW) = IP(JNEW-1) + IW(JNEW-1)
  250 CONTINUE
      DO 290 I = 1,NE
        I0 = I
        I1 = DEST(I)
        IF (I1.EQ.0) GO TO 290
        IF (I1.EQ.I) GO TO 290
        IOLD = IRN(I0)
        AOLD = A(I0)
        DO 270 IDUM = 1,NE
          IF (I1.EQ.I) THEN
            IRN(I1) = IOLD
            A(I1) = AOLD
            GO TO 290
          ELSE
            INEW = IRN(I1)
            ANEW = A(I1)
            IRN(I1) = IOLD
            A(I1) = AOLD
            IOLD = INEW
            AOLD = ANEW
            I0 = I1
            I1 = DEST(I1)
            DEST(I0) = 0
          ENDIF
  270   CONTINUE
  290 CONTINUE
  299 RETURN
      END
#endif
#endif
      SUBROUTINE ZMUMPS_450(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA
      NVAL = 0 
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END
      SUBROUTINE ZMUMPS_451(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
      INTEGER TODO(TDLEN)
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2
  425   CONTINUE
        IF (TD.EQ.0) GO TO 400 
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
        TD = TD - 2
        GO TO 425
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200 
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE ZMUMPS_452(M,N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4,RLX,RINF)
      IMPLICIT NONE
      INTEGER M,N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(M),IW4(3*N+M)
      DOUBLE PRECISION A(NE),RLX,RINF
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX
      EXTERNAL ZMUMPS_450,ZMUMPS_453,ZMUMPS_455
      DO 20 J = 1,N
        FC(J) = J
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      DO 21 I = 1,M
        IW(I) = 0
   21 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL ZMUMPS_453(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF
      BVAL = 0.0
      BMIN = 0.0
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,M
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX-BMIN .LE. RLX) GO TO 1000
            CALL ZMUMPS_450(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 1000
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60 
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = II - IP(J) + 1
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          IF (BMAX-BMIN .LE. RLX) GO TO 1000
          CALL ZMUMPS_450(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 1000
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL ZMUMPS_453(CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &              IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
   90 CONTINUE 
 1000 IF (M.EQ.N .and. NUMX.EQ.N) GO TO 2000
      CALL ZMUMPS_455(M,N,IPERM,IW,W)
 2000 RETURN
      END
      SUBROUTINE ZMUMPS_453
     &           (ID,MOD,M,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER ID,MOD,M,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),
     &        FC(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1
      IF (ID.EQ.1) THEN
        DO 5 I = 1,M
          CV(I) = 0
    5   CONTINUE
        DO 6 J = 1,N
          ARP(J) = 0
    6   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 J = 1,N
            ARP(J) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N 
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END
      SUBROUTINE ZMUMPS_454(M,N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),OUT(N),PR(N),Q(M),L(M)
      DOUBLE PRECISION A(NE),U(M),D(M),RINF,RINF3
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ,RLX
      LOGICAL LORD
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL ZMUMPS_445,ZMUMPS_446,ZMUMPS_447,ZMUMPS_455
      RLX = U(1)
      RINF3 = U(2)
      LORD = (JPERM(1).EQ.6)
      NUM = 0
      DO 10 K = 1,N
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = RINF
   10 CONTINUE
      DO 15 K = 1,M
        U(K) = RINF3
        IPERM(K) = 0
        L(K) = 0
   15 CONTINUE
      DO 30 J = 1,N
         IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 30
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,M
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IF (JPERM(J).EQ.0) THEN
          JPERM(J) = L(I)
          D(J) = U(I)
          NUM = NUM + 1
        ELSEIF (D(J).GT.U(I)) THEN
          K = JPERM(J)
          II = IRN(K)
          IPERM(II) = 0
          JPERM(J) = L(I)
          D(J) = U(I)
        ELSE
          IPERM(I) = 0
        ENDIF
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 45 K = 1,M
        D(K) = ZERO
   45 CONTINUE
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2) GO TO 95
        VJ = RINF
        DO 50 K = K1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60 
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,M
        D(I) = RINF
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = M + 1
        UP = M + 1
        CSP = RINF
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL ZMUMPS_445(I,M,Q,D,L,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).LT.RINF) DMIN = D(I)*(1+RLX)
            IF (DMIN.GE.CSP) GO TO 160
  152       CALL ZMUMPS_446(QLEN,M,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          IF (DMIN.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          K1 = IP(J+1)-1
          IF (LORD) THEN
            IF (CSP.NE.RINF) THEN
              DI = CSP - VJ
              IF (A(K1).GE.DI) THEN
                K0 = JPERM(J)
                IF (K0.GE.K1-6) GO TO 178
  177           CONTINUE
                  K = (K0+K1)/2
                  IF (A(K).GE.DI) THEN 
                    K1 = K
                  ELSE 
                    K0 = K
                  ENDIF
                  IF (K0.GE.K1-6) GO TO 178
                GO TO 177
  178           DO 179 K = K0+1,K1          
                  IF (A(K).LT.DI) GO TO 179
                  K1 = K - 1
                  GO TO 181
  179           CONTINUE
              ENDIF
            ENDIF
  181       IF (K1.EQ.JPERM(J)) K1 = K1 - 1
          ENDIF
          K0 = IP(J)
          DI = CSP - VJ
          DO 155 K = K0,K1
            I = IRN(K)
            IF (L(I).GE.LOW) GO TO 155
            DNEW = A(K) - U(I)
            IF (DNEW.GE.DI) GO TO 155
            DNEW = DNEW + VJ
            IF (DNEW.GT.D(I)) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
              DI = CSP - VJ
            ELSE
              IF (DNEW.GE.D(I)) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                IF (L(I).NE.0) THEN
                  CALL ZMUMPS_447(L(I),QLEN,M,Q,D,L,2)
                ENDIF
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE   
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL ZMUMPS_445(I,M,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        J = JSP
        IPERM(I) = J
        JPERM(J) = ISP
        DO 170 JDUM = 1,NUM
          JJ = PR(J) 
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 182 KK = UP,M
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  182   CONTINUE 
  190   DO 191 KK = UP,M
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  192   CONTINUE 
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 CONTINUE
      DO 1200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
 1200 CONTINUE
      DO 1201 I = 1,M
        IF (IPERM(I).EQ.0) U(I) = ZERO
 1201 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL ZMUMPS_455(M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE ZMUMPS_457
     &           (M,N,IRN,LIRN,IP,LENC,IPERM,NUM,PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER LIRN,M,N,NUM
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK
      EXTERNAL ZMUMPS_455
      DO 10 I = 1,M
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      DO 12 J = 1,N
        ARP(J) = LENC(J) - 1
   12 CONTINUE
      NUM = 0
      DO 1000 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENC(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 1000
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 1000
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
 1000 CONTINUE
      IF (M.EQ.N .and. NUM.EQ.N) GO TO 2000
      CALL ZMUMPS_455(M,N,IPERM,CV,ARP)
 2000 RETURN
      END
      SUBROUTINE ZMUMPS_455(M,N,IPERM,RW,CW)
      IMPLICIT NONE
      INTEGER M,N
      INTEGER RW(M),CW(N),IPERM(M)
      INTEGER I,J,K
      DO 10 J = 1,N
        CW(J) = 0
   10 CONTINUE
      K = 0
      DO 20 I = 1,M
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          RW(K) = I
        ELSE
          J = IPERM(I)
          CW(J) = I
        ENDIF
   20 CONTINUE
      K = 0
      DO 30 J = 1,N
        IF (CW(J).NE.0) GO TO 30
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   30 CONTINUE
      DO 40 J = N+1,M
        K = K + 1
        I = RW(K)
        IPERM(I) = -J
   40 CONTINUE
      RETURN
      END
      INTEGER FUNCTION ZMUMPS_IZAMAX(n,zx,incx)
      double complex zx(*)
      double precision smax
      double complex zdum
      integer i,incx,ix,n
      ZMUMPS_IZAMAX = 0
      if(n.lt.1)return
      ZMUMPS_IZAMAX = 1
      if (n.eq.1) return
      if(incx.eq.1)go to 20
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      smax = abs(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(zx(ix)).le.smax) go to 5
         ZMUMPS_IZAMAX = i
         smax = abs(zx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 smax = abs(zx(1))
      do 30 i = 2,n
         if(abs(zx(i)).le.smax) go to 30
         ZMUMPS_IZAMAX = i
         smax = abs(zx(i))
   30 continue
      return
      END FUNCTION ZMUMPS_IZAMAX
      SUBROUTINE ZMUMPS_ZSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
      CHARACTER          UPLO
      INTEGER            INCX, LDA, N
      COMPLEX*16         ALPHA
      COMPLEX*16         A( LDA, * ), X( * )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            I, INFO, IX, J, JX, KX
      COMPLEX*16         TEMP
      LOGICAL            ZMUMPS_LSAME
      EXTERNAL           ZMUMPS_LSAME
      INTRINSIC          MAX
      INFO = 0
      IF( .NOT.ZMUMPS_LSAME( UPLO, 'U' ) .AND.
     *    .NOT.ZMUMPS_LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) "Internal error in ZMUMPS_ZSYR"
         CALL ZMUMPS_ABORT()
         RETURN
      END IF
      IF( ( N.EQ.0 ) .OR. ( ALPHA.EQ.ZERO ) )
     $   RETURN
      IF( INCX.LE.0 ) THEN
         KX = 1 - ( N-1 )*INCX
      ELSE IF( INCX.NE.1 ) THEN
         KX = 1
      END IF
      IF( ZMUMPS_LSAME( UPLO, 'U' ) ) THEN
         IF( INCX.EQ.1 ) THEN
            DO 20 J = 1, N
               IF( X( J ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO 10 I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            JX = KX
            DO 40 J = 1, N
               IF( X( JX ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX = KX
                  DO 30 I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
   30             CONTINUE
               END IF
               JX = JX + INCX
   40       CONTINUE
         END IF
      ELSE
         IF( INCX.EQ.1 ) THEN
            DO 60 J = 1, N
               IF( X( J ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO 50 I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            JX = KX
            DO 80 J = 1, N
               IF( X( JX ).NE.ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX = JX
                  DO 70 I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX = IX + INCX
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_ZSYR
      LOGICAL FUNCTION ZMUMPS_LSAME( CA, CB )
      CHARACTER          CA, CB
      INTRINSIC          ICHAR
      INTEGER            INTA, INTB, ZCODE
      ZMUMPS_LSAME = CA.EQ.CB
      IF( ZMUMPS_LSAME )
     $   RETURN
      ZCODE = ICHAR( 'Z' )
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      ZMUMPS_LSAME = INTA.EQ.INTB
      RETURN
      END FUNCTION ZMUMPS_LSAME
      SUBROUTINE ZMUMPS_464( K34, K35 )
      INTEGER K34, K35
#if defined(t3e)
      K34 = 8
      K35 = 32
#else
      K34 = 4
      K35 = 16
#endif
      END SUBROUTINE ZMUMPS_464
      SUBROUTINE ZMUMPS_546( K34, K35 )
      INTEGER K34, K35
#if defined(t3e)
      K34 = 8
      K35 = 16
#else
      K34 = 4
      K35 = 8
#endif
      END SUBROUTINE ZMUMPS_546
