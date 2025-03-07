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
      SUBROUTINE CMUMPS_448 (ICNTL,CNTL)
      IMPLICIT NONE
      INTEGER NICNTL, NCNTL
      PARAMETER (NICNTL=10, NCNTL=10)
      INTEGER ICNTL(NICNTL)
      REAL CNTL(NCNTL)
      INTEGER I
      ICNTL(1) =  6
      ICNTL(2) =  6
      ICNTL(3) = -1
      ICNTL(4) = -1
      ICNTL(5) =  0
      DO 10 I = 6,NICNTL
        ICNTL(I) = 0
   10 CONTINUE
      CNTL(1) = 0.0
      CNTL(2) = 0.0
      DO 20 I = 3,NCNTL
        CNTL(I) = 0.0
   20 CONTINUE
      RETURN
      END
      SUBROUTINE CMUMPS_444
     &           (M,N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),PR(N),Q(M),L(M)
      REAL A(NE)
      REAL D(M), RINF
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW
      REAL             CSP,DI,DNEW,DQ0,AI,A0,BV,TBV,RLX
      REAL             ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC ABS,MIN
      EXTERNAL CMUMPS_445, CMUMPS_446, CMUMPS_447, CMUMPS_455
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
              CALL CMUMPS_445 (I,M,Q,D,L,1)
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
              CALL CMUMPS_446 (QLEN,M,Q,D,L,1)
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
                  CALL CMUMPS_447 (L(I),QLEN,M,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL CMUMPS_445 (I,M,Q,D,L,1)
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
      CALL CMUMPS_455 (M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE CMUMPS_445 (I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DI
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
      SUBROUTINE CMUMPS_446 (QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      REAL DK,DR,DI
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
      SUBROUTINE CMUMPS_447 (POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      REAL D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      REAL DK,DR,DI
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
      SUBROUTINE CMUMPS_458 (M,N,NE,IP,IRN,A,MODE,PERM,IW,INFO)
      IMPLICIT NONE
      INTEGER M,N,NE,MODE,INFO
      REAL A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(*),IW(*)
      INTEGER I, K
      EXTERNAL CMUMPS_456,CMUMPS_449
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
        CALL CMUMPS_456 (M,PERM,IW,INFO)
      ELSE
        CALL CMUMPS_456 (N,PERM,IW,INFO)
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
      CALL CMUMPS_449 (M,N,NE,IP,IRN,A,PERM,IW(1),IW(N+1))
  299 RETURN
      END
      SUBROUTINE CMUMPS_449 (M,N,NE,IP,IRN,A,PERM,IW,DEST)
      IMPLICIT NONE
      INTEGER M,N,NE
      REAL A(NE)
      INTEGER IP(N+1),IRN(NE),PERM(N),IW(N),DEST(NE)
      REAL AOLD,ANEW
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
      SUBROUTINE CMUMPS_450 (IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      REAL A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      REAL SPLIT(XX),HA
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
      SUBROUTINE CMUMPS_451 (N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      REAL A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      REAL HA,KEY
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
      SUBROUTINE CMUMPS_452 (M,N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4,RLX,RINF)
      IMPLICIT NONE
      INTEGER M,N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N), 
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(M),IW4(3*N+M)
      REAL A(NE),RLX,RINF
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      REAL             BVAL,BMIN,BMAX
      EXTERNAL CMUMPS_450,CMUMPS_453,CMUMPS_455
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
      CALL CMUMPS_453 (CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
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
            CALL CMUMPS_450 (IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
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
          CALL CMUMPS_450 (IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
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
        CALL CMUMPS_453 (CNT,MOD,M,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &              IW4(1),IW4(N+1),IW4(2*N+1),IW4(2*N+M+1))
   90 CONTINUE 
 1000 IF (M.EQ.N .and. NUMX.EQ.N) GO TO 2000
      CALL CMUMPS_455 (M,N,IPERM,IW,W)
 2000 RETURN
      END
      SUBROUTINE CMUMPS_453
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
      SUBROUTINE CMUMPS_454 (M,N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D,RINF) 
      IMPLICIT NONE
      INTEGER M,N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(M),JPERM(N),OUT(N),PR(N),Q(M),L(M)
      REAL A(NE),U(M),D(M),RINF,RINF3
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW
      REAL             CSP,DI,DMIN,DNEW,DQ0,VJ,RLX
      LOGICAL LORD
      REAL             ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL CMUMPS_445, CMUMPS_446, CMUMPS_447, CMUMPS_455
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
            CALL CMUMPS_445 (I,M,Q,D,L,2)
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
  152       CALL CMUMPS_446 (QLEN,M,Q,D,L,2)
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
                  CALL CMUMPS_447 (L(I),QLEN,M,Q,D,L,2)
                ENDIF
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE   
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL CMUMPS_445 (I,M,Q,D,L,2)
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
      CALL CMUMPS_455 (M,N,IPERM,L,JPERM)
 2000 RETURN
      END
      SUBROUTINE CMUMPS_457
     &           (M,N,IRN,LIRN,IP,LENC,IPERM,NUM,PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER LIRN,M,N,NUM
      INTEGER ARP(N),CV(M),IRN(LIRN),IP(N),IPERM(M),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK
      EXTERNAL CMUMPS_455
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
      CALL CMUMPS_455 (M,N,IPERM,CV,ARP)
 2000 RETURN
      END
      SUBROUTINE CMUMPS_455 (M,N,IPERM,RW,CW)
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
      INTEGER FUNCTION CMUMPS_ICAMAX(n,cx,incx)
      complex cx(*)
      real smax
      integer i,incx,ix,n
      complex zdum
      CMUMPS_ICAMAX = 0
      if( n .lt. 1 ) return
      CMUMPS_ICAMAX = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
      ix = 1
      smax = abs(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(cx(ix)).le.smax) go to 5
         CMUMPS_ICAMAX = i
         smax = abs(cx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 smax = abs(cx(1))
      do 30 i = 2,n
         if(abs(cx(i)).le.smax) go to 30
         CMUMPS_ICAMAX = i
         smax = abs(cx(i))
   30 continue
      return
      END FUNCTION CMUMPS_ICAMAX
      SUBROUTINE CMUMPS_CSYR( UPLO, N, ALPHA, X, INCX, A, LDA )
      CHARACTER          UPLO
      INTEGER            INCX, LDA, N
      COMPLEX            ALPHA
      COMPLEX            A( LDA, * ), X( * )
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            I, INFO, IX, J, JX, KX
      COMPLEX            TEMP
      LOGICAL            CMUMPS_LSAME
      EXTERNAL           CMUMPS_LSAME
      INTRINSIC          MAX
      INFO = 0
      IF( .NOT.CMUMPS_LSAME( UPLO, 'U' ) .AND.
     *    .NOT.CMUMPS_LSAME( UPLO, 'L' ) ) THEN
         INFO = 1
      ELSE IF( N.LT.0 ) THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = 7
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(*,*) "Internal error in CMUMPS_CSYR"
         CALL CMUMPS_ABORT()
         RETURN
      END IF
      IF( ( N.EQ.0 ) .OR. ( ALPHA.EQ.ZERO ) )
     $   RETURN
      IF( INCX.LE.0 ) THEN
         KX = 1 - ( N-1 )*INCX
      ELSE IF( INCX.NE.1 ) THEN
         KX = 1
      END IF
      IF( CMUMPS_LSAME( UPLO, 'U' ) ) THEN
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
      END SUBROUTINE CMUMPS_CSYR
      LOGICAL FUNCTION CMUMPS_LSAME( CA, CB )
      CHARACTER          CA, CB
      INTRINSIC          ICHAR
      INTEGER            INTA, INTB, ZCODE
      CMUMPS_LSAME = CA.EQ.CB
      IF( CMUMPS_LSAME )
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
      CMUMPS_LSAME = INTA.EQ.INTB
      RETURN
      END FUNCTION CMUMPS_LSAME
      SUBROUTINE CMUMPS_464( K34, K35 )
      INTEGER K34, K35
#if defined(t3e)
      K34 = 8
      K35 = 16
#else
      K34 = 4
      K35 = 8
#endif
      END SUBROUTINE CMUMPS_464
      SUBROUTINE CMUMPS_546( K34, K35 )
      INTEGER K34, K35
#if defined(t3e)
      K34 = 8
      K35 = 16
#else
      K34 = 4
      K35 = 8
#endif
      END SUBROUTINE CMUMPS_546
