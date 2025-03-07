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
      SUBROUTINE SMUMPS_246(MYID, N, STEP, FRERE, FILS,
     *     NA, LNA, NE, DAD, ND, PROCNODE, SLAVEF,
     *     NRLADU, NIRADU, NIRNEC, NRLNEC,
     *     NRLNEC_ACTIVE,
     *     MAXFR, OPSA,
     *     KEEP,KEEP8, LOCAL_M, LOCAL_N,
     *     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     *     I_AM_CAND,NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     *     IFLAG, IERROR
     $     ,MAX_FRONT_SURFACE_LOCAL
     *     ,MAX_SIZE_FACTOR
     $     )
      IMPLICIT NONE
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NRLADU, NIRADU, NRLNEC, NIRNEC
      INTEGER NRLADU_CURRENT, NRLADU_ROOT_3, NRLNEC_ACTIVE
      INTEGER MAXFR, NSTEPS
      INTEGER MAX_FRONT_SURFACE_LOCAL
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), NA(LNA), NE(NSTEPS),
     *        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER*8 KEEP8(150)
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      LOGICAL  I_AM_CAND(NMB_PAR2)
      INTEGER  CANDIDATES( SLAVEF+1, NMB_PAR2 )
      INTEGER  PHASE
      PARAMETER (PHASE=0)
      REAL OPSA
      DOUBLE PRECISION OPSA_LOC !Local variable
      INTEGER MAX_SIZE_FACTOR
      REAL OPS_SUBTREE
      DOUBLE PRECISION OPS_SBTR_LOC ! Local variable
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSTKR, TNSTK, IPOOL, 
     *                                      LSTKI 
      INTEGER BLOCKING_RHS
      INTEGER I,ISTKR,ITOP,NELIM,NFR
      INTEGER ISTKI, ISTKIM, STKI
      INTEGER K,LSTK,NSTK, IFATH
      INTEGER INODE, LEAF, NBLEAF, NBROOT, IN
      INTEGER LEVEL, MAXTEMPCB, MAXITEMPCB
      INTEGER CURRENT_ACTIVE_MEM
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF 
      INTEGER LEVELF, NCB, SIZECB, SIZECBI, SIZECBINFR
      INTEGER SIZECB_SLAVE, SIZEHEADER
      INTEGER NBROWMAX, NSLAVES,  NSLAVES_LOC, NSLAVES_PASSED,
     *         NELIMF, NFRF, NCBF, CBMAXR, CBMAXS,
     *         NSLAVESF, NBROWMAXF, LKJIB,
     *         LKJIBT, NBR, NBCOLFAC
      INTEGER LEV3MAXREC
      INTEGER LWK_RR, LIWK_RR
      INTEGER IROOT, SIZE_ROOT
      INTEGER ALLOCOK
      LOGICAL ROOT_OWNER, COMPRESSCB
      DOUBLE PRECISION OPS_NODE, OPS_NODE_MASTER, OPS_NODE_SLAVE
      INCLUDE 'mumps_headers.h'
      INTEGER WHAT
      INTEGER IDUMMY
      INTRINSIC MIN, INT, REAL
      INTEGER SMUMPS_275, SMUMPS_330
      LOGICAL SMUMPS_170
      INTEGER SMUMPS_52
      EXTERNAL SMUMPS_503, SMUMPS_52
      EXTERNAL SMUMPS_275, SMUMPS_330, 
     *         SMUMPS_170
      logical :: FORCE_CAND, CONCERNED, UPDATES, STACKCB, MASTERSON
      integer :: istat, IFSON, LEVELSON
      COMPRESSCB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      MAX_FRONT_SURFACE_LOCAL=0
      MAX_SIZE_FACTOR=0
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), IPOOL(NSTEPS),
     *          LSTKI(NSTEPS) , stat=ALLOCOK)
      if (ALLOCOK .GT. 0) THEN
        IFLAG  =-7
        IERROR = 4*NSTEPS
        RETURN
      endif
      LKJIB = MAX(KEEP(5),KEEP(6))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
      SIZEHEADER = XSIZE+6  ! default size of header
      ISTKR  = 0
      ISTKI  = 0
      ISTKIM = 0
      OPSA_LOC   = DBLE(0.0E0)
      OPS_SBTR_LOC = DBLE(0.0E0)
      NRLADU = 0
      NIRADU = 0
      NRLADU_CURRENT = 0
      NRLADU_ROOT_3 = 0
      NRLNEC_ACTIVE = 0
      NRLNEC = 0
      NIRNEC = 0
      MAXFR  = 0
      ITOP = 0
      MAXTEMPCB  = 0
      MAXITEMPCB = 0
      SBUF_SEND  = 1
      SBUF_REC   = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU = LOCAL_M*LOCAL_N
        NRLADU_ROOT_3 = LOCAL_M*LOCAL_N
        NRLNEC_ACTIVE = NRLADU_CURRENT
        MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,LOCAL_M*LOCAL_N)
        NRLNEC = NRLADU
        IF (SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *                                       .EQ. MYID) THEN
          NIRADU = SIZEHEADER+2*ND(STEP(INODE))
        ELSE
          NIRADU = SIZEHEADER
        ENDIF
        NIRNEC = NIRADU
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         FORCE_CAND=.FALSE.           
      ELSE
         FORCE_CAND=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE 
         WRITE(MYID+6,*) ' ERROR 1 in SMUMPS_246 '
         CALL SMUMPS_ABORT()
      ENDIF
 95   CONTINUE 
      NFR    = ND(STEP(INODE))
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IFATH = DAD(STEP(INODE))
      MASTER = SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
     *           .EQ. MYID
      LEVEL  = SMUMPS_330(STEP(INODE),PROCNODE,SLAVEF)
      UPDATE=.FALSE.
       if(.NOT.FORCE_CAND) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
       else
         if(MASTER.and.(LEVEL.ne.3)) then
            UPDATE = .TRUE.
         else if(LEVEL.eq.2) then
            if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(INODE)))) THEN
              UPDATE = .TRUE.
            end if
         end if
       end if
      NCB      = NFR-NELIM
      SIZECBINFR = NCB*NCB
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IFATH = DAD(STEP(INODE))
        IF ( IFATH.NE.KEEP(38) .AND. COMPRESSCB ) THEN
          SIZECB    = (NCB*(NCB+1))/2
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      SIZECBI  = 2* NCB  + SIZEHEADER
      IF (LEVEL.NE.2) THEN
        NSLAVES_LOC     = -99999999
        SIZECB_SLAVE = -99999997
        NBROWMAX        = NCB
      ELSE
        IF (KEEP(48) .EQ. 5) THEN
          WHAT = 5 ! Compute both NBROWMAX and SIZECB_SLAVE
          IF (FORCE_CAND) THEN
            NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $                    ISTEP_TO_INIV2(STEP(INODE)))
          ELSE
            NSLAVES_LOC=SLAVEF-1
          ENDIF
          NSLAVES_PASSED=NSLAVES_LOC
        ELSE
          WHAT = 2 ! Compute both NBROWMAX and SIZECB_SLAVE
          NSLAVES_PASSED=SLAVEF
          NSLAVES_LOC   =SLAVEF-1
        ENDIF
         CALL SMUMPS_503(WHAT, KEEP,KEEP8,
     &     NCB, NFR, NSLAVES_PASSED, NBROWMAX, SIZECB_SLAVE)
      ENDIF
      IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU = NIRADU+SIZEHEADER+2*ND(STEP(INODE))
         ENDIF
      ENDIF
      IF (LEVEL.EQ.3) THEN
         IF ( 
     *     KEEP(60).LE.1
     *      ) THEN
           NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+
     *                 LOCAL_M*LOCAL_N)
           NRLADU_CURRENT = LOCAL_M*LOCAL_N
           NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     *                        NRLADU_CURRENT+ISTKR)
         ENDIF
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.SMUMPS_170(STEP(INODE),
     $        PROCNODE,SLAVEF)))THEN
            IF(LEVEL.EQ.1)THEN
               MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $              NFR*NFR)
            ELSEIF(LEVEL.EQ.2)THEN
               IF(KEEP(50).EQ.0)THEN
                 MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                 NFR*NELIM)
               ELSE
                 MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                 NELIM*NELIM)
                 IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                  MAX_FRONT_SURFACE_LOCAL=MAX(MAX_FRONT_SURFACE_LOCAL,
     $                  NELIM*(NELIM+1))
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF (LEVEL.EQ.2) THEN
        IF (MASTER) THEN
          IF (KEEP(50).EQ.0) THEN
             SBUF_SEND = MAX(SBUF_SEND, NFR*LKJIB+LKJIB+4)
          ELSE
             SBUF_SEND = MAX(SBUF_SEND, NELIM*LKJIB+NELIM+6)
          ENDIF
        ELSEIF (UPDATE) THEN
            if (KEEP(50).EQ.0) THEN
              SBUF_REC   = MAX(SBUF_REC, NFR*LKJIB+LKJIB+4)
            else
              SBUF_REC = MAX( SBUF_REC, NELIM*LKJIB+NELIM+6 )
              IF (KEEP(50).EQ.1) THEN
                LKJIBT  = LKJIB
              ELSE
                LKJIBT  = MIN( NELIM, LKJIB * 2 )
              ENDIF
              SBUF_SEND = MAX(SBUF_SEND,
     *                        LKJIBT*NBROWMAX+6)
              SBUF_REC = MAX( SBUF_REC, NBROWMAX*LKJIBT+6 )
            endif
        ENDIF
      ENDIF
      IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            NIRADU = NIRADU + 2*NFR + SIZEHEADER
            IF (KEEP(50).EQ.0) THEN
             NRLADU = NRLADU + NELIM*(2*NFR-NELIM)
             NRLADU_CURRENT = NELIM*(2*NFR-NELIM)
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*(2*NFR-NELIM))
            ELSE
             NRLADU = NRLADU + NELIM*NFR
             NRLADU_CURRENT = NELIM*NFR
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*NFR)
            ENDIF
            SIZECBI    = 2* NCB  + 6 + 3
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
              NIRADU = NIRADU+SIZEHEADER +SLAVEF-1+2*NFR 
              IF (KEEP(50).EQ.0) THEN
                NBCOLFAC=NFR
              ELSE
                NBCOLFAC=NELIM
              ENDIF
              NRLADU = NRLADU + NBCOLFAC*NELIM
              NRLADU_CURRENT = NBCOLFAC*NELIM
              MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NBCOLFAC*NELIM)
               SIZECB     = 0
               SIZECBINFR = 0
               SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
             SIZECB=SIZECB_SLAVE
             SIZECBINFR = SIZECB
             NIRADU   = NIRADU+4+NELIM+NBROWMAX
             NRLADU   = NRLADU + NELIM*NBROWMAX
             NRLADU_CURRENT = NELIM*NBROWMAX
             MAX_SIZE_FACTOR=MAX(MAX_SIZE_FACTOR,NELIM*NBROWMAX)
             SIZECBI      = 4 + NBROWMAX + NCB
             IF (KEEP(50).NE.0) SIZECBI=SIZECBI+NSLAVES_LOC+
     *                                  XTRA_SLAVES_SYM
            ENDIF
         ENDIF
         NIRNEC = MAX0(NIRNEC,
     *             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
         IF (KEEP(50).NE.0.AND.UPDATE.AND.LEVEL.EQ.1) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + NELIM*NCB
         ENDIF
         IF (MASTER .AND.  KEEP(219).NE.0.AND.
     *       KEEP(50).EQ.2.AND.LEVEL.EQ.2) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + NELIM
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC = MAX0(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
           NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_CURRENT+
     *             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ELSE
           NRLNEC = MAX0(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_CURRENT+
     *             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+MAXTEMPCB)
         ENDIF
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     *            ': ERROR 2 in SMUMPS_246. ITOP = ',ITOP
                  CALL SMUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF ! NSTK .GT. 0
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = SMUMPS_275(STEP(IFSON),PROCNODE,SLAVEF)
     *                  .EQ.MYID
            LEVELSON  = SMUMPS_330(STEP(IFSON),PROCNODE,SLAVEF)
            if(.NOT.FORCE_CAND) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     *                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFSON)))) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,
     *          ': ERROR 2 in SMUMPS_246. ITOP = ',ITOP
                CALL SMUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) ! process next son
         END DO
      ENDIF
      IF (UPDATE.OR.LEVEL.EQ.3) THEN
         IF (
     *        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     *       .AND.
     *        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     *      )
     *   THEN
            CALL SMUMPS_511(NFR, NELIM, NELIM,KEEP(50),
     *           1,OPS_NODE)
            IF (LEVEL.EQ.2) THEN
              CALL SMUMPS_511(NFR, NELIM, NELIM,KEEP(50),
     *           2,OPS_NODE_MASTER)
              OPS_NODE_SLAVE=OPS_NODE-OPS_NODE_MASTER
            ENDIF
         ELSE
           OPS_NODE = 0.0E0
         ENDIF
         IF ( LEVEL .EQ. 3 ) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE / DBLE( SLAVEF )
         ELSE IF (MASTER .AND. LEVEL.EQ.2) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE_MASTER
         ELSE IF (MASTER .AND. LEVEL.EQ.1) THEN
            OPSA_LOC = OPSA_LOC + DBLE(OPS_NODE)
         ELSE IF (UPDATE) THEN ! Slave task
            OPSA_LOC = OPSA_LOC + 
     &            DBLE(OPS_NODE_SLAVE)/DBLE(NSLAVES_LOC)
         ENDIF
         IF (SMUMPS_170(STEP(INODE),
     *   PROCNODE, SLAVEF) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ELSE
             CALL SMUMPS_511(NFR, NELIM, NELIM,KEEP(50),
     *           1,OPS_NODE)
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ENDIF
         ENDIF
        ENDIF
      IF (IFATH .EQ. 0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
      ELSE
         NFRF = ND(STEP(IFATH))
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         NCBF = NFRF - NELIMF
         LEVELF =SMUMPS_330(STEP(IFATH),PROCNODE,SLAVEF)
         MASTERF=SMUMPS_275(STEP(IFATH),PROCNODE,SLAVEF).EQ.MYID
         UPDATEF=.FALSE.
         if(.NOT.FORCE_CAND) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(MASTERF.and.(LEVELF.ne.3)) then
               UPDATEF = .TRUE.
            else if (LEVELF.eq.2) then
               if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(ifath)))) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         CONCERNED  = UPDATEF .OR. UPDATE
         IF (LEVELF .NE. 2) THEN
           NBROWMAXF = -999999
         ELSE
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     $               ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
           ELSE
               WHAT = 1 ! NBROWMAX only
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL SMUMPS_503( WHAT, KEEP, KEEP8,
     &     NCBF, NFRF, NSLAVES_LOC, NBROWMAXF, IDUMMY )
         ENDIF
         IF(LEVEL.EQ.1.AND.UPDATE.AND.
     *      (UPDATEF.OR.LEVELF.EQ.2)
     *      .AND.LEVELF.NE.3) THEN
             NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,NRLADU_CURRENT+
     *         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
             NRLNEC = MAX0(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
         ENDIF
         IF (UPDATE .AND. LEVEL.EQ.2 .AND. .NOT. MASTER) THEN
             NRLNEC =
     *         MAX0(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
             NRLNEC_ACTIVE = MAX0(NRLNEC_ACTIVE,2*NRLADU_CURRENT+
     *         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ENDIF
        IF (LEVELF.EQ.3) THEN
          IF (LEVEL.EQ.1) THEN
            LEV3MAXREC = MIN(NCB,LOCAL_M) * MIN(NCB,LOCAL_N)
          ELSE
            LEV3MAXREC = MIN(SIZECB,
     *                 MIN(NBROWMAX,LOCAL_M)*MIN(NCB,LOCAL_N)) 
          ENDIF
          MAXTEMPCB  = MAX(MAXTEMPCB, LEV3MAXREC)
          MAXITEMPCB = MAX(MAXITEMPCB,SIZECBI+SIZEHEADER)
          SBUF_REC   = MAX(SBUF_REC, LEV3MAXREC+SIZECBI)
        ENDIF
        IF (CONCERNED) THEN
         IF (LEVELF.EQ.2) THEN
           IF (UPDATE.AND.(LEVEL.NE.2.OR..NOT.MASTER)) THEN
             IF(MASTERF)THEN
                 NBR = MIN(NBROWMAXF,NBROWMAX)
             ELSE
                 NBR = MIN(MAX(NELIMF,NBROWMAXF),NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               CBMAXS = NBR*NCB
             ELSE
               CBMAXS = NBR*NCB - (NBR*(NBR-1))/2
             ENDIF
           ELSE
              CBMAXS = 0
           END IF
           IF (MASTERF) THEN
             IF (LEVEL.EQ.1) THEN
                IF (.NOT.UPDATE) THEN
                  NBR = MIN(NELIMF, NCB)
                ELSE
                  NBR = 0
                ENDIF
             ELSE
                NBR = MIN(NELIMF, NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
                CBMAXR = NBR*NCB
             ELSE
                CBMAXR = NBR*MIN(NCB,NELIMF)-(NBR*(NBR-1))/2
                CBMAXR = MIN(CBMAXR, NELIMF*(NELIMF+1)/2 )
                CBMAXR = MIN(CBMAXR, SIZECB)
                IF ((LEVEL.EQ.1).AND.(.NOT. COMPRESSCB)) THEN
                  CBMAXR = MIN(CBMAXR,(NCB*(NCB+1))/2)
                ENDIF
             ENDIF
           ELSE IF (UPDATEF) THEN
              NBR = MIN(NBROWMAXF,NBROWMAX)
              CBMAXR = NBR * NCB
              IF (KEEP(50).NE.0) THEN
                CBMAXR = CBMAXR - (NBR*(NBR-1))/2
              ENDIF
           ELSE
              CBMAXR = 0
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LEV3MAXREC
           IF (UPDATE.AND. .NOT. (MASTER.AND.LEVEL.EQ.2)) THEN
             CBMAXS = LEV3MAXREC
           ELSE
             CBMAXS = 0
           ENDIF
         ELSE
           IF (MASTERF) THEN
             CBMAXS = 0
             NBR = MIN(NFRF,NBROWMAX)
             IF ((LEVEL.EQ.1).AND.UPDATE) THEN
                NBR = 0
             ENDIF
             CBMAXR = NBR*MIN(NFRF,NCB)
             IF (LEVEL.EQ.2)
     *       CBMAXR = MIN(CBMAXR, SIZECB_SLAVE)
             IF ( KEEP(50).NE.0 )  THEN
              CBMAXR = MIN(CBMAXR,(NFRF*(NFRF+1))/2)
             ELSE
              CBMAXR = MIN(CBMAXR,NFRF*NFRF)
             ENDIF
           ELSE
             CBMAXR = 0
             CBMAXS = SIZECB
           ENDIF
         ENDIF
         IF (UPDATE) THEN
           CBMAXS = MIN(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUF_SEND = MAX(SBUF_SEND, CBMAXS+SIZECBI)
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN
          STACKCB = .TRUE.
          SIZECBI = 2 * NFR + SIZEHEADER
          IF (LEVEL.EQ.1) THEN
             IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     *           .AND.COMPRESSCB) THEN
                 SIZECB = (NCB*(NCB+1))/2
             ELSE
                 SIZECB = NCB*NCB
             ENDIF
             IF (MASTER) THEN
               SIZECBI = 2+ XSIZE
             ELSE IF (LEVELF.EQ.1) THEN
               SIZECB  = MIN(CBMAXR,SIZECB)
               SIZECBI    = 2 * NCB +  9 ! what is sent
               SBUF_REC   = MAX(SBUF_REC, SIZECBI+SIZECB)
               SIZECBI    =  2 * NCB + SIZEHEADER ! what is stacked
             ELSE 
               SIZECBI    = 2 * NCB +  9 ! what is sent
               SBUF_REC   = MAX(SBUF_REC, 
     &                      MIN(SIZECB,CBMAXR) + SIZECBI)
               MAXTEMPCB  = MAX(MAXTEMPCB, MIN(SIZECB,CBMAXR)) 
               SIZECBI    =  2 * NCB + SIZEHEADER ! what is stacked
               MAXITEMPCB = MAX(MAXITEMPCB, SIZECBI)
               SIZECBI= 0
               SIZECB = 0   
             ENDIF
          ELSE ! (of LEVEL.EQ.1)
             SIZECB = SIZECB_SLAVE
             MAXTEMPCB  = MAX(MAXTEMPCB, MIN(CBMAXR,SIZECB) )
             MAXITEMPCB = MAX(MAXITEMPCB,NBROWMAX+NCB+SIZEHEADER)
             IF (.NOT. 
     &        (UPDATE.AND.(.NOT.MASTER).AND.(NSLAVES_LOC.EQ.1))
     &          ) 
     &       SBUF_REC = MAX(SBUF_REC, 
     &            MIN(CBMAXR,SIZECB) +NBROWMAX + NCB + 6)
             IF (MASTER) THEN
              SIZECBI =  NCB + 5 +  SLAVEF - 1 + XSIZE
              SIZECB  = 0
             ELSE IF (UPDATE) THEN
              SIZECBI =  NFR + 6 + SLAVEF - 1 + XSIZE
              IF (KEEP(50).EQ.0) THEN
                SIZECBI = SIZECBI + NBROWMAX + NFR + 
     &                    SIZEHEADER
              ELSE
                SIZECBI = SIZECBI + NBROWMAX + NFR +
     &                    SIZEHEADER+ NSLAVES_LOC
              ENDIF
             ELSE
              SIZECB  = 0 
              SIZECBI = 0
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB = .TRUE.
               SIZECB  = 0 
               SIZECBI = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (COMPRESSCB) THEN 
                      SIZECB  = (NCB*(NCB+1))/2
                  ELSE
                      SIZECB  = NCB*NCB
                  ENDIF
                  SIZECBI = 2 * NCB + SIZEHEADER
               ELSE IF (LEVEL.EQ.2) THEN
                 IF (MASTER) THEN
                   SIZECBI =  NCB + 5 +  SLAVEF - 1 + XSIZE
                 ELSE 
                   SIZECB  = SIZECB_SLAVE
                   SIZECBI = SIZECBI + NBROWMAX + NFR + SIZEHEADER
                 ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
           IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in SMUMPS_246'
                  CALL SMUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in SMUMPS_246 '
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI+LSTKI(ITOP)
           ISTKIM = MAX0(ISTKI,ISTKIM)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = MAX0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
         ENDIF ! end of (stackcb)
        ENDIF ! end of (CONCERNED) 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF ! IFATH .EQ. 0
 115  CONTINUE
#if defined(try_null_space)
      IF ( KEEP(53) .NE. 0 ) THEN
        IF ( KEEP(38) .ne. 0 ) THEN
          IROOT = KEEP( 38 )
        ELSE
          IROOT = KEEP( 20 )
        END IF
        ROOT_OWNER = ( MYID .eq.
     *  SMUMPS_275( STEP(IROOT), PROCNODE, SLAVEF ) )
        SIZE_ROOT = ND(STEP(IROOT))
        CALL SMUMPS_279( PHASE, 
     *       KEEP(51), KEEP(51), SIZE_ROOT,
     *       LOCAL_M, LOCAL_N, ROOT_OWNER, KEEP,KEEP8,
     *       LIWK_RR, LWK_RR )
        IF ( NRLNEC-NRLADU .LT. LWK_RR ) THEN
          NRLNEC = NRLADU + LWK_RR
        END IF
        IF ( NIRNEC-NIRADU .LT. LIWK_RR ) THEN
          NIRNEC = NIRADU + LIWK_RR
        END IF
      END IF
#endif
      NRLNEC = MAX(NRLNEC, NRLADU+4*KEEP(127)*ABS(KEEP(84)))
      IF (KEEP(84) .LT. 0) THEN
        BLOCKING_RHS = - 2 * KEEP(84)
      ELSE
        BLOCKING_RHS = KEEP(84)
      ENDIF
      NRLNEC_ACTIVE = MAX(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     *                    4*KEEP(127)*BLOCKING_RHS)
      SBUF_REC = MAX(SBUF_REC, MAXTEMPCB+MAXITEMPCB )
      SBUF_REC = SBUF_REC   + 10
      SBUF_SEND = SBUF_SEND + 10
      IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
         SBUF_REC = SBUF_REC+KEEP(108)+1
         SBUF_SEND = SBUF_SEND+KEEP(108)+1
      ENDIF
      IF (SLAVEF.EQ.1) THEN 
         SBUF_REC = 1
         SBUF_SEND= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     *          LSTKI )
      OPS_SUBTREE = REAL(OPS_SBTR_LOC)
      OPSA        = REAL(OPSA_LOC)
      RETURN
      END SUBROUTINE SMUMPS_246
      SUBROUTINE SMUMPS_271( COMM_LOAD, ASS_IRECV, 
     *    INODE, NELIM_ROOT, ROOT, 
     *
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
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (SMUMPS_ROOT_STRUC) :: ROOT
      INTEGER KEEP(500), ICNTL( 40 )
      INTEGER*8 KEEP8(150)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER INODE, NELIM_ROOT
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER POSFAC,IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28)),
     *PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS(KEEP(28))
      INTEGER IFLAG, IERROR, COMM
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC( N ), FILS( N )
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) ) 
      INTEGER INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      INCLUDE 'mumps_tags.h'
      INTEGER I, J, OPSFAC, APOS, LCONT, NCOL_TO_SEND, LDA
      INTEGER FPERE, IOLDPS, NFRONT, NPIV, NASS, NSLAVES,
     *        H_INODE, NELIM, NBCOL, LIST_NELIM_ROW, 
     *        LIST_NELIM_COL, NELIM_LOCAL, TYPE_SON, 
     *        POSELT, NROW, NCOL, NBROW, SHIFT_LIST_ROW_SON,
     *        SHIFT_LIST_COL_SON, SHIFT_VAL_SON,LDAFS, IERR,
     *        STATUS( MPI_STATUS_SIZE ), ISON, PDEST_MASTER_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER MSGSOU, MSGTAG
      LOGICAL INVERT, FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER  SMUMPS_275, SMUMPS_330
      EXTERNAL SMUMPS_275, SMUMPS_330
      FPERE = KEEP(38)
      TYPE_SON = SMUMPS_330(STEP(INODE),PROCNODE_STEPS,SLAVEF)
      IF ( SMUMPS_275( STEP(INODE), PROCNODE_STEPS,
     *     SLAVEF ).EQ.MYID) THEN
       IOLDPS   = PTLUST_S(STEP(INODE))
       NFRONT   = IW(IOLDPS+XSIZE)
       NPIV     = IW(IOLDPS+1+XSIZE)
       NASS     = IABS(IW(IOLDPS + 2+XSIZE))
       NSLAVES  =  IW(IOLDPS+5+XSIZE)
       H_INODE  = 6 + NSLAVES + XSIZE
       NELIM    = NASS - NPIV
       NBCOL = NFRONT - NPIV
       LIST_NELIM_ROW = IOLDPS + H_INODE + NPIV
       LIST_NELIM_COL = LIST_NELIM_ROW + NFRONT
           IF (NELIM.LE.0) THEN
            write(6,*) ' ERROR 1 in SMUMPS_271 ', NELIM
            write(6,*) MYID,':Process root2son: INODE=',INODE,
     * 'Header=',IW(PTLUST_S(STEP(INODE)):PTLUST_S(STEP(INODE))+5+XSIZE)
            CALL SMUMPS_ABORT()
           ENDIF
       NELIM_LOCAL = NELIM_ROOT
       DO I=1, NELIM
        root%RG2L_ROW(IW(LIST_NELIM_ROW)) = NELIM_LOCAL
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_ROW = LIST_NELIM_ROW + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       NBROW = NFRONT - NPIV
       NROW = NELIM
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         NCOL = NFRONT - NPIV
       ELSE
         NCOL = NELIM
       END IF
       SHIFT_LIST_ROW_SON = H_INODE + NPIV
       SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
       IF ( KEEP(50).eq.0 .OR. TYPE_SON .eq. 1 ) THEN
         LDAFS = NFRONT
       ELSE
         LDAFS = NASS
       END IF
       SHIFT_VAL_SON      = NPIV * LDAFS + NPIV
       CALL SMUMPS_80( COMM_LOAD,
     *   ASS_IRECV, 
     *   N, INODE, FPERE,
     *   PTLUST_S(1), PTRAST(1),
     *   ROOT, NROW, NCOL, SHIFT_LIST_ROW_SON,
     *   SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDAFS,
     *   ROOT_NON_ELIM_CB, MYID, COMM,
     *   BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *   IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *   NIRBDU, PTRIST, PTLUST_S(1), PTRFAC(1), PTRAST(1),
     *   STEP, PIMASTER, PAMASTER,
     *   NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *   IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *   OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *   INTARR, DBLARR, ICNTL, KEEP,KEEP8, .FALSE., ND, FRERE,
     *   LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
       IF (IFLAG.LT.0 ) RETURN
       IF (TYPE_SON.EQ.1) THEN
        NROW = NFRONT - NASS
        NCOL = NELIM
        SHIFT_LIST_ROW_SON = H_INODE + NASS
        SHIFT_LIST_COL_SON = H_INODE + NFRONT + NPIV
        SHIFT_VAL_SON      = NASS * NFRONT + NPIV
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          INVERT = .FALSE.
        ELSE
          INVERT = .TRUE.
        END IF
        CALL SMUMPS_80( COMM_LOAD, ASS_IRECV,
     *    N, INODE, FPERE,
     *    PTLUST_S, PTRAST,
     *    ROOT, NROW, NCOL, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, NFRONT,
     *    ROOT_NON_ELIM_CB, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
       ENDIF
       IOLDPS = PTLUST_S(STEP(INODE))
       POSELT = PTRAST(STEP(INODE))
       IW(IOLDPS + 4+XSIZE) = STEP(INODE)
       PTRFAC(STEP(INODE))=POSELT
       IF ( TYPE_SON .eq. 1 ) THEN
         NBROW = NFRONT - NPIV
       ELSE
         NBROW = NELIM
       END IF
       IF (KEEP(50).EQ.0) THEN
        IF ((NPIV .NE. 0) .AND. (NBROW .GT. 1)) THEN
           IF ( IOLDPS .NE. PTLUST_S(STEP( INODE ) ) ) THEN
             WRITE(*,*) ' ERROR 5 in process_root2son : ',
     *       ' IOLDPS,PTLUST_S(STEP(INODE))=',
     *       IOLDPS,PTLUST_S(STEP(INODE))
             CALL SMUMPS_ABORT()
           END IF
         OPSFAC = POSELT + NPIV * NFRONT + NPIV
         APOS   = POSELT + (NPIV + 1) * NFRONT
         DO I = 1, NBROW - 1
           DO J = 1, NPIV
             A(OPSFAC) = A(APOS)
             OPSFAC = OPSFAC + 1
             APOS = APOS + 1
           END DO
           APOS = APOS + NBCOL
         END DO
        ENDIF
       ELSE
         IF ( TYPE_SON .eq. 1 ) THEN
           CALL SMUMPS_324(A(POSELT), NFRONT,
     *          NPIV, NBROW,KEEP(59))
         ELSE
           CALL SMUMPS_325(A(POSELT), NPIV, NBROW)
         END IF
       ENDIF
       IW(IOLDPS + XSIZE)     = NBCOL
       IW(IOLDPS + 1 +XSIZE) = NASS - NPIV
       IF (TYPE_SON.EQ.2) THEN
        IW(IOLDPS + 2 +XSIZE) = NASS
       ELSE
        IW(IOLDPS + 2 +XSIZE) = NFRONT
       ENDIF
       IW(IOLDPS + 3 +XSIZE) = NPIV
      CALL SMUMPS_93(MYID,N,IOLDPS,TYPE_SON,IW,LIW,
     *    A, LA, POSFAC, LRLU, LRLUS,
     *    IWPOS, PTRAST,PTRFAC,STEP, KEEP,KEEP8, .FALSE.,INODE,IERR)
      IF(IERR.LT.0)THEN
         IFLAG=IERR
         IERROR=0
         RETURN
      ENDIF
      ELSE 
        ISON = INODE
        PDEST_MASTER_ISON = SMUMPS_275(STEP(ISON),
     *      PROCNODE_STEPS,SLAVEF)
        DO WHILE (
     *     ( IW( PTRIST(STEP(ISON)) + 1  +XSIZE) .NE.
     *       IW( PTRIST(STEP(ISON)) + 3  +XSIZE) ) .OR.
     *     ( KEEP(50) .NE. 0 .AND.
     *       IW( PTRIST(STEP(ISON)) + 6  +XSIZE) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1  +XSIZE) .NE.
     *           IW( PTRIST(STEP(ISON)) + 3  +XSIZE) ) THEN
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING  = .TRUE.
          SET_IRECV = .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL SMUMPS_329( COMM_LOAD, ASS_IRECV,
     *    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *    MSGSOU, MSGTAG,
     *    STATUS,
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
     *    root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE, LPTRAR,
     *    NELT, FRTPTR, FRTELT,
     *    ISTEP_TO_INIV2, TAB_POS_IN_PERE )
          IF ( IFLAG .LT. 0 ) RETURN
        END DO
       IOLDPS = PTRIST(STEP(INODE))
       LCONT  = IW(IOLDPS+XSIZE)
       NROW   = IW(IOLDPS+2+XSIZE)
       NPIV   = IW(IOLDPS+3+XSIZE)
       NASS   = IW(IOLDPS+4+XSIZE)
       NELIM  = NASS-NPIV
           IF (NELIM.LE.0) THEN
            write(6,*) MYID,': ERROR 2 in SMUMPS_271 '
            write(6,*) MYID,': INODE,LCONT, NROW, NPIV, NASS, NELIM=',
     *      INODE,LCONT, NROW, NPIV, NASS, NELIM
            write(6,*) MYID,': IOLDPS=',IOLDPS
            CALL SMUMPS_ABORT()
           ENDIF
       NSLAVES= IW(IOLDPS+5+XSIZE)
       H_INODE = 6 + NSLAVES + XSIZE
       LIST_NELIM_COL = IOLDPS + H_INODE + NROW + NPIV
       NELIM_LOCAL = NELIM_ROOT
       DO I = 1, NELIM
        root%RG2L_COL(IW(LIST_NELIM_COL)) = NELIM_LOCAL
        NELIM_LOCAL = NELIM_LOCAL + 1
        LIST_NELIM_COL = LIST_NELIM_COL + 1
       ENDDO
       SHIFT_LIST_ROW_SON = 6 + IW(IOLDPS+5+XSIZE) + XSIZE
       SHIFT_LIST_COL_SON = SHIFT_LIST_ROW_SON + NROW + NPIV
       NCOL_TO_SEND       = NELIM
       IF (IW(IOLDPS+XXS).EQ.S_NOLCBNOCONTIG38.OR.
     *     IW(IOLDPS+XXS).EQ.S_ALL) THEN
         SHIFT_VAL_SON      = NPIV
         LDA                = LCONT + NPIV
       ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCBCONTIG38) THEN
         SHIFT_VAL_SON=NROW*(LCONT+NPIV-NELIM)
         LDA         =NELIM
       ELSE IF (IW(IOLDPS+XXS).EQ.S_NOLCLEANED38) THEN
         SHIFT_VAL_SON=0
         LDA = NELIM
       ELSE
         write(*,*) MYID,": internal error in SMUMPS_271",
     *   IW(IOLDPS+XXS), "INODE=",INODE
         CALL SMUMPS_ABORT()
       ENDIF
       IF ( KEEP( 50 ) .eq. 0 ) THEN
         INVERT = .FALSE.
       ELSE
         INVERT = .TRUE.
       END IF
       CALL SMUMPS_80( COMM_LOAD, ASS_IRECV, 
     *    N, INODE, FPERE,
     *    PTRIST, PTRAST,
     *    ROOT, NROW, NCOL_TO_SEND, SHIFT_LIST_ROW_SON,
     *    SHIFT_LIST_COL_SON , SHIFT_VAL_SON, LDA,
     *    ROOT_NON_ELIM_CB, MYID, COMM,
     *
     *    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     *    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     *    NIRBDU, PTRIST, PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER,
     *    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     *    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     *    OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *    INTARR, DBLARR, ICNTL, KEEP,KEEP8, INVERT, ND, FRERE,
     *    LPTRAR, NELT, FRTPTR, FRTELT, 
     *   ISTEP_TO_INIV2, TAB_POS_IN_PERE  )
        IF (IFLAG.LT.0 ) RETURN
       IF (KEEP(214).EQ.2) THEN
        CALL SMUMPS_314( N, INODE,
     *      PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA,
     *      LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     *      NIRBDU, IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER, ITLOC,
     *      IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, TYPE_SON
     $      )
       ENDIF
        IF (IFLAG.LT.0) THEN
           CALL SMUMPS_44( MYID, SLAVEF, COMM )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_271
       SUBROUTINE SMUMPS_217(N, NZ, NSCA, 
     *      ASPK, IRN, ICN, COLSCA, ROWSCA, S, MAXS, 
     *      ICNTL, INFO)
      INTEGER N, NZ, NSCA, MAXS
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(40), INFO(40)
      REAL    ASPK(NZ), COLSCA(*), ROWSCA(*)
      REAL    S(MAXS)
      INTEGER MPG,LP
      INTEGER ISPW1, IWNOR
      INTEGER I, K, ITOT
      LOGICAL PROK
      REAL ONE
      PARAMETER( ONE = 1.0E0 )
      LP      = ICNTL(1)
      MPG     = ICNTL(2)
      MPG    = ICNTL(3)
      PROK   = (MPG.GT.0)
      IF (PROK) WRITE(MPG,101)
 101    FORMAT(/' ****** SCALING OF ORIGINAL MATRIX '/)
        IF (NSCA.EQ.1) THEN
         IF (PROK)
     *    WRITE (MPG,*) ' DIAGONAL SCALING '
        ELSEIF (NSCA.EQ.2) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' SCALING BASED ON (MC29)'
        ELSEIF (NSCA.EQ.3) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' COLUMN SCALING'
        ELSEIF (NSCA.EQ.4) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' ROW AND COLUMN SCALING (1 Pass)'
        ELSEIF (NSCA.EQ.5) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' MC29 FOLLOWED BY ROW &COL SCALING'
        ELSEIF (NSCA.EQ.6) THEN
         IF (PROK)
     *   WRITE (MPG,*) ' MC29 FOLLOWED BY COLUMN SCALING'
        ENDIF
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
        IF ((NSCA.EQ.5).OR.
     &      (NSCA.EQ.6))                   THEN
          ITOT = 5*N + NZ 
          IF (ITOT.GT.MAXS) GOTO 400
          ISPW1 = MAXS - NZ + 1
          DO 15 K=1,NZ
           S(ISPW1+K-1) = ASPK(K)
  15      CONTINUE
        ELSE
          ISPW1 = MAXS + 1
          ITOT  = 5*N
          IF (ITOT.GT.MAXS) GOTO 400
        ENDIF
        IWNOR = ISPW1 - 5*N
          IF (NSCA.EQ.1) THEN
            CALL SMUMPS_238(N,NZ,ASPK,IRN,ICN,S(IWNOR),
     *        COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.2) THEN
            CALL SMUMPS_239(N,NZ,ASPK,IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
          ELSEIF (NSCA.EQ.3) THEN
            CALL SMUMPS_241(N,NZ,ASPK,IRN,ICN,S(IWNOR),COLSCA,
     *      MPG)
          ELSEIF (NSCA.EQ.4) THEN
            CALL SMUMPS_287(N,NZ,IRN,ICN,ASPK,
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.5) THEN
            CALL SMUMPS_239(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
            CALL SMUMPS_241(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),
     *          COLSCA, MPG)
          ELSEIF (NSCA.EQ.6) THEN
            CALL SMUMPS_239(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPG,MPG,NSCA)
            CALL SMUMPS_240(NSCA,N,NZ,IRN,ICN,S(ISPW1),
     *            S(IWNOR+N),ROWSCA,MPG)
            CALL SMUMPS_241(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),
     *          COLSCA, MPG)
          ENDIF
      GOTO 500
 400  INFO(1) = -5
      INFO(2) = ITOT-MAXS
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     * WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
 500  RETURN
      END SUBROUTINE SMUMPS_217
      SUBROUTINE SMUMPS_287(N,NZ,IRN,ICN,VAL,
     *    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N)
      REAL    COLSCA(N),ROWSCA(N)
      REAL    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      REAL    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO, ONE
      PARAMETER(ZERO=0.0E0, ONE=1.0E0)
      DO 50 J=1,N
       CNOR(J)   = REAL(ZERO)
       RNOR(J)   = REAL(ZERO)
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     *        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ABS(CNOR(J))) THEN
              CNOR(J) =     VDIAG
            ENDIF
            IF (VDIAG.GT.ABS(RNOR(I))) THEN
              RNOR(I) =     VDIAG
            ENDIF
 100   CONTINUE
      IF (MPRINT.GT.0) THEN
       CMIN = ABS( CNOR(1) )
       CMAX = ABS( CNOR(1) )
       RMIN = ABS( RNOR(1) )
       DO 111 I=1,N
        ARNOR = ABS(RNOR(I))
        ACNOR = ABS(CNOR(I))
        IF (ACNOR.GT.CMAX) CMAX=ACNOR
        IF (ACNOR.LT.CMIN) CMIN=ACNOR
        IF (ARNOR.LT.RMIN) RMIN=ARNOR
 111   CONTINUE
       WRITE(MPRINT,*) '**** STAT. OF MATRIX PRIOR ROW&COL SCALING'
       WRITE(MPRINT,*) ' MAXIMUM NORM-MAX OF COLUMNS:',CMAX
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF COLUMNS:',CMIN
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF ROWS   :',RMIN
      ENDIF
      DO 120 J=1,N
       IF (ABS(CNOR(J)).LE.ZERO) THEN
         CNOR(J)   = REAL(ONE)
       ELSE
         CNOR(J)   = REAL(ONE)/CNOR(J)
       ENDIF
 120  CONTINUE
      DO 130 J=1,N
       IF (ABS(RNOR(J)).LE.ZERO) THEN
         RNOR(J)   = REAL(ONE)
       ELSE
         RNOR(J)   = REAL(ONE)/RNOR(J)
       ENDIF
 130  CONTINUE
       DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
        COLSCA(I) = COLSCA(I) * CNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0)
     *  WRITE(MPRINT,*) ' END OF SCALING BY MAX IN ROW AND COL'
      RETURN
      END SUBROUTINE SMUMPS_287
      SUBROUTINE SMUMPS_239(N,NZ,VAL,ROWIND,COLIND,
     *       RNOR,CNOR,WNOR,MPRINT,MP,
     *       NSCA)
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N),
     *        WNOR(5*N)
      INTEGER COLIND(NZ),ROWIND(NZ)
      INTEGER J,I,K
      INTEGER MPRINT,MP,NSCA
      INTEGER IFAIL9
      REAL ZERO, ONE
      PARAMETER( ZERO = 0.0E0, ONE = 1.0E0 )
      DO 15 I=1,N
       RNOR(I)   = REAL(ZERO)
       CNOR(I)   = REAL(ZERO)
  15  CONTINUE
      CALL SMUMPS_216(N,N,NZ,VAL,ROWIND,COLIND,
     *   RNOR,CNOR,WNOR, MP,IFAIL9)
      DO 30 I=1,N
       CNOR(I) = EXP(CNOR(I))
       RNOR(I) = EXP(RNOR(I))
  30  CONTINUE
      IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
        DO 100 K=1,NZ
          I   = ROWIND(K)
          J   = COLIND(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 100
          VAL(K) = VAL(K) * CNOR(J) * RNOR(I)
 100    CONTINUE
      ENDIF
      IF (MPRINT.GT.0) 
     *   WRITE(MPRINT,*) ' END OF SCALING USING MC29'
      RETURN
      END SUBROUTINE SMUMPS_239
      SUBROUTINE SMUMPS_241(N,NZ,VAL,IRN,ICN,
     *       CNOR,COLSCA,MPRINT)
      INTEGER N,NZ
      REAL VAL(NZ),CNOR(N),COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      REAL VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO, ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      DO 10 J=1,N
       CNOR(J)   = REAL(ZERO)
  10  CONTINUE
      DO 100 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     *      (J.LE.0).OR.(J.GT.N)) GOTO 100
        VDIAG = ABS(VAL(K))
        IF (VDIAG.GT.ABS(CNOR(J))) THEN
           CNOR(J) =     VDIAG
        ENDIF
 100  CONTINUE
      DO 110 J=1,N
       IF (ABS(CNOR(J)).LE.ZERO) THEN
         CNOR(J)   = REAL(ONE)
       ELSE
         CNOR(J)   = REAL(ONE)/CNOR(J)
       ENDIF
 110  CONTINUE
       DO 215 I=1,N
        COLSCA(I) = COLSCA(I) * CNOR(I)
 215   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF COLUMN SCALING'
      RETURN
      END SUBROUTINE SMUMPS_241
      SUBROUTINE SMUMPS_238(N,NZ,VAL,IRN,ICN,RNOR,
     *      COLSCA,ROWSCA,MPRINT)
      INTEGER   N, NZ
      REAL RNOR(N)
      REAL  VAL(NZ), ROWSCA(N),COLSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      REAL      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC SQRT
      REAL ZERO, ONE
      PARAMETER(ZERO=0.0E0, ONE=1.0E0)
      DO 10 I=1,N
       RNOR(I)   = REAL(ONE)
  10  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          IF ((I.GT.N).OR.(I.LE.0)) GOTO 100
          J = ICN(K)
          IF (I.EQ.J) THEN
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ZERO) THEN
              RNOR(J) = REAL(ONE/(SQRT(VDIAG)))
            ENDIF
          ENDIF
 100   CONTINUE
       DO 110 I=1,N
        COLSCA(I) = RNOR(I)
        ROWSCA(I) = RNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF DIAGONAL SCALING'
      RETURN
      END SUBROUTINE SMUMPS_238
      SUBROUTINE SMUMPS_240(NSCA,N,NZ,IRN,ICN,VAL,
     *    RNOR,ROWSCA,MPRINT)
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      REAL VAL(NZ),RNOR(N), ROWSCA(N)
      REAL VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO,ONE
      PARAMETER (ZERO=0.0E0, ONE=1.0E0)
      DO 50 J=1,N
       RNOR(J)   = REAL(ZERO)
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     *        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ABS(RNOR(I))) THEN
              RNOR(I) =  VDIAG
            ENDIF
 100   CONTINUE
      DO 130 J=1,N
       IF (REAL(RNOR(J)).LE.REAL(ZERO)) THEN
         RNOR(J)   = REAL(ONE)
       ELSE
         RNOR(J)   = REAL(ONE)/RNOR(J)
       ENDIF
 130  CONTINUE
      DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
 110  CONTINUE
      IF ( (NSCA.EQ.4) .OR. (NSCA.EQ.6) ) THEN
        DO 150 K=1,NZ
          I   = IRN(K)
          J   = ICN(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 150
          VAL(K) = VAL(K) * RNOR(I)
 150    CONTINUE
      ENDIF
      IF (MPRINT.GT.0)
     *  WRITE(MPRINT,'(A)') '  END OF ROW SCALING'
      RETURN
      END
      SUBROUTINE SMUMPS_216(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      REAL A(NE)
      INTEGER IRN(NE),ICN(NE)
      REAL R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MIN
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      REAL ONE
      REAL SMIN,ZERO
      PARAMETER (ONE=1.0E0,SMIN=0.1E0,ZERO=0.0E0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      REAL E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3
      DO 10 I = 1,M
         R(I) = REAL(ZERO)
         W(I1+I) = REAL(ZERO)
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = REAL(ZERO)
         W(I2+J) = REAL(ZERO)
         W(I3+J) = REAL(ZERO)
         W(I4+J) = REAL(ZERO)
   20 CONTINUE
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE
      DO 60 K = 1,NE
         IF (ABS(A(K)).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = REAL(ZERO)
      Q = ONE
      S = REAL(ZERO)
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (ABS(S).LE.ABS(SM)) GO TO 160
      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (ABS(A(K)).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = REAL(ZERO)
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (ABS(S).LE.ABS(SM)) E = REAL(ZERO)
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (ABS(S).LE.ABS(SM)) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (ABS(A(K)).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = REAL(ZERO)
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (ABS(S).LE.ABS(SM)) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (ABS(S).LE.ABS(SM)) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (ABS(A(K)).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from SMUMPS_216 ****',' IFAIL =',IFAIL
      END SUBROUTINE SMUMPS_216
      SUBROUTINE SMUMPS_89(N,KEEP28, IW,LIW,A,LA,NIRBDU,
     *       LRLU,IPTRLU,IWPOS,
     *       IWPOSCB,PTRIST,PTRAST,STEP, PIMASTER,PAMASTER)
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU,KEEP28,
     &        IPTRLU,IWPOS,IWPOSCB
      INTEGER IW(LIW),PTRIST(KEEP28),PTRAST(KEEP28),
     &        STEP(N), 
     * PIMASTER(KEEP28),
     * PAMASTER(KEEP28)
      REAL A(LA)
      LOGICAL FREE, SAME_PROC
      INTEGER IPTIW,IPTA, SIZFI,SIZFR,LONGI,LONGR
      INTEGER NCOL,NROW,NSLSON,HF  ! Not needed if XSIZE>=3
      INTEGER IPTSHIFT
      INCLUDE 'mumps_headers.h'
      INTEGER I
      IF (IWPOSCB.EQ.NIRBDU) GOTO 100
      IPTIW = IWPOSCB
      IPTA  = IPTRLU
      LONGI = 0
      LONGR = 0
 10   FREE = IW(IPTIW+1+XXS).EQ.S_FREE
      SIZFR   = IW(IPTIW+1+XXR)
      SIZFI   = IW(IPTIW+1+XXI)
      IF (FREE) THEN
        IF (LONGI.NE.0) THEN
         DO I=0,LONGI-1
           IW(IPTIW + SIZFI - I) = IW (IPTIW - I )
         ENDDO
         DO I=0,LONGR-1
           A(IPTA + SIZFR - I)   = A(IPTA - I )
         ENDDO
          DO I=1,KEEP28
           IF  ((PTRIST(I).LE.IPTIW).AND.
     &       (PTRIST(I).GT.IWPOSCB) ) THEN
            PTRIST(I) = PTRIST(I) + SIZFI
            IF (IW(PTRIST(I)+XSIZE).LT.0) THEN
             PAMASTER(I) = PAMASTER(I) + SIZFR
            ELSE
             PTRAST(I) = PTRAST(I) + SIZFR
            ENDIF
           ENDIF 
           IF ((PIMASTER(I).LE.IPTIW).AND.
     &       (PIMASTER(I).GT.IWPOSCB) ) THEN
            PIMASTER(I) = PIMASTER(I) + SIZFI
            PAMASTER(I) = PAMASTER(I) + SIZFR
           ENDIF
          ENDDO 
        ENDIF
        IWPOSCB = IWPOSCB + SIZFI
        IPTRLU = IPTRLU + SIZFR
        LRLU   = LRLU + SIZFR
      ELSE
       LONGI = LONGI + SIZFI
       LONGR = LONGR + SIZFR
      ENDIF
      IPTA  = IPTA  + SIZFR
      IPTIW = IPTIW + SIZFI
      IF (IPTIW.NE.NIRBDU) GOTO 10
 100  CONTINUE
       RETURN
       END SUBROUTINE SMUMPS_89
      SUBROUTINE SMUMPS_221(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,
     *     DKEEP,PIVNUL_LIST,LPN_LIST)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      REAL A(LA) 
      REAL UU, SEUIL
      INTEGER IW(LIW) 
      INTEGER  IOLDPS, POSELT
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      INCLUDE 'mumps_headers.h'
      REAL SWOP
      INTEGER APOS
      REAL AMROW
      REAL ZERO,RMAX,ONE
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INTRINSIC MAX
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        NASSW   = IABS(IW(IOLDPS+3+XSIZE))
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT*(NPIVP1-1) + NPIV
           IDIAG = APOS
           IF(ABS(A(APOS)).LT.SEUIL) THEN
              IF(REAL(A(APOS)) .GE. ZERO) THEN
                 A(APOS) = SEUIL
              ELSE
                 A(APOS) = -SEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (ABS(A(APOS)).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = ZERO
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
             J3    = NASS -NPIV
             JMAX  = ISAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + IPIV - NPIVP1
            IF (RMAX.LE.DKEEP(1)) THEN
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE+
     &                      IW(IOLDPS+5+XSIZE)+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(DKEEP(2).GT.ZERO) THEN
                  IF(REAL(A(IDIAG)) .GE. ZERO) THEN
                     A(IDIAG) = DKEEP(2)
                  ELSE
                     A(IDIAG) = -DKEEP(2)
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS - NPIV + NFRONT - 1
                 DO JJ=J1,J2
                   A(JJ)= REAL(ZERO)
                 ENDDO
                 A(IDIAG) = REAL(ONE)
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 380   
            ENDIF
            IF (ABS(A(IDIAG)).GT. MAX(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LE. MAX(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV*NFRONT
            J2 = J1 + NFRONT - 1
            J3 = POSELT + (IPIV-1)*NFRONT
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV
            J2 = POSELT + NPIV + JMAX - 1
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT
              J2 = J2 + NFRONT
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NFRONT + NPIV + 1 +XSIZE
            ISWPS2 = IOLDPS + 5 + NFRONT + NPIV + JMAX +XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'Detected a null pivot, INODE/NPIV=',INODE,NPIV
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_221
      SUBROUTINE SMUMPS_220(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *   INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,LA,INODE,INOPV
      REAL UU, SEUIL
      REAL A(LA)
      INTEGER IW(LIW)
      REAL AMROW
      REAL ZERO,RMAX
      REAL  SWOP
      INTEGER APOS, POSELT, IOLDPS
      INTEGER NOFFW,NPIV,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INCLUDE 'mumps_headers.h'
      INTRINSIC MAX
      DATA ZERO /0.0E0/
        INOPV   = 0
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
          DO 460 IPIV=NPIVP1,NASS
            APOS = POSELT + NFRONT*NPIV + (IPIV-1)
            JMAX = 1
            AMROW = ZERO
            J1 = APOS
            J3    = NASS -NPIV
            JMAX  = ISAMAX(J3,A(J1),NFRONT)
            JJ    = J1 + (JMAX-1)*NFRONT
            AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = APOS +  (NASS-NPIV) * NFRONT
            J3 = NFRONT - NASS
            IF (J3.EQ.0) GOTO 370
            DO 360 JJ=1,J3
              RMAX = MAX(ABS(A(J1)),RMAX)
              J1 = J1 + NFRONT
  360       CONTINUE
  370       IF (RMAX.EQ.ZERO) GO TO 460
            IDIAG = APOS + (IPIV - NPIVP1)*NFRONT
            IF (ABS(A(IDIAG)).GE.MAX(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LT.MAX(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV
            J3 = POSELT + (IPIV-1)
            DO 390 JJ= 1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J3)
              A(J3) = SWOP
              J1 = J1 + NFRONT
              J3 = J3 + NFRONT
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV*NFRONT
            J2 = POSELT + (NPIV + JMAX - 1)*NFRONT
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1
              J2 = J2 + 1
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIV + 1 + XSIZE
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
       INOPV = 1
  420 RETURN
      END SUBROUTINE SMUMPS_220
      SUBROUTINE SMUMPS_225(IBEG_BLOCK,
     *     NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *     IOLDPS,POSELT,IFINB,LKJIB,LKJIT)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB,LKJIB,IBEG_BLOCK
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    VALPIV
      INTEGER APOS, UUPOS, IOLDPS, POSELT
      INTEGER LKJIT
      REAL ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,LPOS,NEL
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3+XSIZE) = NASS
          ELSE
           IW(IOLDPS+3+XSIZE) = MIN0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3+XSIZE)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3+XSIZE) = MIN0(JROW2+LKJIB,NASS)
          IBEG_BLOCK = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + NPIV*(NFRONT + 1)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541     CONTINUE
         LPOS   = APOS + NFRONT
         UUPOS  = APOS+1
         CALL SGER(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     *              A(LPOS+1),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE SMUMPS_225
      SUBROUTINE SMUMPS_229(NFRONT,N,INODE,IW,LIW,A,LA,IOLDPS,
     *          POSELT)
      IMPLICIT NONE
      INTEGER NFRONT,N,INODE,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT, UUPOS
      INTEGER IOLDPS,NPIV,NEL
      INTEGER LPOS,JROW,IRWPOS
      INCLUDE 'mumps_headers.h'
      REAL ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 650
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 340 JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
  340   CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 JROW = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
  650   RETURN
        END SUBROUTINE SMUMPS_229
      SUBROUTINE SMUMPS_228(NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *       IOLDPS,POSELT,IFINB)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFINB
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT,UUPOS
      INTEGER IOLDPS,NPIV,KROW
      INTEGER NEL,LPOS,ICOL,NEL2,IRWPOS
      INTEGER NPIVP1
      REAL ONE
      DATA ONE /1.0E0/
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + NPIV*(NFRONT + 1)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 541 KROW = 1,NEL
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541    CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 ICOL = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
        RETURN
        END SUBROUTINE SMUMPS_228
      SUBROUTINE SMUMPS_231(A,LA,NFRONT,
     *       NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER LA,POSELT
      REAL    A(LA)
      INTEGER NFRONT, NPIV, NASS
      INTEGER NEL1,NEL11,LPOS2,LPOS1,LPOS
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + NASS*NFRONT
        CALL STRSM('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     *              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + NPIV
        LPOS1  = POSELT + NPIV
        CALL SGEMM('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_231
      SUBROUTINE SMUMPS_232(A,LA,NFRONT,NPIV,NASS,POSELT,LKJIB)
      INTEGER LA, NFRONT, NPIV, NASS, LKJIB
      REAL    A(LA)
      INTEGER POSELT
      INTEGER POSELT_LOCAL
      INTEGER NEL1, NEL11, NPBEG, LPOS, LPOS1, LPOS2
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        POSELT_LOCAL = POSELT
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT_LOCAL + NPIV*NFRONT + NPBEG - 1
        POSELT_LOCAL = POSELT_LOCAL + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL STRSM('L','L','N','N',LKJIB,NEL1,ONE,A(POSELT_LOCAL),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT_LOCAL + LKJIB
        CALL SGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *       NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_232
      SUBROUTINE SMUMPS_233(IBEG_BLOCK,
     *    NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     *    IOLDPS,POSELT,LKJIB,LKJIT )
      IMPLICIT NONE
      INTEGER NFRONT, NASS,N,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW) 
      INTEGER LKJIB, INODE, IBEG_BLOCK
      INTEGER POSELT
      INTEGER IOLDPS, NPIV, JROW2, NPBEG
      INTEGER NONEL, LKJIW, NEL1, NEL11
      INTEGER LBP, IPOS, KPOS, LPOS2, HF
      INTEGER LPOS1,LPOS,LBPT,I1,K1,II,ISWOP,LBP1
      INTEGER LKJIT, POSLOCAL
      INCLUDE 'mumps_headers.h'
      REAL ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NPIV   = IW(IOLDPS+1+XSIZE)
        JROW2  = IABS(IW(IOLDPS+3+XSIZE))
        NPBEG  = IBEG_BLOCK
        HF     = 6 + IW(IOLDPS+5+XSIZE) +XSIZE
        NONEL         = JROW2 - NPIV + 1
        IF ((NASS-NPIV).GE.LKJIT) THEN
         LKJIB       = LKJIB + NONEL
         IW(IOLDPS+3+XSIZE)= MIN0(NPIV+LKJIB,NASS)
        ELSE
          IW(IOLDPS+3+XSIZE) = NASS
        ENDIF
        IBEG_BLOCK = NPIV + 1
        NEL1   = NASS - JROW2
        LKJIW  = NPIV - NPBEG + 1
        NEL11  = NFRONT - NPIV
        IF ((NEL1.EQ.0).OR.(LKJIW.EQ.0)) GO TO 500
        LPOS2  = POSELT + JROW2*NFRONT + NPBEG - 1
         POSLOCAL = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
         CALL STRSM('L','L','N','N',LKJIW,NEL1,ONE,
     *               A(POSLOCAL),NFRONT,
     *               A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIW
        LPOS1  = POSLOCAL + LKJIW
        CALL SGEMM('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   LBP = JROW2 - NPIV
         LBPT  = LBP + LBP
        IF  ((NEL1.GE.LBPT).AND.(NEL1.GE.LKJIB)) THEN
         I1 = IOLDPS + HF + NPIV
         K1 = IOLDPS + HF + NASS - LBP
         DO 10 II=1,LBP
          ISWOP  = IW(I1)
          IW(I1) = IW(K1)
          IW(K1) = ISWOP
          I1     = I1 +1
          K1     = K1 + 1
  10     CONTINUE
         IPOS = POSELT + NPIV*NFRONT
         KPOS = POSELT + (NASS-LBP)*NFRONT
         LBP1 = LBP * NFRONT
         CALL SSWAP(LBP1,A(IPOS),1,A(KPOS),1)
        ENDIF
        END SUBROUTINE SMUMPS_233
      SUBROUTINE SMUMPS_236(A,LA,NPIVB,NFRONT,
     *                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS,LA
      REAL    A(LA)
      INTEGER APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER LPOS, LPOS1, LPOS2, NEL1, NEL11, NPIVE
      REAL    ALPHA, ONE
      PARAMETER(ONE=1.0E0, ALPHA=-1.0E0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + NPIVB*NFRONT + NPIVB
        LPOS2  = APOS + NASSL
        CALL STRSM('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     *              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + NFRONT*NPIVE
        LPOS1  = APOS + NFRONT*NPIVE
        CALL SGEMM('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     *          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        END SUBROUTINE SMUMPS_236
      SUBROUTINE SMUMPS_27( id,  ANORMINF, LSCAL )
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE(SMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      INTEGER, DIMENSION (:), POINTER :: KEEP,INFO
      INTEGER*8, DIMENSION (:), POINTER :: KEEP8
      LOGICAL :: I_AM_SLAVE
      REAL DUMMY(1)
      REAL ZERO
      PARAMETER( ZERO = 0.0E0)
      REAL, ALLOCATABLE :: SUMR(:), SUMR_LOC(:)
      INTEGER :: ALLOCOK, MTYPE, I
      INFO =>id%INFO
      KEEP =>id%KEEP
      KEEP8 =>id%KEEP8
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     *             ( id%MYID .eq. MASTER .AND.
     *               KEEP(46) .eq. 1 ) )
      IF (id%MYID .EQ. MASTER) THEN
       ALLOCATE( SUMR( id%N ), stat =allocok )
       IF (allocok .GT.0 ) THEN
        id%INFO(1)=-13
        id%INFO(2)=id%N
        RETURN
       ENDIF
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
             IF (.NOT.LSCAL) THEN
              CALL SMUMPS_207(id%A(1),
     *          id%NZ, id%N,
     *          id%IRN(1), id%JCN(1),
     *          SUMR, KEEP,KEEP8 )
             ELSE
              CALL SMUMPS_289(id%A(1),
     *          id%NZ, id%N,
     *          id%IRN(1), id%JCN(1), 
     *          SUMR, KEEP, KEEP8,
     *          id%COLSCA(1))
             ENDIF
            ELSE
             MTYPE = 1
             IF (.NOT.LSCAL) THEN
              CALL SMUMPS_119(MTYPE, id%N,
     *           id%NELT, id%ELTPTR,
     *           id%LELTVAR, id%ELTVAR,
     *           id%NA_ELT, id%A_ELT(1),
     *           SUMR, KEEP,KEEP8 )
             ELSE
              CALL SMUMPS_135(MTYPE, id%N,
     *           id%NELT, id%ELTPTR(1),
     *           id%LELTVAR, id%ELTVAR(1),
     *           id%NA_ELT, id%A_ELT(1),
     *           SUMR, KEEP,KEEP8, id%COLSCA(1))
             ENDIF
            ENDIF
          ENDIF
      ELSE
          LSCAL = .FALSE.
          ALLOCATE( SUMR_LOC( id%N ), stat =allocok )
          IF (allocok .GT.0 ) THEN
             id%INFO(1)=-13
             id%INFO(2)=id%N
             RETURN
          ENDIF
          IF ( I_AM_SLAVE .and.
     *           id%NZ_loc .NE. 0 ) THEN
           IF (.NOT.LSCAL) THEN
              CALL SMUMPS_207(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc, 
     *          SUMR_LOC, id%KEEP,id%KEEP8 )
           ELSE
              CALL SMUMPS_289(id%A_loc,
     *          id%NZ_loc, id%N,
     *          id%IRN_loc, id%JCN_loc, 
     *          SUMR_LOC, id%KEEP,id%KEEP8,
     *          id%COLSCA)
           ENDIF
          ELSE
           SUMR_LOC = REAL(ZERO)
          ENDIF
          IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SUMR_LOC, SUMR,
     *        id%N, MPI_REAL,
     *        MPI_SUM,MASTER,id%COMM, IERR)
          ELSE
              CALL MPI_REDUCE( SUMR_LOC, DUMMY,
     *        id%N, MPI_REAL,
     *        MPI_SUM,MASTER,id%COMM, IERR)
          END IF
        DEALLOCATE (SUMR_LOC)
      ENDIF
      IF ( id%MYID .eq. MASTER ) THEN
       ANORMINF = REAL(ZERO)
        IF (LSCAL) THEN
         DO I = 1, id%N
          ANORMINF = MAX(ABS(id%ROWSCA(I) * SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ELSE
         DO I = 1, id%N
          ANORMINF = MAX(ABS(SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ENDIF
      ENDIF
      CALL MPI_BCAST(ANORMINF, 1,
     *              MPI_REAL, MASTER,
     *              id%COMM, IERR )
      IF (id%MYID .eq. MASTER) DEALLOCATE (SUMR)
      RETURN
      END SUBROUTINE SMUMPS_27
      SUBROUTINE SMUMPS_628(IW,LREC,SIZE_FREE)
      INTEGER, intent(in) :: LREC
      INTEGER, intent(in) :: IW(LREC)
      INTEGER, intent(out):: SIZE_FREE
      INCLUDE 'mumps_headers.h'
      IF (IW(1+XXS).EQ.S_NOLCBCONTIG .OR.
     *    IW(1+XXS).EQ.S_NOLCBNOCONTIG) THEN
        SIZE_FREE=IW(1+XSIZE+2)*IW(1+XSIZE+3)
      ELSE IF (IW(1+XXS).EQ.S_NOLCBCONTIG38 .OR.
     *         IW(1+XXS).EQ.S_NOLCBNOCONTIG38) THEN
        SIZE_FREE=IW(1+XSIZE+2)*(IW(1+XSIZE)+
     *            IW(1+XSIZE + 3) -
     *          ( IW(1+XSIZE + 4)
     *          - IW(1+XSIZE + 3) ) )
      ELSE
        SIZE_FREE=0
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_628
      SUBROUTINE SMUMPS_629
     *(IW,LIW,IXXP,ICURRENT,NEXT, RCURRENT,ISIZE2SHIFT)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER LIW,IXXP,ICURRENT,NEXT,RCURRENT,ISIZE2SHIFT
      INTEGER IW(LIW)
      ICURRENT=NEXT
      RCURRENT=RCURRENT-IW(ICURRENT+XXR)
      NEXT=IW(ICURRENT+XXP)
      IW(IXXP)=ICURRENT+ISIZE2SHIFT
      IXXP=ICURRENT+XXP
      RETURN
      END SUBROUTINE SMUMPS_629
      SUBROUTINE SMUMPS_630(IW,LIW,BEG2SHIFT,END2SHIFT,ISIZE2SHIFT)
      IMPLICIT NONE
      INTEGER LIW, BEG2SHIFT, END2SHIFT, ISIZE2SHIFT
      INTEGER IW(LIW)
      INTEGER I
      IF (ISIZE2SHIFT.GT.0) THEN
        DO I=END2SHIFT,BEG2SHIFT,-1
          IW(I+ISIZE2SHIFT)=IW(I)
        ENDDO
      ELSE IF (ISIZE2SHIFT.LT.0) THEN
        DO I=BEG2SHIFT,END2SHIFT
          IW(I+ISIZE2SHIFT)=IW(I)
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_630
      SUBROUTINE SMUMPS_631(A, LA, BEG2SHIFT, END2SHIFT, RSIZE2SHIFT)
      IMPLICIT NONE
      INTEGER LA, BEG2SHIFT, END2SHIFT, RSIZE2SHIFT
      REAL A(LA)
      INTEGER I
      IF (RSIZE2SHIFT.GT.0) THEN
        DO I=END2SHIFT,BEG2SHIFT,-1
          A(I+RSIZE2SHIFT)=A(I)
        ENDDO
      ELSE IF (RSIZE2SHIFT.LT.0) THEN
        DO I=BEG2SHIFT,END2SHIFT
          A(I+RSIZE2SHIFT)=A(I)
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_631
      SUBROUTINE SMUMPS_94(N,KEEP28,IW,LIW,A,LA,NIRBDU,
     *       LRLU,IPTRLU,IWPOS,
     *       IWPOSCB,PTRIST,PTRAST,STEP,PIMASTER,PAMASTER, ITLOC,
     *       KEEP216,LRLUS)
      IMPLICIT NONE
      INTEGER N,LIW,LA,NIRBDU,LRLU,KEEP28,
     &        IPTRLU,IWPOS,IWPOSCB,KEEP216
      INTEGER, intent(IN):: LRLUS
      INTEGER IW(LIW),PTRIST(KEEP28),PTRAST(KEEP28),
     &        STEP(N),
     * PIMASTER(KEEP28),
     * PAMASTER(KEEP28), ITLOC(N)
      REAL A(LA)
      INCLUDE 'mumps_headers.h' ! headers and status of records
      INTEGER ICURRENT, NEXT, RCURRENT, STATE_NEXT
      INTEGER ISIZE2SHIFT, RSIZE2SHIFT
      INTEGER IBEGCONTIG, RBEGCONTIG
      INTEGER INODE
      INTEGER FREE_IN_REC
      INTEGER IXXP
      ISIZE2SHIFT=0
      RSIZE2SHIFT=0
      ICURRENT  = NIRBDU-XSIZE+1
      RCURRENT = LA+1
      IBEGCONTIG = -999999 ! only used for compress
      RBEGCONTIG = -999999 ! only used for compress
      NEXT = IW(ICURRENT+XXP)
      STATE_NEXT = IW(NEXT+XXS)
      IXXP = ICURRENT+XXP
  10     CONTINUE
         IF ( STATE_NEXT .NE. S_FREE .AND.
     *        (KEEP216.EQ.3.OR.
     *         (STATE_NEXT .NE. S_NOLCBNOCONTIG .AND.
     *          STATE_NEXT .NE. S_NOLCBCONTIG .AND.
     *          STATE_NEXT .NE. S_NOLCBNOCONTIG38 .AND.
     *          STATE_NEXT .NE. S_NOLCBCONTIG38))) THEN
            CALL SMUMPS_629(IW,LIW,
     *           IXXP, ICURRENT, NEXT, RCURRENT, ISIZE2SHIFT)
            IF (IBEGCONTIG < 0) THEN
              IBEGCONTIG=ICURRENT+IW(ICURRENT+XXI)-1
            ENDIF
            IF (RBEGCONTIG < 0) THEN
              RBEGCONTIG=RCURRENT+IW(ICURRENT+XXR)-1
            ENDIF
            INODE=IW(ICURRENT+XXN)
            IF (RSIZE2SHIFT .NE. 0) THEN
                IF (PTRAST(STEP(INODE)).EQ.RCURRENT)
     *            PTRAST(STEP(INODE))=
     *            PTRAST(STEP(INODE))+RSIZE2SHIFT
                IF (PAMASTER(STEP(INODE)).EQ.RCURRENT)
     *            PAMASTER(STEP(INODE))=
     *            PAMASTER(STEP(INODE))+RSIZE2SHIFT
            ENDIF
            IF (ISIZE2SHIFT .NE. 0) THEN
                IF (PTRIST(STEP(INODE)).EQ.ICURRENT)
     *            PTRIST(STEP(INODE))=
     *            PTRIST(STEP(INODE))+ISIZE2SHIFT
                IF (PIMASTER(STEP(INODE)).EQ.ICURRENT)
     *            PIMASTER(STEP(INODE))=
     *            PIMASTER(STEP(INODE))+ISIZE2SHIFT
            ENDIF
            IF (NEXT .NE. TOP_OF_STACK) THEN
              STATE_NEXT=IW(NEXT+XXS)
              GOTO 10
            ENDIF
         ENDIF
  20     CONTINUE
         IF (IBEGCONTIG.NE.0 .AND. ISIZE2SHIFT .NE. 0) THEN
           CALL SMUMPS_630(IW,LIW,ICURRENT,IBEGCONTIG,ISIZE2SHIFT)
           IF (IXXP .LE.IBEGCONTIG) THEN
           IXXP=IXXP+ISIZE2SHIFT
           ENDIF
         ENDIF
         IBEGCONTIG=-9999
  25     CONTINUE
         IF (RBEGCONTIG .GT.0 .AND. RSIZE2SHIFT .NE. 0) THEN
           CALL SMUMPS_631(A,LA,RCURRENT,RBEGCONTIG,RSIZE2SHIFT)
         ENDIF
         RBEGCONTIG=-99999
  30     CONTINUE
         IF (NEXT.EQ. TOP_OF_STACK) GOTO 100
         IF (STATE_NEXT .EQ. S_NOLCBCONTIG .OR.
     *       STATE_NEXT .EQ. S_NOLCBNOCONTIG .OR.
     *       STATE_NEXT .EQ. S_NOLCBCONTIG38 .OR.
     *       STATE_NEXT .EQ. S_NOLCBNOCONTIG38) THEN
           IF ( KEEP216.eq.3) THEN
             WRITE(*,*) "Internal error 2 in SMUMPS_94"
           ENDIF
           IF (RBEGCONTIG > 0) GOTO 25
           CALL SMUMPS_629
     *       (IW,LIW,IXXP,ICURRENT,NEXT, RCURRENT,ISIZE2SHIFT)
           IF (IBEGCONTIG < 0 ) THEN
             IBEGCONTIG=ICURRENT+IW(ICURRENT+XXI)-1
           ENDIF
           CALL SMUMPS_628(IW(ICURRENT),
     *                              LIW-ICURRENT+1,
     *                              FREE_IN_REC)
           IF (STATE_NEXT .EQ. S_NOLCBNOCONTIG) THEN
             CALL SMUMPS_627(A,LA,RCURRENT,
     *            IW(ICURRENT+XSIZE+2),
     *            IW(ICURRENT+XSIZE),
     *            IW(ICURRENT+XSIZE)+IW(ICURRENT+XSIZE+3), 0,
     *            IW(ICURRENT+XXS),RSIZE2SHIFT) ! State=S_NOLCBCONTIG
           ELSE IF (STATE_NEXT .EQ. S_NOLCBNOCONTIG38) THEN
             CALL SMUMPS_627(A,LA,RCURRENT,
     *            IW(ICURRENT+XSIZE+2),
     *            IW(ICURRENT+XSIZE),
     *            IW(ICURRENT+XSIZE)+IW(ICURRENT+XSIZE+3),
     *            IW(ICURRENT+XSIZE+4)-IW(ICURRENT+XSIZE+3), ! NELIM
     *            IW(ICURRENT+XXS),RSIZE2SHIFT) ! State=S_NOLCBCONTIG
           ELSE IF (RSIZE2SHIFT .GT.0) THEN
           CALL SMUMPS_631(A, LA, RCURRENT+FREE_IN_REC,
     *                       RCURRENT+IW(ICURRENT+XXR)-1,
     *                       RSIZE2SHIFT)
           ENDIF
           INODE=IW(ICURRENT+XXN)
           IF (ISIZE2SHIFT.NE.0) THEN
             PTRIST(STEP(INODE))=PTRIST(STEP(INODE))+ISIZE2SHIFT
           ENDIF
           PTRAST(STEP(INODE))=PTRAST(STEP(INODE))+RSIZE2SHIFT+
     *                         FREE_IN_REC
           IW(ICURRENT+XXR)=IW(ICURRENT+XXR)-FREE_IN_REC
           IF (STATE_NEXT.EQ.S_NOLCBCONTIG.OR.
     *         STATE_NEXT.EQ.S_NOLCBNOCONTIG) THEN
             IW(ICURRENT+XXS)=S_NOLCLEANED
           ELSE
             IW(ICURRENT+XXS)=S_NOLCLEANED38
           ENDIF
           RSIZE2SHIFT=RSIZE2SHIFT+FREE_IN_REC
           RBEGCONTIG=-9999
           IF (NEXT.EQ.TOP_OF_STACK) THEN
             GOTO 20
           ELSE
             STATE_NEXT=IW(NEXT+XXS)
           ENDIF
           GOTO 30
         ENDIF
         IF (IBEGCONTIG.GT.0) THEN
           GOTO 20
         ENDIF
  40     CONTINUE
         IF (STATE_NEXT == S_FREE) THEN
            ICURRENT=NEXT
            RCURRENT=RCURRENT-IW(ICURRENT+XXR)
            ISIZE2SHIFT = ISIZE2SHIFT+IW(ICURRENT+XXI)
            RSIZE2SHIFT = RSIZE2SHIFT+IW(ICURRENT+XXR)
            NEXT=IW(ICURRENT+XXP)
            IF (NEXT.EQ.TOP_OF_STACK) THEN
              WRITE(*,*) "Internal error 1 in SMUMPS_94"
              CALL SMUMPS_ABORT()
            ENDIF
            STATE_NEXT  = IW(NEXT+XXS)
            GOTO 40
         ENDIF
      GOTO 10
 100  CONTINUE
      IWPOSCB = IWPOSCB + ISIZE2SHIFT
      LRLU    = LRLU    + RSIZE2SHIFT
      IPTRLU  = IPTRLU  + RSIZE2SHIFT
      RETURN
      END SUBROUTINE SMUMPS_94
      SUBROUTINE SMUMPS_632(IREC, IW, NIRBDU,
     *            ISIZEHOLE, RSIZEHOLE)
      IMPLICIT NONE
      INTEGER, intent(in) :: IREC, NIRBDU
      INTEGER, intent(in) :: IW(NIRBDU)
      INTEGER, intent(out):: ISIZEHOLE, RSIZEHOLE
      INTEGER IRECLOC
      INCLUDE 'mumps_headers.h'
      ISIZEHOLE=0
      RSIZEHOLE=0
      IRECLOC = IREC + IW( IREC+XXI )
 10   CONTINUE
      IF (IW(IRECLOC+XXS).EQ.S_FREE) THEN
        ISIZEHOLE=ISIZEHOLE+IW(IRECLOC+XXI)
        RSIZEHOLE=RSIZEHOLE+IW(IRECLOC+XXR)
        IRECLOC=IRECLOC+IW(IRECLOC+XXI)
        GOTO 10
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_632
      SUBROUTINE SMUMPS_627(A, LA, RCURRENT,
     *           NROW, NCB, LD, NELIM, NODESTATE, ISHIFT)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER LD, NROW, NCB, NELIM, NODESTATE
      INTEGER ISHIFT, LA, RCURRENT
      REAL A(LA)
      INTEGER I,J,IOLD,INEW
      LOGICAL NELIM_ROOT
      NELIM_ROOT=.TRUE.
      IF (NODESTATE.EQ. S_NOLCBNOCONTIG) THEN
         NELIM_ROOT=.FALSE.
         IF (NELIM.NE.0)  THEN
           WRITE(*,*) "Internal error 1 IN SMUMPS_627"
           CALL SMUMPS_ABORT()
         ENDIF
      ELSE IF (NODESTATE .NE. S_NOLCBNOCONTIG38) THEN
           WRITE(*,*) "Internal error 2 in SMUMPS_627"
     *                ,NODESTATE
           CALL SMUMPS_ABORT()
      ENDIF
      IF (ISHIFT .LT.0) THEN
        WRITE(*,*) "Internal error 3 in SMUMPS_627",ISHIFT
        CALL SMUMPS_ABORT()
      ENDIF
      IF (NELIM_ROOT) THEN
        IOLD=RCURRENT+LD*NROW-1-NCB+NELIM
      ELSE
        IOLD = RCURRENT+LD*NROW-1
      ENDIF
      INEW = RCURRENT+LD*NROW+ISHIFT-1
      DO I = NROW, 1, -1
        IF (I.EQ.NROW .AND. ISHIFT.EQ.0.AND.
     *    .NOT. NELIM_ROOT) THEN
          IOLD=IOLD-LD
          INEW=INEW-NCB
          CYCLE
        ENDIF
        IF (NELIM_ROOT) THEN
          DO J=1,NELIM
            A( INEW ) = A( IOLD - J + 1)
            INEW = INEW - 1
          ENDDO
        ELSE
          DO J=1, NCB
            A( INEW ) = A( IOLD - J + 1)
            INEW = INEW - 1
          ENDDO
        ENDIF
        IOLD = IOLD - LD
      ENDDO
      IF (NELIM_ROOT) THEN
        NODESTATE=S_NOLCBCONTIG38
      ELSE
        NODESTATE=S_NOLCBCONTIG
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_627
      SUBROUTINE SMUMPS_272(BUFR,LBUFR,
     *     LBUFR_BYTES,
     *     root, N, IW, LIW, A, LA,
     *     NIRBDU, NBPROCFILS, LRLU, IPTRLU, IWPOS, IWPOSCB,
     *     PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     *     COMP, LRLUS, IPOOL, LPOOL, LEAF,
     *     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *     KEEP,KEEP8, IFLAG, IERROR, COMM, COMM_LOAD, ITLOC,
     *     ND,PROCNODE_STEPS,SLAVEF )
      USE SMUMPS_LOAD
      USE SMUMPS_OOC        
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC ) :: ROOT
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER LBUFR, LBUFR_BYTES, N, LIW, LA, NIRBDU, LRLU, IPTRLU,
     *        IWPOS, IWPOSCB, COMP, COMM, COMM_LOAD, IFLAG,
     *        IERROR, LRLUS
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LEAF )
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), ITLOC( N )
      INTEGER BUFR( LBUFR_BYTES ), NBPROCFILS( KEEP(28) )
      INTEGER IW( LIW )
      INTEGER ND(KEEP(28)), PROCNODE_STEPS(KEEP(28)),SLAVEF
      REAL A( LA )
      INTEGER   MYID
      INTEGER FILS( N ), PTRAIW(N), PTRARW( N )
      INTEGER INTARR(MAX(1,KEEP(14)))
      REAL DBLARR(MAX(1,KEEP(13)))
        INCLUDE 'mpif.h'
        INTEGER IERR
        INTEGER POSITION, LOCAL_M, LOCAL_N, POS_ROOT, LREQI, LREQA
        INTEGER NROW_SON, NCOL_SON, IROOT, ISON
        INCLUDE 'mumps_headers.h'
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   ISON, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NROW_SON, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   NCOL_SON, 1, MPI_INTEGER, COMM, IERR )
        IROOT = KEEP( 38 )
        IF ( PTRIST( STEP(IROOT) ) .NE. 0 .OR.
     *       PTLUST_S( STEP(IROOT)) .NE. 0 ) THEN
            NBPROCFILS(STEP(IROOT)) = NBPROCFILS(STEP(IROOT))-1
            IF ( NBPROCFILS( STEP(IROOT) ) .eq. 0 ) THEN
              IF (KEEP(201).NE.0) THEN
              CALL SMUMPS_580(IERR)
              ENDIF
              CALL SMUMPS_507( N, IPOOL, LPOOL,
     *             PROCNODE_STEPS, SLAVEF, KEEP(28), KEEP(76),
     *             KEEP(80), KEEP(47),
     *             STEP, IROOT + N)
              IF (KEEP(47) .GE. 3) THEN
                 CALL SMUMPS_500(
     $                IPOOL, LPOOL, 
     *                PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *                MYID, STEP, N, ND, FILS )
              ENDIF
            END IF
        ELSE
            NBPROCFILS(STEP( IROOT ) ) = -1
           IF (KEEP(60) == 0) THEN
            CALL SMUMPS_284( root, IROOT, N,
     *                     IW, LIW, A, LA,
     *                     NIRBDU,
     *                     FILS, MYID, PTRAIW, PTRARW, INTARR, DBLARR,
     *                     LRLU, IPTRLU,
     *                     IWPOS, IWPOSCB, PTRIST, PTRAST,
     *                     STEP, PIMASTER, PAMASTER, ITLOC,
     *                     COMP, LRLUS, IFLAG, KEEP,KEEP8, IERROR )
            IF ( IFLAG .LT. 0 ) RETURN
           ELSE
             PTRIST(STEP(IROOT)) = -55555
           ENDIF
        END IF
        LREQI = NROW_SON + NCOL_SON
        LREQA = NROW_SON * NCOL_SON
        IF ( (LREQA.NE.0) .AND.
     *       (PTRIST(STEP(IROOT)).LT.0).AND.
     *       KEEP(60)==0) THEN
         WRITE(*,*) ' Error in SMUMPS_272 '
         CALL SMUMPS_ABORT()
        ENDIF
        IF (LREQA.NE.0) THEN
          CALL SMUMPS_22(.FALSE.,.FALSE.,
     *     MYID,N,KEEP,KEEP8,IW,LIW,A, LA, NIRBDU,
     *     LRLU, IPTRLU, IWPOS, IWPOSCB, PTRIST,
     *     PTRAST, STEP, PIMASTER, PAMASTER, ITLOC,
     *     LREQI, LREQA, -1234, S_NOTFREE, .FALSE.,
     *     COMP, LRLUS, IFLAG, IERROR )
          IF ( IFLAG .LT. 0 ) RETURN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   IW( IWPOSCB + 1 ), LREQI,
     *                   MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     *                   A( IPTRLU + 1 ), LREQA,
     *                   MPI_REAL, COMM, IERR )
          IF (KEEP(60) .EQ.0) THEN
          IF ( PTRIST(STEP(IROOT)) .NE. 0 ) THEN
                 LOCAL_N  = -IW( PTRIST(STEP( IROOT )) + XSIZE    )
                 LOCAL_M  =  IW( PTRIST(STEP( IROOT )) + 1 + XSIZE)
                 POS_ROOT = PAMASTER(STEP( IROOT ))
          ELSE
                 LOCAL_N = IW( PTLUST_S(STEP( IROOT ) ) + 1 + XSIZE)
                 LOCAL_M = IW( PTLUST_S(STEP( IROOT ) ) + 2 + XSIZE)
                 POS_ROOT = PTRFAC(IW(PTLUST_S(STEP(IROOT))+4+ XSIZE))
          END IF
          CALL SMUMPS_38( NROW_SON, NCOL_SON, IW( IWPOSCB + 1 ),
     *                     IW( IWPOSCB + NROW_SON + 1 ),
     *                     A( IPTRLU + 1 ),
     *                     A( POS_ROOT ), LOCAL_M, LOCAL_N )
          ELSE
          CALL SMUMPS_38( NROW_SON, NCOL_SON, IW( IWPOSCB + 1 ),
     *                     IW( IWPOSCB + NROW_SON + 1 ),
     *                     A( IPTRLU + 1 ),
     *                     root%SCHUR_POINTER(1),
     *                     root%SCHUR_LLD , root%SCHUR_NLOC)
          ENDIF
          IWPOSCB = IWPOSCB + LREQI
          IPTRLU  = IPTRLU  + LREQA
          LRLU    = LRLU    + LREQA
          LRLUS   = LRLUS   + LREQA
          CALL SMUMPS_471(.FALSE.,.FALSE.,
     *                    LA-LRLUS,0,-LREQA,KEEP,KEEP8,LRLU)
        ENDIF
      RETURN
      END SUBROUTINE SMUMPS_272
      SUBROUTINE SMUMPS_224(NFRONT,NASS,IBEGKJI, LPIV, TIPIV,
     *    N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU,SEUIL,KEEP,KEEP8,
     *     DKEEP,PIVNUL_LIST,LPN_LIST)
      IMPLICIT NONE
      INTEGER IBEGKJI, LPIV 
      INTEGER TIPIV(LPIV)
      INTEGER NFRONT,NASS,N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      REAL A(LA) 
      REAL UU, SEUIL
      INTEGER IW(LIW) 
      INTEGER IOLDPS, POSELT
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(30)
      REAL SWOP
      INTEGER APOS,ILOC
      REAL ZERO, RMAX, AMROW, ONE
      INTEGER NPIV,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF
      INCLUDE 'mumps_headers.h'
      INTEGER ISAMAX
      INTRINSIC MAX
      DATA ZERO /0.0E0/
      DATA ONE /1.0E0/
        NPIV    = IW(IOLDPS+1+XSIZE)
        HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
        NPIVP1  = NPIV + 1
        ILOC    = NPIVP1 - IBEGKJI + 1
        TIPIV(ILOC) = ILOC
        NASSW   = IABS(IW(IOLDPS+3+XSIZE))
        IF(INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT*(NPIVP1-1) + NPIV
           IDIAG = APOS
           IF(ABS(A(APOS)).LT.SEUIL) THEN
              IF(REAL(A(APOS)) .GE. ZERO) THEN
                 A(APOS) = SEUIL
              ELSE
                 A(APOS) = -SEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = ZERO
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
             J3    = NASS -NPIV
             JMAX  = ISAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + IPIV - NPIVP1
            IF (RMAX.LE.DKEEP(1)) THEN
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE+
     &                      IW(IOLDPS+5+XSIZE)+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(DKEEP(2).GT.ZERO) THEN
                  IF(REAL(A(IDIAG)) .GE. ZERO) THEN
                     A(IDIAG) = DKEEP(2)
                  ELSE
                     A(IDIAG) = -DKEEP(2)
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS - NPIV + NFRONT - 1
                 DO JJ=J1,J2
                   A(JJ)= REAL(ZERO)
                 ENDDO
                 A(IDIAG) = REAL(ONE)
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 380   
            ENDIF
            IF (ABS(A(IDIAG)).GT.MAX(UU*RMAX,SEUIL)) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF (AMROW.LE.MAX(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV*NFRONT
            J2 = J1 + NFRONT - 1
            J3 = POSELT + (IPIV-1)*NFRONT
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            TIPIV(ILOC) = ILOC + JMAX - 1
            J1 = POSELT + NPIV
            J2 = POSELT + NPIV + JMAX - 1
            DO 410 KSW=1,NASS
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT
              J2 = J2 + NFRONT
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      IFLAG = -10
      WRITE(*,*) 'NIV2:Detected 0 pivot, INODE,NPIV=',INODE,NPIV
  420 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_224
      SUBROUTINE  SMUMPS_294( COMM_LOAD, ASS_IRECV, 
     *             N, INODE, FPERE,
     *             IW, LIW, 
     *             IOLDPS, POSELT, A, LA, LDA_FS, 
     *             IBEGKJI, IEND, TIPIV, LPIV, LASTBL, NB_BLOC_FAC,
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
      USE SMUMPS_BUFFER
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mpif.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, IBEGKJI, IEND, LPIV, 
     *        IOLDPS, POSELT, LA, LDA_FS, NB_BLOC_FAC
      INTEGER IW(LIW), TIPIV(LPIV)
      LOGICAL LASTBL
      REAL A(LA)
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT ) 
      INTEGER KEEP(500)
      INTEGER*8 KEEP8(150)
      INTEGER NBFIN, IFLAG, IERROR, LEAF, LPOOL,
     *        SLAVEF, ICNTL(40)
      INTEGER POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, 
     *        NIRBDU, COMP
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     *        ITLOC(N), FILS(N),
     *        PTRARW(LPTRAR), PTRAIW(LPTRAR), 
     *        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER INTARR(MAX(1,KEEP(14)))
      INTEGER PTRIST(KEEP(28)), PTRAST(KEEP(28)), PTLUST_S(KEEP(28)),
     *        PTRFAC(KEEP(28)),
     *        STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28)), NSTK_S(KEEP(28)),
     *        NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     *        TAB_POS_IN_PERE(SLAVEF+2,MAX(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      REAL DBLARR(MAX(1,KEEP(13)))
      EXTERNAL  SMUMPS_329
      INCLUDE 'mumps_headers.h'
      INTEGER NPIV, NCOL, APOS, PDEST, NSLAVES
      INTEGER IERR, IERR_MPI, LREQA, LREQI
      INTEGER STATUS( MPI_STATUS_SIZE )
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      DOUBLE PRECISION FLOP1,FLOP2
      NSLAVES= IW(IOLDPS+5+XSIZE)
          IF (NSLAVES.EQ.0) THEN
           WRITE(6,*) ' ERROR 1 in SMUMPS_294 '
           CALL SMUMPS_ABORT()
          ENDIF
      NPIV   = IEND - IBEGKJI + 1
      NCOL   = LDA_FS - IBEGKJI + 1
      APOS   = POSELT + LDA_FS*(IBEGKJI-1) + IBEGKJI - 1
      IF (IBEGKJI > 0) THEN
       CALL SMUMPS_511( LDA_FS, IBEGKJI-1, LPIV,
     *                            KEEP(50),2,FLOP1)
      ELSE
        FLOP1=0
      ENDIF
      CALL SMUMPS_511( LDA_FS, IEND, LPIV,
     *                           KEEP(50),2,FLOP2)
      FLOP2 = FLOP1 - FLOP2
      CALL SMUMPS_190(1, .FALSE., FLOP2, KEEP,KEEP8)
      IF ((NPIV.GT.0) .OR. 
     *    ((NPIV.EQ.0).AND.(LASTBL)) ) THEN
        PDEST  = IOLDPS + 6 + XSIZE
        IERR = -1
        IF ( NPIV .NE. 0 ) THEN
          NB_BLOC_FAC = NB_BLOC_FAC + 1
        END IF
        DO WHILE (IERR .EQ.-1)
          CALL SMUMPS_65( INODE, LDA_FS, NCOL, 
     *               NPIV, FPERE, LASTBL, TIPIV, A(APOS),
     *               IW(PDEST), NSLAVES, KEEP(50), NB_BLOC_FAC,
     *               COMM, IERR )
        IF (IERR.EQ.-1) THEN
           BLOCKING  = .FALSE.
           SET_IRECV = .TRUE.
           MESSAGE_RECEIVED = .FALSE.
           CALL SMUMPS_329( COMM_LOAD, ASS_IRECV, 
     *      BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     *      MPI_ANY_SOURCE, MPI_ANY_TAG,
     *      STATUS, BUFR, LBUFR,
     *      LBUFR_BYTES,
     *      PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     *      LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *      PTLUST_S, PTRFAC,
     *      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, IFLAG,
     *      IERROR, COMM,
     *      NBPROCFILS,
     *      IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     *      root, OPASSW, OPELIW, ITLOC, FILS, PTRARW, PTRAIW,
     *      INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     *      LPTRAR, NELT, FRTPTR, FRTELT, 
     *      ISTEP_TO_INIV2, TAB_POS_IN_PERE )
           IF (MESSAGE_RECEIVED) POSELT = PTRAST(STEP(INODE))
           IF ( IFLAG .LT. 0 ) GOTO 500
         ENDIF
        ENDDO
        IF (IERR .EQ. -2 .OR. IERR.EQ.-3 ) THEN
          IF (IERR.EQ.-2) IFLAG = -17
          IF (IERR.EQ.-3) IFLAG = -20
          LREQA = NCOL*NPIV  
          LREQI = NPIV + 6 + 2*NSLAVES
          IERROR =  LREQI  * KEEP( 34 )+ LREQA * KEEP( 35 )
          GOTO 300
        ENDIF
      ENDIF
      GOTO 500
  300 CONTINUE
      CALL SMUMPS_44( MYID, SLAVEF, COMM )
 500  RETURN
      END SUBROUTINE  SMUMPS_294
      SUBROUTINE SMUMPS_273( ROOT, 
     *    INODE, NELIM, NSLAVES, ROW_LIST,
     *    COL_LIST, SLAVE_LIST, 
     *
     *    PROCNODE_STEPS, IWPOS, IWPOSCB, IPTRLU,
     *    LRLU, LRLUS, N, IW, LIW, A, LA, NIRBDU, PTRIST,
     *    PTLUST_S, PTRFAC,
     *    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, ITLOC, COMP,
     *    IFLAG, IERROR, 
     *    IPOOL, LPOOL, LEAF, MYID, SLAVEF, KEEP,KEEP8,
     *    COMM,COMM_LOAD,FILS,ND )
      USE SMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      TYPE (SMUMPS_ROOT_STRUC) :: ROOT
      INTEGER INODE, NELIM, NSLAVES 
      INTEGER KEEP( 500 )
      INTEGER*8 KEEP8(150)
      INTEGER ROW_LIST(*), COL_LIST(*), 
     *        SLAVE_LIST(*)
      INTEGER IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, NIRBDU
      INTEGER N, LIW, LA
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER PTRIST( KEEP(28) ), PTLUST_S(KEEP(28)), PTRFAC(KEEP(28))
      INTEGER PTRAST(KEEP(28))
      INTEGER STEP(N), 
     * PIMASTER(KEEP(28)),
     *  PAMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S(KEEP(28)), ITLOC( N ), PROCNODE_STEPS( KEEP(28) )
      INTEGER IFLAG, IERROR
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF
      INTEGER COMM,COMM_LOAD,ND(KEEP(28)),FILS(N)
      INTEGER IROOT, TYPE_INODE, DEB_ROW, DEB_COL,
     *        NOINT, NOREAL
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mumps_tags.h'
      INTEGER SMUMPS_330
      EXTERNAL SMUMPS_330
      IROOT        = KEEP(38)
      NSTK_S(STEP(IROOT))= NSTK_S(STEP(IROOT)) - 1
      KEEP(42) = KEEP(42) + NELIM
      TYPE_INODE= SMUMPS_330(STEP(INODE), PROCNODE_STEPS, SLAVEF)
      IF (TYPE_INODE.EQ.1) THEN 
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + 1
        ELSE 
         KEEP(41) = KEEP(41) + 3
        ENDIF
      ELSE
        IF (NELIM.EQ.0) THEN
         KEEP(41) = KEEP(41) + NSLAVES
        ELSE 
         KEEP(41) = KEEP(41) + 2*NSLAVES + 1
        ENDIF
      ENDIF
      IF  (NELIM.EQ.0) THEN
        PIMASTER(STEP(INODE)) = 0 
      ELSE
       NOINT = 6 + NSLAVES + NELIM  + NELIM + XSIZE
       NOREAL= 0
       CALL SMUMPS_22(.FALSE.,.FALSE.,
     *   MYID,N,KEEP,KEEP8,IW,LIW, A, LA, NIRBDU,
     *   LRLU, IPTRLU,IWPOS,IWPOSCB,
     *   PTRIST,PTRAST,STEP, PIMASTER, PAMASTER, ITLOC,
     *   NOINT, NOREAL, INODE, S_NOTFREE, .TRUE.,
     *   COMP, LRLUS, IFLAG, IERROR )
       IF ( IFLAG .LT. 0 ) THEN
         WRITE(*,*) ' Failure in int space allocation in CB area ',
     *    ' during assembly of root : SMUMPS_273',
     *    ' size required was :', NOINT,
     *    'INODE=',INODE,' NELIM=',NELIM, ' NSLAVES=', NSLAVES
         RETURN
        ENDIF
        PIMASTER(STEP( INODE )) = IWPOSCB + 1
        PAMASTER(STEP( INODE )) = IPTRLU  + 1
        IW( IWPOSCB + 1+XSIZE ) = 2*NELIM
        IW( IWPOSCB + 2+XSIZE ) = NELIM
        IW( IWPOSCB + 3+XSIZE ) = 0
        IW( IWPOSCB + 4+XSIZE ) = 0
        IW( IWPOSCB + 5+XSIZE ) = 1
        IW( IWPOSCB + 6+XSIZE ) = NSLAVES
        IF (NSLAVES.GT.0) THEN
         IW( IWPOSCB+7+XSIZE: IWPOSCB+7+XSIZE+NSLAVES-1) = 
     *                   SLAVE_LIST(1:NSLAVES)
        ENDIF
        DEB_ROW = IWPOSCB+7+NSLAVES+XSIZE
        IW(DEB_ROW : DEB_ROW+NELIM -1) = ROW_LIST(1:NELIM)
        DEB_COL = DEB_ROW + NELIM
        IW(DEB_COL : DEB_COL+NELIM -1) = COL_LIST(1:NELIM)
      ENDIF
      IF (NSTK_S(STEP(IROOT)) .EQ. 0 ) THEN
          CALL SMUMPS_507(N, IPOOL, LPOOL, PROCNODE_STEPS,
     *         SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     *         STEP, IROOT )
          IF (KEEP(47) .GE. 3) THEN
             CALL SMUMPS_500(
     $            IPOOL, LPOOL, 
     *            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     *            MYID, STEP, N, ND, FILS )
          ENDIF
      END IF
      RETURN
      END SUBROUTINE SMUMPS_273
      SUBROUTINE SMUMPS_534( N,FRERE, FILS,
     *     NA,NE,ND,PERM,SYM,INFO,LP,K47,K81,K215,
     *     PROCNODE,SLAVEF,PEAK
     $     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, LP, SIZE_MEM_SBTR
      INTEGER FRERE(N), FILS(N)
      INTEGER NA(N), NE(N), ND(N),K47,K81,K215
      INTEGER INFO(40)
      INTEGER SLAVEF,PROCNODE(N)
      REAL PEAK
      INTEGER NBROOT, NBLEAF, LNA, allocok, LEAF, I, NSTEPS,
     *        K47_LOC, K81_LOC
      INTEGER, ALLOCATABLE, DIMENSION (:) :: NEW_NA, STEP
      INTEGER TEMP_MEM(1),SBTR_WHICH_M
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: MEM_SUBTREE,MY_ROOT,
     $     MY_SIZE,MY_LEAF
      INTEGER, ALLOCATABLE, DIMENSION (:) ::      DEPTH_FIRST
      REAL, ALLOCATABLE, DIMENSION (:) :: COST_TRAV
      INTEGER DUMMY_DAD(1), DUMMY_DAD_LENGTH
      PARAMETER (DUMMY_DAD_LENGTH=1)
      LOGICAL USE_DAD
      PARAMETER (USE_DAD=.FALSE.)
      INCLUDE 'mumps_headers.h'
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF = 1 
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
      ENDIF
      LNA = NBROOT + NBLEAF + 2
      ALLOCATE (NEW_NA (LNA), STEP(N), stat=allocok)
      IF (allocok > 0 ) THEN
        INFO(1) = -7
        INFO(2) = LNA + N
        RETURN
      ENDIF
      NEW_NA(1)=NBLEAF
      NEW_NA(2)=NBROOT
      LEAF = 3
      IF ( N == 1 ) THEN
          NEW_NA(LEAF) = 1
          LEAF = LEAF + 1
      ELSE IF (NA(N) < 0) THEN
          NEW_NA(LEAF) = - NA(N)-1
          LEAF = LEAF + 1
          DO I = 1, NBLEAF - 1
            NEW_NA(LEAF) = NA(I)
            LEAF = LEAF + 1
          ENDDO
      ELSE IF (NA(N-1) < 0 ) THEN
          NEW_NA(LEAF) = - NA(N-1) - 1
          LEAF =LEAF + 1
          IF ( NBLEAF > 1 ) THEN
            DO I = 1, NBLEAF - 1
              NEW_NA(LEAF) = NA(I)
              LEAF = LEAF + 1
            ENDDO
          ENDIF
      ELSE
          DO I = 1, NBLEAF
            NEW_NA(LEAF) = NA(I)
            LEAF = LEAF + 1
          ENDDO
      END IF
      SIZE_MEM_SBTR=NEW_NA(2)
      ALLOCATE(MEM_SUBTREE(SIZE_MEM_SBTR,SLAVEF),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= SIZE_MEM_SBTR*SLAVEF
         RETURN
      ENDIF
      ALLOCATE(MY_ROOT(SIZE_MEM_SBTR,SLAVEF),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= SIZE_MEM_SBTR*SLAVEF
         RETURN
      ENDIF
      ALLOCATE(MY_LEAF(SIZE_MEM_SBTR,SLAVEF),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= SIZE_MEM_SBTR*SLAVEF
         RETURN
      ENDIF
      ALLOCATE(MY_SIZE(SIZE_MEM_SBTR,SLAVEF),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= SIZE_MEM_SBTR*SLAVEF
         RETURN
      ENDIF
      DO I = 1, N
        STEP(I)=I
        IF ( FRERE(I) .EQ. 0 ) THEN
          NEW_NA(LEAF) = I
          LEAF = LEAF + 1
        END IF
      END DO
      NSTEPS = N
      K47_LOC = 0
      K81_LOC = 0
      SBTR_WHICH_M=0
      ALLOCATE(DEPTH_FIRST(1),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= 1
         RETURN
      ENDIF
      ALLOCATE(COST_TRAV(1),stat=allocok) 
      IF (allocok .ne.0) THEN
         info(1)= -7
         info(2)= 1
         RETURN
      ENDIF
      CALL SMUMPS_363(N,FRERE, STEP, FILS,
     *     NEW_NA,LNA,NE,ND,
     *     DUMMY_DAD, DUMMY_DAD_LENGTH, USE_DAD,
     *     NSTEPS,PERM,SYM,INFO,LP,K47_LOC,K81_LOC,0,K215,
     *     PROCNODE,MEM_SUBTREE,SLAVEF, SIZE_MEM_SBTR, PEAK,SBTR_WHICH_M
     $     ,1,1,DEPTH_FIRST,COST_TRAV,MY_LEAF,MY_SIZE,MY_ROOT
     *)
      NA(1:NBLEAF)=NEW_NA(3:2+NBLEAF)
      NA(N)=NBROOT
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
      DEALLOCATE(MEM_SUBTREE)
      DEALLOCATE(NEW_NA,STEP)
      DEALLOCATE(DEPTH_FIRST,COST_TRAV)
      DEALLOCATE(MY_LEAF,MY_SIZE,MY_ROOT)
      RETURN
      END SUBROUTINE SMUMPS_534
      SUBROUTINE SMUMPS_363(N,FRERE, STEP, FILS,
     *     NA,LNA,NE,ND, DAD, LDAD, USE_DAD,
     *     NSTEPS,PERM,SYM,INFO,LP,K47,K81,K76,K215,
     *     PROCNODE,MEM_SUBTREE,SLAVEF, SIZE_MEM_SBTR, PEAK
     $     ,SBTR_WHICH_M,SIZE_DEPTH_FIRST,SIZE_COST_TRAV,
     $     DEPTH_FIRST_TRAV,COST_TRAV,MY_FIRST_LEAF,
     $     MY_NB_LEAF,MY_ROOT_SBTR
     $     )
      IMPLICIT NONE
      INTEGER N,PERM,SYM, NSTEPS, LNA, LP, SIZE_MEM_SBTR,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS),K47,K81,K76,K215
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(40)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      INTEGER MEM_SUBTREE(SIZE_MEM_SBTR,SLAVEF),SBTR_WHICH_M
      INTEGER MY_FIRST_LEAF(SIZE_MEM_SBTR,SLAVEF),
     $     MY_ROOT_SBTR(SIZE_MEM_SBTR,SLAVEF),
     $     MY_NB_LEAF(SIZE_MEM_SBTR,SLAVEF)
      EXTERNAL SMUMPS_283,SMUMPS_275
      LOGICAL SMUMPS_283
      INTEGER SMUMPS_275
      REAL PEAK
      INTEGER SIZE_DEPTH_FIRST,DEPTH_FIRST_TRAV(SIZE_DEPTH_FIRST)
      INTEGER SIZE_COST_TRAV
      REAL COST_TRAV(SIZE_COST_TRAV)
      INTEGER, DIMENSION (:), ALLOCATABLE :: DEPTH
      INTEGER ITOP,IFATH,IN,LSTK,NSTK,INODE,K,I,allocok,LOCAL_PERM
      INTEGER*8 NCB
      INTEGER*8 NELIM,NFR
      INTEGER NFR4,NELIM4
      INTEGER LEAF,NBLEAF,NBROOT, SIZE_TAB
      INTEGER, DIMENSION (:), ALLOCATABLE :: ROOT,IPOOL,TNSTK
      INTEGER, DIMENSION (:), ALLOCATABLE,TARGET :: SON,TEMP
      INTEGER*8, DIMENSION (:), ALLOCATABLE :: M,M_TOTAL, FACT
      INTEGER*8, DIMENSION (:), ALLOCATABLE :: TAB1,TAB2
      INTEGER, DIMENSION (:), POINTER :: TAB
      INTEGER J,x,y,z,dernier,fin,RANK_TRAV
      INTEGER cour,t,II,temporary
      INTEGER actuel,CB_current,CB_MAX,ROOT_OF_CUR_SBTR
      INTEGER*8, DIMENSION (:), ALLOCATABLE :: T1,T2
      INTEGER, DIMENSION (:), ALLOCATABLE :: RESULT
      INTEGER*8 MEM_SIZE,FACT_SIZE,SUM,MEM_SEC_PERM,FACT_SIZE_T,
     $     MEM_SIZE_T,TOTAL_MEM_SIZE,TMP_TOTAL_MEM_SIZE,TMP_SUM,
     $     SIZECB
      INTEGER*8 TMP4
      LOGICAL   SBTR_M
      INTEGER INDICE(SLAVEF),ID,FIRST_LEAF,SIZE_SBTR
      EXTERNAL SMUMPS_170,SMUMPS_167
      LOGICAL SMUMPS_170,SMUMPS_167
      DOUBLE PRECISION COST_NODE
      INCLUDE 'mumps_headers.h'
      TOTAL_MEM_SIZE=0
      ROOT_OF_CUR_SBTR=0
      IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.
     $     (PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4).OR.
     $     (PERM.EQ.5).OR.(PERM.EQ.6))THEN
         LOCAL_PERM=0
      ENDIF
      IF (K47 == 4 .OR. ((K47.GE.2).AND.(K81.GE. 1))) THEN
        DO I=1,SLAVEF
          INDICE(I)=1
        ENDDO
        DO I=1,SLAVEF
          DO x=1,SIZE_MEM_SBTR
            MEM_SUBTREE(x,I)=-1
          ENDDO
        ENDDO
      ENDIF
      SBTR_M=((K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))
      MEM_SIZE=0
      FACT_SIZE=0
      IF ((PERM.GT.7).AND.
     * (.NOT.(K47 == 4 .OR. ((K47.GE.2).AND.(K81 .GE. 1))))) THEN
         WRITE(*,*) "Internal Error in SMUMPS_363",PERM
         CALL SMUMPS_ABORT()
      END IF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      IF((PERM.EQ.0).AND.(NBROOT.EQ.NBLEAF)) RETURN
      IF (SBTR_M.OR.(PERM.EQ.2))  THEN
         IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
            ALLOCATE(M_TOTAL(NSTEPS), stat=allocok )
            IF (allocok > 0) THEN
               IF ( LP .GT. 0 )
     *              WRITE(LP,*)'Memory allocation error in
     $              SMUMPS_363'
               INFO(1)=-7
               INFO(2)=NSTEPS
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ALLOCATE( IPOOL(NBLEAF), M(NSTEPS), fact(NSTEPS),
     *          TNSTK(NSTEPS), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     *    WRITE(LP,*)'Memory allocation error in SMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      II=0
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
         IF(NE(I).GE.II) II=NE(I)
      ENDDO
      SIZE_TAB=MAX(II,NBROOT)
      ALLOCATE(SON(II), TEMP(II),
     *         TAB1(SIZE_TAB), TAB2(SIZE_TAB), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     *    WRITE(LP,*)'Memory allocation error in SMUMPS_363'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      ALLOCATE(T1(SIZE_TAB),T2(SIZE_TAB),
     *         RESULT(SIZE_TAB),stat=allocok)
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     *    WRITE(LP,*)'Memory allocation error in SMUMPS_363'
        INFO(1)=-7
        INFO(2)=SIZE_TAB
        RETURN
      ENDIF
      IF(NBROOT.EQ.NBLEAF)THEN
        IF((PERM.NE.1).OR.(PERM.EQ.4).OR.(PERM.EQ.6))THEN
          WRITE(*,*)'Internal Error in reordertree:'
          WRITE(*,*)'  problem with perm parameter in reordertree'
          CALL SMUMPS_ABORT()
        ENDIF
        DO I=1,NBROOT
          TAB1(I)=ND(STEP(NA(I+2+NBLEAF)))
          IPOOL(I)=NA(I+2+NBLEAF)
          M(STEP(IPOOL(I)))=TAB1(I)*TAB1(I)
        ENDDO
        CALL SMUMPS_462(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &    RESULT,T1,T2)
        GOTO 789
      ENDIF
      IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
         ALLOCATE(DEPTH(NSTEPS),stat=allocok)
         IF (allocok > 0) THEN
            IF ( LP .GT. 0 )
     *           WRITE(LP,*)'Memory allocation error in
     $           SMUMPS_363'
            INFO(1)=-7
            INFO(2)=NSTEPS
            RETURN
         ENDIF
         DEPTH=0
         NBROOT = NA(2)
         IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
         fin=NBROOT
         LEAF=NA(1)
 499     CONTINUE
         INODE=IPOOL(fin)
         IF(INODE.LT.0)THEN
            WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
            CALL SMUMPS_ABORT()
         ENDIF
         IN=INODE
 4602    IN = FILS(IN)
         IF (IN .GT. 0 ) THEN
            GOTO 4602
         ENDIF
         IN=-IN
         DO I=1,NE(STEP(INODE))
            SON(I)=IN
            IN=FRERE(STEP(IN))
         ENDDO
         DO I=1,NE(STEP(INODE))
            IPOOL(fin)=SON(I)
            DEPTH(STEP(SON(I)))=DEPTH(STEP(INODE))+1
            SON(I)=0
            fin=fin+1
         ENDDO
         IF(NE(STEP(INODE)).EQ.0)THEN
            LEAF=LEAF-1
         ELSE
            fin=fin-1
            GOTO 499
         ENDIF
         fin=fin-1
         IF(fin.EQ.0) GOTO 489
         GOTO 499
 489     CONTINUE
      ENDIF
      IF(K76.EQ.4)THEN
         RANK_TRAV=NSTEPS
         DEPTH_FIRST_TRAV=0
      ENDIF
      IF((K76.EQ.5).OR.(PERM.EQ.5).OR.(PERM.EQ.6))THEN
         COST_TRAV=0.0E0
         COST_NODE=0.0E0
      ENDIF        
      DO I=1,NSTEPS
         M(I)=0
         IF (SBTR_M.OR.(PERM.EQ.2))  THEN
            IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
               M_TOTAL(I)=0
            ENDIF
         ENDIF
      ENDDO
      DO I=1,NSTEPS
         fact(I)=0
      ENDDO
      IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
      LEAF = NBLEAF + 1
 91   CONTINUE
        IF (LEAF.NE.1) THEN
           LEAF = LEAF -1
           INODE = IPOOL(LEAF)
        ENDIF
 96     CONTINUE
        NFR    = ND(STEP(INODE))
        NSTK   = NE(STEP(INODE))
        NELIM4 = 0
        IN = INODE
 101    NELIM4 = NELIM4 + 1
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 101
        NELIM=NELIM4
        IF(NE(STEP(INODE)).EQ.0) THEN
           M(STEP(INODE))=NFR*NFR
           IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                 M_TOTAL(STEP(INODE))=NFR*NFR
           ENDIF
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.3))THEN
           IF(SMUMPS_170(STEP(INODE),
     $PROCNODE,SLAVEF))THEN
              DEPTH(STEP(INODE))=0
           ENDIF
        ENDIF
        IF ( SYM .eq. 0 ) THEN
          FACT(STEP(INODE))=FACT(STEP(INODE))+
     &      (2*NFR*NELIM)-(NELIM*NELIM)
        ELSE
          FACT(STEP(INODE))=FACT(STEP(INODE))+NFR*NELIM
        ENDIF
        IF (USE_DAD) THEN
          IFATH = DAD( STEP(INODE) )
        ELSE
          IN = INODE
 113      IN = FRERE(IN)
          IF (IN.GT.0) GO TO 113
          IFATH = -IN
        ENDIF
        IF (IFATH.EQ.0) THEN
           NBROOT = NBROOT - 1
           IF (NBROOT.EQ.0) GOTO 116
           GOTO 91
        ELSE
           FACT(STEP(IFATH))=FACT(STEP(IFATH))+FACT(STEP(INODE))
           IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
              DEPTH(STEP(IFATH))=MAX(DEPTH(STEP(INODE)),
     $             DEPTH(STEP(IFATH)))
           ENDIF
        ENDIF
        TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
        IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
           INODE = IFATH        
           IN=INODE
           dernier=IN
           I=1
 5700      IN = FILS(IN)
           IF (IN .GT. 0 ) THEN
             dernier=IN
             I=I+1
             GOTO 5700
           ENDIF
           NCB=ND(STEP(INODE))-I
           IN=-IN
           IF(PERM.NE.7)THEN
              DO I=1,NE(STEP(INODE))
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ELSE
              DO I=NE(STEP(INODE)),1,-1
                 SON(I)=IN
                 TEMP(I)=IN
                 IF(IN.GT.0) IN=FRERE(STEP(IN))
              ENDDO
           ENDIF
           IF(PERM.EQ.7) GOTO 213
           DO II=1,NE(STEP(INODE))
             TAB1(II)=0
             TAB2(II)=0
             cour=SON(II)
             NELIM4=1
 151         cour=FILS(cour)
             IF(cour.GT.0) THEN
                NELIM4=NELIM4+1
                GOTO 151
             ENDIF
             NELIM=NELIM4
             IF((SYM.EQ.0).OR.(K215.NE.0)) THEN
                SIZECB=(ND(STEP(SON(II)))-NELIM)*
     $               (ND(STEP(SON(II)))-NELIM)
             ELSE
                SIZECB=(ND(STEP(SON(II)))-NELIM)*(ND(STEP(SON(II)))-
     $               NELIM+1)/2
             ENDIF
             IF((PERM.EQ.0).OR.(PERM.EQ.5))THEN
               TAB1(II)=M(STEP(SON(II)))-
     &           SIZECB
               TAB2(II)=SIZECB
             ENDIF
             IF((PERM.EQ.1).OR.(PERM.EQ.6)) THEN
                TAB1(II)=M(STEP(SON(II)))-SIZECB
                TAB1(II)=TAB1(II)-FACT(STEP(SON(II)))
                TAB2(II)=SIZECB+FACT(STEP(SON(II)))
             ENDIF
             IF(PERM.EQ.2)THEN
                IF (SMUMPS_170(STEP(INODE),
     $               PROCNODE,SLAVEF))THEN
                   TAB1(II)=M_TOTAL(STEP(SON(II)))-SIZECB
     $                  -FACT(STEP(SON(II)))
                   TAB2(II)=SIZECB
                ELSE
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB             
                ENDIF
             ENDIF
             IF(PERM.EQ.3)THEN
                IF (SMUMPS_170(STEP(INODE),
     $               PROCNODE,SLAVEF))THEN
                   TAB1(II)=M(STEP(SON(II)))-SIZECB
                   TAB2(II)=SIZECB               
                ELSE
                   TAB1(II)=DEPTH(STEP(SON(II)))
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
             IF(PERM.EQ.4)THEN
                IF (SMUMPS_170(STEP(INODE),
     $               PROCNODE,SLAVEF))THEN
                   TAB1(II)=M(STEP(SON(II)))-
     &                  SIZECB-FACT(STEP(SON(II)))
                   TAB2(II)=SIZECB             
                ELSE
                   TAB1(II)=DEPTH(STEP(SON(II)))
                   TAB2(II)=M(STEP(SON(II)))
                ENDIF
             ENDIF
          ENDDO
          CALL SMUMPS_462(SON,NE(STEP(INODE)),TAB1,TAB2,
     $         LOCAL_PERM
     &           ,RESULT,T1,T2)
          IF(PERM.EQ.0) THEN
             DO II=1,NE(STEP(INODE))
               cour=TEMP(II)
               NELIM4=1
 153           cour=FILS(cour)
               IF(cour.GT.0) THEN
                  NELIM4=NELIM4+1
                  GOTO 153
               ENDIF
               NELIM=NELIM4
               IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                  SIZECB=(ND(STEP(TEMP(II)))-NELIM)*
     $                 (ND(STEP(TEMP(II)))-NELIM)
               ELSE
                  SIZECB=(ND(STEP(TEMP(II)))-NELIM)*
     $                 (ND(STEP(TEMP(II)))-NELIM+1)/2
               ENDIF
               TAB1(II)=SIZECB
             ENDDO
             CALL SMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)             
           ENDIF
           IF(PERM.EQ.1) THEN
              DO II=1,NE(STEP(INODE))
                cour=TEMP(II)
                NELIM4=1
 187            cour=FILS(cour)
                IF(cour.GT.0) THEN
                   NELIM4=NELIM4+1
                   GOTO 187
                ENDIF    
                NELIM=NELIM4   
                IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                   SIZECB=(ND(STEP(TEMP(II)))-NELIM)*
     $                  (ND(STEP(TEMP(II)))-NELIM)
                ELSE
                   SIZECB=(ND(STEP(TEMP(II)))-NELIM)*
     $                  (ND(STEP(TEMP(II)))-NELIM+1)/2
                ENDIF
                TAB1(II)=SIZECB+FACT(STEP(TEMP(II)))
             ENDDO
             CALL SMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,3,
     &         RESULT,T1,T2)
           ENDIF
 213       CONTINUE
           IFATH=INODE
           DO II=1,2
              SUM=0
              FACT_SIZE=0
              FACT_SIZE_T=0
              MEM_SIZE=0
              MEM_SIZE_T=0
              CB_MAX=0
              CB_CURRENT=0
              TMP_SUM=0
              IF(II.EQ.1) TAB=>SON 
              IF(II.EQ.2) TAB=>TEMP
              DO I=1,NE(STEP(INODE))
                 cour=TAB(I)
                 NELIM4=1
 149             cour=FILS(cour)
                 IF(cour.GT.0) THEN
                    NELIM4=NELIM4+1
                    GOTO 149
                 ENDIF    
                 NELIM=NELIM4   
                 NFR=ND(STEP(TAB(I)))
                 IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                    SIZECB=(NFR-NELIM)*(NFR-NELIM)
                 ELSE
                    SIZECB=(NFR-NELIM)*(NFR-NELIM+1)/2
                 ENDIF
                 MEM_SIZE=MAX(MEM_SIZE,(M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 IF (SBTR_M.OR.(PERM.EQ.2)) THEN
                       MEM_SIZE_T=MAX(MEM_SIZE_T,(M_TOTAL(STEP(TAB(I)))+
     $                      SUM+
     $                      FACT_SIZE_T))
                       FACT_SIZE_T=FACT_SIZE_T+fact(STEP(TAB(I)))
                 ENDIF
                 TOTAL_MEM_SIZE=MAX(TOTAL_MEM_SIZE,
     $                (M(STEP(TAB(I)))+SUM+FACT_SIZE))
                 TMP_SUM=TMP_SUM+fact(STEP(TAB(I)))
                 SUM=SUM+SIZECB
                 IF((PERM.EQ.1).OR.(PERM.EQ.4))THEN
                    FACT_SIZE=FACT_SIZE+fact(STEP(TAB(I)))
                 ENDIF
              ENDDO
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=NCB*NCB
              ELSE
                 SIZECB=NCB*(NCB+1)/2
              ENDIF
              TOTAL_MEM_SIZE=MAX(TOTAL_MEM_SIZE,((ND(STEP(IFATH))
     $             *ND(STEP(IFATH)))+MAX(SUM,SIZECB)+TMP_SUM))
              IF(II.EQ.1)THEN
                 TMP_TOTAL_MEM_SIZE=TOTAL_MEM_SIZE
              ENDIF
              IF((II.EQ.1).OR.(PERM.EQ.7)) THEN
                 M(STEP(IFATH))=MAX(MEM_SIZE,((ND(STEP(IFATH))
     $             *ND(STEP(IFATH)))+MAX(SUM,SIZECB)+FACT_SIZE))
                 IF (SBTR_M.OR.(PERM.EQ.2))  THEN
                       M_TOTAL(STEP(IFATH))=MAX(MEM_SIZE_T,
     $                      ((ND(STEP(IFATH))
     $                      *ND(STEP(IFATH)))+MAX(SUM,SIZECB)+
     $                      FACT_SIZE_T))
                 ENDIF
              ENDIF
              IF((II.EQ.2).AND.(PERM.EQ.1).OR.(PERM.EQ.0).OR.
     $             (PERM.EQ.5).OR.(PERM.EQ.6).OR.
     $             (.NOT.SBTR_M.OR.(SBTR_WHICH_M.NE.1)))THEN
                 MEM_SEC_PERM=MAX(MEM_SIZE,((ND(STEP(IFATH))
     $             *ND(STEP(IFATH)))+MAX(SUM,SIZECB)+FACT_SIZE))
              ENDIF
              IF((PERM.EQ.2).OR.(PERM.EQ.3).OR.(PERM.EQ.4))THEN
                 MEM_SEC_PERM=HUGE(MEM_SEC_PERM)
              ENDIF
           ENDDO
           IF(PERM.EQ.7) GOTO 96
           IF(MEM_SEC_PERM.EQ.M(STEP(IFATH))) THEN
              TAB=>TEMP
           ELSE IF (MEM_SEC_PERM.LT.M(STEP(IFATH))) THEN
              WRITE(*,*)'Probleme dans reorder!!!!'
              CALL SMUMPS_ABORT()
           ELSE 
              TOTAL_MEM_SIZE=TMP_TOTAL_MEM_SIZE
              TAB=>SON
           ENDIF
           DO I=NE(STEP(INODE)),1,-1
              IF(I.EQ.NE(STEP(INODE))) THEN
                 FILS(dernier)=-TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
              IF(I.EQ.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 FRERE(STEP(TAB(I)))=-INODE
                 GOTO 222
              ENDIF
              IF(I.GT.1) THEN
                 FRERE(STEP(dernier))=TAB(I)
                 dernier=TAB(I)
                 GOTO 222
              ENDIF
 222          CONTINUE
           ENDDO
           GOTO 96
        ELSE
           GOTO 91
        ENDIF
 116    CONTINUE
        NBROOT = NA(2)
        IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        IF(PERM.EQ.7) GOTO 001
        IF (PERM.eq.1) THEN
          DO I=1,NBROOT
            TAB1(I)=M(STEP(NA(I+2+NBLEAF)))-FACT(STEP(NA(I+2+NBLEAF)))
            TAB1(I)=-TAB1(I)
          ENDDO
          CALL SMUMPS_462(NA(2+NBLEAF+1),NBROOT,TAB1,TAB2,4,
     &      RESULT,T1,T2)
          IPOOL(1:NBROOT) = NA(3+NBLEAF:2+NBLEAF+NBROOT)
        ENDIF
 001    CONTINUE
        fin=NBROOT
        LEAF=NA(1)
        FIRST_LEAF=-9999
        SIZE_SBTR=0
 999    CONTINUE
        INODE=IPOOL(fin)
        IF(INODE.LT.0)THEN
           WRITE(*,*)'Internal Error in reordertree INODE < 0 !'
           CALL SMUMPS_ABORT()
        ENDIF
        IF(SIZE_SBTR.NE.0)THEN 
           IF(.NOT.SMUMPS_167(STEP(INODE),PROCNODE,SLAVEF))THEN
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
        IF(SMUMPS_283(STEP(INODE),PROCNODE,SLAVEF))THEN
           ROOT_OF_CUR_SBTR=INODE
        ENDIF
        IF (K76.EQ.4)THEN
           IF(SLAVEF.NE.1)THEN
              WRITE(*,*)'INODE=',INODE,'RANK',RANK_TRAV
              IF(SMUMPS_167(STEP(INODE),PROCNODE,SLAVEF))THEN
                 DEPTH_FIRST_TRAV(STEP(INODE))=DEPTH_FIRST_TRAV(STEP(
     $                ROOT_OF_CUR_SBTR))
              ELSE
                 DEPTH_FIRST_TRAV(STEP(INODE))=RANK_TRAV
              ENDIF
              RANK_TRAV=RANK_TRAV-1
           ENDIF
        ENDIF
        IF (K76.EQ.5)THEN
           IF(SLAVEF.NE.1)THEN
              IF (USE_DAD) THEN
                IFATH=DAD(INODE)
              ELSE
                IN = INODE
 395            IN = FRERE(IN)
                IF (IN.GT.0) GO TO 395
                IFATH = -IN
              ENDIF
              NFR    = ND(STEP(INODE))
              NELIM4 = 0
              IN = INODE
 396          NELIM4 = NELIM4 + 1
              IN = FILS(IN)
              IF (IN .GT. 0 ) GOTO 396
              NELIM=NELIM4
              IF((SYM.EQ.0).OR.(K215.NE.0))THEN
                 SIZECB=(NFR-NELIM)*(NFR-NELIM)
              ELSE
                 SIZECB=(NFR-NELIM)*(NFR-NELIM+1)/2
              ENDIF
              NFR4=NFR
              CALL SMUMPS_511(NFR4,NELIM4,NELIM4,
     *             SYM,1,COST_NODE)
              IF(IFATH.NE.0)THEN
                 IF(SMUMPS_167(STEP(INODE),PROCNODE,SLAVEF))THEN
                    COST_TRAV(STEP(INODE))=COST_TRAV(STEP(
     $                   ROOT_OF_CUR_SBTR))
                 ELSE
                    COST_TRAV(STEP(INODE))=COST_NODE+
     $                   COST_TRAV(STEP(IFATH))+
     $                   (SIZECB*18)
                 ENDIF
              ELSE
                 COST_TRAV(STEP(INODE))=COST_NODE
              ENDIF
              IF(K76.EQ.5)THEN
                 WRITE(*,*)'INODE=',INODE,'COST=',COST_TRAV(STEP(INODE))
              ENDIF
           ENDIF
        ENDIF
        IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF((SLAVEF.NE.1).AND.
     $          SMUMPS_283(STEP(INODE),PROCNODE,SLAVEF))THEN
                IF (NE(STEP(INODE)).NE.0) THEN
                   ID=SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
                   IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=
     $                     M_TOTAL(STEP(INODE))
                   ELSE
                      MEM_SUBTREE(INDICE(ID+1),ID+1)=M(STEP(INODE))
                   ENDIF
                   MY_ROOT_SBTR(INDICE(ID+1),ID+1)=INODE
                  INDICE(ID+1)=INDICE(ID+1)+1
                ENDIF
              ENDIF
              IF((SLAVEF.EQ.1).AND.FRERE(STEP(INODE)).EQ.0)THEN
                 ID=SMUMPS_275(STEP(INODE),PROCNODE,SLAVEF)
                 IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1))THEN
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=
     $                   M_TOTAL(STEP(INODE))
                 ELSE
                    MEM_SUBTREE(INDICE(ID+1),ID+1)=M(STEP(INODE))
                 ENDIF
                 INDICE(ID+1)=INDICE(ID+1)+1
              ENDIF
        ENDIF
        IN=INODE
 5602   IN = FILS(IN)
        IF (IN .GT. 0 ) THEN
           dernier=IN
           GOTO 5602
        ENDIF
        IN=-IN
        DO I=1,NE(STEP(INODE))
           TEMP(I)=IN
           IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
              IF((SLAVEF.NE.1).AND.(.NOT.SMUMPS_170(
     $             STEP(INODE),PROCNODE,SLAVEF)))THEN
                 NFR    = ND(STEP(INODE))
                 NELIM4 = 0
                 II = TEMP(I)
 845             NELIM4 = NELIM4 + 1
                 II = FILS(II)
                 IF (II .GT. 0 ) GOTO 845
                 NELIM=NELIM4
                 NFR4 = NFR
                 CALL SMUMPS_511(NFR4,NELIM4,NELIM4,
     *                SYM,1,COST_NODE)
                 TAB1(I)=COST_NODE+
     $                COST_TRAV(STEP(INODE))
                 TAB2(I)=0.0E0
              ELSE
                 SON(I)=IN
              ENDIF
           ELSE
              SON(I)=IN
           ENDIF
           IN=FRERE(STEP(IN))
        ENDDO
        IF((PERM.EQ.5).OR.(PERM.EQ.6))THEN
           IF((SLAVEF.NE.1).AND.(.NOT.SMUMPS_170(
     $          STEP(INODE),PROCNODE,SLAVEF)))THEN
              CALL SMUMPS_462(TEMP,NE(STEP(INODE)),TAB1,TAB2,
     $             LOCAL_PERM
     &             ,RESULT,T1,T2)
              TAB=>TEMP
              DO I=NE(STEP(INODE)),1,-1
                 IF(I.EQ.NE(STEP(INODE))) THEN
                    FILS(dernier)=-TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
                 IF(I.EQ.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    FRERE(STEP(TAB(I)))=-INODE
                    GOTO 221
                 ENDIF
                 IF(I.GT.1) THEN
                    FRERE(STEP(dernier))=TAB(I)
                    dernier=TAB(I)
                    GOTO 221
                 ENDIF
 221             CONTINUE
                 SON(NE(STEP(INODE))-I+1)=TAB(I)
              ENDDO
           ENDIF
        ENDIF
        DO I=1,NE(STEP(INODE))
           IPOOL(fin)=SON(I)
           SON(I)=0
           fin=fin+1
        ENDDO
        IF(NE(STEP(INODE)).EQ.0)THEN
           IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
              IF(SLAVEF.NE.1)THEN
                 IF(SMUMPS_167(STEP(INODE),PROCNODE,SLAVEF))THEN
                    IF(FIRST_LEAF.EQ.-9999)THEN
                       FIRST_LEAF=INODE
                    ENDIF
                    SIZE_SBTR=SIZE_SBTR+1
                 ENDIF
              ENDIF
           ENDIF
           IF(PERM.NE.7)THEN
              NA(LEAF+2)=INODE
           ENDIF
           LEAF=LEAF-1
        ELSE
           fin=fin-1
           GOTO 999
        ENDIF
        fin=fin-1
        IF(fin.EQ.0) THEN
           IF(SIZE_SBTR.NE.0)THEN 
              IF ( K47 == 4 .OR. ((K81.GE.1).AND.(K47.GE.2))) THEN
                 IF((SLAVEF.NE.1))THEN
                    MY_FIRST_LEAF(INDICE(ID+1)-1,ID+1)=FIRST_LEAF
                    MY_NB_LEAF(INDICE(ID+1)-1,ID+1)=SIZE_SBTR
                    FIRST_LEAF=-9999
                    SIZE_SBTR=0
                 ENDIF
              ENDIF
           ENDIF
           GOTO 789
        ENDIF
        GOTO 999
 789    CONTINUE
        NBROOT=NA(2)
        NBLEAF=NA(1)
        PEAK=0.0E0
        FACT_SIZE=0
        DO I=1,NBROOT
           PEAK=MAX(PEAK,REAL(M(STEP(NA(2+NBLEAF+I)))))
           FACT_SIZE=FACT_SIZE+FACT(STEP(NA(2+NBLEAF+I)))
        ENDDO
 5483   CONTINUE
        DEALLOCATE(IPOOL)
        DEALLOCATE(M)
        DEALLOCATE(fact)
        DEALLOCATE(TNSTK)
        DEALLOCATE(SON)
        DEALLOCATE(TAB2)
        DEALLOCATE(TAB1)
        DEALLOCATE(T1)
        DEALLOCATE(T2)
        DEALLOCATE(RESULT)
        DEALLOCATE(TEMP)
        IF((PERM.EQ.3).OR.(PERM.EQ.4))THEN
           DEALLOCATE(DEPTH)
        ENDIF
        IF (SBTR_M.OR.(PERM.EQ.2))  THEN
           IF((SBTR_WHICH_M.EQ.1).AND.(PERM.NE.1).OR.(PERM.EQ.2))THEN
              DEALLOCATE(M_TOTAL)
           ENDIF
        ENDIF
      RETURN
      END SUBROUTINE SMUMPS_363
      RECURSIVE SUBROUTINE SMUMPS_462(TAB,DIM,TAB1,TAB2,PERM,
     &  RESULT,TEMP1,TEMP2)
      IMPLICIT NONE
      INTEGER DIM
      INTEGER*8 TAB1(DIM),TAB2(DIM)
      INTEGER*8 TEMP1(DIM),TEMP2(DIM)
      INTEGER TAB(DIM), PERM,RESULT(DIM)
      INTEGER I,J,K,C,I1,I2,COR1,COR2
      IF(DIM.EQ.1) THEN
        RESULT(1)=TAB(1)
        TEMP1(1)=TAB1(1)
        TEMP2(1)=TAB2(1)
        RETURN
      ENDIF
      I=DIM/2
      CALL SMUMPS_462(TAB(1),I,TAB1(1),TAB2(1),PERM,
     &  RESULT(1),TEMP1(1),TEMP2(1))
      CALL SMUMPS_462(TAB(I+1),DIM-I,TAB1(I+1),TAB2(I+1),
     &  PERM,RESULT(I+1),TEMP1(I+1),TEMP2(I+1))
      I1=1
      I2=I+1
      J=1
      DO WHILE ((I1.LE.I).AND.(I2.LE.DIM))
        IF((PERM.EQ.3))THEN
          IF(TEMP1(I1).LE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.4).OR.(PERM.EQ.5))THEN
          IF (TEMP1(I1).GE.TEMP1(I2))THEN
            TAB(J)=RESULT(I1)
            TAB1(J)=TEMP1(I1)
            J=J+1
            I1=I1+1
          ELSE
            TAB(J)=RESULT(I2)
            TAB1(J)=TEMP1(I2)
            J=J+1
            I2=I2+1          
          ENDIF
          GOTO 3
        ENDIF
        IF((PERM.EQ.0).OR.(PERM.EQ.1).OR.(PERM.EQ.2)) THEN
          IF(TEMP1(I1).GT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
            GOTO 3
          ENDIF
          IF(TEMP1(I1).LT.TEMP1(I2))THEN
            TAB1(J)=TEMP1(I2)
            TAB2(J)=TEMP2(I2)
            TAB(J)=RESULT(I2)
            J=J+1
            I2=I2+1
            GOTO 3
          ENDIF        
          IF((TEMP1(I1).EQ.TEMP1(I2)))THEN
            IF(TEMP2(I1).LE.TEMP2(I2))THEN
              TAB1(J)=TEMP1(I1)
              TAB2(J)=TEMP2(I1)
              TAB(J)=RESULT(I1)
              J=J+1
              I1=I1+1
            ELSE
              TAB1(J)=TEMP1(I2)
              TAB2(J)=TEMP2(I2)
              TAB(J)=RESULT(I2)
              J=J+1
              I2=I2+1
            ENDIF
          ENDIF
        ENDIF
  3   CONTINUE    
      ENDDO
      IF(I1.GT.I)THEN
        DO WHILE(I2.LE.DIM)
          TAB(J)=RESULT(I2)
          TAB1(J)=TEMP1(I2)
          TAB2(J)=TEMP2(I2)
          J=J+1
          I2=I2+1
        ENDDO
      ELSE
        IF(I2.GT.DIM)THEN
          DO WHILE(I1.LE.I)
            TAB1(J)=TEMP1(I1)
            TAB2(J)=TEMP2(I1)
            TAB(J)=RESULT(I1)
            J=J+1
            I1=I1+1
          ENDDO
        ENDIF
      ENDIF
      DO I=1,DIM
        TEMP1(I)=TAB1(I)
        TEMP2(I)=TAB2(I)
        RESULT(I)=TAB(I)
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_462
