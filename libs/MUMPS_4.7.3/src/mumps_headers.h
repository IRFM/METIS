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
C     Common header positions:
C
C     XXI    ->  size of integer record
C     XXR    ->  size of real record
C     XXS    ->  state of the node
C     XXN    ->  node number
C     XXP    ->  pointer to previous record
C
C  During Solve phase
C     XXS --> state of the master of a level 2 node:
C             -9999  =  undefined
C             C_FINI (no more contribution need be received)
C
      INTEGER XXI, XXR, XXS, XXN, XXP, XSIZE
      PARAMETER(XXI=0,XXR=1,XXS=2,XXN=3,XXP=4,XSIZE=5)
C     ------------------------
C     Compressed CB of level 1
C     ------------------------
      INTEGER S_CB1COMP
      PARAMETER (S_CB1COMP=314)
C     ----------------------
C     State of type 2 strips
C     ----------------------
      INTEGER S_ACTIVE, S_ALL, S_NOLCBCONTIG,
     *        S_NOLCBNOCONTIG, S_NOLCLEANED,
     *        S_NOLCBNOCONTIG38, S_NOLCBCONTIG38,
     *        S_NOLCLEANED38, C_FINI
      PARAMETER(S_ACTIVE=400, S_ALL=401, S_NOLCBCONTIG=402,
     *          S_NOLCBNOCONTIG=403, S_NOLCLEANED=404,
     *          S_NOLCBNOCONTIG38=405, S_NOLCBCONTIG38=406,
     *          S_NOLCLEANED38=407,C_FINI=1)
C     (note: S_NOLCBNOCONTIG38 corresponds to nodes
C            that only have NELIMxNBROW significant
C            entries and that can be compacted accordingly)
C
C     Free : the record can be compressed/removed
      INTEGER S_FREE, S_NOTFREE
      PARAMETER(S_FREE=54321,S_NOTFREE=-123456)
C
C
      INTEGER TOP_OF_STACK
      PARAMETER(TOP_OF_STACK=-999999)
C
C     -----------------------
C     XTRA_SLAVES_SYM: extra dummy slaves added
C     for strips in the symmetric case to hold:
C       - counter of BLFAC_SLAVE messages received
C       - NFRONT
C       - information about recursivity and reception
C         of ROOT2SON during recursivity
C         and strip cannot be freed yet.
      INTEGER XTRA_SLAVES_SYM, S_ROOT2SON_CALLED
      PARAMETER(XTRA_SLAVES_SYM=3, S_ROOT2SON_CALLED=-341)

C
C FIXME, later:
C   New state: S_ONDISK ? Maybe not.
C   IW(XXSIZEONDISK) could provide size written on disk.
C                   = 0 => in-core
C                   = IW(XXR) => all OOC ?
C         pb with type 2 bands, etc to be discussed.
C
