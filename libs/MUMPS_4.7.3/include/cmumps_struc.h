!
!   THIS FILE IS PART OF CMUMPS VERSION 4.7.3
!   This Version was built on Fri May  4 15:54:01 2007
!
!
!  This version of CMUMPS is provided to you free of charge. It is public
!  domain, based on public domain software developed during the Esprit IV
!  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
!  Since this first public domain version in 1999, the developments are
!  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
!  INRIA.
!
!  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
!  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.
!
!  Up-to-date copies of the CMUMPS package can be obtained
!  from the Web pages http://mumps.enseeiht.fr/
!  or http://graal.ens-lyon.fr/CMUMPS
!
!
!   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
!   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
!
!
!  User documentation of any code that uses this software can
!  include this complete notice. You can acknowledge (using
!  references [1], [2], and [3] the contribution of this package
!  in any scientific publication dependent upon the use of the
!  package. You shall use reasonable endeavours to notify
!  the authors of the package of this publication.
!
!   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
!   Multifrontal parallel distributed symmetric and unsymmetric solvers,
!   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).
!
!   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
!   A fully asynchronous multifrontal solver using distributed dynamic
!   scheduling, SIAM Journal of Matrix Analysis and Applications,
!   Vol 23, No 1, pp 15-41 (2001).
!
!   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
!   S. Pralet, Hybrid scheduling for the parallel solution of linear
!   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
!
!     $Id: cmumps_struc.h,v 1.72 2007/03/25 19:43:31 jylexcel Exp $
      INCLUDE 'cmumps_root.h'
      TYPE CMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters 
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite, 
!         SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
         INTEGER SYM, PAR
         INTEGER JOB 
!    --------------------
!    Order of Input matrix 
!    --------------------
         INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
         INTEGER NZ
         COMPLEX, DIMENSION(:), POINTER :: A
         INTEGER, DIMENSION(:), POINTER :: IRN, JCN
         COMPLEX, DIMENSION(:), POINTER :: COLSCA, ROWSCA
!
!        ------------------------------------
!        Case of distributed assembled matrix
!        matrix on entry:
!        ------------------------------------
         INTEGER NZ_loc
         INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
         COMPLEX, DIMENSION(:), POINTER :: A_loc
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
         INTEGER NELT
         INTEGER, DIMENSION(:), POINTER :: ELTPTR
         INTEGER, DIMENSION(:), POINTER :: ELTVAR
         COMPLEX, DIMENSION(:), POINTER :: A_ELT
!
!    -----------------
!    MPI Communicator
!    -----------------
         INTEGER COMM
!
!    ---------------------------------------------
!    Symmetric permutation : 
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data 
! ******************
!    --------------------------------------------------------
!    RHS / SOL_LOC
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
         COMPLEX, DIMENSION(:), POINTER :: RHS, REDRHS
         COMPLEX, DIMENSION(:), POINTER :: RHS_SPARSE
         COMPLEX, DIMENSION(:), POINTER :: SOL_LOC
         INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
         INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
         INTEGER, DIMENSION(:), POINTER :: ISOL_LOC
         INTEGER LRHS, NRHS, NZ_RHS, LSOL_LOC, LREDRHS
         INTEGER PADDING
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
         INTEGER ICNTL(40)
         INTEGER INFO(40) 
         INTEGER INFOG(40)
         REAL CNTL(15)
         REAL RINFO(20)
         REAL RINFOG(20)
!        Cost (flops) of subtrees on local process
         REAL COST_SUBTREES
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!        SYM_PERM: Symmetric permutation 
!        UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
         INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
! 
!    -------------------------------------
!    Case of distributed matrix on entry:
!    CMUMPS potentially provides mapping
!    -------------------------------------
         INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -------------------------------
!    Deficiency and null space basis
!    -------------------------------
         INTEGER Deficiency
         COMPLEX, DIMENSION(:,:), POINTER :: NULL_SPACE
!    -----
!    Schur
!    -----
         INTEGER NPROW, NPCOL, MBLOCK, NBLOCK
         INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
         INTEGER SIZE_SCHUR
         INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
         COMPLEX, DIMENSION(:), POINTER :: SCHUR
         COMPLEX, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
!    -----------
!    Out-of-core
!    -----------
         CHARACTER(LEN=150) :: OOC_TMPDIR, OOC_PREFIX
!
!
! **********************
! INTERNAL Working data
! **********************
!        For MPI
         INTEGER COMM_NODES, MYID_NODES, COMM_LOAD
         INTEGER  MYID, NPROCS, NSLAVES
         INTEGER ASS_IRECV
         INTEGER, DIMENSION(:), POINTER :: POIDS
         INTEGER LBUFR
         INTEGER LBUFR_BYTES
         INTEGER, DIMENSION(:), POINTER ::  BUFR
         INTEGER INST_Number
!        for analysis/facto/solve phases
         INTEGER MAXIS, MAXS
         INTEGER MAXIS1
         INTEGER KEEP(500)
         INTEGER*8 KEEP8(150)
!        IS is used for the factors + workspace for contrib. blocks
         INTEGER, DIMENSION(:), POINTER :: IS
!        is1 (maxis1) contains working arrays computed 
!        and used only during analysis
         INTEGER, DIMENSION(:), POINTER :: IS1
!        The following data/arrays are computed during the analysis
!        phase and used during the factorization and solve phases.
         INTEGER LNA
         INTEGER,POINTER,DIMENSION(:)::STEP, NE_STEPS, ND_STEPS
         INTEGER,POINTER,DIMENSION(:)::FRERE_STEPS, DAD_STEPS
         INTEGER,POINTER,DIMENSION(:)::FILS, PTRAR, FRTPTR, FRTELT
         INTEGER,POINTER,DIMENSION(:)::NA, PROCNODE_STEPS
!        The two pointer arrays computed in facto and used by the solve
!           (except the factors) are PTLUST_S and PTRFAC. 
         INTEGER, DIMENSION(:), POINTER :: PTLUST_S, PTRFAC
!        main real working arrays for factorization/solve phases
         COMPLEX, DIMENSION(:), POINTER :: S
!        Information on mapping
         INTEGER, DIMENSION(:), POINTER :: PROCNODE
         INTEGER nbsa
!        Input matrix ready for numerical assembly 
!            -arrowhead format in case of assembled matrix
!            -element format otherwise
         INTEGER, DIMENSION(:), POINTER :: INTARR
         COMPLEX, DIMENSION(:), POINTER :: DBLARR
!         Element entry: internal data
         INTEGER NELT_LOC, LELTVAR, NA_ELT
         INTEGER, DIMENSION(:), POINTER :: ELTPROC
!         Candidates and node partitionning
         INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
         INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
         INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
         INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE 
	 LOGICAL, DIMENSION(:), POINTER :: I_AM_CAND
!        For heterogeneous architecture
         INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!        Compressed RHS
         INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP
         COMPLEX, DIMENSION(:), POINTER :: RHSCOMP
!        To save the matrix in a simple format
         CHARACTER(LEN=80) WRITE_PROBLEM
!   ------------------------
!   Root structure(internal)
!   ------------------------
         TYPE (CMUMPS_ROOT_STRUC) :: root
!        For C interface
!   Instance number used/managed by the C/F77 interface
         INTEGER INSTANCE_NUMBER
!   Info on the subtrees to be used during facto
         INTEGER, DIMENSION(:),   POINTER :: MEM_SUBTREE
         INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
         INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
         INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
         INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
         DOUBLE PRECISION, DIMENSION(:),   POINTER :: COST_TRAV
         INTEGER NBSA_LOCAL
         INTEGER MAX_SURF_MASTER
!    For simulating parallel out-of-core stack.
         INTEGER, DIMENSION(:),POINTER ::CB_SON_SIZE
!    OOC management data that must persist from factorization to solve.
         INTEGER, DIMENSION(:),   POINTER :: OOC_INODE_SEQUENCE
         INTEGER, DIMENSION(:),   POINTER :: OOC_NUM_FILE
         INTEGER, DIMENSION(:),   POINTER :: OOC_POS_IN_FILE
         INTEGER, DIMENSION(:),   POINTER :: OOC_SIZE_OF_BLOCK
         INTEGER OOC_MAX_NB_NODES_FOR_ZONE,OOC_TOTAL_NB_NODES
         INTEGER OOC_NB_FILES
         CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES  
         INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
         INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Internal control array
         REAL DKEEP(30)
!    Array needed to manage additionnal candidate processor 
!    when using -DNEW_MAPPING
         INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC
!    Version number
         CHARACTER(LEN=80) VERSION_NUMBER
      END TYPE CMUMPS_STRUC
