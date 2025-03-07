/*

   THIS FILE IS PART OF MUMPS VERSION 4.7.3
   This Version was built on Fri May  4 15:54:01 2007


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://mumps.enseeiht.fr/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent,
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet, Hybrid scheduling for the parallel solution of linear
   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).

*/
/* $Id: smumps_c.h,v 1.35 2007/03/25 19:43:31 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined SMUMPS_C_H )
#define SMUMPS_C_H

/* Complex datatypes */

/* Next line defines SMUMPS_INT, SMUMPS_DOUBLE and SMUMPS_DOUBLE2 */
#include "smumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    SMUMPS_INT sym, par, job;
    SMUMPS_INT comm_fortran;    /* Fortran communicator */
    SMUMPS_INT icntl[40];
    SMUMPS_DOUBLE2 cntl[15];
    SMUMPS_INT n;
   
    SMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    SMUMPS_INT nz; SMUMPS_INT *irn; SMUMPS_INT *jcn; SMUMPS_DOUBLE *a;
    /* Distributed entry */
    SMUMPS_INT nz_loc; SMUMPS_INT *irn_loc; SMUMPS_INT *jcn_loc; SMUMPS_DOUBLE *a_loc;
    /* Element entry */
    SMUMPS_INT nelt; SMUMPS_INT *eltptr; SMUMPS_INT *eltvar; SMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    SMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    SMUMPS_INT *sym_perm;
    /* column permutation */
    SMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    SMUMPS_DOUBLE *colsca; SMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    SMUMPS_DOUBLE *rhs, *redrhs, *rhs_sparse, *sol_loc;
    SMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    SMUMPS_INT nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
  SMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  SMUMPS_INT mblock, nblock, nprow, npcol;
    SMUMPS_INT info[40],infog[40];
    SMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    SMUMPS_INT deficiency; SMUMPS_INT * pivnul_list; SMUMPS_INT * mapping;
    /* Schur */
    SMUMPS_INT size_schur; SMUMPS_INT *listvar_schur; SMUMPS_DOUBLE *schur;
    /* Internal parameters */
    SMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
    char version_number[80];
  } SMUMPS_STRUC_C;


#if defined(UPPER) || defined(_WIN32)
#define smumps_f77_ SMUMPS_F77
#define smumps_affect_mapping_ SMUMPS_AFFECT_MAPPING
#define smumps_affect_pivnul_list_ SMUMPS_AFFECT_PIVNUL_LIST
#define smumps_affect_colsca_ SMUMPS_AFFECT_COLSCA
#define smumps_affect_rowsca_ SMUMPS_AFFECT_ROWSCA 
#define smumps_affect_uns_perm_     SMUMPS_AFFECT_UNS_PERM
#define smumps_affect_sym_perm_     SMUMPS_AFFECT_SYM_PERM
#define smumps_nullify_c_mapping_   SMUMPS_NULLIFY_C_MAPPING
#define smumps_nullify_c_pivnul_list_ SMUMPS_NULLIFY_C_PIVNUL_LIST
#define smumps_nullify_c_sym_perm_  SMUMPS_NULLIFY_C_SYM_PERM
#define smumps_nullify_c_uns_perm_  SMUMPS_NULLIFY_C_UNS_PERM
#define smumps_nullify_c_colsca_    SMUMPS_NULLIFY_C_COLSCA
#define smumps_nullify_c_rowsca_    SMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
#define smumps_f77_ smumps_f77__
#define smumps_affect_mapping_ smumps_affect_mapping__
#define smumps_affect_pivnul_list_ smumps_affect_pivnul_list__
#define smumps_affect_colsca_ smumps_affect_colsca__
#define smumps_affect_rowsca_ smumps_affect_rowsca__
#define smumps_affect_uns_perm_     smumps_affect_uns_perm__     
#define smumps_affect_sym_perm_     smumps_affect_sym_perm__     
#define smumps_nullify_c_mapping_   smumps_nullify_c_mapping__    
#define smumps_nullify_c_pivnul_list_ smumps_nullify_c_pivnul_list__  
#define smumps_nullify_c_sym_perm_  smumps_nullify_c_sym_perm__   
#define smumps_nullify_c_uns_perm_  smumps_nullify_c_uns_perm__   
#define smumps_nullify_c_colsca_    smumps_nullify_c_colsca__     
#define smumps_nullify_c_rowsca_    smumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
#define smumps_f77_ smumps_f77
#define smumps_affect_mapping_ smumps_affect_mapping
#define smumps_affect_pivnul_list_ smumps_affect_pivnul_list
#define smumps_affect_colsca_ smumps_affect_colsca
#define smumps_affect_rowsca_ smumps_affect_rowsca
#define smumps_affect_uns_perm_     smumps_affect_uns_perm     
#define smumps_affect_sym_perm_     smumps_affect_sym_perm     
#define smumps_nullify_c_mapping_   smumps_nullify_c_mapping    
#define smumps_nullify_c_pivnul_list_ smumps_nullify_c_pivnul_list
#define smumps_nullify_c_sym_perm_  smumps_nullify_c_sym_perm   
#define smumps_nullify_c_uns_perm_  smumps_nullify_c_uns_perm   
#define smumps_nullify_c_colsca_    smumps_nullify_c_colsca     
#define smumps_nullify_c_rowsca_    smumps_nullify_c_rowsca
#endif

#ifndef MUMPS_CALL
#if defined(_WIN32)
/* Modify/choose between next 2 lines depending
 * on your Windows calling conventions */
/* #define MUMPS_CALL __stdcall */
#define MUMPS_CALL
#else
#define MUMPS_CALL
#endif
#endif

void MUMPS_CALL smumps_c(SMUMPS_STRUC_C * smumps_par);
void MUMPS_CALL smumps_affect_mapping_(SMUMPS_INT * f77mapping);
void MUMPS_CALL smumps_affect_pivnul_list_(SMUMPS_INT * f77pivnul_list);
void MUMPS_CALL smumps_affect_uns_perm_(SMUMPS_INT * f77sym_perm);
void MUMPS_CALL smumps_affect_sym_perm_(SMUMPS_INT * f77uns_perm);
void MUMPS_CALL smumps_nullify_c_mapping_();
void MUMPS_CALL smumps_nullify_c_pivnul_list_();
void MUMPS_CALL smumps_nullify_c_sym_perm_();
void MUMPS_CALL smumps_nullify_c_uns_perm_();
#ifdef return_scaling
void MUMPS_CALL smumps_affect_colsca_(SMUMPS_DOUBLE * f77colsca);
void MUMPS_CALL smumps_affect_rowsca_(SMUMPS_DOUBLE * f77rowsca);
void MUMPS_CALL smumps_nullify_c_colsca_();
void MUMPS_CALL smumps_nullify_c_rowsca_();
#endif

void MUMPS_CALL smumps_f77_(SMUMPS_INT *job, SMUMPS_INT *sym, SMUMPS_INT *par, SMUMPS_INT * comm_fortran,
SMUMPS_INT *n, SMUMPS_INT *icntl,
SMUMPS_DOUBLE2 *cntl, SMUMPS_INT *nz, SMUMPS_INT *irn, SMUMPS_INT *irn_avail,
SMUMPS_INT *jcn, SMUMPS_INT *jcn_avail, SMUMPS_DOUBLE *a, SMUMPS_INT *a_avail,
  SMUMPS_INT *nz_loc, SMUMPS_INT * irn_loc, SMUMPS_INT * irn_loc_avail,
  SMUMPS_INT * jcn_loc, SMUMPS_INT * jcn_loc_avail, SMUMPS_DOUBLE *a_loc,
  SMUMPS_INT * a_loc_avail, SMUMPS_INT * nelt, SMUMPS_INT * eltptr,
  SMUMPS_INT * eltptr_avail, SMUMPS_INT * eltvar, SMUMPS_INT * eltvar_avail,
  SMUMPS_DOUBLE * a_elt, SMUMPS_INT * a_elt_avail, SMUMPS_INT * perm_in,
  SMUMPS_INT * perm_in_avail, SMUMPS_DOUBLE * rhs, SMUMPS_INT * rhs_avail,
  SMUMPS_DOUBLE * redrhs, SMUMPS_INT * redrhs_avail,
  SMUMPS_INT * info, SMUMPS_DOUBLE2 * rinfo, SMUMPS_INT * infog, SMUMPS_DOUBLE2 * rinfog,
  SMUMPS_INT * deficiency, SMUMPS_INT * size_schur, SMUMPS_INT * listvar_schur,
  SMUMPS_INT * listvar_schur_avail, SMUMPS_DOUBLE * schur,
  SMUMPS_INT * schur_avail, SMUMPS_DOUBLE * colsca, SMUMPS_INT * colsca_avail,
  SMUMPS_DOUBLE * rowsca, SMUMPS_INT * rowsca_avail, SMUMPS_INT * instance_number,
  SMUMPS_INT * nrhs, SMUMPS_INT * lrhs, SMUMPS_INT * lredrhs, SMUMPS_DOUBLE * rhs_sparse, SMUMPS_INT * rhs_sparse_avail,
  SMUMPS_DOUBLE * sol_loc, SMUMPS_INT * sol_loc_avail, SMUMPS_INT * irhs_sparse,
  SMUMPS_INT * irhs_sparse_avail, SMUMPS_INT * irhs_ptr, SMUMPS_INT * irhs_ptr_avail,
  SMUMPS_INT * isol_loc, SMUMPS_INT * isol_loc_avail, SMUMPS_INT * nz_rhs, SMUMPS_INT * lsol_loc
  ,SMUMPS_INT * schur_mloc, SMUMPS_INT *schur_nloc, SMUMPS_INT * schur_lld,
  SMUMPS_INT * schur_mblock, SMUMPS_INT * schur_nblock, SMUMPS_INT * schur_nprow,
  SMUMPS_INT * schur_npcol,
  SMUMPS_INT * ooc_tmpdir,
  SMUMPS_INT * ooc_prefix,
  SMUMPS_INT * ooc_tmpdirlen,
  SMUMPS_INT * ooc_prefixlen
  );
#endif

