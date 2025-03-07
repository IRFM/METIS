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
/* $Id: dmumps_c.h,v 1.35 2007/03/25 19:43:31 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined DMUMPS_C_H )
#define DMUMPS_C_H

/* Complex datatypes */

/* Next line defines DMUMPS_INT, DMUMPS_DOUBLE and DMUMPS_DOUBLE2 */
#include "dmumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    DMUMPS_INT sym, par, job;
    DMUMPS_INT comm_fortran;    /* Fortran communicator */
    DMUMPS_INT icntl[40];
    DMUMPS_DOUBLE2 cntl[15];
    DMUMPS_INT n;
   
    DMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    DMUMPS_INT nz; DMUMPS_INT *irn; DMUMPS_INT *jcn; DMUMPS_DOUBLE *a;
    /* Distributed entry */
    DMUMPS_INT nz_loc; DMUMPS_INT *irn_loc; DMUMPS_INT *jcn_loc; DMUMPS_DOUBLE *a_loc;
    /* Element entry */
    DMUMPS_INT nelt; DMUMPS_INT *eltptr; DMUMPS_INT *eltvar; DMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    DMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    DMUMPS_INT *sym_perm;
    /* column permutation */
    DMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    DMUMPS_DOUBLE *colsca; DMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    DMUMPS_DOUBLE *rhs, *redrhs, *rhs_sparse, *sol_loc;
    DMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    DMUMPS_INT nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
  DMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  DMUMPS_INT mblock, nblock, nprow, npcol;
    DMUMPS_INT info[40],infog[40];
    DMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    DMUMPS_INT deficiency; DMUMPS_INT * pivnul_list; DMUMPS_INT * mapping;
    /* Schur */
    DMUMPS_INT size_schur; DMUMPS_INT *listvar_schur; DMUMPS_DOUBLE *schur;
    /* Internal parameters */
    DMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
    char version_number[80];
  } DMUMPS_STRUC_C;


#if defined(UPPER) || defined(_WIN32)
#define dmumps_f77_ DMUMPS_F77
#define dmumps_affect_mapping_ DMUMPS_AFFECT_MAPPING
#define dmumps_affect_pivnul_list_ DMUMPS_AFFECT_PIVNUL_LIST
#define dmumps_affect_colsca_ DMUMPS_AFFECT_COLSCA
#define dmumps_affect_rowsca_ DMUMPS_AFFECT_ROWSCA 
#define dmumps_affect_uns_perm_     DMUMPS_AFFECT_UNS_PERM
#define dmumps_affect_sym_perm_     DMUMPS_AFFECT_SYM_PERM
#define dmumps_nullify_c_mapping_   DMUMPS_NULLIFY_C_MAPPING
#define dmumps_nullify_c_pivnul_list_ DMUMPS_NULLIFY_C_PIVNUL_LIST
#define dmumps_nullify_c_sym_perm_  DMUMPS_NULLIFY_C_SYM_PERM
#define dmumps_nullify_c_uns_perm_  DMUMPS_NULLIFY_C_UNS_PERM
#define dmumps_nullify_c_colsca_    DMUMPS_NULLIFY_C_COLSCA
#define dmumps_nullify_c_rowsca_    DMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
#define dmumps_f77_ dmumps_f77__
#define dmumps_affect_mapping_ dmumps_affect_mapping__
#define dmumps_affect_pivnul_list_ dmumps_affect_pivnul_list__
#define dmumps_affect_colsca_ dmumps_affect_colsca__
#define dmumps_affect_rowsca_ dmumps_affect_rowsca__
#define dmumps_affect_uns_perm_     dmumps_affect_uns_perm__     
#define dmumps_affect_sym_perm_     dmumps_affect_sym_perm__     
#define dmumps_nullify_c_mapping_   dmumps_nullify_c_mapping__    
#define dmumps_nullify_c_pivnul_list_ dmumps_nullify_c_pivnul_list__  
#define dmumps_nullify_c_sym_perm_  dmumps_nullify_c_sym_perm__   
#define dmumps_nullify_c_uns_perm_  dmumps_nullify_c_uns_perm__   
#define dmumps_nullify_c_colsca_    dmumps_nullify_c_colsca__     
#define dmumps_nullify_c_rowsca_    dmumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
#define dmumps_f77_ dmumps_f77
#define dmumps_affect_mapping_ dmumps_affect_mapping
#define dmumps_affect_pivnul_list_ dmumps_affect_pivnul_list
#define dmumps_affect_colsca_ dmumps_affect_colsca
#define dmumps_affect_rowsca_ dmumps_affect_rowsca
#define dmumps_affect_uns_perm_     dmumps_affect_uns_perm     
#define dmumps_affect_sym_perm_     dmumps_affect_sym_perm     
#define dmumps_nullify_c_mapping_   dmumps_nullify_c_mapping    
#define dmumps_nullify_c_pivnul_list_ dmumps_nullify_c_pivnul_list
#define dmumps_nullify_c_sym_perm_  dmumps_nullify_c_sym_perm   
#define dmumps_nullify_c_uns_perm_  dmumps_nullify_c_uns_perm   
#define dmumps_nullify_c_colsca_    dmumps_nullify_c_colsca     
#define dmumps_nullify_c_rowsca_    dmumps_nullify_c_rowsca
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

void MUMPS_CALL dmumps_c(DMUMPS_STRUC_C * dmumps_par);
void MUMPS_CALL dmumps_affect_mapping_(DMUMPS_INT * f77mapping);
void MUMPS_CALL dmumps_affect_pivnul_list_(DMUMPS_INT * f77pivnul_list);
void MUMPS_CALL dmumps_affect_uns_perm_(DMUMPS_INT * f77sym_perm);
void MUMPS_CALL dmumps_affect_sym_perm_(DMUMPS_INT * f77uns_perm);
void MUMPS_CALL dmumps_nullify_c_mapping_();
void MUMPS_CALL dmumps_nullify_c_pivnul_list_();
void MUMPS_CALL dmumps_nullify_c_sym_perm_();
void MUMPS_CALL dmumps_nullify_c_uns_perm_();
#ifdef return_scaling
void MUMPS_CALL dmumps_affect_colsca_(DMUMPS_DOUBLE * f77colsca);
void MUMPS_CALL dmumps_affect_rowsca_(DMUMPS_DOUBLE * f77rowsca);
void MUMPS_CALL dmumps_nullify_c_colsca_();
void MUMPS_CALL dmumps_nullify_c_rowsca_();
#endif

void MUMPS_CALL dmumps_f77_(DMUMPS_INT *job, DMUMPS_INT *sym, DMUMPS_INT *par, DMUMPS_INT * comm_fortran,
DMUMPS_INT *n, DMUMPS_INT *icntl,
DMUMPS_DOUBLE2 *cntl, DMUMPS_INT *nz, DMUMPS_INT *irn, DMUMPS_INT *irn_avail,
DMUMPS_INT *jcn, DMUMPS_INT *jcn_avail, DMUMPS_DOUBLE *a, DMUMPS_INT *a_avail,
  DMUMPS_INT *nz_loc, DMUMPS_INT * irn_loc, DMUMPS_INT * irn_loc_avail,
  DMUMPS_INT * jcn_loc, DMUMPS_INT * jcn_loc_avail, DMUMPS_DOUBLE *a_loc,
  DMUMPS_INT * a_loc_avail, DMUMPS_INT * nelt, DMUMPS_INT * eltptr,
  DMUMPS_INT * eltptr_avail, DMUMPS_INT * eltvar, DMUMPS_INT * eltvar_avail,
  DMUMPS_DOUBLE * a_elt, DMUMPS_INT * a_elt_avail, DMUMPS_INT * perm_in,
  DMUMPS_INT * perm_in_avail, DMUMPS_DOUBLE * rhs, DMUMPS_INT * rhs_avail,
  DMUMPS_DOUBLE * redrhs, DMUMPS_INT * redrhs_avail,
  DMUMPS_INT * info, DMUMPS_DOUBLE2 * rinfo, DMUMPS_INT * infog, DMUMPS_DOUBLE2 * rinfog,
  DMUMPS_INT * deficiency, DMUMPS_INT * size_schur, DMUMPS_INT * listvar_schur,
  DMUMPS_INT * listvar_schur_avail, DMUMPS_DOUBLE * schur,
  DMUMPS_INT * schur_avail, DMUMPS_DOUBLE * colsca, DMUMPS_INT * colsca_avail,
  DMUMPS_DOUBLE * rowsca, DMUMPS_INT * rowsca_avail, DMUMPS_INT * instance_number,
  DMUMPS_INT * nrhs, DMUMPS_INT * lrhs, DMUMPS_INT * lredrhs, DMUMPS_DOUBLE * rhs_sparse, DMUMPS_INT * rhs_sparse_avail,
  DMUMPS_DOUBLE * sol_loc, DMUMPS_INT * sol_loc_avail, DMUMPS_INT * irhs_sparse,
  DMUMPS_INT * irhs_sparse_avail, DMUMPS_INT * irhs_ptr, DMUMPS_INT * irhs_ptr_avail,
  DMUMPS_INT * isol_loc, DMUMPS_INT * isol_loc_avail, DMUMPS_INT * nz_rhs, DMUMPS_INT * lsol_loc
  ,DMUMPS_INT * schur_mloc, DMUMPS_INT *schur_nloc, DMUMPS_INT * schur_lld,
  DMUMPS_INT * schur_mblock, DMUMPS_INT * schur_nblock, DMUMPS_INT * schur_nprow,
  DMUMPS_INT * schur_npcol,
  DMUMPS_INT * ooc_tmpdir,
  DMUMPS_INT * ooc_prefix,
  DMUMPS_INT * ooc_tmpdirlen,
  DMUMPS_INT * ooc_prefixlen
  );
#endif

