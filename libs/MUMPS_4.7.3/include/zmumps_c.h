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
/* $Id: zmumps_c.h,v 1.35 2007/03/25 19:43:31 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined ZMUMPS_C_H )
#define ZMUMPS_C_H

/* Complex datatypes */
typedef struct {double r,i;} mumps_double_complex;

/* Next line defines ZMUMPS_INT, ZMUMPS_DOUBLE and ZMUMPS_DOUBLE2 */
#include "zmumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    ZMUMPS_INT sym, par, job;
    ZMUMPS_INT comm_fortran;    /* Fortran communicator */
    ZMUMPS_INT icntl[40];
    ZMUMPS_DOUBLE2 cntl[15];
    ZMUMPS_INT n;
   
    ZMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    ZMUMPS_INT nz; ZMUMPS_INT *irn; ZMUMPS_INT *jcn; ZMUMPS_DOUBLE *a;
    /* Distributed entry */
    ZMUMPS_INT nz_loc; ZMUMPS_INT *irn_loc; ZMUMPS_INT *jcn_loc; ZMUMPS_DOUBLE *a_loc;
    /* Element entry */
    ZMUMPS_INT nelt; ZMUMPS_INT *eltptr; ZMUMPS_INT *eltvar; ZMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    ZMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    ZMUMPS_INT *sym_perm;
    /* column permutation */
    ZMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    ZMUMPS_DOUBLE *colsca; ZMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    ZMUMPS_DOUBLE *rhs, *redrhs, *rhs_sparse, *sol_loc;
    ZMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    ZMUMPS_INT nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
  ZMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  ZMUMPS_INT mblock, nblock, nprow, npcol;
    ZMUMPS_INT info[40],infog[40];
    ZMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    ZMUMPS_INT deficiency; ZMUMPS_INT * pivnul_list; ZMUMPS_INT * mapping;
    /* Schur */
    ZMUMPS_INT size_schur; ZMUMPS_INT *listvar_schur; ZMUMPS_DOUBLE *schur;
    /* Internal parameters */
    ZMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
    char version_number[80];
  } ZMUMPS_STRUC_C;


#if defined(UPPER) || defined(_WIN32)
#define zmumps_f77_ ZMUMPS_F77
#define zmumps_affect_mapping_ ZMUMPS_AFFECT_MAPPING
#define zmumps_affect_pivnul_list_ ZMUMPS_AFFECT_PIVNUL_LIST
#define zmumps_affect_colsca_ ZMUMPS_AFFECT_COLSCA
#define zmumps_affect_rowsca_ ZMUMPS_AFFECT_ROWSCA 
#define zmumps_affect_uns_perm_     ZMUMPS_AFFECT_UNS_PERM
#define zmumps_affect_sym_perm_     ZMUMPS_AFFECT_SYM_PERM
#define zmumps_nullify_c_mapping_   ZMUMPS_NULLIFY_C_MAPPING
#define zmumps_nullify_c_pivnul_list_ ZMUMPS_NULLIFY_C_PIVNUL_LIST
#define zmumps_nullify_c_sym_perm_  ZMUMPS_NULLIFY_C_SYM_PERM
#define zmumps_nullify_c_uns_perm_  ZMUMPS_NULLIFY_C_UNS_PERM
#define zmumps_nullify_c_colsca_    ZMUMPS_NULLIFY_C_COLSCA
#define zmumps_nullify_c_rowsca_    ZMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
#define zmumps_f77_ zmumps_f77__
#define zmumps_affect_mapping_ zmumps_affect_mapping__
#define zmumps_affect_pivnul_list_ zmumps_affect_pivnul_list__
#define zmumps_affect_colsca_ zmumps_affect_colsca__
#define zmumps_affect_rowsca_ zmumps_affect_rowsca__
#define zmumps_affect_uns_perm_     zmumps_affect_uns_perm__     
#define zmumps_affect_sym_perm_     zmumps_affect_sym_perm__     
#define zmumps_nullify_c_mapping_   zmumps_nullify_c_mapping__    
#define zmumps_nullify_c_pivnul_list_ zmumps_nullify_c_pivnul_list__  
#define zmumps_nullify_c_sym_perm_  zmumps_nullify_c_sym_perm__   
#define zmumps_nullify_c_uns_perm_  zmumps_nullify_c_uns_perm__   
#define zmumps_nullify_c_colsca_    zmumps_nullify_c_colsca__     
#define zmumps_nullify_c_rowsca_    zmumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
#define zmumps_f77_ zmumps_f77
#define zmumps_affect_mapping_ zmumps_affect_mapping
#define zmumps_affect_pivnul_list_ zmumps_affect_pivnul_list
#define zmumps_affect_colsca_ zmumps_affect_colsca
#define zmumps_affect_rowsca_ zmumps_affect_rowsca
#define zmumps_affect_uns_perm_     zmumps_affect_uns_perm     
#define zmumps_affect_sym_perm_     zmumps_affect_sym_perm     
#define zmumps_nullify_c_mapping_   zmumps_nullify_c_mapping    
#define zmumps_nullify_c_pivnul_list_ zmumps_nullify_c_pivnul_list
#define zmumps_nullify_c_sym_perm_  zmumps_nullify_c_sym_perm   
#define zmumps_nullify_c_uns_perm_  zmumps_nullify_c_uns_perm   
#define zmumps_nullify_c_colsca_    zmumps_nullify_c_colsca     
#define zmumps_nullify_c_rowsca_    zmumps_nullify_c_rowsca
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

void MUMPS_CALL zmumps_c(ZMUMPS_STRUC_C * zmumps_par);
void MUMPS_CALL zmumps_affect_mapping_(ZMUMPS_INT * f77mapping);
void MUMPS_CALL zmumps_affect_pivnul_list_(ZMUMPS_INT * f77pivnul_list);
void MUMPS_CALL zmumps_affect_uns_perm_(ZMUMPS_INT * f77sym_perm);
void MUMPS_CALL zmumps_affect_sym_perm_(ZMUMPS_INT * f77uns_perm);
void MUMPS_CALL zmumps_nullify_c_mapping_();
void MUMPS_CALL zmumps_nullify_c_pivnul_list_();
void MUMPS_CALL zmumps_nullify_c_sym_perm_();
void MUMPS_CALL zmumps_nullify_c_uns_perm_();
#ifdef return_scaling
void MUMPS_CALL zmumps_affect_colsca_(ZMUMPS_DOUBLE * f77colsca);
void MUMPS_CALL zmumps_affect_rowsca_(ZMUMPS_DOUBLE * f77rowsca);
void MUMPS_CALL zmumps_nullify_c_colsca_();
void MUMPS_CALL zmumps_nullify_c_rowsca_();
#endif

void MUMPS_CALL zmumps_f77_(ZMUMPS_INT *job, ZMUMPS_INT *sym, ZMUMPS_INT *par, ZMUMPS_INT * comm_fortran,
ZMUMPS_INT *n, ZMUMPS_INT *icntl,
ZMUMPS_DOUBLE2 *cntl, ZMUMPS_INT *nz, ZMUMPS_INT *irn, ZMUMPS_INT *irn_avail,
ZMUMPS_INT *jcn, ZMUMPS_INT *jcn_avail, ZMUMPS_DOUBLE *a, ZMUMPS_INT *a_avail,
  ZMUMPS_INT *nz_loc, ZMUMPS_INT * irn_loc, ZMUMPS_INT * irn_loc_avail,
  ZMUMPS_INT * jcn_loc, ZMUMPS_INT * jcn_loc_avail, ZMUMPS_DOUBLE *a_loc,
  ZMUMPS_INT * a_loc_avail, ZMUMPS_INT * nelt, ZMUMPS_INT * eltptr,
  ZMUMPS_INT * eltptr_avail, ZMUMPS_INT * eltvar, ZMUMPS_INT * eltvar_avail,
  ZMUMPS_DOUBLE * a_elt, ZMUMPS_INT * a_elt_avail, ZMUMPS_INT * perm_in,
  ZMUMPS_INT * perm_in_avail, ZMUMPS_DOUBLE * rhs, ZMUMPS_INT * rhs_avail,
  ZMUMPS_DOUBLE * redrhs, ZMUMPS_INT * redrhs_avail,
  ZMUMPS_INT * info, ZMUMPS_DOUBLE2 * rinfo, ZMUMPS_INT * infog, ZMUMPS_DOUBLE2 * rinfog,
  ZMUMPS_INT * deficiency, ZMUMPS_INT * size_schur, ZMUMPS_INT * listvar_schur,
  ZMUMPS_INT * listvar_schur_avail, ZMUMPS_DOUBLE * schur,
  ZMUMPS_INT * schur_avail, ZMUMPS_DOUBLE * colsca, ZMUMPS_INT * colsca_avail,
  ZMUMPS_DOUBLE * rowsca, ZMUMPS_INT * rowsca_avail, ZMUMPS_INT * instance_number,
  ZMUMPS_INT * nrhs, ZMUMPS_INT * lrhs, ZMUMPS_INT * lredrhs, ZMUMPS_DOUBLE * rhs_sparse, ZMUMPS_INT * rhs_sparse_avail,
  ZMUMPS_DOUBLE * sol_loc, ZMUMPS_INT * sol_loc_avail, ZMUMPS_INT * irhs_sparse,
  ZMUMPS_INT * irhs_sparse_avail, ZMUMPS_INT * irhs_ptr, ZMUMPS_INT * irhs_ptr_avail,
  ZMUMPS_INT * isol_loc, ZMUMPS_INT * isol_loc_avail, ZMUMPS_INT * nz_rhs, ZMUMPS_INT * lsol_loc
  ,ZMUMPS_INT * schur_mloc, ZMUMPS_INT *schur_nloc, ZMUMPS_INT * schur_lld,
  ZMUMPS_INT * schur_mblock, ZMUMPS_INT * schur_nblock, ZMUMPS_INT * schur_nprow,
  ZMUMPS_INT * schur_npcol,
  ZMUMPS_INT * ooc_tmpdir,
  ZMUMPS_INT * ooc_prefix,
  ZMUMPS_INT * ooc_tmpdirlen,
  ZMUMPS_INT * ooc_prefixlen
  );
#endif

