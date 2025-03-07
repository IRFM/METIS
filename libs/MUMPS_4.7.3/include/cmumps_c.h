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
/* $Id: cmumps_c.h,v 1.35 2007/03/25 19:43:31 jylexcel Exp $ */
/* Mostly written in march 2002 (JYL) */

#if ( ! defined CMUMPS_C_H )
#define CMUMPS_C_H

/* Complex datatypes */
typedef struct {float r,i;} mumps_complex;

/* Next line defines CMUMPS_INT, CMUMPS_DOUBLE and CMUMPS_DOUBLE2 */
#include "cmumps_prec.h"
/*
 * Definition of the (simplified)
 * MUMPS C structure
 */
typedef struct
  {
    CMUMPS_INT sym, par, job;
    CMUMPS_INT comm_fortran;    /* Fortran communicator */
    CMUMPS_INT icntl[40];
    CMUMPS_DOUBLE2 cntl[15];
    CMUMPS_INT n;
   
    CMUMPS_INT nz_alloc; /* used in matlab interface to decide if
                       we free + malloc when we have large variation */

    /* Assembled entry */
    CMUMPS_INT nz; CMUMPS_INT *irn; CMUMPS_INT *jcn; CMUMPS_DOUBLE *a;
    /* Distributed entry */
    CMUMPS_INT nz_loc; CMUMPS_INT *irn_loc; CMUMPS_INT *jcn_loc; CMUMPS_DOUBLE *a_loc;
    /* Element entry */
    CMUMPS_INT nelt; CMUMPS_INT *eltptr; CMUMPS_INT *eltvar; CMUMPS_DOUBLE *a_elt;

    /* Ordering, if given by user */
    CMUMPS_INT *perm_in;

    /* Orderings returned to user */
    /* symmetric permutation */
    CMUMPS_INT *sym_perm;
    /* column permutation */
    CMUMPS_INT *uns_perm;

    /* Scaling (input only in this version) */
    CMUMPS_DOUBLE *colsca; CMUMPS_DOUBLE *rowsca;
    /* RHS, solution, ouptput data and statistics */
    CMUMPS_DOUBLE *rhs, *redrhs, *rhs_sparse, *sol_loc;
    CMUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc;
    CMUMPS_INT nrhs, lrhs, lredrhs, nz_rhs, lsol_loc;
  CMUMPS_INT schur_mloc, schur_nloc, schur_lld;
  CMUMPS_INT mblock, nblock, nprow, npcol;
    CMUMPS_INT info[40],infog[40];
    CMUMPS_DOUBLE2 rinfo[20], rinfog[20];
    /* Null space */
    CMUMPS_INT deficiency; CMUMPS_INT * pivnul_list; CMUMPS_INT * mapping;
    /* Schur */
    CMUMPS_INT size_schur; CMUMPS_INT *listvar_schur; CMUMPS_DOUBLE *schur;
    /* Internal parameters */
    CMUMPS_INT instance_number;
    /* For out-of-core */
    char ooc_tmpdir[151];
    char ooc_prefix[151];
    char version_number[80];
  } CMUMPS_STRUC_C;


#if defined(UPPER) || defined(_WIN32)
#define cmumps_f77_ CMUMPS_F77
#define cmumps_affect_mapping_ CMUMPS_AFFECT_MAPPING
#define cmumps_affect_pivnul_list_ CMUMPS_AFFECT_PIVNUL_LIST
#define cmumps_affect_colsca_ CMUMPS_AFFECT_COLSCA
#define cmumps_affect_rowsca_ CMUMPS_AFFECT_ROWSCA 
#define cmumps_affect_uns_perm_     CMUMPS_AFFECT_UNS_PERM
#define cmumps_affect_sym_perm_     CMUMPS_AFFECT_SYM_PERM
#define cmumps_nullify_c_mapping_   CMUMPS_NULLIFY_C_MAPPING
#define cmumps_nullify_c_pivnul_list_ CMUMPS_NULLIFY_C_PIVNUL_LIST
#define cmumps_nullify_c_sym_perm_  CMUMPS_NULLIFY_C_SYM_PERM
#define cmumps_nullify_c_uns_perm_  CMUMPS_NULLIFY_C_UNS_PERM
#define cmumps_nullify_c_colsca_    CMUMPS_NULLIFY_C_COLSCA
#define cmumps_nullify_c_rowsca_    CMUMPS_NULLIFY_C_ROWSCA
#elif defined(Add__)
#define cmumps_f77_ cmumps_f77__
#define cmumps_affect_mapping_ cmumps_affect_mapping__
#define cmumps_affect_pivnul_list_ cmumps_affect_pivnul_list__
#define cmumps_affect_colsca_ cmumps_affect_colsca__
#define cmumps_affect_rowsca_ cmumps_affect_rowsca__
#define cmumps_affect_uns_perm_     cmumps_affect_uns_perm__     
#define cmumps_affect_sym_perm_     cmumps_affect_sym_perm__     
#define cmumps_nullify_c_mapping_   cmumps_nullify_c_mapping__    
#define cmumps_nullify_c_pivnul_list_ cmumps_nullify_c_pivnul_list__  
#define cmumps_nullify_c_sym_perm_  cmumps_nullify_c_sym_perm__   
#define cmumps_nullify_c_uns_perm_  cmumps_nullify_c_uns_perm__   
#define cmumps_nullify_c_colsca_    cmumps_nullify_c_colsca__     
#define cmumps_nullify_c_rowsca_    cmumps_nullify_c_rowsca__    
#elif defined(Add_)
/* Default. Nothing to do */
#else
/* Name without underscore is used */
#define cmumps_f77_ cmumps_f77
#define cmumps_affect_mapping_ cmumps_affect_mapping
#define cmumps_affect_pivnul_list_ cmumps_affect_pivnul_list
#define cmumps_affect_colsca_ cmumps_affect_colsca
#define cmumps_affect_rowsca_ cmumps_affect_rowsca
#define cmumps_affect_uns_perm_     cmumps_affect_uns_perm     
#define cmumps_affect_sym_perm_     cmumps_affect_sym_perm     
#define cmumps_nullify_c_mapping_   cmumps_nullify_c_mapping    
#define cmumps_nullify_c_pivnul_list_ cmumps_nullify_c_pivnul_list
#define cmumps_nullify_c_sym_perm_  cmumps_nullify_c_sym_perm   
#define cmumps_nullify_c_uns_perm_  cmumps_nullify_c_uns_perm   
#define cmumps_nullify_c_colsca_    cmumps_nullify_c_colsca     
#define cmumps_nullify_c_rowsca_    cmumps_nullify_c_rowsca
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

void MUMPS_CALL cmumps_c(CMUMPS_STRUC_C * cmumps_par);
void MUMPS_CALL cmumps_affect_mapping_(CMUMPS_INT * f77mapping);
void MUMPS_CALL cmumps_affect_pivnul_list_(CMUMPS_INT * f77pivnul_list);
void MUMPS_CALL cmumps_affect_uns_perm_(CMUMPS_INT * f77sym_perm);
void MUMPS_CALL cmumps_affect_sym_perm_(CMUMPS_INT * f77uns_perm);
void MUMPS_CALL cmumps_nullify_c_mapping_();
void MUMPS_CALL cmumps_nullify_c_pivnul_list_();
void MUMPS_CALL cmumps_nullify_c_sym_perm_();
void MUMPS_CALL cmumps_nullify_c_uns_perm_();
#ifdef return_scaling
void MUMPS_CALL cmumps_affect_colsca_(CMUMPS_DOUBLE * f77colsca);
void MUMPS_CALL cmumps_affect_rowsca_(CMUMPS_DOUBLE * f77rowsca);
void MUMPS_CALL cmumps_nullify_c_colsca_();
void MUMPS_CALL cmumps_nullify_c_rowsca_();
#endif

void MUMPS_CALL cmumps_f77_(CMUMPS_INT *job, CMUMPS_INT *sym, CMUMPS_INT *par, CMUMPS_INT * comm_fortran,
CMUMPS_INT *n, CMUMPS_INT *icntl,
CMUMPS_DOUBLE2 *cntl, CMUMPS_INT *nz, CMUMPS_INT *irn, CMUMPS_INT *irn_avail,
CMUMPS_INT *jcn, CMUMPS_INT *jcn_avail, CMUMPS_DOUBLE *a, CMUMPS_INT *a_avail,
  CMUMPS_INT *nz_loc, CMUMPS_INT * irn_loc, CMUMPS_INT * irn_loc_avail,
  CMUMPS_INT * jcn_loc, CMUMPS_INT * jcn_loc_avail, CMUMPS_DOUBLE *a_loc,
  CMUMPS_INT * a_loc_avail, CMUMPS_INT * nelt, CMUMPS_INT * eltptr,
  CMUMPS_INT * eltptr_avail, CMUMPS_INT * eltvar, CMUMPS_INT * eltvar_avail,
  CMUMPS_DOUBLE * a_elt, CMUMPS_INT * a_elt_avail, CMUMPS_INT * perm_in,
  CMUMPS_INT * perm_in_avail, CMUMPS_DOUBLE * rhs, CMUMPS_INT * rhs_avail,
  CMUMPS_DOUBLE * redrhs, CMUMPS_INT * redrhs_avail,
  CMUMPS_INT * info, CMUMPS_DOUBLE2 * rinfo, CMUMPS_INT * infog, CMUMPS_DOUBLE2 * rinfog,
  CMUMPS_INT * deficiency, CMUMPS_INT * size_schur, CMUMPS_INT * listvar_schur,
  CMUMPS_INT * listvar_schur_avail, CMUMPS_DOUBLE * schur,
  CMUMPS_INT * schur_avail, CMUMPS_DOUBLE * colsca, CMUMPS_INT * colsca_avail,
  CMUMPS_DOUBLE * rowsca, CMUMPS_INT * rowsca_avail, CMUMPS_INT * instance_number,
  CMUMPS_INT * nrhs, CMUMPS_INT * lrhs, CMUMPS_INT * lredrhs, CMUMPS_DOUBLE * rhs_sparse, CMUMPS_INT * rhs_sparse_avail,
  CMUMPS_DOUBLE * sol_loc, CMUMPS_INT * sol_loc_avail, CMUMPS_INT * irhs_sparse,
  CMUMPS_INT * irhs_sparse_avail, CMUMPS_INT * irhs_ptr, CMUMPS_INT * irhs_ptr_avail,
  CMUMPS_INT * isol_loc, CMUMPS_INT * isol_loc_avail, CMUMPS_INT * nz_rhs, CMUMPS_INT * lsol_loc
  ,CMUMPS_INT * schur_mloc, CMUMPS_INT *schur_nloc, CMUMPS_INT * schur_lld,
  CMUMPS_INT * schur_mblock, CMUMPS_INT * schur_nblock, CMUMPS_INT * schur_nprow,
  CMUMPS_INT * schur_npcol,
  CMUMPS_INT * ooc_tmpdir,
  CMUMPS_INT * ooc_prefix,
  CMUMPS_INT * ooc_tmpdirlen,
  CMUMPS_INT * ooc_prefixlen
  );
#endif

