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
/*    $Id: zmumps_io.h,v 1.30 2006/12/06 12:05:49 jylexcel Exp $ */

#if defined(_WIN32) || defined (UPPER)
#define zmumps_is_there_finished_request ZMUMPS_IS_THERE_FINISHED_REQUEST
#define zmumps_clean_request ZMUMPS_CLEAN_REQUEST
#define zmumps_test_request ZMUMPS_TEST_REQUEST
#define zmumps_wait_request ZMUMPS_WAIT_REQUEST
#define zmumps_wait_all_requests ZMUMPS_WAIT_ALL_REQUESTS
#define zmumps_low_level_init_prefix ZMUMPS_LOW_LEVEL_INIT_PREFIX
#define zmumps_low_level_init_tmpdir ZMUMPS_LOW_LEVEL_INIT_TMPDIR
#define zmumps_low_level_init_ooc_c ZMUMPS_LOW_LEVEL_INIT_OOC_C
#define zmumps_low_level_write_ooc_c ZMUMPS_LOW_LEVEL_WRITE_OOC_C
#define zmumps_low_level_read_ooc_c ZMUMPS_LOW_LEVEL_READ_OOC_C
#define zmumps_low_level_direct_read ZMUMPS_LOW_LEVEL_DIRECT_READ
#define zmumps_clean_io_data_c ZMUMPS_CLEAN_IO_DATA_C
#define zmumps_get_max_nb_req ZMUMPS_GET_MAX_NB_REQ
#define zmumps_get_max_file_size ZMUMPS_GET_MAX_FILE_SIZE
#define zmumps_ooc_get_nb_files ZMUMPS_OOC_GET_NB_FILES
#define zmumps_ooc_get_file_name ZMUMPS_OOC_GET_FILE_NAME
#define zmumps_ooc_set_file_name ZMUMPS_OOC_SET_FILE_NAME
#define zmumps_ooc_start_low_level ZMUMPS_OOC_START_LOW_LEVEL
#define zmumps_ooc_alloc_pointers ZMUMPS_OOC_ALLOC_POINTERS
#define zmumps_ooc_print_stats ZMUMPS_OOC_PRINT_STATS
#define zmumps_ooc_remove_file ZMUMPS_OOC_REMOVE_FILE
#define zmumps_ooc_init_vars ZMUMPS_OOC_INIT_VARS
#define zmumps_ooc_end_write ZMUMPS_OOC_END_WRITE
#elif defined(Add__)
#define zmumps_is_there_finished_request zmumps_is_there_finished_request__
#define zmumps_clean_request zmumps_clean_request__
#define zmumps_test_request zmumps_test_request__
#define zmumps_wait_request zmumps_wait_request__
#define zmumps_wait_all_requests zmumps_wait_all_requests__
#define zmumps_low_level_init_prefix zmumps_low_level_init_prefix__
#define zmumps_low_level_init_tmpdir zmumps_low_level_init_tmpdir__
#define zmumps_low_level_init_ooc_c zmumps_low_level_init_ooc_c__
#define zmumps_low_level_write_ooc_c zmumps_low_level_write_ooc_c__
#define zmumps_low_level_read_ooc_c zmumps_low_level_read_ooc_c__
#define zmumps_low_level_direct_read zmumps_low_level_direct_read__
#define zmumps_clean_io_data_c zmumps_clean_io_data_c__
#define zmumps_get_max_nb_req zmumps_get_max_nb_req__
#define zmumps_get_max_file_size zmumps_get_max_file_size__ 
#define zmumps_ooc_get_nb_files zmumps_ooc_get_nb_files__
#define zmumps_ooc_get_file_name zmumps_ooc_get_file_name__
#define zmumps_ooc_set_file_name zmumps_ooc_set_file_name__
#define zmumps_ooc_start_low_level zmumps_ooc_start_low_level__
#define zmumps_ooc_alloc_pointers zmumps_ooc_alloc_pointers__
#define zmumps_ooc_print_stats zmumps_ooc_print_stats__
#define zmumps_ooc_remove_file zmumps_ooc_remove_file__
#define zmumps_ooc_init_vars zmumps_ooc_init_vars__
#define zmumps_ooc_end_write zmumps_ooc_end_write__
#elif defined(Add_)
#define zmumps_is_there_finished_request zmumps_is_there_finished_request_
#define zmumps_clean_request zmumps_clean_request_
#define zmumps_test_request zmumps_test_request_
#define zmumps_wait_request zmumps_wait_request_
#define zmumps_wait_all_requests zmumps_wait_all_requests_
#define zmumps_low_level_init_prefix zmumps_low_level_init_prefix_
#define zmumps_low_level_init_tmpdir zmumps_low_level_init_tmpdir_
#define zmumps_low_level_init_ooc_c zmumps_low_level_init_ooc_c_
#define zmumps_low_level_write_ooc_c zmumps_low_level_write_ooc_c_
#define zmumps_low_level_read_ooc_c zmumps_low_level_read_ooc_c_
#define zmumps_low_level_direct_read zmumps_low_level_direct_read_
#define zmumps_clean_io_data_c zmumps_clean_io_data_c_
#define zmumps_get_max_nb_req zmumps_get_max_nb_req_
#define zmumps_get_max_file_size zmumps_get_max_file_size_ 
#define zmumps_ooc_get_nb_files zmumps_ooc_get_nb_files_
#define zmumps_ooc_get_file_name zmumps_ooc_get_file_name_
#define zmumps_ooc_set_file_name zmumps_ooc_set_file_name_
#define zmumps_ooc_start_low_level zmumps_ooc_start_low_level_
#define zmumps_ooc_alloc_pointers zmumps_ooc_alloc_pointers_
#define zmumps_ooc_print_stats zmumps_ooc_print_stats_
#define zmumps_ooc_remove_file zmumps_ooc_remove_file_
#define zmumps_ooc_init_vars zmumps_ooc_init_vars_
#define zmumps_ooc_end_write zmumps_ooc_end_write_
#endif

#include "zmumps_compat.h"

/*
 *  Two character arrays that are used by low_level_init_prefix
 *  and low_level_init_tmpdir to store intermediate file names.
 *  They are passed to zmumps_io_basic.c inside the routine
 *  zmumps_low_level_init_ooc_c.
 *  Note that both low_level_init_prefix and low_level_init_tmpdir
 *  MUST be called before low_level_init_ooc_c.
 * 
 */
static char zmumps_ooc_store_prefix[150]; static int zmumps_ooc_store_prefixlen=-1;
static char zmumps_ooc_store_tmpdir[150]; static int zmumps_ooc_store_tmpdirlen=-1;

void MUMPS_CALL zmumps_low_level_init_prefix(int * dim, char * str, zmumps_ftnlen l1);
void MUMPS_CALL zmumps_low_level_init_tmpdir(int * dim, char * str, zmumps_ftnlen l1);
void MUMPS_CALL zmumps_low_level_init_ooc_c(int* _myid, int* total_size_io,int* size_element,
                               int* async, int* k211, int* ierr);


void MUMPS_CALL zmumps_is_there_finished_request(int* flag,int* ierr);

void MUMPS_CALL zmumps_clean_request(int* request_id,int* ierr);

void MUMPS_CALL zmumps_test_request(int* request_id,int* flag,int* ierr);

void MUMPS_CALL zmumps_wait_request(int* request_id,int* ierr);

void MUMPS_CALL zmumps_wait_all_requests(int* ierr);

void MUMPS_CALL zmumps_low_level_write_ooc_c( const int * strat_IO, 
                                 void * address_block,
                                 int * block_size,
                                 int * pos_in_file,
                                 int * file_number,
                                 int * inode,
                                 int * request_arg,
                                 int * type,
                                 int * ierr);

void MUMPS_CALL zmumps_low_level_read_ooc_c( const int * strat_IO, 
                                 void * address_block,
                                 int * block_size,
                                 int * from_where,
                                 int * file_number,
                                 int * inode,
                                 int * request_arg,
                                 int * type,
                                 int * ierr);

void MUMPS_CALL zmumps_low_level_direct_read(void * address_block,
                                 int * block_size,
                                 int * from_where,
                                 int * file_number,
                                 int * type,
                                 int * ierr);

void MUMPS_CALL zmumps_clean_io_data_c(int* myid,int* step,int* ierr);

void MUMPS_CALL zmumps_get_max_nb_req(int *max,int* ierr);

void MUMPS_CALL zmumps_get_max_file_size(double * max_ooc_file_size);

void MUMPS_CALL zmumps_ooc_get_nb_files(int* nb_files);

void MUMPS_CALL zmumps_ooc_get_file_name(int* indice, int * length, char* name, zmumps_ftnlen l1);

void MUMPS_CALL zmumps_ooc_set_file_name(int* indice, int* length, int* ierr, char* name, zmumps_ftnlen l1);

void MUMPS_CALL zmumps_ooc_alloc_pointers(int* dim,int* ierr);

void MUMPS_CALL zmumps_ooc_init_vars(int* myid_arg, int* nb_file_arg,
                                   int* size_element,int* async, int* k211,
                                   int *ierr);

void MUMPS_CALL zmumps_ooc_start_low_level(int* ierr);

void MUMPS_CALL zmumps_ooc_print_stats();

void MUMPS_CALL zmumps_ooc_remove_file(int *ierr, char *name, zmumps_ftnlen l1);

void MUMPS_CALL zmumps_ooc_end_write(int *ierr);
