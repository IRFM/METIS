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
/* $Id: smumps_io_basic_var.h,v 1.16 2006/12/20 09:41:25 aguermou Exp $ */

#include "smumps_io_basic.h"


#if ! defined (_WIN32) 
#if defined(WITH_PFUNC) && ! defined (WITHOUT_PTHREAD)
#include <pthread.h>
pthread_mutex_t smumps_io_pwrite_mutex;
#endif
/* int* smumps_io_pfile_pointer_array; */
/* int* smumps_io_current_file; */
/* #else /\*_WIN32*\/ */
/* FILE** smumps_io_current_file; */
/* FILE** smumps_io_pfile_pointer_array; */
#endif /*_WIN32*/

/* smumps_file_struct* smumps_io_pfile_pointer_array;
   smumps_file_struct* smumps_io_current_file; */

smumps_file_type* smumps_files;

/* int smumps_io_current_file_number; */
char* smumps_ooc_file_prefix;
/* char** smumps_io_pfile_name; */
/* int smumps_io_current_file_position; */
/* int smumps_io_write_pos; */
/* int smumps_io_last_file_opened; */
int smumps_elementary_data_size;
int smumps_io_is_init_called;
int smumps_io_myid;
/* int smumps_io_nb_file; */
int smumps_io_flag_async;
int smumps_io_k211;
/* int smumps_flag_open;*/
int smumps_directio_flag;
