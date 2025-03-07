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
/* $Id: dmumps_io_basic_extern.h,v 1.14 2006/12/20 09:41:25 aguermou Exp $ */
#include "dmumps_io_basic.h"



#if ! defined (_WIN32) 
#if defined (WITH_PFUNC) && ! defined (WITHOUT_PTHREAD)
#include <pthread.h>
extern pthread_mutex_t dmumps_io_pwrite_mutex;
#endif
/* extern int* dmumps_io_pfile_pointer_array; */
/* extern int* dmumps_io_current_file; */
/* #else /\*_WIN32*\/ */
/* extern FILE** dmumps_io_current_file; */
/* extern FILE** dmumps_io_pfile_pointer_array; */
#endif /*_WIN32*/

/*extern dmumps_file_struct* dmumps_io_pfile_pointer_array;
  extern dmumps_file_struct* dmumps_io_current_file;*/
extern dmumps_file_type* dmumps_files;

/* extern int dmumps_io_current_file_number; */
extern char* dmumps_ooc_file_prefix;
/* extern char** dmumps_io_pfile_name; */
/* extern int dmumps_io_current_file_position; */
/* extern int dmumps_io_write_pos; */
/* extern int dmumps_io_last_file_opened; */
extern int dmumps_elementary_data_size;
extern int dmumps_io_is_init_called;
extern int dmumps_io_myid;
/* extern int dmumps_io_nb_file; */
extern int dmumps_io_flag_async;
extern int dmumps_io_k211;
/* extern int dmumps_flag_open; */
extern int directio_flag;
extern int dmumps_directio_flag;
