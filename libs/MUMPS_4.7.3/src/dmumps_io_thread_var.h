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
/* $Id: dmumps_io_thread_var.h,v 1.4 2006/07/21 17:30:28 jylexcel Exp $ */

#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#include "dmumps_io_thread.h"

int io_flag_stop,current_req_num;
pthread_t io_thread,main_thread;
pthread_mutex_t io_mutex;


pthread_cond_t cond_io,cond_nb_free_finished_requests,cond_nb_free_active_requests,cond_stop;
pthread_mutex_t io_mutex_cond;
int int_sem_io,int_sem_nb_free_finished_requests,int_sem_nb_free_active_requests,int_sem_stop;


int with_sem;


struct request_io *io_queue;
int first_active,last_active,nb_active;
int *finished_requests_inode,*finished_requests_id,first_finished_requests,
  last_finished_requests,nb_finished_requests,smallest_request_id;
int dmumps_owns_mutex;
int test_request_called_from_mumps;

#endif
