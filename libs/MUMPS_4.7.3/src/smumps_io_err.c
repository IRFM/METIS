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
/*    $Id: smumps_io_err.c,v 1.9 2006/07/21 16:02:25 aguermou Exp $  */

#include "smumps_io_err_var.h"
#include "smumps_io_basic_extern.h"

  
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)

int smumps_io_protect_err(){
  if(smumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_lock(&err_mutex);
  }
  return 0;
}
int smumps_io_unprotect_err(){
  if(smumps_io_flag_async==IO_ASYNC_TH){
    pthread_mutex_unlock(&err_mutex);
  }
  return 0;
}

int smumps_io_init_err_lock(){
  pthread_mutex_init(&err_mutex,NULL);
  return 0;
}
int smumps_io_destroy_err_lock(){
  pthread_mutex_destroy(&err_mutex);
  return 0;
}

int smumps_check_error_th(){
  return err_flag;
}

#endif /* _WIN32 && WITHOUT_PTHREAD */

int smumps_io_prop_err_info(int ierr){
  /* Copies the error description string in a fortran character
     array. */
  int i;
  for(i=0;i<(int)strlen(error_str);i++){
    smumps_err[i]=error_str[i];
  }
  *dim_mumps_err=(int)strlen(smumps_err);
  return 0;
}

int smumps_io_build_err_str(int errnum, int smumps_err,const char* desc,char* buf,int size){
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
  smumps_io_protect_err();
#endif
  if(err_flag==0){
    sprintf(buf,"%s: %s",desc,strerror(errnum));
    err_flag=smumps_err;
  }
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)  
  smumps_io_unprotect_err();
#endif
  return 0;
}

void MUMPS_CALL smumps_low_level_init_err_str(int * dim, char* err_str, smumps_ftnlen l1){
  /* Keeps a C pointer to store error description string that will be
     displayed by the fortran layers*/
  smumps_err=err_str;
  dim_mumps_err=dim;
  err_flag=0;
  return;
}

