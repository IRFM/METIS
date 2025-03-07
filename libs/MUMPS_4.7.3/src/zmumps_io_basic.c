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
/*    $Id: zmumps_io_basic.c,v 1.92 2007/02/28 14:52:30 jylexcel Exp $  */

#include "zmumps_io_basic_var.h"
#include "zmumps_io_err_extern.h"

int _zmumps_next_file(int type){
  /* Defines the pattern for the file name. The last 6 'X' will be replaced
     so as to name were unique */
  char name[150];
  int ret_code, fd;  
  zmumps_file_struct  *zmumps_io_pfile_pointer_array;

#ifndef _WIN32
  strcpy(name,zmumps_ooc_file_prefix);
  fd=mkstemp(name);  
  /* Note that a file name is built by mkstemp and that the file is 
     opened. fd hold the file descriptor to access it.
     We want to close the file that will be opened later
     and might be removed before the end of the processus.
    */
    if(fd<0){
      ret_code=-99;
      return ret_code;
    }
    else{ 
      close(fd); 
    }

#else
  sprintf(name,"%s_%d",zmumps_ooc_file_prefix,((zmumps_files+type)->zmumps_io_current_file_number)+1);
#endif
  
  if ((zmumps_files+type)->zmumps_io_current_file_number >= ((zmumps_files+type)->zmumps_io_nb_file)-1){
    /* Exception : probably thrown because of a bad estimation
       of number of files. */
    /* We increase the number of file needed and then realloc. */
    ((zmumps_files+type)->zmumps_io_nb_file)++;
/* #ifndef _WIN32 */
/*     zmumps_io_pfile_pointer_array=realloc((void *)zmumps_io_pfile_pointer_array,zmumps_io_nb_file*sizeof(int)); */
/* #else */
/*     zmumps_io_pfile_pointer_array=realloc((void *)zmumps_io_pfile_pointer_array,zmumps_io_nb_file*sizeof(FILE*)); */
/* #endif */
    (zmumps_files+type)->zmumps_io_pfile_pointer_array=realloc((void *)(zmumps_files+type)->zmumps_io_pfile_pointer_array,((zmumps_files+type)->zmumps_io_nb_file)*sizeof(zmumps_file_struct));
    if((zmumps_files+type)->zmumps_io_pfile_pointer_array==NULL){
      sprintf(error_str,"Allocation problem in low-level OOC layer\n");
      return -13;
    }
/*     zmumps_io_pfile_name=realloc((void*)zmumps_io_pfile_name,(zmumps_io_nb_file)*sizeof(char *)); */
/*     if(zmumps_io_pfile_name==NULL){ */
/*       sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/*       return -13; */
/*     } */
  }
  zmumps_io_pfile_pointer_array=(zmumps_files+type)->zmumps_io_pfile_pointer_array;


  /* Do change the current file */
  ((zmumps_files+type)->zmumps_io_current_file_number)++;
  /*  *(zmumps_io_pfile_pointer_array+zmumps_io_current_file_number)=fopen(name,"w+"); */
/*   *(zmumps_io_pfile_name+zmumps_io_current_file_number)=(char *)malloc((strlen(name)+1)*sizeof(char)); */

/*   if(*(zmumps_io_pfile_name+zmumps_io_current_file_number)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/*     return -13; */
/*   } */

  strcpy((zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number)->name,name);
  /* See zmumps_io_basic.h for comments on the I/O flags passed to open */
#ifndef _WIN32
   (zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number)->file=open(name,(zmumps_files+type)->zmumps_flag_open); 
  /* 
CPA: for LU factor file: 
(dmumps_io_pfile_pointer_array+(dmumps_files+type)->dmumps_io_current_file_number)->file= open(name, O_WRONLY | O_CREAT | O_TRUNC, 0666); */


  if((zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number)->file==-1){
    zmumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200);
    return -90;
  }
#else
  (zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number)->file=fopen(name,(zmumps_files+type)->zmumps_flag_open);
  if((zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number)->file==NULL){
    /*    zmumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200); */
    sprintf(error_str,"Problem while opening OOC file");
    return -90;
   }
#endif

  
  (zmumps_files+type)->zmumps_io_current_file=(zmumps_io_pfile_pointer_array+(zmumps_files+type)->zmumps_io_current_file_number);
  ((zmumps_files+type)->zmumps_io_last_file_opened)++;
  /*  if(*(zmumps_io_pfile_pointer_array+zmumps_io_current_file_number)==NULL) */
  ((zmumps_files+type)->zmumps_io_current_file)->write_pos=0;
  /*  printf("new file created %d\n",zmumps_io_current_file_number);*/
  return 0;
}

void _zmumps_update_current_file_position(zmumps_file_struct* file_arg){
  file_arg->current_pos=file_arg->write_pos;
/*   zmumps_io_current_file_position=zmumps_io_write_pos; */
}


int _zmumps_compute_where_to_write(const double to_be_written,const int type){
  /* Check if the current file has enough memory to receive the whole block*/
  double size;
  int ret_code;
  zmumps_file_struct *current_file;
  current_file=(zmumps_files+type)->zmumps_io_current_file;
  _zmumps_update_current_file_position(current_file);
  size=(double)MAX_FILE_SIZE-((double)current_file->current_pos+to_be_written);
  if (size > 0){
  }else{
    ret_code=_zmumps_next_file(type);
    if(ret_code<0){
      return ret_code;
    }
    current_file=(zmumps_files+type)->zmumps_io_current_file;
    _zmumps_update_current_file_position(current_file);
  }

  return 0;
}

int _zmumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number,const int type){
  int ret_code;
  ret_code=_zmumps_compute_where_to_write(to_be_written,type);
  if(ret_code<0){
    return ret_code;
  }
  *pos_in_file=((zmumps_files+type)->zmumps_io_current_file)->current_pos;
  /* should be modified to take into account the file arg */
  *file_number=(zmumps_files+type)->zmumps_io_current_file_number;
  return 0;
}

int _zmumps_compute_nb_concerned_files(int * block_size, int * nb_concerned_files){
  *nb_concerned_files=(int)ceil((((double)(*block_size))*((double)(zmumps_elementary_data_size)))/((double)MAX_FILE_SIZE));
  return 0;
}

int zmumps_io_do_write_block(void * address_block,
		     int * block_size,
		     int * pos_in_file,
		     int * file_number,		     
         	     int * type_arg,
		     int * ierr){   
  /* Type of fwrite : size_t fwrite(const void *ptr, size_t size, 
                                    *size_t nmemb, FILE *stream); */
  size_t write_size;
  int i;
  int nb_concerned_files=0;
  int ret_code,file_number_loc,pos_in_file_loc;
  double to_be_written;
#ifndef _WIN32
  int* file;
#else
  FILE** file;
#endif
  int where;
  void* loc_addr;
  int type;
  type=*type_arg;

  loc_addr=address_block;
  _zmumps_compute_nb_concerned_files(block_size,&nb_concerned_files);
  to_be_written=((double)zmumps_elementary_data_size)*((double)(*block_size));
  if((nb_concerned_files==0)&&(to_be_written==0)){
    *file_number=(zmumps_files+type)->zmumps_io_current_file_number;
    *pos_in_file=((zmumps_files+type)->zmumps_io_current_file)->write_pos;
    return 0;
  }

  for(i=0;i<nb_concerned_files;i++){

#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(zmumps_io_flag_async==IO_ASYNC_TH){
      zmumps_io_protect_pointers();
    }
#endif
#endif
    ret_code=_zmumps_prepare_pointers_for_write(to_be_written,&pos_in_file_loc,&file_number_loc,type);
    
    if(ret_code<0){
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(zmumps_io_flag_async==IO_ASYNC_TH){
      zmumps_io_unprotect_pointers();
    }
#endif
#endif
      return ret_code;
    }

    if(i==0){
      *pos_in_file=pos_in_file_loc;
      if(nb_concerned_files==1){
	*file_number=file_number_loc;
      }else{
	*file_number=-file_number_loc;
      }
    }
    if((double)(MAX_FILE_SIZE-((zmumps_files+type)->zmumps_io_current_file)->write_pos)>to_be_written){
      write_size=(size_t)to_be_written;
    }else{
      write_size=(size_t)((double)(MAX_FILE_SIZE-((zmumps_files+type)->zmumps_io_current_file)->write_pos));
    }
#if defined (_WIN32)
    write_size=(size_t)(int)((write_size)/zmumps_elementary_data_size);
#endif
    file=&(((zmumps_files+type)->zmumps_io_current_file)->file);
    where=((zmumps_files+type)->zmumps_io_current_file)->write_pos;
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(zmumps_io_flag_async==IO_ASYNC_TH){
      zmumps_io_unprotect_pointers();
    }
#endif
#endif
    /* printf("1 write -> size = %d | off = %d | file = %d\n",(int)write_size,where,*file_number);*/
    ret_code=zmumps_io_write__(file,loc_addr,write_size,where);
    
    if(ret_code<0){
      return ret_code;
    }
    
#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(zmumps_io_flag_async==IO_ASYNC_TH){
      zmumps_io_protect_pointers();
    }
#endif
#endif

#ifndef _WIN32
    ((zmumps_files+type)->zmumps_io_current_file)->write_pos=((zmumps_files+type)->zmumps_io_current_file)->write_pos+((int)write_size);
    to_be_written=to_be_written-((int)write_size);
    loc_addr=(void*)((size_t)loc_addr+write_size);
/*     zmumps_io_write_pos=zmumps_io_write_pos+((int)write_size); */
/*     to_be_written=to_be_written-((int)write_size); */
/*     loc_addr=(void*)((size_t)loc_addr+write_size); */
#else
    /* fread and write */
    ((zmumps_files+type)->zmumps_io_current_file)->write_pos=((zmumps_files+type)->zmumps_io_current_file)->write_pos+((int)write_size*zmumps_elementary_data_size);
    to_be_written=to_be_written-((int)write_size*zmumps_elementary_data_size);
    loc_addr=(void*)((size_t)loc_addr+(size_t)((int)write_size*zmumps_elementary_data_size));
/*     zmumps_io_write_pos=zmumps_io_write_pos+((int)write_size*zmumps_elementary_data_size); */
/*     to_be_written=to_be_written-((int)write_size*zmumps_elementary_data_size); */
/*     loc_addr=(void*)((size_t)loc_addr+(size_t)((int)write_size*zmumps_elementary_data_size)); */
#endif

#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC
    if(zmumps_io_flag_async==IO_ASYNC_TH){
      zmumps_io_unprotect_pointers();
    }
#endif
#endif
  }
  if(to_be_written!=0){
    sprintf(error_str,"Internal (1) error in low-level I/O operation %lf",to_be_written);
    return -90;
  }
  
  return 0;
}

int zmumps_io_do_read_block(void * address_block,
	            int * block_size,
                    int * from_where,
		    int * file_number,
                    int * type_arg,
                    int * ierr){
  int ret_code;
  size_t size;
#ifndef _WIN32
  int* file;
#else
  FILE** file;
#endif
  double read_size;
  int local_fnum,local_offset;
  void *loc_addr;
  int type;
  type=*type_arg;
  /*  if(((double)(*block_size))*((double)(zmumps_elementary_data_size))>(double)MAX_FILE_SIZE){
    sprintf(error_str,"Internal error in low-level I/O operation (requested size too big for file system) \n");
    return -90;
    }*/

  if(*block_size==0){
    return 0;
  }
  read_size=(double)zmumps_elementary_data_size*(double)(*block_size);
  local_fnum=*file_number;
  if(local_fnum<0)
    local_fnum=-local_fnum;
  /*  if((*file_number<0)&&(read_size<(double)MAX_FILE_SIZE)){
    sprintf(error_str,"Internal error (1) in low level read op\n");
    return -90;
    }*/
  local_offset=*from_where;
  loc_addr=address_block;

  while(read_size>0){
    file=&((((zmumps_files+type)->zmumps_io_pfile_pointer_array)+local_fnum)->file);
#ifndef _WIN32
    if(read_size+(double)local_offset>(double)MAX_FILE_SIZE){
      size=(size_t)MAX_FILE_SIZE-(size_t)local_offset;
    }else{
      size=(size_t)read_size;
    }
#else
    if(read_size+(double)local_offset>(double)MAX_FILE_SIZE){
      size=((size_t)MAX_FILE_SIZE-(size_t)local_offset)/(size_t)zmumps_elementary_data_size;
    }else{
      size=(size_t)(read_size/zmumps_elementary_data_size);
    }
#endif
    ret_code=zmumps_io_read__(file,loc_addr,size,local_offset);

    if(ret_code<0){
      return ret_code;
    }

#ifdef _WIN32
    size=size*zmumps_elementary_data_size;
#endif

    read_size=read_size-(double)size;
    loc_addr=(void*)((size_t)loc_addr+size);
    local_fnum++;
    local_offset=0;
    if(local_fnum>(zmumps_files+type)->zmumps_io_nb_file){
      sprintf(error_str,"Internal error (2) in low level read op\n");
      return -90;
    }
  }
  return 0;
}
int zmumps_free_file_pointers(int *step){
  int i,j,bound,ierr;
  if(*step==0){
    /* factorization */
    bound=NB_FILE_TYPE_FACTO;
  }else{
    /* solve */
    bound=NB_FILE_TYPE_SOLVE;
  }
  for(j=0;j<bound;j++){
    for(i=0;i<=(zmumps_files+j)->zmumps_io_last_file_opened;i++){
#ifndef _WIN32
#ifdef SOLVE_CRAY
#endif
      ierr=close((((zmumps_files+j)->zmumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
	zmumps_io_build_err_str(errno,-90,"Problem while closing OOC file",error_str,200);
	return -90;
      }
#else
      ierr=fclose((((zmumps_files+j)->zmumps_io_pfile_pointer_array)+i)->file);
      if(ierr==-1){
	sprintf(error_str,"Problem while opening OOC file\n");
	return -90;
      }    
#endif
      /*     free(*(zmumps_io_pfile_name+i)); */
    } 
    free((zmumps_files+j)->zmumps_io_pfile_pointer_array);
  }
/*   free(zmumps_io_pfile_name); */
  free(zmumps_files);
/*   Free prefix only for facto  */
  if (*step == 0) free(zmumps_ooc_file_prefix);
  return 0;
}

int zmumps_init_file_structure(int* _myid, int* total_size_io,int* size_element){
  /* Computes the number of files needed. Uses ceil value. */
  int ierr;
#ifndef _WIN32
  int k211_loc;
  int zmumps_flag_open;
#endif
  int i,nb;
  int zmumps_io_nb_file;


  zmumps_io_nb_file=(int)((((double)(*total_size_io))*((double)(*size_element)))/(double)MAX_FILE_SIZE)+1;

  zmumps_directio_flag=0;

#ifndef _WIN32 
  zmumps_flag_open=0;
#endif



  zmumps_io_myid=*_myid;

  zmumps_elementary_data_size=*size_element;
  /* Allocates the memory necessary to handle the file pointer array.*/
  zmumps_files=(zmumps_file_type *)malloc(NB_FILE_TYPE_FACTO*sizeof(zmumps_file_type));
  if(zmumps_files==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  nb=zmumps_io_nb_file;
  for(i=0;i<NB_FILE_TYPE_FACTO;i++){
    (zmumps_files+i)->zmumps_io_nb_file=nb;
#ifndef _WIN32
    (zmumps_files+i)->zmumps_flag_open=zmumps_flag_open|O_WRONLY;
#else
    strcpy((zmumps_files+i)->zmumps_flag_open,"wb");
#endif
    ierr=zmumps_io_alloc_file_struct(&nb,i);
    if(ierr<0){
      return ierr;
    }
    ierr=_zmumps_next_file(i);
    if(ierr<0){
      return ierr;
    }
    nb=1;
  }
  /* Init the current file.*/
  return 0;
}
int zmumps_init_file_name(char* zmumps_dir,char* zmumps_file,
			 int* zmumps_dim_dir,int* zmumps_dim_file,int* _myid){
  int i;
  char *tmp_dir,*tmp_fname;
  char base_name[20];
  int dir_flag=0,file_flag=0;
  char zmumps_base[10]="zmumps_";
  tmp_dir=(char *)malloc(((*zmumps_dim_dir)+1)*sizeof(char));
  if(tmp_dir==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  tmp_fname=(char *)malloc(((*zmumps_dim_file)+1)*sizeof(char));
  if(tmp_fname==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  for(i=0;i<*zmumps_dim_dir;i++){
    tmp_dir[i]=zmumps_dir[i];
  }
  tmp_dir[i]=0;
  for(i=0;i<*zmumps_dim_file;i++){
    tmp_fname[i]=zmumps_file[i];
  }
  tmp_fname[i]=0;  
  if(strcmp(tmp_dir,UNITIALIZED)==0){
    dir_flag=1;
    free(tmp_dir);
    tmp_dir=getenv("MUMPS_OOC_TMPDIR");
    if(tmp_dir==NULL){
#ifdef SP_
#ifndef CINES_
      tmp_dir=getenv("TMPDIR");
      if(tmp_dir==NULL){
	tmp_dir=ZMUMPS_OOC_DEFAULT_DIR;
      }
#else
      tmp_dir=ZMUMPS_OOC_DEFAULT_DIR;       
#endif       
#else
      tmp_dir=ZMUMPS_OOC_DEFAULT_DIR;
#endif      
    }
  }
  if(strcmp(tmp_fname,UNITIALIZED)==0){
    free(tmp_fname);
    tmp_fname=getenv("MUMPS_OOC_PREFIX");
    file_flag=1;
  }

  if(tmp_fname!=NULL){
#ifndef _WIN32
      sprintf(base_name,"_%s%d_XXXXXX",zmumps_base,*_myid);
#else
      sprintf(base_name,"_%s%d",zmumps_base,*_myid);
#endif
      zmumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(tmp_fname)+strlen(base_name)+1+1)*sizeof(char));
      if(zmumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(zmumps_ooc_file_prefix,"%s%s%s%s",tmp_dir,SEPARATOR,tmp_fname,base_name);
  }else{
#ifndef _WIN32
    sprintf(base_name,"%s%s%d_XXXXXX",SEPARATOR,zmumps_base,*_myid);
#else
    sprintf(base_name,"%s%s%d",SEPARATOR,zmumps_base,*_myid);
#endif
      zmumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(base_name)+1)*sizeof(char));
      if(zmumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(zmumps_ooc_file_prefix,"%s%s%s",tmp_dir,SEPARATOR,base_name);
  }  
  if(!dir_flag){
    free(tmp_dir);
  }
  if(!file_flag){
    free(tmp_fname);
  }
  return 0;
}

int zmumps_io_alloc_file_struct(int* nb,int which){
  /* allocate structure for files corresponding to factors */
  (zmumps_files+which)->zmumps_io_pfile_pointer_array=(zmumps_file_struct *)malloc((*nb)*sizeof(zmumps_file_struct));
  if((zmumps_files+which)->zmumps_io_pfile_pointer_array==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  


/*   zmumps_io_pfile_name=(char **)malloc((*nb)*sizeof(char *)); */
/*   if(zmumps_io_pfile_name==NULL){ */
/* #ifndef _WIN32       */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/* #else */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer\n"); */
/* #endif     */
/*     return -13; */
/*     /\*    return -1;*\/ */
/*   } */
  (zmumps_files+which)->zmumps_io_current_file_number = -1;
  (zmumps_files+which)->zmumps_io_last_file_opened=-1;
  return 0;
}

int zmumps_io_get_nb_files(int* nb_files,int* type){
  *nb_files=((zmumps_files+*type)->zmumps_io_last_file_opened)+1;
  return 0;
}

int zmumps_io_get_file_name(int* indice,char* name,int* length,int* type){
  int i;
  i=(*indice)-1;
  strcpy(name,(((zmumps_files+*type)->zmumps_io_pfile_pointer_array)+i)->name);
  *length=strlen(name)+1;
  return 0;  
}

int zmumps_io_alloc_pointers(int * dim){
  int ierr;
  int i;
  /* This is called by solve step, we have only one type of files */
  zmumps_files=(zmumps_file_type *)malloc(sizeof(zmumps_file_type));
  zmumps_files->zmumps_io_nb_file=*dim;  
  if(zmumps_files==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  /* FIXME : This loop must be modified when NB_FILE_TYPE_SOLVE is
     greater than 1 (the dim parameter is not good) */
  for(i=0;i<NB_FILE_TYPE_SOLVE;i++){
    ierr=zmumps_io_alloc_file_struct(dim,i);
    if(ierr<0){
      return ierr;
    }
  }
  return 0;
}

int zmumps_io_init_vars(int* myid_arg, int* nb_file_arg,int* size_element,int* async_arg){

#ifndef _WIN32
  int k211_loc;
  int zmumps_flag_open;
#endif
  int i;
  
  zmumps_directio_flag=0;

#ifndef _WIN32 
  zmumps_flag_open=0;
#endif


  /* must be changed when we will have more than one file type during solve step */

  for(i=0;i<NB_FILE_TYPE_FACTO;i++){
#ifndef _WIN32
    (zmumps_files+i)->zmumps_flag_open=zmumps_flag_open|O_RDONLY;
#else
    strcpy((zmumps_files+i)->zmumps_flag_open,"rb");
#endif
  }

  zmumps_io_myid=*myid_arg;
  zmumps_elementary_data_size=*size_element;
  zmumps_io_flag_async=*async_arg;

  return 0;
}

int zmumps_io_set_file_name(int* indice,char* name,int* length,int* type){
  int i;
  i=(*indice)-1;
/*   *(zmumps_io_pfile_name+i)=(char *) malloc((*length)*strlen(name)); */
/*   if(*(zmumps_io_pfile_name+i)==NULL){ */
/*     sprintf(error_str,"Allocation problem in low-level OOC layer"); */
/*     return -13; */
/*   } */
  strcpy((((zmumps_files+*type)->zmumps_io_pfile_pointer_array)+i)->name,name);
  return 0;  
}

int zmumps_io_open_files_for_read(){
  int i,j;
  zmumps_file_struct  *zmumps_io_pfile_pointer_array;
#ifdef IRIX64_
  struct dioattr dio;
#endif
  for(j=0;j<NB_FILE_TYPE_SOLVE;j++){
    zmumps_io_pfile_pointer_array=(zmumps_files+j)->zmumps_io_pfile_pointer_array;  
    for(i=0;i<(zmumps_files+j)->zmumps_io_nb_file;i++){
#ifndef _WIN32
      (zmumps_io_pfile_pointer_array+i)->file=open((zmumps_io_pfile_pointer_array+i)->name,(zmumps_files+j)->zmumps_flag_open);
      
      if((zmumps_io_pfile_pointer_array+i)->file==-1){
	zmumps_io_build_err_str(errno,-90,"Problem while opening OOC file",error_str,200);
	return -90;
      }
#else
      (zmumps_io_pfile_pointer_array+i)->file=fopen((zmumps_io_pfile_pointer_array+i)->name,(zmumps_files+j)->zmumps_flag_open);
      if((zmumps_io_pfile_pointer_array+i)->file==NULL){
	sprintf(error_str,"Problem while opening OOC file");
	return -90;
      }
#endif
    }
  }
  return 0;
}

int zmumps_io_set_last_file(int* dim,int* type){
  (zmumps_files+*type)->zmumps_io_last_file_opened=*dim-1;
  return 0;
}

#if ! defined (_WIN32) && ! defined (WITHOUT_PTHREAD)
#ifdef WITH_PFUNC

int zmumps_io_protect_pointers(){
  pthread_mutex_lock(&zmumps_io_pwrite_mutex);
  return 0;
}
int zmumps_io_unprotect_pointers(){
  pthread_mutex_unlock(&zmumps_io_pwrite_mutex);
  return 0;
}

int zmumps_io_init_pointers_lock(){
  pthread_mutex_init(&zmumps_io_pwrite_mutex,NULL);
  return 0;
}
int zmumps_io_destroy_pointers_lock(){
  pthread_mutex_destroy(&zmumps_io_pwrite_mutex);
  return 0;
}

#endif /*WITH_PFUNC*/
#endif /* _WIN32 && WITHOUT_PTHREAD */

int zmumps_io_read__(void * file,void * loc_addr,size_t size,int local_offset){
  int ret_code;
#ifndef _WIN32
  if(!zmumps_directio_flag){
    ret_code=zmumps_io_read_os_buff__(file,loc_addr, size,local_offset);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=zmumps_io_read_win32__(file,loc_addr, size,local_offset);
  if(ret_code<0){
    return ret_code;
  }
#endif  
  return 0;
}

#ifndef _WIN32
int zmumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  /*  printf("Read with buff %d %d %d\n",(int) size, local_offset,*((int *)file));*/
#ifdef WITH_PFUNC
  ret_code=pread(*(int *)file,loc_addr,size,local_offset);
#else
  lseek(*(int *)file,(long) local_offset,SEEK_SET);
  ret_code=read(*(int *)file,loc_addr,size);
#endif
  if((int) ret_code==-1){
    zmumps_io_build_err_str(errno,-90,"Problem with low level read",error_str,200);
    return -90;
  }
  return 0;
}

#endif

#ifdef _WIN32
int zmumps_io_read_win32__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  fseek(*(FILE **)file,(long) local_offset,SEEK_SET);
  ret_code=fread(loc_addr,zmumps_elementary_data_size,size,*(FILE **)file);
  if(ret_code!=size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif

int zmumps_io_write__(void *file, void *loc_addr, size_t write_size, int where){
  int ret_code;
#ifndef _WIN32
  if(!zmumps_directio_flag){
    ret_code=zmumps_io_write_os_buff__(file,loc_addr, write_size,where);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=zmumps_io_write_win32__(file,loc_addr, write_size,where);
  if(ret_code<0){
    return ret_code;
  }
#endif
  return 0;
}

#ifndef _WIN32
int zmumps_io_write_os_buff__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
#ifdef WITH_PFUNC
  ret_code=pwrite(*(int *)file,loc_addr,write_size,where);
#else
  /*in this case all the I/O's are made by the I/O thread => we don't
    need to protect the file pointer.*/
  lseek(*(int *)file,(long)where,SEEK_SET); 
  ret_code=write(*(int *)file,loc_addr,write_size);
#endif
  if(((int)ret_code==-1)||(ret_code!=write_size)){
    zmumps_io_build_err_str(errno,-90,"Problem with low level write",error_str,200);
    return -90;
  }
  return 0;
}

#endif

#ifdef _WIN32
int zmumps_io_write_win32__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
  fseek(*(FILE **)file,(long)where,SEEK_SET);  
  ret_code=fwrite(loc_addr,zmumps_elementary_data_size, write_size,*(FILE**)file);
  if(ret_code!=write_size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif


