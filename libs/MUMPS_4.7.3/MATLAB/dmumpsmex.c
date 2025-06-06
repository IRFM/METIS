#include "mex.h"
#include "dmumps_c.h"

#ifdef MATLAB73
long longone=1;
#else
int longone=1;
#endif

#define double2 double
#define mxREAL2 mxREAL
#define EXTRACT_CMPLX_FROM_C_TO_MATLAB EXTRACT_FROM_C_TO_MATLAB
#define EXTRACT_CMPLX_FROM_MATLAB_TOPTR EXTRACT_FROM_MATLAB_TOPTR

#define SYM        (prhs[0])
#define JOB        (prhs[1])
#define ICNTL      (prhs[2])
#define CNTL       (prhs[3])
#define PERM_IN    (prhs[4])
#define COLSCA     (prhs[5])
#define ROWSCA     (prhs[6])
#define RHS        (prhs[7])
#define VAR_SCHUR  (prhs[8])
#define INST       (prhs[9])
#define REDRHS_IN  (prhs[10])
#define A_IN       (prhs[11])

#define INFO_OUT   (plhs[0])
#define RINFO_OUT  (plhs[1])
#define RHS_OUT    (plhs[2])
#define INST_OUT   (plhs[3])
#define SCHUR_OUT  (plhs[4])
#define REDRHS_OUT (plhs[5])
#define PIVNUL_LIST (plhs[6])
#define PERM_OUT   (plhs[7])
#define UNS_PERM   (plhs[8])
#define ICNTL_OUT  (plhs[9])
#define CNTL_OUT   (plhs[10])


#define MYMALLOC(ptr,l,type)                      \
  if(!(ptr = (type *) malloc(l*sizeof(type)))){   \
    mexErrMsgTxt ("Malloc failed in mumpsmex.c");	\
  }                                               



#define MYFREE(ptr) \
  if(ptr){   \
    free(ptr); \
		ptr = 0;	 \
 }



#define EXTRACT_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)          \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  MYFREE(mumpspointer);			 		         	         \
  if(ptr_matlab[0] != -9999){                                                    \
    MYMALLOC(mumpspointer,length,type);                               		 \
    for(i=0;i<length;i++){                                                       \
      mumpspointer[i] = ptr_matlab[i];                                           \
    }                                                                            \
  }


#define EXTRACT_FROM_MATLAB_TOARR(mxcomponent,mumpsarray,type,length)            \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  if(ptr_matlab[0] != -9999){                                                    \
    for(i=0;i<length;i++){                                                       \
      if(ptr_matlab[i] != -9998){                                                \
        mumpsarray[i] = ptr_matlab[i];                                           \
      }                                                                          \
    }                                                                            \
  }                                                                              \

#define EXTRACT_FROM_MATLAB_TOVAL(mxcomponent,mumpsvalue)                        \
  ptr_matlab = mxGetPr(mxcomponent);                                             \
  if(ptr_matlab[0] != -9999){                                                    \
      mumpsvalue = ptr_matlab[0];                                                \
  }                                                                              \


#ifdef MATLAB73
#define EXTRACT_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)                \
  if(mumpspointer == 0){                                                         \
    mxcomponent = mxCreateDoubleMatrix (longone,longone, mxREAL);                \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    ptr_matlab[0] = -9999;                                                       \
  }else{                                                                         \
    mxcomponent = mxCreateDoubleMatrix (longone,(long) length,mxREAL);           \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    for(i=0;i<length;i++){                                                       \
      ptr_matlab[i]=(double)(mumpspointer)[i];                                   \
    }                                                                            \
  }                                                                              
#else
#define EXTRACT_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)                \
  if(mumpspointer == 0){                                                         \
    mxcomponent = mxCreateDoubleMatrix (longone,longone, mxREAL);                \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    ptr_matlab[0] = -9999;                                                       \
  }else{                                                                         \
    mxcomponent = mxCreateDoubleMatrix (longone, length,mxREAL);                 \
    ptr_matlab = mxGetPr (mxcomponent);                                          \
    for(i=0;i<length;i++){                                                       \
      ptr_matlab[i]=(double)(mumpspointer)[i];                                   \
    }                                                                            \
  }                                                                              
#endif

void DMUMPS_free(DMUMPS_STRUC_C **dmumps_par){
  if(*dmumps_par){
  MYFREE( (*dmumps_par)->irn );
  MYFREE( (*dmumps_par)->jcn  );
  MYFREE( (*dmumps_par)->a );
  MYFREE( (*dmumps_par)->irn_loc );
  MYFREE( (*dmumps_par)->jcn_loc );
  MYFREE( (*dmumps_par)->a_loc );
  MYFREE( (*dmumps_par)->eltptr );
  MYFREE( (*dmumps_par)->eltvar );
  MYFREE( (*dmumps_par)->a_elt );
  MYFREE( (*dmumps_par)->perm_in );
  MYFREE( (*dmumps_par)->colsca );
  MYFREE( (*dmumps_par)->rowsca  );
  MYFREE( (*dmumps_par)->pivnul_list );
  MYFREE( (*dmumps_par)->listvar_schur );
  MYFREE( (*dmumps_par)->sym_perm );
  MYFREE( (*dmumps_par)->uns_perm );
  MYFREE( (*dmumps_par)->irhs_ptr);
  MYFREE( (*dmumps_par)->irhs_sparse);
  MYFREE( (*dmumps_par)->rhs_sparse);
  MYFREE( (*dmumps_par)->rhs);
  MYFREE( (*dmumps_par)->redrhs);
  MYFREE(*dmumps_par);
  }
}

void DMUMPS_alloc(DMUMPS_STRUC_C **dmumps_par){

  MYMALLOC((*dmumps_par),1,DMUMPS_STRUC_C);
  (*dmumps_par)->irn  = NULL;
  (*dmumps_par)->jcn  = NULL;
  (*dmumps_par)->a  = NULL;
  (*dmumps_par)->irn_loc  = NULL;
  (*dmumps_par)->jcn_loc  = NULL;
  (*dmumps_par)->a_loc  = NULL;
  (*dmumps_par)->eltptr  = NULL;
  (*dmumps_par)->eltvar  = NULL;
  (*dmumps_par)->a_elt  = NULL;
  (*dmumps_par)->perm_in  = NULL;
  (*dmumps_par)->colsca  = NULL;
  (*dmumps_par)->rowsca  = NULL;
  (*dmumps_par)->rhs  = NULL;
  (*dmumps_par)->redrhs  = NULL;
  (*dmumps_par)->rhs_sparse = NULL;
  (*dmumps_par)->irhs_sparse = NULL;
  (*dmumps_par)->irhs_ptr = NULL;
  (*dmumps_par)->pivnul_list  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}



void mexFunction(int nlhs, mxArray *plhs[ ],
		 int nrhs, const mxArray *prhs[ ]) { 
  
  int i,j,pos;
  int *ptr_int;
  double *ptr_double;
  double *ptr_matlab;
  int tmp_m,tmp_n;

  /* C pointer for input parameters */
  size_t inst_address;
  int n,m,ne,netrue,inst,job;
  #ifdef MATLAB73
  long *irn_in,*jcn_in;
  #else
  int *irn_in,*jcn_in; 
  #endif  

  /* variable for multiple and sparse rhs */
  int nbrhs,ldrhs,nz_rhs,posrhs;
  #ifdef MATLAB73
  long *irhs_ptr, *irhs_sparse;
  #else
  int *irhs_ptr, *irhs_sparse;
  #endif
  double *rhs_sparse;

  DMUMPS_STRUC_C *dmumps_par;
  int dosolve = 0;
  int donullspace = 0;
  int doanal = 0;
  
  
  EXTRACT_FROM_MATLAB_TOVAL(JOB,job);


  dosolve = (job == 3 || job == 5 || job == 6);
  doanal = (job == 1 || job == 4 || job == 6);

  if(job == -1){
    DMUMPS_alloc(&dmumps_par);
    EXTRACT_FROM_MATLAB_TOVAL(SYM,dmumps_par->sym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc = -1;
  }else{
    EXTRACT_FROM_MATLAB_TOVAL(INST,inst_address);
    ptr_int = (int *) inst_address;
    
    dmumps_par = (DMUMPS_STRUC_C *) ptr_int;
    
    if(job == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      DMUMPS_free(&dmumps_par);
    }else{
      
      /* check of input arguments */
      n = mxGetN(A_IN);
      m = mxGetM(A_IN);
			
      if (!mxIsSparse(A_IN) || n != m ){
	mexErrMsgTxt ("Input matrix must be a sparse square matrix");
      }
      jcn_in = mxGetJc(A_IN);
      ne = jcn_in[n];
      irn_in = mxGetIr(A_IN);
      dmumps_par->n = n;
      
      if(dmumps_par->sym != 0){
	netrue = (n+ne)/2;
      }else{
	netrue = ne;
      }
      if(dmumps_par->nz_alloc < netrue || dmumps_par->nz_alloc >= 2*netrue){  
	MYFREE(dmumps_par->jcn);
	MYFREE(dmumps_par->irn);
	MYFREE(dmumps_par->a);
	MYMALLOC((dmumps_par->jcn),netrue,int);
	MYMALLOC((dmumps_par->irn),netrue,int);
	MYMALLOC((dmumps_par->a),netrue,double2);
	dmumps_par->nz_alloc = netrue;
      }


      if(dmumps_par->sym == 0){
	/* if analysis already performed then we only need to read
	   numerical values
	   Note that we suppose that matlab did not change the internal
	   format of the matrix between the 2 calls */
	if(doanal){ 
	  /* || dmumps_par->info[22] == 0 */
	  for(i=0;i<n;i++){
	    for(j=jcn_in[i];j<jcn_in[i+1];j++){
	      (dmumps_par->jcn)[j] = i+1;
	      (dmumps_par->irn)[j] = irn_in[j]+1;
	    }
	  }
	}
	ptr_matlab = mxGetPr(A_IN);
	for(i=0;i<ne;i++){                                                   
	  (dmumps_par->a)[i] = ptr_matlab[i];
	}
	dmumps_par->nz = ne;
      }else{
	/* in the symmetric case we do not need to check doanal */
	pos = 0;
	ptr_matlab = mxGetPr(A_IN);
	for(i=0;i<n;i++){
	  for(j=jcn_in[i];j<jcn_in[i+1];j++){
	    if(irn_in[j] >= i){
	      if(pos >= netrue){
		mexErrMsgTxt ("Input matrix must be symmetric");}
	      (dmumps_par->jcn)[pos] = i+1;
	      (dmumps_par->irn)[pos] = irn_in[j]+1;
	      (dmumps_par->a)[pos] = ptr_matlab[j];
	      pos++;
	    }
	  }
	}
	dmumps_par->nz = pos;
      }
    

      EXTRACT_FROM_MATLAB_TOVAL(JOB,dmumps_par->job);
      EXTRACT_FROM_MATLAB_TOARR(ICNTL,dmumps_par->icntl,int,40);
      EXTRACT_FROM_MATLAB_TOARR(CNTL,dmumps_par->cntl,double,15);
      EXTRACT_FROM_MATLAB_TOPTR(PERM_IN,(dmumps_par->perm_in),int,n);

      EXTRACT_CMPLX_FROM_MATLAB_TOPTR(COLSCA,(dmumps_par->colsca),double,n);
      EXTRACT_CMPLX_FROM_MATLAB_TOPTR(ROWSCA,(dmumps_par->rowsca),double,n);

      dmumps_par->size_schur = mxGetN(VAR_SCHUR);
      EXTRACT_FROM_MATLAB_TOPTR(VAR_SCHUR,(dmumps_par->listvar_schur),int,dmumps_par->size_schur);
      if(!dmumps_par->listvar_schur) dmumps_par->size_schur = 0;

      ptr_matlab = mxGetPr (RHS);

/*
 * To follow the "spirit" of the matlab/scilab interfaces, treat case of null
 * space separately. In that case, we initialize lrhs and nrhs automatically
 * allocate the space needed, and do not rely on what is provided by the user
 * in component RHS, that is not touched.
 *
 * Note that at the moment the user should not call the solution step combined
 * with the factorization step when he/she sets icntl[25-1] to a non-zero value.
 * Hence we suppose infog[28-1] is available and we can use it.
 * 
 * For users of scilab/matlab, it would still be nice to be able to set ICNTL(25)=-1,
 * and use JOB=6. If we want to make this functionality available, we should
 * call separately job=2 and job=3 even if job=5 or 6 and set nbrhs (and allocate
 * space correctly) between job=2 and job=3 calls to MUMPS.
 *
 */
      if ( dmumps_par->icntl[25-1] == -1 && dmumps_par->infog[28-1] > 0 ) {
          dmumps_par->nrhs=dmumps_par->infog[28-1];
          donullspace = dosolve;
         }
      else if ( dmumps_par->icntl[25-1] > 0 && dmumps_par->icntl[25-1] <= dmumps_par->infog[28-1] ) {
           dmumps_par->nrhs=1;
           donullspace = dosolve;
         }
      else {
	    donullspace=0;
	 }
      if (donullspace) {
        nbrhs=dmumps_par->nrhs; ldrhs=n;
	dmumps_par->lrhs=n;
	MYMALLOC((dmumps_par->rhs),((dmumps_par->n)*(dmumps_par->nrhs)),double2);
         }
      else if((!dosolve) || ptr_matlab[0] == -9999 ) { /* rhs not already provided, or not used */
/*JY:  Case where dosolve is true and ptr_matlab[0]=-9999, this could cause problems:
 *        1/ RHS was not initialized while it should have been
 *        2/ RHS was explicitely initialized to -9999 but is not allocated of the right size
 */
        EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS,(dmumps_par->rhs),double,1);
      }else{
	nbrhs = mxGetN(RHS);
	ldrhs = mxGetM(RHS);
	dmumps_par->nrhs = nbrhs;
	dmumps_par->lrhs = ldrhs;
	if(ldrhs != n){
	  mexErrMsgTxt ("Incompatible number of rows in RHS");
	}
	if (!mxIsSparse(RHS)){ /* full rhs */
	  dmumps_par->icntl[19] = 0;
	  EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS,(dmumps_par->rhs),double,(nbrhs*ldrhs));
	}else{ /* sparse rhs */
	  /* printf("sparse RHS ldrhs = %d nrhs = %d\n",ldrhs,nbrhs); */
	  dmumps_par->icntl[19] = 1;
          irhs_ptr = mxGetJc(RHS);
	  irhs_sparse = mxGetIr(RHS);
	  rhs_sparse = mxGetPr(RHS);

	  nz_rhs = irhs_ptr[nbrhs];

	  MYMALLOC((dmumps_par->irhs_ptr),(nbrhs+1),int);
	  MYMALLOC((dmumps_par->irhs_sparse),nz_rhs,int);
	  MYMALLOC((dmumps_par->rhs_sparse),nz_rhs,double2);
	  MYMALLOC((dmumps_par->rhs),(nbrhs*ldrhs),double2);

	  for(i=0;i<nbrhs;i++){
	    for(j=irhs_ptr[i];j<irhs_ptr[i+1];j++){
	      (dmumps_par->irhs_sparse)[j] = irhs_sparse[j]+1;
	    }
	    (dmumps_par->irhs_ptr)[i] = irhs_ptr[i]+1;
	  }
	  (dmumps_par->irhs_ptr)[nbrhs] = nz_rhs+1;
	  for(i=0;i<nz_rhs;i++){                                                   
	    (dmumps_par->rhs_sparse)[i] = rhs_sparse[i];
	  }
	}
      }

      if(dmumps_par->size_schur > 0){
	MYMALLOC((dmumps_par->schur),((dmumps_par->size_schur)*(dmumps_par->size_schur)),double2);
	dmumps_par->icntl[18] = 1;
      }else{
	dmumps_par->icntl[18] = 0;
      }
       /* Reduced RHS */
       if ( dmumps_par->size_schur > 0 && dosolve ) {
          if ( dmumps_par->icntl[26-1] == 2 ) {
            /* REDRHS is on input */
            tmp_m= mxGetM(REDRHS_IN);
            tmp_n= mxGetN(REDRHS_IN);
            if (tmp_m != dmumps_par->size_schur || tmp_n != dmumps_par->nrhs) {
              mexErrMsgTxt ("bad dimensions for REDRHS in mumpsmex.c");
            }
            EXTRACT_CMPLX_FROM_MATLAB_TOPTR(REDRHS_IN,(dmumps_par->redrhs),double,tmp_m*tmp_n);
            dmumps_par->lredrhs=dmumps_par->size_schur;
          }
          if ( dmumps_par->icntl[26-1] == 1 ) {
            /* REDRHS on output. Must be allocated before the call */
            MYFREE(dmumps_par->redrhs);
            if(!(dmumps_par->redrhs=(double2 *)malloc((dmumps_par->size_schur*dmumps_par->nrhs)*sizeof(double2)))){
              mexErrMsgTxt("malloc redrhs failed in intmumpsc.c");
            }
          }
       }
      dmumps_c(dmumps_par);
    }
  }
  if(nlhs > 0){
    EXTRACT_FROM_C_TO_MATLAB( INFO_OUT  ,(dmumps_par->infog),40);
    EXTRACT_FROM_C_TO_MATLAB( RINFO_OUT ,(dmumps_par->rinfog),20);
    if(dmumps_par->rhs && dosolve){
      /* nbrhs may not have been set (case of null space) */
      nbrhs=dmumps_par->nrhs;
      #ifdef MATLAB73
      RHS_OUT = mxCreateDoubleMatrix ((long) n,(long) nbrhs,mxREAL2);
      #else
       RHS_OUT = mxCreateDoubleMatrix ( n, nbrhs,mxREAL2);
      #endif
      ptr_matlab = mxGetPr (RHS_OUT);
      for(j=0;j<nbrhs;j++){
	posrhs = j*n;
	for(i=0;i<n;i++){
	  ptr_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i];	
	}              
      }
    }else{
      EXTRACT_CMPLX_FROM_C_TO_MATLAB( RHS_OUT,(dmumps_par->rhs),1);
    }

    ptr_int = (int *)dmumps_par;
    inst_address = (size_t) ptr_int;
    EXTRACT_FROM_C_TO_MATLAB( INST_OUT   ,&inst_address,1); 
    EXTRACT_FROM_C_TO_MATLAB( PIVNUL_LIST ,dmumps_par->pivnul_list,dmumps_par->infog[27]);
    EXTRACT_FROM_C_TO_MATLAB( PERM_OUT   ,dmumps_par->sym_perm,n);
    EXTRACT_FROM_C_TO_MATLAB( UNS_PERM   ,dmumps_par->uns_perm,n);
    EXTRACT_FROM_C_TO_MATLAB( ICNTL_OUT   ,dmumps_par->icntl,40);
    EXTRACT_FROM_C_TO_MATLAB( CNTL_OUT   ,dmumps_par->cntl,15);

    if(dmumps_par->size_schur > 0){
      #ifdef MATLAB73
      SCHUR_OUT = mxCreateDoubleMatrix((long) dmumps_par->size_schur,(long) dmumps_par->size_schur,mxREAL2);
      #else
      SCHUR_OUT = mxCreateDoubleMatrix( dmumps_par->size_schur, dmumps_par->size_schur,mxREAL2);
      #endif
      ptr_matlab = mxGetPr (SCHUR_OUT);
      for(i=0;i<dmumps_par->size_schur;i++){
	pos = i*(dmumps_par->size_schur);
	for(j=0;j<dmumps_par->size_schur;j++){
	  ptr_matlab[j+pos] = (dmumps_par->schur)[j+pos];
	}
      }
    }else{
      SCHUR_OUT = mxCreateDoubleMatrix(longone,longone,mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
      ptr_matlab[0] = -9999; 
    }
    /* REDRHS on output */
    if ( dmumps_par->icntl[26-1]==1 && dmumps_par->size_schur > 0 && dosolve ) {
      #ifdef MATLAB73
      REDRHS_OUT = mxCreateDoubleMatrix((long) dmumps_par->size_schur,(long) dmumps_par->nrhs,mxREAL2);
      #else
      REDRHS_OUT = mxCreateDoubleMatrix( dmumps_par->size_schur, dmumps_par->nrhs,mxREAL2);
      #endif
      ptr_matlab = mxGetPr(REDRHS_OUT);
      for(i=0;i<dmumps_par->nrhs*dmumps_par->size_schur;i++){
	  ptr_matlab[i] = ((dmumps_par->redrhs)[i]);
      }
    }else{
      REDRHS_OUT = mxCreateDoubleMatrix(longone,longone,mxREAL2);
      ptr_matlab = mxGetPr (REDRHS_OUT);
      ptr_matlab[0] = -9999; 
    }


    MYFREE(dmumps_par->redrhs);
    MYFREE(dmumps_par->schur);
    MYFREE(dmumps_par->irhs_ptr);
    MYFREE(dmumps_par->irhs_sparse);
    MYFREE(dmumps_par->rhs_sparse);
    MYFREE(dmumps_par->rhs);
  }
}
