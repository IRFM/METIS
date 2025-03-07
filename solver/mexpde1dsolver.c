#include "mex.h"
#include <stdio.h>
/*
 pde1solver(f_in,fp_in,dfpdx_in,dfpdx2_in,        &
                      mata_in,matb_in,matc_in,matd_in,      &
                      matap_in,matbp_in,matcp_in,matdp_in,  &
                      dimk_in,dimn_in,dimm_in,dx_in,dt_in,sca_f_in, &
                      T0_in,T1_in,V0_in,V1_in,mode_in)
*/
extern void pde1solver_(double *f_in,double *fp,
		  double *mata_in,double *matb_in,double *matc_in,double *matd_in,    
		  double *matap_in,double *matbp_in,double *matcp_in,double *matdp_in,
		  int *dimk_in,int *dimn_in,int *dimm_in,double *dx_in,double *dt_in,double *sca_f_in,
		  double *T0_in,double *T1_in,double *V0_in,double *V1_in,double *mode_in,
                  double *dfpdx_in,double *dfpdx2_in);

void mexFunction(int nlhs,       mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  /* Declare variables */ 
  mwSize number_of_dims ;
  const mwSize  *dim_array;
  
  mwSize dimk_in,dimm_in,dimn_in;
  int dimk_4in,dimm_4in,dimn_4in;
  double *dt_in,*dx_in,*sca_f_in;
    
  double *f_in,*fp_in,*fp;
  double *dfpdx_in,*dfpdx2_in,*matd_in,*matdp_in;
  double *mata_in,*matb_in,*matc_in,*matap_in,*matbp_in,*matcp_in;
  double *T0_in,*T1_in;
  double *V0_in,*V1_in;
  double *mode_in;
  int i,j,k;
  double dfpdx_aux=1, dfpdx2_aux=1;
 
  mata_in    =mxGetPr(prhs[0]);
  matb_in    =mxGetPr(prhs[1]);
  matc_in    =mxGetPr(prhs[2]);
  matd_in    =mxGetPr(prhs[3]);
  matap_in   =mxGetPr(prhs[4]);
  matbp_in   =mxGetPr(prhs[5]);
  matcp_in   =mxGetPr(prhs[6]);
  matdp_in   =mxGetPr(prhs[7]);
  f_in       =mxGetPr(prhs[8]);
  fp_in      =mxGetPr(prhs[9]);
  V0_in      =mxGetPr(prhs[10]);
  T0_in      =mxGetPr(prhs[11]);
  V1_in      =mxGetPr(prhs[12]);
  T1_in      =mxGetPr(prhs[13]);
  mode_in    =mxGetPr(prhs[14]);
  sca_f_in   =mxGetPr(prhs[15]);
  dx_in      =mxGetPr(prhs[16]);
  dt_in      =mxGetPr(prhs[17]);

  number_of_dims=mxGetNumberOfDimensions(prhs[0]);
  dim_array=mxGetDimensions(prhs[0]);
  dimk_in=dim_array[0];
  dimm_in=dim_array[1];
    
  dimn_in=0;
  for(i=0;i<dimm_in;i++)
    {
      if (mode_in[i] == 0) 
	{
	  dimn_in=dimn_in+1;
	}
    }

  plhs[0]=mxCreateDoubleMatrix(dimk_in,dimm_in,mxREAL);
  fp=mxGetPr(plhs[0]);
  for(i=0;i<dimk_in*dimm_in;i++)
    {
      fp[i]=fp_in[i];
    }

  if (nlhs > 1) 
   {
    plhs[1]=mxCreateDoubleMatrix(dimk_in,dimm_in,mxREAL);
    dfpdx_in=mxGetPr(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(dimk_in,dimm_in,mxREAL);
    dfpdx2_in=mxGetPr(plhs[2]);
   } else {
     dfpdx_in  = &dfpdx_aux;
     dfpdx2_in = &dfpdx2_aux;
   }

  dimk_4in= (int) dimk_in; 
  dimm_4in= (int) dimm_in; 
  dimn_4in= (int) dimn_in;
 
  /*printf("START MEXPDE1SOLVER\n");
	mexPrintf("START MEXPDE1SOLVER\n");*/
  pde1solver_(f_in,fp,
	      mata_in,matb_in,matc_in,matd_in,    
	      matap_in,matbp_in,matcp_in,matdp_in,
	      &dimk_4in,&dimn_4in,&dimm_4in,dx_in,dt_in,sca_f_in,
	      T0_in,T1_in,V0_in,V1_in,mode_in,dfpdx_in,dfpdx2_in);
  /*printf("%f %f %f %f \n",fp[0],fp[1],fp[2],fp[3]) ;
  printf("%f %f\n",fp[dimk_in],fp[dimk_in+1]) ;
  printf("MEXPDE1SOLVER: OK\n");
	mexPrintf(" MEXPDE1SOLVER OK\n");*/
}
