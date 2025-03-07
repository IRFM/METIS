#include <mex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

/* definitions pour la bibliotheque NCR */
#define NR_END 1
#define FREE_ARG char*

/* definitions pour le mexfile */
#define	X	prhs[0]  
#define Y	prhs[1]
#define XI	prhs[2] 
#define YI	plhs[0]

/* missing function in visual studio libs */
#ifdef _WIN64

    #define fmax(A,B) ((A)>(B) ? (A):(B))
    #define fmin(A,B) ((A)<(B) ? (A):(B))
    #define fabs(A) ((A)<(0) ? (-A):(A))
#endif

double *dvector(unsigned int nl, unsigned int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)mxCalloc(nh-nl+1+NR_END,sizeof(double));
	if (!v) mexErrMsgTxt("allocation failure in dvector()");
	return v-nl+NR_END;
}

void free_dvector(double *v, unsigned int nl, unsigned int nh)
/* free a double vector allocated with dvector() */
{
	mxFree((FREE_ARG) (v+nl-NR_END));
}

/* debut du mexfile */
/*function d = pchipslopes(x,y,del)
%PCHIPSLOPES  Derivative values for shape-preserving Piecewise Cubic Hermite
% Interpolation.
% d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k)).*/


/* pour mac */
/*void mexFunction (int nargout, Matrix *plhs[] , int nargin , Matrix *prhs[])*/

/* pour station */
void mexFunction(int nargout, mxArray *plhs[], int nargin, const mxArray *prhs[])
{  
	/* declaration des variables */
	double *matx,*maty,*matxi,*matyi,*h;
	double hs,w1,w2,dmin,dmax,s1,s2;
	int n,o,k,l;
	int mx,my,mxi,nx,ny,nxi;
	
/*	test du nombre d'argument en E/S 
	
	if (nargout !=1)
		mexErrMsgTxt("this function required one output parameter!");
		
	if (nargin !=3)		
		mexErrMsgTxt("this function required 3 input parameters !");
	recupere les matrice d'entree 
*/	
	matx=mxGetPr(X);
	mx=mxGetM(X);
	nx=mxGetN(X);
		
	maty=mxGetPr(Y);
	my=mxGetM(Y);
	ny=mxGetN(Y);
	
	matxi=mxGetPr(XI);
	mxi=mxGetM(XI);
	nxi=mxGetN(XI);

/*	
       affichage des dimensions
*/       
/*
	mexPrintf("taille de x = %d %d\n",mx,nx);
	mexPrintf("taille de y = %d %d\n",my,ny);
	mexPrintf("taille de xi = %d %d\n",mxi,nxi);
*/	
/* 
        size of ouput vector 
*/	
	o = mx * nx -1;
/* 
        size of input verctor 
*/	
	n = mx * nx;   
	
/*	mexPrintf("tailles = %d %d\n",n,o); */
			
/*      
	creation de la matrice de sortie 
*/
/*	mexPrintf("allocation matice de sortie\n"); */
	YI=mxCreateDoubleMatrix(1,n,mxREAL);
/*	mexPrintf("apres allocation matice de sortie\n");	*/
	if (YI)
	
	{
		matyi=mxGetPr(YI);
				
/*      
        %  Special case n=2, use linear interpolation. 

	n = length(x);
	if n==2  
	d = repmat(del(1),size(y));
	return
	end 
*/		
	      if (n == 2)
	      { 
/*	        mexPrintf("cas a 2 points\n");	*/
		for(k=0;k<=o;k++)
		{		
		    matyi[k] = matxi[0];
/*	            mexPrintf("yi[%d] = %g\n",k,matyi[k]);	*/
		}
		
	      }
	      else
	      {

/* %  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
%  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

   d = zeros(size(y));
  
   k = find(sign(del(1:n-2)).*sign(del(2:n-1)) > 0);
  
   h = diff(x);
   hs = h(k)+h(k+1);
   w1 = (h(k)+hs)./(3*hs);
   w2 = (hs+h(k+1))./(3*hs);
   dmax = max(abs(del(k)), abs(del(k+1)));
   dmin = min(abs(del(k)), abs(del(k+1)));
   d(k+1) = dmin./conj(w1.*(del(k)./dmax) + w2.*(del(k+1)./dmax));
		
*/
/*		mexPrintf("allocation h\n");*/
                h=dvector((unsigned int)0,(unsigned int) n+1);
/*		mexPrintf("boucle h\n");*/
		l = n - 2;
		for(k=0;k<=l;k++)
		{
		 h[k] =  matx[k+1] - matx[k];
/*		 mexPrintf("h = %g\n",h[k]); */
		}

/*		mexPrintf("boucle calcul\n");*/
/* 
la mise a zero est incluse dans mxCreateDoubleMatrix
		for(k=0;k<=o;k++)
		{
		   mexPrintf("mise a 0\n"); 
		   matyi[k] = 0.0;
		}
*/		
/*		   mexPrintf("%d %d %d\n",k,l,o); */
		for(k=0;k<=o;k++)
		{
		   
		   if (k < l)
		   {
/*		      mexPrintf("test signe %g\n",(matxi[k]*matxi[k+1])); */
		      s1 = matxi[k]*matxi[k+1];
		      if (s1 > 0.0)
		      {
/*		        mexPrintf("calcul pout k = %d\n",k+1); */
			hs = h[k] + h[k+1];
			w1 = (h[k]+hs)/(3*hs);
			w2 = (hs+h[k+1])/(3*hs);
			dmax = fmax(fabs(matxi[k]), fabs(matxi[k+1]));
			dmin = fmin(fabs(matxi[k]), fabs(matxi[k+1]));
			matyi[k+1] = dmin/(w1*(matxi[k]/dmax) + w2*(matxi[k+1]/dmax));	
/*		        mexPrintf("hs = %g, w1 = %g, w2 = %g, dmax = %g, dmin = %g\n",hs,w1,w2,dmax,dmin); */
			
		      }
		   }
		}
/*
%  Slopes at end points.
%  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

   d(1) = ((2*h(1)+h(2))*del(1) - h(1)*del(2))/(h(1)+h(2));
   if sign(d(1)) ~= sign(del(1))
      d(1) = 0;
   elseif (sign(del(1)) ~= sign(del(2))) && (abs(d(1)) > abs(3*del(1)))
      d(1) = 3*del(1);
   end
   d(n) = ((2*h(n-1)+h(n-2))*del(n-1) - h(n-1)*del(n-2))/(h(n-1)+h(n-2));
   if sign(d(n)) ~= sign(del(n-1))
      d(n) = 0;
   elseif (sign(del(n-1)) ~= sign(del(n-2))) && (abs(d(n)) > abs(3*del(n-1)))
      d(n) = 3*del(n-1);
   end
*/
/*		mexPrintf("limite 0\n");*/
		matyi[0] = ((2*h[0]+h[1])*matxi[0] - h[0]*matxi[1])/(h[0]+h[1]); 
		s1 = fabs(matyi[0]);
		s2 = fabs(3*matxi[0]);
/*		mexPrintf("matyi[0] = %g\n",matyi[0]);*/
		if (((matyi[0] > 0.0) && (matxi[0] <= 0.0)) || ((matyi[0] <= 0.0) && (matxi[0] > 0.0)))
		{
		  matyi[0] = 0.0;
/*		  mexPrintf("matyi[0] = %g\n",matyi[0]);*/
		}
		else if ((((matxi[0] > 0.0) && (matxi[1] <= 0.0)) || ((matxi[0] <= 0.0) && (matxi[1] > 0.0))) && (s1 > s2))
		{
		  matyi[0] = 3*matxi[0];
/*		  mexPrintf("matyi[0] = %g\n",matyi[0]);*/
		}
		
/*		mexPrintf("limite n\n"); */
		matyi[o] =  ((2*h[o-1]+h[o-2])*matxi[o-1] - h[o-1]*matxi[o-2])/(h[o-1]+h[o-2]);
		s1 = fabs(matyi[o]);
		s2 = fabs(3*matxi[o-1]);
/*		mexPrintf("matyi[o] = %g\n",matyi[o]);*/
		if (((matyi[o] > 0.0) && (matxi[o-1] <= 0.0)) || ((matyi[o] <= 0.0) && (matxi[o-1] > 0.0)))
		{
		  matyi[o] = 0.0;		  
/*		  mexPrintf("matyi[o] = %g\n",matyi[o]);*/
		}
		else if((((matxi[o-1] > 0.0) && (matxi[o-2] <= 0.0)) || ((matxi[o-1] <= 0.0) && (matxi[o-2] > 0.0))) && (s1 > s2))
		{
		  matyi[o] = 3*matxi[o-1];
/*		  mexPrintf("matyi[o] = %g\n",matyi[o]);*/
		}
/*	      	mexPrintf("liberation h\n"); */
		free_dvector(h,(unsigned int)0,(unsigned int)n);
	      }
  	
	}
	
	else 
		mexErrMsgTxt("unable to allocate memory for output variable");
			
		
		
} 
/* 
fin de la fonction mex 
*/
