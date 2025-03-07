#include <mex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/* definitions pour la bibliotheque NCR */
#define NR_END 1
#define FREE_ARG char*

/* definitions pour le mexfile */
#define	X	prhs[0]  
#define Y	prhs[1]
#define XI	prhs[2] 
#define YI	plhs[0]


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

void spline( double x[], double y[], int n, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

	u=dvector((unsigned int)1,(unsigned int) (n-1));
	
	y2[1]=u[1]=0.0;

	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	qn=un=0.0;

	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector(u,(unsigned int)1,(unsigned int)(n-1));
}

void splint( double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0)
	{ 
		mexPrintf("Bad xa input to routine splint");
		*y=(ya[khi]+ya[klo])/2;
	}
	else
	{
		a=(xa[khi]-x)/h;
		b=(x-xa[klo])/h;
		*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	}
}


/* debut du mexfile */

/* pour mac */
/*void mexFunction (int nargout, Matrix *plhs[] , int nargin , Matrix *prhs[])*/

/* pour station */
void mexFunction(int nargout, mxArray *plhs[], int nargin, const mxArray *prhs[])
{  
	/* declaration des variables */
	double *matx,*maty,*matxi,*matyi,*x,*y,*y2;
	double xs,ys;
	int m,n,o,ll,k,adr;
	int mx,my,mxi,nx,ny,nxi;
	
	/* test du nombre d'argument en E/S */
	
	if (nargout !=1)
		mexErrMsgTxt("Il faut un argument en sortie !");
		
	if (nargin !=3)		
		mexErrMsgTxt("Il faut 3 arguments en entree !");
		
	/* recupere les matrice d'entree */
	
	matx=mxGetPr(X);
	mx=mxGetM(X);
	nx=mxGetN(X);
		
	maty=mxGetPr(Y);
	my=mxGetM(Y);
	ny=mxGetN(Y);
	
	matxi=mxGetPr(XI);
	mxi=mxGetM(XI);
	nxi=mxGetN(XI);
	
	/* affichage des dimensions
	mexPrintf("taille de x = %d %d",mx,nx);
	mexPrintf("taille de y = %d %d",my,ny);
	mexPrintf("taille de xi = %d %d",mxi,nxi); */
	
	
	/* verification des dimension des matrices */
	
	if ((mx!=my)||(nx!=ny))
	{
		mexPrintf("les deux premieres matrices doivent avoir memes dimensions !");
		YI=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}	
	if (mx!=mxi)
	{
		mexPrintf("toutes les matrices doivent avoir la meme premiere dimension!");
		YI=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	/* assignation des varibles */
	
	o=mx-1;
	n=nx;
	m=nxi;
			
	/* creation de la matrice de sortie */
	
	YI=mxCreateDoubleMatrix(mx,nxi,mxREAL);
	
	if (YI)
	
	{
		matyi=mxGetPr(YI);
				
		y2=dvector((unsigned int)0,(unsigned int)(n+1));
		x=dvector((unsigned int)0,(unsigned int)(n+1));
		y=dvector((unsigned int)0,(unsigned int)(n+1));
	
		
		for ( ll=0 ; ll<=o ; ll++)
		{
	
			for(k=1;k<=n;k++)
			{		
				adr=(k-1)*(o+1)+ll;
				x[k]=matx[adr];
				y[k]=maty[adr];	
				
				/* control des entrees
				mexPrintf("%g	%g\n", x[k],y[k]); */
					
			}

			spline(x,y,n,y2);
	
			for(k=1;k<=m;k++)
			{					
				adr=(k-1)*(o+1)+ll;	
				xs=matxi[adr];
				splint( x ,y , y2 , n , xs , &ys);
				matyi[adr]=ys;
			
				/* control des sorties		
				mexPrintf("%g	%g\n", xs , ys); */
			}

		}
		free_dvector(y2,(unsigned int)0,(unsigned int)(n+1));
		free_dvector(y,(unsigned int)0,(unsigned int)(n+1));
		free_dvector(x,(unsigned int)0,(unsigned int)(n+1));
	
	}
	
	else 
		mexErrMsgTxt("erreur d'allocation memoire pour la matrice de sortie");
			
		
		
} /* fin de la fonction mex */
