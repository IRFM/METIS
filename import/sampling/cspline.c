/* @(#)cspline.c	 */
/* Entrees :
x : vecteur ou matrice des temps x(temps) ou x(temps,espace)
y : matrice des valeurs y(temps,espace)
xx : vecteur de reechantillonage des donnees xx(temps_nouveau)
Sortie :
yy : matrice de sortie 

dimensions(x)=[nt_entree,1] ou [nt_entree,nsig]
dimensions(y)=[nt_entree,nsig]
dismensions(xx)=[nt_sortie,1]
dimensions(yy)=[nt_sortie,nsig];

si x est un vecteur, ncolx =1
si x est une matrice, ncolx=nsig

Cette fonction suppose que la matrice yy est allouee par l'appelant */

/* La gestion des Not A Number n'est prise en compte
qu'avec la norme IEEE virgule flottante */

# define _IEEE 1
#include <mex.h>
#include <math.h>

#ifdef _MSC_VER
#include <float.h>  /* for _isnan() on VC++*/
#define isnan(x) _isnan(x)  // VC++ uses _isnan() instead of isnan()
#else
#include <math.h>  /* for isnan() everywhere else*/
#endif

/*
 * These define statements allow the Matlab arrays to be indexed correctly.
 * Matlab
 indexes into its arrays column-wise, not row-wise like C. Also, C
 * is zero-based, 
 not one-based like Matlab. These defined macros perform
 * the necessary accessing
 of the matlab arrays. 
 */

#define XX(i,j) ( *(xx+((int)(j)-1)*nt_sortie+(int)(i)-1) )
#define X(i,j) ( *(x+((int)(j)-1)*nt_entree+(int)(i)-1) )
#define Y(i,j) ( *(y+((int)(j)-1)*nt_entree+(int)(i)-1) )
#define YY(i,j) ( *(yy+((int)(j)-1)*nt_sortie+(int)(i)-1) )
#define Y2(i,j) ( *(y2Pr+((int)(j)-1)*y2M+(int)(i)-1) )
#define U(i,j) ( *(uPr+((int)(j)-1)*uM+(int)(i)-1) )



/*#ifdef _WIN32
#define isnan(x) ((x) != (x))
#endif*/

int cspline(double *x,double *y,int nt_entree,int nsig,int ncolx,double *xx,int nt_sortie,double *yy)

{

 /* Scalar variables */
 double          a, b, sig, h, uu, p, xc;

 /* Matrix Variables */
 mxArray         *y2, *u;
 double         *y2Pr, *uPr;
 int             y2M, uM;
#ifdef _MSC_VER
 unsigned long nan[2]={0xffffffff, 0x7fffffff};
#endif
 /* varibale de boucle */
 int	j,k,l,m,n,klo,khi;
 
 /* Not a Number */
 double NotANum;
 
#ifdef _MSC_VER
  NotANum = *( double* )nan;
#else
 NotANum = 0.0/0.0;
#endif
/*((dnan *)&(NotANum))->nan_parts.exponent = 0x7ff;*/
/*NotANum = NAN(NotANum);*/
/* Allocate Memory for y2 */
 y2 = mxCreateDoubleMatrix((int) nt_entree, (int) 1, mxREAL);
 y2Pr = mxGetPr(y2);
 y2M = mxGetM(y2);

 /* Allocate Memory for u */
 u = mxCreateDoubleMatrix((int) nt_entree, (int) 1, mxREAL);
 uPr = mxGetPr(u);
 uM = mxGetM(u);

 /* affichage de test */
 /* printf("nsig = %8i\n",nsig); */
 /* printf("ncolx = %8i\n",ncolx); */
 /* printf("nt_entree = %8i\n",nt_entree); */
 /* printf("nt_sortie = %8i\n",nt_sortie); */
  
  
  
  /* debut du calcul */
  
  
   /* boucle sur les ligne */
   for (k = 1; k <= nsig; k ++) {


     /* Detection de Not a Number dans les donnees */
     for (j=0;j<nt_entree;j++)
       {
        if (isnan(Y(j,k))) 
          {
           mexPrintf("Not a Number detecte dans les mesures \n");
           mxDestroyArray(y2);
           mxDestroyArray(u);
           return -1;
          }
       }

	  /* selon la dimension de x*/
	  if (ncolx > 1){
	  	n=k;
	  } else {
	  	n=1;
	  }
	  
      /* condition en 1 */
      Y2(1, 1) = 0.0;
      U(1, 1) = 0.0;

      /* calcul de u et y2 */
       for (l = 2; l <= (nt_entree - 1); l ++) {

         sig = (X(l, n) - X(l - 1, n)) / (X(l + 1, n) - X(l - 1, n));
         p = sig * Y2(l - 1, 1) + 2.0;
         Y2(l, 1) = (sig - 1.0) / p;
         uu = (Y(l + 1, k) - Y(l, k)) / (X(l + 1, n) - X(l, n));
         uu = uu - (Y(l, k) - Y(l - 1, k)) / (X(l, n) - X(l - 1, n));
         U(l, 1) = (6.0 * uu / (X(l + 1, n) - X(l - 1, n)) - sig * U(l - 1, 1)) / p;

      }

      /* condition en nt_entree */
      Y2(nt_entree, 1) = 0.0;


      /* calcul final de y2 */
      for (l = (nt_entree - 1); l >= 1; l -- ) {

         Y2(l, 1) = Y2(l, 1) * Y2(l + 1, 1) + U(l, 1);
      }

      /* calcul de yy */
      for (l = 1; l <= nt_sortie; l ++ ) {

         /* recherche de l'intervalle */
         klo = 1;
         khi = nt_entree;
         xc = XX(l,1);

         /* dans l'intervalle */
         if ((xc >= X(1,n)) && (xc <= X(nt_entree,n))) {

            /* boucle de recherche */
            while ((khi - klo) > 1) {
               m = ((khi + klo) >> 1);

               if (X(m, n) > xc) {
                  khi = m;

               } else {
                  klo = m;
               }
            }

            /* calcul de yy(l) */
            h = X(khi, n) - X(klo, n);

            if (h == 0) {
               YY(l, k) = NotANum;

            } else {
               a = (X(khi, n) - xc) / h;
               b = (xc - X(klo, n)) / h;
               uu =((a * a * a - a) * Y2(klo, 1) + (b * b * b - b) * Y2(khi, 1));
               YY(l, k) = a * Y(klo, k) + b * Y(khi, k)+ uu * (h * h) / 6.0 ;
            }

         } else {
 
             /* hors intervalle */
            YY(l, k) = NotANum;
         }
 
       /* fin boucle sur les temps */
       }
 
   /* fin de la boucle sur les signaux */
   }

   /* liberation de la memoire */
 
   mxDestroyArray(y2);
   mxDestroyArray(u);
 
  
  /* fin de la fonction */
  return 0;
}
