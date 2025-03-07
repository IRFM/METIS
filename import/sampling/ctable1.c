/* @(#)ctable1.c	 */
/* Entrees :
X0 : matrice de taille [lX0,ntg] (X0 monotone croissant)
lX0 : nombre de lignes de X0
ntg : nombre de colonnes de X0
Y0 : matrice de taille [lX0,nsg] 
nsg : nombre de colonnes de Y0
X : nouveau vecteur de taille lX (X monotone croissant )
(Si les bornes de X debordent de X0, les valeurs de Y correspondantes seront des NaNs)
lX : nombre de lignes de X
Sortie :
Y : matrice interpolee de taille [lX,nsg] 
Cette fonction suppose que le vecteur Y est alloue par l'appelant 

Remarque : 
L'utilisateur devra s'assurer que les temps utilises soient bien monotones croissants */

/* La gestion des Not A Number n'est prise en compte
qu'avec la norme IEEE virgule flottante */
#include <mex.h>
#include <stdio.h>
# define _IEEE 1
#include <math.h>

#ifdef _MSC_VER
#include <float.h>  /* for _isnan() on VC++*/
#define isnan(x) _isnan(x)  // VC++ uses _isnan() instead of isnan()
#else
#include <math.h>  /* for isnan() everywhere else*/
#endif

int ctable1(double *X0,double *Y0,int lX0,int nsg,int ntg,double *X,int lX,double *Y)

{
 int	i,j,k;
 int m,m0,n0;
 double NotANum;
 
#ifdef _MSC_VER
 unsigned long nan[2]={0xffffffff, 0x7fffffff};
#endif
         
 /* Construction d'un Not A Number */
#ifdef _MSC_VER
  NotANum = *( double* )nan;
#else
 NotANum =  0.0/0.0; /*1e1*/;
#endif
/* ((dnan *)&(NotANum))->nan_parts.exponent = 0x7ff;*/

 n0=0;
 /* si la donnee est un groupe, chaque signal est traite separement */
 for (k=0;k<nsg;k++)
   {
    m = k * lX;    /* deplacement dans Y */
    m0 = k * lX0;  /* deplacement dans Y0 */
    /* si chaque signal du groupe a son propre vecteur temps,
       deplacement dans X0 = deplacement dans Y0 */
    if (ntg==nsg) n0=m0;
    /* Verification de l'inclusion du vecteur X dans le vecteur X0 
    if (X[0] < X0[0+n0] || X[lX-1] > X0[lX0-1+n0])
      {
       printf("ctable1 : valeurs du temps interpolant hors bornes du temps de la donnee ! \n");
      }*/
    /* Detection de Not a Number dans les donnees */
    for (j=0;j<lX0;j++)
      {
       if (isnan(Y0[j+m0])) 
         {
          mexPrintf("Not a Number detecte dans les mesures \n");
          return -1;
         }
      }
    /* pour un X[i] donne, recherche des valeurs encadrantes de X0,
       puis interpolation lineaire */
    j=0;
    for (i=0;i<lX;i++)
      {
       if (X[i] < X0[n0] || X[i] > X0[lX0-1+n0])
         {
          Y[i+m]=NotANum; /* Les valeurs de Y valent NaN pour les temps interpolant hors bornes*/
         }
       else
         {  
          for (;j<lX0-1;j++)
            {
	     if (X[i] >= X0[j+n0] && X[i] <= X0[j+1+n0]) 
               {
                Y[i+m] = Y0[j+m0]+(Y0[j+1+m0]-Y0[j+m0])*(X[i]-X0[j+n0])/(X0[j+1+n0]-X0[j+n0]);
               /* l'indice est trouve, on sort de la boucle sur X0, sans incrementer j */
	        break;
	       }
	    }
         }
      }
   }
  return 0;
}
