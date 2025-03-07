/* @(#)cround.c	 */
/* Entrees :
X0 : matrice de taille [lX0,ntg] (X0 monotone croissant)
lX0 : nombre de lignes de X0
ntg : nombre de colonnes de X0
Y0 : matrice de taille [lX0,nsg] 
nsg : nombre de colonnes de Y0
X : nouveau vecteur de taille lX (X monotone croissant)
(Si les bornes de X debordent de X0, les valeurs de Y correspondantes seront des NaNs)
lX : nombre de lignes de X
Sortie :
Y : nouveau vecteur interpole de taille lX 
Cette fonction suppose que le vecteur Y est alloue par l'appelant 

Remarque : 
L'utilisateur devra s'assurer que les temps utilises soient bien monotones croissants */

/*#include </usr/include/math.h>*/
#include <math.h>


/* La gestion des Not A Number n'est prise en compte
qu'avec la norme IEEE virgule flottante */
# define _IEEE 1

int cround(double *X0,double *Y0,int lX0,int nsg,int ntg,double *X,int lX,double *Y)

{
 double x,y,NotANum;
 int	i,j,k;
 int m,m0,n0,ix;
#ifdef _WIN32
 unsigned long nan[2]={0xffffffff, 0x7fffffff};
#endif

 /* Construction d'un Not A Number */
#ifdef _WIN32
  NotANum = *( double* )nan;
#else
 NotANum = 0.0/0.0;
#endif
 /* ((dnan *)&(NotANum))->nan_parts.exponent = 0x7ff;*/
/*NotANum = NAN(NotANum);*/
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
       printf("cround : valeurs de X a l'exterieur des bornes de X0! \n");
      }
    /* pour un X[i] donne, recherche des valeurs encadrantes de X0,
       puis affectation a la valeur la plus proche */
    ix = 0;
    j = 1;
    for (i=0;i<lX;i++)
      {
       Y[i+m]=NotANum; /* A priori, les valeurs de Y valent NaN */
       x = fabs(X[i]-X0[ix+n0]);

       for (;j<lX0;j++)
         {
          y = fabs(X[i]-X0[j+n0]);/* Valeur absolue du delta t */
          /* on retient l'indice pour lequel delta t est minimum */
          if (y < x)
           {
            x = y;
            ix = j;
           }
          else
           {
            /* l'indice est trouve, on sort de la boucle sur X0, sans incrementer j */
            break;
           }
	 }

       Y[i+m] = Y0[ix+m0];
      }
   }
  return 0;
}
