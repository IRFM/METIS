/* TSAMPLE	Reechantillonnage de donnees selon un nouveau temps.
[A,B,C,D,E,F,G,H,I,J,K,L,M,N] = tsample(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)

 [r1,r2,...,t] = tsample(x1,t1,x2,t2,...);
 [r1,r2,...,t] = tsample(x1,t1,x2,t2,...[,T|t][,mode]);

 r1,... donnees reechantillonnees, un signal par colonne
 t      temps reechantillonne

 x1,... signal ou groupe (un signal par colonne) complete eventuellement de NaN's
 t1,... temps correspondant

Options :
  contrainte temps pour le reechantillonnage 
  T      periode de reechantillonnage desiree 
   ou
  t      temps desire
 mode de reechantillonnage 'abc' ('fel' est le defaut)
  a :	Contrainte periode
      f  selon la periode la plus rapide (fast)
      s  selon la periode la plus lente (slow)
  b :  Contrainte vecteur temps
      r  restriction au vecteur temps inclus dans tous
      e  extension au plus grand des vecteurs temps 
       (les mesures sont completees eventuellement de NaN's)
  c :  Type reechantillonnage
      l  interpolation lineaire
      c  interpolation par spline cubique
      n  valeur du plus proche echantillon

Fonctionnement et restrictions :
1) Si des valeurs du temps de reechantillonnage sont a l'exterieur du temps a interpoler,
%les mesures correspondantes seront des "Not A Number" et il apparaitra le message suivant :
% "valeurs du temps interpolant hors bornes du temps de la donnee"
2) Cette routine a pour objet d'etre plus rapide que son correspondant en Matlab.
Elle est prevue uniquement pour une abscisse temps, qui par definition est monotone croissant.
3) Les donnees a reechantillonner ne doivent pas contenir de Not A Number.
Dans ce cas, la routine s'interrompt et affiche "Utilisation incorrecte de TSAMPLE".

Routines appelees : cround, ctable1 et cspline

 R.Masset  - CEA DRFC Tore Supra - Juin 1995 
 JF.Artaud - CEA DRFC Tore Supra - Juillet 1995 */

#include <mex.h>
#include <stdio.h>
#ifdef _WIN32
  #include <string.h>
  #include <float.h>
#else
  #include <strings.h>
#endif
#include <memory.h>
#include <math.h>
#include <limits.h>


/* plus petit nombre de type "double" */
#define	DBL_MIN		(-DBL_MAX) 
#define max(A,B) ((A)>(B) ? (A):(B))
#define min(A,B) ((A)<(B) ? (A):(B))

/* routine mexFunction obligatoire pour creer le mex-file : tsample.mex.
nrhs arguments d'entree et nlhs arguments de sortie.
prhs : tableau des nrhs pointeurs des arguments d'entree.
plhs : tableau des nlhs pointeurs des arguments de sortie.
*/

extern int ctable1(double *,double *,int ,int ,int ,double *,int ,double *);
extern int cround(double *,double *,int ,int ,int ,double *,int ,double *);
extern int cspline(double *,double *,int ,int ,int ,double *,int ,double *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
double *tabin[26],*tabout[14],*time;
unsigned int m,n;
double tf,tl,T,Tf,Ts,tpf,tpl;
int nbmes,nbt,icr,npair,reste,ntp[14],nSigGrp[14],nTpsGrp[14];
int i,j,k,cr;
char mode[4];

if (nlhs==0) mexErrMsgTxt("Il faut au moins 1 argument de sortie !");
/* Assignation des pointeurs aux variables d'entree */
for (i=0;i<nrhs;i++)     tabin[i]=mxGetPr(prhs[i]);

/* Nombre de paires d'arguments d'entree */
npair = nrhs/2 ;
/*printf("nombre de paires = %d \n",npair);*/
/* Reste de la division par 2 */
reste = nrhs%2 ;

/* si le dernier argument en entree est une chaine de caracteres */
if (mxIsChar(prhs[nrhs-1])) 
  {
   /*printf("dernier argument = chaine de caracteres \n");*/
   mxGetString(prhs[nrhs-1],mode,4);
   mode[3]='\0';
   /*printf("mode = %s \n",mode);*/
   /* si le nombre d'arguments d'entree est impair */
   if (reste == 1) 
     {   
      printf("pas de contrainte temps pour le reechantillonnage demandee \n");
     }
   /* si le nombre d'arguments d'entree est pair */
   else
     {
      npair = npair - 1;
      /*printf("nombre de paires = %d \n",npair);*/
      /* la contrainte temps, demandee pour le reechantillonnage, est l'avant dernier argument */
      nbt = mxGetM(prhs[nrhs-2])*mxGetN(prhs[nrhs-2]);
      /* si la contrainte temps est une periode 
         une periode T est imposee ( mode(0)='i' )*/
      if (nbt==1) 
        {
         mode[0]='i';
         T = tabin[nrhs-2][0];
        }
      /*  i la contrainte temps est un vecteur temps */
      else
        {
         mode[1]='i';
         time=tabin[nrhs-2];
        }
     }
  }
 
/* si le dernier argument en entree n'est pas une chaine de caracteres */
else
  {
   /* si le nombre d'arguments d'entree est impair,
      seule une contrainte temps pour le reechantillonnage est demandee */
   if (reste == 1) 
     {   
      /* la contrainte temps, demandee pour le reechantillonnage, est le dernier argument */
      nbt = mxGetM(prhs[nrhs-1])*mxGetN(prhs[nrhs-1]);
      /* si la contrainte temps est une periode  
         une periode T est imposee ( mode(0)='i' )*/
      if (nbt==1) 
        {
         strcpy(mode,"iel");
         T = tabin[nrhs-1][0];
        }
      /* si la contrainte temps est un vecteur temps */
      else  
        {
         strcpy(mode,"fil");
         time = tabin [nrhs-1];
        }   
     }   
   /* si le nombre d'arguments d'entree est pair,
      pas de contrainte temps pour le reechantillonnage demandee  */  
   else    
   /* mode par defaut : 
      periode la plus rapide, extension au plus grand temps, interpolation lineaire */
     {   
      strcpy(mode,"fel");
     }   
  }

/* Evaluation des caracterististiques temporelles de chaque signal */
 Tf = DBL_MAX; Ts =DBL_MIN;

/* pour chaque temps de chaque paire (signal,temps) */
for (i=0;i<npair;i++)
  {
   /* Nombre de signaux du groupe */
   nSigGrp[i] = mxGetN(prhs[2*i]);
   k = 2*i + 1;
   /* Nombre de vecteurs temps du groupe */
   nTpsGrp[i] = mxGetN(prhs[k]);
   /* Nombre de temps du groupe */
   ntp[i] = mxGetM(prhs[k]);
   /* Le nombre de mesures doit etre egal au nombre de temps du groupe */
   if (ntp[i] !=  mxGetM(prhs[2*i])) mexErrMsgTxt("Pour une donnee, le nombre de temps est different du nombre de mesures");
   /* Le nombre de temps doit etre au moins egal a 2 */
   if (ntp[i] < 2) mexErrMsgTxt("Pour une donnee, le nombre de temps est inferieur a 2 ");
   tpf=tf; tpl=tl;
   tf = tabin[k][0]; tl = tabin[k][0];
   for (j=1;j<ntp[i];j++)
     {
      tf = min(tf,tabin[k][j]);
      tl = max(tl,tabin[k][j]);
      Tf = min(Tf,(tabin[k][j]-tabin[k][j-1]));
      Ts = max(Ts,(tabin[k][j]-tabin[k][j-1]));
     }
   if (mode[1] == 'r' && i>0)
     {
      tf = max(tf,tpf);
      tl = min(tl,tpl);
     }
   else if (mode[1] == 'e' && i>0)
     {
      tf = min(tf,tpf);
      tl = max(tl,tpf);
     }
  }

/* choix de la periode rapide ou lente */
if (mode[0] == 'f')   T = Tf;
else if (mode[0] == 's')  T = Ts;

/* Constitution du vecteur temps pour interpolation, 
   si celui-ci n'est pas impose */
if (mode[1] != 'i')
  {
   nbt = 1 + (tl - tf)/T;
   /*printf("nbt = %d \n",nbt);*/
   time = (double *) mxCalloc(nbt,sizeof(double));
   for (i=0;i<nbt;i++) time[i] = tf + i * T;
  }

/* le vecteur temps pour interpolation a une seule colonne */
nSigGrp[npair]=1;
 
/* Creation des matrices de sortie et Assignation des pointeurs aux variables */
for (i=0;i<min(npair,nlhs);i++)
  {
   plhs[i] = mxCreateDoubleMatrix(nbt,nSigGrp[i],mxREAL);
   tabout[i] = mxGetPr(plhs[i]);
  }
if (nlhs>npair)
  {
   plhs[npair] = mxCreateDoubleMatrix(nbt,1,mxREAL);
   tabout[npair] = mxGetPr(plhs[npair]);
   for (i=0;i<nbt;i++) tabout[npair][i] = time[i];
  }
   
/* interpolation
   pour chaque paire (signal, temps correspondant) */
for (i=0;i<min(npair,nlhs);i++)
  {
   /* interpolation lineaire --> resultats A,B,C,D,.....*/
   if (mode[2] == 'l') 
     {
      cr = ctable1(tabin[2*i+1],tabin[2*i],ntp[i],nSigGrp[i],nTpsGrp[i],time,nbt,tabout[i]);
      if (cr==-1) mexErrMsgTxt("Utilisation incorrecte de TSAMPLE");
     }

   /* arrondi a la valeur la plus proche --> resultats A,B,C,D,.....*/
   else if (mode[2] == 'n')
     {
      cr = cround(tabin[2*i+1],tabin[2*i],ntp[i],nSigGrp[i],nTpsGrp[i],time,nbt,tabout[i]);
     }

   /* interpolation par spline cubique  --> resultats A,B,C,D,.....*/
   else if (mode[2] == 'c')
     {
      cr = cspline(tabin[2*i+1],tabin[2*i],ntp[i],nSigGrp[i],nTpsGrp[i],time,nbt,tabout[i]);
      if (cr==-1) mexErrMsgTxt("Utilisation incorrecte de TSAMPLE");
     }
  }
}
