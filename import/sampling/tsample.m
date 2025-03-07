% TSAMPLE	Reechantillonnage de donnees selon un nouveau temps.
%[A,B,C,D,E,F,G,H,I,J,K,L,M,N] = tsample(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
%
% [r1,r2,...,t] = tsample(x1,t1,x2,t2,...);
% [r1,r2,...,t] = tsample(x1,t1,x2,t2,...[,T|t][,mode]);
%
% r1,... donnees reechantillonnees, un signal par colonne
% t      temps reechantillonne
%
% x1,... signal ou groupe (un signal par colonne) complete eventuellement de NaN's
% t1,... temps correspondant
%
% Options :
% contrainte temps pour le reechantillonnage 
%  T      periode de reechantillonnage desiree 
%   ou
%  t      temps desire
% mode de reechantillonnage 'abc' ('fel' est le defaut)
%  a :	Contrainte periode
%      f  selon la periode la plus rapide (fast)
%      s  selon la periode la plus lente (slow)
%  b :  Contrainte vecteur temps
%      r  restriction au vecteur temps inclus dans tous
%      e  extension au plus grand des vecteurs temps 
%       (les mesures sont completees eventuellement de NaN's)
%  c :  Type reechantillonnage
%      l  interpolation lineaire
%      c  interpolation par spline cubique
%      n  valeur du plus proche echantillon
%
%Fonctionnement et restrictions :
%1) Si des valeurs du temps de reechantillonnage sont a l'exterieur du temps a interpoler,
%les mesures correspondantes seront des "Not A Number" et il apparaitra le message suivant :
% "valeurs du temps interpolant hors bornes du temps de la donnee"
%2) Cette routine a pour objet d'etre plus rapide que son correspondant en Matlab.
%Elle est prevue uniquement pour une abscisse temps, qui par definition est monotone croissant.
%3) Les donnees a reechantillonner ne doivent pas contenir de Not A Number.
%Dans ce cas, la routine s'interrompt et affiche "Utilisation incorrecte de TSAMPLE".
%
% Routines appelees : cround, ctable1 et cspline

% R.Masset - CEA DRFC Tore Supra  - Juin 1995
% JF.Artaud - CEA DRFC Tore Supra - Juillet 1995 

