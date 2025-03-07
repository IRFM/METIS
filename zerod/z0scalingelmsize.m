% calcul de la taille des elms a partir des donnees de la figure 11 de la reference :
% A. Loarte et al , PPCF  45 (2003) p 1549-1569
function dwelmswped = z0scalingelmsize(te95,ti95,dens95,q95,R95,a95,zeff95,zi)

% calcul de nustare
[nustare,nustari,tei] = z0fnustar(te95,ti95,dens95,q95,R95,a95,zeff95,zi);
% donnees pour le fit (juste pour avoir la bonne dependance)
nu_star  = [eps,0.01,0.3,20,1/eps];  
elm_size = [1,0.25,0.13,0.025,0]; 

% fit 
dwelmswped = interp1(log(nu_star),elm_size,log(nustare),'pchip','extrap');


%  figure(22);clf
%  semilogx(nustare,dwelmswped,'b',nu_star,elm_size,'or');
%  set(gca,'Xlim',[min(0.01,min(nustare)),max(20,max(nustare))]);


