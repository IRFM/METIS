% function AAleg = equicoronal_leger(iz,te,display_fig) 
%
%	M. Mattioli
%	Calcul des abondances fractionnelles d'une espèce légère (jusqu'à Z=10)
%	dans le cadre du modèle coronal.
%
%***********  /mattioli/mac/equicorleg.m
%*****
%*****  clear,  iz=10,  equicorleg, 
%
%	Cette routine appelle:
%		sion_leger.m
%		alfarec_leger.m
%
%	Arguments d'entrée:
%		iz = numéro atomique de l'espèce
%		te = vecteur des températures électroniques (eV)
%		display_fig =	1: affichage des abondances fractionnelles
%						0: pas d'affichage
%
%	Argument de sortie:
%			AAleg = abondances fractionnelles, matrice  nbre de Te x nbre d'états d'ionisation
%
%	Modif R. Guirlet le  2/08/2006 : ajout de l'argument display_fig
%	Modif R. Guirlet le 14/09/2006 : ajout de l'aluminium
%   Modif R. Guirlet le  7/01/2011 : tous types de symboles acceptés (ex carbone: C, c, CA, ca, Ca)
%
function AAleg = equicoronal_leger(iz,te,display_fig) 
%
if nargin==2
	display_fig = 0;
end
%clear te
%if iz==10, te=logspace(-0.,3.6,37); fmin=20.;  fmax=30.; end
%if iz==8,  te=logspace(-0.,3.6,37); fmin=20.;  fmax=35.; end
%if iz==7,  te=logspace(-0.,3.2,33); fmin=20.;  fmax=35.; end
%if iz==6,  te=logspace(-0.,3.2,33); fmin=20.;  fmax=40.; end
%if iz==5,  te=logspace(-0.,3.0,31); fmin=15.;  fmax=40.; end
%if iz==4,  te=logspace(-0.,3.0,31); fmin=15.;  fmax=40.; end
if iz==10
	fmin=20.;  fmax=30.; 
elseif iz==8
	fmin=20.;  fmax=35.;
elseif iz==7
	fmin=20.;  fmax=35.;
elseif iz==6
	fmin=20.;  fmax=40.;
elseif iz==5
	fmin=15.;  fmax=40.;
elseif iz==4
	fmin=15.;  fmax=40.; 
elseif iz==2
	disp('Sorry - no coronal computation available for He - Assume He 2+ only')
	AAleg = [zeros(length(te),2) ones(length(te),1)];
	return
end
if iz>13
	disp('Sorry - no coronal computation available for this element')
	AAleg = 0;
	return
end

Amin=1e-30;    X=1.;  NETAT=iz+1;
%clear AAleg
pion4=[9.3 18.2 153.9 217.7];
pion5=[ 8.30 25.15 37.93 259.37 340.22];
%pion6=[ 11.26 24.38 41.38 64.49 392.08 489.98]; erreur potentiel  d'ionisation C III
% d'apres Kelly
pion6=[ 11.26 24.38 47.88 64.49 392.08 489.98];
%pion7=[ 14.53 29.60 47.45 69.13 97.89 552.06 667.03];erreur NIV
pion7=[ 14.53 29.60 47.45 77.5 97.89 552.06 667.03];
%pion8=[ 13.62 35.12 54.93 68.6 103.68 138.12 739.32 871.39]; erreur O IV V
pion8=[ 13.62 35.12 54.93 77.4 113.90 138.12 739.32 871.39];
%pion10=[ 21.6 41.1 63.5 92.5 126.2 157.9 207.5 239.1 1196. 1360.];
pion10=[ 21.6 41.1 63.5 97.1 126.2 157.9 207.5 239.1 1196. 1360.];

if iz==4 | iz==5 | iz==6 | iz==7 | iz==8 | iz==10
%	for k=1:length(te)
%		for i=1:NETAT 
%			AAleg(k,i)=0.; 
%		end 
%	end
	AAleg = zeros(length(te),NETAT);
	pion = eval(['pion',int2str(iz)]);
	% Evaluation des coefficients d'ionisation (cleg) et de recombinaison (rleg)
	cleg = sion_leger(iz,te);
	[rleg,rradleg,rdielleg] = alfarec_leger(iz,te);
	%
	for k=1:length(te) 
		% Définition des états d'ionisation min et max à prendre en compte à cette Te
		mina=1;maxa=2; 
		for i=2:iz
 			if (pion(i)-te(k)/fmin)<=0 
				mina=mina+1; 
			end
 			if (te(k)*fmax-pion(i))>=0 
				maxa=maxa+1; 
			end
		end
		%
		minb=mina+1; 
		P(mina)=1.; 
		SOMME=P(mina); 
		for i=mina:maxa-1    
			i1=i+1; 
			SALF=cleg(i,k)/rleg(i,k);
   			P(i1)=P(i)*SALF; 
			SOMME=SOMME+P(i1);  
		end
		AAleg(k,mina)=X/SOMME;
		for i=minb:maxa
			AAleg(k,i)=P(i)*AAleg(k,mina); 
			if AAleg(k,i)<=Amin
				AAleg(k,i)=Amin; 
			end
		end
		if AAleg(k,mina)<=Amin 
			AAleg(k,mina)=Amin; 
		end
		S(k)=sum(AAleg(k,:));
	end
%S 
%max(S), min(S)
	if iz==10 | iz==7 | iz==5 
		AAlegold=AAleg; 
	end 
%
elseif iz == 3
elseif iz == 9
elseif iz == 11
elseif iz == 12
elseif iz == 13
	log10te_ref = 4:.1:8;		% log10(Te_reference), Te_reference en kelvins
	log10te_ref = log10te_ref(:);	% log10te en colonne
	frac_ab_ref = [-1.997  -0.004 -10     -10	  -10	  -10	  -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -2.183  -0.003 -4.078  -10	  -10	  -10	  -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -2.277  -0.003 -2.788  -10	  -10	  -10	  -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -2.378  -0.010 -1.745  -10	  -10	  -10	  -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -2.533  -0.056 -0.926  -3.393  -10     -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -2.846  -0.246 -0.394  -1.567  -10     -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -3.518  -0.782 -0.323  -0.445  -10     -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -4.671  -1.790 -0.815  -0.081  -10     -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -10     -2.872 -1.454  -0.016  -10     -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -10     -3.821 -2.019  -0.004  -4.121  -10     -10	  -10	  -10	   -10    -10     -10     -10     -10
                   -10     -4.617 -2.477  -0.003  -2.467  -10     -10     -10	  -10	   -10    -10     -10	  -10	  -10
                   -10     -10    -2.819  -0.027  -1.227  -10	  -10     -10	  -10	   -10    -10     -10	  -10	  -10
                   -10     -10    -3.129  -0.186  -0.459  -3.010  -10     -10	  -10      -10 	  -10	  -10	  -10 	  -10
                   -10     -10    -3.539  -0.561  -0.151  -1.728  -4.725  -10	  -10      -10    -10     -10	  -10	  -10
                   -10     -10    -4.031  -1.051  -0.109  -0.880  -2.864  -10	  -10      -10    -10     -10	  -10	  -10
                   -10     -10    -4.649  -1.655  -0.277  -0.377  -1.531  -4.361  -10      -10    -10     -10	  -10	  -10
                   -10     -10    -10     -2.438  -0.690  -0.226  -0.709  -2.684  -10      -10    -10     -10	  -10	  -10
                   -10     -10    -10     -3.373  -1.304  -0.366  -0.308  -1.583  -3.280  -10     -10	  -10	  -10 	  -10
                   -10     -10    -10     -4.378  -2.026  -0.684  -0.192  -0.892  -1.923  -3.537  -10	  -10	  -10 	  -10
                   -10     -10    -10     -10     -2.829  -1.142  -0.292  -0.518  -0.998  -1.936  -3.365  -4.739  -10	  -10
                   -10     -10    -10     -10     -3.813  -1.825  -0.673  -0.502  -0.519  -0.890  -1.657  -2.341  -10     -10
                   -10     -10    -10     -10     -10     -2.928  -1.512  -1.003  -0.624  -0.513  -0.737  -0.851  -10     -10
                   -10     -10    -10     -10     -10     -4.536  -2.884  -2.080  -1.356  -0.836  -0.614  -0.254  -10     -10
                   -10     -10    -10     -10     -10     -10	  -4.396  -3.330  -2.302  -1.429  -0.843  -0.090  -4.710  -10
                   -10     -10    -10     -10     -10     -10	  -10	  -4.480  -3.181  -2.001  -1.113  -0.040  -3.461  -10
                   -10     -10    -10     -10     -10     -10	  -10	  -10	  -3.931  -2.482  -1.338  -0.023  -2.476  -10
                   -10     -10    -10     -10     -10     -10	  -10	  -10	  -4.560  -2.872  -1.504  -0.024  -1.699  -4.151
                   -10     -10    -10     -10     -10     -10	  -10	  -10	  -10	  -3.206  -1.636  -0.048  -1.104  -2.840
                   -10     -10    -10     -10     -10     -10	  -10	  -10	  -10	  -3.544  -1.787  -0.119  -0.680  -1.823
                   -10     -10    -10     -10     -10     -10	  -10	  -10	  -10	  -3.944  -2.014  -0.272  -0.429  -1.078
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -4.451  -2.354  -0.535  -0.353  -0.586
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -2.805  -0.898  -0.425  -0.305
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -3.319  -1.314  -0.584  -0.161
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -3.847  -1.738  -0.776  -0.089
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -4.368  -2.143  -0.970  -0.053
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -4.866  -2.524  -1.156  -0.033
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -10	  -2.879  -1.329  -0.021
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -10	  -3.211  -1.490  -0.015
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -10	  -3.522  -1.640  -0.010
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -10	  -3.815  -1.780  -0.007
                   -10     -10    -10     -10     -10 	  -10	  -10	  -10	  -10	  -10     -10	  -4.093  -1.910  -0.005];
	te_K = te(:)*1.2e4;	% Vecteur des températures en kelvins
	log10te_K = log10(te_K);
	frac_ab = tsample(frac_ab_ref,log10te_ref,log10te_K);
	AAleg = exp(frac_ab);
elseif iz == 14
elseif iz == 15
elseif iz == 16
end



if display_fig


figure
if iz==13
	loglog(te,AAleg)
	axis([10 2000 1e-2 1.05])
	grid
	title(' Aluminium equilibre coronal ')
elseif iz==10
	loglog(te,AAlegold)
	axis([10 2000 1e-2 1.05])
	grid
	title(' NEON  equilibre coronal ')
elseif iz==8
	loglog(te,AAleg)
	axis([10 1000 1e-2 1.05])
	grid
	title(' OXYGENE  equilibre coronal ')
elseif iz==7
	loglog(te,AAlegold)
	axis([10 700 1e-2 1.05])
	grid
	title('AZOTE equilibre coronal ')
elseif iz==6
	loglog(te,AAleg)
	%axis([10 400 1e-2 1.05])
	grid
	title(' CARBONE  equilibre coronal ')
elseif iz==5
	loglog(te,AAlegold)
	axis([10 300 1e-2 1.05])
	grid
	title('BORE equilibre coronal ')
elseif iz==4
	loglog(te,AAleg)
	axis([10 100 1e-2 1.05])
	grid
	title(' BERYLLIUM equilibre coronal ')
end
xlabel('Température (eV)')

end
