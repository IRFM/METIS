% Calcul de l'intervalle de Te qui contient 90% du Fe XV à l'équilibre coronal
%
%	R.? Guirlet, 12/10/05
%
% Vecteur de temérature électronique
Te = 1:3000;
% Abondance fractionnelle du Fe XV d'après equicoronal_lourd.m
load Equi_coronal_Fe
Abondance_FeXV = Abondance_Fe(:,15);
%Abondance_FeXV_tot = sum(Abondance_FeXV);
%
%[max_Abondance_FeXV, i_max_Abondance_FeXV] = max(Abondance_FeXV);
%
% PEC du Fe XV 284.1 A (1 colonne = 1 température, toutes les densités)
ne_PEC = [1.00E+08 2.34E+08 5.46E+08 1.27E+09 2.98E+09 6.95E+09 1.62E+10 3.79E+10 8.86E+10 2.07E+11 4.83E+11 1.13E+12 2.64E+12 6.16E+12 1.44E+13 3.36E+13 7.85E+13 1.83E+14 4.28E+14 1.00E+15]';
Te_PEC = [1 9.69E+00 1.36E+01 1.94E+01 2.91E+01 3.88E+01 5.82E+01 9.69E+01 1.36E+02 1.94E+02 2.91E+02 3.88E+02 5.82E+02 9.69E+02 1.36E+03 1.94E+03 2.91E+03 3.88E+03 5.82E+03 9.69E+03 1.36E+04 1.94E+04 2.91E+04 3.88E+04 5.82E+04];
PEC_FeXV284 = [7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09 ...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09 ;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09 ...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09 ;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09 ...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09 ;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.39E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.30E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.38E-08 1.39E-08...
 			1.36E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.29E-09 4.97E-09 8.46E-09 1.06E-08 1.27E-08 1.38E-08 1.39E-08...
 			1.35E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.99E-09 6.20E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.29E-09 4.97E-09 8.45E-09 1.06E-08 1.27E-08 1.38E-08 1.39E-08...
 			1.35E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.98E-09 6.19E-09 5.31E-09 4.79E-09 4.29E-09 3.79E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.29E-09 4.96E-09 8.44E-09 1.06E-08 1.27E-08 1.38E-08 1.39E-08...
 			1.35E-08 1.29E-08 1.23E-08 1.14E-08 1.02E-08 9.35E-09 8.50E-09 7.59E-09...
 			6.98E-09 6.19E-09 5.31E-09 4.79E-09 4.29E-09 3.78E-09 3.46E-09 3.05E-09;
 			7.54E-10 7.54E-10 2.29E-09 4.95E-09 8.41E-09 1.05E-08 1.26E-08 1.38E-08 1.38E-08...
 			1.35E-08 1.28E-08 1.23E-08 1.14E-08 1.01E-08 9.34E-09 8.49E-09 7.58E-09...
 			6.98E-09 6.19E-09 5.31E-09 4.79E-09 4.29E-09 3.78E-09 3.46E-09 3.05E-09;
 			7.55E-10 7.55E-10 2.29E-09 4.92E-09 8.35E-09 1.05E-08 1.25E-08 1.37E-08 1.37E-08...
 			1.34E-08 1.28E-08 1.22E-08 1.13E-08 1.01E-08 9.32E-09 8.48E-09 7.57E-09...
 			6.97E-09 6.18E-09 5.30E-09 4.78E-09 4.29E-09 3.78E-09 3.46E-09 3.05E-09;
 			7.57E-10 7.57E-10 2.27E-09 4.87E-09 8.21E-09 1.03E-08 1.23E-08 1.34E-08 1.35E-08...
 			1.32E-08 1.26E-08 1.21E-08 1.12E-08 1.01E-08 9.27E-09 8.44E-09 7.55E-09...
 			6.95E-09 6.17E-09 5.29E-09 4.78E-09 4.28E-09 3.78E-09 3.46E-09 3.05E-09;
 			7.60E-10 7.60E-10 2.25E-09 4.75E-09 7.93E-09 9.88E-09 1.18E-08 1.29E-08 1.31E-08...
 			1.28E-08 1.23E-08 1.18E-08 1.10E-08 9.93E-09 9.17E-09 8.36E-09 7.49E-09...
 			6.90E-09 6.14E-09 5.27E-09 4.76E-09 4.27E-09 3.77E-09 3.45E-09 3.04E-09];
% La dépendance en densité est très faible --> on moyenne sur ne
PEC_FeXV284 = mean(PEC_FeXV284);
% Rééchantillonage sur la température du calcul de l'équilibre coronal
PEC_FeXV284 = tsample(PEC_FeXV284',Te_PEC',te);
%
Emiss_FeXV_284 = Abondance_FeXV.*PEC_FeXV284;
Emiss_FeXV_284_tot = sum(Emiss_FeXV_284);
[max_Emiss_FeXV_284, i_max_Emiss_FeXV_284] = max(Emiss_FeXV_284);

%
% 1e méthode: on intègre l'abondance fractionnelle autour du maximum et on
% en déduit la fraction de la densité totale dans la bande de Te correspondante
% Te_min et Te_max sont les limites du domaine d'intégration. 
% On fait varier la largeur de ce domaine en intégrant symétriquement autour du max:
% à + ou - un point de Te, 2 points,... 100 points.
% Inconvénient: la fraction non prise en compte à plus haute température est
% beaucoup plus grande que celle à plus basse température
i_max = 100;
for i=1:i_max
	fraction_Ab1(i) = sum(Emiss_FeXV_284((i_max_Emiss_FeXV_284-i):(i_max_Emiss_FeXV_284+i))) / Emiss_FeXV_284_tot;
	fraction_Ab1_dessous(i) = sum(Emiss_FeXV_284(1:(i_max_Emiss_FeXV_284-i-1))) / Emiss_FeXV_284_tot;
	fraction_Ab1_dessus(i)  = sum(Emiss_FeXV_284((i_max_Emiss_FeXV_284+i+1):length(Emiss_FeXV_284))) / Emiss_FeXV_284_tot;
	Te_min1(i) = te(i_max_Emiss_FeXV_284-i);
	Te_max1(i) = te(i_max_Emiss_FeXV_284+i);
end
figure
subplot(311)
	plot((1:i_max)*mean(diff(te)),fraction_Ab1,'r',...
		 (1:i_max)*mean(diff(te)),fraction_Ab1_dessous,'r--',...
		 (1:i_max)*mean(diff(te)),fraction_Ab1_dessus,'r-.',...
		 (1:i_max)*mean(diff(te)),fraction_Ab1+fraction_Ab1_dessous+fraction_Ab1_dessus,'r:')
	set(gca,'ylim',[0 1.05])
	grid
	title('Fraction de l''émissivité totale de FeXV 28.41 nm - méthode 1')
subplot(312)
	plot((1:i_max)*mean(diff(te)),Te_min1,'r')
	ylabel('eV')
	grid
	title('T_e min')
subplot(313)
	plot((1:i_max)*mean(diff(te)),Te_max1,'r')
	xlabel('1/2 largeur de l''intervalle d''intégration (eV)')
	ylabel('eV')
	grid
	title('T_e max')
%
% 2e méthode: on cherche la position des points de part et d'autre du maximum
% où l'abondance est une fraction alpha du maximum. On repère les températures 
% correspondantes et on calcule la fraction du total comprise dans cette bande.
% La dissymétrie est plus faible qu'avec la méthode précédente mais encore pas
% nulle.
i_Ab0 = min(find(Emiss_FeXV_284>=1e-4*max_Emiss_FeXV_284));
for i=1:19
	alpha(i) = 1 - i*.05;
	i_min = iround(Emiss_FeXV_284((i_Ab0+1):i_max_Emiss_FeXV_284)/max_Emiss_FeXV_284,alpha(i));
	i_max = iround(Emiss_FeXV_284(i_max_Emiss_FeXV_284:length(Emiss_FeXV_284))/max_Emiss_FeXV_284,alpha(i));
	fraction_Ab2(i) = sum(Emiss_FeXV_284((i_Ab0+i_min):(i_max_Emiss_FeXV_284+i_max))) / Emiss_FeXV_284_tot;
	fraction_Ab2_dessous(i) = sum(Emiss_FeXV_284(1:(i_Ab0+i_min-1))) / Emiss_FeXV_284_tot;
	fraction_Ab2_dessus(i)  = sum(Emiss_FeXV_284((i_max_Emiss_FeXV_284+i_max+1):length(Emiss_FeXV_284))) / Emiss_FeXV_284_tot;
	Te_min2(i) = te(i_Ab0+i_min);
	Te_max2(i) = te(i_max_Emiss_FeXV_284+i_max-1);
end
figure
subplot(311)
	plot(alpha,fraction_Ab2,'b',...
		 alpha,fraction_Ab2_dessous,'b--',...
		 alpha,fraction_Ab2_dessus,'b-.',...
		 alpha,fraction_Ab2+fraction_Ab2_dessous+fraction_Ab2_dessus,'b:')
	set(gca,'ylim',[0 1.05])
	grid
	title('Fraction de l''émissivité totale de FeXV 28.41 nm - méthode 2')
subplot(312)
	plot(alpha,Te_min2,'b')
	ylabel('eV')
	grid
	title('T_e min')
subplot(313)
	plot(alpha,Te_max2,'b')
	xlabel('Fraction du maximum d''émissivité de Fe XV 28.41 nm')
	ylabel('eV')
	grid
	title('T_e max')
