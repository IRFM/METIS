% Calcul de l'intervalle de Te qui contient 90% du Fe XV à l'équilibre coronal
%
%	R.? Guirlet, 12/10/05
%
% Vecteur de temérature électronique
Te = 1:3000;
% Abondance fractionnelle du Fe XV d'après equicornal_lourd.m
load Equi_coronal_Fe
Abondance_FeXV = Abondance_Fe(:,15);
Abondance_FeXV_tot = sum(Abondance_FeXV);
%
[max_Abondance_FeXV, i_max_Abondance_FeXV] = max(Abondance_FeXV);
%
%
% PEC de la raie Fe XV 28.41 nm
Te_PEC = [		50       100       150       200       250       350		450    550];	% en eV
ne_PEC = [1e12	3e12	6e12	1e13	3e13	6e13	1e14	3e14];						% en cm-3
ne_PEC = ne_PEC * 1e6;																		% en m-3
% PEC: une colonne = une température, une ligne = une densité;	cm3.s-1
PEC = [		1.62e+04  1.92e+04  1.91e+04  1.86e+04  1.81e+04  1.71e+04	1.63e+04  1.56e+04;
		    4.86e+04  5.75e+04  5.73e+04  5.59e+04  5.43e+04  5.13e+04	4.90e+04  4.69e+04;
			9.71e+04  1.15e+05  1.15e+05  1.12e+05  1.08e+05  1.03e+05	9.79e+04  9.37e+04;
			1.62e+05  1.92e+05  1.91e+05  1.86e+05  1.81e+05  1.71e+05	1.63e+05  1.56e+05;
			4.85e+05  5.74e+05  5.73e+05  5.59e+05  5.42e+05  5.12e+05	4.89e+05  4.68e+05;
			9.70e+05  1.15e+06  1.14e+06  1.12e+06  1.08e+06  1.02e+06	9.78e+05  9.36e+05;
			1.61e+06  1.91e+06  1.90e+06  1.86e+06  1.80e+06  1.70e+06	1.63e+06  1.56e+06;
			4.82e+06  5.70e+06  5.68e+06  5.55e+06  5.39e+06  5.09e+06	4.86e+06  4.66e+06];
PEC = PEC * 1e-6;			%	m3.s-1

%
% 1e méthode: on intègre l'abondance fractionnelle autour du maximum et on
% en déduit la fraction de la densité totale dans la bande de Te correspondante
% Te_min et Te_max sont les limites du domaine d'intégration. 
% On fait varier la largeur de ce domaine en intégrant à + ou - un point de Te,
% 2 points,... 100 points.
% Inconvénient: la fraction non prise en compte à plus haute température est
% beaucoup plus grande que celle à plus basse température
for i=1:100
	fraction_Ab1(i) = sum(Abondance_FeXV((i_max_Abondance_FeXV-i):(i_max_Abondance_FeXV+i))) / Abondance_FeXV_tot;
	Te_min1(i) = te(i_max_Abondance_FeXV-i);
	Te_max1(i) = te(i_max_Abondance_FeXV+i);
end
figure
subplot(311)
	plot(fraction_Ab1)
	grid
	title('Fraction de la densité totale de FeXV - méthode 1')
subplot(312)
	plot(Te_min1)
	ylabel('eV')
	grid
	title('T_e min')
subplot(313)
	plot(Te_max1)
	xlabel('Fraction du maximum d''abondance')
	ylabel('eV')
	grid
	title('T_e max')
%
% 2e méthode: on cherche la position des points de part et d'autre du maximum
% où l'abondance est une fraction alpha du maximum. On repère les températures 
% correspondantes et on calcule la fraction du total comprise dans cette bande.
% La dissymétrie est plus faible qu'avec la méthode précédente mais encore pas
% nulle.
i_Ab0 = min(find(Abondance_FeXV>=1e-4));
for i=1:19
	alpha(i) = 1 - i*.05;
	i_min = iround(Abondance_FeXV((i_Ab0+1):i_max_Abondance_FeXV)/max_Abondance_FeXV,alpha(i));
	i_max = iround(Abondance_FeXV(i_max_Abondance_FeXV:length(Abondance_FeXV))/max_Abondance_FeXV,alpha(i));
	fraction_Ab2(i) = sum(Abondance_FeXV(i_min:(i_max_Abondance_FeXV+i_max))) / Abondance_FeXV_tot;
	Te_min2(i) = te(i_Ab0+i_min);
	Te_max2(i) = te(i_max_Abondance_FeXV+i_max-1);
end
figure
subplot(311)
	plot(alpha,fraction_Ab2)
	grid
	title('Fraction de la densité totale de FeXV - méthode 2')
subplot(312)
	plot(alpha,Te_min2)
	ylabel('eV')
	grid
	title('T_e min')
subplot(313)
	plot(alpha,Te_max2)
	xlabel('Fraction du maximum d''abondance')
	ylabel('eV')
	grid
	title('T_e max')
