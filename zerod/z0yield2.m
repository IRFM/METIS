function [ya0,Limp,Lt,Eth,Q,Mimp,Mt]  = z0yield2(Te,Ti,Zimp,Zt,mode)

% list of elements
persistent name_liste
persistent num_liste

% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)




ya0   = NaN .* ones(size(Te));
Limp = '';
Lt  = '';
Eth = NaN;
Q   = NaN;

if nargin < 5
	mode = 1;
end

% reference pour les donnees :
%  N. Matsunami et al
% Energy Dependence of the Yields of Ion-Induced Sputtering of Monatomic Solids
% IPPJ-AM-32 (Institute of Plasma Physics, Nagoya University, Japan, 1983).
% repris de  http://www.iap.tuwien.ac.at/www/surface/sputteryield
% Michael Schmid, IAP/TU Wien Surface Physics Group 2006-2009.
% les donn?es proviennent de cette page web  
% name,Z,Mass,Ex,Q
liste = {
{'Ag',47,107.8682, 2.95, 1.21},
{'Al',13,26.98154, 3.39, 1.09},
{'Au',79,196.96654,3.81, 1.04},
{'B',  5,10.81,    5.77, 4.6},
{'Be', 4,9.01,     3.32, 2.17},
{'C', 6,12.011,    7.37, 3.1},
{'Cd',48,112.41,   1.16, 1},
{'Co',27,58.93,    4.39, 1.00},
{'Cr',24,51.996,   4.10, 1.23},
{'Cu',29,63.546,   3.49, 1.30},
{'Fe',26,55.845,   4.28, 1.06},
{'Ge',32,72.61,    3.85, 0.83},
{'Hf',72,178.49,   6.44, 0.75},
{'Ir',77,192.217,  6.94, 1.37},
{'Mn',25,54.94,    2.92, 1.13},
{'Mo',42,95.94,    6.82, 0.84},
{'Nb',41,92.90638, 7.57, 1.02},
{'Ni',28,58.6934,  4.44, 1.06},
{'Os',76,190.2,    8.17, 1.47},
{'Pb',82,207.2,    2.03, 1},
{'Pd',46,106.42,   3.89, 1.10},
{'Pt',78,195.08,   5.84, 1.13},
{'Re',75,186.21,   8.03, 1.27},
{'Rh',45,102.9055, 5.75, 1.26},
{'Ru',44,101.07,   6.74, 1.52},
{'Si',14,28.0855,  4.63, 0.78},
{'Sn',50,118.69,   3.14, 0.47},
{'Th',90,232.04,   6.20, 0.9},
{'Ti',22,47.9,     4.85, 0.58},
{'Ta',73,180.9479, 8.10, 0.78},
{'U', 92,238.03,   5.55, 0.81},
{'V', 23,50.9415,  5.31, 0.9},
{'W', 74,183.84,   8.90, 1.10},
{'Zn',30,65.38,    1.35, 1},
{'Zr',40,91.22,    6.25, 0.70},
{'Ar', 18,39.95},
{'D',   1,2},
{'T',   1,3},
{'Ga', 31,69.75},
{'H',   1,1},
{'He',  2,4},
{'He3', 2,3},
{'Kr', 36,83.8},
{'Li',  3,6.94},
{'Na', 11,22.99},
{'Ne', 10,20.18},
{'O',   8,16.00},
{'Xe', 54,131.3},
{'N',7,14},
{'Cl',17,34.453},
};

if isempty(name_liste)
    for k = 1:length(liste)
        name_liste{k} = liste{k}{1};
        num_liste(k)  = liste{k}{2};
    end
end
% recherche avec le nom ou la charge
if ischar(Zimp)
	 ind_impur = strmatch(Zimp,name_liste,'exact');
else
	 ind_impur = find(Zimp == num_liste,1);
end
if ischar(Zt)
	 ind_target = strmatch(Zt,name_liste,'exact');
else
	 ind_target = find(Zt == num_liste,1);
end

% % recherche avec le nom ou la charge
% ind_impur =  NaN;
% for k = 1:length(liste)
%     if ischar(Zimp)
% 	if strmatch(Zimp,liste{k}{1},'exact')
% 	      ind_impur = k;
% 	      break
% 	end
%     else
%        if liste{k}{2} == Zimp
% 	    ind_impur = k;
% 	    break
%        end
%     end
% end
% if ~isfinite(ind_impur)
% 	return
% end
% ind_target =  NaN;
% for k = 1:length(liste)
%     if ischar(Zt)
% 	if strmatch(Zt,liste{k}{1},'exact')
% 	      ind_target = k;
% 	      break
% 	end
%     else
%        if liste{k}{2} == Zt
% 	    ind_target = k;
% 	    break
%        end
%     end
% end
if ~isfinite(ind_target)
	return
elseif length(liste{ind_target}) < 5
	return
end

Limp = liste{ind_impur}{1};
Zimp = liste{ind_impur}{2};
Mimp = liste{ind_impur}{3};
Lt   = liste{ind_target}{1};
Zt   = liste{ind_target}{2};
Mt   = liste{ind_target}{3};
Ex   = liste{ind_target}{4};
Q    = liste{ind_target}{5};


% correction de Zimp et Zt pour W
if Zimp == 74
    Zimp = z0wavez(Te);
end
if Zt == 74
    Zt = z0wavez(Te);
end

% correction de Zimp et Zt pour Sn
if Zimp == 50
    Zimp = z0snavez(Te);
end
if Zt == 50
    Zt = z0snavez(Te);
end

% lien entre la temperature et l'energie incidente 
% Material migration in divertor tokamaks
% G.F. Matthews / Journal of Nuclear Materials 337-339 (2005)
%Ei  =  3 .* Zimp .* Te + 2 .* Ti;
% new more complete expression (thank to Rajiv Goswami).
% M. Warrier et al. / Computer Physics Communications 160 (2004) 46?68
phis = 0.5 .* Te .* log(2.* pi .* phys.me ./ phys.ua ./ Mimp .* (1 + Ti ./ Te));
gamma_heat = 1 ; %  1 = isothermal ; 5/3 = adiabatic.
Ei = Te + (2 + gamma_heat) .* Ti - Zimp .* phis;

%  figure(21);
%  subplot(3,1,1)
%  plot(Te,Ei,'.',Te, 3 .* Zimp .* Te + 2 .* Ti);
%  subplot(3,1,2);
%  plot(Te,Zimp,'.');
%  subplot(3,1,3);
%  plot(Te,Zt,'.');
%  drawnow


if mode == 0
    % formule :N. Matsunami et al
    %
    alpha = 2.8;
    epsi  = 0.03255 .* Ei .* Mt ./ (Mt + Mimp) ./ Zimp ./ Zt ./ sqrt(Zimp .^ (2/3) + Zt .^ (2/3));   
    K     = 8.478 .* Zimp .* Zt ./ sqrt(Zimp .^ (2/3) + Zt .^ (2/3)) .* Mimp ./ (Mimp + Mt);
    Eth   = Ex .* (1.9 + 3.8 .* Mimp ./ Mt + 0.134 .* (Mt ./ Mimp) .^ 1.24);
    a     = 0.08 + 0.164 .*    (Mt ./ Mimp) .^ 0.4 + 0.0145 .* (Mt ./ Mimp) .^ 1.29;
    b     = 3.441 .* sqrt(epsi ) .* log(epsi +exp(1)) ./ (1 + 6.355 .* sqrt(epsi) + epsi .* (6.882 .* sqrt(epsi) -1.708));
    c     = sqrt(epsi) .* 0.079 .* (Mimp + Mt) .^ (3/2) ./ Mimp ./ sqrt(Mimp) ./ sqrt(Mt) .* Zimp .^ (2/3) .* sqrt(Zt) ./ (Zimp .^(2/3) + Zt .^(2/3)) .^ (3/4);
    ya0   = 0.42 .* a .* Q .* K .* b ./ Ex ./ (1 + 0.35 .* Ex .* c) .*  (1 - sqrt(Eth ./ Ei)) .^ alpha;
    ya0(Eth > Ei) = 0;
    
else
    % formule :
    % Unified analytic representation of physical sputtering yield
    % Journal of Nuclear Materials 290??293 (2001) 104-106
    % R.K. Janev a,b, Yu.V. Ralchenko a,c, T. Kenmotsu a,*, K. Hosaka a
    alpha = 3;
    Apar  = 0.436;
    Bpar  = 0.212;
    
    Eth   = Ex .* (1.9 + 3.8 .* Mimp ./ Mt + 0.134 .* (Mt ./ Mimp) .^ 1.24);
    Etf   = 30.74 .* (Mimp + Mt) ./ Mt .* Zimp .* Zt .* sqrt(Zimp .^(2/3) + Zt .^ (2/3));
    
    epsi  = Ei ./ Etf;
    delta = Eth ./ Etf;
    a     = 1.265 .* delta ./ (0.18 + delta .^ (2/3));
    b     = 20.5 .* delta .^ (2/5) ./ (1 + 112 .* delta);
    gamma = 0.81 .* (0.0051 + delta .^ (4/5)) ./ (0.013 + delta .^ (3/5));
    g     = 0.85 + 4.0 .* exp(-2.94 .* delta .^(3/5));
    
    eta   = a .* (epsi ./ delta - 1) + b .* ((epsi ./ delta) .^ gamma -1) + 1;
    yt    = (1 - 1 ./ eta) .^ alpha .* (Apar .* log(eta) ./ eta  + Bpar ./ eta .^ 2);
    ya0   = yt .* Q .* g;
     
end


ya0(imag(ya0) ~= 0) = 0;
ya0(ya0 < 0) = 0;
