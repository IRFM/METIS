% ref : R. S. Sutherland, Mon. Not. R. Astron. Soc. 300, 321-330 (1998).
function [gg,Cbrem,gamma2,tab_gaunt] = z0gaunt_brem(te,Z)

% te en keV
Ry     = 13.605698; % eV
gamma2 = Z .^ 2 .* Ry ./ (te .* 1000);


tab_gaunt =[
-20.00   2.*sqrt(3)/pi
-10.00   2.*sqrt(3)/pi
-4.00   1.11388
-3.80   1.11698
-3.60   1.12089
-3.40   1.12581
-3.20   1.13200
-3.00   1.13975
-2.80   1.14945
-2.60   1.16149
-2.40   1.17635
-2.20   1.19447
-2.00   1.21622
-1.80   1.24182
-1.60   1.27104
-1.40   1.30328
-1.20   1.33711
-1.00   1.37040
-0.80   1.40029
-0.60   1.42365
-0.40   1.43768
-0.20   1.44060
 0.00   1.43220 
 0.20   1.41391
 0.40   1.38830
 0.60   1.35832
 0.80   1.32658
 1.00   1.29496
 1.20   1.26462
 1.40   1.23618
 1.60   1.20993
 1.80   1.18594
 2.00   1.16421
 2.20   1.14464
 2.40   1.12711
 2.60   1.11147
 2.80   1.09757
 3.00   1.08526
 3.20   1.07438
 3.40   1.06481
 3.60   1.05640
 3.80   1.04904
 4.00   1.04264
 10.00   1.0
 20.00   1.0
 ];



gg = interp1(tab_gaunt(:,1),tab_gaunt(:,2),log10(gamma2),'pchip',1);


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

% Te ev keV 
Cbrem = phys.e .^ 6 ./ ( 6 .* sqrt(3/2) .* pi .^ (3/2) .* phys.epsi0 .^ 3 .* phys.c .^ 3 .* phys.h .* phys.me .^ (3/2)).* sqrt(phys.e.*1000);











