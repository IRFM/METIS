% fit de Peres de la section efficace DT
% utilisee pour les calculs ITER 
% unite : systeme SI sauf T en eV
% section efficace integree de la reaction D-T 
function sigmavdt = zperes(T,phys)

% en kev 
T = 1e-3 .* T;

% la fontion est autonome pour usage externe a cronos
if nargin < 2
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
end

% contsantes du fit
p1 = 9.494748e-16;
p2 = 2.818421e-2;
p3 = 6.116184e-2;
p4 = 2.834474e-3;
p5 = 8.955113e-3;
p6 = -5.734052e-5;

% rapport de masse 
muc2 = 1.124656e6;

% B
B = sqrt(2).* phys.e .^ 2 ./ 4 ./ phys.epsi0 ./ ...
    (phys.h /2 / pi) ./ phys.c .* sqrt(muc2);

% fit
u = 1 - T .* (p2 + (p4 + p6 .* T) .* T) ./ (1 + (p3 + p5 .* T) .* T);
b = 3 .* (B ./ 2) .^ (2/3);
a = p1 .* (B ./2) .^ (1/3) ./ sqrt(muc2);

%la section efficace
sigmavdt = a ./ T .^ (2/3) ./ u .^ (5/6) .* exp (- b .* (u ./ T) .^ (1/3));

