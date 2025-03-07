% averaged energy of He4 in TT reaction
% data 

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

dE = 11.3e6; %eV
c2 = 20;
en = 2 .* dE .* phys.e ./ phys.mp;

% Monte Carlo loop
res =[];
for k=1:1e5
  %v1 = rand(1) .* sqrt(en);
  v1 = randn(1) .* sqrt(en);
  cu = 1 - 2 .* rand(1);
  c1 = 4 .* v1 .* cu;
  c0 = 2 * v1 .^ 2 - en;
  [x1,x2]= zpolyroot(c0,c1,c2);
  if ~iscomplex(x1)
      res(end + 1) = abs(x1);
  end
  if ~iscomplex(x2)
      res(end + 1) = abs(x2);
  end
end
figure
hist(res);
ehe4_mean = mean(4 .* phys.mp .* res .^ 2 ./ 2) ./ phys.e
ehe4_std = std(4 .* phys.mp .* res .^ 2 ./ 2) ./ phys.e
