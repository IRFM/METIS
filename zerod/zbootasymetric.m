% this function compute the asymetric counter part of bootstrap 
% ref : Yu. V. Gott and E.I. Yurchenko, PoP 16 (2009) p 112502
% x = flux surface label
% q(t,x) = safety factor
% epsi(t,x) = inverse aspect ratio
% Raxe(t,x) = flux surface radius (m)
% btot(t,x) = total magnetic field (m)
% dens(t,x) = density (m^-3)
% temp(t,x) = temperature (eV)
% mass = number of mass (if charge >0)
% charge = number of charge 
% b0     = normalisation magnetic field (T)
%
% ja = asymetric current
% jb = bootstrap evaluation.
function [ja,jb] = zbootasymetric(x,q,epsi,Raxe,btot,dens,temp,mass,charge,b0)


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
%
if charge < 0
  z = phys.e;
  m = phys.me;
else
  z = charge .* phys.e;
  m = mass .* phys.mp;
end
%
if all(size(b0) ==1)
  b0 = b0 .* ones(size(q));
elseif size(b0,2) == 1
  b0 = b0 * ones(size(x));
end
%
if all(size(z) ==1)
  z = z .* ones(size(q));
elseif size(z,2) == 1
  z = z * ones(size(x));
end
%
if all(size(m) ==1)
  m = m .* ones(size(q));
elseif size(m,2) == 1
  m = m * ones(size(x));
end

% larmor radius
vth  = sqrt(2 .* phys.e .* temp ./ m);
wc   = z .* btot ./ m;
rho  = vth ./wc;
% drift parameter
seta = 2 .* q .* rho ./ Raxe;

% current (// to B)
ct = 1;
ftr = sqrt(epsi + seta .^ (2/3));
p   = phys.e .* temp .* dens;
warning off
dpdepsi = pdederive(x,p,0,2,2,1) ./  pdederive(x,epsi,0,2,2,1);
dpdepsi(:,1) = 0;
%
ja = 2 .* ct .* p .* q ./ ftr ./ Raxe  ./ b0;
jb = -2.44 .* ftr .* ct .* dpdepsi .* q ./ epsi ./ Raxe ./ b0;
jb(:,1) = 0;
warning on
