% for testing:
% eta_t = z0eta_turb(post.profil0d.nep,post.profil0d.tep,post.profil0d.epsi,post.profil0d.zeff)

% turbulent resistivity 
function eta_t = z0eta_turb(ne,te,epsi,zeff,Diff)

% phyical constants 
phys = cphys;

% model sourced from :
% L. Colas and G. Giruzzi 1993 Nucl. Fus. 33 156
% I. Chadavdarovski and R. Gatto PoP 24 072512 2017
%vth     = sqrt(phys.e .* te ./ phys.me);
%nuevth3 = 4 .* pi .* phys.e .^ 4 .* ne .* lambda_ln_e(ne,te) ./ phys.me .^ 2;
mu_t    = sqrt(2 .* epsi ./ (1 + epsi));
Z       = max(1,min(5,zeff));
delta   = (0.15 + 0.04 .* Z) .* mu_t + (0.9 - 0.15 .* Z) .* mu_t .^ 2;
I_mu    = (7 + 2 .* zeff) ./ (3 + zeff) - delta;
%sigma_t = 2 ^ 6 .* sqrt(2 .* pi) .* (phys.e ./ nuevth3) .^ 2 .*  ...
%          I_mu ./ (5 + Z) ./ phys.me .* vth .^ 7  .* ne .* D_tild;
% for simplification, used Spitzer * factor to get right parameters dependances
sigma_t = 1.9012e4 .* te .^ (3/2) .* (I_mu ./ (5 + Z))./ lambda_ln_e(ne,te)  .* ...
          (te ./ (max(te,[],2) * ones(1,size(te,2)))) .^ 2  ./  ...
          (ne ./ (max(ne,[],2) * ones(1,size(ne,2)))) .* Diff; 
eta_t   = 1 ./ max(eps,sigma_t); % normally decrease of plasma resistivity with turbulence !

function out = lambda_ln_e(ne,te)

out = 31.3 - log(sqrt(ne)./te); % OK


function phys = cphys

% constante physique (phys)
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.s3igma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)