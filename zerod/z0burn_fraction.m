function [burn_fraction,averaged_burn_fraction,source_T_plasma,source_T_pellet,source_T_gaz,source_T_nbi,source_leakage,source_recycling]=z0burn_fraction(post,eta_gaz,eta_pellet,pellet_fraction)

% physical constant
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
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3

%all in Pa * m^3 /s
% usefull vector 
vt = ones(size(post.z0dinput.cons.temps));

% pellet fraction
if (nargin < 4) || all(pellet_fraction == 0);
    frac_pellet = post.zerod.frac_pellet .* (post.zerod.frac_pellet > 1e-2);
elseif length(pellet_fraction) == 1
    frac_pellet = vt * pellet_fraction;
else
    frac_pellet = pellet_fraction;
end 
% NBI contribution
% case of Hydrogen NBI in DT plasma
if isfield(post.z0dinput.option,'forced_H_NBI') && (post.z0dinput.option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = post.z0dinput.option.gaz; 
end

% totale source ?
switch gas_nbi
    case -1
        source_T_nbi  = 0 .* vt;
        perte_T_nbi   = 0 .* vt;
        
    case 3
        source_T_nbi  = real(post.z0dinput.cons.ftnbi) .* real(post.zerod.pnbi)  ./ post.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo) + ...
            imag(post.z0dinput.cons.ftnbi) .* imag(post.zerod.pnbi)  ./ post.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo);
        
        perte_T_nbi   = real(post.z0dinput.cons.ftnbi) .* real(post.zerod.pnbi)  ./ post.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(post.zerod.frnbi)) .* (1 - real(post.zerod.frnbi)) + ...
            imag(post.z0dinput.cons.ftnbi) .* imag(post.zerod.pnbi)  ./ post.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(post.zerod.frnbi)) .* (1 - imag(post.zerod.frnbi));
        
        
    otherwise
        source_T_nbi  = 0 .* vt;
        perte_T_nbi   = 0 .* vt;
end

% time derivative plasma content:
N_T = post.zerod.vp .* post.zerod.nTm ./ (4.41e-4 .* phys.avo);
dN_Tdt = z0dxdt(N_T,post.z0dinput.cons.temps);

% internal term 
internal = post.zerod.salpha ./ (4.41e-4 .* phys.avo) + dN_Tdt;

% Losses of Tritium from central plasma to the egde (not edge fast losses linked to tau_p) : tau_He without rdivertor recycling
if (post.z0dinput.option.tauhemul < 0)  && (post.z0dinput.option.Recycling < 1)
        % modele qui prend en compte le confinement reel et le recyclage dans le divertor
	% ref Stangeby section 6.7 
        % ref originale : D. Reiter et al, PPCF vol 33 (1991) p 1579-1600
	tauhe    = post.zerod.tauhe - post.z0dinput.option.Recycling ./ (1 - post.z0dinput.option.Recycling) .* post.zerod.taup;

else
	tauhe    = post.zerod.tauhe;
end
source_leakage   = N_T ./ tauhe;
% fraction qui retourne dans le plasma sous forme de neutres
switch post.z0dinput.option.configuration
case {0,1}
    source_recycling = post.z0dinput.option.fn0a .* source_leakage;
case {2,3}
    source_recycling    = (post.zerod.xpoint .* post.z0dinput.option.fn0a_div + (~post.zerod.xpoint) .* post.z0dinput.option.fn0a) .* source_leakage;
otherwise
    source_recycling = post.z0dinput.option.fn0a_div .* source_leakage;
end
if (nargin >= 4) && any(pellet_fraction ~= 0)
  fpe_old          = min(1,post.zerod.taup ./ max(post.zerod.taup,post.zerod.taue)) .* post.zerod.frac_pellet;
  fpe              = min(1,post.zerod.taup ./ max(post.zerod.taup,post.zerod.taue)) .* frac_pellet;
  source_recycling = source_recycling.* max(0.1,1 - fpe) ./ max(0.1,1 - fpe_old);
end
% input effective source in plasma
source_T_plasma = internal + source_leakage - source_recycling;

% source_T_pellet
source_T_pellet = frac_pellet .* (source_T_plasma - source_T_nbi) ./ max(eps,eta_pellet);

% source_T_gaz
source_T_gaz    = (1 - frac_pellet) .* (source_T_plasma - source_T_nbi) ./ eta_gaz;
source_T_gaz(~isfinite(source_T_gaz)) = 0;

% burn fraction of tritium
burn_fraction   = internal ./ (source_leakage - source_recycling + (1 - eta_gaz) .* source_T_gaz + (1 - eta_pellet) .* source_T_pellet + perte_T_nbi);

% global shot performance:
averaged_burn_fraction = trapz(post.z0dinput.cons.temps, burn_fraction .* N_T,1) ./ max(eps,trapz(post.z0dinput.cons.temps, N_T,1));



