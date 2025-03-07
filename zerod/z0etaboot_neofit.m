% formulation from:
% "A new set of analytical formulae for the computation of the bootstrap current and the neoclassical conductivity in tokamaks"
% A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664
%
function [eta,jboot_b0,l31] = z0etaboot_neofit(x,te,ti,ne,ni,q,zeff,zion,Raxe,ft,epsi,dpsidx,F,modeh,fgg)

% physical constants 
persistent phys
phys = cphys;

% le facteur est un fit avec NClass (l'integrale n'est pas precise au bord)
if nargin < 15
    fgg = 1;
end



% common quantities
sp    = sigma_spitzer(ne,te,zeff); % remark in Sauter paper erratum even if formula cited in Redl paper point to original Sauter paper.
nue   = nue_star(q,Raxe,epsi,ne,te,zion);
nui   = nui_star(q,Raxe,epsi,ni,ti,zion);

% resistivity
x33   = ft_eff_33(zeff,ft,nue);
sr    = sigma_ratio(x33,zeff);
sigma = sp .*  sr;
eta   = 1 ./ sigma;
if nargout == 1
    return;
end

% variables for LXX function
x31    = ft_eff_31(zeff,ft,nue);
x32e   = ft_eff_32_e(zeff,ft,nue);
x32ei  = ft_eff_32_ei(zeff,ft,nue);
x34    = x31;
% alpha
alpha = falpha(zeff,ft,nui);
% LXX 
l31   =  L31(x31,zeff);
if nargout == 3
    jboot_b0   = [];
    return;
end
l32   =  L32(x32e,x32ei,zeff);
l34   =  l31; % new choice in this paper; different from previous publication

% pressure
pel  = phys.e .* ne .* te;
pion = phys.e .* ni .* ti;
ptot = pel + pion;

% computation off gradient
dnedx        = pdederive(x,ne,0,2,2,1);
dtedx        = pdederive(x,te,0,2,2,1);
dtidx        = pdederive(x,ti,0,2,2,1);
% correction piedestal (derivee 2 points + air du rectangle / aire triangle)
indh = find(modeh);
if ~isempty(indh) && (abs(fgg-1) > sqrt(eps)) && (fgg ~= 0)  % this formulation has the correct formulation for pedestal (so it is not useful to activate this mechanism 
	dx                = x(end) - x(end-1);
	dnedx(indh,end-1) = fgg .* (ne(indh,end) - ne(indh,end-1)) ./ dx;
	dtedx(indh,end-1) = fgg .* (te(indh,end) - te(indh,end-1)) ./ dx;
	dtidx(indh,end-1) = fgg .* (ti(indh,end) - ti(indh,end-1)) ./ dx;
end
dln_ne_dx = dnedx ./ ne;
dln_te_dx = dtedx ./ te;
dln_ti_dx = dtidx ./ ti;

% boostrap current <J_BS.B>
jboot_b0 = F ./ dpsidx .*  (ptot .* l31 .* dln_ne_dx +  ...
                             pel .* (l31 + l32) .* dln_te_dx  + ...
                             pion .* (l31 + alpha .* l34) .* dln_ti_dx);

% les gradients s'annulent au centre
jboot_b0(:,1)     = 0;                         
                        
% rep = whos;
% if any([rep.complex]) || any(~isfinite(jboot_b0(:)))
%      keyboard
% end

function out = L31(x31,zeff)

fzeff = zeff .^ 1.2 - 0.71; % OK
out = (1 + 0.15 ./ fzeff) .* x31 - 0.22 ./ fzeff .* x31 .^ 2 +  ...
       0.01 ./ fzeff .* x31 .^ 3 + 0.06 ./ fzeff .* x31 .^ 4; %  OK
      

function out = L32(x32e,x32ei,zeff)

out = F32_ei(x32ei,zeff) + F32_ee(x32e,zeff); % OK


function out = F32_ee(x32e,zeff)

out = (0.1 + 0.6 .* zeff) ./ (zeff .* (0.77 + 0.63 .* (1 + (zeff -1) .^ 1.1))) .* (x32e - x32e .^ 4) + ...
      0.7 ./ (1 + 0.2 .* zeff) .* (x32e .^ 2 - x32e .^ 4  - 1.2 .* (x32e .^ 3 - x32e .^ 4)) + ...
      1.3 ./ (1 + 0.5 .* zeff) .* x32e .^ 4; % OK


function out = F32_ei(x32ei,zeff)

out = - (0.4 + 1.93 .* zeff) ./ zeff ./ (0.8 + 0.6 .* zeff) .*  (x32ei - x32ei .^ 4) + ...
        5.5 ./ (1.5 + 2 .* zeff) .* (x32ei .^ 2 - x32ei .^ 4  - 0.8 .* (x32ei .^ 3 - x32ei .^ 4)) - ...
        1.3 ./ (1 + 0.5 .* zeff) .* x32ei .^ 4; % OK


function out = F34(x34,zeff)

out =  L31(x34,zeff); % OK


function out = ft_eff_31(zeff,ft,nue)

term1  = 0.67 .*(1 - 0.7 .* ft) .* sqrt(nue) ./ (0.56 + 0.44 .* zeff);  % OK
term2  = (0.52 + 0.086 .* sqrt(nue)) .* (1+ 0.87 .* ft) .* nue ./ ...
         (1+1.13 .* sqrt(zeff -1)); % OK
out    = ft ./ (1 + term1 + term2); % OK


function out = ft_eff_32_e(zeff,ft,nue)

term1  = 0.23 .* (1 - 0.96 .* ft ) .* sqrt(nue) ./ sqrt(zeff); % OK
term21 = sqrt(1 + 2 .* sqrt(zeff - 1));  % OK 
term22 = ft .^ 2 .* sqrt((0.075 + 0.25 .* (zeff -1) .^ 2) .* nue); % OK
term2  = 0.13 .* (1 - 0.38 .* ft) .* nue ./ zeff .^ 2 .* (term21 + term22); % OK
out    = ft ./ (1 + term1 + term2); % OK


function out = ft_eff_32_ei(zeff,ft,nue)

term1  = 0.87 .* (1 + 0.39 .* ft) .* sqrt(nue) ./ (1+ 2.95 .* (zeff -1) .^ 2); %OK
term2  = 1.53 .* (1 - 0.37 .* ft) .* nue .^ 2 .* (2 + 0.375 .* (zeff -1)); % OK
out    = ft ./ (1 + term1 + term2); % OK 


function out = ft_eff_33(zeff,ft,nue)

term1  = 0.25 .* (1 - 0.7 .* ft) .* sqrt(nue) .* (1+ 0.45 .* sqrt(zeff -1)); % OK
term2  = 0.61 .* (1 - 0.41 .* ft) .* nue ./ sqrt(zeff); %OK 
out    = ft ./ (1 + term1 + term2); % OK 


function out = falpha(zeff,ft,nui)

out = (falpha0(zeff,ft) + 0.7 .* zeff .* sqrt(ft .* nui)) ./  ...
      (1 + 0.18 .* sqrt(nui)) - 0.002 .* nui .^ 2 .* ft .^ 6;  % OK
out = out ./ (1 + 0.004 .* nui .^ 2 .* ft .^ 6); % OK


function out = falpha0(zeff,ft)

out = - (0.62  + 0.055 .* (zeff -1)) ./ (0.53 + 0.17 .* (zeff -1)) .* ...
        (1 - ft) ./ (1 - (0.31 - 0.065 .* (zeff -1)) .* ft - 0.25 .* ft .^ 2); %OK
    

function out = sigma_ratio(x33,zeff)

out = 1 - (1 + 0.21 ./ zeff) .* x33 + 0.54 ./ zeff .* x33 .^ 2 - ...
      0.33 ./ zeff .* x33 .^ 3; %OK
  

% from O. Sauter original paper + erratum
function out = sigma_spitzer(ne,te,zion)

nz  = 0.58 + 0.74 ./ (0.76 + zion);   % OK
out = 1.9012e4 .* te .^ (3/2) ./ zion ./ nz ./ lambda_ln_e(ne,te); % OK


function out = nue_star(q,R,epsi,ne,te,zion)

out = 6.921e-18 .* q .* R .* ne .* zion .* lambda_ln_e(ne,te) ./ ...
      te .^ 2 ./ epsi .^(3/2); % OK


function out = nui_star(q,R,epsi,ni,ti,zion)

out = 4.90e-18 .* q .* R .* ni .* zion .^ 4 .* lambda_ln_ii(zion,ni,ti) ./ ...
      ti .^ 2 ./ epsi .^ (3/2); % OK
  

function out = lambda_ln_e(ne,te)

out = 31.3 - log(sqrt(ne)./te); % OK


function out = lambda_ln_ii(zion,ni,ti)

out = 30 - log(zion .^ 3 .* sqrt(ni)./ ti .^ (3/2)); % OK



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
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
