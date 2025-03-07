% function computing prediction from simple model for Te Ti depending of the heating scheme
% all input are in SI unit but Te and Ti in eV
%
% syntax:
%  [Te_x1,Ti_x1,Pei,error_best] = analytical_teti_model(scaling,Padd,fpadd_ion,eta,R,a,K,Bt,Ip,ne,Zeff,Zmain,Zimp,Amain,Aimp,x1,[expo_tite,expo_fgr,aolti_max,x2])
%
% input:
%    scaling   = energuy confinement time scaling law name  ('ITERH-96P(th)','ITERH-98P(y,2)','Sauter-Martin','ITPA20','neo Alcator','Wesson','Lukash')
%    Padd      = additional power (W); for 'Lukash' = P_aux + sqrt(-1) * P_OH
%    fpadd_ion = fraction of additional power heating directly the ion
%    eta       = Chi_i / Chi_e  or factor in front of (Ti/Te)^expo_tite * f_gr ^ expo_fgr
%    R         = major radius (m)
%    a         = minor radius (m)
%    K         = plasam elongation
%    Bt        = magnetic field (T)
%    Ip        = plasma current (A)
%    ne        = electon density (m^-3) or Greenwald fraction
%    Zeff      = effective charge
%    Zmain     = main ion number of charge
%    Zimp      = impurities number of charge
%    Amain     = main ion number of mass
%    Aimp      = impurities number of mass
%    x1        = normalized Lao radius inside which additional power is deposited
%    expo_tite = exponent of denpendance of Chi_i/Chi_e in Ti/Te:
%                 Chi_i/Chi_e = (Ti/Te)^expo_tite
%                 default = 0
%    expo_fgr  = exponent of dependance of Chi_i/Chi_e in Greenwald fraction (default = 0); 
%    aolti_max = maximal normalized Ti gradient (~2, default = Inf)
%    x2        = optionnal integration radius for full analytical computation (only for test)
%                if NaN, compute numericaly the integral.
%
% output:
%   Te_x1      = predicted electron temperature in x1 (eV)
%   Ti_x1      = predicted ion temperature in x1 (eV)
%   Pei        = power transfered from electron to ion by collision (W)
%   error_best = model convergence error
%   
% test:
%
% padd =  linspace(0e6,10e6);
% [Te,Ti,Pei,error_best] = analytical_teti_model('ITERH-96P(th)',padd,0,NaN,2.5,0.5,1.5,3.7,5e5,0.5,2,1,7,2,14,0.2);
% figure;
% subplot(2,2,1)
% plot(padd/1e6,Te,'r',padd/1e6,Ti,'b')
% subplot(2,2,2)
% plot(padd/1e6,Ti./ Te)
% subplot(2,2,3)
% plot(padd/1e6,Pei/1e6)
% subplot(2,2,4)
% plot(padd/1e6,error_best)

function [Te_out,Ti_out,Pei_out,error_best] = analytical_teti_model(scaling,Padd,fpadd_ion,eta,R,a,K,Bt,Ip,ne,Zeff,Zmain,Zimp,Amain,Aimp,x1,expo_tite,expo_fgr,aolti_max,x2)

% test input
if nargin < 17
    expo_tite = 0;
else
    expo_tite = real(expo_tite);
end
if nargin < 18
    expo_fgr = 0;
end
if nargin < 19
    aolti_max = Inf;
end
if nargin < 20
    x2 = NaN;
end


% constants
phys = cphys;

% test input 
if ne < 10
    disp('Greenwald fraction is get from input')
    % this is a fraction of Greenwald density
    ne = ne .* 1e20 .* (Ip ./1e6) ./ pi ./ a.^2
end
% Volume element
Cv = 2*pi^2.*R.*a.^2.*K;

% plasma composition
Cmain = (Zimp - Zeff) ./ (Zimp - Zmain) ./ Zmain;
Cimp  = (Zeff - Zmain) ./ (Zimp - Zmain) ./ Zimp;
Ci    = (Zimp + Zmain -Zeff) ./ Zmain ./ Zimp;

% decode
Pohm = imag(Padd);
Padd = real(Padd) + imag(Padd);

% scaling law (term beta et alpha)
qcyl  = 5 .* K .* a .^ 2 .* Bt ./ R ./ (Ip / 1e6);
qcyl_star  = 5 .* (1 + K .^ 2) ./ 2 .* a .^ 2 .* Bt ./ R ./ (Ip / 1e6);
Hg    = (ne/1e20) .* pi .* K .* a .^ 2 ./ (Ip / 1e6);
nsat  = 0.06 .* (Ip / 1e6) .* R .* sqrt(Amain) ./ K ./ a .^ 2.5;
nstar = min( nsat,(ne/1e20)); 
% transition d'apres J.G. Cordey rapport JET-P(85)28
rap   = min(1,max(0,erfc(max(0,Padd) ./ max(1,Pohm) - 1)));



d95   = 0;
Vp    = Cv;
if isnumeric(scaling)
    alpha = 2/3;
    beta  = scaling ./ max(1,Padd)  .^ (1 - alpha);
else    
    switch scaling
        case 'ITERH-96P(th)'
            alpha = 0.73;
            beta  = 23e-3  .* (Ip ./ 1e6) .^ 0.96 .* Bt .^ 0.03 .* (ne/1e19) .^ 0.4 .* (1e-6) .^ -0.73 .* ...
                R .^ 1.83 .* K .^ 0.64 .* (a./R) .^ -0.06 .* Amain .^ 0.2;
            
        case 'ITERH-98P(y,2)'
            alpha  = 0.69;
            beta   = 56.2e-3  .* (Ip /1e6) .^ 0.93 .* Bt .^ 0.15 .* (ne/1e19) .^ 0.41 .* (1e-6) .^ -0.69 .* ...
                R .^ 1.97 .* K .^ 0.78 .* (a/R) .^ 0.58 .* Amain .^ 0.19;
        case 'Sauter-Martin'
            alpha  = (2/3);
            beta   = 0.0307 .* a .* K .* Bt .* sqrt(ne/1e19) ./ qcyl .* (1e-6 ./ Vp) .^ -(2/3);
        case 'ITPA20'
            alpha  = 0.669;
            beta   = 0.053 .* (Ip/1e6) .^ 0.98 .* Bt .^ 0.22 .* R .^ (1.71-0.35) .* (ne / 1e19) .^ 0.24 .* ...
                a .^ 0.35 .* K .^ 0.8 .* Aimp .^ 0.2 .* (1e-6) .^  (-0.669) .* (1 + abs(d95)) .^ 0.36;
            
        case 'neo Alcator'
            alpha = 0;
            beta  = 0.07 .* (ne/1e20)  .* a .* R .^ 2 .* qcyl_star .* sqrt(Amain);
            
        case 'Wesson'
            alpha = 0;
            beta  = 0.07 .* nstar .* a .* R .^ 2 .* qcyl;
            
        case 'Lukash'
            
            alpha = 0.73;
            beta_l  = 23e-3  .* (Ip ./ 1e6) .^ 0.96 .* Bt .^ 0.03 .* (ne/1e19) .^ 0.4 .* (1e-6) .^ -0.73 .* ...
                      R .^ 1.83 .* K .^ 0.64 .* (a./R) .^ -0.06 .* Amain .^ 0.2;
            
            beta_oh = 0.14 .* a .* R .* Hg .* Bt .^ 1.1 ./ (1 + 0.63 .* Bt .* Hg .^ 2) .* Padd .^ alpha;
            beta = rap .* beta_oh + (1-rap) .* beta_l;
    end
end

% loop of convergence on mu
if isnumeric(scaling)
    mu_in = linspace(0.1,2,10001);
else
    switch scaling
        case 'Lukash'
            mu_in = linspace(0.1,2,10001);
        otherwise
            mu_in = logspace(-1,1,100001);
    end
end
x     = linspace(1,x1);
ve    = ones(size(x));
if ~isfinite(eta)
    eta_ref = 1;
else
    eta_ref = eta;
end
fprintf('Chi_i / Chi_e = %g * (Ti/Te)^{%g} * f_Gr^%g\n',eta_ref,expo_tite,expo_fgr);

if ~isfinite(x2)
    x2 = 0.7;
    mode_full = true;
    disp('full integration');
else
    mode_full = false;
end
if isfinite(aolti_max)
    fprintf('Clamping a/L_{Ti} at %g\n',aolti_max)
end

for k=1:length(mu_in)
    %%fprintf('.');
    % loop on mu
    mu_loc = mu_in(k);
    
    % compuation of Te in x1
    Te = 2 .* beta .* Padd .^ (1 - alpha) .* (1 -x1) ./ (phys.e .* Cv .* ne .* (1 + mu_loc .* Ci) .* (1 - x1 .^ 3));
    
    % Equiparttion constant term
    warning off
    lnei          =  14.9 - 0.5.*log(ne ./ 1e20) + log(Te ./ 1e3);
    warning on
    ind = find(~isfinite(lnei) | (lnei <10));
    if ~isempty(ind)
        lnei(ind) = 10 .* ones(1,length(ind));
    end
    Aei  = (Cmain .* Zmain .^ 2 ./ Amain + Cimp .* Zimp .^ 2 ./ Aimp) .* phys.e .^ (7/2) .* sqrt(2 .* phys.me) .* lnei ./ ...
           4 ./ pi .^ (3/2) ./ phys.mp ./ phys.epsi0 .^ 2;

    % ion temperature
    Pei_x2 = ne .^ 2 .* Aei .* (1 - mu_loc) ./ sqrt(Te) .* Cv .* (x1 .^ 2 + 4/3 .* ((x1 + 2) .* sqrt(1-x1) - (x2 + 2) .* sqrt(1-x2)));
    error_borne = 0 + 1e38 .* double(Pei_x2 <= (- 0.99 * fpadd_ion .* Padd)) + 1e38 .* double(Pei_x2 >= (0.99 .* (1 -fpadd_ion) .* Padd));
    Pei_x2 = max(- 0.99 * fpadd_ion .* Padd,min(0.99 .* (1 -fpadd_ion) .* Padd, Pei_x2));
    
    % main dependance from Asp paper (doi:10.1088/0741-3335/47/3/007)
    fgr = ne ./( 1e20 .* (Ip ./1e6) ./ pi ./ a.^2);
    eta = eta_ref * mu_loc^expo_tite .* fgr .^ expo_fgr;
    
    Ti  = max(13.6,Te ./ eta ./ Ci .* (fpadd_ion .* Padd + Pei_x2) ./ ( (1 -fpadd_ion) .* Padd - Pei_x2)); 
    if isfinite(aolti_max)
       for l=1:11
           dtidx = min(Ti ./ (1-x1),aolti_max.* Ti);
           Ti    = dtidx .* (1-x1);
       end
    end
    
    if mode_full
        xx = ones(size(Ti(:)))*x;
        v1 = - Te ./ eta ./ Ci ./ (1-x1);
        if length(v1) > 1
            v1 = v1(:) * ve;
        end
        v2 = ne .^ 2 .* Aei .* (1 - mu_loc) ./ sqrt(Te) .* Cv;
        if length(v2) > 1
            v2 = v2(:) * ve;
        end
        v3 = fpadd_ion .* Padd;
        if length(v3) > 1
            v3 = v3(:) * ve;
        end
        v4 = (1 - fpadd_ion) .* Padd;
        if length(v4) > 1
            v4 = v4(:) * ve;
        end
        v5 = x1;
        if length(v5) > 1
            v5 = v5(:) * ve;
        end
        v6 = Padd;
        if length(v6) > 1
            v6 = v6(:) * ve;
        end        
        dtidx = v1 .* max(v6/1e3,v3 + v2 .* (v5 .^ 2 + 4/3 .* ((v5 + 2) .* sqrt(1-v5) - (xx + 2) .* sqrt(1-xx)))) ./ max(v6/1e3,v4 - v2 .* (v5 .^ 2 + 4/3 .* ((v5 + 2) .* sqrt(1-v5) - (xx + 2) .* sqrt(1-xx))));
        Ti    = max(13.6,trapz(x,dtidx,2));
        if isfinite(aolti_max)
            for l=1:11
                if length(Ti) > 1
                    dtidx = max(dtidx,-aolti_max.* (Ti*ve));
                else
                    dtidx = max(dtidx,-aolti_max.* Ti);
                end
                Ti    = max(13.6,trapz(x,dtidx,2));
            end
        end
    end
    Pei = ne .^ 2 .* Aei .* (1 - mu_loc) ./ sqrt(Te) .* Cv .* (x1 .^ 2 + 4/3 .* (x1 + 2) .* sqrt(1-x1));
    mu_rep  = Ti(:) ./ Te(:);
    %mu_rep(~isfinite(mu_rep)) = 1e38;
    if length(Ti) ~= length(Te)
            Te  = Te * ones(size(Ti));
    end
    if length(Ti) ~= length(Pei)
            Pei  = Pei * ones(size(Ti));
    end
    if k== 1
        Te_out  = Te(:);
        Ti_out  = Ti(:);
        Pei_out = Pei(:);
        error_best = abs(mu_rep(:) - mu_loc) + error_borne(:);
        Te_out(error_borne > 0)  = NaN;
        Ti_out(error_borne > 0)  = NaN;
        Pei_out(error_borne > 0) = NaN;
    else
        mu_rep = mu_rep(:);
        Te  = Te(:);
        Ti  = Ti(:);
        Pei = Pei(:);
        error = abs(mu_rep - mu_loc) + error_borne(:);
        mask  = error < error_best;
        Te_out(mask)     = Te(mask);
        Ti_out(mask)     = Ti(mask);
        Pei_out(mask)    = Pei(mask);
        error_best(mask) = error(mask);
    end  
end



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
