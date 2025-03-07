% fusion rate for p + B11 ->  3 * alpha
% reference: The thermonuclear fusion rate coefficient for p-11B reactions
% W.M. Nevins and R. Swain 2000 Nucl. Fusion 40 865, 
function [rate,rate_NRL,rate_NRH,rate_R] = pb11_fusion_rate_nevins(T_H_eV,T_B_eV)

% physical constants
phys = cphys;

% reduce masse in eV
% for pure B11 and proton
mu = 11.0093054 * 1.00782503207 / (11.0093054 + 1.00782503207) * phys.ua * phys.c ^ 2  / phys.e;
% for natural Boron and natural hydrogen 
%mu = 859.526e6;

% Gamov energy
%E_G = 22.589e6; %eV
% to be coherent
E_G = 2 .* mu .* (pi .* phys.alpha .* 1 .* 5) ^ 2;

% effective temperature
if nargin > 1
    T_eff_eV = (1.00782503207 .* T_H_eV + 11.0093054 .* T_B_eV) ./  ...
               (1.00782503207 + 11.0093054);
else
    T_eff_eV = T_H_eV;
end
if any(T_eff_eV(:) > 500e3)
    warning('pb11_fusion_rate_nevins: some T_eff are higher than 500 keV, fusion rate are maybe inaccurate!');
end

% low temperature (T <= 70 keV)
% non resonant part
E0       = (E_G/4) ^ (1/3) .* T_eff_eV .^ (2/3);
delta_E0 = 4 .* sqrt(T_eff_eV .* E0 ./ 3);
tau      = 3 .* E0 ./ T_eff_eV;
%
C0 = 197e6; % eV
C1 = 0.240e6;
C2 = 2.31e2;
Al = 1.82e10;
El = 148e3;
delta_El = 2.35e3;
%
S_eff = C0 .* (1 + 5 ./ 12 ./ tau) + C1 .* ((E0 + 35/36 .* T_eff_eV) ./ 1e3) + ...
        C2 .* ((E0 .^ 2 + 89/36 .* E0 .* T_eff_eV) ./ 1e6);
S_eff = 1e-28 .* S_eff; % go from barn to m ^2 
%
rate_NRL = phys.c .* sqrt(2 .* T_eff_eV ./ mu) .* delta_E0 .* S_eff ./ T_eff_eV .^ 2 .* exp(-tau);%
%
% resonant part
rate_R  = phys.c .* 1e-28 .* sqrt(8 .* pi .* T_eff_eV ./ mu) .* ...
          (1e3 ./ T_eff_eV) .^ 2 .*  Al ./ delta_El .* ...
          exp(-(sqrt(E_G ./ El) + El ./ T_eff_eV));
%rate_R_alt  = 5.41e-21 .* (1e3 ./ T_eff_eV) .^(3/2) .* exp(-148e3 ./ T_eff_eV);
rate_low = rate_NRL + rate_R;

% high temperature
%  mu is already defined
P1 = 4.4467e-11; % in eV m^3 / s
P2 = -5.9357e-5; % in eV^-1
P3 = 2.0165e-4;  % in eV^-1
P4 = 1.0404e-9;  % in eV^-2
P5 = 2.7621e-9;  % in eV^-2
P6 = -9.1653e-15; % in eV^-3
P7 = 9.8305e-16; % in eV^-3
%
nume  = P2 + T_eff_eV .* (P4 + T_eff_eV .* P6);
deno  = 1 + T_eff_eV .* (P3 + T_eff_eV .* (P5 + T_eff_eV .* P7));
theta = T_eff_eV ./ (1 - T_eff_eV .* nume ./ deno);
%
Xi = (E_G ./ 4 ./ theta) .^ (1/3);
rate_NRH = P1 .* theta .* sqrt(Xi ./ mu ./ T_eff_eV .^ 3) .* exp(- 3 .* Xi);
%
rate_high = rate_NRH + rate_R;
%
% raccord
%
f = min(1,max(0,(T_eff_eV - 50e3) / (70e3 -50e3)));
rate = (1 - f) .* rate_low + f .* rate_high;



