% cross section for p + B11 ->  3 * alpha
% reference: The thermonuclear fusion rate coefficient for p-11B reactions
% W.M. Nevins and R. Swain 2000 Nucl. Fusion 40 865, 
function [sigma,S1,S2,S3,S] = pb11_cross_section_nevins(E_eV)

% physical constants
phys = cphys;

% formula in keV
E_keV = E_eV / 1000;

% Gamov energy
E_G = 22.589e6; %eV

% separation in 3 energy interval for S(E)
% from  near 0 to 400 keV
C0 = 197e6; % eV
C1 = 0.240e6;
C2 = 2.31e2;
Al = 1.82e10;
El = 148e3;
delta_El = 2.35e3;
%
S1 = C0 + C1 .* E_keV + C2 .* E_keV .^ 2 + ...
     Al ./ ((E_keV - El ./ 1e3) .^ 2 + (delta_El ./ 1e3) .^ 2);

% from 400 keV to 642 keV
D0 = 330e6; %eV
D1 = 66.1e6;
D2 = -20.3e6;
D5 = -1.58e6;
XE = (E_eV - 400e3) ./ 1e5;
S2 = D0 + D1 .* XE + D2 .* XE .^ 2 + D5 .* XE .^ 5;

% from 642 keV to Inf 
Ai  = [2.57e12   5.67e11     1.34e11     5.68e11]; % eV
Ei  = [581.3e3   1083e3      2405e3      3344e3]; % eV
dEi = [85.7e3    234e3       138e3       309e3]; % eV
B   = 4.38e6;
S3  = zeros(size(E_eV));
for k=1:length(Ai)
   S3 = S3 +  Ai(k) ./ (((E_eV -Ei(k)) / 1e3) .^ 2 + (dEi(k) ./ 1e3) .^2);
end
S3 = S3 + B;

% selection
S = double(E_eV < 400e3) .* S1 +  ...
    double((E_eV >= 400e3) & (E_eV < 642e3)) .* S2 + ...
    double(E_eV >= 642e3) .* S3;

% cross section (m^2)
sigma = 1e-28 .* S ./ E_eV .* exp(-sqrt(E_G ./ E_eV));



