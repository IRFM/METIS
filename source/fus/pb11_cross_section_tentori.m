% cross section for p + B11 ->  3 * alpha
% reference: Revisiting p-11B fusion cross section and reactivity, and their analytic approximations
% Alessandro Tentori and Fabio Belloni 2023 Nucl. Fusion 63 086001
%
% and
% Low-Energy Cross Sections for 11B(p, 3e)*
% H.W. Becker, C. Rolfs, and H.P. Trautvetter
% Z. Phys. A - Atomic Nuclei 327, 341-355 (1987)
%
function [sigma,S1,S2,S3,S] = pb11_cross_section_tentori(E_eV)

% physical constants
phys = cphys;

% formula in keV
E_keV = E_eV / 1000;
E_MeV = E_keV / 1000;

% Gamov energy
E_G = 22.589e6; %eV

% separation in 3 energy interval for S(E)
% from  near 0 to 400 keV
C0 = 197e6; % eV
C1 = 0.269e6;
C2 = 2.54e2;
Al = 1.82e10;
El = 148e3;
delta_El = 2.35e3;
%
S1 = C0 + C1 .* E_keV + C2 .* E_keV .^ 2 + ...
     Al ./ ((E_keV - El ./ 1e3) .^ 2 + (delta_El ./ 1e3) .^ 2);
% original fit from measure for the resonance at 148 keV for testing
% Sal0 = 2.1e6 - 1.26e6 .* E_MeV  - 0.14e6  .* E_MeV .^ 2 + ...
%     0.69e3 ./ ((E_MeV -0.1483) .^ 2 + 7.13e-6);
% Sal1 = 195e6 + 241e6 .* E_MeV  + 231e6  .* E_MeV .^ 2 + ...
%     1.76e4 ./ ((E_MeV -0.148) .^ 2 + 5.52e-6);
% S1low = Sal0 + Sal1;    

% from 400 keV to 668 keV
D0 = 346e6; %eV
D1 = 150e6;
D2 = -59.9e6;
D5 = -0.46e6;
XE = (E_eV - 400e3) ./ 1e5;
S2 = D0 + D1 .* XE + D2 .* XE .^ 2 + D5 .* XE .^ 5;

% from 668 keV to 9760 keV
Ai  = [1.98e12   3.89e12     1.36e12     3.71e12]; % eV
Ei  = [640.9e3   1211e3      2340e3      3294e3]; % eV
dEi = [85.5e3    414e3       221e3       351e3]; % eV
B   = 0.381e6;
S3  = zeros(size(E_eV));
for k=1:length(Ai)
   S3 = S3 +  Ai(k) ./ (((E_eV -Ei(k)) / 1e3) .^ 2 + (dEi(k) ./ 1e3) .^2);
end
S3 = S3 + B;

% selection
S = double(E_eV < 400e3) .* S1 +  ...
    double((E_eV >= 400e3) & (E_eV < 668e3)) .* S2 + ...
    double(E_eV >= 668e3) .* S3;
S(E_eV > 9.76e6) = NaN;

% cross section (m^2)
sigma = 1e-28 .* S ./ E_eV .* exp(-sqrt(E_G ./ E_eV));



