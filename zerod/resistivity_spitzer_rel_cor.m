function cor = resistivity_spitzer_rel_cor(ne,te,Zeff)
% this function comptute the correction to Spitzer resistitivity due to relativistic efects.
% reference:Transport coefficients of a relativistic plasma, O. J. Pike and S. J. Rose1, PHYSICAL REVIEW E 93, 053208 (2016)
% equation 57.
% Relativistic transport theory for a two-temperature magnetized plasma,
% T. Metens; R. Balescu, Phys. Fluids B 2, 2076-2090 (1990)
% https://doi.org/10.1063/1.859428


% with normalisation hidden in formula
cor = eta_norm_rel(ne,te,Zeff) ./ z0eta_spitzer(ne,te,Zeff) .*  ...
      z0eta_spitzer(ne,1,Zeff) ./ eta_norm_rel(ne,1,Zeff);
cor = cor ./ min(cor);


function out = eta_norm_rel(ne,te,Zeff)
% phyical constants
phys = cphys;

% normalized inverse energy
t = (phys.e .* te) ./ (phys.me .* phys.c .^ 2);

% coef 
Z  = [1 2 3 4 5 6 7 8 10 12 14 20 30 60 114];
a1 = [316 159 112 89.5 76.6 68.2 62.3 58.0 51.9 48.0 45.2 40.1 36.3 32.4 28.5];
%a2 = [12.0 10.2 9.36 8.88 8.57 8.35 8.18 8.04 7.85 7.71 7.61 7.42 7.25 7.08 6.89];
a3 = [109 66.3 51.7 44.3 39.9 36.9 34.7 33.0 30.7 29.2 28.0 26.0 24.4 22.7 21.0];
a4 = [639 279 185 143 120 105 95.2 87.8 77.7 71.2 66.6 58.6 52.4 46.4 40.4];
%a5 = [142 86.4 67.7 58.2 52.5 48.7 45.9 43.9 40.9 38.9 37.5 34.9 32.8 30.7 28.4];
a6 = [1390 573 373 287 240 210 190 175 156 143 134 118 107 95.3 83.9];
%a7 = [229 146 117 102 93.3 87.3 82.9 79.6 75.0 71.9 69.7 65.6 62.4 59.1 55.6];

% interpolation
a1i = interp1(Z,a1,Zeff,'pchip','extrap');
a3i = interp1(Z,a3,Zeff,'pchip','extrap');
a4i = interp1(Z,a4,Zeff,'pchip','extrap');
a6i = interp1(Z,a6,Zeff,'pchip','extrap');

% formula 57
alpha_c = 1 - (a1i + a3i .* t .* (1 + t)) ./ (a4i + a6i .* t .* (1+t));
%alpha_c = alpha_c ./ (1 -a1i ./ a4i);

% from alpha_c to alpha (equation 35-38)
warning off
lnei = 31.3 - log(sqrt(ne)./te); % OK
%lnei          =  15.2 - 0.5 .* log(ne ./ 1e20) + log(te ./ 1e3);
warning on
fei    = Zeff .* ne .* phys.e .^ 4 .* lnei./ (4 .* pi .* phys.epsi0 .^ 2 .* phys.me .^ 2);
K1     = besselk(1,1./t);
K2     = besselk(2,1./t);
x      = 1./t;
% eulergamma = 0.5772156649015328606;
%K2_low = 2./x.^2 -1/2 + (x.^2 .* (3 - 4 .* eulergamma + log(16) - 4 .*  log(x))) ./ 32 + (x .^ 4 .* (17 - 12 .*  eulergamma + log(4096) - 12 .* log(x))) ./ 1152 + ...
%         (x .^ 6 .* (43 - 24 .* eulergamma + 24 .* log(2) - 24 .* log(x))) ./ 73728; % + O[x]^7
%K2_high =  exp(-x) .* (sqrt(pi/2) .* sqrt(1./x) + 15/8  .* sqrt(pi/2) .* (1./x).^(3/2) + 105/128 .*  sqrt(pi/2) .* (1./x).^(5/2) - ...
%           (315 .* sqrt(pi/2) .* (1./x).^(7/2)) ./ 1024 + (10395 .* sqrt(pi/2) .* (1./x).^(9/2)) ./ 32768 -  ...
%           (135135 .* sqrt(pi/2) .* (1./x) .^ (11/2)) ./ 262144); %+ O((1/x)^(13/2)))
gav   = 3 .* t + K1 ./ K2;
rel_flag = (gav -1) > 1e-2;
gav(~rel_flag) = 1 + 1.5125 .* t(~rel_flag);
%tstar = fei .* gav .* t .* exp(-1./t) .* K2 ./ (1 + 2 .* t + 2 .* t .^2);
%out   = alpha_c ./ tstar .* gav
tstar = 3 .* phys.c .^ 3 .* gav .* t .* exp(1./t) .* K2 ./ (1 + 2 .* t + 2 .* t .^2) ./ fei;
tau = 3 .* sqrt(pi) .* phys.c .^ 3 .* t .^ (3/2) ./ sqrt(2) ./ fei;
tstar(~isfinite(tstar)) = 1.0021 .* tau(~isfinite(tstar)); % for low temperature asymtotic behaviour
out   = alpha_c ./  tstar .* gav .* ne .* phys.me; 




