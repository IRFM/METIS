function rap = ratio_sigma_stop_B_vs_H(E_eV)

% no stop section formula exist for neutral boron beam
% some data are available in:
% Ionization, stopping, and thermalization of hydrogen and boron beams injected in fusion plasmas
% Agustin F. Lifschitz; Ricardo Farengo; Nestor R. Arista
% Phys. Plasmas 7, 30363041 (2000)  https://doi.org/10.1063/1.874156

% Assuming other dependance are the same, use data from this paper to
% compute sigma_B / sigma_H = f(E_ev) from attenuation length of H^0 and
% B^0 from 1 tp 0.2 of state fraction computed in 20 % concentration of
% boron in hydrogen plasma of 10^21 m^-3 density. The dependance in plasma
% temperature is small

% input
lE_eV = log(max(1,E_eV));

% data
E_in = [200e3,  640e3]; % eV
L_H  = [7.8e-2, 3.8e-2]; %m
L_B  = [0.5,   0.15]; %m

% fit
pp_H = polyfit(log(E_in),log(L_H),1);
pp_B = polyfit(log(E_in),log(L_B),1);

% ratio
LBE = exp(polyval(pp_B,lE_eV));
LHE = exp(polyval(pp_H,lE_eV));
% factor by which stop cross section should be multiplied:
rap = min(1,max(0.1,LHE ./ LBE));


