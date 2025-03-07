% section efficace T(d,n)He4 toutes energies
% utilisee pour l'interraction faisceau plasma
function cs = z0csdhe3(ecm)
%ref : improved formulas for fusion cross-sctions and thermal reactivities,
% H-S Bosch and G.M. Hale, NF vol 32, p 611- , 1992.
% ecm en keV et cs en 1e-3 barn , soit 1e-3 * 1e-4 * 1e-24 m ^ 2 = 1e-31 m ^ 2

% ecm  0.3 a 900 keV
a1 = 5.7501e6;    a2 = 2.5226e3;    a3 = 4.5566e1;   a4 = 0.0;  a5 = 0.0;
b1 = -3.1995e-3;  b2 = -8.5530e-6;  b3 = 5.9014e-8;  b4 = 0.0;
sse1 = (a1 + ecm .* (a2 + ecm .* (a3 + ecm .* (a4 + ecm .* a5)))) ./ ...
       (1 + ecm .* (b1 + ecm .* (b2 + ecm .* (b3 + ecm .* b4))));
% de 900 a 4800
sse2 = -8.3993e5 ./(1 + ecm .* (-2.6830e-3 + ecm .* (1.1633e-6 + ecm .* (-2.1332e-10 +ecm .* 1.4250e-14))));

sse = sse1 .* (ecm <= 900) + sse2 .* (ecm >900);
cs = sse ./ ecm ./exp(68.7508./sqrt(ecm));


