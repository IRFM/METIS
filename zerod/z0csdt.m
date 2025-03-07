% section efficace T(d,n)He4 toutes energies
% utilisee pour l'interraction faisceau plasma
function cs = z0csdt(ecm)
%ref : improved formulas for fusion cross-sctions and thermal reactivities,
% H-S Bosch and G.M. Hale, NF vol 32, p 611- , 1992.
% ecm en keV et cs en 1e-3 barn , soit 1e-3 * 1e-4 * 1e-24 m ^ 2 = 1e-31 m ^ 2

% ecm  <550 keV
a1 = 6.927e4; a2 = 7.454e8;  a3 = 2.05e6;   a4 = 5.2002e4; a5 = 0;
b1 = 6.38e1;  b2 = -9.95e-1; b3 = 6.981e-5; b4 = 1.728e-4;
sse1 = (a1 + ecm .* (a2 + ecm .* (a3 + ecm .* (a4 + ecm .* a5)))) ./ ...
       (1 + ecm .* (b1 + ecm .* (b2 + ecm .* (b3 + ecm .* b4))));
% de 550 a 4700
sse2 = -1.4714e6 ./(1 + ecm .* (-8.4127e-3 + ecm .* (4.7983e-6 + ecm .* (-1.0748e-9 +ecm .* 8.5184e-14))));
sse = sse1 .* (ecm <= 550) + sse2 .* (ecm >550);
cs = sse ./ ecm ./exp(34.3827./sqrt(ecm));


