% utilisee pour l'interraction faisceau plasma D(d,n)He3
function cs = z0csdd(ecm)
%ref : improved formulas for fusion cross-sctions and thermal reactivities,
% H-S Bosch and G.M. Hale, NF vol 32, p 611- , 1992.
% ecm en keV et cs en 1e-3 barn , soit 1e-3 * 1e-4 * 1e-24 m ^ 2 = 1e-31 m ^ 2

% ecm  0.5 a 4900 keV
a1 = 5.3701e4; a2 = 3.3027e2;  a3 = -1.2706e-1;   a4 = 2.9327e-5; a5 = -2.5151e-9;
b1 = 0;  b2 = 0; b3 = 0; b4 = 0;
sse = (a1 + ecm .* (a2 + ecm .* (a3 + ecm .* (a4 + ecm .* a5)))) ./ ...
       (1 + ecm .* (b1 + ecm .* (b2 + ecm .* (b3 + ecm .* b4))));
cs = sse ./ max(eps,ecm) ./exp(31.3970./sqrt(max(eps,ecm)));
cs(ecm < 0.5) = 0;

