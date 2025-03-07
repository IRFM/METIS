% calcul de la mutuel pour le flux du champ vertical et la self externe
% model ref :Hirshman and Neilson Phys. of Fluids 29 (3) 1986, p 792
function [Lext,MBv] = z0lextmutbv(R0,amin,K)

e = max(eps, amin ./ R0);
a = (1 + 1.81 .* sqrt(e) + 2.05 .* e) .* log(8 ./ e) - ...
    (2 + 9.25 .* sqrt(e) - 1.21 .* e);
b = 0.73 .* sqrt(e) .* (1 + 2 .* e .^ 4 - 6 .* e .^ 5 + 3.7 .* e .^ 6);
c = 1 + 0.98 .* e .^ 2 + 0.49 .* e .^ 4 + 1.47 .* e .^ 6;
d = 0.25 .* e .* (1 + 0.84 .* e - 1.44 .* e .^ 2);

Lext = a .* (1 - e) ./ (1 - e + b .* K);
MBv  = (1 - e) .^ 2 ./ ((1 - e) .^ 2 .* c + d .* sqrt(K));

mu0 = 4*pi*10^(-7);
Lext = Lext .* mu0 .* R0;

