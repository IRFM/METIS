% test correction equipartition
te = logspace(0,6,1001);
ti = te;
Ai = 1;
Zi = 1;

% physical constantes
phys = cphys;

tstar = phys.e .* te ./ phys.me ./ phys.c .^ 2; 
%
%epsi_low = 3/2 + 15/8 .* tstar - 15/8 .* tstar .^ 2;
%epsi_high= 3 - 1./ tstar + 1 ./ 2 ./ tstar .^ 2;
epsi  = (besselk(1, 1./ tstar) ./ besselk(2, 1 ./ tstar) + 3 .* tstar - 1) ./ tstar;

% normalisation to be used with qei
% il y a sans doute une inversion
% le temps d'equipartition augmente avec te relativist dans le papier 
% donc la puissance echangée diminue avec ce facteur.
rc = (3./2) ./ epsi;


ten = phys.e .* te ./ phys.me ./ phys.c .^ 2;
tin = phys.e .* ti ./ (Ai .* phys.ua) ./ phys.c .^ 2;

ex_sp  = (ten + tin) .^ -(3/2);
ex_rel = exp(-1 ./ ten) ./ besselk(2,1 ./ ten) .* sqrt(ten ./ (ten + tin) .^ 3) .* ...
         (2 .* (ten+tin) .^ 2 + 2 .* (ten + tin) +1);
     
ec = sqrt(2 .* pi) ./ 2 .* ex_rel ./ ex_sp; 


cordey = 1 + 0.3 .* ten;

figure;
subplot(2,2,1)
plot(te,rc);

subplot(2,2,2)
plot(te,ec);

subplot(2,2,3)
plot(te,rc .* ec,'b',te,cordey,'r',te,ec);






