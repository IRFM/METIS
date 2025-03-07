% elagrissement du depot du fait de la largeure d'orbite
%  x        = coordonnees
%  rhomax   = equi.rhomax
%  energie  = energie des particule (eV)
%  aj       = masse des particule
%  zj       = charge des particule
%  rmoy     = rayon moyen
%  q        = facteur de securite
%  btot     = champ magnetique total
%  ftrap    = fraction de trapper
%  vpr      = element de volume
%  source   = depot sans elargissement
%  phys     = structure phys

function so = zorbitel(x,rhomax,energie,aj,zj,rmoy,bpol,ftrap,vpr,source,phys)

phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)

vth      = sqrt(2 .* phys.e .* energie ./ aj ./ phys.mp );
dr       = aj .* phys.mp ./ zj ./ phys.e .* rmoy .* vth ./ 4 ./ rmoy(end) ./bpol; 
dr(:,1)    = dr(:,2);

w        = min(1/4,dr ./ rhomax);

ve       = ones(size(x));
dd       = ((1 ./ sqrt(2 .* pi .* w' .^ 2)) * ve)' .* ...
           exp(- ( (x' *  ve - ve' * x) .^ 2 ./ (w' * ve)' .^ 2)); 

norme   = ve' * trapz(x',dd,1);
so      = trapz(x,dd ./ norme .* (ve' * source),2)';
so      = ftrap .* so + (1 - ftrap) .* source; 


eout    = rhomax .* trapz(x,vpr .* so);
ein     = rhomax .* trapz(x,vpr .* source);
so      = so .* ein ./ eout;


