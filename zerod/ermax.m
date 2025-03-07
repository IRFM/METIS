function [er,tause] = ermax(te,ne,zeff,opt)

% constante physique (phys)
c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
epsi0       =   1./c.^2./mu0;             % permitivite du vide (F/m)  (definition)
me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)

ne = max(1e13,ne);
te = max(13.6,te);

lb = 15.2 - 0.5 .* log(ne./1e20) + log(te ./1e3);
er = e .^ 3 .* ne .*  lb  ./ (4 .* pi .* epsi0 .^ 2 .*  me .* c .^ 2);
er = er .* (zeff + 2);

if opt == 2
   er = er .* (me .* c .^ 2) ./ (e .* te);
end

vthe = sqrt(2.*e .*te ./me);
tause =  (4.*pi.*epsi0 .^ 2.*me .^ 2 .* vthe .^ 3 )./ (ne .* e .^ 4 .* lb);
