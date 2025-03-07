% calcul d'une aproximation de bpol
Sp = data.equi.spr(:,end);
K  = data.equi.e(:,end);
ka = data.gene.volume ./ (2 .* pi .^ 2 .* data.equi.rmoy(:,end) .* a.^ 2);
a  = data.equi.a(:,end);
x  = linspace(0,1,101);
jli = data.equi.jmoy;
spr   = (2 .* Sp) * x;
bpol       = cumtrapz(x,4.* pi .* 1e-7 .* jli .*spr,2) ./ (2.*pi.* (a*x) + eps);
bpol(:,1) = 0;
% correction de l'elongation
%lini       = 2 .* trapz(x,bpol.^ 2 .* (a *x),2) ./ a  ./ bpol(:,end) .^ 2;
lini  =  trapz(x,bpol.^ 2 .* (a *x),2) ./ bpol(:,end) .^ 2 ./ a  .* (1 + K .^ 2)  ./ ka .^ (3/2) ;

figure(18);plot(data.gene.temps,lini,data.gene.temps,data.gene.li);
