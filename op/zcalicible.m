function [ilh,fgre] = zcalicible(data) 
%
% efficacite LH
%--------------
cst    = 0.75*sqrt(2.1+3)/4*7;
nbar   = data.gene.nbar;
b0     = data.geo.b0;
zeffm  = data.gene.zeffm;
cons   = data.cons.hyb(:,1);
r0     = data.geo.r0;
a      = data.geo.a;
iboot  = data.gene.iboot;
eta    = 6.5e18;
ilh    = eta .* abs(cons) ./ r0 ./ nbar;
fgre   = (nbar/1e20)./((iboot+ilh)/1e6)*3.1415.*(a.^2);
