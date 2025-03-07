function [psi,dpsidt,q,epar,F,bpol,li] = z0j2q(t,x,jmoy,ptot,spr,r0,b0,rm,peri,ip,vloop)

% on travail a ip donne pour le 0D
% definition
mu0 = 4 .* pi .* 1e-7;
ve = ones(size(x));
vt = ones(size(t));

% calcul de psi
psi   =  j2psi(r0,rm,x,jmoy); 

% valeur au bord
dpsidt_edge  = - vloop ./ 2 ./ pi;
psi_edge     = cumtrapz(t,dpsidt_edge,1);
psi  = (psi - psi(:,end) * ve) +  psi_edge * ve;
dpsidt = gradient(psi')' ./ (gradient(t) * ve);

% calcul de F
dpsidx   = pdederive(x,psi,0,2,2,1);
df2dpsi  = ((2.* mu0 .* r0 ) * ve) .* (jmoy - (r0*ve) .* pdederive(x,ptot,0,2,2,1) ./ dpsidx);
df2dx    = df2dpsi .* dpsidx;
F2       = (r0 .^ 2 .* b0 .^ 2) *ve + cumtrapz(x(end:-1:1),df2dx(:,end:-1:1),2);
F        = sqrt(F2(:,end:-1:1));

% calcul de epar
epar = - F ./ ((b0 .* r0 .^ 2) *ve) .* dpsidt;

% calcul de q
q = - F .* spr ./ (r0 *ve) ./ 2 ./ pi ./ dpsidx .* (rm * ve); 

% calcul de li
bpol  = abs(dpsidx) ./ ((r0 .* rm) * ve);
li  =  2 .* rm .*  trapz(x,bpol.^ 2 .* spr .* ((2 .* pi .* r0)*ve),2) ./ ((mu0 .*ip ) .^ 2  .* r0);





function psi = j2psi(r0,rm,x,j)

mu0 = 4 .* pi .* 1e-7;
vt  = ones(size(r0));
ve  = ones(size(x));
inte1 = ((mu0 .* r0  .* rm .^ 2) * x ) .* j;
inte2 = cumtrapz(x,inte1,2) ./ (vt *x + eps);
inte2(:,1) = 0;
psi   = cumtrapz(x(end:-1:1),inte2(:,end:-1:1),2);
psi   = - psi(:,end:-1:1);
psi   = psi - (psi(:,end) * ve); 
 
