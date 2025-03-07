function [val,val0,val1]=coucherip(R,Z,B0,freq,n,ep,A);
%
% function rip=coucherip(R,Z,B0,freq);
% minimum d'une couche FCI avec le ripple
%
% freq en MHz
% n  harmonique
% ep = +- 1 (1+ep*ripple)
% A masse atomique
% 

mp   = 1.6726e-27;
el   = 1.6022e-19;
me   = 9.1095e-31;
cst  = el*2.37*B0/2/pi/mp/A;
freq = freq*1e6;


X   = R-2.04;
Y   = 0.52*X+1;
b   = 0.26;
b2  = 2*b^2;
b4  = 4*b^2;
rc  = sqrt( ( Y - sqrt( Y.*Y - b4*(X.*X+Z.*Z) ) ) /b2);
rip = (1+ep*2.2e-4*exp(rc.*(5+rc*1.6)));

val = cst*rip./R-freq/n;


B    = B0*2.37./R;
fp   = el*B/2/pi/mp/1e6/A;
ind  = find(freq/n/1e6 >= fp);
if ~isempty(ind) 
  val1 = R(ind(1));
end

if sign(min(val)) == -1
  ind  = find(val <=0);
  ind  = ind(1);
  val0 = R(ind);
else
  rip  = 2.2e-4*exp(rc.*(5+rc*1.6));
  while sign(min(val)) == 1
    rip = rip/1.1;
    val = cst*(1+ep*rip)./R-freq/n;
  end
  ind  = find(val <=0);
  ind  = ind(1);
  val0 = R(ind);
end
