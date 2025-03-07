function [piqa,piqq]=piquage(x,y,val,tot,v1max,v1min)
%____________________________________________________________________________|
% piquage (1-x2^piqa)^piqq, [piqa,piqq]=piquage(x,y,val,tot,v1max,v1min)     |
% x : rho normalise                                                          |
% y : profil                                                                 |
% piqq [0.7 5.5]                                                             |
% val = 1 -> piqa = 1, sinon piqa [0.4,2.5]                                  |
% v1max: valeur maximum de piqq [10.5], v1min: valeur minimum de piqq [0.7]  |
% tot : nombre de points dans chaque intervalle de recherche                 |
% V. Basiuk, A. Becoulet, 8 fevrier 1996                                     |
%____________________________________________________________________________|
[n,m]  = size(y);
l      = length(x);
if nargin < 6
   v1min = 0.7;
end
if nargin < 5 
   v1max = 10.5;
end
if nargin < 4
  tot = 50;
end
nb     = tot;
if nargin < 3
  val = 1;
end
if val == 1
  pasa = 1;
else
  pasa = linspace(0.4,3.5,nb);
end
pasq   = linspace(v1min,v1max,nb);

for kpiq=1:m
  fit     = y(:,kpiq);
  rexp    = x;
  profexp = fit/max(fit);
  psi     = rexp.^2;
  if val == 1
    piqj = zeros(1,nb);
    ni   = 1;
  else
    piqj = zeros(nb,nb);
    ni   = nb;
  end
  nfin = length(profexp);
  for it=1:ni
     alpha        = pasa(it);
     piqj(it,:)   = pasq;
     alj(it,:)    = alpha*ones(1,nb);
     val1         = (1-(psi'*ones(1,nb)).^alpha).^(ones(l,1)*pasq);
     prof         = (profexp(1)-profexp(nfin))*val1+profexp(nfin);
     meansq(it,:) = sum((prof-(profexp*ones(1,nb))).^2);
  end
  kk=find(meansq==min(min(meansq)));
  if isempty(kk) == 0
    kk=kk(1);
    piqq(kpiq) = piqj(kk);
    piqa(kpiq) = alj(kk);
  else
    piqq(kpiq) = 0;
    piqa(kpiq) = 0;
  end
end
