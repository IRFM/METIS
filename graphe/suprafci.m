function [sup,frsup]=suprafci(data,param)
% -- calcul de la fraction de suprathermique
% [sup,frsup]=suprafci(data,param)
%
x = param.gene.x;

for k=1:length(data.gene.temps)
   Psup      = data.source.fci.psupra(k,:);
   sup(k)    = zintvol(Psup,x,data.equi.vpr(k,:),data.equi.rhomax(k))/1e6;
   frsup(k)  = sup(k)/(data.gene.paddfci(k)/1e6)* 100; 
end


function s=zintvol(e,x,vpr,rhomax)   

  s = rhomax.*trapz(x,vpr .* e,2);
