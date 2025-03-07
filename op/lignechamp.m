function [nnRRref,nnZZref,nBRref,nBZref,nBPref]=lignechamp(data,angf);   
% [nRRref,nZZref,BRref,BZref,BPref]=lignechamp(data,angf);
%
% inputs
% - angf : angle toroidal (en degre) sur lequel on veut les nouvelles surfaces de flux)
% angf > 0
%
% outputs
% - nnRRref,nnZZref : nouvelles surface de flux
% - BRref, BZref, BPref : nouvelles composantes du champ magnetique, 
%   prise en compte de l'effet du ripple sur les positions R Z et les composantes du champ
%   a theta et rho fixe.
%
geom(1) = data.geo.r0;
geom(2) = data.geo.a;
geom(3) = data.equi.ip;
geom(4) = data.equi.d(1);
geom(5) =  (data.geo.b0*data.geo.r0)/7.3008E-3;
%
% geom(9) = 1, champ poloidal analytique, 0, pas de champ poloidal
%
geom(9) = 0;
%
% piquage du courant 
%
x=data.equi.a/max(data.equi.a);
[pa,pq]=piquage(x,data.equi.jmoy',1);
geom(10) = pq;
disp('calcul de l''effet du ripple, patience ...')
RRref = squeeze(double(data.equi.R));
ZZref = squeeze(double(data.equi.Z));
BRref = squeeze(double(data.equi.BR));
BZref = squeeze(double(data.equi.BZ));
BPref = squeeze(double(data.equi.BPHI));
nRRref   = NaN .* ones(size(RRref));
nZZref   = NaN .* ones(size(RRref));
nnRRref  = NaN .* ones(size(RRref));
nnZZref  = NaN .* ones(size(RRref));
dBRref   = NaN .* ones(size(RRref));
dBZref   = NaN .* ones(size(RRref));
dBPref   = NaN .* ones(size(RRref));
ndBRref  = NaN .* ones(size(RRref));
ndBZref  = NaN .* ones(size(RRref));
ndBPref  = NaN .* ones(size(RRref));

[n,m] = size(RRref);
fprintf('\n');
for k=1:n
    fprintf('.');
   for l=1:m
     geom(6)  = RRref(k,l);
     geom(7)  = ZZref(k,l);
     geom(8)  = 0;
     geom(11) = angf*pi/180;
    [Rlig,Zlig,Plig,BRa,BZa,BPa]=lignemex(geom);
     nRRref(k,l)  = Rlig(end-1);
     nZZref(k,l)  = Zlig(end-1);
     dBPref(k,l)  = BPa(end-1)-BPa(1);
     dBRref(k,l)  = BRa(end-1)-BRa(1);
     dBZref(k,l)  = BZa(end-1)-BZa(1);
  end
  % reclacul des angles
  if k>1
    inda   = 1:(m-1);
    R0ref  = 0.5 .* (max(RRref(k,inda)) + min(RRref(k,inda)));
    Z0ref  = 0.5 .* (max(ZZref(k,inda)) + min(ZZref(k,inda)));
    angref = unwrap(angle((RRref(k,inda)-R0ref)+i*(ZZref(k,inda)-Z0ref)));


    nR0ref  = 0.5 .* (max(nRRref(k,inda)) + min(nRRref(k,inda)));
    nZ0ref  = 0.5 .* (max(nZZref(k,inda)) + min(nZZref(k,inda)));
    nangref = unwrap(angle((nRRref(k,inda)-nR0ref)+i*(nZZref(k,inda)-nZ0ref)));

    lnangref = unwrap(cat(2,nangref,nangref,nangref))-2*pi;
    lnRRref = cat(2,nRRref(k,inda),nRRref(k,inda),nRRref(k,inda));
    lnZZref = cat(2,nZZref(k,inda),nZZref(k,inda),nZZref(k,inda));
    lndBRref = cat(2,dBRref(k,inda),dBRref(k,inda),dBRref(k,inda));
    lndBZref = cat(2,dBZref(k,inda),dBZref(k,inda),dBZref(k,inda));
    lndBPref = cat(2,dBPref(k,inda),dBPref(k,inda),dBPref(k,inda));

    nnRRref(k,inda)=interp1(lnangref,lnRRref,angref,'cubic');
    nnZZref(k,inda)=interp1(lnangref,lnZZref,angref,'cubic');
    ndBRref(k,inda)=interp1(lnangref,lndBRref,angref,'cubic');
    ndBZref(k,inda)=interp1(lnangref,lndBZref,angref,'cubic');
    ndBPref(k,inda)=interp1(lnangref,lndBPref,angref,'cubic');

    nnRRref(k,m)= nnRRref(k,1);
    nnZZref(k,m)= nnZZref(k,1);
    ndBRref(k,m)= ndBRref(k,1);
    ndBZref(k,m)= ndBZref(k,1);
    ndBPref(k,m)= ndBPref(k,1);

  else
    nnRRref(k,:)=nRRref(k,:);
    nnZZref(k,:)=nZZref(k,:);
    ndBRref(k,:)=dBRref(k,:);
    ndBZref(k,:)=dBZref(k,:);
    ndBPref(k,:)=dBPref(k,:);
  end

end  
fprintf('\n');


nBRref = BRref + ndBRref;
nBZref = BZref + ndBZref;
nBPref = BPref + ndBPref;
