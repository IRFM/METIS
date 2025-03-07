function [nRRref,nZZref]=lignemex1(data,itor,RRref,ZZref,angd,angf);   
% [nRRref,nZZref]=lignemex1(data,itor,RRref,ZZref,angd,angf);
%
% inputs
% - itor : courant toroidal
% - RRref et ZZref, : position des surfaces de flux en R et Z pour l'angle
% toroidal angd (en degre) donne
% - angf : angle toroidal (en degre) sur lequel on veut les nouvelles surfaces de flux)
% angf > angc
%
% outputs
% - nRRref,nZZref : nouvelles surface de flux
%
%
geom(1) = data.geo.r0;
geom(2) = data.geo.a;
geom(3) = data.equi.ip;
geom(4) = data.equi.d(1);
geom(5) = itor;
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
disp('prise en compte du ripple, patience !!!')
[n,m]=size(RRref);
fprintf('\n');
for k=1:n
    fprintf('.');
   for l=1:m
     geom(6) = RRref(k,l);
     geom(7) = ZZref(k,l);
     geom(8) = angd*pi/180;
     geom(11) = angf*pi/180;
    [Rlig,Zlig,Plig]=lignemex(geom);
%     ind          = min(find(Plig*180/pi>=angf));
     nRRref(k,l)  = Rlig(end-1);
     nZZref(k,l)  = Zlig(end-1);
  end
end  
fprintf('\n');
