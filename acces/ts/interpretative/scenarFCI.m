function [scenar,scenstr]=scenarFCI(geo,phys,freq,pos,composition)
%
%  [scenar,scenstr]=scenarFCI(geo,phys,freq,pos,composition);
%  
%  scenario FCI
% -- entrees --
%  freq          : frequence (MHz) 
%  pos           : position des antennes par rapport a la derniere surface magnetique
%  geo           : structure contenant les donnees geometriques du plasma
%  composition   : structure contenant la composition du plasma
%  phys          : grandeur physique (me, el, mp)
%
% -- sorties --
% scenar :  	1  -> minoritaire ou second harmonique 	-> utilisation de pion
% 		0  -> FWEH 				-> utilisation d'absor
% 		-1 -> probleme
% scenstr : 	HMIN_H
%		HMIN_D 
%		HMIN_He3 
%		HMIN_He
%               FWEH
%               HARM_2H 
%
%  Création : 11 février 2002
%  Auteur   : V. Basiuk
%


mp    = phys.mp;
el    = phys.e;
me    = phys.me;
%
% champ B
%
B0    = mean(geo.b0);
R     = linspace(mean(geo.r0)-mean(geo.a)-0.7,mean(geo.r0)+mean(geo.a)+0.7,500);
B     = B0*2.37./R;
fp    = el*B/2/pi/mp/1e6;
fd    = el*B/2/pi/mp/2/1e6;
fhe3  = 2*el*B/2/pi/mp/3/1e6;

f    = freq;
nh   = 1;
nd   = 1;
nhe3 = 1;
rH   = [0 0 0];
rD   = [0 0 0];
rHe3 = [0 0 0];
%
% position grossiere des couches
%
for k=1:3
%
% hydrogene
%
   ind            = find(f/k>=fp);  
   if isempty(ind) == 0
     ind         = ind(1);
     rH(nh)      = R(ind);
     nh          = nh+1;
   end
%
% deuterium ou helium 4
%
   ind           = find(f/k>=fd); 
   if isempty(ind) == 0 
      if ind(1) > 1
        ind      = ind(1); 
        rD(nd)   = R(ind);
        nd       = nd+1;
      end
   end
%
% helium3
%
   ind           = find(f/k>=fhe3); 
   if isempty(ind) == 0 
      if ind(1) > 1
        ind        = ind(1); 
        rHe3(nhe3) = R(ind);
        nhe3       = nhe3+1;
      end
   end

end
%
% choix du scenario
% 1  : minoritaire ou second harmonique 	-> utilisation de pion
% 0  : FWEH 				-> utilisation d'absor
% -1 : probleme
scenar  = -1;
scenstr = '????????';
%
% nature du minoritaire
%
if composition.z(2) == 1 &  composition.a(2) == 1
  gaz = 'H';
end
if composition.z(2) == 1 &  composition.a(2) == 2
  gaz = 'D';
end
if composition.z(2) == 2 &  composition.a(2) == 3
  gaz = 'He3';
end
if composition.z(2) == 2 &  composition.a(2) == 4
  gaz = 'He';
end
%
% position des couches par rapport a l'axe magnetique, normalise au petit rayon
%
for k=1:3
  distH(k)   = abs((rH(k)-mean(geo.r0)))./mean(geo.a);
  distD(k)   = abs((rD(k)-mean(geo.r0)))./mean(geo.a);
  distHe3(k) = abs((rHe3(k)-mean(geo.r0)))./mean(geo.a);
end
%
% grand rayon max du plasma [distp] et grand rayon min [distm]
%
distp  = mean(geo.r0+geo.a);
distm  = mean(geo.r0-geo.a);
centre = mean(geo.r0);
%
% determination du scenario
%
if strcmp(gaz,'H') & distH(1) > 0.5 & rH(1) < centre & B0 < 3

  scenar = 0;
  scenstr = 'FWEH    ';

              
end
	     	
if strcmp(gaz,'H') & distH(1) < 0.7

  scenar = 1;
  scenstr = 'HMIN_H  ';

              
end
if strcmp(gaz,'D') & distD(1) < 0.7

  scenar = 1;
  scenstr = 'HMIN_D  ';

              
end
if strcmp(gaz,'He') & distD(1) < 0.7 

  scenar = 1;
  scenstr = 'HMIN_He ';

end
if strcmp(gaz,'He3') & distHe3(1) < 0.7  

  scenar = 1;
  scenstr = 'HMIN_He3';

end
if ~strcmp(gaz,'He3') & distH(1) > 0.7 & distH(2) > 0.7

  scenar  = 0;
  scenstr = 'FWEH    ';

end
if strcmp(gaz,'H') & distH(1) > 0.7 & distH(2) < 0.7

  scenar  = 1;
  scenstr = 'HARM_2H ';

end

	   




