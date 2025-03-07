function [scenar,scenstr,Rres]=scenarzerodfci(geo,freq,pos,gaz)
%
%  [scenar,scenstr,Rres]=scenarzerodfci(geo,freq,pos,composition);
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
%  Cr�tion : 11 f�rier 2002
%  Auteur   : V. Basiuk
%
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)


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
% position des couches par rapport a l'axe magnetique, normalise au petit rayon
%
for k=1:3
  distH(k)   = abs((rH(k)-mean(geo.r0+geo.d0)))./mean(geo.a);
  distD(k)   = abs((rD(k)-mean(geo.r0+geo.d0)))./mean(geo.a);
  distHe3(k) = abs((rHe3(k)-mean(geo.r0+geo.d0)))./mean(geo.a);
end
%
% grand rayon max du plasma [distp] et grand rayon min [distm]
%
distp  = mean(geo.r0+geo.a);
distm  = mean(geo.r0-geo.a);
centre = mean(geo.r0+geo.d0);
%
% determination du scenario
%
Rres = NaN;
if strcmp(gaz,'H') & distH(1) > 0.5 & rH(1) < centre & B0 < 3

  scenar = 0;
  scenstr = 'FWEH    ';
  return
              
end
	     	
if strcmp(gaz,'H') & distH(1) < 0.7

  scenar = 1;
  scenstr = 'HMIN_H  ';
  Rres    = rH(1);
  return
              
end
if strcmp(gaz,'D') & distD(1) < 0.7

  scenar = 1;
  scenstr = 'HMIN_D  ';
  return
              
end
if strcmp(gaz,'He') & distD(1) < 0.7 

  scenar = 1;
  scenstr = 'HMIN_He ';
  return
              
end
if strcmp(gaz,'He3') & distHe3(1) < 0.7  

  scenar = 1;
  scenstr = 'HMIN_He3';
  return
              
end
if ~strcmp(gaz,'He3') & distH(1) > 0.7 & distH(2) > 0.7

  scenar  = 0;
  scenstr = 'FWEH    ';
  return
            
end
if strcmp(gaz,'H') & distH(1) > 0.7 & distH(2) < 0.7

  scenar  = 1;
  scenstr = 'HARM_2H ';
  return

end

	   




