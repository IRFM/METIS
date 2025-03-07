% calcul de la vitesse de rotation moyenne
%  pour test : wrad = z0rot(zs,z0dinput.option,z0dinput.cons,z0dinput.geo,1);
function [wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,fact] = z0rot(zs,option,cons,geo,test)

% constante physique (phys)
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


% calcul des sources
% idn
rnbi      = sin(option.angle_nbi./180*pi);
snbi      = geo.R .* (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi) .* phys.ua .* (zs.pnbi_th ./ option.einj./phys.e)  .* sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi)) .* rnbi;

% icrh (constante a determiner)
% largeur orbite wesson p 132 (patato)
ral    = 7e-3.* sqrt(zs.einj_icrh ./ 1e3) ./ (geo.b0 ./ geo.R); % en m
ds     = geo.R .* (2 .* ral ./ geo.R ) .^ 2/3 .* (zs.qa .^ (2/3) +1);
% facteur de largeur d'orbite (formule de  NF 42 , p  959  (2002)  L.-G. Eriksson and F. Porcelli su la quantite de mouvement)
dp      = (zs.q0 ./ geo.b0) .^ (2/3) .* (geo.R .* zs.einj_icrh).^ (1/3);
% puissance au ions
pion  = zs.pion_icrh;
pion  = (pion > 1e5) .* pion;

switch option.mino
case 'He3'
   ag = 3;
   zg = 2;
case 'T'
   ag = 3;
   zg = 1;
case 'He4'
   ag = 4;   
   zg = 2;
case 'D'
   ag = 2;   
   zg = 1;
otherwise
   ag = 1;   
   zg = 1;
end
% facteur libre
gg        = 0.015;
% la couche est supposee centrale ou cote LFS  => rotation co-courant 
% de  NF 42 , p  959  (2002)  L.-G. Eriksson and F. Porcelli
%sicrh     = gg .* sqrt(ag.*phys.ua).* pion ./ geo.a .^ 2 ./ sqrt(max(1e3,zs.einj_icrh).*phys.e) .* sqrt(d./geo.R) .* d .^ 2 ./ geo.R ./ geo.a ./ option.fmin;

% recherche de la postion de la resonnance par rapport au centre plasma 
rc     = zs.d0 + geo.R;
bres   =  2 .* pi .* option.freq .* 1e6 ./ (95.5e6 .*zg ./ ag);
% utilisation d'un forme de profil
x   = linspace(0,1,21);
ux  = (1  - x .^ 2);
ve  = ones(size(x));
vt  = ones(size(geo.R));
%qp  = min(zs.q0,zs.qmin) * ve + (zs.qa - min(zs.q0,zs.qmin)) * (x .^ 2);
qp  = z0qp(x,min(zs.q0,zs.qmin),zs.qa);
rh  = geo.R * ve + zs.d0 * ux - geo.a * x;
rl  = geo.R * ve + zs.d0 * ux + geo.a * x;
rr  = cat(2,rh(:,end:-1:2),rl);
xa  = cat(2,-geo.a *x(end:-1:2),geo.a *x);
qpp = cat(2,qp(:,end:-1:2),qp);
btot = (geo.b0 * ones(1,size(xa,2))) .* sqrt((xa ./ qpp ./ (geo.R *ones(1,size(xa,2)))) .^ 2 + ((geo.R *ones(1,size(xa,2))) ./ rr) .^ 2);
dd  = (btot - bres) .^ 2;
mask = (dd == (min(dd,[],2)*ones(1,size(xa,2))));
res   = sum(rr .*mask,2);
sicrh     = gg .* sqrt(ag.*phys.ua).* pion ./ geo.a .^ 2 ./ sqrt(max(1e3,zs.einj_icrh).*phys.e) .* sqrt(dp./geo.R);

% selon la position
pos   = (res + ral + ds - rc);
test  = tanh(pos./zs.d0);
sicrh = test .* sicrh;

% terme de fusion
% largeur orbite
%ral    = 4.55e-3.* sqrt(3.56e6 ./ 1e3) ./ (geo.b0 ./ geo.R); % en m
%d      = geo.R .* (2 .* (min(30,zs.qa) + 1) ./2 .* ral ./ geo.R ) .^ 2/3;
qeff = 1./(1./zs.qa + 1./zs.qmin);
dp      = (qeff ./ geo.b0) .^ (2/3) .* (geo.R .* 3.56e6).^ (1/3);
sfus   = gg .* sqrt(4.*phys.ua) .* zs.pfus .* zfract0(zs.ecrit_he,3.56e6) ./ geo.a .^ 2 ./ sqrt(3.56e6.*phys.e) .* sqrt(dp./geo.R);

% facteur pour passer de vtor au moment total vtor *fact = M
arot = min(1,max(0,(1 - zs.tite))) .* zs.ane + min(1,max(0,zs.tite)) .* zs.ate; % selon le chauffage des ions
mass = (zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm  + 3 .* zs.nTm + 4 .* zs.nhem +  2 .* option.zmax .* zs.nimpm;
fact = phys.ua .* mass .* zs.vp .* geo.R .* (1 + zs.ane) .* (1 + arot) ./ ( 1 + zs.ane + arot);
% equation resolue : dV/dt = -V/tau + S;


% effet du ripple :
warning off
z1   =  (1+(zs.meff >3));
qeff = 1./(1./zs.qa + 1./zs.qmin);
% Tokamaks, Wesson p 661
ras    = geo.a ./ geo.R;
lnii   = 17.3 - 0.5.*log(zs.nem ./ 1e20) + 3/2 .* log(zs.tem .* zs.tite ./ 1e3); % pour H, D et T uniquement
% pour l'espece principale
% Tokamaks, Wesson p 663
taui   =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* zs.meff)  .* ...
	                   (phys.e .* zs.tem .* zs.tite) .^ (3/2) ./ zs.nim ./ lnii ./ z1 .^ 4 ;
% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
nui    = 1./ min(taui,1);
vthi   = sqrt(2 .* phys.e .*zs.tem .* zs.tite ./ zs.meff ./ phys.mp);
nuis   = nui .* geo.R .* qeff .* ras .^ (-3/2) ./ vthi;
fk     = (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* ras .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* ras .^ 3);
fk     = max(-2.1,min(1.7,fk));
frot   = (fk + (7/2) - 1) ;
warning on

% terme ohmique du au champ electrique // (C. Bourdelle these  5-15)
% EXPRESSION THESE C. BOURDELLE 5-19 (il y a toujours un peu de ripple dans un tokamak)
sripth    =  - frot   .* zs.tem .* zs.tite ./ geo.a ./ z1 ./ ( phys.mu0 .* max(1e3,zs.ip ./ zs.peri)) .* fact ./ zs.taue ;

if option.rip == 1
   priplh   = max(0,cons.plh - zs.plh_th);
   pripicrh = max(0,cons.picrh - zs.picrh_th);
   sriplh   =   sripth .* priplh .*0.000015;  % le facteur est empirique 
   sripicrh =   test .* sripth .* pripicrh .* 0.000015; % Le signe n'est pas connu

else
   sripth    =  1e-5 .* sripth; % cas faible ripple ~ 11/1.9
   sriplh    =  0 .* sripth;
   sripicrh  =  0 .* sripth;
end



% terme de friction modifiant le terme en taue
% c'est une correction sur le volume du plasma en rotation
% le calcul est complique, on neglige pour le moment
% le facteur devant l'expression est a adapter
%cxv = 7e-14;
%disp('il faut mettre le vrai volume concerne par le neutres') 
%friction = phys.ua .* zs.meff .* (zs.nim ./ zs.nem .* zs.nebord) .^ 2 .* cxv .* geo.R .* zs.vp;
%friction = (zs.tite .* zs.tebord ./ zs.tem) .*phys.ua .* zs.meff .* (zs.nim ./ zs.nem .* zs.nebord) .^ 2 .* cxv .* geo.R .* zs.vp;
% (zs.tite .* zs.tebord ./ zs.tem) = rapport rotation moyenne a la rotation de bord vue pour la friction
%taueff         = (friction + 1./ max(1e-6,zs.taue .* (1 + zs.modeh))).^ -1;


% equa diff
srot           = (snbi + sicrh + sfus + sripth+sriplh+sripicrh) ./ fact;
vini           = zs.wrad(1).* geo.R(1);
if ~isfinite(vini)
   vini = 0;
end
%[tvtorm,vtorm] = ode23('ztf0',cons.temps,vini,[],cons.temps,srot,taueff);
[tvtorm,vtorm] = z0ode(cons.temps,srot,zs.taue,vini);
wrad           = vtorm ./ geo.R;

if nargin > 4
   t = cons.temps;
   wexp = evalin('base','post.z0dinput.exp0d.wrad');
   vref   =(2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi) .* phys.ua .* (zs.pnbi_th ./ option.einj./phys.e)  .* ...
            sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi)) .* rnbi;
   vref   = vref .* zs.taue ./ ( zs.meff .* phys.ua.* (zs.nem./ (1+(zs.meff == 4))) .* zs.vp);
   figure(51);clf;plot(t,vref,'b',t,wrad.*geo.R,'r');
   figure(52);clf;plot(t,wexp,'b',t,wrad .* (1+arot),'r');
   keyboard
end
