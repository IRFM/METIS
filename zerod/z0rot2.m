% dans metis la dirction toroidal est dans le sens du courant plasma, 
% de telle sorte que Btheta est toujours positif
% le moment injecte par l'IDN est posistif si l'injection est co courant
% le champs toroidal est positif s'il est dans le sens du courant
% calcul de la vitesse de rotation moyenne
function [wrad,snbi,sicrh,sfus,sripth,sriplh,sripicrh,sturb,fact,profli] = z0rot2(zs,option,cons,geo,profli)


error('This function is obsolette ...')

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
%rnbi      = sin(option.angle_nbi./180*pi);
rnbi      = zs.mu0_nbi;
snbi      = geo.R .* (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi) .* phys.ua .* (zs.pnbi_th ./ option.einj./phys.e)  .* sqrt(2 .*option.einj .* phys.e ./ phys.ua ./ (2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi)) .* rnbi;

% icrh (constante a determiner)
% largeur orbite wesson p 132 (potato)
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

% facteur libre (calibre sur le calcul du courant genere par les alphas de spot)
% valide sur le chocs JET tel 53521 ...
gg        = 0.15;
% la couche est supposee centrale ou cote LFS  => rotation co-courant 
% de  NF 42 , p  959  (2002)  L.-G. Eriksson and F. Porcelli
%sicrh     = gg .* sqrt(ag.*phys.ua).* pion ./ geo.a .^ 2 ./ sqrt(max(1e3,zs.einj_icrh).*phys.e) .* sqrt(d./geo.R) .* d .^ 2 ./ geo.R ./ geo.a ./ option.fmin;

% recherche de la postion de la resonnance par rapport au centre plasma 
rc     = zs.d0 + geo.R;
bres   =  2 .* pi .* option.freq .* 1e6 ./ (95.5e6 .*zg ./ ag);
% utilisation d'un forme de profil
x   = profli.xli;
ux  = (1  - x .^ 2);
ve  = ones(size(x));
vt  = ones(size(geo.R));
qp  = profli.qjli;
rh  = profli.Raxe - geo.a * x;
rl  = profli.Raxe + geo.a * x;
rr  = cat(2,rh(:,end:-1:2),rl);
xa  = cat(2,-geo.a *x(end:-1:2),geo.a *x);
qpp = max(0.5,cat(2,qp(:,end:-1:2),qp));
btot = (geo.b0 * ones(1,size(xa,2))) .* sqrt((xa ./ qpp ./ (geo.R *ones(1,size(xa,2)))) .^ 2 + ((geo.R *ones(1,size(xa,2))) ./ rr) .^ 2);

% choix de l'harmonique
maskh = (bres > max(btot,[],2));
bres  = maskh .* bres ./ 2 + (~maskh) .* bres;
harm  = round(mean(maskh + 1));
zg    = zg .* harm;

% position
dd  = (btot - bres * ones(1,size(xa,2))) .^ 2;
mask = (dd == (min(dd,[],2)*ones(1,size(xa,2))));
res   = sum(rr .*mask,2);
sicrh     = gg .* sqrt(ag.*phys.ua).* pion ./ geo.a .^ 2 ./  ...
            sqrt(max(1e3,zs.einj_icrh).*phys.e) .* sqrt(dp./geo.R);

% selon la position
pos   = (res + ral + ds - rc);
test  = tanh(pos./zs.d0);
sicrh = test .* sicrh;

% terme de fusion
% largeur orbite
%ral    = 4.55e-3.* sqrt(3.56e6 ./ 1e3) ./ (geo.b0 ./ geo.R); % en m
%d      = geo.R .* (2 .* (min(30,zs.qa) + 1) ./2 .* ral ./ geo.R ) .^ 2/3;
% facteur libre (calibre sur le calcul du courant genere par les alphas de spot)
gfus        = 0.15;
qeff = 1./(1./zs.qa + 1./zs.qmin);
dp      = (qeff ./ geo.b0) .^ (2/3) .* (geo.R .* 3.56e6).^ (1/3);
sfus   = gfus .* sqrt(4.*phys.ua) .* zs.pfus .* zfract0(zs.ecrit_he,3.56e6) ./ geo.a .^ 2 ./ sqrt(3.56e6.*phys.e) .* sqrt(dp./geo.R);

% facteur pour passer de omega au moment total omega  * fact = M
%arot = min(1,max(0,(1 - zs.tite))) .* zs.ane + min(1,max(0,zs.tite)) .* zs.ate; % selon le chauffage des ions
%mass = (zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm  + 3 .* zs.nTm + 4 .* zs.nhem +  2 .* option.zmax .* zs.nimpm;
%fact = phys.ua .* mass .* zs.vp .* geo.R .* (1 + zs.ane) .* (1 + arot) ./ ( 1 + zs.ane + arot);
% calcul plus precis
mass = ((((zs.n1m - zs.nDm - zs.nTm) + 2 .* zs.nDm  + 3 .* zs.nTm) ./ zs.n1m) * ve)  .* profli.n1p + ...
       4 .* profli.nhep  + ((option.zimp+1) * 2 + option.rimp .* (option.zmax+1) .* 2) .* profli.nzp;
fact = phys.ua .* trapz(x,profli.vpr .* mass .* profli.tip .* profli.Raxe,2) ./  ...
       trapz(x,profli.vpr .* profli.tip,2) .* trapz(x,profli.vpr,2) ;

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
sripth    =  - frot   .* zs.tem .* zs.tite ./ geo.a ./ z1 ./ ...
              ( phys.mu0 .* max(1e3,zs.ip ./ zs.peri)) .* fact ./ max(1e-6,zs.taue) .* 3 ;
calrip   = sripth  ./ (zs.ip .^ 2 .* zs.RR);

if option.rip == 1
   % on suppose la proportionnalite a la puissance (il faut sans doute un facteur devant)
   priplh   = max(0,cons.plh - zs.plh_th);
   pripicrh = max(0,cons.picrh - zs.picrh_th);
   % le terme source est suppose poportionnele au  courant perdu
   % seul le courant ionique compte (courant en retour pour LH)
   sriplh   =    0.3 .*   (zs.plh_th > 1e4)  .* priplh .* calrip;   % a 3 vth 
   sripicrh =    - (pion > 1e4) .* sqrt(zs.tem ./ max(1e3,zs.einj_icrh))  .* pripicrh .* calrip; % Le signe n'est pas connu
   sripth    =   zs.pohm .* calrip;
else
   sriplh    =  0 .* sripth;
   sripicrh  =  0 .* sripth;
   sripth    =  zs.pohm .* calrip;  % effet  minime en ripple faible 
end

% il faut encore regle les gain du ripple icrh et lh et le seuil 

% source de rotation spontannee (expose ttf 2006 Marseille, Rice , ...)
% utilisation de wion uniquement pour TS et TCV
%wion   = (3/2) .* trapz(x,profli.nip .* profli.tip .* phys.e  .* profli.vpr,2);
%betan  = zs.wth .* (1.6.*pi./3) .* geo.a ./ zs.vp ./ geo.b0 ./ zs.ip;
%betan  = 2 .* wion .* (1.6.*pi./3) .* geo.a ./ zs.vp ./ geo.b0 ./ zs.ip;
%ca     = sqrt(geo.b0 .^ 2 ./ phys.mu0 ./ max(1e13,zs.nem) ./ zs.meff ./ phys.mp);
%vturb2  = betan  .* ca ;
% nouveau scaling Rice NF  47, (2007) p 1618-1624
% utilisation de Pion uniquement pour TS et TCV
press  = 2 .* trapz(x,profli.nip .* profli.tip .* phys.e  .* profli.vpr,2) ./ zs.vp;
vturb  = 1e11 .* geo.b0 .^ 1.1 .* press .* zs.ip .^ (-1.9) .* geo.R .^2.2;
%figure(2100);clf;plot(cons.temps,vturb,'r',cons.temps,vturb2,'b')
sturb  = fact  ./ max(1e-6,zs.taue .* option.taurotmul) .* vturb;

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

%  figure(151);clf
%  plot(cons.temps,sturb,'b',cons.temps,snbi,'r',cons.temps,sicrh,'m',cons.temps,sfus,'c', ...
%       cons.temps,sripth,'g',cons.temps,sriplh,'c:',cons.temps,sripicrh,'m:')
%  drawnow


% equa diff
srot           = (sturb + snbi + sicrh + sfus + sripth + sriplh + sripicrh) ./ fact;
srot(~isfinite(srot)) = 0;

if option.evolution == 1
	vini = zs.wrad(1) .* geo.R(1);
else
	vini           = srot(1) .* zs.taue(1) .* option.taurotmul .* geo.R(1);
	if ~isfinite(vini)
   		vini = 0;
	end
end

%[tvtorm,vtorm] = ode23('ztf0',cons.temps,vini,[],cons.temps,srot,taueff);
% le facteur 3 vient des donnees de Jet (Mantica TTF Marseille 2006)
% l'effet est visible sur TS
% il faut sans doute utiliser la diffusivite et la convection de la matiere
[tvtorm,vtorm] = z0ode(cons.temps,srot,zs.taue.* option.taurotmul,vini);
wrad           = vtorm ./ geo.R;

% calcul du profil de rotation toroidal
% palsma de fond
switch option.gaz
case 1
   zj = 1;
   aj = 1;
case 2
   zj = 1;
   aj = 2;
case 3
   zj = 1;
   aj = mean((2  + 3 .* cons.iso)  ./  (1+ cons.iso));
case 4
   zj = 2;
   aj = 4;
end

% impurete principale
zimp = option.zimp;
aimp = ceil(zimp .* (7/3));

% 2ieme impurete
zmax = option.zmax;
amax = ceil(zmax .* (7/3));


% pour chaque espece d'ions
nDp   = max(1e13,profli.n1p .* ((zs.nDm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nTp   = max(1e13,profli.n1p .* ((zs.nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nHp   = max(1e13,profli.n1p - nTp - nDp);
nhep  = max(1e13,profli.nhep);
nz1p  = max(1e13,profli.nzp);
nz2p  = max(1e11,profli.nzp .* option.rimp);   
% masse
Mtor    = phys.mp .*  (nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p);
% omega homothetic a Ti (cas Jet avec NBI)
omega  = profli.tip;
if option.qdds > 0
	mask1  = (profli.qjli <= 1);
	omega1 = max((~mask1) .* omega,[],2) * ones(size(x));
	omega  = (~mask1) .* omega + mask1 .* omega1;
end

% renormalisation de la quantite de mouvement
normv = trapz(x,profli.vpr .* Mtor .* omega .* profli.Raxe,2);
normf = trapz(x,profli.vpr .* Mtor .* profli.Raxe,2);
profli.omega  = omega .* (((wrad .* normf) ./ max(eps,normv)) * ve);
% calcul de la rottaion poloidal 
warning off


% formulaire ORNL
lnii   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(1); 
lnhe   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2); 
lnz1   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zimp); 
lnz2   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zmax); 

% pour l'espece principale
% Tokamaks, Wesson p 663
taui_h     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 1)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nHp ./ lnii ./ 1 .^ 4 ;
taui_d     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 2)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nDp ./ lnii ./ 1 .^ 4 ;
taui_t     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nTp ./ lnii ./ 1 .^ 4 ;
taui_he     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 4)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe ./ 2 .^ 4 ;
taui_z1     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* aimp)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz1p ./ lnz1 ./ zimp .^ 4 ;
taui_z2     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* amax)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz2p ./ lnz2 ./ zmax .^ 4 ;

% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 1 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_h;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkh     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 2 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi./ taui_d;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkd     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_t;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkt     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 4 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkhe     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ aimp ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz1     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ amax ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz2     = max(-2.1,min(1.7,fk));

warning on



% changement de repere
dpsidrho = abs(pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve));
dpsidrho(:,1) = NaN;

% vitesse poloidal 
b2 = (geo.b0 .^ 2) * ve; 
% standard  + orbit squeezing
gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .*nHp,0,2,2,1) ./ nHp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_h  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_h(:,1) = 2 .* utheta_h(:,2) - utheta_h(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nDp,0,2,2,1) ./ nDp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_d  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_d(:,1) = 2 .* utheta_d(:,2) - utheta_d(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nTp,0,2,2,1) ./ nTp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_t  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_t(:,1) = 2 .* utheta_t(:,2) - utheta_t(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep,0,2,2,1) ./ nhep;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_he  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho; 
utheta_he(:,1) = 2 .* utheta_he(:,2) - utheta_he(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nz1p,0,2,2,1) ./ nz1p;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_z1  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zimp .* gtheta ./ dpsidrho; 
utheta_z1(:,1) = 2 .* utheta_z1(:,2) - utheta_z1(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nz2p,0,2,2,1) ./ nz2p;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_z2  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zmax .* gtheta ./ dpsidrho; 
utheta_z2(:,1) = 2 .* utheta_z2(:,2) - utheta_z2(:,3);

% calcul du moment poloidal injecte par IdN
% dans le sens de la vistesse diamagnetique electronique (couplage par les collisions sur les electrons)
pfast      = profli.pnbi .* ((zs.esup_nbi ./ max(1,trapz(x,profli.vpr .* profli.pnbi,2))) * ve);
nfast      = max(1,pfast .* 2 ./ (phys.e .* option.einj)); 
btot       = sqrt(profli.fdia .^2 .* profli.r2i + profli.bpol .^ 2);
utheta_nbi = profli.jnbicd ./ nfast ./ btot ./ phys.e;

% cacul du moment total
Rtor    = Mtor    .* profli.omega ./ profli.r2i;
dpsidrho(:,1) = 0;
% calcul du champ electrique radial (Er gradient(rho))
Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
          2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
          3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
          2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
          aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
          amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1);
Ptor    = Ptor ./ (profli.rmx(:,end) * ve);
			
			
Unbi   = phys.mp .*  ((2 .* (1-cons.ftnbi) + 3 .* cons.ftnbi) * ve) .* nfast .* utheta_nbi;
Utor    = phys.mp .*  (      nHp .* utheta_h + 2 .* nDp  .* utheta_d  +  ...
                        3 .* nTp .* utheta_t + 4 .* nhep .* utheta_he +  ...
			aimp .* nz1p .* utheta_z1  + amax .* nz2p .* utheta_z2);
						
Utor    = option.signe .* profli.fdia .* profli.r2i .* dpsidrho .* (Utor + Unbi);			

			
profli.er  = (Ptor + Rtor .* profli.r2i .* dpsidrho  - Utor)  ./ Mtor;

% les champ en Rmax
btor         = option.signe .* (profli.fdia ./ (profli.Raxe + geo.a * x));
grho         = abs((profli.rmx(:,end) * ve) ./ max(eps,pdederive(x,profli.Raxe + geo.a * x,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profli.psi,0,2,2,1)./ (profli.Raxe + geo.a * x) .* grho ./ (profli.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);

% calul de la vitessse toroidal en Rmax pour l'impurete principale
profli.utheta  = utheta_z1 + Unbi ./ Mtor;
profli.vtheta  = profli.utheta .* bpol;
profli.vtor    = (profli.er .* grho - grho ./ phys.e ./ zimp ./ nz1p .* ...
                 pdederive(x,phys.e .* profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) + ...
		 profli.vtheta .* btor) ./ max(eps,bpol);
profli.vtor(:,1) = 2 .*  profli.vtor(:,2) - profli.vtor(:,3);
		 