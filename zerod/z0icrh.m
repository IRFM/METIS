% cette fonction calcul les effets du chauffage minoriatire a l'harmonique 1 de ICRH
% pour les test : [pth,pel,pion,esupra,einj,ecrit,teff0,taus,frloss] = z0icrh(zs,z0dinput.geo,z0dinput.cons,z0dinput.option);
function [pth,pel,pion,esupra,einj,ecrit,teff0,taus,frloss,rres,xres,frac,harm,nmino] = z0icrh(zs,geo,cons,option,profli)

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

% isotopic composition for option.gaz == 5
if option.gaz == 5
    nHe3onD = real(cons.iso);
    nTonD   = imag(cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(cons.iso));
    nTonD   = real(cons.iso);
end
cons.iso = real(cons.iso);

% choix du minoritaire
switch option.mino
    case 'He3'
        ag = 3;
        zg = 2;
        lg = 7.92e-3;
    case 'T'
        ag = 3;
        zg = 1;
        lg = 7.92e-3;
    case 'He4'
        ag = 4;
        zg = 2;
        lg = 4.55e-3;
    case 'D'
        ag = 2;
        zg = 1;
        lg = 6.46e-3;
    case 'B'
        ag = 11;
        zg = 5;
        lg = 4.576e-3 .* sqrt(11) / 5;
    otherwise
        ag = 1;
        zg = 1;
        lg = 4.576e-3;
end

% 1ere partie : la resonance
% recherche de la postion de la resonnance par rapport au centre plasma
rc     = zs.d0 + geo.R;
bres   =  2 .* pi .* option.freq .* 1e6 ./ (95.5e6 .* zg ./ ag);
% utilisation d'un forme de profil
if isfield(profli,'xli')
	x = profli.xli;
else
	x   = linspace(0,1,21);
end
ux  = (1  - x .^ 2);
ve  = ones(size(x));
vt  = ones(size(geo.R));
%qp  = min(zs.q0,zs.qmin) * ve + (zs.qa - min(zs.q0,zs.qmin)) * (x .^ 2);
if isfield(profli,'qjli')
	qp = profli.qjli;
else
	qp = z0qp(x,max(1,min(zs.q0,zs.qmin)),zs.qa);
end
if isfield(profli,'Raxe')
	rh  = profli.Raxe - geo.a * x;
	rl  = profli.Raxe + geo.a * x;
else
	rh  = geo.R * ve + zs.d0 * ux - geo.a * x;
	rl  = geo.R * ve + zs.d0 * ux + geo.a * x;
end
rr  = cat(2,rh(:,end:-1:2),rl);
xa  = cat(2,-geo.a *x(end:-1:2),geo.a *x);
qpp = max(0.5,cat(2,qp(:,end:-1:2),qp));
btot = (geo.b0 * ones(1,size(xa,2))) .* sqrt((xa ./ qpp ./ (geo.R *ones(1,size(xa,2)))) .^ 2 + ((geo.R *ones(1,size(xa,2))) ./ rr) .^ 2);

% choix de l'harmonique
switch option.mino
    case 'D'
        if (option.gaz == 2)  || (option.gaz == 3)
            % D dans D, recherche de l'harmonique superieure
            dd2 = abs(btot(:,21) - bres ./ 2);
            dd3 = abs(btot(:,21) - bres ./ 3);
            maskh = double(dd3 < dd2);
            bres  = maskh .* bres ./ 3 + (~maskh) .* bres ./ 2;
            harm  = 2 + maskh;
            zg    = zg .* harm;
        else
            maskh = (bres > max(btot,[],2));
            bres  = maskh .* bres ./ 2 + (~maskh) .* bres;
            harm  = round(mean(maskh + 1));
            zg    = zg .* harm;
        end
    case 'B'
            dd2 = abs(btot(:,end) - bres ./ 2);
            dd3 = abs(btot(:,end) - bres ./ 3);
            maskh = double(dd3 < dd2);
            bres  = maskh .* bres ./ 3 + (~maskh) .* bres ./ 2;
            harm  = 2 + maskh;
            zg    = zg .* harm;        
    otherwise
        maskh = (bres > max(btot,[],2));
        bres  = maskh .* bres ./ 2 + (~maskh) .* bres;
        harm  = round(mean(maskh + 1));
        zg    = zg .* harm;
end
dd  = (btot - bres * ones(1,size(btot,2))) .^ 2;
mask = (dd == (min(dd,[],2)*ones(1,size(xa,2))));
res   = sum(rr .*mask,2);
rres = res; % pour la sortie
xres  = min(0.95,abs(sum(xa .* mask,2)) ./ geo.a);
% position de la resonance pour le calcul de perte
aperte = sum(xa .* mask,2);
% facteur de securite sur la resonance
qperte = z0qp(xres,max(1,min(zs.q0,zs.qmin)),zs.qa);


% si les profils sont disponibles
if isfield(profli,'tep') && isfield(profli,'picrh')
	mask      = 0.1 .* exp(-(vt * profli.xli - xres * ve) .^ 2 ./ 0.003) + profli.picrh ./ max(eps,max(profli.picrh,[],2) *ve);
        %figure(24);clf
	%plot(profli.xli,mask)
    norme_res = trapz(profli.xli,max(eps,mask) .* profli.vpr,2);
	teres = trapz(profli.xli,max(eps,mask) .* profli.vpr .* profli.tep,2) ./ norme_res;
	neres = trapz(profli.xli,max(eps,mask) .* profli.vpr .* profli.nep,2) ./ norme_res;
	nheres = trapz(profli.xli,max(eps,mask) .* profli.vpr .* profli.nhep,2) ./ norme_res;
	tjres = trapz(profli.xli,max(eps,mask) .* profli.vpr .* profli.tip,2) ./ norme_res;
    switch option.gaz
        case 1
            njres  = neres .* max(1e13,zs.n1m - zs.nDm -zs.nTm) ./ zs.nem;
        case {2,5}
            njres  = neres .* zs.nDm ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
        case 3
            njres  = neres .* zs.nDm ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
        case 11
            njres  = neres .* max(1e13,zs.n1m - zs.nDm) ./ zs.nem;            
        otherwise
            njres  = neres .* zs.nhem ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
    end
    qperte = trapz(profli.xli,max(eps,mask) .* profli.vpr .* profli.qjli,2) ./ norme_res;

else
	% donnees associees (attention ce n'est pas les valeur local !)
	teres  = (zs.tem .* (1 + zs.ate) .* (1 - xres .^ 2) .^ zs.ate + zs.tem) ./ 2;
	neres  = ((zs.nem .* (1 + zs.ane) - zs.nebord) .* (1 - xres .^ 2) .^ zs.ane + zs.nebord + zs.nem) ./ 2;
    nheres = neres .* zs.nhem ./ zs.nem;
    tjres  = teres .* zs.tite;
    switch option.gaz
        case 1
            njres  = neres .* max(1e13,zs.n1m - zs.nDm -zs.nTm) ./ zs.nem;
        case {2,5}
            njres  = neres .* zs.nDm ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
        case 3
            njres  = neres .* zs.nDm ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
        case 11
            njres  = neres .* max(1e13,zs.n1m - zs.nDm) ./ zs.nem;            
        otherwise
            njres  = neres .* zs.nhem ./ zs.nem;   % attention on suppose un  plasma pure avec un minoritaire sans impuretees
    end
end	


% choix du majoritaire
switch option.gaz
    case 1
        nref = max(1e13,zs.n1m - zs.nDm -zs.nTm);
    case {2,3}
        nref = zs.nDm;
    case 5
        nref = zs.nDm;        
    case 11
        nref = zs.n1m - zs.nDm;
    otherwise
        nref = zs.nhem;
end
% densité minoritaire 
switch option.gaz
    case {1,2,3,4}
        switch option.mino
            case 'T'
                nmino  = neres .* nref ./ zs.nem .* max(1e-2,cons.iso .* option.cmin ./ harm .^ 2);
            case 'He'
                nmino  = nheres  .* option.cmin ./ harm .^ 2;
            otherwise
                nmino  = neres .* nref ./ zs.nem .* option.cmin ./ harm .^ 2;
        end
    case 5
        switch option.mino
            case 'He3'
                nmino  = nheres  .* option.cmin ./ harm .^ 2;
            case 'He'
                nmino  = neres .* nref ./ zs.nem .* option.frhe0 ./ harm .^ 2;
             otherwise
                nmino  = neres .* nref ./ zs.nem .* option.cmin ./ harm .^ 2;
        end
    case 11
       switch option.mino
            case 'He'
                nmino  = nheres  .* option.cmin ./ harm .^ 2;
            case 'B'
                nmino  = neres .* nref ./ zs.nem .* cons.iso  ./ harm .^ 2;
             otherwise
                nmino  = neres .* nref ./ zs.nem .* option.cmin ./ harm .^ 2;
        end
    otherwise
        nref = zs.nhem;
end

% plasma de fond
switch option.gaz
    case 1
        zj = 1;
        aj = 1;
    case 2
        zj = 1;
        aj = 2;
    case 3
        zj = 1;
        aj = (2  + 3 .* cons.iso)  ./ (1 + cons.iso);
    case 4
        zj = 2;
        aj = 4;
    case 5
        zj = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
        aj = mean((2  + 3 .* cons.iso)   ./  (1 + cons.iso));
    case 11
        zj = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        aj = mean((1  + 11 .* cons.iso)  ./  (1 + cons.iso));
end

% volume de plasma pour la resonance
% fact est determiner �l'aide de PION
% ne fonctionne pas tres bien
% choix du minoritaire
switch option.mino
    case {'T','B'}
        fact   = 1;
    otherwise
        fact   = 3.2;
end
if option.fact_mino > 0
    fact =  option.fact_mino;
end
  
%ebeta = phys.e .* 14.8 .*teres .^ (1/3) .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* njres ./ neres .* tjres) .^ (2/3);
%ecrit = phys.e .* 14.8 .*teres  .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* njres ./ neres) .^ (2/3);
% attention ecrit doit etre en eV 
ebeta = 14.8 .*teres .^ (1/3) .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* njres ./ neres .* tjres) .^ (2/3);
ecrit = 14.8 .*teres  .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* njres ./ neres) .^ (2/3);

egamma = 14.8 .* teres  .* (2.* sqrt(ag) .* zs.zeff) .^ (2/3);
tpara  = egamma ./ 8;
kpara  = option.nphi ./ geo.R;
dr     = geo.R .* kpara .* sqrt(2 .* tpara .* phys.e ./ phys.ua ./ ag) ./ 2 ./ pi ./ option.freq ./ 1e6;
frac   = min(1,max(0.05,fact .* 2 .* geo.a .*geo.K .* dr ./ zs.sp));


% en attendant la mise au point (depot dans un rayon de 40 cm pour TS et Jet)
%frac   = 0.25;
vmino  = zs.vp .* frac;   % fraction  effective du volume de minoritaire accelere par ICRH

% puissance par unite de volume delivree aux ions minoritaire
pm = cons.picrh .* (1 - zs.frloss_icrh) ./ vmino;

if all(pm <= 1e3) 
   % gain de temps
   pth = pm;
   pel = pm;
   pion = pm;
   esupra = pm;
   einj  = 2 .* zs.te0;
   teff0 = zs.te0 .* zs.tite;
   taus =  6.27e8 .* ag ./ zg .^ 2 ./ 17 ./ (neres ./ 1e6).* teres .^ (3/2);
   frloss = 0 .* zs.te0;
   harm = harm .* vt;
   return
end

% 2eme partie : la fonction de distribution (j == f)
lnldei  = 15.2 - 0.5 .* log(neres./1e20) + log(teres ./1e3);
taus  =  6.27e8 .* ag ./ zg .^ 2 ./ lnldei ./ (neres ./ 1e6).* teres .^ (3/2);
taus  = taus .* (pm >= 1e3) + 1e-6 .* (pm <1e3); 

%ecgs  = phys.e ./ sqrt(4 .* pi .* phys.epsi0);
%zetab  = (phys.ua .* ag) .* pm  ./ ( 8 .* sqrt(pi) .* neres .* nmino .* zg .^ 2 .* ecgs .^ 4 .* 17) .* sqrt(2 .* phys.e .* teres ./ phys.me);
zeta  = pm .* taus ./ 3 ./ nmino ./ phys.e ./ teres;
% vecteur des vitesses
vmax   = phys.c;
%v = logspace(0,log(vmax) / log(10),1001);
% optimisation temps de calcul a verifier
v = logspace(0,log(vmax) / log(10),101);
% suite
ep    = 2/3/sqrt(pi);
lj    = sqrt(aj .* phys.ua ./ 2 ./ phys.e ./ tjres);
le    = sqrt(phys.me ./ 2 ./ phys.e ./ teres);
rj    = njres .* zj  .^ 2 .* lj ./ neres ./ le;
ev    = ag .* phys.ua .* v .^ 2 ./2;
ej    = ag ./ aj .* phys.e .* tjres  .* ((2 + 2 .* rj + 3 .* zeta) ./ 2 ./ ep ./ (2 + 3 .* zeta)) .^ (2/3);
vv    = ones(size(ev));
vt    = ones(size(teres));
inter = rj .* ((2 .* ag + aj) .* (2 + 3 .* zeta) .* teres -  4 .* ag .* tjres)./  ...
        (2 .* ag  .* tjres .* (2 +  2 .* rj + 3 .* zeta));
lnfv  = - 2 .* (vt * ev) ./ (phys.e .* (teres * vv)) ./ (2  + 3 .* (zeta *vv)) .* (1 + (inter * vv) .* hh((vt * ev) ./ (ej*vv)));
lnfv  = lnfv .* ((pm >= 1e3) *vv);

% temperature effective
%teff  = (teres * vv).* (1+ zeta * vv) ./ (1 + ((rj .* (teres - tjres + zeta .* teres) ./ tjres ./ ( 1 + rj + zeta)) * vv) ./ (1 + ( (vt * ev)  ./ (ej * vv)) .^ (3/2)));
%teff0  = teres.* (1+ zeta) ./ (1 + inter);

% normalisation de la fonction de distribution
fv    = exp(lnfv);
fvn   = (nmino * vv) .* fv./ (2 .* pi .* trapz(v,(vt * v) .* fv,2) * vv) .* ((pm >= 1e3) *vv);


% matrice des pentes
pentes = ((log(max(eps,fvn(:,1))) * vv) - log(max(eps,fvn))) ./ max(eps,vt * (ev-ev(1)));
pentes = max(pentes(:,3:end),[],2);
% thermique a soustraire
warning off
fth = min(fvn,exp(log(fvn(:,1)) * vv  - pentes * (ev-ev(1))));
warning on
fth(~isfinite(fth)) = 0;
teff0  =  pi .* ag .* phys.ua .* trapz(v,(vt * v) .^ 3 .* fth,2) ./ nmino  ./ phys.e;


% thermique a soustraire
%fth   = (((ag .* phys.ua ./ 2 ./ pi ./ teff0 ./ phys.e) .^ (3/2)) * vv) .* exp(-ag .* phys.ua .* (vt * v).^ 2 ./(teff0 * vv) ./ 2 ./ phys.e);  
%fth   = (((ag .* phys.ua ./ 2 ./ pi ./ (zs.tite .* teres) ./ phys.e) .^ (3/2)) * vv) .* exp(-ag .* phys.ua .* (vt * v).^ 2 ./((zs.tite .* teres) * vv) ./ 2 ./ phys.e);  
% normalisation a v = 0 
%fth   = fth ./ (2 .* pi .* trapz(v,(vt * v)  .* fth,2) * vv);
%fth   = fth .* ((fvn(:,1) ./ fth(:,1)) * vv); 


% fonction de distribution des  vitesses perpendiculaires (deja integree sur v//)
% calcul de l'energie suprathermique contenu dans le plasma :    (fvn -fth) * (1/2*m*v^2) * (2*pi*v)dv
% avec ou sans effet du ripple
esupra_st = max(0,pi .* ag .* phys.ua .* trapz(v,(vt * v) .^ 3 .* (fvn - fth),2) .* vmino);
etot_st = max(0,pi .* ag .* phys.ua .* trapz(v,(vt * v) .^ 3 .* fvn,2) .* vmino);

% calcul de l'energie moyenne des ions
einj  = max(tjres,pi .* ag .* phys.ua .* trapz(v,(vt * v) .^ 3 .* fvn,2) ./ nmino  ./ phys.e);

% puissance sur les electrons (staydy state, L-G Eriksson)
%pel_st   = 2 .* esupra_st ./ max(eps,taus) .* (pm >= 1e3); 
%pel_st   = max(0,2 .* etot_st ./ max(eps,taus) .* (pm >= 1e3)); 
pel_st   = max(0,2 .* esupra_st ./ max(eps,taus) .* (pm >= 1e3)); 

% cas du ripple
if option.rip == 1
   prip   = cons.picrh  .* (1 - zs.frloss_icrh) - zs.picrh;
   pel_st = max(pel_st - prip,0.1 .* pel_st);
   pion   = max(0,cons.picrh .* (1 - zs.frloss_icrh) - pel_st - prip);
else
   % conservation du total
   pion     = max(0,cons.picrh  .* (1 - zs.frloss_icrh) - pel_st);       % le chauffage des ions est instantanne
end

if option.transitoire == 1
   % temps de variation sur les electron
   % calcul de taueff pour
   tauseff = taus ./ 2;
   tauseff(tauseff <= 0) = zs.taue(tauseff <=0)./1e3;
   tauseff(pm<1e3) = mean(tauseff(pm >=1e3));
   tauseff(1)      = mean(tauseff(pm >=1e3));
   % calcul de la puissance thermique  
   if option.evolution == 1
   	esupra_ini = zs.esup_icrh(1);	        
   else
   	esupra_ini = 0.01 .* pel_st(1) .* mean(tauseff);
   end
   [te,esupra] = z0ode(cons.temps,pel_st,tauseff,esupra_ini);
   esupra = max(0,esupra);
   pel   = max(0,esupra ./ tauseff);
   pth   = max(0,pion + pel);
else
   % temps de variation sur les electron
   % calcul de taueff pour
   tauseff = taus ./ 2;
   tauseff(tauseff == 0) = zs.taue(tauseff ==0)./1e3;
   % calcul de la puissance thermique  
   esupra=esupra_st;
   pel   = pel_st;
   pth   = max(0,pion + pel);
end

%  % Calcul alternatif
%  [esupra_alt,palt_th]   = zsupra0(cons.temps,cons.picrh,taus,0.1,ecrit,einj,ag);
%  figure(21)
%  clf
%  subplot(2,2,1)
%  plot(cons.temps,frac)
%  
%  subplot(2,2,2)
%  plot(cons.temps,esupra,'b',cons.temps,esupra_st,'r',cons.temps,esupra_alt,'k.');
%  
%  subplot(2,2,3)
%  loglog(0.5 .* ag .* phys.ua .* v .^ 2 ./ phys.e,fvn)
%  axis([1,2e6,1,Inf])
%  
%  subplot(2,2,4)
%  plot(cons.temps,zfract0(ecrit,einj),'r',cons.temps,pion ./ max(eps,pth),'b');
%  %semilogy(cons.temps,einj,cons.temps,teff0)
%  drawnow
% securite simulation courte
pel  = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pel);
pth  = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pth);
pion = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pion);


% estimation de la fraction perdue lors de la 1ere orbite
rloc   = (geo.R  + aperte) * vv;
ral    = lg .* sqrt((vt * ev) ./ 1e3 ./ phys.e) ./ (((geo.b0 .* geo.R)*vv) ./ rloc); % en m
%dp_     = (geo.R*vv) .* (2 .* (qperte * vv) .* ral ./ (geo.R *vv)) .^ 2/3;
% si la resonnance est loin de l'axe magnetique
%  % patato
dp1     = (rres * vv) .* (2 .* (qperte * vv) .* ral ./ (rres *vv)) .^ 2/3;
%  % banana
dp2     = sqrt(abs(aperte * vv) ./ (rres * vv)) .* ral .* (qperte * vv);
dp      = dp2 .* (dp2 < abs(aperte * vv)) + dp1 .* (dp2 >= abs(aperte * vv));
posext = (ral + dp + (zs.d0 .* (1 - xres .^ 2)) * vv  + aperte * vv  - geo.a  * vv);
%mask   = (ral + dp + (zs.d0 .* (1 - xres .^ 2)) * vv  + aperte * vv  - geo.a  * vv) > 0;
% le depot a une certaine largeur
poids  = min(1,max(0,(posext + (dr*vv))  ./ max(2 .* dr*vv,1e-3)));
eloss  = 4 .* pi .* trapz(v,(vt * v) .^ 2 .* fvn .* poids ,2);
eref   = 4 .* pi .* trapz(v,(vt * v) .^ 2 .* fvn ,2);
frloss = max(0,min(1, eloss ./ (eref + eps)));
% pour eviter de propager des erreur numerique dans les tests de non regression
frloss(frloss < (3 .* eps)) = 0;

% mise a dimension
harm = harm .* vt;

%  disp('in z0icrh')
%  keyboard

% fonction H p 750 NF 1975 (15), Stix.
function s= hh(x)

warning off
s = (2./x) .* ((1/6) .* log((1 - sqrt(x) + x)./(1 + 2 .* sqrt(x) + x)) +  ...
      sqrt(1/3).* ((pi ./ 6) + atan((2 .* sqrt(x) - 1) ./ sqrt(3)))); 
s(x==0) = 1;
% control
%ss= cumtrapz(x,1./(1+x .^(3/2)))./x

warning on


% fonction G p 749 NF 1975 (15), Stix.
function s=gg(x)

ep =2/3/sqrt(pi);
s = ep .* x ./ (1+ 2 .* ep .* x .^ 3);


% fonction PHI p 749 NF 1975 (15), Stix.
function s=phi(x)

ep =2/3/sqrt(pi);
s = ep .* (3 .* x + 2 .* x .^ 3) ./ (1+ 2 .* ep .* x .^ 3);

