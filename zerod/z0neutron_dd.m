% calcul le flux de neutron dd
% le calcul a ete valide pour TS , JET sur le choc 38437 et pour les prevision CDR de JT60-SA
function [neutron_total,neutron_th,neutron_nbi_th,neutron_nbi_nbi,pddfus,proton_dd, picrh_nbi,einj] = ...
    z0neutron_dd(option,cons,zs,profli,pnorme)

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

% case of Hydrogen NBI in DT plasma
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    gas_nbi = -1;
else
    gas_nbi = option.gaz;
end

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



% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - option.zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* option.zimp;
end
dd   = abs(Z_el - option.zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* option.zmax;
end


x = profli.xli;
temps = cons.temps;
ux  = 1 - x .^ 2;
ve  = ones(size(x));
vt  = ones(size(temps));
spr = profli.spr;
vpr = profli.vpr;
nep = profli.nep;
tep = profli.tep;
tip = profli.tip;
%nip = profli.nip;
n1p = profli.n1p;
nD  = profli.n1p .* ((zs.nDm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
%nT  = profli.n1p .* ((zs.nTm ./ max(1,trapz(x,vpr .* abs(profli.n1p),2)) .* trapz(x,vpr,2)) * ve);
%pnbi = profli.pnbi;
zeffp = profli.zeff;

% cas with rwo NBI
pnbi_in = profli.nbishape_el + profli.nbishape_ion;
pnbi_in_tot  = trapz(profli.xli,(real(pnbi_in) + imag(pnbi_in)) .* profli.vpr,2);
pnbi_in = pnbi_in .* (((real(zs.pnbi_th) + imag(zs.pnbi_th)) ./ max(1,pnbi_in_tot)) * ve);

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
    case 5
        zj = mean((1  + 4 .* cons.iso)   ./  (1 + 2 .* cons.iso));
        aj = mean((2  + 3 .* cons.iso)   ./  (1 + cons.iso));
    case 11
        zj = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        aj = mean((1  + 11 .* cons.iso)  ./  (1 + cons.iso));
end

% gestion de la puissance NBI
if isfield(option,'nb_nbi')
    nb_nbi = option.nb_nbi;
else
    nb_nbi = 1;
end
if nargin > 4
    rapnbi  =  real(pnorme) ./ max(1,trapz(x,vpr .* real(pnbi_in),2));
    pnbi  = real(pnbi_in) ./ (max(1,trapz(x,vpr .* real(pnbi_in),2)) * ve) .* (real(pnorme) * ve);
    ftnbi = zeros(size(cons.ftnbi));
    pnbi_th = real(zs.pnbi_th) .* real(rapnbi);
    if nb_nbi == 2
        rapnbi2  =  imag(pnorme) ./ max(1,trapz(x,vpr .* imag(pnbi_in),2));
        pnbi2  = imag(pnbi_in) ./ (max(1,trapz(x,vpr .* imag(pnbi_in),2)) * ve) .* (imag(pnorme) * ve);
        ftnbi2 = zeros(size(cons.ftnbi));
        pnbi_th2 = imag(zs.pnbi_th) .* rapnbi2;
        pnbi_th = pnbi_th + sqrt(-1) .* pnbi_th2;
        pnbi    = pnbi  + sqrt(-1) .* pnbi2;
        ftnbi   = ftnbi + sqrt(-1) .* ftnbi2;
    end
else
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        ftnbi = zeros(size(cons.ftnbi));
    else
        ftnbi = real(cons.ftnbi);
    end
    rapnbi = vt;
    pnbi_th = real(zs.pnbi_th);
    pnbi  = real(pnbi_in) ./ (max(1,trapz(x,vpr .* real(pnbi_in),2)) * ve) .* (real(zs.pnbi_th) * ve);
    if nb_nbi == 2
        if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
            ftnbi2 = zeros(size(cons.ftnbi));
        else
            ftnbi2 = imag(cons.ftnbi);
        end
        rapnbi2 = vt;
        pnbi_th2 = imag(zs.pnbi_th);
        pnbi2  = imag(pnbi_in) ./ (max(1,trapz(x,vpr .* imag(pnbi_in),2)) * ve) .* (imag(zs.pnbi_th) * ve);
        pnbi_th = pnbi_th + sqrt(-1) .* pnbi_th2;
        pnbi    = pnbi  + sqrt(-1) .* pnbi2;
        ftnbi   = ftnbi + sqrt(-1) .* ftnbi2;
    end
end

% choix du minoritaire
% switch option.mino
%     case 'He3'
%         ag = 3;
%         zg = 2;
%         lg = 7.92e-3;
%     case 'T'
%         ag = 3;
%         zg = 1;
%         lg = 7.92e-3;
%     case 'He4'
%         ag = 4;
%         zg = 2;
%         lg = 4.55e-3;
%     case 'D'
%         ag = 2;
%         zg = 1;
%         lg = 6.46e-3;
%     case 'B'
%         ag = 11;
%         zg = 5;
%         lg = 4.576e-3 .* sqrt(11) / 5;
%     otherwise
%         ag = 1;
%         zg = 1;
%         lg = 4.576e-3;
% end
%correction de einj pour tenir compte de l'interaction avec icrh
%fdicrh    = 1 .* zs.nem; %?
fdicrh    = 0.25e-18 .* zs.nem.^2; %?
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    Wdt = zeros(size(zs.esup_nbi));
else
    Wdt = (1 - real(cons.ftnbi)) .* real(zs.esup_nbi) ./ zs.vp;
end
wcdt      = 15.2e6 .* zj ./ aj .* profli.fdia(:,1) ./ profli.Raxe(:,1);
if option.cmin > 0
    ficrh_nbi = fdicrh .* Wdt ./ zs.nmino ./ wcdt;
else
    ficrh_nbi = 0 .* vt;
end
ficrh_nbi = max(0,min(10,ficrh_nbi));
picrh_nbi = ficrh_nbi ./ (1 + ficrh_nbi) .* zs.picrh;
einj      = option.einj .* max(1,min(10,(picrh_nbi + real(pnbi_th)) ./ max(1,real(pnbi_th))));
pnbi_th   = picrh_nbi + pnbi_th;

if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        Wdt2 = zeros(size(zs.esup_nbi));
    else
        Wdt2  = (1 - imag(cons.ftnbi)) .* imag(zs.esup_nbi) ./ zs.vp;
    end
    if option.cmin > 0
        ficrh_nbi2 = fdicrh .* Wdt2 ./ zs.nmino ./ wcdt;
    else
        ficrh_nbi2 = 0 .* vt;
    end
    ficrh_nbi2 = max(0,min(10,ficrh_nbi2));
    picrh_nbi2 = ficrh_nbi2 ./ (1 + ficrh_nbi2) .* zs.picrh;
    einj2      = option.einj2 .* max(1,min(10,(picrh_nbi2 + pnbi_th2) ./ max(1,pnbi_th2)));
    pnbi_th    = pnbi_th + sqrt(-1) .* picrh_nbi2;
    ficrh_nbi  = ficrh_nbi + sqrt(-1) .* ficrh_nbi2;
    picrh_nbi  = picrh_nbi + sqrt(-1) .* picrh_nbi2;
    einj       = einj + sqrt(-1) .* einj2;
end
%  figure(21);clf
%  subplot(3,1,1)
%  plot(cons.temps,ficrh_nbi);
%  subplot(3,1,2)
%  plot(cons.temps,einj);
%  subplot(3,1,3)
%  plot(cons.temps,picrh_nbi,'b',cons.temps,pnbi_th,'r');
%  drawnow
%keyboard

%  figure(21)
%  clf;
%  plot(cons.temps,zs.nmino)
%  drawnow

% fusion DD du plasma thermique (correction ,doit etre petit, l'enrichissement du milieu en T, H et He3 n'est pas prise en compte)
[dd_p,dd_n,dt,dhe3,tt,the3_pn,the3_d]=zsigmavfusion(max(tip(:),13.6));
sp       = 0.5 .* nD .^ 2 .* reshape(dd_p.sv,size(tip)) .* (tip >= dd_p.timin) .*  (tip <= dd_p.timax);
sn       = 0.5 .* nD .^ 2 .* reshape(dd_n.sv,size(tip)) .* (tip >= dd_n.timin) .*  (tip <= dd_n.timax);
proton_dd_th     = trapz(x,vpr .* sp,2);
neutron_th    = trapz(x,vpr .* sn,2);

% calcul de l'interraction faisceau-plasma
pmod        = max(1,real(pnbi) + imag(pnbi));
nDi         = max(1e13,trapz(x,nD .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
tii         = max(30,trapz(x,tip .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
tei         = max(30,trapz(x,tep .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
nei         = max(1e13,trapz(x,nep .* pmod .* vpr,2) ./ max(1,trapz(x,pmod .* vpr,2)));
%taus_nbi    = 6.27e8 .* 2 .* tei .^ (3/2) ./ (nei./ 1e6) ./ 17;
switch gas_nbi
    case 5
        fact        = (zs.nDm ./2 + zs.nTm ./ 3 + (zs.n1m - zs.nTm - zs.nDm) + (4/3) .* zs.nhem + option.frhe0 .* zs.nem +  ...
                       zs.nimpm  .* (option.zimp .^ 2 ./ aimp + option.rimp .* option.zmax .^ 2 ./ amax)) ./ zs.nem;
    case 11
        fact        = (zs.nDm ./2 + 25/11 .* zs.nTm + (zs.n1m - zs.nDm) + zs.nhem  + ...
                       zs.nimpm  .* (option.zimp .^ 2 ./ aimp + option.rimp .* option.zmax .^ 2 ./ amax)) ./ zs.nem;
   otherwise
        fact        = (zs.nDm ./2 + zs.nTm ./ 3 + (zs.n1m - zs.nTm - zs.nDm) + zs.nhem  + ...
                      zs.nimpm  .* (option.zimp .^ 2 ./ aimp + option.rimp .* option.zmax .^ 2 ./ amax)) ./ zs.nem;
end
lnldei      = 15.2 - 0.5 .* log(nei./1e20) + log(tei ./1e3);

%ecrit_nbi   = max(30,14.8 .* tei .* (2 .^ (3/2)  .* fact) .^ (2/3));   % c'est l'energie liee a vc
switch gas_nbi
    case -1
        % c'est l'energie liee a vc
        ecrit_nbi = max(30,14.8 .* tei .* (1 .^ (3/2)  .* fact) .^ (2/3));
        taus_nbi  = 6.27e8  .* 1 .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
    case 11
        z_ave_2  =  1 .* (1-real(ftnbi)) +  25 .* real(ftnbi);
        z_ave_43 =  1 .* (1-real(ftnbi)) +  5 ^(4/3) .* real(ftnbi);
        % c'est l'energie liee a vc
        ecrit_nbi= max(30,14.8 .* tei .* ((1 .* (1-real(ftnbi)) + 11 .* real(ftnbi)).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
        taus_nbi       = 6.27e8 .* (1 .* (1-real(ftnbi)) + 11 .* real(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei ./ z_ave_2;
    case 5
        z_ave_2  =  1 .* (1-real(ftnbi)) +  4 .* real(ftnbi);
        z_ave_43 =  1 .* (1-real(ftnbi)) +  2 ^(4/3) .* real(ftnbi);
        % c'est l'energie liee a vc
        ecrit_nbi= max(30,14.8 .* tei .* ((2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
        taus_nbi       = 6.27e8 .* (2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei ./ z_ave_2;
    case 3
        % c'est l'energie liee a vc
        ecrit_nbi = max(30,14.8 .* tei .* ((2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
        taus_nbi       = 6.27e8 .* (2 .* (1-real(ftnbi)) + 3 .* real(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
    otherwise
        % c'est l'energie liee a vc
        ecrit_nbi = max(30,14.8 .* tei .* ((2 .* (1-real(ftnbi)) + 1 .* real(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
        taus_nbi  = 6.27e8 .* (2 .* (1-real(ftnbi)) + 1 .* real(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
end



if nb_nbi == 2
    switch gas_nbi
        case -1
            % c'est l'energie liee a vc
            ecrit_nbi2 = max(30,14.8 .* tei .* (1 .^ (3/2)  .* fact) .^ (2/3));
            taus_nbi2  = 6.27e8 .* 1  .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
        case 11
            z_ave_2  =  1 .* (1-imag(ftnbi)) +  25 .* imag(ftnbi);
            z_ave_43 =  1 .* (1-imag(ftnbi)) +  5 ^(4/3) .* imag(ftnbi);
            % c'est l'energie liee a vc
            ecrit_nbi2 = max(30,14.8 .* tei .* ((1 .* (1-imag(ftnbi)) + 11 .* imag(ftnbi)).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
            taus_nbi2       = 6.27e8 .* (1 .* (1-imag(ftnbi)) + 11 .* imag(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei ./ z_ave_2;
        case 5
            z_ave_2  =  1 .* (1-imag(ftnbi)) +  4 .* imag(ftnbi);
            z_ave_43 =  1 .* (1-imag(ftnbi)) +  2 ^(4/3) .* imag(ftnbi);
            % c'est l'energie liee a vc
            ecrit_nbi2 = max(30,14.8 .* tei .* ((2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)).^ (3/2)  .* fact) .^ (2/3)) .* z_ave_43;
            taus_nbi2       = 6.27e8 .* (2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei ./ z_ave_2;
        case 3
            % c'est l'energie liee a vc
            ecrit_nbi2 = max(30,14.8 .* tei .* ((2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
            taus_nbi2       = 6.27e8 .* (2 .* (1-imag(ftnbi)) + 3 .* imag(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
        otherwise
            % c'est l'energie liee a vc
            ecrit_nbi2 = max(30,14.8 .* tei .* ((2 .* (1-imag(ftnbi)) + 1 .* imag(ftnbi)).^ (3/2)  .* fact) .^ (2/3));
            taus_nbi2  = 6.27e8 .* (2 .* (1-imag(ftnbi)) + 1 .* imag(ftnbi)) .* tei .^ (3/2) ./ (nei./ 1e6) ./ lnldei;
    end
    ecrit_nbi  = ecrit_nbi + sqrt(-1) .* ecrit_nbi2;
    taus_nbi   = taus_nbi  + sqrt(-1) .* taus_nbi2;
end
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    [neutron_nbi_th,emean] = znbi_dd(nDi,tii,real(pnbi_th) .* 0,real(taus_nbi),real(einj),real(ecrit_nbi));
else
    [neutron_nbi_th,emean] = znbi_dd(nDi,tii,real(pnbi_th) .* (1 - real(ftnbi)),real(taus_nbi),real(einj),real(ecrit_nbi));
end
emean                  = max(zs.tem,emean);
if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        [neutron_nbi_th2,emean2] = znbi_dd(nDi,tii,imag(pnbi_th) .* 0,imag(taus_nbi),imag(einj),imag(ecrit_nbi));
    else
        [neutron_nbi_th2,emean2] = znbi_dd(nDi,tii,imag(pnbi_th) .* (1 - imag(ftnbi)),imag(taus_nbi),imag(einj),imag(ecrit_nbi));
        
    end
    emean2                   = max(zs.tem,emean2);
    neutron_nbi_th           = neutron_nbi_th + neutron_nbi_th2;
    emean                    = emean + sqrt(-1) .* emean2;
end

% approximation faisceau-faisceau avec la meme formule (c'est faux)
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    [esupranbi,pDth,tauseff] = zsupra0(temps,real(pnbi_th) .* 0,real(taus_nbi),1e-6*ones(size(taus_nbi)),real(ecrit_nbi),real(einj),2);
else
    [esupranbi,pDth,tauseff] = zsupra0(temps,real(pnbi_th) .* (1 - real(ftnbi)),real(taus_nbi),1e-6*ones(size(taus_nbi)),real(ecrit_nbi),real(einj),2);
end
fon  = ((real(pnbi_th) - real(picrh_nbi)) ./ zs.vp) > 1e3;
if any(fon)
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        ctr  =  0 .* real(tauseff) ./ 2 ./ zs.vp ./ (1.602176462e-19 .* real(emean));
    else
        ctr  = (1 - real(ftnbi)).* real(tauseff) ./ 2 ./ zs.vp ./ (1.602176462e-19 .* real(emean));
    end
    ctr  = fon .* ctr + (~fon) .* mean(ctr(find(fon)));
    ttr  = fon .* real(taus_nbi) + (~fon) .* mean(real(taus_nbi(find(fon))));
    neutron_nbi_nbi   = 0.5 .* znbi_dd(real(pnbi_th) .* ctr, real(emean),pDth,ttr,real(einj),real(ecrit_nbi));
else
    neutron_nbi_nbi = zeros(size(cons.temps));
end

neutron_nbi_nbi_1 = neutron_nbi_nbi;
if nb_nbi == 2
    if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
        [esupranbi2,pDth2,tauseff2] = zsupra0(temps,imag(pnbi_th) .* 0,imag(taus_nbi),1e-6*ones(size(taus_nbi)),imag(ecrit_nbi),imag(einj),2);
    else
        [esupranbi2,pDth2,tauseff2] = zsupra0(temps,imag(pnbi_th) .* (1 - imag(ftnbi)),imag(taus_nbi),1e-6*ones(size(taus_nbi)),imag(ecrit_nbi),imag(einj),2);
    end
    fon2  = ((imag(pnbi_th) - imag(picrh_nbi)) ./ zs.vp) > 1e3;
    if any(fon2)
        if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
            ctr2  = 0.* tauseff2 ./ 2 ./ zs.vp ./ (1.602176462e-19 .* imag(emean));
        else
            ctr2  = (1 - imag(ftnbi)).* tauseff2 ./ 2 ./ zs.vp ./ (1.602176462e-19 .* imag(emean));
        end
        ctr2  = fon2 .* ctr2 + (~fon2) .* mean(ctr2(find(fon2)));
        ttr2  = fon2 .* imag(taus_nbi) + (~fon2) .* mean(imag(taus_nbi(find(fon2))));
    end
    % 2 sur 2
    if any(fon2)
        neutron_nbi_nbi_2   = 0.5 .* znbi_dd(imag(pnbi_th) .* ctr2,imag(emean),pDth2,ttr2,imag(einj),imag(ecrit_nbi));
    else
        neutron_nbi_nbi_2   = 0 * neutron_nbi_nbi_1;
    end
    
    if any(fon2) && any(fon)
        % 1 sur 2
        neutron_nbi_nbi_3   = 0.5 .* znbi_dd(real(pnbi_th) .* ctr,real(emean),pDth2,ttr2,imag(einj),imag(ecrit_nbi));
        % 2 sur 1
        neutron_nbi_nbi_4   = 0.5 .* znbi_dd(imag(pnbi_th) .* ctr2,imag(emean),pDth,ttr,real(einj),real(ecrit_nbi));
    else
        neutron_nbi_nbi_3   = 0 * neutron_nbi_nbi_1;
        neutron_nbi_nbi_4   = 0 * neutron_nbi_nbi_1;
    end
    esupranbi    = esupranbi + sqrt(-1) .* esupranbi2;
    pDth         = pDth + sqrt(-1) .* pDth2;
    tauseff      = tauseff  + sqrt(-1) .* tauseff2;
    neutron_nbi_nbi = neutron_nbi_nbi_1 + neutron_nbi_nbi_2 + neutron_nbi_nbi_3 + neutron_nbi_nbi_4;
    
end

%  figure(21);clf
%  subplot(3,1,1)
%  plot(cons.temps,ctr);
%  subplot(3,1,2)
%  plot(cons.temps,ttr);
%  subplot(3,1,3)
%  plot(cons.temps,neutron_nbi_nbi);
%  drawnow



% calcul exact pour un temps
if isfield(option,'forced_H_NBI') && (option.forced_H_NBI ~= 0)
    % not needed
elseif (max((real(pnbi_th) + imag(pnbi_th)) ./ max(1,zs.vp)) > 1000) && (any(pnbi(:) ~= 0))
    %indr		     = find(pnbi_th == max(pnbi_th),1);
    pnbi_th_sum = (real(pnbi_th) .* sqrt(real(einj)) + imag(pnbi_th) .* sqrt(imag(einj))) ./ ...
        max(1,sqrt(real(einj)) + sqrt(imag(einj)));
    pnbi_sum    = (real(pnbi) .* (sqrt(real(einj)) * ve) + imag(pnbi) .* (sqrt(imag(einj)) * ve)) ./ ...
        (max(eps,sqrt(real(einj)) + sqrt(imag(einj))) * ve);
    vpr_m = zpmv(cons.temps,pnbi_th_sum,vpr);
    nep_m = zpmv(cons.temps,pnbi_th_sum,nep);
    tep_m = zpmv(cons.temps,pnbi_th_sum,tep);
    n1p_m = zpmv(cons.temps,pnbi_th_sum,n1p);
    tip_m = zpmv(cons.temps,pnbi_th_sum,tip);
    zeffp_m = zpmv(cons.temps,pnbi_th_sum,zeffp);
    pnbi_m  = zpmv(cons.temps,pnbi_th_sum,pnbi_sum);
    nDm_m   = zpmean(cons.temps,pnbi_th_sum,zs.nDm);
    nTm_m   = zpmean(cons.temps,pnbi_th_sum,zs.nTm);
    n1m_m   = zpmean(cons.temps,pnbi_th_sum,zs.n1m);
    nhem_m  = zpmean(cons.temps,pnbi_th_sum,zs.nhem);
    nimpm_m = zpmean(cons.temps,pnbi_th_sum,zs.nimpm);
    nem_m   = zpmean(cons.temps,pnbi_th_sum,zs.nem);
    neutron_ref  = z0beambeam(1,x,vpr_m,nep_m,tep_m,n1p_m,tip_m,zeffp_m,pnbi_m, ....
        nDm_m,nTm_m,n1m_m,nhem_m,nimpm_m,nem_m, ...
        zpmean(cons.temps,real(pnbi_th),real(einj)), ...
        zpmean(cons.temps,real(pnbi_th),real(einj)), ...
        zpmean(cons.temps,real(pnbi_th),real(pnbi_th) .* (1 - real(ftnbi))), ...
        zpmean(cons.temps,real(pnbi_th),real(pnbi_th) .* (1 - real(ftnbi))), ...
        zpmean(cons.temps,real(pnbi_th),real(zs.mu0_nbi)), ...
        zpmean(cons.temps,real(pnbi_th),real(zs.mu0_nbi)),1);
    
    
    if nb_nbi == 2
        neutron_ref_1 = neutron_ref;
        neutron_nbi_nbi   = neutron_nbi_nbi_1   ./ max(1,zpmean(cons.temps,real(pnbi_th),neutron_nbi_nbi_1)) .* ...
            neutron_ref_1;
        if any(neutron_nbi_nbi_2 ~= 0)
            neutron_ref_2  = z0beambeam(1,x,vpr_m,nep_m,tep_m,n1p_m,tip_m,zeffp_m,pnbi_m, ....
                nDm_m,nTm_m,n1m_m,nhem_m,nimpm_m,nem_m, ...
                zpmean(cons.temps,imag(pnbi_th),imag(einj)), ...
                zpmean(cons.temps,imag(pnbi_th),imag(einj)), ...
                zpmean(cons.temps,imag(pnbi_th),imag(pnbi_th) .* (1 - imag(ftnbi))), ...
                zpmean(cons.temps,imag(pnbi_th),imag(pnbi_th) .* (1 - imag(ftnbi))), ...
                zpmean(cons.temps,imag(pnbi_th),imag(zs.mu0_nbi)), ...
                zpmean(cons.temps,imag(pnbi_th),imag(zs.mu0_nbi)),1);
            neutron_nbi_nbi = neutron_nbi_nbi + neutron_nbi_nbi_2   ./...
                max(1,zpmean(cons.temps,imag(pnbi_th),neutron_nbi_nbi_2)) .* neutron_ref_2;
            
        end
        if any(neutron_nbi_nbi_3 ~= 0) || any(neutron_nbi_nbi_4 ~= 0)
            neutron_ref_34  = z0beambeam(1,x,vpr_m,nep_m,tep_m,n1p_m,tip_m,zeffp_m,pnbi_m, ....
                nDm_m,nTm_m,n1m_m,nhem_m,nimpm_m,nem_m, ...
                zpmean(cons.temps,real(pnbi_th_sum),real(einj)), ...
                zpmean(cons.temps,imag(pnbi_th_sum),imag(einj)), ...
                zpmean(cons.temps,real(pnbi_th_sum),real(pnbi_th) .* (1 - real(ftnbi))), ...
                zpmean(cons.temps,imag(pnbi_th_sum),imag(pnbi_th) .* (1 - imag(ftnbi))), ...
                zpmean(cons.temps,real(pnbi_th_sum),real(zs.mu0_nbi)), ...
                zpmean(cons.temps,imag(pnbi_th_sum),imag(zs.mu0_nbi)),0);
            neutron_nbi_nbi = neutron_nbi_nbi + (neutron_nbi_nbi_3 + neutron_nbi_nbi_4)   ./...
                max(1,zpmean(cons.temps,imag(pnbi_th_sum),neutron_nbi_nbi_3 + neutron_nbi_nbi_4)) .* neutron_ref_34;
            
        end
    else
        neutron_nbi_nbi   = neutron_nbi_nbi   ./ max(1,zpmean(cons.temps,pnbi_th,neutron_nbi_nbi)) .* neutron_ref;
    end
end

% autre calcul a partir de J.D Strachan NF 33,7 (1993) p 991-
% normalisation beam beam sur beam-plasma
%sbt  = 0.167e-9 .* nDi .* pnbi_th.* (1 -ftnbi) .* taus_nbi;
%sbb  = 1.3e13 .* (pnbi_th.* (1 -ftnbi)) .^ 2 .* taus_nbi .^ 2 ./ option.einj .^2;
%disp([mean(neutron_nbi_th),mean(neutron_nbi_nbi)]);
%disp([mean(sbt),mean(sbb)]);
%disp([mean(sbb ./ sbt .* neutron_nbi_th),mean(sbb),mean(neutron_nbi_nbi)]);
%neutron_nbi_nbi = sbb ./ sbt .* neutron_nbi_th;


% nombre de neutron total
neutron_th    = max(0,neutron_th);
neutron_nbi_th    = max(0,neutron_nbi_th);
neutron_nbi_nbi   = max(0,neutron_nbi_nbi);
neutron_total =neutron_th + neutron_nbi_th +  neutron_nbi_nbi;

% puissance de fusion complementaire
% le rapport entre les voies n et p est celui du thermique (en gros 2 * nombre de neutrons)
pddfus    = 1.602176462e-19 .* ((dd_p.p + dd_p.t) .* proton_dd_th ./ max(1,neutron_th)  + dd_n.he3) .* neutron_total;
proton_dd.th      = proton_dd_th;
proton_dd.nbi_th  = proton_dd_th ./ max(1,neutron_th) .* neutron_nbi_th;
proton_dd.nbi_nbi = proton_dd_th ./ max(1,neutron_th) .* neutron_nbi_nbi;
proton_dd.total   = proton_dd_th ./ max(1,neutron_th) .* neutron_total;


%    figure(21);clf
%    plot(cons.temps,neutron_th,'g',cons.temps,neutron_nbi_th,'b',...
%    cons.temps,neutron_nbi_nbi,'m', cons.temps,neutron_total,'c', ...
%    evalin('base','z0dinput.cons.temps'),evalin('base','z0dinput.exp0d.ndd'),'r');
%    drawnow
%keyboard



function  s  = zpmv(t,p,e)

%  figure(21);clf
%  subplot(2,1,1)
%  plot(t,p)
%  subplot(2,1,2);
%  plot(t,e);
%  drawnow
ve = ones(1,size(e,2));

indok = find(isfinite(p) & all(isfinite(e),2));
if length(indok) < 3
    indok = find(all(isfinite(e),2));
    if isempty(indok)
        s = NaN*ve;
        disp('zpmean : invalid data ...');
    else
        s = mean(e(indok,:),1);
    end
else
    t = t(indok);
    p = p(indok);
    e = e(indok,:);
    indp = find(p > 0);
    %disp(length(indp))
    if length(indp) > 2
        s = trapz(t(indp),(p(indp) * ve) .* e(indp,:),1) ./ (trapz(t(indp),eps + p(indp),1) * ve);
    else
        s = mean(e,1);
    end
end


function  s  = zpmean(t,p,e)

%  figure(21);clf
%  subplot(2,1,1)
%  plot(t,p)
%  subplot(2,1,2);
%  plot(t,e);
%  drawnow

indok = find(isfinite(p) & isfinite(e));
if length(indok) < 3
    indok = find(isfinite(e));
    if isempty(indok)
        s = NaN;
        disp('zpmean : invalid data ...');
    else
        s = mean(e(indok));
    end
else
    t = t(indok);
    p = p(indok);
    e = e(indok);
    indp = find(p > 0);
    %disp(length(indp))
    if length(indp) > 2
        s = trapz(t(indp),p(indp) .* e(indp)) ./ trapz(t(indp),eps + p(indp));
    else
        s = mean(e);
    end
end



