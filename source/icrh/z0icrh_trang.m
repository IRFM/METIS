% cette fonction calcul les effets du chauffage minoriatire a l'harmonique 1 de ICRH
% pour les test : [pth,pel,pion,esupra,einj,ecrit,teff0,taus,frloss,rres,xres,frac,harm,nmino,Pe,PM,esupra_par,jmino,imino,pfw] = z0icrh_trang(post.zerod,post.z0dinput.geo,post.z0dinput.cons,post.z0dinput.option,post.profil0d);
function [pth,pel,pion,esupra,einj,ecrit,teff0,taus,frloss,rres,xres,frac,harm,nmino,Pe,PM,esupra_par,jmino,imino,pfw,pabs] = z0icrh_trang(zs,geo,cons,option,profli)

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
dd  = (btot - bres * ones(size(btot))) .^ 2;
mask = (dd == (min(dd,[],2)*ones(1,size(xa,2))));
res   = sum(rr .* mask,2);
rres = res; % pour la sortie
xres  = min(0.95,abs(sum(xa .* mask,2)) ./ geo.a);

% choix de l'harmonique
switch option.mino
    case 'D'
        
        if (option.gaz == 2)  || (option.gaz == 3)
            % D dans D, recherche de l'harmonique superieure
            dd2 = abs(btot(:,21) - bres ./ 2);
            dd3 = abs(btot(:,21) - bres ./ 3);
            maskh = (dd3 < dd2);
            harm  =  round(mean(cons.picrh .* (2 + maskh))./ max(1,mean(cons.picrh)));
        else
            bmax = trapz(cons.temps, cons.picrh .* max(btot,[],2)) ./ trapz(cons.temps, cons.picrh);
            bmin = trapz(cons.temps, cons.picrh .* min(btot,[],2)) ./ trapz(cons.temps, cons.picrh);
            if (bres < bmax)  &&  (bres > bmin)
                harm = 1;
            elseif ((bres/2) < bmax)  &&  ((bres/2) > bmin)
                harm = 2;
            elseif((bres/3) < bmax)  &&  ((bres/3) > bmin)
                harm = 3;
            else
                harm = 0;
            end
        end        
    otherwise
        bmax = trapz(cons.temps, cons.picrh .* max(btot,[],2)) ./ trapz(cons.temps, cons.picrh);
        bmin = trapz(cons.temps, cons.picrh .* min(btot,[],2)) ./ trapz(cons.temps, cons.picrh);
        if (bres < bmax)  &&  (bres > bmin)
            harm = 1;
        elseif ((bres/2) < bmax)  &&  ((bres/2) > bmin)
            harm = 2;
        elseif((bres/3) < bmax)  &&  ((bres/3) > bmin)
            harm = 3;
        else
            harm = 0;
        end
end

% choix du minoritaire
% compatibility with z0ich.m
switch option.gaz
    case 1
        nref = max(1e13,zs.n1m - zs.nDm -zs.nTm);
    case {2,3}
        nref = zs.nDm;
    case 5
        nref = zs.nhem;        
    case 11
        nref = max(1e13,zs.n1m - zs.nDm);
    otherwise
        nref = zs.nhem;
end


switch option.mino
    case 'T'
        nmino  = nref  .* max(1e-2,cons.iso .* option.cmin);
    case 'He'
        nmino  = nref  .* option.cmin;
    otherwise
        nmino  = nref .* option.cmin;
end
cmino = nmino ./ zs.nim;

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



% puissance
% cas du ripple
prip   = max(0,(option.rip == 1) .* (cons.picrh  .* (1 - zs.frloss_icrh) - zs.picrh));
pm     = max(0,cons.picrh .* (1 - zs.frloss_icrh) - prip);
pfw    = pm .* abs(option.fabs_fw);
pm     = pm - pfw;

if all(pm <= 1e3) || ~isfield(profli,'qjli') || (harm == 0)
    % if harm == 0, no ressonnace in plasma
    if harm == 0
        pm = eps .* ones(size(pm));
        fprintf('o');
    end
    % gain de temps
    pth = pm;
    pel = pm ./ 2;
    pion = pm ./ 2;
    esupra = zeros(size(pm));
    esupra_par  = zeros(size(pm));
    einj  = 2 .* zs.te0;
    teff0 = zs.te0 .* zs.tite;
    taus =  6.27e8 .* ag ./ zg .^ 2 ./ 17 ./ (zs.nem ./ 1e6).* zs.tem .^ (3/2);
    ecrit = 14.8 .* zs.tem .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* nref ./ zs.nem) .^ (2/3);
    frloss = 0 .* zs.te0;
    harm = harm .* vt;
    frac = zeros(size(pm));
    shape = ones(size(pm)) * (1 - x .^ 2) .^ 4;
    vpr   = (4 .* pi .* geo.a .^ 2 .* geo.R .* geo.K) * x;
    shape = shape .* ((pm ./ trapz(x,vpr .* shape,2)) * ones(size(x)));
    Pe    = shape ./ 2;
    PM    = shape ./ 2;
    jmino = 0 .* Pe;
    imino = zeros(size(pm));
    return
end

% appel du code de Trang
% % preparation
% tune.ichcd.active = 1;
% tune.ichcd.uindices = 2;
% % ICRH_config
% ipa.active = tune.ichcd.active;
% ipa.uindices = tune.ichcd.uindices;
% ipa.check = false;

%% METIS Test Find all in  post.z0dinput.option
ipa.eps_Xi = 1e-3;
ipa.upshift = 0;
ipa.p  = harm;
ipa.f  = option.freq*1e6; % Hz 5.5e10; 5.5e7;
ipa.fz = option.icrh_width;
ipa.rate_ni =  cmino; % nm/nM
ipa.n = option.nphi;
% ipa.widthtau = 1;2; 1.22; %*delR == post.z0dinput.option.icrh_width
ipa.P_lossTau = 0;
ipa.Nu = 1000;
ipa.Nz = 100;
% In model.atom
ipa.Zm = zg;% H in METIS;
ipa.ZM = zj;% Deuteron 1 % charge of Majority with electron
ipa.Am = ag;% number of mass of minority
ipa.AM = aj;% number of mass of Majority

ichcd_params = ipa;
g = 0; % geometric params like C2, C3,...
model.geom.rhogauss = profli.xli';
% rho = model.geom.rhogauss;
Pe = zeros(size(profli.xli,2),size(cons.temps,1));
PM = Pe;
% Pm = Pe;
Wperp = Pe;
Wpar = Pe;
Tpar = Pe;
pabs = Pe;
rres = geo.R;
P_loss_conv = zeros(size(cons.temps));
flag_non_conv = zeros(size(cons.temps));

% boucle sur les temps
for it = 1:length(cons.temps)
    U = pm(it);
    
    if U > 1e3 % for icrh
        q  = profli.qjli(it,:)';
        te = profli.tep(it,:)';
        ne = profli.nep(it,:)';
        ti = profli.tip(it,:)';
        ni = profli.nip(it,:)';
        % to take into the shafranov shift we use Raxe(xres) instead of R0
        rigidity = geo.b0(it) .* geo.R(it);
        model.equi.R0 = geo.R(it);
        model.equi.epsilon = geo.a(it)/model.equi.R0;
        model.equi.B0 = rigidity ./ model.equi.R0;
        model.equi.kappa  = geo.K(it);
        model.equi.Raxe   = profli.Raxe(it,:)';
        model.equi.xRaxe  = profli.xli';
        ichcd_params.picrh = U;
        ichcd_params.Ploss = 0; % already acounted u=in pm computation
        if length(cmino) > 1
            ichcd_params.rate_ni =  cmino(it);
        end
        if length(aj) > 1
            ichcd_params.AM =  aj(it);
        end
        [Pe(:,it),PM(:,it),pabs(:,it),rres(it),Wperp(:,it),Wpar(:,it),Tpar(:,it),P_loss_conv(it),flag_non_conv(it),ichcd_params.Xi] = ...
            ichcd_colli_test(q,te,ne,ti,ni,U,model,ichcd_params,g);
        if flag_non_conv(it) == 1
            ichcd_params.Xi = 0.02 .* ones(size(ichcd_params.Xi));
        end
        
        % Metis have only two channels for heat power transport: ions and electrons. heat that goes directly to minority ions must be account in ion channel
        PM(:,it) = PM(:,it) + max(0,pabs(:,it) - Pe(:,it) - PM(:,it));
        pabs(:,it) = max(pabs(:,it),PM(:,it) + Pe(:,it));
        
        %TO BROADEN THE FAST ION SOURCE PROFILE ACCORDING TO ORBIT  WIDTH EFFECTS
        if option.orbit_width == 1
            einj_loc   = (2/3 .* (Wperp(:,it) +  Wpar(:,it))' + (3/2) .* phys.e .* ti' .* ni'  .* ichcd_params.rate_ni) ...
                ./ max(1,ni' .* ichcd_params.rate_ni) ./ phys.e;
            equi.vpr    = profli.vpr(it,:);      % VOLUME ELEMENT PROFILE                        (M2)
            equi.q      = profli.qjli(it,:);        % SECURITY FACTOR PROFILE                       (-)
            equi.ftrap  = profli.ftrap(it,:);    % PROFILE OF FRACTION OF TRAPPED PARTICLES      (-)
            equi.a      = geo.a(it) .* profli.xli;        % PROFILE OF MINOR RADIUS, (RMAX-RMIN)/2        (M)
            equi.raxe   = profli.Raxe(it,:);     % PROFILE OF FLUX SURFACE CENTRE, (RMAX+RMIN)/2 (M)
            equi.b2     = (profli.fdia(it,:) ./ equi.raxe) .^ 2 + profli.bpol(it,:) .^ 2; % SQUARE OFMAGNETIC FIELD PROFILE                        (T^2)
            %
            Pe(:,it) = zorbitel2(profli.xli,einj_loc,ag,zg,Pe(:,it)',equi,phys)';
            PM(:,it) = zorbitel2(profli.xli,einj_loc,ag,zg,PM(:,it)',equi,phys)';
            %figure(22);plot(profli.xli,einj_loc);ichcd_params.rate_ni
            %figure(21);plot(profli.xli',Pe_mem,'c',profli.xli',Pe(:,it),'b',profli.xli',PM_mem,'m',profli.xli',PM(:,it),'r');drawnow
        end
        
    end
    
end
if any(flag_non_conv == 1)
    nc = fix(10 .* sum(flag_non_conv) ./ length(flag_non_conv));
    if nc > 9
        fprintf('Y');
    elseif nc == 0
        fprintf('y');
    else
        fprintf('%d',nc);
    end
end

if length(cmino) > 1
    cmino = cmino * ve;
end
% energy suprathermique
esupra_perp = trapz(profli.xli,profli.vpr .* Wperp',2);
esupra_par  = trapz(profli.xli,profli.vpr .* Wpar',2);
esupra_st   = esupra_perp + esupra_par;
etot_st     = esupra_st + (3/2) .* trapz(profli.xli,profli.vpr .* pabs' .* profli.tip .*  profli.nip .* cmino,2) ./  ...
    max(1,trapz(profli.xli,profli.vpr .* pabs',2)) .* phys.e;

% calcul de l'energie moyenne des ions
einj   = trapz(profli.xli,profli.vpr .* (Wperp' + Wpar' + (3/2) .* phys.e .* profli.tip .* profli.nip  .* cmino) .* pabs',2) ./ ...
    max(1,trapz(profli.xli,profli.vpr .* pabs' .* profli.nip .* cmino,2)) ./ phys.e;
%  ecrit  = trapz(profli.xli,profli.vpr .* (Wpar' + (3/2) .* phys.e .* profli.tip .* profli.nip .* cmino) .* pabs',2) ./ ...
%           max(1,trapz(profli.xli,profli.vpr .* pabs' .* profli.nip .* cmino,2)) ./ phys.e;
if length(aj) > 1
    ecrit_mat = 14.8 .* profli.tep  .* (ag .^ (3/2) ./ (aj*ve) .* zj .^ 2 .* profli.nip .* cmino ./ profli.nep) .^ (2/3);
else
    ecrit_mat = 14.8 .* profli.tep  .* (ag .^ (3/2) ./ aj .* zj .^ 2 .* profli.nip .* cmino ./ profli.nep) .^ (2/3);
end
ecrit =trapz(profli.xli,profli.vpr .* ecrit_mat  .* pabs',2) ./ max(1,trapz(profli.xli,profli.vpr .* pabs',2));
%figure(21);plot(cons.temps,ecrit_alt,'b',cons.temps,ecrit,'r');drawnow

%% effective temperature (averaged on profile deposition
teff0  = trapz(profli.xli,profli.vpr .* Tpar' .* pabs' .* profli.nip .* cmino,2) ./ ...
    max(1,trapz(profli.xli,profli.vpr .* pabs' .* profli.nip .* cmino,2));

%% minority density
nmino  = trapz(profli.xli,profli.vpr .* pabs' .* profli.nip .* cmino,2) ./ ...
    max(1,trapz(profli.xli,profli.vpr .* pabs',2));
%% fraction of accelerated minirity ion
%  poids  = pabs' ./ (max(1,trapz(profli.xli,profli.vpr .* pabs',2)) * ve);
%  frac   = trapz(profli.xli,profli.vpr .* poids .* profli.nip .* cmino,2) ./ ...
%           max(1,trapz(profli.xli,profli.vpr .* profli.nip .* cmino,2));
%% alternative computation of fraction
frac = (etot_st ./ max(eps,einj) ./ phys.e) ./ max(1,trapz(profli.xli,profli.vpr .* profli.nip .* cmino,2));

% puissances totales
pel_st = trapz(profli.xli,profli.vpr .* Pe',2);
pion   = trapz(profli.xli,profli.vpr .* PM',2);
pabs_int  = trapz(profli.xli,profli.vpr .* pabs',2);
% nous supposons que la puissance non absorbee va directement aux electrons
%figure(21);clf;plot(cons.temps,pel_st,'r',cons.temps,pion,'b',cons.temps,pabs_int,'c',cons.temps,pfw,'g',cons.temps,pm,'m');drawnow
pfw    = pfw + max(0,pm - pabs_int);
%  zzz = flag_non_conv;
%  zzz(flag_non_conv ==  1) = NaN;
%  figure(21);
%  subplot(2,1,1)
%  plot(cons.temps,pel_st + pion + pfw,'b',cons.temps,pm,'-.r',cons.temps,pabs_int,'g',cons.temps,zzz,'ok');
%  subplot(2,1,2)
%  plot(cons.temps,fact_power);
%  drawnow
%figure(21);plot(profli.xli,Pe','b',profli.xli,PM','r');drawnow

% taus par inversion de la formule
% pel_st   = max(0,2 .* esupra_st ./ max(eps,taus) .* (pm >= 1e3));
taus = 2 .* esupra_st ./ max(eps,pel_st);
taus(flag_non_conv == 1) = min(taus(flag_non_conv == 1),zs.taue(flag_non_conv == 1));
taus = min(10 .* zs.taue,taus);

if isfield(option,'transitoire') && (option.transitoire == 1)
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
    
    %     dbstop if all error
    %     figure(22);
    %     clf
    %     plot(cons.temps,cumtrapz(cons.temps,pel_st),'r',cons.temps,cumtrapz(cons.temps,pel),'b');
    %     drawnow
    
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
% Part of IBW from MC that heat directly the electrons
if isfield(option,'MC_onoff') 
    switch option.MC_onoff
        case 'off'
            P_loss_conv(:) = 0;
    end
end
pconv = min(pion,(1 - option.fMC_loss) .* P_loss_conv);
P_loss_conv = P_loss_conv - pconv;
Pe =  Pe   .* (ones(size(Pe,1),1) *(pel ./ max(eps,trapz(profli.xli,profli.vpr .* Pe',2)))') +  ...
      pabs .* (ones(size(pabs,1),1) *(pconv ./ max(eps,trapz(profli.xli,profli.vpr .* pabs',2)))');
pel  = pel  + pconv;
pion = pion - pconv;

% securite simulation courte
pel  = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pel);
pth  = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pth);
pion = min(max(cons.picrh .* (1 - zs.frloss_icrh)),pion);


% estimation des pertes de premiere orbite
if isfield(profli,'Raxe')
    rloc     = profli.Raxe + geo.a * profli.xli;
    d0ux     = profli.Raxe - geo.R * ve;
else
    rloc   = geo.R * ve + geo.a * profli.xli;
    d0ux   = 0;
end
btot = (geo.b0 * ones(1,size(profli.xli,2))) .* sqrt(((geo.a * profli.xli) ./ qp ./ rloc) .^ 2 +  ...
    ((geo.R *ones(1,size(profli.xli,2))) ./ rloc) .^ 2);

% only isotrope part
einj_loc   = ((2/3) .* (Wperp + Wpar)' + (3/2) .* phys.e .* profli.tip .* profli.nip  .* cmino) ./ max(1,profli.nip .* cmino) ./ phys.e;
ral    = lg .* sqrt(einj_loc ./ 1e3) ./ btot; % en m
% seul l'injection contre courant envoie les particules vers l'exterieur
% patato
dp1     = rloc .* (2 .* qp .* ral ./ rloc) .^ (2/3);
% banana
dp2     = sqrt(geo.a * profli.xli ./ rloc) .* ral .* qp;
dp      = dp2 .* (dp2 < (geo.a * profli.xli)) + dp1 .* (dp2 >= (geo.a * profli.xli));
mask   = (ral + dp + d0ux + geo.a * profli.xli - geo.a  * ve) >0;
%frloss = max(0,min(1,(P_loss_conv + trapz(profli.xli,profli.vpr .* pabs' .* mask,2)) ./ (trapz(profli.xli,profli.vpr .* pabs',2)+eps)));
frloss = max(0,min(1,(P_loss_conv + trapz(profli.xli,profli.vpr .* pabs' .* mask,2)) ./ max(pm,eps)));

% fast particles induced current
n_fast = zerod_fast_ions_density(profli.nep,profli.tep,profli.zeff,zs.meff * ve,ag,zg,einj_loc,Wperp' + Wpar');
% calcul du courant genere par les alpha (ref L-G Eriksson and F. Porcelli, Plasma Phys. Controlled Fus. 43 (2001) R145-)
% formule p R175, 9.4
dp0R    = dp ./ rloc;
vmino   = sqrt(2 .*  einj_loc .* phys.e ./  phys.ua ./ ag);
jmino   = (1/2) .* abs(option.ifast_icrh) .* zg .* phys.e .* vmino .* sqrt(dp0R) .* n_fast;

% calcul de l'effet du courant de retour
switch option.e_shielding
    case 'Honda-NEO'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        % ref A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664
        %GZ = z0sauterL31(xp,tep,nep,qpr,zeffp,Raxe,ftrap,epsi);
        [~,~,GZ] = z0etaboot_neofit(profli.xli,profli.tep,profli.tep,profli.nep,profli.nep,profli.qjli,profli.zeff,zj,profli.Raxe,profli.ftrap,profli.epsi);
    case 'Honda-Sauter'
        % ref : M. Honda et al, Nucl. Fus. 52 (2012) p 023021
        GZ = z0sauterL31(profli.xli,profli.tep,profli.nep,profli.qjli,profli.zeff,profli.Raxe,profli.ftrap,profli.epsi);
    otherwise
        
        xt = profli.ftrap ./ (1 - profli.ftrap);
        D  = 1.414 .* profli.zeff + profli.zeff .^ 2 + xt .* (0.754 + 2.657 .* profli.zeff + 2 .* profli.zeff .^ 2) + ...
            xt .^2 .* ( 0.348 + 1.243 .* profli.zeff + profli.zeff .^ 2);
        GZ  = xt .* ((0.754 + 2.21 .* profli.zeff + profli.zeff .^2) + xt .* (0.348 + 1.243 .* profli.zeff + profli.zeff .^2)) ./ D;
end
jmino    = (1 - (1 - GZ) .* zg ./ profli.zeff) .* jmino;
imino    = trapz(profli.xli,profli.spr .* jmino,2);
%figure(21);clf;plot(profli.xli,jmino);drawnow


% computation of xres
if isfield(profli,'Raxe')
    rres(~isfinite(rres)) = profli.Raxe(~isfinite(rres),1);
else
    rres(~isfinite(rres)) = geo.R(~isfinite(rres));
end
xx_loc  = linspace(-1,1,size(xa,2));
xx_fine = linspace(-1,1,2001);
xa_fine = pchip(xx_loc,xa,xx_fine);
rr_fine = pchip(xx_loc,rr,xx_fine);
dd  = (rr_fine - rres * ones(1,size(rr_fine,2))) .^ 2;
mask = (dd == (min(dd,[],2)*ones(1,size(rr_fine,2))));
xres  = min(0.95,abs(sum(xa_fine .* mask,2)) ./ geo.a);

% mise a dimension
harm = harm .* vt;

% normalisation of output
pe_int   = trapz(profli.xli,profli.vpr .* Pe',2);
pM_int   = trapz(profli.xli,profli.vpr .* PM',2);
% flip profile output
Pe = Pe' .* ((pel  ./ max(1,pe_int)) * ones(1,size(Pe,1)));
PM = PM' .* ((pion ./ max(1,pM_int)) * ones(1,size(PM,1)));
pabs = max(pabs',Pe+PM);



% adding FW contribution
pth = pth + pfw;
pel = pel + pfw;
%disp('in z0icrh_TRANG')
%  keyboard
%  pe_int   = trapz(profli.xli,profli.vpr .* Pe,2);
%  pM_int   = trapz(profli.xli,profli.vpr .* PM,2);
%  pverif   = pM_int + pe_int + pfw;
%  figure(21);
%  clf
%  plot(cons.temps,cons.picrh,'r',cons.temps,cons.picrh .* (1 - frloss),'m',cons.temps,pth,'-.b',cons.temps,pel+pion,'k.',cons.temps,pverif,'og')
%  drawnow