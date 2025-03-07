% script to see effect of edge density at LCFS on other parameters near LCFS
% Power
Pnbi = 3e6;                   % NBI power in W
Pin_other = 0;                % other heating power in W
Prad_core = 0.5e6;              % all power lossed in the core plasma in W
fpower    = 0.6;              % fraction of power that goes to outer target
% plasma geometrie
R         = 3;                % major radius in m
a         = 1;                % minor radius in m
K         = 1.6;              % elongation
q         = 5;                % edge safety factor
ip        = 2e6;              % plasma current in A
Bt        = 3.1;              % vacuum magentic field in T

% other parameter
mode      = 1;                % 1 = divertor; 0 = poloidal limiter
Recycle   = 0.995;            % recycling coefficient
gaz_puff  = 1e21;             % core gaz-puff (so including efficiency effect) in e/s
% gestion du recyclage (default pour un divertor etanche)
Recycle_lim   = Recycle .* (1 - mode);     % direct recycling at LCFS (limiter = 1 ; divertor = 0)
Recycle_div   = Recycle .* mode;           % recycling near target (limiter = 0 ; divertor = 1)
% sion fournir les valeurs
% Recycle_lim   = 0.1;          % direct recycling at LCFS (limiter = 1 ; divertor = 0)
% Recycle_div   = 1;            % recycling near target (limiter = 0 ; divertor = 1)
%
Chi       = 0.1;              % heat transport coefficient near LCFS in m^2/s (L-mode = 10; H-mode = 0.1)
vrod      = 3;                % ratio between convective and diffusive tranport coefficient for particle near LCFS; >0 for inward pinch
dped      = 0.025 .* a;       % pedestal width or width of high gradient zone in density profile in m
Einj      = 100e3;            % NBI energie injection in eV
mass      = 2.5;              % number of mass of main ion
charge    = 1;                % number of charge of main ion 
gamma     = 7;                % sheath heat transmission coefficient
fcond     = 1;                % conduction factor
fpe       = 0.5;              % electron heat faction
fmom      = 1;                % momentum loss factor
fie       = 3;                % 1+ Tia/Tea

% control parameter
nea       = ip/1e6/pi/a^2.*1e20 .* logspace(-2,1,301)';  % edge density in m^-3

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

% start on computation
Sext = (4 * pi ^ 2) .* a .* R .* sqrt(K);       % surface of the plasma at LCFS
Peri = 2 .* pi .* a .* sqrt(K);                 % perimeter of the plasma at LCFS
%lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* R;
lcx   = sqrt(Peri .^ 2  + (pi .* R .* q) .^ 2);  
lc    = mode .* lcx + (~mode) .* lcpol;

% 2 points model
Ploss = Pnbi + Pin_other - Prad_core;           % Power crossing the LCFS
lambda    = 5671 .* real(Ploss) .^ (1/8) .* (1 + K .^ 2) .^ (5/8) .*  ...
                  a .^ (17/8) .* Bt .^ (1/4) ./ ip .^(9/8) ./ R .* ...
                  (2.* mass ./ charge);

qpar  = fpower .* Ploss ./ (2 .* pi .* a .* lambda);
vt = ones(size(nea));
qpl_corr = qpar .* vt;

figure(21):clf;hold on
plot(nea,qpar * vt ,'r');
for k = 1:31
    [tea,tet,net,indbad,noconv] = z0twopoints(qpl_corr,nea,lc .* vt ,vt ,mass,fcond .* fpe .* vt,fmom .* fie .* vt,vt,gamma .* vt);
    qpl_corr_new = min(qpar * vt, (qpar * vt)./ (1 + Recycle_div .* z0eioniz_div(tet,net) ./ gamma ./ max(eps,tet)));
    qpl_corr = 0.7 .* qpl_corr + 0.3 .* qpl_corr_new;
    plot(nea,qpl_corr)
    drawnow
end
% pedestal zone
Snbi = Pnbi ./ phys.e ./ Einj;                  % particle source due to NBI
S_R  = Recycle_lim .* nea .* Peri .* lambda .* sqrt(2 .* phys.e ./ phys.ua ./ mass .* tea);
Sn   = Snbi + S_R + gaz_puff;
nep  = (Sn ./ Sext + nea .* Chi ./ 2 ./ q ./ dped .* ( 1 + dped .* vrod ./ 2 ./ R)) ./ (  Chi ./ 2 ./ q ./ dped .* ( 1 - dped .* vrod ./ 2 ./ R));
ne   = (nea + nep)  ./ 2;
tep  = (Ploss ./ phys.e ./ Sext + tea .* (2 .* ne .* Chi ./ dped - (5/2) .* Sn ./ Sext))  ./ (2 .* ne .* Chi ./ dped + (5/2) .* Sn ./ Sext);


% script pour le plot des donnees HRTS de JEt dans les simulation avec METIS
h = findobj(0,'type','figure','tag','z0plot_jet_edge1');
if isempty(h)
       h=figure('tag','z0plot_jet_edge1');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')

k = 3;
subplot(k,1,1)
loglog(nea ./ 1e19,tet,nea ./ 1e19,tea,nea ./ 1e19,tep);
legend('T_{e,target}','T_{e,LCFS}','T_{e,Top pedestal}','location','best');
ylabel('Temperature (eV)');
%xlabel('n_{e,a} 1e19 m^{-3}');
title('effect of LCFS density on pedestal and target parameters');
z0loglin(gca);
subplot(k,1,2)
loglog(nea ./ 1e19,net./ 1e19,nea ./ 1e19,nea./ 1e19,nea ./ 1e19,nep./ 1e19);
legend('n_{e,target}','n_{e,LCFS}','n_{e,Top pedestal}','location','best');
ylabel('Density (1e19 m^{-3})');
set(gca,'ylim',[0 1e4]);
z0loglin(gca);

subplot(k,1,3)
semilogx(nea ./ 1e19,2 .* tet .* net .* phys.e ,nea ./ 1e19,2 .* tea .* nea .* phys.e ,nea ./ 1e19,2 .* tep .* nep .* phys.e );
legend('P_{target}','P_{LCFS}','P_{Top pedestal}','location','best');
ylabel('Pressure (Pa)');
xlabel('n_{e,a} 1e19 m^{-3}');
z0loglin(gca);
joint_axes(h,k);
edition2


