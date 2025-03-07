% cette fonction estime le champ electrique radial a partir des donnees du 0D
% les vitesses sont donnees pour l'espece principale
% pour les test : [t,x,er,vtheta,vtor] = z0ploter(zs,z0dinput.geo,z0dinput.cons,z0dinput.option,profli);
function [t,x,er,vtheta,vtor,wradp] = z0ploter(zs,geo,cons,option,profli,noplot)

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
    case 5
        zj = mean((1  + 4 .* real(cons.iso))   ./  (1 + 2 .* real(cons.iso)));
        aj = mean((2  + 3 .* real(cons.iso))   ./  (1 + real(cons.iso)));
        warning('nHe3onD & nTonD not yet implemented !');

    case 11
        zj = mean((1  + 25 .* cons.iso)   ./ (1 + 5.* cons.iso));
        aj = mean((1  + 11 .* cons.iso)  ./  (1 + cons.iso));
end

% impurete principale
zimp = option.zimp;
%aimp = ceil(zimp.* (7/3));

% 2ieme impurete
zmax = option.zmax;
%amax = ceil(zmax.* (7/3));

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - zimp);
mask = (dd == min(dd));
aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(aimp)
    aimp = 7/3 .* zimp;
end
dd   = abs(Z_el - zmax);
mask = (dd == min(dd));
amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(amax)
    amax = 7/3 .* zmax;
end


t = profli.temps;
x = profli.xli;
ve = ones(size(x));
vt = ones(size(t));
a = interp1(zs.temps,geo.a,t,'linear');
rext         = (profli.Raxe + a * x);

btor         = option.signe .* (profli.fdia./ rext);
grho         = abs((profli.rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profli.psi,0,2,2,1)./ rext .* grho ./ (profli.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);
erp          = grho .* pdederive(x,phys.e .* profli.tip .* profli.nzp,0,2,2,1) ./ ...
                     profli.nzp ./ phys.e ./ zimp ./ (profli.rmx(:,end) * ve);

wradp = profli.omega;		     
		     
% le plot
hz =findobj(0,'type','figure','tag','0der');
if isempty(hz)
  	  hz=figure('tag','0der','name','Radial eletric field');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1);
zplotprof(gca,t*ve,rext,btor,'color','r');
zplotprof(gca,t*ve,rext,bpol,'color','b');
ylabel('Btor (r) & Bpol (b)  ( in T)');
xlabel('Rext (m)')
title('Bphi  & Ip in the same direction');
subplot(2,2,2);
zplotprof(gca,t*ve,rext,profli.vtheta,'color','b');
zplotprof(gca,t*ve,rext,profli.vtor,'color','r');
zplotprof(gca,t*ve,rext,profli.omega .* rext,'color','g');
ylabel('Vtheta (b)  , Vtor (r) & bulk rotation (g) (in m/s)');
xlabel('Rext (m)')
title('Rotation velocity of main impurity');
subplot(2,2,3);
zplotprof(gca,t*ve,rext,erp,'color','c');
zplotprof(gca,t*ve,rext,- profli.vtheta .* btor,'color','b');
zplotprof(gca,t*ve,rext,- profli.vtheta .* btor + erp,'color','m');
zplotprof(gca,t*ve,rext, profli.vtor .* bpol,'color','r');
ylabel('component : pressure (c), vtheta (b), pressure + vtheta (m) & vtor (r)  (in V/m)');
xlabel('Rext (m)')

subplot(2,2,4);
zplotprof(gca,t*ve,rext,profli.er .* grho,'color','r');
ylabel('Er (in V/m)');
xlabel('Rext (m)')
