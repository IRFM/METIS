liste          = {};
liste{end+1}   = 'ppf/@shot/BOLO/TOPI';    % puissance rayonnee totale improved
liste{end+1}   = 'ppf/@shot/BOLO/TOPO';    % puissance rayonnee totale old
liste{end+1}   = 'ppf/@shot/BOLO/TOBU';    % puissance rayonnee coeur  upper
liste{end+1}   = 'ppf/@shot/BOLO/TOBH';    % puissance rayonnee coeur  horizontal
liste{end+1}   = 'ppf/@shot/BOLO/TXPN';    % puissance rayonnee xpoint and divertor   new
liste{end+1}   = 'ppf/@shot/BOLO/TOXP';    % puissance rayonnee xpoint and divertor old
bolo         = cgcgetjet(post.z0dinput.shot,liste,'','');
bolo = bolo.ppf.BOLO;

% donnees du modele a 2 points
zs      = post.zerod;
option  = post.z0dinput.option; 
geo     = post.z0dinput.geo; 
cons    = post.z0dinput.cons; 
profli  = post.profil0d;

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
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3

% compatibilite
if ~isfield(option,'yield_model')
 	option.yield_model      = 'Javev';
end

% puissance conduite a la separatrice
pl        = max(zs.pin .* sqrt(eps),zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz);
% fraction perdue en volume dans la sol par rayonnement:
%fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
switch option.sol_rad
case 'coupled'
	fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
        pradsol = zs.pradsol + max(0,1 - option.fprad) .* zs.prad;
otherwise
	fesol   = max(0,min(1, zs.pradsol ./ max(1,pl)));
        pradsol = zs.pradsol;
end

% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* option.lcx .* zs.qa) .^ 2);  
switch option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lcpol;
case 3
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lclim;
otherwise
	lc  = lcx;
end
%lc  = zs.modeh .* lcx + (~zs.modeh) .* lclim;

if isfield(zs,'dsol')
	dsol = zs.dsol;
elseif option.sol_lscale  == 0
    dsol        = geo.a ./ 100;
elseif option.sol_lscale  > 0
    dsol        = geo.a .* option.sol_lscale;
else
    dsol        = - geo.R .* option.sol_lscale;
end


% flux //  (formula 5.64 et 5.75 Stangeby)
x      = profli.xli;
ve     =  ones(size(x));
Raxea  = interp1(profli.temps,profli.Raxe,zs.temps,'pchip','extrap');
Fa     = interp1(profli.temps,profli.fdia,zs.temps,'pchip','extrap');
psi    = interp1(profli.temps,profli.psi,zs.temps,'pchip','extrap');
rmx    = interp1(profli.temps,profli.rmx,zs.temps,'pchip','extrap');
zeffp  = interp1(profli.temps,profli.zeff,zs.temps,'pchip','extrap');
nzp    = interp1(profli.temps,profli.nzp,zs.temps,'pchip','extrap');
nwp    = interp1(profli.temps,profli.nwp,zs.temps,'pchip','extrap');
nhep   = interp1(profli.temps,profli.nhep,zs.temps,'pchip','extrap');
n1p    = interp1(profli.temps,profli.n1p,zs.temps,'pchip','extrap');
nep    = interp1(profli.temps,profli.nep,zs.temps,'pchip','extrap');
Qe     = interp1(profli.temps,profli.qe(:,end),zs.temps,'pchip','extrap');
Qi     = interp1(profli.temps,profli.qi(:,end),zs.temps,'pchip','extrap');
rext         = Raxea + geo.a * x;
btor         = Fa ./ rext;
grho         = abs((rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,psi,0,2,2,1) ./ rext .* grho ./ (rmx(:,end) * ve);
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));
% flux //  (formula 5.64 et 5.75 Stangeby)
%bpola = interp1(profli.temps,profli.bpol(:,end),zs.temps,'pchip','extrap');
%Fa = interp1(profli.temps,profli.fdia(:,end),zs.temps,'pchip','extrap');
%Raxea = interp1(profli.temps,profli.Raxe(:,end),zs.temps,'pchip','extrap');

%ut = atan(bpola ./  Fa .* Raxea);
%qpl_tot     = pl  ./ (4 .* pi  .* Raxea(:,end) .* dsol .* sin(ut));
Asol_para = 4 .* pi  .* Raxea(:,end) .* dsol .* sin(ut);
switch option.sol_model
case '2_points';
	% rien
otherwise
	warndlg('the 2 points model is not used in this simulation','2 points model');
	option.sol_model = '2_points';
end
option_mem = option;
option.plot2points = 'Yes';
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot_in,pl_in,zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond] = ...
    z0convergence_2points_dic(option,post.z0dinput.cons,post.z0dinput.geo,post.zerod,post.profil0d);

prad_div_x = (qpl_rad_div + qpl_neutral_div) .* (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut)) ./ option.fpower;
prad_core  = post.zerod.prad  + post.zerod.pbrem + post.zerod.pcyclo;
prad_tot   = prad_core + prad_div_x + pradsol;

h = findobj(0,'type','figure','tag','z0plot_jet_bolo');
if isempty(h)
       h=figure('tag','z0plot_jet_bolo');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'defaultaxeslinewidth',3,'color',[1 1 1])
colormap('hot')

k = 3;
l = 1;
w = 301;
subplot(k,1,l)
hold on
leg = {};
if ~isempty(bolo.TOBH.data)
  plot(bolo.TOBH.t,bolo.TOBH.data ./ 1e6,'m');
  leg{end+1} = 'Horizontal camera';
end
if ~isempty(bolo.TOBU.data)
  plot(bolo.TOBU.t,bolo.TOBU.data ./ 1e6,'c');
  leg{end+1} = 'Upper vertical camera';
end
plot(post.zerod.temps,prad_core ./ 1e6,'k');
leg{end+1} = 'METIS';
if ~isempty(bolo.TOBH.data)
  plot(bolo.TOBH.t,sgolayfilt(bolo.TOBH.data,1,w) ./ 1e6,'r');
  leg{end+1} = 'Horizontal camera filtered';
end
if ~isempty(bolo.TOBU.data)
  plot(bolo.TOBU.t,sgolayfilt(bolo.TOBU.data,1,w) ./ 1e6,'b');
  leg{end+1} = 'Upper vertical camera filtered';
end
title(sprintf('METIS : %s@%d / JET bolometry ', ...
	  post.z0dinput.machine,post.z0dinput.shot));
ylabel('Core radiative power (MW)');
legend(leg);
z0loglin(gca);


l = l+1;
subplot(k,1,l)
hold on
leg = {};
if ~isempty(bolo.TXPN.data)
  plot(bolo.TXPN.t,bolo.TXPN.data ./ 1e6,'m');
  leg{end+1} = 'new method';
end
if ~isempty(bolo.TOXP.data)
  plot(bolo.TOXP.t,bolo.TOXP.data ./ 1e6,'c');
  leg{end+1} = 'old method';
end
plot(post.zerod.temps,prad_div_x ./ 1e6,'k');
leg{end+1} = 'METIS';
if ~isempty(bolo.TXPN.data)
 plot(bolo.TXPN.t,sgolayfilt(bolo.TXPN.data,1,w) ./ 1e6,'r');
 leg{end+1} = 'new method filtered'; 
end
if ~isempty(bolo.TOXP.data)
  plot(bolo.TOXP.t,sgolayfilt(bolo.TOXP.data,1,w) ./ 1e6,'b');
  leg{end+1} = 'old method filtered';
end
ylabel('X-point and divertor (MW)');
legend(leg);
z0loglin(gca);


l = l+1;
subplot(k,1,l)
leg = {};
hold on
if ~isempty(bolo.TOPI.data)
  plot(bolo.TOPI.t,bolo.TOPI.data ./ 1e6,'m');
  leg{end+1} = 'improved';
end
if ~isempty(bolo.TOPO.data)
  plot(bolo.TOPO.t,bolo.TOPO.data ./ 1e6,'c');
  leg{end+1} = 'old method';
end
plot(post.zerod.temps,prad_tot ./ 1e6,'k');
leg{end+1} = 'METIS';
if ~isempty(bolo.TOPI.data)
  plot(bolo.TOPI.t,sgolayfilt(bolo.TOPI.data,1,w) ./ 1e6,'r');
  leg{end+1} =  'improved filtered';
end
if ~isempty(bolo.TOPO.data)
  plot(bolo.TOPO.t,sgolayfilt(bolo.TOPO.data,1,w) ./ 1e6,'b');
  leg{end+1} = 'old method filtered';
end
ylabel('total radiative power (MW)');
legend(leg);
z0loglin(gca);
xlabel('time (s)');
joint_axes(h,k);
edition2


