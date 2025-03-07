% reference  : Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER, 
% T. Eich et al, Nucl. Fusion 53 (2013) 093031 (7pp); doi:10.1088/0029-5515/53/9/093031
% plot des donnees du modele a 2 points
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
pl        = max(zs.pin ./ 100,zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz - zs.pradsol);
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
switch option.sol_model
case '2_points';
	% rien
otherwise
	warndlg('the 2 points model is not used in this simulation','2 points model');
	option.sol_model = '2_points';
end
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot,pl,zeff_div,gamma,mach_target] = ...
    z0convergence_2points_dic(option,post.z0dinput.cons,post.z0dinput.geo,post.zerod,post.profil0d);
Asol_para = 2 .* pi  .* Raxea(:,end) .* dsol .* sin(ut);
% 50 % du rayonnement du divertor retourne sur les plaques.
pref = (qpl_target  + 0.5 .* qpl_rad_div) .* Asol_para;

maskx = ones(size(zs.temps));
maskx(zs.xpoint == 0) = NaN;

% configuring the display of spread and flux expansion
spread_interval=0.2;
fl_exp_interval=0.2;


prompt1={sprintf('Scanning interval of flux spread (in mm):')};
prompt2={sprintf('Scanning interval of flux expansion from LCFS to Divertor (DIV/LCFS, starting from 1):')};
name1='Spread interval';
name2='Flux expansion interval';
numlines=2;
defaultanswer1={sprintf('%g',spread_interval)};
defaultanswer2={sprintf('%g',fl_exp_interval)};
answer1=inputdlg(prompt1,name1,numlines,defaultanswer1);
answer2=inputdlg(prompt2,name2,numlines,defaultanswer2);

if ~isempty(answer1)
    spread_interval = min(1,max(eps,str2num(answer1{1})));
end
if ~isempty(answer2)
    fl_exp_interval = min(1,max(eps,str2num(answer2{1})));
end
    
% estimation parametrique du la puissance deposee
% pour 5 angles pour les plaques de divertor : 90, 25,10, 5 et 3???
% pour 3 valeurs de l'expension du flux 1,  3 et 10,
% pour 3 valeurs de  S : 0.5,1,2 mm
% soit 6 graphes de 5 courbes
fullscreen = get(0,'ScreenSize');
h = findobj(0,'type','figure','tag','z02p_peak_power');
if isempty(h)
       h=figure('tag','z02p_peak_power');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times','toolbar','figure', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen,'name', ...
         sprintf('METIS : %s@%d / Sketch of outer divertor target power deposition (%g %% of power crossing LCFS + 1/2 of divertor radiative power, indicative data !)',post.z0dinput.machine,post.z0dinput.shot,100 .* option.fpower))
m = 1;
% the possible poloidal angles on divertor are machine-dependant
switch z0dinput.option.machine
    case 'COMPASS'
        angle    = [90 45 35 20 5];
    case 'COMPASS-U'
        angle    = [45 25 10 5 3];
    otherwise
        angle    = [90 25 10 5 3];
end
flux_exp = 1+[0 1 2 3] * fl_exp_interval;
S        = ([1 2 3 4] * spread_interval ).* 1e-3;

indx = find(post.zerod.xpoint);
if length(indx) > 100
	pas = max(1,fix(length(indx)/100));
	indx = indx(1:pas:length(indx));
elseif isempty(indx)
    indx = find(geo.d == max(geo.d));
end
for k=1:length(S)
	for l=1:length(flux_exp)
		subplot(length(S),length(flux_exp),m)
		coul      = get(gca,'colororder');
        for n=1:length(angle)
            [sb,qpdep,fnorm] = z0div_power_dep(angle(n),S(k),flux_exp(l),dsol(indx),pref(indx), ...
                abs(option.fR_target) .* geo.R(indx),max(geo.a(indx) ./ 2));
            %disp([mean(fnorm),max(fnorm),min(fnorm),median(fnorm),std(fnorm)]);
            zplotprof(gca,post.zerod.temps(indx),sb ,qpdep ./ 1e6,'color',coul(n,:));
            drawnow
        end
        if m == 1
            % made the legend description a function of the array
            % values in angle!
            leg = {};
            for lk = 1:length(angle)
                leg{lk} = sprintf('\\alpha_t = %2.0f (deg.)',angle(lk));
                legend(leg)
            end
        end
        xlabel('s (m, 0 = separatrix)')
        ylabel('MW/m^2');
        z0loglin(gca);
        title(sprintf('Spread = %g (mm) & Flux_{exp} = %g',S(k).*1e3,flux_exp(l)));
        m = m +1;
        drawnow
    end
end





