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

if ~isfield(option,'Sn_fraction')
 	option.Sn_fraction      = 0;
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

% recycling control
if option.Recycling_target == 0
      % use global recycling coeficient
      Recycling = option.Recycling;
else
      % use taget recycling coeficient
      Recycling = option.Recycling_target;
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

if isfield(option,'Sq') && (option.Sq > 0)
  f_Sq   = 1.64 .* option.Sq ./ (zs.dsol + 1.64 .* option.Sq);
  fpower = option.fpower.* f_Sq;
elseif isfield(option,'Sq') && (option.Sq < 0)
  %A. Scarabosio paper (Journal of Nuclear Materials 463 (2015) 49-54)
  Sq_mm = abs(option.Sq) .* (double(zs.modeh) .* (0.12 .* (zs.plim ./ 1e6) .^ 0.21 .*  max(1e13/1e19,zs.nebord ./ 1e19) .^ -0.02) .* interp1(profli.temps,profli.bpol(:,end),cons.temps,'linear','extrap') .^ -0.82 .* geo.R .^ 0.71 +  ...
	  (1- double(zs.modeh)) .* (0.13 .* max(1e13/1e19,zs.nebord ./ 1e19) .^ 1.1) .* interp1(profli.temps,profli.bpol(:,end),cons.temps,'linear','extrap') .^ -0.84);
  %figure(21);plot(cons.temps,Sq_mm);drawnow
  f_Sq   = 1.64 .* Sq_mm ./ (zs.dsol + 1.64 .* Sq_mm);
  fpower = option.fpower .* f_Sq;
else
  fpower = option.fpower; 
end



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
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot_in,pl_in,zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond,profli] = ...
    z0convergence_2points_dic(option,post.z0dinput.cons,post.z0dinput.geo,post.zerod,post.profil0d);
option = option_mem;
qpl_rad_sol = fpower .* pradsol ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_tot   = fpower .* pl ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_in_max = fpower .* zs.pin ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
%(4*pi*rp*sol_width/q95)
%qpar=p_sep*1e6/(4*pi*rp*sol_width/q95); % W/m2 estimate of parallel q
%L=pi*q95*rp;
qpl_target  = qpl_tot -  qpl_rad_sol - qpl_rad_div - qpl_neutral_div;

fie = 1 + zs.tibord ./ zs.tebord;
tite_loc = zs.tibord ./ zs.tebord;
%fpe = min(1,max(0.1,zs.pel ./ (zs.pel + zs.pion)));
fpe = min(1,max(0.1,Qe ./ (Qe + Qi)));



if option.fmom == 0
	% longueur de ionisation pres des plaques
	[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(zs.telim,zs.telim .* (zs.tibord ./ zs.tebord));  % attention ici le rapport ti/te est calculer a l'exterieur, tebord ne doit pas etre mis a jour
	%% equilibre entre 1s et 2s pour le neutres de centre 
	alphas = nelim .* sss ./ (nelim .* sss + Ass);
	% etat d'equilibre 1s/2s
	sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
	alpha  = (sviss  + sii)  ./ (sviss  + sii+ svcx);
	%flimx  = min(1,exp(telim - option.eioniz));
	fmom_corr    = max(0.1,min(1,0.17 + 2 .* (alpha ./ (alpha + 1)) .^ ((alpha + 1) ./ 2)));
else
	fmom_corr = option.fmom .* ones(size(zs.telim));		
end


% flux de matiere dans la sol
switch option.gaz
    case 1
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1);
    case 11
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
                            (1 + 11 .* cons.iso) .* (1 + cons.iso));
    case 2
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2);
    case {3,5}
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
                             (2 + 3.* cons.iso) .* (1 + cons.iso));
    case 4
        cs0          = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4);
end
flux_target =  mach_target .* cs0 .* sqrt(zs.telim) .* zs.nelim;
neu         =  (4 .* zs.telim .* zs.nelim) ./ (fmom .* fie) ./ zs.tebord;

if option.mach_corr == 1
	ftm = (gamma .* (1 + mach_target .^ 2) ./ (2 .* abs(mach_target) .* (gamma - 1 + mach_target .^ 2))) .^ 2;
	fnm = 2 ./ (ftm .* (1 + mach_target .^ 2));
	neu = neu ./ ftm ./ fnm;
        flux_target = flux_target .* mach_target;
end
flux_output = interp1(post.profil0d.temps,post.profil0d.ge(:,end) .* post.profil0d.grho2(:,end) .* post.profil0d.vpr_tor(:,end), ...
			   post.zerod.temps,'pchip','extrap');

maskx = ones(size(zs.temps));
maskx(zs.xpoint == 0) = NaN;

h = findobj(0,'type','figure','tag','z02p');
if isempty(h)
       h=figure('tag','z02p');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1],'toolbar','figure')
k    = 6;
subplot(k,1,1)
plot(zs.temps,lc,zs.temps,zs.sext ./ 100,zs.temps,10 .* Raxea(:,end),zs.temps,1000 .* dsol,zs.temps,10 .* zs.sext .* dsol,zs.temps,400 .* pi  .* Raxea(:,end) .* dsol .* sin(ut));
xlabel('time (s)');
legend('Lc (m)','Area (10^2 m^2)','R_0 * 10 (m)','Width (mm)','V_{SOL} * 10 (m^3)','Asol * 100 (m^2)','Location','north','orientation','horizontal')
switch post.z0dinput.option.sol_model
case '2_points'
	title(sprintf('METIS : %s@%d/ Two Point model', ...
          post.z0dinput.machine,post.z0dinput.shot));
otherwise
	title(sprintf('METIS : %s@%d/ Edge scaling', ...
          post.z0dinput.machine,post.z0dinput.shot));
end
subplot(k,1,2)
plot(zs.temps,sin(ut),zs.temps,fie,zs.temps,fpe,zs.temps,fcond ./ fpe ,zs.temps,zs.zeff,zs.temps,zeff_div,zs.temps,zs.tauhe ./ max(eps,zs.taue),zs.temps,mach_target,zs.temps,fmom,'-.');
%xlabel('time (s)');
z0loglin(gca);
legend('sin(atan(Bp/Bt))','F_{ie}','F_{pe}','F_{cond}','Zeff_{core}','Zeff_{div}','tau_{he}/tau_{E}','Mach_T','f_{mom}','Location','north','orientation','horizontal')
xx=[zs.temps(:);zs.temps(end:-1:1)];
yy=[max(fmom,fmom_corr);min(fmom(end:-1:1),fmom_corr(end:-1:1))];
hl=patch(xx,yy,[0.7,0.7,0.7],'erasemode','xor','edgecolor','none');

subplot(k,1,3)
plot(zs.temps,maskx .* qpl_tot./1e6,'r',zs.temps,maskx .* qpl_rad_sol./1e6,'b',zs.temps,maskx .* qpl_rad_div./1e6,'c',zs.temps,maskx .* qpl_neutral_div./1e6,'m',zs.temps,maskx .* qpl_target./1e6,'k');
%xlabel('time (s)');
ylabel('MW/m^2');
legend('qpl_{total}','qpl_{rad}^{SOL}','qpl_{rad}^{DIV}','qpl_{neutral}','qpl_{DIV}','Location','north','orientation','horizontal')

subplot(k,1,4)
semilogy(zs.temps,maskx .* flux_target./1e22,zs.temps,zs.n0a./1e22,zs.temps,flux_output./1e22);
xpoint = zs.xpoint;
xpoint(xpoint < 1) = NaN;
hold on
plot(zs.temps,xpoint,'og')
z0loglin(gca);
%xlabel('time (s)');
legend('Target flux (10^{22} e/s)','Recycling LCFS (10^{22} e/s)','Output plasma flux (10^{22} e/s)','Xpoint','Location','north','orientation','horizontal')

subplot(k,1,5)
plot(zs.temps,zs.telim,'r',zs.temps,zs.tebord,'b');
%xlabel('time (s)');
legend('T_{e,plate} (eV)','T_{e,LCFS} (eV)','Location','north','orientation','horizontal')
xx=[zs.temps(:);zs.temps(end:-1:1)];
yy=[max(zs.telim,telim);min(zs.telim(end:-1:1),telim(end:-1:1))];
hl=patch(xx,yy,[1,0.7,0.7],'erasemode','xor','edgecolor','none');
yy=[max(zs.tebord,tebord);min(zs.tebord(end:-1:1),tebord(end:-1:1))];
hl=patch(xx,yy,[0.7,0.7,1],'erasemode','xor','edgecolor','none');

subplot(k,1,6)
semilogy(zs.temps,zs.nelim/1e19,'r',zs.temps,zs.nebord/1e19,'b',zs.temps,maskx .* neu ./ 1e19,':g');
z0loglin(gca);
xlabel('time (s)');
legend('n_{e,plate} (10^{19} m^{-3})','n_{e,LCFS} (10^{19} m^{-3})','n_{e,u} (10^{19} m^{-3})','Location','north','orientation','horizontal')
xx=[zs.temps(:);zs.temps(end:-1:1)];
yy=[max(zs.nelim,nelim);min(zs.nelim(end:-1:1),nelim(end:-1:1))] ./ 1e19;
hl=patch(xx,yy,[1,0.7,0.7],'erasemode','xor','edgecolor','none');

joint_axes(h,k);



% fuite des impuretes depuis le divertor
k = 0;	
ya0_imp = [];
ya0_max = [];
ya0_w = [];

switch option.yield_model
case 'Matsunami'
	modely = 0;
otherwise
	modely = 1;
end

try
    
    % carbon feed back on plate
    if isfield(option,'carbonblow') && (option.carbonblow ~= 0) && (option.zimp == 6)
        [ya0_imp,zavez_C]    = z0carbonblow(option,zs,profli);
    else
        czimp = nzp(:,end) ./ nep(:,end);
        % concentration de l'imurete zmax
        czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
        % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
        cws =  nwp(:,end) ./ nep(:,end);
        % concentration en He
        switch option.gaz
            case 5
                 che  =  option.frhe0 .* profil0d.nep(:,end);
                 che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            otherwise
                 che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
                 che3 =  zeros(size(profil0d.nep(:,end)));
        end
        % le reste pour HDT
        if option.Sn_fraction > 0
            zwsn   = (1 - option.Sn_fraction) .* z0wavez(zs.tebord) + option.Sn_fraction .* z0snavez(zs.tebord);
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn - 2 .* che - 2 .* che3);
        else
            creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zs.tebord) - 2 .* che - 2 .* che3);
        end
        % H D T
        switch option.gaz
            case 11
                cB  = zerod.nTm ./zerod.nem;
                creste = max(0,crest - 5 .* cB);
                cT  = zeros(size(zerod.n1m));
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
            otherwise
                cT  = creste .* zerod.nTm ./zerod.n1m;
                cD  = creste .* zerod.nDm ./zerod.n1m;
                cB  = zeros(size(zerod.n1m));
                cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
        end
        
        % le yield effectif est la somme des yields ponderee :
        % il y une seule temperature sur le limiteur
        ya0_imp  = czimp .* z0yield2(zs.telim,zs.telim,option.zimp,option.zimp,modely) + czmax .* z0yield2(zs.telim,zs.telim,option.zmax,option.zimp,modely) + ...
            cws   .* z0yield2(zs.telim,zs.telim,'W',option.zimp,modely) + che .* z0yield2(zs.telim,zs.telim,'He',option.zimp,modely) + ...
            cH    .* z0yield2(zs.telim,zs.telim,'H',option.zimp,modely) + cD  .* z0yield2(zs.telim,zs.telim,'D',option.zimp,modely)  + ...
            cT    .* z0yield2(zs.telim,zs.telim,'T',option.zimp,modely) + che3 .* z0yield2(zerod.telim,zerod.telim,'He3',option.zimp,modely) + ...
            cB    .* z0yield2(zerod.telim,zerod.telim,'B',option.zimp,modely);
    end
    if all(isfinite(ya0_imp))
        k       = k + 1;
        [void,Lvoid,leg_imp]  = z0yield2(1000,1000,option.zimp,option.zimp,1);
    else
        ya0_imp = [];
        leg_imp = '';
    end
catch
    ya0_imp = [];
    leg_imp = '';
end

try
    czimp = nzp(:,end) ./ nep(:,end);
    % concentration de l'imurete zmax
    czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
    % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
    cws =  nwp(:,end) ./nep(:,end);
    % concentration en He
    switch option.gaz
        case 5
            che  =  option.frhe0 .* profil0d.nep(:,end);
            che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
        otherwise
            che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
            che3 =  zeros(size(profil0d.nep(:,end)));
    end
    % le reste pour HDT
    if option.Sn_fraction > 0
        zwsn   = (1 - option.Sn_fraction) .* z0wavez(zs.tebord) + option.Sn_fraction .* z0snavez(zs.tebord);
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn - 2 .* che - 2 .* che3);
    else
        
        creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zs.tebord) - 2 .* che - 2 .* che3);
    end
    % H D T
    switch option.gaz
        case 11
            cB  = zerod.nTm ./zerod.nem;
            creste = max(0,crest - 5 .* cB);
            cT  = zeros(size(zerod.n1m));
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
        otherwise
            cT  = creste .* zerod.nTm ./zerod.n1m;
            cD  = creste .* zerod.nDm ./zerod.n1m;
            cB  = zeros(size(zerod.n1m));
            cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
    end
    
    % le yield effectif est la somme des yields ponderee :
    % il y une seule temperature sur le limiteur
    ya0_max  = czimp .* z0yield2(zs.telim,zs.telim,option.zimp,option.zmax,modely) + czmax .* z0yield2(zs.telim,zs.telim,option.zmax,option.zmax,modely) + ...
        cws   .* z0yield2(zs.telim,zs.telim,'W',option.zmax,modely) + che .* z0yield2(zs.telim,zs.telim,'He',option.zmax,modely) + ...
        cH    .* z0yield2(zs.telim,zs.telim,'H',option.zmax,modely) + cD  .* z0yield2(zs.telim,zs.telim,'D',option.zmax,modely)  + ...
        cT    .* z0yield2(zs.telim,zs.telim,'T',option.zmax,modely)+ che3 .* z0yield2(zerod.telim,zerod.telim,'He3',option.zmax,modely) + ...
        cB    .* z0yield2(zerod.telim,zerod.telim,'B',option.zmax,modely);
    if all(isfinite(ya0_max))
        k       = k + 1;
        [void,Lvoid,leg_max]  = z0yield2(1000,1000,option.zmax,option.zmax,1);
    else
        ya0_max = [];
        leg_max = '';
    end
catch
    ya0_max = [];
    leg_max = '';
end

switch option.gaz
    case 11
        try
            czimp = nzp(:,end) ./ nep(:,end);
            % concentration de l'imurete zmax
            czmax = 0.01 .* abs(option.fzmax_div) + option.rimp .* czimp;
            % la concentration pour le self sputtering : prise en compte de ce qui n'est pas redepose sur place.
            cws =  nwp(:,end) ./nep(:,end);
            % concentration en He
            switch option.gaz
                case 5
                    che  =  option.frhe0 .* profil0d.nep(:,end);
                    che3 =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
                otherwise
                    che  =  profil0d.nhep(:,end) ./ profil0d.nep(:,end);
                    che3 =  zeros(size(profil0d.nep(:,end)));
            end
            % le reste pour HDT
            if option.Sn_fraction > 0
                zwsn   = (1 - option.Sn_fraction) .* z0wavez(zs.tebord) + option.Sn_fraction .* z0snavez(zs.tebord);
                creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* zwsn - 2 .* che - 2 .* che3);
            else
                
                creste = max(0,1 - czimp .* option.zimp  - czmax .*  option.zmax - cws .* z0wavez(zs.tebord) - 2 .* che - 2 .* che3);
            end
            % H D T
            switch option.gaz
                case 11
                    cB  = zerod.nTm ./zerod.nem;
                    creste = max(0,crest - 5 .* cB);
                    cT  = zeros(size(zerod.n1m));
                    cD  = creste .* zerod.nDm ./zerod.n1m;
                    cH  = creste .* max(0,zerod.n1m - zerod.nDm) ./zerod.n1m;
                otherwise
                    cT  = creste .* zerod.nTm ./zerod.n1m;
                    cD  = creste .* zerod.nDm ./zerod.n1m;
                    cB  = zeros(size(zerod.n1m));
                    cH  = creste .* max(0,zerod.n1m - zerod.nDm - zerod.nTm) ./zerod.n1m;
            end
            
            % le yield effectif est la somme des yields ponderee :
            % il y une seule temperature sur le limiteur
            ya0_B  = czimp .* z0yield2(zs.telim,zs.telim,option.zimp,'B',modely) + czmax .* z0yield2(zs.telim,zs.telim,5,'B',modely) + ...
                cws   .* z0yield2(zs.telim,zs.telim,'W','B',modely) + che .* z0yield2(zs.telim,zs.telim,'He','B',modely) + ...
                cH    .* z0yield2(zs.telim,zs.telim,'H','B',modely) + cD  .* z0yield2(zs.telim,zs.telim,'D','B',modely)  + ...
                cT    .* z0yield2(zs.telim,zs.telim,'T','B',modely)+ che3 .* z0yield2(zerod.telim,zerod.telim,'He3','B',modely) + ...
                cB    .* z0yield2(zerod.telim,zerod.telim,'B','B',modely);
            if all(isfinite(ya0_B))
                k       = k + 1;
                [void,Lvoid,leg_B]  = z0yield2(1000,1000,'B','B',1);
            else
                ya0_B = [];
                leg_B = '';
            end
        catch
            ya0_B = [];
            leg_B = '';
        end
    otherwise
        ya0_B = [];
        leg_B = '';
end

% calcul des donnees relatives a W
[void_nwp,void_nwm,void_zu1w,void_zu2w,tleak_w,fw,fwth_fit,ya0_w,cw_edge,fraction_prompt,zwavep] = ...
    z0acctungsten(option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.profil0d);
zwavep          = interp1(post.profil0d.temps,zwavep,post.z0dinput.cons.temps,'linear','extrap');
fraction_prompt = interp1(post.profil0d.temps,fraction_prompt,post.z0dinput.cons.temps,'linear','extrap');
cw_edge         = interp1(post.profil0d.temps,cw_edge,post.z0dinput.cons.temps,'linear','extrap');
ya0_w           = interp1(post.profil0d.temps,ya0_w,post.z0dinput.cons.temps,'linear','extrap');
fwth_fit        = interp1(post.profil0d.temps,fwth_fit,post.z0dinput.cons.temps,'linear','extrap');
fw              = interp1(post.profil0d.temps,fw,post.z0dinput.cons.temps,'linear','extrap');
tleak_w         = interp1(post.profil0d.temps,tleak_w,post.z0dinput.cons.temps,'linear','extrap');
if option.W_effect > 0
    
    if ~isempty(ya0_w)
        k       = k + 1;
    end
    
else
    ya0_w = [];
end


% estimation de delta_s
% Stangeby 4.55 :
[svi1s,svi2s,svcx,svrec,sii,sss,Ass] = z0sectionh(post.zerod.telim,post.zerod.telim .* post.zerod.tibord ./ post.zerod.tebord);

%% equilibre entre 1s et 2s pour le neutres de centre 
alphas = post.zerod.nelim .* sss ./ ( post.zerod.nelim .* sss + Ass);
% etat d'equilibre 1s/2s
sviss  = svi2s .* alphas + svi1s .* (1 - alphas);
s3  = sviss  + sii + svcx;
delta_s = sqrt(2 .* phys.e .* post.zerod.telim  ./ post.zerod.meff ./ phys.mp) ./ post.zerod.nelim ./  s3 ./ ut;

% temperature de fuite
% stangeby 6.87 p 317
%tleak_w       = 3.8e-9 .* (abs(post.z0dinput.option.ftwleak) .* z0wavez(post.zerod.tebord) + (1 - abs(post.z0dinput.option.ftwleak)) .* z0wavez(post.zerod.telim)) .* sqrt(post.zerod.nelim .* delta_s);
%fw         = max(1e-308,exp(- (tleak_w ./ post.zerod.telim) .^ 2));
%  if option.ftwleak < 0
%  	% from ref :DIVIMP simulation ..., A. Jarvinen et al, Physica Sripta T145 (2011) p 014013-
%  	% etalonage de la fuite
%  	fw_fit    = [10 1 0.3253 1.5101e-16 1e-45];
%  	fw_leak   = [1 1 0.035  4.0000e-06 1e-308];
%  	fwth_fit  = min(1,exp(interp1(log(fw_fit),log(fw_leak),log(fw),'pchip','extrap')));
%  	fwth_fit(~isfinite(fwth_fit)) = 0;
%  	fwth_fit(fw < 1e-45) = 0;
%  else
%  	fwth_fit = fw;
%  end
tleak_zmax = 3.8e-9 .* post.z0dinput.option.zmax .* sqrt(post.zerod.nelim .* delta_s);
tleak_zimp = 3.8e-9 .* post.z0dinput.option.zimp .* sqrt(post.zerod.nelim .* delta_s);
tleak_B = 3.8e-9 .* 5 .* sqrt(post.zerod.nelim .* delta_s);



%  % model de redepostion 
%  % ref : Stangeby NF 51 (2011) 063001
%  % ref D. Naujoks et al, NF 36 (1996)  p671-687
%  % la section efficace de ionisation de W est complexe  a calculer avec ions + electrons + effet cinetique
%  tp   = [0 	1 	10 	100 	1000 	1e5];
%  svp  = [1e-10	1e-8	1.5e-7  5e-7	2e-7	1e-7];
%  svwi = pchip(tp,svp,zs.telim) .* 1e-6; % m^-3  * s^-1
%  % self sputtering dominant ou par impurete legere : E0 = telim
%  lioniz   =  sqrt(2 .* zs.telim .* phys.e ./ (183.84 .* phys.ua)) ./ zs.nelim ./  svwi;
%  % 1 fois ionise
%  larmor_w =  4.57e-3 .* sqrt(183.84) .* sqrt(zs.telim ./ 1e3) ./ geo.b0;
%  p        = lioniz ./ larmor_w;
%  % from ERO with all effect
%  pp       = [0 	0.2 	0.6 	2.2 	3.2 	6	1000];
%  fp       = [1	0.85	0.65	0.3	0.25	0.15	0.1];   
%  fraction_prompt =  min(1-eps,max(0.1,pchip(pp,fp,max(0,min(p,1000)))));



h = findobj(0,'type','figure','tag','z0leak');
if isempty(h)
    h=figure('tag','z0leak');
else
    figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
    'defaultlinelinewidth',1,'color',[1 1 1],'toolbar','figure')
l = 1;
subplot(k+2,1,l)
% plot de Tleak
semilogy(post.zerod.temps,tleak_w,'b',post.zerod.temps,tleak_zmax,'c',post.zerod.temps,tleak_zimp,'m',post.zerod.temps,tleak_B,'g',post.zerod.temps,post.zerod.telim,'r')
title(sprintf('METIS : %s@%d / Divertor leak ', ...
    post.z0dinput.machine,post.z0dinput.shot));
ylabel('eV');
legend('T_e^{W,leak}',sprintf('T_e^{Z=%d,leak}',post.z0dinput.option.zmax),sprintf('T_e^{Z=%d,leak}',post.z0dinput.option.zimp),'T_e^{B,leak}','T_e^{target}');
z0loglin(gca);

l = l+1;
subplot(k+2,1,l)
semilogy(post.zerod.temps,fwth_fit,'b',post.zerod.temps,exp(- (tleak_zmax ./ post.zerod.telim) .^ 2),'c', ...
    post.zerod.temps,exp(- (tleak_zimp ./ post.zerod.telim) .^ 2),'m',post.zerod.temps,exp(- (tleak_B ./ post.zerod.telim) .^ 2),'g')
set(gca,'ylim',[1e-6,1]);
ylabel('leakage factor');
legend('f_{W,leak}',sprintf('f_{Z=%d,leak}',post.z0dinput.option.zmax),sprintf('f_{Z=%d,leak}',post.z0dinput.option.zimp),'f_{B,leak}');
z0loglin(gca);

l = l+1;

% pour les limites : REF Stangeby NF 51 (2011) p 063001
if ~isempty(ya0_imp)
    if k > 0
        subplot(k+2,1,l)
    end
    if option.zimp  ~= 74
        semilogy(post.zerod.temps,ya0_imp,'b')
        z0loglin(gca);
        ylabel('Sputtering yield in divertor');
        legend(leg_imp)
    else
        semilogy(post.zerod.temps,ya0_imp,'b', ...
            post.zerod.temps,0.00045 .* ones(size(post.zerod.temps)),'r',  ...
            post.zerod.temps,1e-8 .* ones(size(post.zerod.temps)),'g')
        z0loglin(gca);
        ylabel('Sputtering yield in divertor');
        legend(leg_imp,'ITER design','Reactor limit')
    end
    set(gca,'ylim',[1e-10,Inf]);
    l = l+1;
end
if ~isempty(ya0_max)
    if k > 1
        subplot(k+2,1,l)
    end
    if (option.zmax  ~= 74)
        semilogy(post.zerod.temps,ya0_max,'b')
        z0loglin(gca);
        ylabel('Sputtering yield in divertor');
        legend(leg_max)
    else
        semilogy(post.zerod.temps,ya0_max,'b', ...
            post.zerod.temps,0.00045 .* ones(size(post.zerod.temps)),'r',  ...
            post.zerod.temps,1e-8 .* ones(size(post.zerod.temps)),'g');
        z0loglin(gca);
        ylabel('Sputtering yield in divertor');
        legend(leg_max,'ITER design','Reactor limit')
    end
    set(gca,'ylim',[1e-10,Inf]);
    l = l+1;
    
end
if ~isempty(ya0_B)
    if k > 1
        subplot(k+2,1,l)
    end
    semilogy(post.zerod.temps,ya0_B,'b', ...
        post.zerod.temps,0.00045 .* ones(size(post.zerod.temps)),'r',  ...
        post.zerod.temps,1e-8 .* ones(size(post.zerod.temps)),'g');
    z0loglin(gca);
    ylabel('Sputtering yield in divertor');
    legend(leg_B,'ITER design','Reactor limit')
    set(gca,'ylim',[1e-10,Inf]);
    l = l+1;
    
end
if option.W_effect > 0
    if k > 1
        subplot(k+2,1,l)
    end
    if ~isempty(ya0_w)
        semilogy(post.zerod.temps,ya0_w,'b',post.zerod.temps,ya0_w .* (1 - fraction_prompt),'c', ...
            post.zerod.temps,0.00045 .* ones(size(post.zerod.temps)),'r',  ...
            post.zerod.temps,1e-8 .* ones(size(post.zerod.temps)),'g');
        z0loglin(gca);
        ylabel('Sputtering yield in divertor');
        legend('W or Sn','W or Sn with prompt','ITER design','Reactor limit')
        set(gca,'ylim',[1e-10,Inf]);
        
    end
end
xlabel('time (s)');
if k > 1
    joint_axes(h,k);
end

%si W
if option.W_effect > 0
    
    zwave = trapz(post.profil0d.xli,z0wavez(post.profil0d.tep) .* post.profil0d.nep .* post.profil0d.vpr,2) ./  ...
        trapz(post.profil0d.xli,post.profil0d.nep .* post.profil0d.vpr,2);
    
    if option.Sn_fraction > 0
        zwave = (1 - option.Sn_fraction) .*  zwave + option.Sn_fraction .* trapz(post.profil0d.xli,z0snavez(post.profil0d.tep) .* post.profil0d.nep .* post.profil0d.vpr,2) ./  ...
            trapz(post.profil0d.xli,post.profil0d.nep .* post.profil0d.vpr,2);
    end
    
    [nwp,nwm,zu1w,zu2w] = z0acctungsten(option,post.z0dinput.geo,post.z0dinput.cons,post.zerod,post.profil0d);
    
    zu1 = (option.zimp + option.rimp .* option.zmax + zu1w);
    zu2 = (option.zimp .^ 2 + option.rimp .* option.zmax .^ 2 + zu2w);
    
    if option.Sn_fraction > 0
        zwloc  = (1 - option.Sn_fraction) .* z0wavez(post.profil0d.tep) + option.Sn_fraction .* z0snavez(post.profil0d.tep);
    else
        zwloc  = z0wavez(post.profil0d.tep);
    end
    nzwloc = nwp .* zwloc;
    dnep_w = nzwloc - nzwloc(:,end) * ve;
    
    
    fullscreen = get(0,'ScreenSize');
    h = findobj(0,'type','figure','tag','z0W');
    if isempty(h)
        h=figure('tag','z0W');
    else
        figure(h);
    end
    clf
    set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',1,'color',[1 1 1],'Position',fullscreen,'toolbar','figure')
    
    subplot(3,1,1)
    semilogy(post.profil0d.temps,zu1,'b',post.profil0d.temps,zu2,'r')
    z0loglin(gca);
    xlabel('time (s)');
    legend('Volume averaged Z  ','Volume averaged Z^2');
    title(sprintf('METIS : %s@%d / Tungsten effects', ...
        post.z0dinput.machine,post.z0dinput.shot));
    
    subplot(3,3,4)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.zeff,'color','b');
    xlabel('r/a');
    ylabel('Zeff');
    
    subplot(3,3,5)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.nzp ./ post.profil0d.nip,'color','b');
    xlabel('r/a');
    ylabel('n_{imp} / n_{i}');
    
    
    subplot(3,3,6)
    if option.Sn_fraction > 0
        zplotprof(gca,post.profil0d.temps,post.profil0d.xli,z0wavez(post.profil0d.tep),'color','b');
        hold on
        zplotprof(gca,post.profil0d.temps,post.profil0d.xli,z0snavez(post.profil0d.tep),'color','r');
        xlabel('r/a');
        ylabel('<Z_W> in b & <Z_{Sn}> in r');
    else
        zplotprof(gca,post.profil0d.temps,post.profil0d.xli,z0wavez(post.profil0d.tep),'color','b');
        xlabel('r/a');
        ylabel('<Z_W>');
    end
    
    subplot(3,3,7)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,dnep_w ./ post.profil0d.nep,'color','b');
    xlabel('r/a');
    ylabel('dn_{e,W} / n_{e}');
    
    subplot(3,3,8)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.prad,'color','b');
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.pbrem,'color','r');
    xlabel('r/a');
    ylabel('W.m^{-3}');
    legend('P_{line}','P_{brem}');
    
    subplot(3,3,9)
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,nwp,'color','b');
    zplotprof(gca,post.profil0d.temps,post.profil0d.xli,post.profil0d.nzp,'color','r');
    set(gca,'yscale','log');
    z0loglin(gca);
    xlabel('r/a');
    ylabel('m^{-3}');
    legend('n_{W}','n_{imp}');
    
    
    
end



