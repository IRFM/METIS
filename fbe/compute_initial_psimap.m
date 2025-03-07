% fonction to initialise poloidal flux map for FBE code : plasma part ( the coil part must be added)
% itime = index of time 
% post  = METIS data structure
% (Rfbe,Zfbe) = coordinate of the points,
% equi_extrap = if = 0, exptrapolation in vaccum use simple interpolation;
% if = 1, use analytical solution of G-S equation in vacuum; if = 2,force  to use analytical solution of G-S equation in vacuum even if not % converged
% 
function [psi2d,BR2d,BZ2d,BPHI2d]= compute_initial_psimap(time,post,Rfbe,Zfbe,equi_extrap)

%% PHYSICS CONSTANTS
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


% normalisation 
factor_two_pi = 1;  
if nargin < 5
  equi_extrap = 0;
end

%% CONSTANT = NUMBER OF POLOIDAL POINTS
if isfield(post.z0dinput.option,'nb_points_pol')
    nbth = post.z0dinput.option.nb_points_pol;
else
    nbth = 65;
end  
rzp = max(1,max(post.z0dinput.geo.K(:)));
% decode
profil0d = post.profil0d;
itime = find(profil0d.temps >= time,1);
if isempty(itime)
  itime = length(profil0d.temps);
end

xpoint = interp1_imas(post.zerod.temps,post.zerod.xpoint,profil0d.temps,'nearest','extrap');
boundary_type = xpoint(itime);   
z0 = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.z0,profil0d.temps(itime),'nearest','extrap');
if isfield(profil0d,'Rsepa') && isfield(profil0d,'Zsepa') && all(isfinite(profil0d.Rsepa(itime,:))) && all(isfinite(profil0d.Zsepa(itime,:)))
     % deduire de la sepa
    control =sepa_moments(profil0d.Rsepa(itime,:),profil0d.Zsepa(itime,:) + z0);
    z0_mag = control.z0;
    outlineisfilled = 1;
else
    control = [];
    z0_mag = z0;
    outlineisfilled = 0;
end

%% METIS DATA
prof.x    = profil0d.xli;
prof.kx   = profil0d.kx(itime,:);     
prof.dx   = profil0d.dx(itime,:);      
prof.rho  = profil0d.rmx(itime,:);     
prof.Raxe = profil0d.Raxe(itime,:);
prof.epsi = profil0d.epsi(itime,:);
prof.psi  = profil0d.psi(itime,:);
prof.phi  = profil0d.phi(itime,:);
prof.dphidx = pdederive(prof.x,prof.phi,0,2,2,1);
prof.q    = profil0d.qjli(itime,:);
prof.fdia = profil0d.fdia(itime,:);
prof.jmoy = profil0d.jli(itime,:);
prof.ptot = profil0d.ptot(itime,:);
prof.vpr  = profil0d.vpr_tor(itime,:);
prof.dvdx = profil0d.vpr(itime,:);
prof.spr  = profil0d.spr(itime,:) ./ pdederive(prof.x,prof.rho,1,2,2,1);
prof.dsdx = profil0d.spr(itime,:);
prof.volume  = cumtrapz(profil0d.xli,profil0d.vpr(itime,:),2);
prof.surface = cumtrapz(profil0d.xli,profil0d.spr(itime,:),2);
%
geo_input.a     = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.a ,profil0d.temps(itime),'pchip','extrap');
geo_input.R     = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.R ,profil0d.temps(itime),'pchip','extrap');
geo_input.K     = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.K ,profil0d.temps(itime),'pchip','extrap');
geo_input.d     = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.d ,profil0d.temps(itime),'pchip','extrap');
geo_input.b0    = post.z0dinput.option.signe .* interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.b0 ,profil0d.temps(itime),'pchip','extrap'); 
z0_offset       = interp1_imas(post.z0dinput.cons.temps,post.z0dinput.geo.z0 ,profil0d.temps(itime),'pchip','extrap');
geo_input.sp    = interp1_imas(post.z0dinput.cons.temps,post.zerod.sp ,profil0d.temps(itime),'pchip','extrap');
geo_input.vp    = interp1_imas(post.z0dinput.cons.temps,post.zerod.vp ,profil0d.temps(itime),'pchip','extrap');
geo_input.sext  = interp1_imas(post.z0dinput.cons.temps,post.zerod.sext ,profil0d.temps(itime),'pchip','extrap');
  
if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa') && all(isfinite(profil0d.Rsepa(itime,:))) && all(isfinite(profil0d.Zsepa(itime,:)))
  geo_input.Rsepa = profil0d.Rsepa(itime,:);
  geo_input.Zsepa = profil0d.Zsepa(itime,:);
  % safety rule: must be convex
  if isfield(post.z0dinput.option,'Convex_LCFS') && (post.z0dinput.option.Convex_LCFS == 1)
      KH = sort(unique(convhull(geo_input.Rsepa,geo_input.Zsepa)));
      if (length(KH) ~= length(geo_input.Rsepa))
	  %figure(21);clf;
	  %plot(geo_input.Rsepa,geo_input.Zsepa,'b',geo_input.Rsepa(KH),geo_input.Zsepa(KH),'.r');
	  %disp([length(KH),length(geo_input.Rsepa)])
	  index_full = 1:length(geo_input.Rsepa);
	  Rsepa = geo_input.Rsepa(KH);
	  Zsepa = geo_input.Zsepa(KH);
	  geo_input.Rsepa = interp1(KH,Rsepa,index_full,'linear');
	  geo_input.Zsepa = interp1(KH,Zsepa,index_full,'linear');
	  indbad_lcfs = find(~isfinite(geo_input.Rsepa) | ~isfinite(geo_input.Zsepa));
	  if ~isempty(indbad_lcfs)
	      geo_input.Rsepa(indbad_lcfs) = [];
	      geo_input.Zsepa(indbad_lcfs) = [];
	  end
	  %hold on;
	  %plot(geo_input.Rsepa,geo_input.Zsepa,'k');
	  %drawnow
	  %geo_input.Rsepa = geo_input.Rsepa(KH);
	  %geo_input.Zsepa = geo_input.Zsepa(KH);
      end
  end
  geo_input.z0    = (max(profil0d.Zsepa(itime,:)) + min((profil0d.Zsepa(itime,:)))) ./ 2;
  geo_input.Zsepa = geo_input.Zsepa - geo_input.z0; 
  box.a    = max(post.z0dinput.geo.a);
  box.z0   = sum(post.z0dinput.geo.z0 .* post.z0dinput.cons.ip) ./ sum(post.z0dinput.cons.ip);
  box.rmin = min(profil0d.Rsepa(:)) - box.a/2;
  box.rmax = max(profil0d.Rsepa(:)) + box.a/2;
  box.zmin = min(profil0d.Zsepa(:)) - box.a/2*rzp + box.z0;
  box.zmax = max(profil0d.Zsepa(:)) + box.a/2*rzp + box.z0;

else
  geo_input.Rsepa = [];
  geo_input.Zsepa = []; 
  post.z0dinput.option.morphing = 0;
  geo_input.z0 = 0;
  box.a    = max(post.z0dinput.geo.a);
  box.z0   = sum(post.z0dinput.geo.z0 .* post.z0dinput.cons.ip) ./ sum(post.z0dinput.cons.ip);
  box.rmin = min(post.z0dinput.geo.R - post.z0dinput.geo.a) - box.a/2;
  box.rmax = max(post.z0dinput.geo.R + post.z0dinput.geo.a) + box.a/2;
  box.zmin = min(- post.z0dinput.geo.a .* post.z0dinput.geo.K) - box.a/2*rzp + box.z0;
  box.zmax = max(post.z0dinput.geo.a .* post.z0dinput.geo.K) + box.a/2*rzp + box.z0;
end

% normalisation volume
prof.volume  = prof.volume ./ prof.volume(end)  .* geo_input.vp;
prof.surface = prof.surface ./ prof.surface(end) .* geo_input.sp;
% derivative of psi
% dpsidx s'annule au centre
prof.psid1    = pdederive(prof.x,prof.psi,0,2,2,1);
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
prof.psid2    = pdederive(prof.x,prof.psi,1,0,2,2);
if post.z0dinput.option.cronos_regul == 4
	  prof.psid2(1) = - geo_input.b0 ./ prof.q(1) .* prof.rho(end) .^ 2;
end
%
[profil,deuxd,moment,Ip_lcfs] = metis2equi1t_imas(prof,geo_input,phys,sign(factor_two_pi) .* nbth,post.z0dinput.option.morphing,factor_two_pi);
deuxd.Z = deuxd.Z + geo_input.z0;
moment.zaxe = moment.zaxe   +  geo_input.z0;  
Iprof = cumtrapz(profil0d.xli,profil0d.spr(itime,:) .* profil0d.jli(itime,:),2);
Iprof_bis = profil0d.rmx(itime,end) .* cumtrapz(profil0d.xli,profil0d.vpr_tor(itime,:) .* profil0d.jli(itime,:) ./ 2 ./ pi .* profil0d.ri(itime,:),2);
ip_loc  = interp1_imas(post.zerod.temps,post.zerod.ip,profil0d.temps(itime),'pchip','extrap');
delta_ip = ip_loc - Ip_lcfs(end);
if delta_ip > 0
  ip_loc_error_upper = delta_ip;
  ip_loc_error_lower = 0;
else
  ip_loc_error_upper = 0;
  ip_loc_error_lower = delta_ip;
end
error_2d = abs(ip_loc - Ip_lcfs(end)) ./ max(1,abs(ip_loc));
fprintf(' error_2D = %g | ',error_2d);
% rescale to preserve plasma current even with high pedestal
Iprof  = abs(Iprof ./ Iprof(end) .* ip_loc);
factor = abs(Iprof ./ max(1,abs(Ip_lcfs)));
factor(1) = 1;
%figure(119);clf;plot(factor);drawnow
vth            = ones(1,size(deuxd.dPSIdx,2));
deuxd.dPSIdx   = (factor' * vth) .* deuxd.dPSIdx;
deuxd.BR       = (factor' * vth) .* deuxd.BR;
deuxd.BZ       = (factor' * vth) .* deuxd.BZ;
deuxd.dPSIdR   = (factor' * vth) .* deuxd.dPSIdR; 
deuxd.dPSIdZ   = (factor' * vth) .* deuxd.dPSIdZ; 

% there is a problem with z0
deuxd.Z     = deuxd.Z + z0_offset;
moment.zaxe = moment.zaxe   + z0_offset;

% remplissage de la structure pour le precalcul de la grille rectangulaire
% remplissage de la structure equivide
% passage en coodornnees de flux
xout = profil.rho ./ max(profil.rho);  
% calcul de rhomax
equicronos.rhomax          = profil.rho(end);
% utilisation des donnnees du mapping pour les indice de temps problematique
equicronos.phi(1,:)             = profil.phi; 
equicronos.psi(1,:)             = profil.psi ./ factor_two_pi; 
% le R et le Z
equicronos.R(1,:,:)          =  shiftdim(deuxd.R,-1);
equicronos.Z(1,:,:)          =  shiftdim(deuxd.Z,-1) - moment.zaxe(1);
equicronos.rhoRZ(1,:)        =  profil.rho';
equicronos.psiRZ(1,:)        =  profil.psi'./ factor_two_pi;
equicronos.df2RZ(1,:)        =  2 .* profil.fdia' .* pdederive(profil.psi,profil.fdia,2,2,2,1)';
equicronos.dprRZ(1,:)        =  pdederive(profil.psi,profil.ptot,2,2,2,1)';
% La carte de champ
equicronos.BR(1,:,:)     = shiftdim(deuxd.BR,-1) .* sign(factor_two_pi);
equicronos.BZ(1,:,:)     = shiftdim(deuxd.BZ,-1) .* sign(factor_two_pi);
equicronos.BPHI(1,:,:)   = shiftdim(deuxd.BPHI,-1);
  
% equilibre etendu hors du plasma
% EVALUATION IN THE VACCUUM
method_extrap = NaN;
if equi_extrap == 0
  equivide = zfitvide_interp(equicronos,0,[]);
  method_extrap = 0;
else
  method_extrap = 1;
  if ~isempty(control) && ~isempty(control.x_point)
    xpoints.R = control.x_point{1}.r;
    xpoints.Z = control.x_point{1}.z;
    if length(control.x_point) > 1
	xpoints.R(2) = control.x_point{2}.r;
	xpoints.Z(2) = control.x_point{2}.z;
    end
    equivide = zfitvide_jorek(equicronos,0,[],xpoints);   
  else
    equivide = zfitvide_jorek(equicronos,0,[]); 
  end
  if ((equivide.dpsi_error > 0.1) || (equivide.dB_error > 0.2)) && (equi_extrap <= 1)
      fprintf('changing to direct methode |');
      method_extrap = 0;
      dpsi_error_mem = equivide.dpsi_error;
      equivide = zfitvide_interp(equicronos,0,[]);
      equivide.dpsi_error = equivide.dpsi_error + dpsi_error_mem;
  end 
end
if method_extrap == 0
  [psi2d,BR2d,BZ2d] = zpsivide_interp(equivide,Rfbe,Zfbe); 
else 
  [psi2d,BR2d,BZ2d] = zpsivide_updown_16_jorek(equivide,Rfbe,Zfbe); 
end
BR2d  = BR2d .* sign(factor_two_pi);
BZ2d  = BZ2d .* sign(factor_two_pi);
psi2d = psi2d .* factor_two_pi;
mask = zinout(equivide.R(end,:),equivide.Z(end,:),Rfbe,Zfbe);
BPHI2d   = profil0d.fdia(itime,end) ./ Rfbe;
BPHI2d(mask)  = interp1(profil0d.psi(itime,:),profil0d.fdia(itime,:),psi2d(mask),'pchip','extrap') ./ Rfbe(mask);
% undefined at R = 0
psi2d(~isfinite(psi2d)) = 0;


function control =sepa_moments(Rsepa,Zsepa)

% calcul des moments
% calcul de R0 et Z0
maskrmax  = (Rsepa == max(Rsepa,[],2));
% recalcul des parametres sur le vecteur final
rmin  = min(Rsepa,[],2);
rmax  = max(Rsepa,[],2);
a = 0.5 .* (rmax - rmin);
R0 = 0.5 .* (rmax + rmin);
control.R0  = R0;
control.a   = a;
zmin  = min(Zsepa,[],2);
zmax  = max(Zsepa,[],2);
control.K    = (zmax - zmin) ./ 2 ./ a;
rzmax = Rsepa(min(find(Zsepa == zmax)));
rzmin = Rsepa(min(find(Zsepa == zmin)));
z0    = Zsepa(min(find(Rsepa == rmax)));
% IMAS definition
control.z0_geo= (zmax + zmin) ./ 2;
control.z0    = z0;
control.Ku    = (zmax - z0)  ./ a;
control.Kl    = (z0 - zmin)  ./ a;
control.d     = abs(rzmax + rzmin -  2 .* R0) ./ 2 ./ a;
control.du    = abs(rzmax - R0) ./ a;
control.dl    = abs(rzmin - R0) ./ a;

% Xpoint detection
indh     = find(Zsepa == zmax,1);
indh     = cat(2,indh - 3, indh - 2,indh - 1,indh,indh + 1,indh + 2,indh + 3); 
indh     = mod(indh-1,size(Zsepa,2))+1;
rh       = Rsepa(indh);
zh       = Zsepa(indh);
ph       = polyfit(rh,zh,2);
eh       = sqrt(mean((zh - polyval(ph,rh)).^ 2)) ./ (max(rh) - min(rh));
indl     = find(Zsepa == zmin,1);
indl     = cat(2,indl - 3, indl - 2,indl - 1,indl,indl + 1,indl + 2,indl + 3);
indl     = mod(indl-1,size(Zsepa,2))+1;
rl       = Rsepa(indl);
zl       = Zsepa(indl);
pl       = polyfit(rl,zl,2);
el       = sqrt(mean((zl - polyval(pl,rl)).^ 2)) ./ (max(rl) - min(rl));
control.x_point = {};
if el > 2e-2
  indlz  = find(zl == min(zl),1);
  control.x_point{end + 1}.r = rl(indlz);
  control.x_point{end}.z = zl(indlz);
end
if eh > 2e-2
  indhz  = find(zh == max(zh),1);
  control.x_point{end + 1}.r = rh(indhz);
  control.x_point{end}.z = zh(indhz);
end


% calcul de la solution dans le vide par decompostion en serie des fonctions de green et fit sur le bord
% permet de prolonge la solution a l'exterieur de la separatrice.
% pour calculer en (R,Z), il faut appele, ensuite : [psi_out,BR_out,BZ_out] = zpsivide_updown_16(equivide,R,Z);
function equivide = zfitvide_jorek(equi,order,plotonoff,xpoints,tolerance,methode,smooth_psi,geo)

if nargin < 2
	order = [];
end
if nargin < 3
	plotonoff = 0;
end
if nargin < 4 
	xpoints = [];
end

if nargin < 5	
	tolerance = 0;
elseif isempty(tolerance)
	tolerance = 0;
end
if nargin < 6
	methode = 0;
elseif isempty(methode)
	methode = 0;
end
if (methode == 0) &&  isempty(order)
	order = 0;
elseif isempty(order)
	order = ceil(size(equi.R,3) ./ 10);
end

if nargin < 7
	smooth_psi = 0;
elseif isempty(smooth_psi)
	smooth_psi = 0;
end


if ~all(isfinite(equi.psiRZ))
   x      = linspace(0,1,length(equi.psi));
   rhoh   = equi.rhoRZ;
   xh     = rhoh ./ rhoh(end);
   equi.psiRZ = pchip(x,equi.psi,xh);
else 
   equi.psiRZ = equi.psiRZ;
end

% nombre de coefficient
equivide = zpsivide_updown_16_jorek(order);
% extraction des infos sur la DSMF
equivide.R     = squeeze(equi.R);
equivide.Z     = squeeze(equi.Z);
equivide.BR    = squeeze(equi.BR);
equivide.BZ    = squeeze(equi.BZ);
equivide.PSI   = equi.psiRZ' * ones(1,size(equivide.R,2));
equivide.signe = 1;

% points de controle
indd = size(equivide.R,1);
indc = 1:(size(equivide.R,2) - 1);
R     = equivide.R(indd,indc);
Z     = equivide.Z(indd,indc);
BR    = equivide.BR(indd,indc);
BZ    = equivide.BZ(indd,indc);
PSI   = equivide.PSI(indd,indc);
% contrainte au centre (pour donnee le sens de psi)
RC    = equivide.R(1,1);
ZC    = equivide.Z(1,1);
PSIC  = equivide.PSI(1,1);

% separatrice
if nargin > 7
    RS    = geo.R(:);
    ZS    = geo.Z(:);
    PSIS  = equi.psiRZ(end) * ones(size(RS));
else
    RS    = R;
    ZS    = Z;
    PSIS  = PSI;       
end
%
%equivide.Rnorm = mean(RS);
%equivide.Zoffset = mean(ZS);
equivide.Rnorm = RC;
equivide.Zoffset = ZC;


% creation de la matrice a inversee
if ~isempty(xpoints)
	mpol = NaN * ones(prod(size(RS)) + 2 .* prod(size(R)) + 2 .* length(xpoints.R) + 3 ,length(equivide.coef));	
else
	mpol = NaN * ones(prod(size(RS)) + 2 .* prod(size(R)) + 3,length(equivide.coef));
	%mpol_test = NaN * ones(prod(size(RS)) + 3,length(equivide.coef));
end
w_axe = 1;
for k = 1:length(equivide.coef)
  equivide.coef = zeros(size(equivide.coef));
  equivide.coef(k) = 1;
  [psix,BRx,BZx] = zpsivide_updown_16_jorek(equivide,R,Z,1); 
  [psic,BRc,BZc] = zpsivide_updown_16_jorek(equivide,RC,ZC,1); 
  [psis,BRs,BZs] = zpsivide_updown_16_jorek(equivide,RS,ZS,1); 
  if ~isempty(xpoints)
	[psi0,BR0,BZ0] = zpsivide_updown_16_jorek(equivide,xpoints.R,xpoints.Z,1);
	mpol(:,k) = cat(1,psis(:),BRx(:),BZx(:),psic .* w_axe,BRc .* w_axe,BZc .* w_axe,equivide.order .* BR0(:), equivide.order .* BZ0(:));
  else
	mpol(:,k) = cat(1,psis(:),BRx(:),BZx(:),psic .* w_axe,BRc .* w_axe,BZc .* w_axe);
	%mpol_test(:,k) = cat(1,psis(:),psic,BRc,BZc);
  end
end
if ~isempty(xpoints)
	v0 = zeros(1,length(xpoints.R));
	ypol = cat(1,PSIS(:),BR(:),BZ(:),PSIC .* w_axe,0,0,v0(:),v0(:));	
else
	ypol = cat(1,PSIS(:),BR(:),BZ(:),PSIC .* w_axe,0,0);
	%ypol_test = cat(1,PSIS(:),PSIC,0,0);
end
% resolution au sens des moindres carres
if tolerance == 0
	equivide.coef =  mpol \ ypol;
else
	equivide.coef =  pinv(mpol,tolerance)  * ypol;
	if plotonoff > 0
		[u,s,v] = svd(mpol,'econ');
		h = findobj(0,'type','figure','tag','equivide3');
		if isempty(h)
		h=figure('tag','equivide3');
		else
		figure(h);
		end   
		clf
		set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
			'defaultlinelinewidth',1,'color',[1 1 1])

		semilogy(diag(s),'o');
	end
end
% check for orientation
[psix,BRx,BZx] = zpsivide_updown_16_jorek(equivide,RC,ZC,1);
[psis,BRs,BZs] = zpsivide_updown_16_jorek(equivide,RS,ZS,1);
if sign(mean(psis) - psic) ~= sign(mean(PSIS) - PSIC)
  error_sign = 1e38;
else
  error_sign = 0;
end

% traitement de l discontinuite
% recalcul 
[psix,BRx,BZx] = zpsivide_updown_16_jorek(equivide,R,Z,1);
[psis,BRs,BZs] = zpsivide_updown_16_jorek(equivide,RS,ZS,1);
equivide.BR_lcfs_delta     = BR - BRx;
equivide.BZ_lcfs_delta     = BZ - BZx;
equivide.psi_lcfs_delta    = PSIS - psis;
equivide.dpsi_error = mean(abs(equivide.psi_lcfs_delta)) ./ abs(PSI(1) -PSIC) + error_sign;
equivide.dB_error   = sqrt(mean((equivide.BR_lcfs_delta .^ 2 + equivide.BZ_lcfs_delta .^ 2) ./ (BR .^ 2 + BZ .^ 2)));
%equivide

fprintf('dpsi_error = %g & dB_error = %g |',equivide.dpsi_error,equivide.dB_error);

%  figure(35);clf
%  plot(1:length(BR),BR,'r',1:length(BRx),BRx,'m',1:length(BZ),BZ,'b',1:length(BZx),BZx,'c')
%  drawnow


if isempty(plotonoff)
    return
elseif plotonoff ~= 1
    return
end


% recalcul 
[psix,BRx,BZx] = zpsivide_updown_16_jorek(equivide,R,Z,1);
% evaluation dans le vide
a   = (max(R(:)) - min(R(:))) ./ 2;
rzp = (max(Z(:)) - min(Z(:))) ./ a ./ 2;
nb  = 51;
[Ra,Za] = meshgrid(linspace(a ./ 2 ,max(R(:)) + a .* 1.5,nb), ...
                   linspace(min(Z(:))-  a   .* rzp,max(Z(:)) + a  .* rzp,ceil(nb .* rzp)));
% avec interpolation ->
[psia,BRa,BZa] =  zpsivide_updown_16_jorek(equivide,Ra,Za);

% sans interpolation ->
%[psia,BRa,BZa] =  zpsivide(equivide,Ra,Za,1);

h = findobj(0,'type','figure','tag','equivide1');
if isempty(h)
       h=figure('tag','equivide1');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(2,2,1)
plot((psix  - PSI) ./abs(equi.psiRZ(1) - equi.psiRZ(end)),'r');
ylabel('flux error @ LCMS')

subplot(2,2,2)
plot(indc,BR,'b',indc,BRx,'r')
ylabel('BR  @ LCMS')
legend('target','ident')
subplot(2,2,3)
plot(indc,BZ,'b',indc,BZx,'r')
ylabel('BZ   @ LCMS')
legend('target','ident')

subplot(2,2,4)
plot(indc,BR .^ 2 + BZ .^ 2,'b',indc,BRx .^ 2 + BZx .^ 2,'r')
ylabel('Bpol^2 @ LCMS')
legend('target','ident')

Rh = squeeze(equi.R);
Zh = squeeze(equi.Z);

h = findobj(0,'type','figure','tag','equivide2');
if isempty(h)
       h=figure('tag','equivide2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
	
x    = linspace(0,1,length(equi.psiRZ));
xx   = linspace(1,3,101);
psiu = cat(2,equi.psiRZ(6:5:end),pchip(x,equi.psiRZ,xx));	
if smooth_psi == 1
	contour(Ra,Za,zsmooth_psi(psia),psiu);
else
	contour(Ra,Za,psia,psiu);
end

hold on
plot(Rh(1:5:end,:)',Zh(1:5:end,:)','linestyle',':');
contour(Ra,Za,zsmooth_psi(psia),[1-eps,1+eps].*equi.psiRZ(end),'color','g');
if smooth_psi == 1
	contour(Ra,Za,zsmooth_psi(psia),[1-eps,1+eps].*equi.psiRZ(end),'color','g');
else
	contour(Ra,Za,psia,[1-eps,1+eps].*equi.psiRZ(end),'color','g');
end
if ~isempty(xpoints)
	plot(xpoints.R,xpoints.Z,'x');
end
axis('equal');
xlabel('R (m)');
xlabel('Z (m)');
title('iso flux');
keyboard

% fonction appele par zfitvide
function [psi_out,BR_out,BZ_out] = zpsivide_updown_16_jorek(equivide,Rin,Zin,fit)


% calcul du  nombre de coef (coef en entree donne l'ordre)
ind = 0;
if nargin <= 1
	psi_out.R   = [];
	psi_out.Z   = [];
	psi_out.BR  = [];
	psi_out.BZ  = [];
	psi_out.PSI = [];
	psi_out.Rnorm = [];
	psi_out.Zoffset = [];	
	psi_out.coef = zeros(1,16 + 2 .* equivide);
	psi_out.order = equivide;
	fprintf(' plasma up/down asymetry | ')
 	return
end


% separation
s1   = equivide.coef(1);
s2   = equivide.coef(2);
A    = equivide.coef(3);
B    = equivide.coef(4);
cn   = equivide.coef(5:16);
if equivide.order > 0
	gn  = equivide.coef(17:end);
	gn1 = gn(1:2:end-1);
	gn2 = gn(2:2:end);
end

% pour test
%  fcn  =zeros(size(cn));
%  ind = 1+fix(rand(1) .* length(fcn));
%  fcn(ind) = 1;
%  cn  = cn .* fcn;

% dimension
nrz = size(Rin);
Rin  = Rin(:) ./ equivide.Rnorm;
Zin  = (Zin(:) - equivide.Zoffset) ./ equivide.Rnorm;

% init out
psi_out = NaN .* ones(size(Rin));
BR_out  = NaN .* ones(size(Rin));  
BZ_out  = NaN .* ones(size(Rin));

% point interieur
if nargin < 4
	notinside = 1;
else
	notinside = 0;
end

mask = zeros(size(Rin));
R    = Rin(~mask);
Z    = Zin(~mask);

% init (la partie en a + b * r^2) est dans le developpment
% solution trivial de Psi
rr   = sqrt(R .^ 2 + Z .^ 2);
psiv = s1 .* rr;
BRv  = - s1 .* Z ./ rr ./ R;
BZv  = s1  ./ rr;

% nouvelle contribution (groupe de symetrie de l'equation)
% ref : Y.E. Litvinenko  PoP 17, 074502 (2010). 
psiv  = psiv + s2 .* (Z .* sqrt(R .^ 2 + Z .^ 2) + R .^ 2 .* asinh(Z./R));
BRv   = BRv  - s2 .* 2 .* sqrt(Z .^ 2 + R .^ 2) ./ R;
BZv   = BZv  + s2 .* 2.* asinh(Z ./ R);


% calcul (partie courant)
% solution soloviev
% ref : P. J. Mc Carthy POP 1999 p 3554-...
psiv = psiv +  A .* R .^ 4 ./ 8 + B .* Z .^2 ./ 2;
BRv  = BRv  - B .* Z ./ R; 
BZv  = BZv +  A .* R .^ 2 ./ 2;


% calcul (partie contribution du vide)
% ref :
% Analytical tokamak equilibrium for shaped plasmas
% S. B. Zheng,a) A. J. Wootton, and Emilia R. Solano
% Phys. Plasmas 3 (3), March 1996

% A.J. Cerfon and J.P. Freidberg, POP 17, 032502 (2010)
% choix R0 = 1
% 1
psi_1 = 1;
BR_1  = 0;
BZ_1  = 0;

% 2 
psi_2 = R .^ 2; 
BR_2  = 0;
BZ_2  = 2;

% 3
psi_3 = Z .^ 2  - R .^ 2 .* log(R);
BR_3  = - 2 .* Z ./ R;
BZ_3  = - 2 .* log(R) - 1;

%4 
psi_4 = R .^ 4 - 4 .* R .^ 2 .* Z .^ 2;
BR_4  = 8 .* R .* Z;
BZ_4  = 4 .* R .^ 2 - 8 .* Z .^ 2;

%5
psi_5 = 2 .* Z .^ 4  - 9 .* R .^ 2 .* Z .^ 2  + 3 .* R .^ 4  .* log(R)  - 12 .* R .^ 2 .* Z .^ 2 .* log(R);
BR_5  = - (8 .* Z .^ 3 + (-24 .* R .^ 2 .* log(R) - 18 .* R .^ 2) .* Z) ./ R;
BZ_5  =  -(24 .* log(R) + 30) .* Z .^ 2 + 12 .* R .^ 2 .* log(R) + 3 .* R .^ 2;

%6
psi_6 = R .^ 6  - 12 .* Z .^ 2 .* R .^ 4 + 8 .* Z .^ 4 .* R .^ 2;
BR_6  = -(32 .* R .* Z .^ 3 - 24 .* R .^ 3 .* Z);
BZ_6  = 16 .* Z .^ 4 - 48 .* R .^ 2 .* Z .^ 2 + 6 .* R .^ 4;

%
psi_7 = 8 .* Z .^ 6  - 140 .* Z .^ 4 .* R .^ 2 + 75 .* Z .^ 2 .* R .^ 4  - 15  .* R .^ 6 .* log(R) + 180 .* R .^ 4 .* Z .* 2 .* log(R) - 120 .* R .^ 2 .* Z .^ 4 .* log(R); 
BR_7  = -(48 .* Z .^ 5 + (-480 .* R .^ 2 .* log(R) - 560 .* R .^ 2) .* Z .^ 3 + 150 .* R .^ 4 .* Z + 360 .* R .^ 4 .* log(R)) ./ R;
BZ_7  = (-240 .* log(R) - 400) .* Z .^ 4 + 300 .* R .^ 2 .* Z .^ 2 + (1440 .* R .^ 2 .* log(R) + 360 .* R .^ 2) .* Z - 90 .* R .^ 4 .* log(R) - 15 .* R .^ 4;


% 8
psi_8 = Z;
BR_8  = - 1./ R;
BZ_8  = 0;

% 9
psi_9 = Z .* R .^ 2;
BR_9  = - R;
BZ_9  = 2.* Z;

% 10
psi_10 = Z .^ 3 - 3 .* Z .* R .^ 2 .* log(R);
BR_10  = -(3 .* Z .^ 2 - 3 .* R .^ 2 .* log(R)) ./ R;
BZ_10  = (-6 .* log(R) - 3) .* Z;

% 11
psi_11 = 3 .* Z .* R .^ 4  - 4 .* Z .^ 3 .* R .^ 2;
BR_11  = - 3 .* R .^ 3 + 12 .* R .* Z .^ 2;
BZ_11 = 12 .* R .^ 2 .* Z - 8 .* Z .^ 3;

%12
psi_12 = 8 .* Z .^ 5 - 45 .* Z .* R .^ 4 - 80 .* Z .^ 3 .* R .^ 2 .* log(R) + 60 .* Z .* R .^ 4 .* log(R); 
BR_12  = - (40 .* Z .^ 4 - 240 .* R .^ 2 .* log(R) .* Z .^ 2 + 60 .* R .^ 4 .* log(R) - 45 .* R .^ 4) ./ R;
BZ_12  = (-160 .* log(R) - 80) .* Z .^ 3 + (240 .* R .^ 2 .* log(R) - 120 .* R .^ 2) .* Z;

% sommation
psiv = psiv + cn(1) .* psi_1 + cn(2)  .* psi_2 +  cn(3)  .* psi_3  + cn(4)  .* psi_4 + ...
              cn(5) .* psi_5 + cn(6)  .* psi_6 +  cn(7)  .* psi_7  + cn(8)  .* psi_8 + ...
              cn(9) .* psi_9 + cn(10) .* psi_10 + cn(11) .* psi_11 + cn(12) .* psi_12;

BRv = BRv + cn(1) .* BR_1 + cn(2)  .* BR_2 +  cn(3)  .* BR_3  + cn(4)  .* BR_4 + ...
            cn(5) .* BR_5 + cn(6)  .* BR_6 +  cn(7)  .* BR_7  + cn(8)  .* BR_8 + ...
            cn(9) .* BR_9 + cn(10) .* BR_10 + cn(11) .* BR_11 + cn(12) .* BR_12;

BZv = BZv + cn(1) .* BZ_1 + cn(2)  .* BZ_2 +  cn(3)  .* BZ_3  + cn(4)  .* BZ_4 + ...
            cn(5) .* BZ_5 + cn(6)  .* BZ_6 +  cn(7)  .* BZ_7  + cn(8)  .* BZ_8 + ...
            cn(9) .* BZ_9 + cn(10) .* BZ_10 + cn(11) .* BZ_11 + cn(12) .* BZ_12;


if equivide.order > 0
	% calcul (partie courant)
	% fonction de green developpee en series pour un filament de courant
	% chaque moment est libre 
	% ref : P. J. Mc Carthy POP 1999 p 3554-...
	for k = 1:length(gn1)
			l = k + 2;
			% 
			[u,upr,upz] = f(l,R,Z,gn1(k),gn2(k));
			%
			psiv = psiv +  R .^ l .* u;
			BRv  = BRv  - (R .^ l .* upz) ./ R;
			BZv  = BZv  + (R .^ l .* upr + l .* R .^ (l - 1) .* u) ./ R;
	end
end


% probleme de convention de signe 
if equivide.signe == 1
  BRv = - BRv;
  BZv = - BZv;
end
% recopie dans les sortie
psi_out(~mask) = psiv;
BR_out(~mask)  = BRv ./ equivide.Rnorm .^ 2;
BZ_out(~mask)  = BZv ./ equivide.Rnorm .^ 2;


% test & calcul dans le plasma

% donnees suplementaires
if notinside == 1

%  % resample on equi grid
%  FPSI = scatteredInterpolant(Rin,Zin,psiv,'natural','linear');
%  FBR  = scatteredInterpolant(Rin,Zin,BRv,'natural','linear');
%  FBZ  = scatteredInterpolant(Rin,Zin,BZv,'natural','linear');
%  PSI_ref = FPSI(equivide.R,equivide.Z);
%  BR_ref  = FBR(equivide.R,equivide.Z);
%  BZ_ref  = FBZ(equivide.R,equivide.Z);
%  
%  figure
%  contour(equivide.R,equivide.Z,equivide.PSI,101,'color','b','linestyle',':')
%  hold on
%  contour(equivide.R,equivide.Z,PSI_ref,101,'color','r','linestyle',':')
%  quiver(equivide.R,equivide.Z,equivide.BR,equivide.BZ,'color','b');
%  quiver(equivide.R,equivide.Z,BR_ref,BZ_ref,'color','r');
%  title('Bpoloidal: blue = METIS & red = Extrapolation');
%  xlabel('R (m)');
%  ylabel('Z (m)');
%  keyboard

      % unscale
      Rin  = Rin(:) .* equivide.Rnorm;
      Zin  = Zin(:) .* equivide.Rnorm + equivide.Zoffset;

      % rescale internal analytical map
      psi_LCFS_map = griddata(Rin,Zin,psi_out,equivide.R(end,:),equivide.Z(end,:),'cubic');
      offset = equivide.PSI(end,:) - psi_LCFS_map;
      if all(~isfinite(offset))
	offset = 0;
      else
	offset = mean(offset);
      end
      psi_out = psi_out + offset;
      fprintf('offset = %g | ',offset);
      
%       	psi_out_alt   = reshape(psi_out,nrz);
%  	BR_out_alt    = reshape(BR_out,nrz);
%  	BZ_out_alt    = reshape(BZ_out,nrz);
%  	Rin_alt      = reshape(Rin,nrz);
%  	Zin_alt       = reshape(Zin,nrz);
%  
%  	figure(34)
%  	clf
%  	contour(equivide.R,equivide.Z,equivide.PSI,equivide.PSI(:,1),'color','b','linestyle','--')
%  	hold on
%  	contour(Rin_alt,Zin_alt,psi_out_alt,equivide.PSI(:,1),'color','r','linestyle','-')
%  	quiver(equivide.R,equivide.Z,equivide.BR,equivide.BZ,'color','c');
%  	quiver(Rin_alt,Zin_alt,BR_out_alt,BZ_out_alt,'color','m');
%  	title('Bpoloidal: blue = METIS & red = Extrapolation');
%  	xlabel('R (m)');
%  	ylabel('Z (m)');
%  	plot(equivide.R(end,:),equivide.Z(end,:),'.k');
%  	drawnow
 
      
      % donnees helena
      Req   = equivide.R(2:end,1:end-1);
      Zeq   = equivide.Z(2:end,1:end-1);
      BReq  = equivide.BR(2:end,1:end-1);
      BZeq  = equivide.BZ(2:end,1:end-1);
      PSIeq = equivide.PSI(2:end,1:end-1);
      pas   = sqrt(prod(size(Rin)));  
      Lplus = sqrt((Req(end,:) - equivide.R(1,1)) .^ 2 + (Zeq(end,:) - equivide.Z(1,1)) .^ 2);
      Rplus = Req(end,:) + (Req(end,:) - equivide.R(1,1)) ./ Lplus ./ pas;
      Zplus = Zeq(end,:) + (Zeq(end,:) - equivide.Z(1,1)) ./ Lplus ./ pas;
      [psi_plus,BR_plus,BZ_plus] = zpsivide_updown_16_jorek(equivide,Rplus,Zplus,1);
      Req   = cat(1,equivide.R(1,1),Req(:),Rplus(:));
      Zeq   = cat(1,equivide.Z(1,1),Zeq(:),Zplus(:));
      BReq  = cat(1,equivide.BR(1,1),BReq(:),BR_plus(:));
      BZeq  = cat(1,equivide.BZ(1,1),BZeq(:),BZ_plus(:));
      PSIeq = cat(1,equivide.PSI(1,1),PSIeq(:),psi_plus(:));

      BR_out_eq  = griddata(Req,Zeq,BReq,Rin,Zin,'cubic');
      BZ_out_eq  = griddata(Req,Zeq,BZeq,Rin,Zin,'cubic');
      psi_out_eq = griddata(Req,Zeq,PSIeq,Rin,Zin,'cubic');
      disp('interpolation');
      if equivide.PSI(1,1) < equivide.PSI(end,1)
            frac = (psi_out - max(PSIeq(:))) ./ (min(PSIeq(:)) - max(PSIeq(:)));
      else
            frac = (psi_out - min(PSIeq(:))) ./ (max(PSIeq(:)) - min(PSIeq(:)));
      end
      frac = max(0,min(1,frac));
      frac(~isfinite(psi_out_eq)) = 0;
      frac = frac ./ max(frac(:));
      psi_out_eq(~isfinite(psi_out_eq)) = 0;
      BR_out_eq(~isfinite(BR_out_eq)) = 0;
      BZ_out_eq(~isfinite(BZ_out_eq)) = 0;     
      psi_out = (1- frac) .* psi_out + frac .* psi_out_eq;
      BR_out  = (1- frac) .* BR_out + frac .* BR_out_eq;
      BZ_out  = (1- frac) .* BZ_out + frac .* BZ_out_eq;

%        if numel(Rin)> 3
%  	% mise en forme
%  	psi_out   = reshape(psi_out,nrz);
%  	BR_out    = reshape(BR_out,nrz);
%  	BZ_out    = reshape(BZ_out,nrz);
%  	Rin       = reshape(Rin,nrz);
%  	Zin       = reshape(Zin,nrz);
%  
%  	figure(31)
%  	clf
%  	contour(equivide.R,equivide.Z,equivide.PSI,equivide.PSI(:,1),'color','b','linestyle','--')
%  	hold on
%  	contour(Rin,Zin,psi_out,equivide.PSI(:,1),'color','r','linestyle','-')
%  	contour(Rin,Zin,psi_out,101,'color','g','linestyle','-')
%  	quiver(equivide.R,equivide.Z,equivide.BR,equivide.BZ,'color','c');
%  	quiver(Rin,Zin,BR_out,BZ_out,'color','m');
%  	title('Bpoloidal: blue = METIS & red = Extrapolation');
%  	xlabel('R (m)');
%  	ylabel('Z (m)');
%  	plot(equivide.R(end,:),equivide.Z(end,:),'.k');
%  	drawnow
%  	keyboard
%        end
end

% mise en forme
psi_out   = reshape(psi_out,nrz);
BR_out    = reshape(BR_out,nrz);
BZ_out    = reshape(BZ_out,nrz);


% fonction appele par zfitvide_interp
function [psi_out,BR_out,BZ_out] = zfitvide_interp(equivide,varargin)
  psi_out = zpsivide_interp(equivide);

% fonction appele par zfitvide_interp
function [psi_out,BR_out,BZ_out] = zpsivide_interp(equivide,Rin,Zin,fit)


if nargin == 1
  equivide.R     = squeeze(equivide.R);
  equivide.Z     = squeeze(equivide.Z);
  equivide.BR    = squeeze(equivide.BR);
  equivide.BZ    = squeeze(equivide.BZ);
  equivide.PSI   = equivide.psiRZ' * ones(1,size(equivide.R,2));
  % adding zero on D far outside
  R   = equivide.R(:);
  Z   = equivide.Z(:);
  PSI = equivide.PSI(:);
  %
  Z0  = equivide.Z(1,1);
  R0  = equivide.R(1,1);
  a   = (max(R) - min(R)) ./ 2;
  Rplus   = zeros(1,size(equivide.R,2));
  Zplus   = zeros(1,size(equivide.R,2));
  PSIplus = zeros(1,size(equivide.R,2));
  for k=1:size(equivide.R,2)
       l = a ./ abs((equivide.R(end,k) - R0) + sqrt(-1) .* (equivide.Z(end,k) - Z0)) ./ size(equivide.R,1);
       Rplus(k) = (1+l) .* (equivide.R(end,k) - R0) + R0;
       Zplus(k) = (1+l) .* (equivide.Z(end,k) - Z0) + Z0;
       PSIplus(k) = equivide.PSI(end,k)  - equivide.R(end,k) .* (equivide.BZ(end,k) .* (Rplus(k) -  equivide.R(end,k)) - ...
                    equivide.BR(end,k) .* (Zplus(k) -  equivide.Z(end,k))); 
  end
  equivide.Rext   = cat(1,R,Rplus');
  equivide.Zext   = cat(1,Z,Zplus');
  equivide.PSIext = cat(1,PSI,PSIplus');
  %figure(117);clf;plot(equivide.Rext,equivide.Zext);drawnow
  warning off
  if exist('scatteredInterpolant')
	equivide.F_PSI = scatteredInterpolant(equivide.Rext,equivide.Zext,equivide.PSIext,'natural','linear');
  else
	equivide.F_PSI = TriScatteredInterp(equivide.Rext,equivide.Zext,equivide.PSIext,'natural');  
  end
  warning on
  
  % traitement de l discontinuite
  % recalcul 
  psis =  equivide.F_PSI(equivide.R(end,:),equivide.Z(end,:));
  equivide.psi_lcfs_delta    = equivide.PSI(end,:) - psis;
  equivide.dpsi_error = mean(abs(equivide.psi_lcfs_delta)) ./ abs( equivide.PSI(1,1) - equivide.PSI(end,1));
  dl      =  sqrt(eps) .* mean(equivide.R(:));
  BRs     =   (equivide.F_PSI(equivide.R(end,:),equivide.Z(end,:) + dl/2) - equivide.F_PSI(equivide.R(end,:),equivide.Z(end,:) - dl/2)) ./ equivide.R(end,:) ./ dl;
  BZs     = - (equivide.F_PSI(equivide.R(end,:) + dl/2,equivide.Z(end,:)) - equivide.F_PSI(equivide.R(end,:) - dl/2,equivide.Z(end,:))) ./ equivide.R(end,:) ./ dl;
  equivide.BR_lcfs_delta     = equivide.BR(end,:) - BRs;
  equivide.BZ_lcfs_delta     = equivide.BZ(end,:) - BZs;
  equivide.dB_error   = sqrt(mean((equivide.BR_lcfs_delta .^ 2 + equivide.BZ_lcfs_delta .^ 2) ./ (equivide.BR(end,:) .^ 2 + equivide.BZ(end,:) .^ 2)));
  fprintf('dpsi_error = %g & dB_error = %g |',equivide.dpsi_error,equivide.dB_error);
  psi_out = equivide;
  
else 
  %
  fprintf('interpolation\n');
  %
  nsin = size(Rin);
  Rin = Rin(:);
  Zin = Zin(:);
  mask_out = zinout(equivide.R(end,:),equivide.Z(end,:),Rin,Zin);
  dl      =  sqrt(eps) .* mean(equivide.R(:));
  warning off
  psi_out =   equivide.F_PSI(Rin,Zin);
  BR_out  =   (equivide.F_PSI(Rin,Zin + dl/2) - equivide.F_PSI(Rin,Zin - dl/2)) ./ Rin ./ dl;
  BZ_out  = - (equivide.F_PSI(Rin + dl/2,Zin) - equivide.F_PSI(Rin - dl/2,Zin)) ./ Rin ./ dl;
  psi_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.PSI(:),Rin(mask_out),Zin(mask_out),'cubic');
  BR_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.BR(:),Rin(mask_out),Zin(mask_out),'cubic');
  BZ_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.BZ(:),Rin(mask_out),Zin(mask_out),'cubic');
  warning on
  psi_out = reshape(psi_out,nsin);
  BR_out  = reshape(BR_out,nsin);
  BZ_out  = reshape(BZ_out,nsin);
  Rin  = reshape(Rin,nsin);
  Zin  = reshape(Zin,nsin);

%    figure(31)
%    clf
%    contour(equivide.R,equivide.Z,equivide.PSI,equivide.PSI(:,1),'color','b','linestyle','--')
%    hold on
%    contour(Rin,Zin,psi_out,equivide.PSI(:,1),'color','r','linestyle','-')
%    contour(Rin,Zin,psi_out,101,'color','g','linestyle','-') 
%    quiver(equivide.R,equivide.Z,equivide.BR,equivide.BZ,'color','c');
%    quiver(Rin,Zin,BR_out,BZ_out,'color','m');
%    title('Bpoloidal: blue = METIS & red = Extrapolation');
%    xlabel('R (m)');
%    ylabel('Z (m)');
%    plot(equivide.R(end,:),equivide.Z(end,:),'.k');
%    drawnow
%    keyboard
%    
%    if any(~isfinite(psi_out(:))) || any(~isfinite(BR_out(:))) || any(~isfinite(BZ_out(:)))
%      keyboard
%    end
end

	      
