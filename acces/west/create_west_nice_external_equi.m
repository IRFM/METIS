function create_west_nice_external_equi(shot)

% constant
mu0 = 4 * pi * 1e-7;

% test if there is some data
[ip,tip]= tsbase(shot,'smag_ip');
if isempty(ip)
    disp('No data')
    return
end
t_igni = tsbase(shot,'rignitron');
tplusdip = ftplusdip(shot);

% read equilibrium
equi = imas_west_get(shot,'equilibrium',0,1);
if isempty(equi) || isempty(equi.time)
  error('No available data for equilibrium');
end

% selection of valide time slices only
indok = find(equi.code.output_flag>=0);
if length(indok)< 3
    error('Not sufficient high number of valid equilibrium time slices available in this shot');
end

% computation of time derivative quantities
tequi   = equi.time(indok)' - t_igni;
psi1d   = - equi.profiles_1d.psi(indok,:)./ 2 ./ pi;
jphi_1d = equi.profiles_1d.j_tor(indok,:);
Rim     = equi.profiles_1d.gm9(indok,:);
Rm      = equi.profiles_1d.gm8(indok,:);
volume  = equi.profiles_1d.volume(indok,:);
area    = equi.profiles_1d.area(indok,:);
ephi    = - Rim .* z0dxdt(psi1d,tequi);
ej      = jphi_1d .* ephi;
pohm_equi1d = sum(diff(volume,1,2) .* (ej(:,1:end-1) + ej(:,2:end)) ./ 2,2);
ip_equi1d   = sum(diff(area,1,2) .* (jphi_1d(:,1:end-1) + jphi_1d(:,2:end)) ./ 2,2);
vres_equi1d = pohm_equi1d ./ ip_equi1d;

% computation of derived quantities
aminor = (equi.profiles_1d.r_outboard(indok,:) - equi.profiles_1d.r_inboard(indok,:)) ./ 2;
raxe   = (equi.profiles_1d.r_outboard(indok,:) + equi.profiles_1d.r_inboard(indok,:)) ./ 2;
epsi   = aminor ./ raxe;
x      = aminor ./ max(eps,aminor(:,end) * ones(1,size(aminor,2)));
%
Bt     = abs(equi.vacuum_toroidal_field.r0 .* equi.vacuum_toroidal_field.b0(indok)' ./ raxe(:,end));
rmx    = sqrt(abs(equi.profiles_1d.phi(indok,:)) ./ pi ./ (Bt * ones(1,size(x,2))));
% dpsidx s'annule au centre
psid1    = pdederive(x,psi1d,0,2,2,1);
% dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
psid2    = pdederive(x,psi1d,1,0,2,2);
psid2(:,1) = - Bt ./ abs(equi.profiles_1d.q(indok,1)) .* rmx(:,end) .^ 2;
%
bpol  = sqrt(equi.profiles_1d.gm2(indok,:) .* equi.profiles_1d.dpsi_drho_tor(indok,:) .^ 2) ./ 2 ./ pi;


% create external data structure 
data.x             = x;
data.time          = tequi;
data.psi           = psi1d ;
data.phi           = abs(equi.profiles_1d.phi(indok,:));
%data.phi_tor       = data.phi;   % to be checked !
data.dphidt        = z0dxdt(data.phi,data.time);
data.psid1         = psid1;
data.psid2         = psid2;
data.dpsidt        = z0dxdt(data.psi,data.time);
data.qjli          = abs(equi.profiles_1d.q(indok,:));
data.jli           = jphi_1d;
data.epar          = ephi;
data.bpol          = bpol;
data.rmx           = rmx;
data.rm            = max(rmx,[],2);
data.drmdt         = z0dxdt(data.rm,data.time);
data.grho2r2       = equi.profiles_1d.gm2(indok,:);
data.r2i           = equi.profiles_1d.gm1(indok,:);
data.ri            = equi.profiles_1d.gm9(indok,:);
data.grho          = abs(equi.profiles_1d.gm7(indok,:));
data.grho2         = abs(equi.profiles_1d.gm3(indok,:));
data.spr_tor       = equi.profiles_1d.darea_drho_tor(indok,:);
if all(data.spr_tor(:) == 0)
    data.spr_tor      =  equi.profiles_1d.gm9(indok,:) .* equi.profiles_1d.dvolume_drho_tor(indok,:) ./ (2 * pi);
end
data.vpr_tor       = equi.profiles_1d.dvolume_drho_tor(indok,:); 
data.volume        = equi.profiles_1d.volume(indok,:);
data.surface_pol   = equi.profiles_1d.area(indok,:);
data.surface_flux  = equi.profiles_1d.surface(indok,:);
data.dphidx        = pdederive(x,data.phi,0,2,2,1);
data.vpr           = pdederive(x,data.volume,0,2,2,1);
data.spr           = pdederive(x,data.surface_pol,0,2,2,1);
data.C2            = equi.profiles_1d.dvolume_drho_tor(indok,:) .* data.grho2r2;
data.C3            = equi.profiles_1d.dvolume_drho_tor(indok,:) .* data.r2i;
data.ej            = ej;
data.jeff          = equi.profiles_1d.j_parallel(indok,:);
data.jres          = NaN * data.jeff; % must be computed with METIS data
data.Raxe          = raxe;
data.epsi          = epsi;
data.a             = aminor;
data.difcurconv    = equi.convergence.iterations_n(indok)';
%data.jgs           = data.jli;
data.df2dpsi       = - (4*pi) .* equi.profiles_1d.f_df_dpsi(indok,:);
data.dptotdpsi     = - (2*pi) .*  equi.profiles_1d.dpressure_dpsi(indok,:);
data.ptot          = abs(equi.profiles_1d.pressure(indok,:));
data.fdia          = abs(equi.profiles_1d.f(indok,:));
data.wbp           = mu0 .* equi.global_quantities.ip(indok)' .^ 2 .* raxe(:,end) ./ 4 .* equi.global_quantities.li_3(indok)';
data.dwbpdt        = z0dxdt(data.wbp,data.time);
data.kx            = equi.profiles_1d.elongation(indok,:);
data.dx            = (equi.profiles_1d.triangularity_lower(indok,:) + equi.profiles_1d.triangularity_upper(indok,:)) ./ 2;
data.ipout         = abs(equi.global_quantities.ip(indok)');
data.lif           = equi.global_quantities.li_3(indok)';
data.betap         = equi.global_quantities.beta_pol(indok)';
data.pohm          = pohm_equi1d;
%%data.b0             = Bt;
data.q95           = abs(equi.global_quantities.q_95(indok)');
data.qmin          = min(data.qjli,[],2);
data.q0            = abs(equi.global_quantities.q_axis(indok)');
data.psin          = (data.psi - data.psi(:,1) * ones(1,size(data.psi,2))) ./ ...
                     ((data.psi(:,end) - data.psi(:,1)) * ones(1,size(data.psi,2)));
dd95               = abs(data.psin - 0.95);
mask95             = double(dd95 == min(dd95,[],2) *  ones(1,size(data.psin,2)));
data.d95           = sum(data.dx .* mask95,2) ./ max(1, sum(mask95,2));                
data.K95           = sum(data.kx .* mask95,2) ./ max(1, sum(mask95,2));   
data.piqj          = max(0.1,min(10,data.jli(:,1) ./ (sum(diff(data.surface_pol,1,2) .*  ...
                     (data.jli(:,1:end-1) + data.jli(:,2:end)) ./ 2,2) ./ data.surface_pol(:,end))));
                 
% LCFS
data.geo.a         = equi.boundary.minor_radius(indok)';
data.geo.R         = equi.boundary.geometric_axis.r(indok)';
data.geo.z0        = equi.boundary.geometric_axis.z(indok)';
data.geo.K         = equi.boundary.elongation(indok)';
data.geo.d         = equi.boundary.triangularity(indok)';
data.geo.b0        = abs(equi.vacuum_toroidal_field.r0 .* equi.vacuum_toroidal_field.b0(indok)' ./ data.geo.R);
data.geo.vp        = equi.profiles_1d.volume(indok,end);
data.geo.sp        = equi.profiles_1d.area(indok,end);
data.geo.peri      = equi.profiles_1d.surface(indok,end) ./ (2*pi*data.geo.R);
data.geo.sext      = equi.profiles_1d.surface(indok,end);
data.geo.xpoint    = double(isfinite(equi.boundary.x_point.r(indok)') & isfinite(equi.boundary.x_point.z(indok)'));
data.geo.Rsepa     = equi.boundary.outline.r(indok);
data.geo.Zsepa     = equi.boundary.outline.z(indok);

for k = 1:length(data.geo.Rsepa)
        rep = z0polygeom(data.geo.Rsepa{k},data.geo.Zsepa{k});
        data.geo.peri(k)      = rep(4);
end

% inversion radius
mask = (data.qjli <= 1);
ggr  = ones(length(data.time),1) * (1:size(data.qjli,2));
data.indice_inv = max(mask .* ggr,[],2);

% Poynting flux
% ref :Ohmic ???ux consumption during initial operation of the NSTX spherical torus
% J.E. Menard et al, NF 2001 41 1197
% ici on prend en compte que la composante est utile pour la consommation de flux poloidal.
data.poynting = NaN * ones(size(data.time));
dWbpdt        = z0dxdt(data.bpol .^ 2 ./ 2 ./ mu0,data.time);
for k=1:length(data.time)
    data.poynting(k) = data.rm(k) .* trapz(x(k,:),data.vpr_tor(k,:) .* (dWbpdt(k,:) + data.ej(k,:)),2);
end

% effective trapped fraction
if  ~isempty(equi.profiles_1d.trapped_fraction)
    data.ftrap = equi.profiles_1d.trapped_fraction(indok,:);
else
    data.ftrap = ftrap(data.x,data.epsi,data.psi,data.bpol,data.fdia,data.Raxe,data.r2i,data.ri,data.grho2r2,data.grho2,data.a);
end

% set same stucture for current diffusion and eqauilibrium
setappdata(0,'CURDIFF_EXP',data);
setappdata(0,'EQUILIBRIUM_EXP',data);


function ft = ftrap(x,epsi,psi,bpol,fdia,Raxe,r2i,ri,grho2r2,grho2,a)

b2m       = bpol .^ 2  + fdia .^ 2 .* r2i;
bm        = fdia .* ri .* (1 + 0.5 .* bpol ./ sqrt(grho2r2) ./ fdia .* sqrt(grho2));
% calcul de bmax
btor         = (fdia ./ (Raxe - a));
grho         = abs(1 ./ max(eps,abs(pdederive(x,Raxe - a,0,2,2,1))));
grho(1)      = 2 .* grho(2) - grho(3);
bpol         = abs(pdederive(x,psi,0,2,2,1))./ (Raxe - a) .* grho;
bmax         = sqrt(btor .^ 2 + bpol .^ 2);

% variable
h  = min(1,bm ./ bmax);
h2 = min(1,b2m ./ bmax .^ 2);

% expression de ftrap  Lin-Liu and Miller Phys. of Plasmas 2 (5) 1995
ftu          = 1 - h2 ./ h .^ 2 .* ( 1 - sqrt(1 - h) .* (1 + 0.5 .* h));
ftl          = 1 - h2 .* (1./ h2 - sqrt(1-sqrt(h2)) ./ h2  - sqrt(1 - sqrt(h2)) ./ 2 ./h);
ftc          = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));  % Wesson , Tokamak
fte          = 0.75 .* ftu + 0.25 .* ftl;
ft           = x .* fte + (1 - x) .* ftc;

ft_alt     = 1 - (1-epsi) .^ 2 ./ sqrt(1 - epsi .^ 2) ./ (1 + 1.46 .* sqrt(epsi));
ft(~isfinite(ft)) = ft_alt(~isfinite(ft));

