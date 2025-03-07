%% FILL CORE_PROFILE IDS WITH ONLY INTERESTING VARIABLES
function [edge_profiles,edge_transport] = mapedge_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,pulse_schedule,sigma_B0_eff)


%% IMAS part start here
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

% backward compatibility
if ~isfield(z0dstruct.z0dinput.option,'Sn_fraction')
    z0dstruct.z0dinput.option.Sn_fraction = 0;
end

% isotopic composition for option.gaz == 5
if z0dstruct.z0dinput.option.gaz == 5
    nHe3onD = real(z0dstruct.z0dinput.cons.iso);
    nTonD   = imag(z0dstruct.z0dinput.cons.iso);
    warning('nHe3onD & nTonD not yet used !');
else
    nHe3onD = zeros(size(z0dstruct.z0dinput.cons.iso));
    nTonD   = real(z0dstruct.z0dinput.cons.iso);
end
z0dstruct.z0dinput.cons.iso = real(z0dstruct.z0dinput.cons.iso);

% precomputation
% donnees du modele a 2 points
zs      = data_zerod;
option  = z0dstruct.z0dinput.option;
geo     = z0dstruct.z0dinput.geo;
cons    = z0dstruct.z0dinput.cons;
profli  = profil0d;
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
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot_in,pl_in,zeff_div,gamma,mach_target,prad_loc,pradsol_loc,fcond,profli] = ...
    z0convergence_2points_dic(option,cons,geo,zs,profli);
option = option_mem;
qpl_rad_sol = option.fpower .* pradsol ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_tot     = option.fpower .* pl ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
qpl_in_max  = option.fpower .* zs.pin ./ (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut));
%(4*pi*rp*sol_width/q95)
%qpar=p_sep*1e6/(4*pi*rp*sol_width/q95); % W/m2 estimate of parallel q
%L=pi*q95*rp;
qpl_target  = qpl_tot -  qpl_rad_sol - qpl_rad_div - qpl_neutral_div;
praddiv = qpl_neutral_div .* (2 .* pi  .* Raxea(:,end) .* zs.dsol .* sin(ut)) ./ option.fpower;

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
    fmom_corr = max(0.1,min(1,0.17 + 2 .* (alpha ./ (alpha + 1)) .^ ((alpha + 1) ./ 2)));
else
    fmom_corr = option.fmom .* ones(size(zs.telim));
end


% flux de matiere dans la sol
switch option.gaz
    case 1
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 1);
    case 2
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 2);
    case 3
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
            (2 + 3.* cons.iso) .* (1 + cons.iso));
    case 4
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ 4);
        
    case 5
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
            (2 + 3.* cons.iso) .* (1 + cons.iso));
    case 11
        cs0 = sqrt(1.602176462e-19 .* (zs.tite + zs.zeff)./1.6726485e-27 ./ ...
              (1 + 11 .* cons.iso) .* (1 + cons.iso));
 
end
flux_target = mach_target .* cs0 .* sqrt(zs.telim) .* zs.nelim;
neu =  (4 .* zs.telim .* zs.nelim) ./ (fmom .* fie) ./ zs.tebord;

if option.mach_corr == 1
    ftm = (gamma .* (1 + mach_target .^ 2) ./ (2 .* abs(mach_target) .* (gamma - 1 + mach_target .^ 2))) .^ 2;
    fnm = 2 ./ (ftm .* (1 + mach_target .^ 2));
    neu = neu ./ ftm ./ fnm;
    flux_target = flux_target .* mach_target;
end
flux_output = interp1(profli.temps,profli.ge(:,end) .* profli.grho2(:,end) .* profli.vpr_tor(:,end),zs.temps,'pchip','extrap');

maskx = ones(size(zs.temps));
maskx(zs.xpoint == 0) = NaN;

% memory allocation
edge_profiles  = ids_gen('edge_profiles');
%edge_sources   = ids_gen('edge_sources');
edge_transport = ids_gen('edge_transport');

% copy of common data with summary
if ~isempty(pulse_schedule)
    edge_profiles.code   = pulse_schedule.code;
    %edge_sources.code    = summary.code;
    edge_transport.code  = pulse_schedule.code;
end
edge_profiles.ids_properties.comment           = 'METIS edge data from 2 points model';
%edge_sources.ids_properties.comment            = 'METIS edge data from 2 points model';
edge_transport.ids_properties.comment          = 'METIS edge data from 2 points model';
edge_profiles.ids_properties.homogeneous_time  = 1;
%edge_sources.ids_properties.homogeneous_time   = 1;
edge_transport.ids_properties.homogeneous_time = 1;

% models description
edge_transport.model{2} = edge_transport.model{1};
edge_transport.model{1}.identifier.name        = 'Heat flux';
edge_transport.model{1}.identifier.index       = 1;
edge_transport.model{1}.identifier.description = 'heat flux from 2 points model (point 1 = LCFS and point 2 = target/limiter)';
edge_transport.model{2}.identifier.name        = 'SOL and divertor radiative and neutral losses';
edge_transport.model{2}.identifier.index       = 2;
edge_transport.model{2}.identifier.description = 'SOL and divertor radiative and neutral  losses from 2 points model (point 1 = SOL or egde losses and point 2 = near target/limiter losses; electrons fluxes are radiative losses and neutral is neutral friction losses)';

% time vectors
edge_profiles.time    = data_zerod.temps;
edge_transport.time   = data_zerod.temps;

% sub structure initialisation
for k=1:length(data_zerod.temps)
    edge_profiles.profiles_1d{k}   = edge_profiles.profiles_1d{1};
    edge_profiles.grid_ggd{k}      = edge_profiles.grid_ggd{1};
    edge_profiles.ggd{k}           = edge_profiles.ggd{1};
    edge_profiles.ggd_fast{k}      = edge_profiles.ggd_fast{1};
    edge_transport.grid_ggd{k}     = edge_transport.grid_ggd{1};
    edge_transport.model{1}.ggd{k} = edge_transport.model{1}.ggd{1};
    edge_transport.model{2}.ggd{k} = edge_transport.model{1}.ggd{1};
    edge_transport.model{1}.ggd_fast{k} = edge_transport.model{1}.ggd_fast{1};
    edge_transport.model{2}.ggd_fast{k} = edge_transport.model{1}.ggd_fast{1};
end
ggd         = edge_transport.grid_ggd{1};
if isfield(ggd,'grid')
    grid        = ggd.grid;
    grid_subset = grid.grid_subset{1};
else
    grid        = ggd;
    grid_subset = grid.grid_subset{1};
end
element     = grid_subset.element{1};
object      = element.object{1};


% start to fill data fields
% boucle sur les temps
for k=1:length(data_zerod.temps)
    % common grid
    % generation of 2 points grid for one time slice
    grid.identifier.name            = 'twopoints';
    grid.identifier.index           = 1;
    grid.identifier.description     = 'grid for two points model';
    space.geometry_type.name        = 'magnetic field line points';
    space.geometry_type.index       = 1;
    space.geometry_type.description = 'magnetic field line points';
    space.coordinates_type          = 0;
    %boundary.index                  =
    %boundary.neighbours             =
    %object.boundary{1}              = boundary;
    object.geometry                 = 0;
    object.nodes                    = 1;
    object.measure                  = 0;
    objects_per_dimension.object{1} = object;
    %boundary.index                  =
    %boundary.neighbours             =
    %object.boundary{1}              = boundary;
    object.geometry                 = 0;
    object.nodes                    = 2;
    object.measure                  = lc(k);
    objects_per_dimension.object{2} = object;
    space.objects_per_dimension{1}  = objects_per_dimension;
    grid.space{1}                   = space;
    grid_subset.identifier.name     = 'magnetic field line';
    grid_subset.identifier.index    = 1;
    grid_subset.identifier.description ='magnetic field line';
    grid_subset.dimension              = 1;
    object.space                       = 1;
    object.dimension                   = 1;
    object.index                       = 1;
    element.object{1}                  = object;
    object.space                       = 1;
    object.dimension                   = 1;
    object.index                       = 2;
    element.object{2}                  = object;
    grid_subset.element{1}             = element;
    %base.jacobian         		   =
    %base.tensor_covariant              =
    %base.tensor_contravariant          =
    %grid_subset.base{1}                = base;
    %grid_subset.metric.jacobian        =
    %grid_subset.metric.tensor_covariant =
    %grid_subset.metric.tensor_contravariant =
    grid.grid_subset{1}             = grid_subset;
    % switch  on version of imas
    if test_imas_ggd_version_new
        ggd.time =  data_zerod.temps(k);
        edge_profiles.grid_ggd{k}  = ggd;
        %edge_sources.grid_ggd{k}   = ggd;
        edge_transport.grid_ggd{k} = ggd;
        %
        edge_profiles.ggd{k}.electrons = edge_profiles.ggd{1}.electrons;
        edge_profiles.ggd{k}.electrons.temperature{1} = edge_profiles.ggd{1}.electrons.temperature{1};
        edge_profiles.ggd{k}.electrons.temperature{2} = edge_profiles.ggd{1}.electrons.temperature{1};
        edge_profiles.ggd{k}.electrons.temperature{1}.values = zs.tebord(k);
        edge_profiles.ggd{k}.electrons.temperature{1}.grid_index = 1;
        edge_profiles.ggd{k}.electrons.temperature{1}.grid_subset_index = 1;
        edge_profiles.ggd{k}.electrons.temperature{2}.values = zs.telim(k);
        edge_profiles.ggd{k}.electrons.temperature{2}.grid_index = 2;
        edge_profiles.ggd{k}.electrons.temperature{2}.grid_subset_index = 1;
        edge_profiles.ggd{k}.electrons.density{1} = edge_profiles.ggd{1}.electrons.density{1};
        edge_profiles.ggd{k}.electrons.density{2} = edge_profiles.ggd{1}.electrons.density{1};
        edge_profiles.ggd{k}.electrons.density{1}.values = zs.nebord(k);
        edge_profiles.ggd{k}.electrons.density{1}.grid_index = 1;
        edge_profiles.ggd{k}.electrons.density{1}.grid_subset_index = 1;
        edge_profiles.ggd{k}.electrons.density{2}.values = zs.nelim(k);
        edge_profiles.ggd{k}.electrons.density{2}.grid_index = 2;
        edge_profiles.ggd{k}.electrons.density{2}.grid_subset_index = 1;
        edge_profiles.ggd{k}.zeff{1} = edge_profiles.ggd{1}.zeff{1};
        edge_profiles.ggd{k}.zeff{2} = edge_profiles.ggd{1}.zeff{1};
        edge_profiles.ggd{k}.zeff{1}.values = profli.zeff(k,end);
        edge_profiles.ggd{k}.zeff{1}.grid_index = 1;
        edge_profiles.ggd{k}.zeff{1}.grid_subset_index = 1;
        edge_profiles.ggd{k}.zeff{2}.values = zeff_div(k);
        edge_profiles.ggd{k}.zeff{2}.grid_index = 2;
        edge_profiles.ggd{k}.zeff{2}.grid_subset_index = 1;
        edge_profiles.ggd{k}.time = data_zerod.temps(k);
        %
        edge_transport.model{1}.ggd{k}.electrons = edge_transport.model{1}.ggd{1}.electrons;
        edge_transport.model{2}.ggd{k}.electrons = edge_transport.model{1}.ggd{1}.electrons;
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1}.values    = qpl_tot(k);
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1}.grid_index = 1;
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1}.grid_subset_index = 1;
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2}.values    = qpl_target(k);
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2}.grid_index = 2;
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2}.grid_subset_index = 1;
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.particles.flux{1} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.particles.flux{2} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1}.values = flux_output(k);
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1}.grid_index = 1;
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1}.grid_subset_index = 1;
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2}.values = flux_target(k);
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2}.grid_index = 2;
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2}.grid_subset_index = 1;
        edge_transport.model{1}.ggd{k}.time = data_zerod.temps(k);
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1}.values    = - qpl_rad_sol(k);
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1}.grid_index = 1;
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1}.grid_subset_index = 1;
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2}.values    = - qpl_target(k);
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2}.grid_index = 2;
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2}.grid_subset_index = 1;
        edge_transport.model{1}.ggd{k}.neutral{1} = edge_transport.model{1}.ggd{1}.neutral{1};
        edge_transport.model{2}.ggd{k}.neutral{1} = edge_transport.model{1}.ggd{1}.neutral{1};
        edge_transport.model{1}.ggd{k}.neutral{1}.energy.flux{1} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{1}.ggd{k}.neutral{1}.energy.flux{2} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1}.values   = -9.0000e+40;
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1}.grid_index = 1;
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1}.grid_subset_index = 1;
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2}.values   = - qpl_neutral_div(k);
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2}.grid_index = 2;
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2}.grid_subset_index = 1;
        edge_transport.model{2}.ggd{k}.time = data_zerod.temps(k);
        
        edge_transport.model{1}.ggd_fast{k}.power{1} = edge_transport.model{1}.ggd_fast{1}.power{1};
        edge_transport.model{1}.ggd_fast{k}.power{2} = edge_transport.model{1}.ggd_fast{1}.power{1};
        edge_transport.model{2}.ggd_fast{k}.power{1} = edge_transport.model{1}.ggd_fast{1}.power{1};
        edge_transport.model{2}.ggd_fast{k}.power{2} = edge_transport.model{1}.ggd_fast{1}.power{1};
        edge_transport.model{1}.ggd_fast{k}.power{1}.value    = - pradsol(k);
        edge_transport.model{1}.ggd_fast{k}.power{1}.grid_index = 1;
        edge_transport.model{1}.ggd_fast{k}.power{1}.grid_subset_index = 1;
        
        edge_transport.model{1}.ggd_fast{k}.power{2}.value    =  -9.0000e+40;
        edge_transport.model{1}.ggd_fast{k}.power{2}.grid_index = 2;
        edge_transport.model{1}.ggd_fast{k}.power{2}.grid_subset_index = 1;
        
        edge_transport.model{2}.ggd_fast{k}.power{1}.value    =  -9.0000e+40;
        edge_transport.model{2}.ggd_fast{k}.power{1}.grid_index = 1;
        edge_transport.model{2}.ggd_fast{k}.power{1}.grid_subset_index = 1;
        
        edge_transport.model{2}.ggd_fast{k}.power{2}.value    = - praddiv(k);
        edge_transport.model{2}.ggd_fast{k}.power{2}.grid_index = 2;
        edge_transport.model{2}.ggd_fast{k}.power{2}.grid_subset_index = 1;
        
        
    else
        ggd.grid              = grid;
        edge_profiles.ggd{k}  = ggd;
        %edge_sources.ggd{k}   = ggd;
        edge_transport.model{1}.ggd{k} = ggd;
        edge_transport.model{2}.ggd{k} = ggd;
        %
        edge_profiles.ggd{k}.electrons.temperature{1} = edge_profiles.ggd{1}.electrons.temperature{1};
        edge_profiles.ggd{k}.electrons.temperature{2} = edge_profiles.ggd{1}.electrons.temperature{1};
        edge_profiles.ggd{k}.electrons.temperature{1}.values = zs.tebord(k);
        edge_profiles.ggd{k}.electrons.temperature{2}.values = zs.telim(k);
        edge_profiles.ggd{k}.electrons.density{1} = edge_profiles.ggd{1}.electrons.density{1};
        edge_profiles.ggd{k}.electrons.density{2} = edge_profiles.ggd{1}.electrons.density{1};
        edge_profiles.ggd{k}.electrons.density{1}.values = zs.nebord(k);
        edge_profiles.ggd{k}.electrons.density{2}.values = zs.nelim(k);
        edge_profiles.ggd{k}.zeff{1} = edge_profiles.ggd{1}.zeff{1};
        edge_profiles.ggd{k}.zeff{2} = edge_profiles.ggd{1}.zeff{1};
        edge_profiles.ggd{k}.zeff{1}.values = profli.zeff(k,end);
        edge_profiles.ggd{k}.zeff{2}.values = zeff_div(k);
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1}; 
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1}; 
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1}; 
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2} = edge_transport.model{1}.ggd{1}.electrons.energy.flux{1}; 
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{1}.values    = qpl_tot(k);
        edge_transport.model{1}.ggd{k}.electrons.energy.flux{2}.values    = qpl_target(k);
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.particles.flux{1} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{2}.ggd{k}.electrons.particles.flux{2} = edge_transport.model{1}.ggd{1}.electrons.particles.flux{1};
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{1}.values = flux_output(k);
        edge_transport.model{1}.ggd{k}.electrons.particles.flux{2}.values = flux_target(k);
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{1}.values    = - qpl_rad_sol(k);
        edge_transport.model{2}.ggd{k}.electrons.energy.flux{2}.values    = - qpl_target(k);
        edge_transport.model{1}.ggd{k}.neutral{1} = edge_transport.model{1}.ggd{1}.neutral{1};
        edge_transport.model{2}.ggd{k}.neutral{1} = edge_transport.model{1}.ggd{1}.neutral{1};
        edge_transport.model{1}.ggd{k}.neutral{1}.energy.flux{1} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{1}.ggd{k}.neutral{1}.energy.flux{2} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2} = edge_transport.model{1}.ggd{1}.neutral{1}.energy.flux{1};
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{1}.values   = -9.0000e+40;
        edge_transport.model{2}.ggd{k}.neutral{1}.energy.flux{2}.values   = - qpl_neutral_div(k);
    end
end


function rep = test_imas_ggd_version_new

ep = struct(ids_gen('edge_profiles'));

if isfield(ep,'grid_ggd') && isfield(ep,'ggd_fast')
    rep = true;
else
    rep = false;
end


