% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function equilibrium = mapequilibrium_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,equilibrium,summary,opt_grille,sigma_B0_eff,equi_extrap,factor_two_pi)

% pour profile2d : opt_grille = 0 => grille en (psi,theta)
%                  opt_grille = 1 => Delaunay
if nargin < 8
  opt_grille = 0;
elseif isempty(opt_grille)
  opt_grille = 0;
end 
if nargin < 9
  sigma_B0_eff = 1;
elseif isempty(sigma_B0_eff)
  sigma_B0_eff = 1;
end 
if nargin < 10
  equi_extrap = 0;
elseif isempty(equi_extrap)
  equi_extrap = 0;
end 
FEEQS_post = false;
alternative_extrapolation = 'Interpolation';
if equi_extrap >= 2
    if equi_extrap == 3
        alternative_extrapolation = 'G-S polynomial';
    end
    FEEQS_post = true;
    equi_extrap = 0;
end

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

%% CONSTANT = NUMBER OF POLOIDAL POINTS
if isfield(z0dstruct.z0dinput.option,'nb_points_pol')
    nbth = z0dstruct.z0dinput.option.nb_points_pol;
else
    nbth = 65;
end

%% NUMBER OF GRID POINTS
if isfield(z0dstruct.z0dinput.option,'nb_points_radial')
    nbeqdsk = z0dstruct.z0dinput.option.nb_points_radial;
else
    nbeqdsk = 51;
end
%% SUMMARY HERITAGE
%  if isfield(summary,'ids_properties')
%    if isfield(summary.ids_properties,'comment')
%      equilibrium.ids_properties.comment = summary.ids_properties.comment;
%    
%    end
%  end
equilibrium.ids_properties.comment = 'METIS equilibrium';


%% GENERAL INFORMATION
equilibrium.code.name = 'METIS';

%% TIME
ntime = length(profil0d.temps);

%% adaptation of output_flag
equilibrium.code.output_flag = zeros(size(profil0d.temps));

for itime=2:ntime
  equilibrium.time_slice{itime} = equilibrium.time_slice{1};
end
for itime=1:ntime
  equilibrium.time_slice{itime}.time = profil0d.temps(itime);
end
equilibrium.time = profil0d.temps;


%% BOUNDARY STRUCTURE
xpoint = interp1_imas(data_zerod.temps,data_zerod.xpoint,profil0d.temps,'nearest','extrap');
for itime=1:ntime

  equilibrium.time_slice{itime}.boundary.type = xpoint(itime);   
  z0 = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.z0,profil0d.temps(itime),'nearest','extrap');
  if isfield(profil0d,'Rsepa') && isfield(profil0d,'Zsepa') && all(isfinite(profil0d.Rsepa(itime,:))) && all(isfinite(profil0d.Zsepa(itime,:)))
      % LFCS must be open in IMAS
      if (profil0d.Rsepa(itime,end) == profil0d.Rsepa(itime,1)) &&  (profil0d.Zsepa(itime,end) == profil0d.Zsepa(itime,1))
          equilibrium.time_slice{itime}.boundary.outline.r	= profil0d.Rsepa(itime,1:end-1);
          equilibrium.time_slice{itime}.boundary.outline.z	= profil0d.Zsepa(itime,1:end-1) + z0;
      else
          equilibrium.time_slice{itime}.boundary.outline.r	= profil0d.Rsepa(itime,:);
          equilibrium.time_slice{itime}.boundary.outline.z	= profil0d.Zsepa(itime,:) + z0;
      end
      % deduire de la sepa
      control =sepa_moments(profil0d.Rsepa(itime,:),profil0d.Zsepa(itime,:) + z0,equilibrium.time_slice{1}.boundary.x_point{1});
      equilibrium.time_slice{itime}.boundary.minor_radius           = control.a;
      equilibrium.time_slice{itime}.boundary.elongation             = control.K;
      equilibrium.time_slice{itime}.boundary.triangularity_upper    = control.du;
      equilibrium.time_slice{itime}.boundary.triangularity_lower    = control.dl;
      equilibrium.time_slice{itime}.boundary.elongation_upper       = control.Ku;
      equilibrium.time_slice{itime}.boundary.elongation_lower       = control.Kl;
      equilibrium.time_slice{itime}.boundary.squareness_lower_inner = control.squareness_lower_inner;
      equilibrium.time_slice{itime}.boundary.squareness_upper_inner = control.squareness_upper_inner;
      equilibrium.time_slice{itime}.boundary.squareness_lower_outer = control.squareness_lower_outer;
      equilibrium.time_slice{itime}.boundary.squareness_upper_outer = control.squareness_upper_outer;
      if isfield(control,'x_point')
          equilibrium.time_slice{itime}.boundary.x_point             = control.x_point;
      end
      equilibrium.time_slice{itime}.boundary.geometric_axis.r    = profil0d.Raxe(itime,end);
      equilibrium.time_slice{itime}.boundary.geometric_axis.z    = control.z0_geo;
      z0_mag = control.z0_geo;
      outlineisfilled = 1;
  else
      control = [];
      equilibrium.time_slice{itime}.boundary.outline.r	= [];
      equilibrium.time_slice{itime}.boundary.outline.z	= [];
      equilibrium.time_slice{itime}.boundary.minor_radius        = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a,profil0d.temps(itime),'nearest','extrap');
      equilibrium.time_slice{itime}.boundary.elongation          = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.K,profil0d.temps(itime),'nearest','extrap');
      equilibrium.time_slice{itime}.boundary.triangularity_upper = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.d,profil0d.temps(itime),'nearest','extrap');
      equilibrium.time_slice{itime}.boundary.triangularity_lower = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.d,profil0d.temps(itime),'nearest','extrap');
      equilibrium.time_slice{itime}.boundary.geometric_axis.r    = profil0d.Raxe(itime,end);
      equilibrium.time_slice{itime}.boundary.geometric_axis.z    = z0;
      z0_mag = z0;
      outlineisfilled = 0;
  end
  equilibrium.time_slice{itime}.boundary.psi = profil0d.psi(itime,end);
  noms_loc = fieldnames( equilibrium.time_slice{itime}.boundary);
  for lkz = 1:length(noms_loc)
      equilibrium.time_slice{itime}.boundary_separatrix.(noms_loc{lkz}) = equilibrium.time_slice{itime}.boundary.(noms_loc{lkz});
  end
  %% SCALAR PARAMETERS
  %% global_quantities
  equilibrium.time_slice{itime}.global_quantities.magnetic_axis.r = profil0d.Raxe(itime,1);
  equilibrium.time_slice{itime}.global_quantities.magnetic_axis.z               = z0_mag;
  equilibrium.time_slice{itime}.global_quantities.magnetic_axis.b_field_tor     = profil0d.fdia(itime,1) ./ profil0d.Raxe(itime,1);
  equilibrium.time_slice{itime}.global_quantities.ip           = interp1_imas(data_zerod.temps,data_zerod.ip,profil0d.temps(itime),'pchip','extrap');
  %
  equilibrium.time_slice{itime}.global_quantities.beta_tor     = interp1_imas(data_zerod.temps,data_zerod.betan,profil0d.temps(itime),'pchip','extrap') .*  ...
      equilibrium.time_slice{itime}.global_quantities.ip ./ (profil0d.fdia(itime,end) ./ profil0d.Raxe(itime,end)) ./  ...
      equilibrium.time_slice{itime}.boundary.minor_radius ./ 1e6;
  equilibrium.time_slice{itime}.global_quantities.beta_normal  =  interp1_imas(data_zerod.temps,100 .* data_zerod.betan,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.li_3         = interp1_imas(data_zerod.temps,data_zerod.li,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.volume       = interp1_imas(data_zerod.temps,data_zerod.vp,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.area         = interp1_imas(data_zerod.temps,data_zerod.sp,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.surface      = interp1_imas(data_zerod.temps,data_zerod.sext,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.length_pol   = interp1_imas(data_zerod.temps,data_zerod.peri,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.global_quantities.psi_axis     = profil0d.psi(itime,1);
  equilibrium.time_slice{itime}.global_quantities.psi_boundary = profil0d.psi(itime,end);
    % not the same definition as in METIS
  % Poloidal beta. Defined as betap = 4 int(p dV) / [R_0 * mu_0 * Ip^2] {dynamic} [-]
  equilibrium.time_slice{itime}.global_quantities.beta_pol     = 4 .* interp1_imas(data_zerod.temps,data_zerod.betaptot,profil0d.temps(itime),'pchip','extrap') *...
                                                                 (profil0d.bpol(itime,end) .^ 2  ./ 2 ./ phys.mu0)  .* equilibrium.time_slice{itime}.global_quantities.volume ./...
                                                                 (phys.mu0 * equilibrium.time_slice{itime}.global_quantities.ip .^ 2 * profil0d.Raxe(itime,end));

  %
  %equilibrium.time_slice{itime}.global_quantities.magnetic_axis.b_field_tor = profil0d.fdia(itime,1) ./ profil0d.Raxe(itime,1);
  equilibrium.time_slice{itime}.global_quantities.q_axis          = profil0d.qjli(itime,1);
  %equilibrium.time_slice{itime}.global_quantities.magnetic_axis.r = interp1_imas(data_zerod.temps,summary.centre.Rmag.value,profil0d.temps(time),'nearest','extrap');
  %equilibrium.time_slice{itime}.global_quantities.magnetic_axis.z = interp1_imas(data_zerod.temps,summary.centre.Zmag.value,profil0d.temps(itime),'nearest','extrap');
  %
  %equilibrium.time_slice{itime}.global_quantities.q_95 = interp1_imas(data_zerod.temps,summary.local.pedestal.q.value,profil0d.temps(itime),'nearest','extrap');
  psin_loc = (profil0d.psi(itime,:) - profil0d.psi(itime,1)) ./ (profil0d.psi(itime,end) - profil0d.psi(itime,1));
  equilibrium.time_slice{itime}.global_quantities.q_95 = pchip(psin_loc,profil0d.qjli(itime,:),0.95);
  [equilibrium.time_slice{itime}.global_quantities.q_min.value,iqmin] = min(profil0d.qjli(itime,:),[],2);
  equilibrium.time_slice{itime}.global_quantities.q_min.rho_tor_norm = profil0d.rmx(itime,iqmin) ./ max(profil0d.rmx(itime,:));
  %
  rb0 =   interp1_imas(data_zerod.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
  r0  =   interp1_imas(data_zerod.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
  
  equilibrium.vacuum_toroidal_field.r0 = mean(z0dstruct.z0dinput.geo.R);
  equilibrium.vacuum_toroidal_field.b0(itime) = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0(itime) ./ equilibrium.vacuum_toroidal_field.r0;
  %
  equilibrium.time_slice{itime}.global_quantities.w_mhd = interp1_imas(data_zerod.temps,data_zerod.w,profil0d.temps(itime),'nearest','extrap');
    
  %% FOR THE MACH NUMBER! (GONDOR TONE)
  ve   = ones(1,size(profil0d.tip,2));
  meff = interp1_imas(data_zerod.temps,data_zerod.meff,profil0d.temps(itime),'nearest','extrap');
  
  %% ALLOCATE MEMORY SPACE
  v1 = NaN .* ones(size(profil0d.psi(1,:)));
  equilibrium.time_slice{itime}.profiles_1d.psi = v1;
  equilibrium.time_slice{itime}.profiles_1d.phi = v1;
  equilibrium.time_slice{itime}.profiles_1d.pressure = v1;
  equilibrium.time_slice{itime}.profiles_1d.f = v1;
  equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi = v1;
  equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi = v1;
  equilibrium.time_slice{itime}.profiles_1d.j_tor = v1;
  equilibrium.time_slice{itime}.profiles_1d.j_parallel = v1;
  equilibrium.time_slice{itime}.profiles_1d.q = v1;
  equilibrium.time_slice{itime}.profiles_1d.magnetic_shear = v1;
  equilibrium.time_slice{itime}.profiles_1d.r_inboard = v1;
  equilibrium.time_slice{itime}.profiles_1d.r_outboard = v1;
  equilibrium.time_slice{itime}.profiles_1d.rho_tor = v1;
  equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm = v1;
  equilibrium.time_slice{itime}.profiles_1d.dpsi_drho_tor = v1;
  equilibrium.time_slice{itime}.profiles_1d.elongation = v1;
  equilibrium.time_slice{itime}.profiles_1d.triangularity_upper = v1;
  equilibrium.time_slice{itime}.profiles_1d.triangularity_lower = v1;
  equilibrium.time_slice{itime}.profiles_1d.volume = v1;
  equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi = v1;
  equilibrium.time_slice{itime}.profiles_1d.dvolume_drho_tor = v1;
  equilibrium.time_slice{itime}.profiles_1d.area = v1;
  equilibrium.time_slice{itime}.profiles_1d.darea_dpsi = v1;
  equilibrium.time_slice{itime}.profiles_1d.darea_drho_tor = v1;
  equilibrium.time_slice{itime}.profiles_1d.surface	= v1;
  equilibrium.time_slice{itime}.profiles_1d.trapped_fraction = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm1 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm2 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm3 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm4 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm5 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm6 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm7 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm8 = v1;
  equilibrium.time_slice{itime}.profiles_1d.gm9 = v1;
  equilibrium.time_slice{itime}.profiles_1d.b_field_average = v1;
  equilibrium.time_slice{itime}.profiles_1d.b_field_min = v1;
  equilibrium.time_slice{itime}.profiles_1d.b_field_max = v1;
  %
  iprof = 1;
  % initialisation sub structure
  equilibrium.time_slice{itime}.profiles_2d{iprof} = equilibrium.time_slice{1}.profiles_2d{1};
  if opt_grille == 0
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim1 = []; 	
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim2 = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.connect = [];
    %equilibrium.time_slice{itime}.profiles_2d{iprof}.grid_type.description = char({'1','inverse','3','polar'});
    equilibrium.time_slice{itime}.profiles_2d{iprof}.psi = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.j_tor = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.j_parallel = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_r = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_z = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_r = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_z = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_tor = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.r = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.z = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.phi = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.theta = [];
    
  else
    v3  = NaN .* ones(size(profil0d.psi,2),nbth);
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim1 = []; 	
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim2 = [];
    equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.connect = [];
    %equilibrium.time_slice{itime}.profiles_2d{iprof}.grid_type.description = char({'3','irregular','4','unstruct'});
    equilibrium.time_slice{itime}.profiles_2d{iprof}.psi = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.j_tor = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.j_parallel = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_r = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_z = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_r = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_z = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_tor = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.r = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.z = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.phi = v3;
    equilibrium.time_slice{itime}.profiles_2d{iprof}.theta = v3;
    
  end

  %
  v3  = NaN .* ones(size(profil0d.psi,2),nbth);
  v3c = NaN .* ones(size(profil0d.psi,2),nbth,3,3);
  equilibrium.time_slice{itime}.coordinate_system.grid_type.name  = [];
  equilibrium.time_slice{itime}.coordinate_system.grid.dim1 = NaN .* ones(length(profil0d.xli));
  equilibrium.time_slice{itime}.coordinate_system.grid.dim2 = NaN .* ones(1,nbth);
  equilibrium.time_slice{itime}.coordinate_system.jacobian  = v3;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant = v3c;
  equilibrium.time_slice{itime}.coordinate_system.r = v3;
  equilibrium.time_slice{itime}.coordinate_system.z = v3;
  
  %end %% END LOOP ITIME=1:NTIME
  
  %% CALCULATION ON A RECTANGULAR GRID
  profiles_2d_2{1}.grid.dim1    = []; 	
  profiles_2d_2{1}.grid.dim2    = [];
  profiles_2d_2{1}.grid.connect = [];
  %profiles_2d_2{1}.grid_type.description = char({'1','rectangular',' ',' '});;
  profiles_2d_2{1}.psi	      = [];
  profiles_2d_2{1}.j_tor        = [];
  profiles_2d_2{1}.j_parallel   = [];
  profiles_2d_2{1}.b_r 	      = [];
  profiles_2d_2{1}.b_z 	      = [];
  profiles_2d_2{1}.b_field_r    = [];
  profiles_2d_2{1}.b_field_z    = [];
  profiles_2d_2{1}.b_field_tor  = [];
  profiles_2d_2{1}.r            = [];
  profiles_2d_2{1}.z            = [];      
  profiles_2d_2{1}.phi          = [];
  profiles_2d_2{1}.theta        = [];
  %profiles_2d_2.pressure        = [];

  %% THIS IS EXTREMELY SLOW IN MATLAB
  %% [equilibrium_cpo.coord_sys.position(1:size(profil0d.psi,1),1:size(profil0d.psi,2),1:nbth).R] = deal(NaN);
  %% [equilibrium_cpo.coord_sys.position(1:size(profil0d.psi,1),1:size(profil0d.psi,2),1:nbth).Z] = deal(NaN);
  %% THE RECTANGULAR GRID SHOULD HAVE A FIXED NUMBER OF POINTS
  rzp = max(1,max(z0dstruct.z0dinput.geo.K(:)));

  %% LOOP FOR CALCULATING MISSING OUTPUTS
  %for k =1:ntime

  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d/%d\t',itime,length(equilibrium.time_slice));
  
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
  geo_input.a     = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a ,profil0d.temps(itime),'pchip','extrap');
  geo_input.R     = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps(itime),'pchip','extrap');
  geo_input.K     = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.K ,profil0d.temps(itime),'pchip','extrap');
  geo_input.d     = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.d ,profil0d.temps(itime),'pchip','extrap');
  geo_input.b0    = z0dstruct.z0dinput.option.signe .* interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.b0 ,profil0d.temps(itime),'pchip','extrap'); 
  z0_offset       = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.z0 ,profil0d.temps(itime),'pchip','extrap');
  geo_input.sp    = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sp ,profil0d.temps(itime),'pchip','extrap');
  geo_input.vp    = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.vp ,profil0d.temps(itime),'pchip','extrap');
  geo_input.sext  = interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sext ,profil0d.temps(itime),'pchip','extrap');
  if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa') && all(isfinite(profil0d.Rsepa(itime,:))) && all(isfinite(profil0d.Zsepa(itime,:)))
      geo_input.Rsepa = profil0d.Rsepa(itime,:);
      geo_input.Zsepa = profil0d.Zsepa(itime,:);
      % safety rule: must be convex
      if isfield(z0dstruct.z0dinput.option,'Convex_LCFS') && (z0dstruct.z0dinput.option.Convex_LCFS == 1)
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
    box.a    = max(z0dstruct.z0dinput.geo.a);
    box.z0   = sum(z0dstruct.z0dinput.geo.z0 .* z0dstruct.z0dinput.cons.ip) ./ sum(z0dstruct.z0dinput.cons.ip);
    box.rmin = max(eps,min(profil0d.Rsepa(:)) - box.a/2);
    box.rmax = max(profil0d.Rsepa(:)) + box.a/2;
    box.zmin = min(profil0d.Zsepa(:)) - box.a/2*rzp + box.z0;
    box.zmax = max(profil0d.Zsepa(:)) + box.a/2*rzp + box.z0;
  else
    geo_input.Rsepa = [];
    geo_input.Zsepa = []; 
    z0dstruct.z0dinput.option.morphing = 0;
    geo_input.z0 = 0;
    box.a    = max(z0dstruct.z0dinput.geo.a);
    box.z0   = sum(z0dstruct.z0dinput.geo.z0 .* z0dstruct.z0dinput.cons.ip) ./ sum(z0dstruct.z0dinput.cons.ip);
    box.rmin = max(eps,min(z0dstruct.z0dinput.geo.R - z0dstruct.z0dinput.geo.a) - box.a/2);
    box.rmax = max(z0dstruct.z0dinput.geo.R + z0dstruct.z0dinput.geo.a) + box.a/2;
    box.zmin = min(- z0dstruct.z0dinput.geo.a .* z0dstruct.z0dinput.geo.K) - box.a/2*rzp + box.z0;
    box.zmax = max(z0dstruct.z0dinput.geo.a .* z0dstruct.z0dinput.geo.K) + box.a/2*rzp + box.z0;
  end
  % normalisation volume
  prof.volume  = prof.volume ./ prof.volume(end)  .* geo_input.vp;
  prof.surface = prof.surface ./ prof.surface(end) .* geo_input.sp;
  % derivative of psi
  % dpsidx s'annule au centre
  prof.psid1    = pdederive(prof.x,prof.psi,0,2,2,1);
  % dspidx = 0 au centre et d2psidx2 doit etre nul au bord pour que ip soit defini precisement
  prof.psid2    = pdederive(prof.x,prof.psi,1,0,2,2);
  if z0dstruct.z0dinput.option.cronos_regul == 4
	    prof.psid2(1) = - geo_input.b0 ./ prof.q(1) .* prof.rho(end) .^ 2;
  end
  %
  [profil,deuxd,moment,Ip_lcfs] = metis2equi1t_imas(prof,geo_input,phys,sign(factor_two_pi) .* nbth,z0dstruct.z0dinput.option.morphing,factor_two_pi);
  deuxd.Z = deuxd.Z + geo_input.z0;
  moment.zaxe = moment.zaxe   +  geo_input.z0;  
  Iprof = cumtrapz(profil0d.xli,profil0d.spr(itime,:) .* profil0d.jli(itime,:),2);
  Iprof_bis = profil0d.rmx(itime,end) .* cumtrapz(profil0d.xli,profil0d.vpr_tor(itime,:) .* profil0d.jli(itime,:) ./ 2 ./ pi .* profil0d.ri(itime,:),2);
%    figure(118);clf
%    plot(prof.x,Ip_lcfs,'b',profil0d.xli,Iprof,'r',1,equilibrium.time_slice{itime}.global_quantities.ip,'ok',profil0d.xli,Iprof_bis,'g');
%    drawnow
  delta_ip = equilibrium.time_slice{itime}.global_quantities.ip - Ip_lcfs(end);
  if delta_ip > 0
    equilibrium.time_slice{itime}.global_quantities.ip_error_upper = delta_ip;
    equilibrium.time_slice{itime}.global_quantities.ip_error_lower = 0;
  else
    equilibrium.time_slice{itime}.global_quantities.ip_error_upper = 0;
    equilibrium.time_slice{itime}.global_quantities.ip_error_lower = abs(delta_ip);
  end
  error_2d = abs(abs(equilibrium.time_slice{itime}.global_quantities.ip) - abs(Ip_lcfs(end))) ./ max(1,abs(equilibrium.time_slice{itime}.global_quantities.ip));
  fprintf(' error_2D = %g | ',error_2d);
  % rescale to preserve plasma current even with high pedestal
  Iprof  = abs(Iprof ./ Iprof(end) .* equilibrium.time_slice{itime}.global_quantities.ip);
  factor = abs(Iprof ./ max(1,abs(Ip_lcfs)));
  factor(1) = 1;
  %figure(119);clf;plot(factor);drawnow
  vth            = ones(1,size(deuxd.dPSIdx,2));
  deuxd.dPSIdx   = (factor' * vth) .* deuxd.dPSIdx;
  deuxd.BR       = (factor' * vth) .* deuxd.BR;
  deuxd.BZ       = (factor' * vth) .* deuxd.BZ;
  deuxd.dPSIdR   = (factor' * vth) .* deuxd.dPSIdR; 
  deuxd.dPSIdZ   = (factor' * vth) .* deuxd.dPSIdZ; 
  
%    figure(37)
%    clf
%    plot(deuxd.R',deuxd.Z' + z0_offset,'r');
%    hold on
%    if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa')
%      plot(profil0d.Rsepa(itime,:),profil0d.Zsepa(itime,:) + z0_offset ,'k');
%    end
%    xlabel('R (m)');
%    ylabel('Z (m)');
%    drawnow
%    % controle
%    noms = fieldnames(profil);
%    lz = 1;
%    figure(21);clf
%    for kz=1:length(noms)
%      if strfind(noms{kz},'_ctr')
%        subplot(4,4,lz)
%        plot(profil0d.xli,profil0d.(strrep(noms{kz},'_ctr',''))(itime,:),'r',profil.x,profil.(noms{kz}),'b')
%        title(noms{kz})
%        lz =lz +1;
%      end
%    end
%    subplot(4,4,lz)
%    plot(profil0d.xli,profil0d.bpol(itime,:),'r',profil.x,profil.bpolm,'b')
%    title('bpolm')
%    drawnow
  
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
  
  % remplissage des profils;
  equilibrium.time_slice{itime}.profiles_1d.psi = profil0d.psi(itime,:);
  switch z0dstruct.z0dinput.option.COCOS
      case {1,2,3,4,11,12,13,14}
           equilibrium.time_slice{itime}.profiles_1d.phi = profil.phi .* sign(mean(profil.phi)) .* sign(profil0d.psi(itime,end) - profil0d.psi(itime,1)) .* sign(mean(profil0d.qjli(itime,:)));
    otherwise
          equilibrium.time_slice{itime}.profiles_1d.phi = profil.phi .* sign(mean(profil.phi)) .* sign(profil0d.psi(itime,1) - profil0d.psi(itime,end)) .* sign(mean(profil0d.qjli(itime,:)));
  end
  equilibrium.time_slice{itime}.profiles_1d.pressure = profil0d.ptot(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.f = mean(sign(profil0d.fdia(itime,:))) .* profil0d.fdia(itime,:) .* sign(equilibrium.vacuum_toroidal_field.b0(itime));
  
  % new computation with improved precision on magnetic axis
  %equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi = pdederive(profil0d.psi(itime,:),profil0d.ptot(itime,:),0,2,2,1);
  %equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi = pdederive(profil0d.psi(itime,:),profil0d.fdia(itime,:) .^ 2,0,2,2,1) ./ 2 ;
  mu0 = 4e-7 .* pi;
  psid1_             = pdederive(profil0d.xli,profil0d.psi(itime,:),0,2,2,1);
  psid1_(psid1_ == 0) = sign(mean(psid1_)) .* eps;
  psid1_(end)      = -(2*pi) .* mu0 .* profil0d.rmx(itime,end) .* equilibrium.time_slice{itime}.global_quantities.ip  ./ profil0d.C2(itime,end) .* factor_two_pi;
  psid2_    = pdederive(profil0d.xli,profil0d.psi(itime,:),1,0,2,2);
  dpdpsi_   = pdederive(profil0d.xli,profil0d.ptot(itime,:),0,1,2,1) ./ psid1_;
  dpdpsi_(abs(psid1_)  == eps) = 0;
  sm = sign(mean(dpdpsi_));
  dpdpsi_(1) = double(sm >= 0) .* max(0,dpdpsi_(1)) + double(sm < 0).* min(0,dpdpsi_(1)); 
  dpdx_1 = (profil0d.ptot(itime,2) - profil0d.ptot(itime,1)) ./ profil0d.xli(2) .^ 2;
  dpdpsi_1 = sm .* abs(2 .* dpdx_1 ./ (profil0d.fdia(itime,1) ./ profil0d.Raxe(itime,1)) .* profil0d.qjli(itime,1) ./ profil0d.rmx(itime,end) .^ 2 ./ factor_two_pi);
  dpdpsi_(1) = double(sm >= 0) .* max(dpdpsi_(1),dpdpsi_1) + double(sm < 0) .* min(dpdpsi_(1),dpdpsi_1);
  df2dpsi_  = 2 .* mu0 .* (max(0,sign(equilibrium.time_slice{itime}.global_quantities.ip) .* profil0d.jli(itime,:)) .* profil0d.ri(itime,:) - ...
              sign(equilibrium.time_slice{itime}.global_quantities.ip) .* dpdpsi_ * factor_two_pi)./ profil0d.r2i(itime,:) / factor_two_pi;
  equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi = dpdpsi_;
  equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi = df2dpsi_ ./ 2;
  
%  figure(21);clf;plot( equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi,'r', ...
%  profil0d.psi(itime,:),pdederive(profil0d.psi(itime,:),profil0d.ptot(itime,:),0,2,2,1),'sm', ...
%  equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi * 1000,'b', ...
%  profil0d.psi(itime,:),pdederive(profil0d.psi(itime,:),profil0d.fdia(itime,:) .^ 2,0,2,2,1) ./ 2 *1000,'oc');drawnow;pause(3);
%  keyboard
  
  
  equilibrium.time_slice{itime}.profiles_1d.j_tor = profil0d.jli(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.j_parallel = profil0d.jeff(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.q = profil0d.qjli(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.magnetic_shear = pdederive(profil0d.rmx(itime,:),profil0d.qjli(itime,:),0,2,2,1) ./ profil0d.qjli(itime,:) .* profil0d.rmx(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.r_inboard = min(deuxd.R,[],2);
  equilibrium.time_slice{itime}.profiles_1d.r_outboard = max(deuxd.R,[],2);

  %equilibrium.time_slice{itime}.profiles_1d.li_3 = 2 .* cumtrapz(profil.rho,profil.bpol .^ 2 .* profil.vpr,2) ./  ...
  %    ((phys.mu0 .* max(1,cumtrapz(profil.rho,profil.jmoy .* profil.spr,2)) ) .^ 2  .* profil.Raxe);
  equilibrium.time_slice{itime}.profiles_1d.rho_tor = profil0d.rmx(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.rho_tor_norm = profil0d.rmx(itime,:) ./ profil0d.rmx(itime,end);
  equilibrium.time_slice{itime}.profiles_1d.dpsi_drho_tor = pdederive(profil0d.rmx(itime,:),profil0d.psi(itime,:),0,2,2,1);
  equilibrium.time_slice{itime}.profiles_1d.elongation = moment.e;
  equilibrium.time_slice{itime}.profiles_1d.triangularity_upper = moment.trh;
  equilibrium.time_slice{itime}.profiles_1d.triangularity_lower = moment.trl;
  equilibrium.time_slice{itime}.profiles_1d.volume = profil.volume;
  equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi = pdederive(profil.psi,profil.volume,1,2,2,1);
  % center value with q_0
  volmem_area = equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi(1);
  equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi(1)  = - equilibrium.time_slice{itime}.profiles_1d.q(1) .* ...
	deuxd.R(1,1) .^ 2 ./ (equilibrium.time_slice{itime}.profiles_1d.f(1))  .*  factor_two_pi;
  fact_area_dpsi = equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi(1) ./ volmem_area;
  % n'augmente pas la precision
  %equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi = equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi .* profil.vpr ./ ...
  %    (equilibrium.time_slice{itime}.profiles_1d.dvolume_dpsi(end) .* equilibrium.time_slice{itime}.profiles_1d.dpsi_drho_tor(end));
  equilibrium.time_slice{itime}.profiles_1d.dvolume_drho_tor = pdederive(profil.rho,profil.volume,0,2,2,1);
  equilibrium.time_slice{itime}.profiles_1d.area = profil.surface;
  equilibrium.time_slice{itime}.profiles_1d.darea_dpsi = pdederive(profil.psi,profil.surface,1,2,2,1);
  equilibrium.time_slice{itime}.profiles_1d.darea_dpsi(1) =  equilibrium.time_slice{itime}.profiles_1d.darea_dpsi(1) .* fact_area_dpsi;  
  % n'augmente pas la precision
  %equilibrium.time_slice{itime}.profiles_1d.darea_dpsi = equilibrium.time_slice{itime}.profiles_1d.darea_dpsi .* profil.spr ./ ...
  %    (equilibrium.time_slice{itime}.profiles_1d.darea_dpsi(end) .* equilibrium.time_slice{itime}.profiles_1d.dpsi_drho_tor(end));
  equilibrium.time_slice{itime}.profiles_1d.darea_drho_tor = pdederive(profil.rho,profil.surface,0,2,2,1);
  sext =  interp1_imas(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sext ,profil0d.temps(itime),'pchip','extrap');
  equilibrium.time_slice{itime}.profiles_1d.surface = profil0d.vpr_tor(itime,:) .* profil0d.grho(itime,:); 
  equilibrium.time_slice{itime}.profiles_1d.surface = (sext ./ equilibrium.time_slice{itime}.profiles_1d.surface(end)) .*  equilibrium.time_slice{itime}.profiles_1d.surface;  

  equilibrium.time_slice{itime}.profiles_1d.trapped_fraction = profil0d.ftrap(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.gm1 = profil0d.r2i(itime,:); 
  equilibrium.time_slice{itime}.profiles_1d.gm2 = profil0d.grho2r2(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.gm3 = profil0d.grho2(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.gm4 = profil.b2i;
  equilibrium.time_slice{itime}.profiles_1d.gm5 = profil.b2;
  equilibrium.time_slice{itime}.profiles_1d.gm6 = profil.grho2b2;
  equilibrium.time_slice{itime}.profiles_1d.gm7 = profil0d.grho(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.gm8 = profil.rmoy;
  equilibrium.time_slice{itime}.profiles_1d.gm9 = profil0d.ri(itime,:);
  equilibrium.time_slice{itime}.profiles_1d.b_field_average = sqrt(profil.b2); % must be check
  
  bmap = sqrt(deuxd.BR .^ 2 +deuxd.BZ .^ 2 + ((profil.fdia' * ones(1,size(deuxd.BR,2))) ./ deuxd.R) .^ 2);
  b_min = min(bmap,[],2);
  b_max = max(bmap,[],2);
  equilibrium.time_slice{itime}.profiles_1d.b_field_min = b_min;
  equilibrium.time_slice{itime}.profiles_1d.b_field_max = b_max;
  
  %equilibrium_cpo.profiles_1d.omega(itime,:)      = profil0d.omega(itime,:);
  %equilibrium_cpo.profiles_1d.omegaprime(itime,:) = pdederive(profil0d.psi(itime,:),profil0d.omega(itime,:),0,2,2,1);
  
  % calcul du nombre de mach
  %rhoi = 4.57e-3 .* sqrt(meff) .* sqrt(profil0d.tip(itime,:) ./ 1e3) ./ sqrt(profil.b2);
  %valf = sqrt(profil.b2) ./ sqrt(phys.mu0 .* rhoi);
  %equilibrium_cpo.profiles_1d.mach_a(itime,:) = sqrt(profil0d.vtor(itime,:) .^ 2  + profil0d.vtheta(itime,:) .^ 2) ./ valf;
  
  %equilibrium_cpo.profiles_1d.phi_flow(itime,:) =  profil0d.rmx(itime,:) .*  profil0d.vtheta(itime,:) .* profil0d.bpol(itime,:);
  % moment pour coherence interne
  equilibrium.time_slice{itime}.boundary.geometric_axis.z = moment.zaxe(1);
  
  % calcul de Jphi et jpar
  pprim  = equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi'  * ones(1,size(deuxd.R,2));
  ffprim = equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi' * ones(1,size(deuxd.R,2));
  pressure = equilibrium.time_slice{itime}.profiles_1d.pressure' * ones(1,size(deuxd.R,2)); 
  phi      = equilibrium.time_slice{itime}.profiles_1d.phi' * ones(1,size(deuxd.R,2)); 
  jphi  = pprim .* deuxd.R + ffprim ./ deuxd.R ./ phys.mu0; 
  jpar  = jphi .* deuxd.BPHI + ffprim ./ phys.mu0 ./ (equilibrium.time_slice{itime}.profiles_1d.f' * ones(1,size(deuxd.R,2))) .*  ...
	  (deuxd.BR .^ 2 +deuxd.BZ .^ 2);	  
  vphi   = (profil0d.omega(itime,:)' * ones(1,size(deuxd.R,2))) .* deuxd.R;
  bpol2d = abs(deuxd.BR + sqrt(-1) .*  deuxd.BZ);
  vtheta =  (profil0d.utheta(itime,:)' * ones(1,size(bpol2d,2))) .* bpol2d;
  b0    = equilibrium.time_slice{itime}.profiles_1d.f(end) ./ profil0d.Raxe(itime,end);
  jphi  = sign(mean(equilibrium.time_slice{itime}.profiles_1d.j_tor)) .* jphi;
  jpar  = sign(mean(equilibrium.time_slice{itime}.profiles_1d.j_parallel)) .* jpar ./ b0;
  RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
  ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
  theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
  equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim1 = deuxd.PSI(:,1);
  equilibrium.time_slice{itime}.profiles_2d{iprof}.grid.dim2 = deuxd.th(1,:);
  equilibrium.time_slice{itime}.profiles_2d{iprof}.grid_type.name = 'inverse rhopolar_polar';
  equilibrium.time_slice{itime}.profiles_2d{iprof}.grid_type.index = 2;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.grid_type.description = '2D polar coordinates (rho, theta) with magnetic axis as center of grid; theta and values following the corresponding COCOS convention, theta=atan2(z-zaxis,r-raxis) polar angle (COCOS=11 is assumed in ITER)';
  equilibrium.time_slice{itime}.profiles_2d{iprof}.r = deuxd.R;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.z = deuxd.Z;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.psi = deuxd.PSI;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.j_tor = jphi;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.j_parallel = jpar;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.b_r = deuxd.BR;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.b_z = deuxd.BZ;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_r = deuxd.BR; 
  equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_z = deuxd.BZ;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.b_field_tor = deuxd.BPHI;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.theta = theta;
  equilibrium.time_slice{itime}.profiles_2d{iprof}.phi = phi;
   
  % filling outline with profiles_2d LCFS if it is empty (no real LCFS provided to METIS)
  if outlineisfilled == 0
      R_LCFS_pd2 = deuxd.R(end,:);
      Z_LCFS_pd2 = deuxd.Z(end,:);
      
      % LFCS must be open in IMAS
      if (R_LCFS_pd2(end) == R_LCFS_pd2(1)) && (Z_LCFS_pd2(end) == Z_LCFS_pd2(1))
          equilibrium.time_slice{itime}.boundary.outline.r = R_LCFS_pd2(1:end-1);
          equilibrium.time_slice{itime}.boundary.outline.z = Z_LCFS_pd2(1:end-1);
      else
          equilibrium.time_slice{itime}.boundary.outline.r = R_LCFS_pd2;
          equilibrium.time_slice{itime}.boundary.outline.z = Z_LCFS_pd2;
      end
      disp('Using profiles_2d data to provide LCFS');
  end
  
  % backward compatibility
  equilibrium.time_slice{itime}.boundary.lcfs.r = equilibrium.time_slice{itime}.boundary.outline.r;
  equilibrium.time_slice{itime}.boundary.lcfs.z = equilibrium.time_slice{itime}.boundary.outline.z;
    
  % initilialisation substructure 
  grid2 = equilibrium.time_slice{1}.coordinate_system.grid;
  % remplissage de coord_sys
  grid2.dim1       = deuxd.PSI(:,1);
  grid2.dim2       = deuxd.th(1,:);
  equilibrium.time_slice{itime}.coordinate_system.grid_type.name = 'inverse rhopolar_polar';
  equilibrium.time_slice{itime}.coordinate_system.grid_type.index = 2;
  equilibrium.time_slice{itime}.coordinate_system.grid_type.description = '2D polar coordinates (rho, theta) with magnetic axis as center of grid; theta and values following the corresponding COCOS convention, theta=atan2(z-zaxis,r-raxis) polar angle (COCOS=11 is assumed in ITER)';
  equilibrium.time_slice{itime}.coordinate_system.grid = grid2;  
  jacobian                                       = abs(2 .* pi .* deuxd.R .* (deuxd.dRdPSI .* deuxd.dZdth - deuxd.dZdPSI .* deuxd.dRdth));
  jacobian(1,:)                                  = 0;
  equilibrium.time_slice{itime}.coordinate_system.jacobian = abs(jacobian);
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,1,1) = (deuxd.dPSIdR .^ 2 + deuxd.dPSIdZ .^ 2);
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,1,2) = deuxd.dthdR .* deuxd.dPSIdR +  deuxd.dthdZ .* deuxd.dPSIdZ;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,1,3) = 0;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,2,1) =  deuxd.dthdR .* deuxd.dPSIdR +  deuxd.dthdZ .* deuxd.dPSIdZ;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,2,2) =  deuxd.dthdR .^2 + deuxd.dthdZ .^ 2;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,2,3) = 0;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,3,1) = 0;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,3,2) = 0;
  equilibrium.time_slice{itime}.coordinate_system.tensor_contravariant(:,:,3,3) = 1./ deuxd.R .^2;
  equilibrium.time_slice{itime}.coordinate_system.r        = deuxd.R;
  equilibrium.time_slice{itime}.coordinate_system.z        = deuxd.Z;
  
  % equilibre etendu hors du plasma
  % EVALUATION IN THE VACCUUM
  method_extrap = NaN;
  % for testing 
  if equi_extrap == 0
    %equivide = zfitvide_interp(equicronos,0,[]);
    equivide = zfit_poly_interp(equicronos,0,[]);
    method_extrap = -1;
  else
    method_extrap = 1;
    if ~isempty(control) && isfield(control,'x_point') &&  ~isempty(control.x_point)
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
    if (equivide.dpsi_error > 0.1) || (equivide.dB_error > 0.2)   
        fprintf('changing to direct methode |');
        method_extrap = 0;
        dpsi_error_mem = equivide.dpsi_error;
        equivide = zfitvide_interp(equicronos,0,[]);
        equivide.dpsi_error = equivide.dpsi_error + dpsi_error_mem;
    end 
  end
  equilibrium.time_slice{itime}.global_quantities.psi_boundary_error_upper = abs(equivide.dpsi_error);
  equilibrium.time_slice{itime}.global_quantities.psi_boundary_error_lower = abs(equivide.dpsi_error);

  % initialisation sub structure grid 
  grid3 = equilibrium.time_slice{1}.profiles_2d{1}.grid;
  % (R,Z) GRID
  a   = (max(equivide.R(:)) - min(equivide.R(:))) ./ 2;
  if isfield(z0dstruct.z0dinput.option,'fixed_grid') && (z0dstruct.z0dinput.option.fixed_grid == 1)
      grid3.dim1 = linspace(box.rmin,box.rmax,nbeqdsk);
      grid3.dim2 = linspace(box.zmin,box.zmax,ceil(nbeqdsk .* rzp)) - moment.zaxe(1); 
  else
      grid3.dim1 = linspace(max(eps,min(equivide.R(:)) - a ./ 4) ,max(equivide.R(:)) + a ./ 4,nbeqdsk);
      grid3.dim2 = linspace(min(equivide.Z(:)) - ( a ./ 4 .* rzp),max(equivide.Z(:)) + ( a ./ 4 .* rzp),ceil(nbeqdsk .* rzp));
  end
    
  [r2d,z2d] = meshgrid(grid3.dim1,grid3.dim2);

  % multi matlab version
%    if 0
%    try
%      tri          = delaunay(r2d',z2d',{'Qt','Qbb','Qc','Qz'});
%    catch
%      tri          = delaunay(r2d',z2d');
%    end
%    if ~isfield(grid3,'connect')
%      grid3.connect(1:size(tri,1),1:3)  = tri;
%    elseif size(tri,1) <= size(grid3.connect,2)
%      grid3.connect(1:size(tri,1),1:3)  = tri;
%    else
%      connect_mem = grid3.connect(1:max(1,itime-1),:,:);
%      grid3.connect    = 0 .* ones(size(profil0d.psi,1),size(tri,1),size(connect_mem,3));
%      for lk =1:max(1,itime-1)
%        grid3.connect(lk,1:size(connect_mem,2),:)  = connect_mem(lk,:,:);
%      end
%      grid3.connect(itime,1:size(tri,1),1:3)  = tri;
%    end
%    end
    if method_extrap == -1
        [psi2d,BR2d,BZ2d] = zpsi_poly_interp(equivide,r2d,z2d);
    elseif method_extrap == 0
        [psi2d,BR2d,BZ2d] = zpsivide_interp(equivide,r2d,z2d);
    else
        [psi2d,BR2d,BZ2d] = zpsivide_updown_16_jorek(equivide,r2d,z2d);
    end
  BR2d  = BR2d .* sign(factor_two_pi);
  BZ2d  = BZ2d .* sign(factor_two_pi);
  psi2d = psi2d .* factor_two_pi;
  mask = zinout(equivide.R(end,:),equivide.Z(end,:),r2d,z2d);
  BPHI2d   = profil0d.fdia(itime,end) ./ r2d;
  warning off
  if all(diff(equilibrium.time_slice{itime}.profiles_1d.psi) > 0) || all(diff(equilibrium.time_slice{itime}.profiles_1d.psi) < 0)
    BPHI2d(mask)  = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.f,psi2d(mask),'pchip','extrap') ./ r2d(mask);
  else
    warning('Problem for B_field_tor interolation in 2D');
    BPHI2d(mask)  = profil0d.fdia(itime,1)./ r2d(mask);
  end
%    bpol_test = griddata(r2d,z2d,abs(BR2d + sqrt(-1) .* BZ2d),equivide.R,equivide.Z);
%    bpol_test = mean(bpol_test,2);
%    bpol_test_alt = mean(abs(deuxd.BR + sqrt(-1) .* deuxd.BZ),2);
%    figure(32);plot(profil0d.xli,profil0d.bpol(itime,:),'r',profil0d.xli,bpol_test,'b',profil0d.xli,bpol_test_alt,'k');drawnow
  %figure(31);clf;
  %contour(r2d,z2d,psi2d,equilibrium.time_slice{itime}.profiles_1d.psi,'r');
  %hold on
  %plot(equilibrium.time_slice{itime}.profiles_2d{iprof}.r',equilibrium.time_slice{itime}.profiles_2d{iprof}.z',':k')
  if all(diff(equilibrium.time_slice{itime}.profiles_1d.psi) > 0) || all(diff(equilibrium.time_slice{itime}.profiles_1d.psi) < 0)
      pprim  = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi, ...
                       psi2d,'pchip','extrap');
      pprim(~mask) = 0;            
      ffprim  = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi, ...
                       psi2d,'pchip','extrap');
      ffprim(~mask) = 0;            
      jphi  = pprim .* r2d + ffprim ./ r2d ./ phys.mu0; 
      jphi(~mask) = 0;
      jpar  = jphi .* BPHI2d + ffprim ./ phys.mu0 ./ (BPHI2d .* r2d) .*  ...
	      (BR2d .^ 2 + BZ2d .^ 2);
      jpar(~mask) = 0;
      b0    = equilibrium.time_slice{itime}.profiles_1d.f(end) ./ profil0d.Raxe(itime,end);
      jpar  = jpar ./ b0;

      psin    =  abs(profil0d.psi(itime,:) - profil0d.psi(itime,end)) ./ abs(profil0d.psi(itime,1) - profil0d.psi(itime,end));
      psin_loc = abs(equilibrium.time_slice{itime}.profiles_1d.psi - equilibrium.time_slice{itime}.profiles_1d.psi(end)) ./ ...
                 abs(equilibrium.time_slice{itime}.profiles_1d.psi(1) - equilibrium.time_slice{itime}.profiles_1d.psi(end));
      omega    = interp1(psin,profil0d.omega(itime,:),psin_loc,'pchip','extrap');
      vphi    = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,omega,psi2d,'pchip','extrap') .* r2d;
      vphi(~mask) = 0;
      utheta    = interp1(psin,profil0d.utheta(itime,:),psin_loc,'pchip','extrap');
      vtheta    = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,utheta,psi2d,'pchip','extrap') .* abs(BR2d + sqrt(-1) .* BZ2d);
      vtheta(~mask) = 0;

      phi  = interp1(equilibrium.time_slice{itime}.profiles_1d.psi,equilibrium.time_slice{itime}.profiles_1d.phi, ...
                       psi2d,'pchip','extrap');  
                       
      %contour(r2d,z2d,phi,equilibrium.time_slice{itime}.profiles_1d.phi,'b-.');
      %drawnow
      phi(~mask) = NaN;
      RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
      ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
      theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
      theta  = griddata(deuxd.R,deuxd.Z,theta,r2d,z2d,'cubic');
      theta(~mask) = 0;
  
  else
      pprim  = equilibrium.time_slice{itime}.profiles_1d.dpressure_dpsi'  * ones(1,size(equivide.R,2));
      pprim  = griddata(equivide.R,equivide.Z,pprim,r2d,z2d,'cubic');
      ffprim = equilibrium.time_slice{itime}.profiles_1d.f_df_dpsi'  * ones(1,size(equivide.R,2));
      ffprim  = griddata(equivide.R,equivide.Z,ffprim,r2d,z2d,'cubic');
      jphi  = pprim .* r2d + ffprim ./ r2d ./ phys.mu0; 
      jphi(~mask) = 0;
      jpar  = jphi .* BPHI2d + ffprim ./ phys.mu0 ./ (BPHI2d .* r2d) .*  ...
	      (BR2d .^ 2 + BZ2d .^ 2);
      jpar(~mask) = 0;
      vphi   = (profil0d.omega(itime,:)' * ones(1,size(deuxd.R,2))) .* deuxd.R;
      bpol2d = abs(deuxd.BR + sqrt(-1) .*  deuxd.BZ);
      vtheta =  (profil0d.utheta(itime,:)' * ones(1,size(bpol2d,2))) .* bpol2d;
      b0    = equilibrium.time_slice{itime}.profiles_1d.f(end) ./ profil0d.Raxe(itime,end);
      jpar  = jpar ./ b0;
      vphi    = griddata(deuxd.R,deuxd.Z,vphi,r2d,z2d,'cubic');
      vphi(~mask) = 0;
      vtheta  = griddata(deuxd.R,deuxd.Z,vtheta,r2d,z2d,'cubic');
      vtheta(~mask) = 0;
      %pressure = equilibrium_cpo.profiles_1d.pressure(itime,:)' * ones(1,size(deuxd.R,2)); 
      %pressure      = griddata(deuxd.R,deuxd.Z,pressure,r2d,z2d,'cubic');
      %pressure(~mask) = 0;
      phi      = equilibrium.time_slice{itime}.profiles_1d.phi' * ones(1,size(deuxd.R,2)); 
      phi      = griddata(deuxd.R,deuxd.Z,phi,r2d,z2d,'cubic');
      phi(~mask) = NaN;
      RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
      ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
      theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
      theta  = griddata(deuxd.R,deuxd.Z,theta,r2d,z2d,'cubic');
      theta(~mask) = 0;
  end
  
  % initialisation sub structure
  profiles_2d_2{1} = equilibrium.time_slice{1}.profiles_2d{1};

  warning on
  profiles_2d_2{1}.psi 	    = psi2d';
  profiles_2d_2{1}.j_tor       = sign(mean(equilibrium.time_slice{itime}.profiles_1d.j_tor)) .* jphi';
  profiles_2d_2{1}.j_parallel  = sign(mean(equilibrium.time_slice{itime}.profiles_1d.j_parallel)) .* jpar';
  profiles_2d_2{1}.b_r	    = BR2d';
  profiles_2d_2{1}.b_z 	    = BZ2d';
  profiles_2d_2{1}.b_field_r   = BR2d';
  profiles_2d_2{1}.b_field_z   = BZ2d';
  profiles_2d_2{1}.b_field_tor = BPHI2d';
  profiles_2d_2{1}.r           = r2d';
  profiles_2d_2{1}.z           = z2d' + moment.zaxe(1);
  %profiles_2d_2.vtheta 	 = vtheta';
  %profiles_2d_2.vphi 	 = vphi';
  profiles_2d_2{1}.theta 	 = theta';
  profiles_2d_2{1}.phi 	 = phi';
  %profiles_2d_2.pressure 	 = pressure';
  profiles_2d_2{1}.grid = grid3;
  profiles_2d_2{1}.grid_type.name = 'rectangular';
  profiles_2d_2{1}.grid_type.index = 1;
  profiles_2d_2{1}.grid_type.description = 'cylindrical R,Z ala eqdsk, within the corresponding COCOS convention (COCOS=11 is assumed in ITER)';

  fprintf('\n');
  
  

  %% ORDER CONVENTION (WPCD, IMAS?)
  profiles_2d_1{1} = equilibrium.time_slice{itime}.profiles_2d{1};
  equilibrium.time_slice{itime} = rmfield(equilibrium.time_slice{itime},'profiles_2d');     
  equilibrium.time_slice{itime}.profiles_2d{1} = profiles_2d_2{1};
  equilibrium.time_slice{itime}.profiles_2d{2} = profiles_2d_1{1};
  

%     %% graphe for control
  if itime ==  -1
	figure(31)
	clf
	plot(profiles_2d_1{1}.r',profiles_2d_1{1}.z','m');
	hold on
	contour(profiles_2d_2{1}.r,profiles_2d_2{1}.z,profiles_2d_2{1}.psi,equilibrium.time_slice{itime}.profiles_1d.psi,'color','c','linestyle','--')
	contour(profiles_2d_2{1}.r,profiles_2d_2{1}.z,profiles_2d_2{1}.psi,101,'color','g','linestyle','-')
	quiver(profiles_2d_2{1}.r,profiles_2d_2{1}.z,profiles_2d_2{1}.b_field_r,profiles_2d_2{1}.b_field_z,'color','r');
	quiver(profiles_2d_1{1}.r,profiles_2d_1{1}.z,profiles_2d_1{1}.b_field_r,profiles_2d_1{1}.b_field_z,'color','b');
	if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa')
	  plot(profil0d.Rsepa(itime,:),profil0d.Zsepa(itime,:) + z0_offset ,'k');
	end	
	
	plot(equilibrium.time_slice{itime}.boundary.outline.r,equilibrium.time_slice{itime}.boundary.outline.z,'k.');
	plot(equilibrium.time_slice{itime}.boundary.geometric_axis.r ,equilibrium.time_slice{itime}.boundary.geometric_axis.z,'or');
	plot(equilibrium.time_slice{itime}.global_quantities.magnetic_axis.r ,equilibrium.time_slice{itime}.global_quantities.magnetic_axis.z ,'r+');
      %    wall = load('jt60sawall');
      %    rwall = wall.wall.data(1:end,1);
      %    zwall = wall.wall.data(1:end,2);
      %    plot(rwall,zwall,'k');
	
      xlabel('R (m)');
      ylabel('Z (m)');
      drawnow
      keyboard
	  pause(1);
  end
  % adding ggd (copy of profiles_2d_2) to be compatible with ProfileMaker and associated tools
  ggd = imas_write_regular_grid(profiles_2d_2{1},equilibrium.time_slice{1}.ggd{1});
  equilibrium.time_slice{itime}.ggd{1} = ggd;
  ggd = imas_write_regular_grid(profiles_2d_1{1},equilibrium.time_slice{1}.ggd{1});
  equilibrium.time_slice{itime}.ggd{2} = ggd;
  % new grid definition
  equilibrium.grids_ggd{itime}.grid{1} = equilibrium.time_slice{itime}.ggd{1}.grid;
  equilibrium.grids_ggd{itime}.grid{2} = equilibrium.time_slice{itime}.ggd{2}.grid;
  equilibrium.grids_ggd{itime}.time    = equilibrium.time(itime);
  % COCOS
  s_ip  = sign(mean(sign(equilibrium.time_slice{itime}.global_quantities.ip)));
  s_b0  = sign(mean(sign(equilibrium.vacuum_toroidal_field.b0)));
  
  s_psi = sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.psi(end) - equilibrium.time_slice{itime}.profiles_1d.psi(1))));
  s_bp  = s_psi .* s_ip;
  s_rho_phi_theta = sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.q(:)))) .* s_ip .* s_b0;
  if sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.f(:)))) ~= s_b0
    warning('mapequilibrium : sign missmatch between B0 and Fdia');
  end
  % this test is false, see reference paper on COCOS
  %if sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.phi(:)))) ~= s_b0
  %   warning('mapequilibrium : sign missmatch between B0 and Phi');
  %end

  if  sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.j_tor(:)))) ~= s_ip
    warning('mapequilibrium : sign missmatch between Ip and jmoy');
  end
  if  sign(mean(sign(equilibrium.time_slice{itime}.profiles_1d.j_parallel(:)))) ~= s_ip
    warning('mapequilibrium : sign missmatch between Ip and jeff');
    
  end

end %% END LOOP ITIME=1:NTIME


% Call of fixed boundary equilibrium solver of FFEQS.M code in post
% processing
if FEEQS_post
    disp('Call of fixed boundary solver of FEEQS.M code in post processing of METIS IMAS interface:')
    try
         ids_equilibrium_out = improve_equilibrium_with_feeqs(equilibrium,false,alternative_extrapolation);
    catch
        disp('Problem during the call of FFEQS code in post processing of METIS IMAS interface:');
        disp(lasterr)
        ids_equilibrium_out = [];
    end
    if ~isempty(ids_equilibrium_out)
        disp('METIS MHD equilibrium will be replaced by the solution computed in post processing using fixed boundary solver of FEEQS.M code.');
        equilibrium = ids_equilibrium_out;
    end
end



%% ----------------------------------------------------------------------------------------------------------


% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_imas(t,y,tt,varargin)


if length(t) <= 1
	rep = y(end,:);
else
	rep = interp1(t,y,tt,varargin{:});
end

% surcouche pour traiter le cas a 1 temps en meme temps
function rep = griddata_imas(t,x,y,tt,xx,methode)


if size(t,1) <= 1
	if all(x == xx)
		rep = y;
	else
		switch methode
		case 'cubic'
			methode = 'pchip';
		end
		rep  = interp1(x,y,xx,methode,'extrap');
	end
else
    try
        rep = griddata(t,x,y,tt,xx,methode);
    catch
        rep = griddata(t,x,y,tt,xx,methode,{'Qt','Qbb','Qc','Qz'});        
    end
end


% smooth
function psi2d = zsmooth_psi(psi2d)

% smooth si necessaire cas difficile
f33 = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
psi_mem = psi2d;
psi2d = conv2(psi2d ,f33,'same');
psi2d(1,:)   = psi_mem(1,:);
psi2d(end,:) = psi_mem(end,:);
psi2d(:,1)   = psi_mem(:,1);
psi2d(:,end) = psi_mem(:,end);




function [u,upr,upz] = f(k,r,z,c1,c2)

u    =  c1 .* sin(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z) +  ...
        c2 .* cos(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z);
upr  = -c1 .* cos(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z) ./ r .^ 2 .* sqrt(k) .* sqrt(k - 0.2e1) .* z +  ...
        c2 .* sin(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z) ./ r .^ 2 .* sqrt(k) .* sqrt(k - 0.2e1) .* z;
upz  =  c1 .* cos(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z) ./ r .* sqrt(k) .* sqrt(k - 0.2e1) -  ...
        c2 .* sin(0.1e1 ./ r .* sqrt(k) .* sqrt(k - 0.2e1) .* z) ./ r .* sqrt(k) .* sqrt(k - 0.2e1);

	      

% ZINOUT determine les points interieurs et exterieurs au plasma
%---------------------------------------------------------------
% fichier zinout.m ->  zinout
%
%
% fonction Matlab 5 :
%
% Cette fonction determine les points interieurs et exterieurs au plasma.
% La frontiere appartient au plasma.
%  
% syntaxe  :
%  
%     mask = zinout(R,Z,r,z,stric);
%    
% entree :
%
%     R, Z =  coordonnees de la derniere surface magnetique (vecteur)
%     r, z =  coordonneees des points a tester [M,N]
%     strict = si = 0 la frontiere est contenu dans la surface;
%             si = 1 que les points strictement dans la surface
%             si = -1 que les points sur le contour
%             si = -2 que les points du contour
%
% sorties :
% 
%     mask =  mask  des points interieurs (1 = interieur, 0 = exterieur) [M,N]
%
% fonction ecrite par J-F Artaud , poste 62-15
% version 4.1 , du 06/10/2008.
% 
% 
%--------------------------------------------------------------
%
function mask_out = zinout(R,Z,r,z,strict)

if nargin < 5
	strict = 1;
elseif isempty(strict)
	strict = 1;
end


% mise en ligne
R    = R(:);
Z    = Z(:);

% securite ITM
u = R + sqrt(-1) * Z;
indbad = find(diff(u) == 0);
while(~isempty(indbad))
  u(indbad) = [];
  indbad = find(diff(u) == 0);
end
R = real(u);
Z = imag(u);

if (R(1) ~= R(end)) | (Z(1) ~= Z(end))
	R(end+1) = R(1);
	Z(end+1) = Z(1);
end

% inpolygon is faster
[IN ON] = inpolygon(r, z, R, Z);
switch strict
  case 0
    mask_out = IN | ON;
  case 1
    mask_out = IN;
  case -1
    mask_out = ON;
  case -2
    [void,indice] = intersect(r + sqrt(-1) * z , R + sqrt(-1) * Z);
    mask_out = zeros(size(r));
    mask_out(indice) = 1;
  otherwise
    error('input strict value not supported')
end



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
%% ref : P. J. Mc Carthy POP 1999 p 3554-...
% fixed boudnary solution
%psiv = psiv +  A .* R .^ 4 ./ 8 + B .* Z .^2 ./ 2;
%BRv  = BRv  - B .* Z ./ R; 
%BZv  = BZv +  A .* R .^ 2 ./ 2;
%  % Antoine J. Cerfon and Jeffrey P. Freidberg, Citation: Physics of Plasmas (1994-present) 17, 032502 (2010); doi: 10.1063/1.3328818
%  % replaced by FBE solution
%  psiv = psiv + A .* R .^ 4 ./ 8 +  0 .* Z .^2 ./ 2  + B .* 0.5 .* R .* log(R);
%  BRv  = BRv  -  0 .* Z ./ R; 
%  BZv  = BZv  + A .* 0.5 .* R .^ 2  +  0 .* Z .^2 ./ 2 + B .* (log(R) + 0.5);
psiv = psiv +  A .* R  .^ 4 ./ 8 + B .* Z  .^ 2 ./ 2 ;
BRv  = BRv  - B .* Z  ./ R ; 
BZv  = BZv +  A .* R  .^ 2 ./ 2 ;

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
psi_7 = 8 .* Z .^ 6  - 140 .* Z .^ 4 .* R .^ 2 + 75 .* Z .^ 2 .* R .^ 4  - 15  .* R .^ 6 .* log(R) + 180 .* R .^ 4 .* Z .^ 2 .* log(R) - 120 .* R .^ 2 .* Z .^ 4 .* log(R);
BR_7  = -(48 * Z .^ 5 - 560 .* Z .^ 3  .* R .^ 2 + 150 .* Z  .* R .^ 4 + 360 .* R .^ 4 .* Z .* log(R) - 480 .* R .^ 2 .* Z .^ 3 .* log(R))./ R;
BZ_7  = - 280 .* Z .^ 4 + 300 .* Z .^ 2 .* R .^ 2 - 90  .* R .^ 4 .* log(R) - 15  .* R .^ 4 + 720 .* R .^ 2 .* Z .^ 2 .* log(R) + 180 .* R .^ 2 .* Z .^ 2 - 240  .* Z .^ 4 .* log(R) - 120  .* Z .^ 4;
%BR_7  = -(48 .* Z .^ 5 + (-480 .* R .^ 2 .* log(R) - 560 .* R .^ 2) .* Z .^ 3 + 150 .* R .^ 4 .* Z + 360 .* R .^ 4 .* log(R)) ./ R;
%BR_7  = -(48 * Z .^ 5 -560 .* Z .^ 3 .* R .^ 2 + 150 .* Z  .* R .^ 4 + 360 .* R .^ 4 .* Z  .* log(R) - 480 .* R .^ 2 .* Z .^ 3 .* log(R)) ./ R;
%BZ_7  = (-240 .* log(R) - 400) .* Z .^ 4 + 300 .* R .^ 2 .* Z .^ 2 + (1440 .* R .^ 2 .* log(R) + 360 .* R .^ 2) .* Z - 90 .* R .^ 4 .* log(R) - 15 .* R .^ 4;
%BZ_7  = - 280 .* Z .^ 4 + 300 .* Z .^ 2 .* R .^ 2  - 90 .* R .^ 4 .* log(R) - 15  .* R .^ 4  + 720 .* R .^ 2 .* Z .* 2 .* log(R) + 180 .* R .^ 2 .* Z .^ 2  - 240 .* Z .^ 4 .* log(R) - 120  .* Z .^ 4;

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

if any(imag(psi_out(:)))
    keyboard
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
pause(1);

%  % function to test alternative to griddata
%  function out = griddata_alt(xin,yin,var,xout,yout,methode)
%  ss = size(var);
%  
%  F = scatteredInterpolant(xin(:),yin(:),var(:),'natural','linear');
%  out = F(xout,yout);
%  out2 = griddata(xin,yin,var,xout,yout,methode);
%  delta = out(:) - out2(:);
%  sqrt(sum((delta(isfinite(delta)) .^ 2))./ sum(out2(isfinite(delta)) .^ 2))

function control =sepa_moments(Rsepa,Zsepa,x_point)

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

% compute squareness
% reference: T. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009
% close contour if needed
if (Rsepa(1) ~= Rsepa(end)) || (Zsepa(1) ~= Zsepa(end))
   Rsepa(end+1) = Rsepa(1);
   Zsepa(end+1) = Zsepa(1);  
end
% angle vector for ellipse
uel = linspace(0,2*pi,401);
% outter part
% more outer point is (rmax,z0)
% more lower point is (rzmin,zmin)
% more upper point is (rzmax,zmax)
R_E_point_lower = rmax;
Z_E_point_lower = zmin;
R_O_point_lower = rzmin;
Z_O_point_lower = z0;
R_E_point_upper = rmax;
Z_E_point_upper = zmax;
R_O_point_upper = rzmax;
Z_O_point_upper = z0;
rep_lower = zinterx(cat(2,Rsepa(:),Zsepa(:))', ...
            cat(2,cat(1,R_O_point_lower,R_E_point_lower), ...
            cat(1,Z_O_point_lower,Z_E_point_lower))');
rep_upper = zinterx(cat(2,Rsepa(:),Zsepa(:))', ...
            cat(2,cat(1,R_O_point_upper,R_E_point_upper), ...
            cat(1,Z_O_point_upper,Z_E_point_upper))');
if isempty(rep_upper)
    control.squareness_upper_outer = NaN;
else
    L_OD = sqrt((R_O_point_upper - rep_upper(1)).^2 + (Z_O_point_upper - rep_upper(2)).^2);
    L_OE = sqrt((R_O_point_upper - R_E_point_upper).^2 + (Z_O_point_upper - Z_E_point_upper).^2);
    cost = (R_O_point_upper - R_E_point_upper) / L_OE;
    sint = (Z_O_point_upper - Z_E_point_upper)  / L_OE;
    l_h  = R_O_point_upper - rmax;
    l_v  = Z_O_point_upper - zmax;
    Rel  = R_O_point_upper + l_h * cos(uel);
    Zel  = Z_O_point_upper + l_v * sin(uel);
    rep_el = zinterx(cat(2,Rel(:),Zel(:))', ...
            cat(2,cat(1,R_O_point_upper,R_E_point_upper), ...
            cat(1,Z_O_point_upper,Z_E_point_upper))');   
    L_OC = sqrt((R_O_point_upper - rep_el(1)).^2 + (Z_O_point_upper - rep_el(2)).^2);
    control.squareness_upper_outer = (L_OD -L_OC) ./ abs(L_OE- L_OC);
end
if isempty(rep_lower)
    control.squareness_lower_outer = NaN;
else
    L_OD = sqrt((R_O_point_lower - rep_lower(1)).^2 + (Z_O_point_lower - rep_lower(2)).^2);
    L_OE = sqrt((R_O_point_lower - R_E_point_lower).^2 + (Z_O_point_lower - Z_E_point_lower).^2);
    cost = (R_O_point_lower - R_E_point_lower) / L_OE;
    sint = (Z_O_point_lower - Z_E_point_lower)  / L_OE;
    l_h  = R_O_point_lower - rmax;
    l_v  = Z_O_point_lower - zmin;
    Rel  = R_O_point_lower + l_h * cos(uel);
    Zel  = Z_O_point_lower + l_v * sin(uel);
    rep_el = zinterx(cat(2,Rel(:),Zel(:))', ...
            cat(2,cat(1,R_O_point_lower,R_E_point_lower), ...
            cat(1,Z_O_point_lower,Z_E_point_lower))');   
    L_OC = sqrt((R_O_point_lower - rep_el(1)).^2 + (Z_O_point_lower - rep_el(2)).^2);
    control.squareness_lower_outer = (L_OD -L_OC) ./ abs(L_OE - L_OC);
end
%inner part
R_E_point_lower = rmin;
Z_E_point_lower = zmin;
R_O_point_lower = rzmin;
Z_O_point_lower = z0;
R_E_point_upper = rmin;
Z_E_point_upper = zmax;
R_O_point_upper = rzmax;
Z_O_point_upper = z0;
rep_lower = zinterx(cat(2,Rsepa(:),Zsepa(:))', ...
            cat(2,cat(1,R_O_point_lower,R_E_point_lower), ...
            cat(1,Z_O_point_lower,Z_E_point_lower))');
rep_upper = zinterx(cat(2,Rsepa(:),Zsepa(:))', ...
            cat(2,cat(1,R_O_point_upper,R_E_point_upper), ...
            cat(1,Z_O_point_upper,Z_E_point_upper))');
if isempty(rep_upper)
    control.squareness_upper_inner = NaN;
else
    L_OD = sqrt((R_O_point_upper - rep_upper(1)).^2 + (Z_O_point_upper - rep_upper(2)).^2);
    L_OE = sqrt((R_O_point_upper - R_E_point_upper).^2 + (Z_O_point_upper - Z_E_point_upper).^2);
    cost = (R_O_point_upper - R_E_point_upper) / L_OE;
    sint = (Z_O_point_upper - Z_E_point_upper)  / L_OE;
    l_h  = R_O_point_upper - rmin;
    l_v  = Z_O_point_upper - zmax;
    Rel  = R_O_point_upper + l_h * cos(uel);
    Zel  = Z_O_point_upper + l_v * sin(uel);
    rep_el = zinterx(cat(2,Rel(:),Zel(:))', ...
            cat(2,cat(1,R_O_point_upper,R_E_point_upper), ...
            cat(1,Z_O_point_upper,Z_E_point_upper))');   
    L_OC = sqrt((R_O_point_upper - rep_el(1)).^2 + (Z_O_point_upper - rep_el(2)).^2);

    control.squareness_upper_inner = (L_OD - L_OC) ./ abs(L_OE - L_OC);
end
if isempty(rep_lower)
    control.squareness_lower_inner = NaN;
else
    L_OD = sqrt((R_O_point_lower - rep_lower(1)).^2 + (Z_O_point_lower - rep_lower(2)).^2);
    L_OE = sqrt((R_O_point_lower - R_E_point_lower).^2 + (Z_O_point_lower - Z_E_point_lower).^2);
    cost = (R_O_point_lower - R_E_point_lower) / L_OE;
    sint = (Z_O_point_lower - Z_E_point_lower)  / L_OE;
    l_h  = R_O_point_lower - rmin;
    l_v  = Z_O_point_lower - zmin;
    Rel  = R_O_point_lower + l_h * cos(uel);
    Zel  = Z_O_point_lower + l_v * sin(uel);
    rep_el = zinterx(cat(2,Rel(:),Zel(:))', ...
            cat(2,cat(1,R_O_point_lower,R_E_point_lower), ...
            cat(1,Z_O_point_lower,Z_E_point_lower))');   
    L_OC = sqrt((R_O_point_lower - rep_el(1)).^2 + (Z_O_point_lower - rep_el(2)).^2);
    control.squareness_lower_inner = (L_OD - L_OC) ./ abs(L_OE - L_OC);
end  


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
nxpts = 1;
if el > 2e-2
  indlz  = find(zl == min(zl),1);
  control.x_point{nxpts} = x_point;
  control.x_point{nxpts}.r = rl(indlz);
  control.x_point{nxpts}.z = zl(indlz);
  nxpts = nxpts + 1;
end
if eh > 2e-2
  indhz  = find(zh == max(zh),1);
  control.x_point{nxpts} = x_point;
  control.x_point{nxpts}.r = rh(indhz);
  control.x_point{nxpts}.z = zh(indhz);
end

% function to map ggd from regular grid store in profiles_2d
% thank to Th. Aniel
function ggd = imas_write_regular_grid(p2d,ggd)

ggd.grid.identifier.name        = p2d.grid_type.name;
ggd.grid.identifier.index       = p2d.grid_type.index;
ggd.grid.identifier.description = p2d.grid_type.description; 
%
ggd.grid.space{1}.geometry_type.name        = [];
ggd.grid.space{1}.geometry_type.index       =  -999999999;
ggd.grid.space{1}.geometry_type.description = [];   
ggd.grid.space{1}.coordinates_type    = [];
%	
n1 = length(p2d.grid.dim1); 
n2 = length(p2d.grid.dim2); 

for k1 = [1:n1]
   % initialisation sustructure
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1} = ggd.grid.space{1}.objects_per_dimension{1}.object{1};
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.boundary{1}.index      = -999999999;
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.boundary{1}.neighbours = [];
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.geometry(1)            = p2d.grid.dim1(k1);
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.nodes                  = [];
   ggd.grid.space{1}.objects_per_dimension{1}.object{k1}.measure                = -9.0e40;
   
end

% initialisation substructure
ggd.grid.space{2} = ggd.grid.space{1};
for k2 = [1:n2]
   % initialisation sustructure
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2} = ggd.grid.space{2}.objects_per_dimension{1}.object{1};		
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.boundary{1}.index      = -999999999;
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.boundary{1}.neighbours = [];
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.geometry(1)            = p2d.grid.dim2(k2);
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.nodes                  = [];
   ggd.grid.space{2}.objects_per_dimension{1}.object{k2}.measure                = -9.0e40;
   
end
   
ggd.grid.grid_subset{1}.identifier.name                = '';
ggd.grid.grid_subset{1}.identifier.index               = -999999999;
ggd.grid.grid_subset{1}.identifier.description         = '';
ggd.grid.grid_subset{1}.dimension                      = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.space     = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.dimension = -999999999;
ggd.grid.grid_subset{1}.element{1}.object{1}.index     = -999999999;
ggd.grid.grid_subset{1}.base{1}.jacobian               = [];
ggd.grid.grid_subset{1}.base{1}.g11_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g12_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g13_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g21_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g22_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g23_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g31_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g32_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g33_covariant          = [];
ggd.grid.grid_subset{1}.base{1}.g11_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g12_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g13_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g21_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g22_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g23_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g31_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g32_contravariant      = [];
ggd.grid.grid_subset{1}.base{1}.g33_contravariant      = [];

ggd.grid.grid_subset{1}.metric.jacobian          = [];
ggd.grid.grid_subset{1}.metric.g11_covariant     = [];
ggd.grid.grid_subset{1}.metric.g12_covariant     = [];
ggd.grid.grid_subset{1}.metric.g13_covariant     = [];
ggd.grid.grid_subset{1}.metric.g21_covariant     = [];
ggd.grid.grid_subset{1}.metric.g22_covariant     = [];
ggd.grid.grid_subset{1}.metric.g23_covariant     = [];
ggd.grid.grid_subset{1}.metric.g31_covariant     = [];
ggd.grid.grid_subset{1}.metric.g32_covariant     = [];
ggd.grid.grid_subset{1}.metric.g33_covariant     = [];
ggd.grid.grid_subset{1}.metric.g11_contravariant = [];
ggd.grid.grid_subset{1}.metric.g12_contravariant = [];
ggd.grid.grid_subset{1}.metric.g13_contravariant = [];
ggd.grid.grid_subset{1}.metric.g21_contravariant = [];
ggd.grid.grid_subset{1}.metric.g22_contravariant = [];
ggd.grid.grid_subset{1}.metric.g23_contravariant = [];
ggd.grid.grid_subset{1}.metric.g31_contravariant = [];
ggd.grid.grid_subset{1}.metric.g32_contravariant = [];
ggd.grid.grid_subset{1}.metric.g33_contravariant = [];

ggd.r{1}.grid_index                     = -999999999;
ggd.r{1}.grid_subset_index              = -999999999;
ggd.r{1}.values                         = reshape(p2d.r.',n1 * n2,1);
ggd.z{1}.grid_index                     = -999999999;
ggd.z{1}.grid_subset_index              = -999999999;
ggd.z{1}.values                         = reshape(p2d.z.',n1 * n2,1);
ggd.theta{1}.grid_index                 = -999999999;
ggd.theta{1}.grid_subset_index          = -999999999;
if isfield(p2d,'theta')
  ggd.theta{1}.values                   = reshape(p2d.theta.',n1 * n2,1);;
else
  ggd.theta{1}.values                   = [];
end
ggd.psi{1}.grid_index                   = -999999999;
ggd.psi{1}.grid_subset_index            = -999999999;
ggd.psi{1}.values                       = reshape(p2d.psi.',n1 * n2,1);
ggd.phi{1}.grid_index                   = -999999999;
ggd.phi{1}.grid_subset_index            = -999999999;
if isfield(p2d,'phi')
  ggd.phi{1}.values                     = reshape(p2d.phi.',n1 * n2,1);;
else
  ggd.phi{1}.values                     = [];
end
ggd.j_tor{1}.grid_index                 = -999999999;
ggd.j_tor{1}.grid_subset_index          = -999999999;
if isfield(p2d,'j_tor')
  ggd.j_tor{1}.values                   = reshape(p2d.j_tor.',n1 * n2,1);;
else
  ggd.j_tor{1}.values                   = [];
end
ggd.j_parallel{1}.grid_index            = -999999999;
ggd.j_parallel{1}.grid_subset_index     = -999999999;
if isfield(p2d,'j_parallel')
  ggd.j_parallel{1} .values             = reshape(p2d.j_parallel.',n1 * n2,1);;
else
  ggd.j_parallel{1}.values              = [];
end
ggd.b_field_r{1}.grid_index             = -999999999;
ggd.b_field_r{1}.grid_subset_index      = -999999999;
ggd.b_field_r{1}.values                 = reshape(p2d.b_field_r.',n1 * n2,1);
ggd.b_field_z{1}.grid_index             = -999999999;
ggd.b_field_z{1}.grid_subset_index      = -999999999;
ggd.b_field_z{1}.values                 = reshape(p2d.b_field_z.',n1 * n2,1);
ggd.b_field_tor{1}.grid_index           = -999999999;
ggd.b_field_tor{1}.grid_subset_index    = -999999999;
ggd.b_field_tor{1}.values               = reshape(p2d.b_field_tor.',n1 * n2,1);

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
  psi_out_mem = psi_out;
  BR_out_mem  = BR_out;
  BZ_out_mem  = BZ_out;
  psi_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.PSI(:),Rin(mask_out),Zin(mask_out),'cubic');
  BR_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.BR(:),Rin(mask_out),Zin(mask_out),'cubic');
  BZ_out(mask_out) = griddata(equivide.R(:),equivide.Z(:),equivide.BZ(:),Rin(mask_out),Zin(mask_out),'cubic');
  warning on
  %
  psi_out(~isfinite(psi_out)) = psi_out_mem(~isfinite(psi_out));
  BR_out(~isfinite(BR_out))   = BR_out_mem(~isfinite(BR_out));
  BZ_out(~isfinite(BZ_out))   = BZ_out_mem(~isfinite(BZ_out));
  %
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

% fonction appele par zfitvide_interp
function [psi_out,BR_out,BZ_out] = zfit_poly_interp(equivide,varargin)
  psi_out = zpsi_poly_interp(equivide);

% fonction appele par zfitvide_interp
function [psi_out,BR_out,BZ_out] = zpsi_poly_interp(equivide,Rin,Zin,fit)


if nargin == 1
    
  equivide.R     = squeeze(equivide.R);
  equivide.Z     = squeeze(equivide.Z);
  equivide.BR    = squeeze(equivide.BR);
  equivide.BZ    = squeeze(equivide.BZ);
  equivide.PSI   = equivide.psiRZ' * ones(1,size(equivide.R,2));
  %
  [equivide.F_PSI,equivide.F_BR,equivide.F_BZ] = extrapolate_from_LCFS(equivide.R(end,:),equivide.Z(end,:),equivide.PSI(end,:),-equivide.BR(end,:),-equivide.BZ(end,:));
  
 
  
  % traitement de l discontinuite
  % recalcul 
  psis =  equivide.F_PSI(equivide.R(end,:),equivide.Z(end,:));
  equivide.psi_lcfs_delta    = equivide.PSI(end,:) - psis;
  equivide.dpsi_error = mean(abs(equivide.psi_lcfs_delta)) ./ abs( equivide.PSI(1,1) - equivide.PSI(end,1));
  equivide.BR_lcfs_delta     = equivide.BR(end,:) + equivide.F_BR(equivide.R(end,:),equivide.Z(end,:));
  equivide.BZ_lcfs_delta     = equivide.BZ(end,:) + equivide.F_BZ(equivide.R(end,:),equivide.Z(end,:));
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
  mask_out = zinout(equivide.R(end,:),equivide.Z(end,:),Rin,Zin,0);
  warning off
  %psi_out =  2 .* equivide.PSI(end,1) - equivide.F_PSI(Rin,Zin);
  psi_out =     equivide.F_PSI(Rin,Zin);
  BR_out  =   - equivide.F_BR(Rin,Zin);
  BZ_out  =   - equivide.F_BZ(Rin,Zin);
  psi_out_mem = psi_out;
  BR_out_mem  = BR_out;
  BZ_out_mem  = BZ_out;
  F_PSI_EXT   = scatteredInterpolant(equivide.R(:),equivide.Z(:),equivide.PSI(:),'natural','linear');
  F_BR_EXT    = scatteredInterpolant(equivide.R(:),equivide.Z(:),equivide.BR(:),'natural','linear');
  F_BZ_EXT    = scatteredInterpolant(equivide.R(:),equivide.Z(:),equivide.BZ(:),'natural','linear');
  psi_out(mask_out) = F_PSI_EXT(Rin(mask_out),Zin(mask_out));
  BR_out(mask_out)  = F_BR_EXT(Rin(mask_out),Zin(mask_out));
  BZ_out(mask_out)  = F_BZ_EXT(Rin(mask_out),Zin(mask_out));
  warning on
  %
  psi_out(~isfinite(psi_out)) = psi_out_mem(~isfinite(psi_out));
  BR_out(~isfinite(BR_out))   = BR_out_mem(~isfinite(BR_out));
  BZ_out(~isfinite(BZ_out))   = BZ_out_mem(~isfinite(BZ_out));
  %
  psi_out = reshape(psi_out,nsin);
  BR_out  = reshape(BR_out,nsin);
  BZ_out  = reshape(BZ_out,nsin);
  Rin  = reshape(Rin,nsin);
  Zin  = reshape(Zin,nsin);
  
  if   false
      figure(31)
      clf
      contour(equivide.R,equivide.Z,equivide.PSI,equivide.PSI(:,1),'color','b','linestyle','--')
      hold on
      contour(Rin,Zin,psi_out,equivide.PSI(:,1),'color','r','linestyle','-')
      lpsi = cat(1,equivide.PSI(:,1),linspace(equivide.PSI(1,1),equivide.PSI(end,1),11).'+22/21 .* (equivide.PSI(end,1) - equivide.PSI(1,1)));
      contour(Rin,Zin,psi_out,lpsi,'color','g','linestyle','-')
      quiver(equivide.R,equivide.Z,equivide.BR,equivide.BZ,'color','c');
      quiver(Rin,Zin,BR_out,BZ_out,'color','m');
      title('Bpoloidal: blue = METIS & red = Extrapolation');
      xlabel('R (m)');
      ylabel('Z (m)');
      plot(equivide.R(end,:),equivide.Z(end,:),'.k');
      drawnow
      figure(33)
      clf
      contour(equivide.R,equivide.Z,equivide.PSI,equivide.PSI(:,1),'color','b','linestyle','--')
      hold on
      contour(Rin,Zin,psi_out,equivide.PSI(:,1),'color','r','linestyle','-')
      contour(Rin,Zin,psi_out,lpsi)
      title('Bpoloidal: blue = METIS & red = Extrapolation');
      xlabel('R (m)');
      ylabel('Z (m)');
      plot(equivide.R(end,:),equivide.Z(end,:),'.k');
      drawnow
  end
  if any(~isfinite(psi_out(:))) || any(~isfinite(BR_out(:))) || any(~isfinite(BZ_out(:)))
      keyboard
  end
end
