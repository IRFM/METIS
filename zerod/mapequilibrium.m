% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function equilibrium = mapequilibrium(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,equilibrium,scenario,opt_grille,sigma_B0_eff,equi_extrap,factor_two_pi)


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
    equi_extrap = 1;
elseif isempty(equi_extrap)
    equi_extrap = 1;
end 

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

% heritage de scenario
equilibrium.datainfo = scenario.datainfo;
% temps
equilibrium.time			= profil0d.temps;


% structure eqgeometry
equilibrium.eqgeometry.source              = 'METIS';
xpoint = interp1_itm(data_zerod.temps,data_zerod.xpoint,profil0d.temps,'nearest','extrap');
%  for k=1:length(profil0d.temps)
%   	if xpoint(k) == 1
%   		equilibrium.eqgeometry.boundarytype{k}  = 'separatrix';
%   	else 
%   		equilibrium.eqgeometry.boundarytype{k}  = 'limiter';
%   	end
%  end
equilibrium.eqgeometry.boundarytype = xpoint;
%
if isfield(profil0d,'Rsepa') && isfield(profil0d,'Zsepa')
	equilibrium.eqgeometry.boundary.r 		= profil0d.Rsepa;
	equilibrium.eqgeometry.boundary.z 		= profil0d.Zsepa + interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.z0,profil0d.temps,'nearest','extrap') * ...
	                                                  ones(1,size(profil0d.Zsepa,2));
else
	equilibrium.eqgeometry.boundary.r 		= [];
	equilibrium.eqgeometry.boundary.z 		= [];

end
equilibrium.eqgeometry.geom_axis.r 		= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.R,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.geom_axis.z 		= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.z0,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.a_minor 	 		= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.a,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.elongation 	 	= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.K,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.tria_upper 	 	= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.d,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.tria_lower 	 	= interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.d,profil0d.temps,'nearest','extrap');
equilibrium.eqgeometry.xpts.r 			= [];
equilibrium.eqgeometry.xpts.z 			= [];
equilibrium.eqgeometry.left_low_st.r 		= [];
equilibrium.eqgeometry.left_low_st.z 		= [];
equilibrium.eqgeometry.right_low_st.r 		= [];
equilibrium.eqgeometry.right_low_st.z 		= [];
equilibrium.eqgeometry.left_up_st.r 		= []; 
equilibrium.eqgeometry.left_up_st.z 		= []; 
equilibrium.eqgeometry.right_up_st.r 		= [];
equilibrium.eqgeometry.right_up_st.z 		= [];
equilibrium.eqgeometry.active_limit.r 		= [];  
equilibrium.eqgeometry.active_limit.z 		= [];

% parametres scalaires
equilibrium.global_param.beta_pol       = interp1_itm(data_zerod.temps,scenario.global_param.beta_pol.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.beta_tor       = interp1_itm(data_zerod.temps,scenario.global_param.beta_tor.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.beta_normal    = interp1_itm(data_zerod.temps,scenario.global_param.beta_normal.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.i_plasma       = interp1_itm(data_zerod.temps,scenario.global_param.ip.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.li             = interp1_itm(data_zerod.temps,scenario.global_param.li.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.volume         = interp1_itm(data_zerod.temps,scenario.global_param.volume.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.area           = interp1_itm(data_zerod.temps,scenario.global_param.area_pol.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.psi_ax         = interp1_itm(data_zerod.temps,scenario.centre.psi0.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.psi_bound      = interp1_itm(data_zerod.temps,scenario.edge.psi_edge.value,profil0d.temps,'nearest','extrap');
%
equilibrium.global_param.mag_axis.bphi  = profil0d.fdia(:,1) ./ profil0d.Raxe(:,1);
equilibrium.global_param.mag_axis.q     = interp1_itm(data_zerod.temps,scenario.centre.q0.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.mag_axis.position.r   = interp1_itm(data_zerod.temps,scenario.centre.Rmag.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.mag_axis.position.z   = interp1_itm(data_zerod.temps,scenario.centre.Zmag.value,profil0d.temps,'nearest','extrap');
%
equilibrium.global_param.q_95           = interp1_itm(data_zerod.temps,scenario.ninety_five.q_95.value,profil0d.temps,'nearest','extrap');
equilibrium.global_param.q_min           = min(profil0d.qjli,[],2);
%
rb0 =   interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
equilibrium.global_param.toroid_field.r0          = mean(z0dstruct.z0dinput.geo.R);
equilibrium.global_param.toroid_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ equilibrium.global_param.toroid_field.r0;
%
equilibrium.global_param.w_mhd           = interp1_itm(data_zerod.temps,data_zerod.w,profil0d.temps,'nearest','extrap');

% pour le nombre de mach
ve   = ones(1,size(profil0d.tip,2));
meff = interp1_itm(data_zerod.temps,data_zerod.meff,profil0d.temps,'nearest','extrap');

% reservation de la place memoire
v1 = NaN .* ones(size(profil0d.psi));
equilibrium.profiles_1d.psi 		= v1;
equilibrium.profiles_1d.phi 		= v1;
equilibrium.profiles_1d.pressure 	= v1;
equilibrium.profiles_1d.F_dia 		= v1;
equilibrium.profiles_1d.pprime 		= v1;
equilibrium.profiles_1d.ffprime 	= v1;
equilibrium.profiles_1d.jphi 		= v1;
equilibrium.profiles_1d.jparallel 	= v1;
equilibrium.profiles_1d.q 		= v1;
equilibrium.profiles_1d.r_inboard 	= v1;
equilibrium.profiles_1d.r_outboard 	= v1;
equilibrium.profiles_1d.beta_pol 	= v1;
equilibrium.profiles_1d.li      	= v1;
equilibrium.profiles_1d.rho_tor 	= v1;
equilibrium.profiles_1d.dpsidrho_tor = v1;
equilibrium.profiles_1d.rho_vol 	= v1;
equilibrium.profiles_1d.elongation 	= v1;
equilibrium.profiles_1d.tria_upper 	= v1;
equilibrium.profiles_1d.tria_lower 	= v1;
equilibrium.profiles_1d.volume 		= v1;
equilibrium.profiles_1d.vprime 		= v1;
equilibrium.profiles_1d.area 		= v1;
equilibrium.profiles_1d.aprime 		= v1;
equilibrium.profiles_1d.surface 	= v1;
equilibrium.profiles_1d.ftrap 		= v1;
equilibrium.profiles_1d.gm1 		= v1;
equilibrium.profiles_1d.gm2 		= v1;
equilibrium.profiles_1d.gm3 		= v1;
equilibrium.profiles_1d.gm4 		= v1;
equilibrium.profiles_1d.gm5 		= v1;
equilibrium.profiles_1d.gm6 		= v1;
equilibrium.profiles_1d.gm7 		= v1;
equilibrium.profiles_1d.gm8 		= v1;
equilibrium.profiles_1d.gm9 		= v1;
equilibrium.profiles_1d.b_av 		= v1;
equilibrium.profiles_1d.b_min 		= v1;
equilibrium.profiles_1d.b_max 		= v1;
equilibrium.profiles_1d.omega 		= v1;
equilibrium.profiles_1d.omegaprime  = v1;
equilibrium.profiles_1d.mach_a      = v1;
%
if opt_grille == 0
      equilibrium.profiles_2d.grid.dim1       = []; 	
      equilibrium.profiles_2d.grid.dim2       = [];
      equilibrium.profiles_2d.grid.connect    = [];
      equilibrium.profiles_2d.grid_type       = char({'2','inverse','3','polar'});
      equilibrium.profiles_2d.psi  	      = [];
      equilibrium.profiles_2d.jphi  	      = [];
      equilibrium.profiles_2d.jpar  	      = [];
      equilibrium.profiles_2d.br 		= [];
      equilibrium.profiles_2d.bz 		= [];
      equilibrium.profiles_2d.bphi 		= [];
      equilibrium.profiles_2d.r                 = [];
      equilibrium.profiles_2d.z                 = [];
      equilibrium.profiles_2d.phi               = [];
      equilibrium.profiles_2d.theta             = [];
      equilibrium.profiles_2d.pressure          = [];
     
else
      v3  = NaN .* ones(size(profil0d.psi,1),size(profil0d.psi,2),nbth);
      equilibrium.profiles_2d.grid.dim1       = []; 	
      equilibrium.profiles_2d.grid.dim2       = [];
      equilibrium.profiles_2d.grid.connect    = [];
      equilibrium.profiles_2d.grid_type       = char({'3','irregular','4','unstruct'});
      equilibrium.profiles_2d.psi    	      = v3;
      equilibrium.profiles_2d.jphi            = v3;
      equilibrium.profiles_2d.jpar            = v3;
      equilibrium.profiles_2d.br 	      = v3;
      equilibrium.profiles_2d.bz 	      = v3;
      equilibrium.profiles_2d.bphi 	      = v3;
      equilibrium.profiles_2d.r               = v3;
      equilibrium.profiles_2d.z               = v3;
      equilibrium.profiles_2d.phi             = v3;
      equilibrium.profiles_2d.theta           = v3;
      equilibrium.profiles_2d.pressure        = v3;

end
%
v3  = NaN .* ones(size(profil0d.psi,1),size(profil0d.psi,2),nbth);
equilibrium.coord_sys.grid_type  = [];
equilibrium.coord_sys.grid.dim1  = NaN .* ones(size(profil0d.psi));
equilibrium.coord_sys.grid.dim2  = NaN .* ones(size(profil0d.psi,1),nbth);
equilibrium.coord_sys.jacobian   = v3;
equilibrium.coord_sys.g_11   	 = v3;
equilibrium.coord_sys.g_12   	 = v3;
equilibrium.coord_sys.g_13   	 = v3;
equilibrium.coord_sys.g_22   	 = v3;
equilibrium.coord_sys.g_23   	 = v3;
equilibrium.coord_sys.g_33   	 = v3;
equilibrium.coord_sys.position.r = v3;
equilibrium.coord_sys.position.z = v3;


% calcul sur une grille rectangulaire
profiles_2d_2.grid.dim1       = []; 	
profiles_2d_2.grid.dim2       = [];
profiles_2d_2.grid.connect    = [];
profiles_2d_2.grid_type       = char({'1','rectangular','0','none'});
profiles_2d_2.psi	      = [];
profiles_2d_2.jphi            = [];
profiles_2d_2.jpar            = [];
profiles_2d_2.br 	      = [];
profiles_2d_2.bz 	      = [];
profiles_2d_2.bphi 	      = [];
profiles_2d_2.r               = [];
profiles_2d_2.z               = [];      
profiles_2d_2.phi             = [];
profiles_2d_2.theta           = [];
profiles_2d_2.pressure        = [];

%


% c'est extremment lent dans matlab
%	[equilibrium.coord_sys.position(1:size(profil0d.psi,1),1:size(profil0d.psi,2),1:nbth).R] = deal(NaN);
%	[equilibrium.coord_sys.position(1:size(profil0d.psi,1),1:size(profil0d.psi,2),1:nbth).Z] = deal(NaN);
% la grille rectangulaire doit avoir un nombre de point fixe 
rzp = max(1,max(z0dstruct.z0dinput.geo.K(:)));

% boucle de calcul des sorties manquantes
for k =1: length(equilibrium.time)

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d/%d\t',k,length(equilibrium.time));

    % donnees de metis 
    prof.x    = profil0d.xli;
    prof.kx   = profil0d.kx(k,:);     
    prof.dx   = profil0d.dx(k,:);      
    prof.rmx  = profil0d.rmx(k,:);     
    prof.Raxe = profil0d.Raxe(k,:);
    prof.psi  = profil0d.psi(k,:);
    prof.fdia = profil0d.fdia(k,:);
    prof.jmoy = profil0d.jli(k,:);
    prof.ptot = profil0d.ptot(k,:);
    %
    geo_input.a     = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.a ,profil0d.temps(k),'pchip','extrap');
    geo_input.R     = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps(k),'pchip','extrap');
    geo_input.K     = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.K ,profil0d.temps(k),'pchip','extrap');
    geo_input.d     = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.d ,profil0d.temps(k),'pchip','extrap');
    geo_input.b0    = z0dstruct.z0dinput.option.signe .* interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.b0 ,profil0d.temps(k),'pchip','extrap'); 
    z0_offset       = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.z0 ,profil0d.temps(k),'pchip','extrap');
    geo_input.sp    = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sp ,profil0d.temps(k),'pchip','extrap');
    geo_input.vp    = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.vp ,profil0d.temps(k),'pchip','extrap');
    geo_input.sext  = interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sext ,profil0d.temps(k),'pchip','extrap');
    if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa') && all(isfinite(profil0d.Rsepa(k,:))) && all(isfinite(profil0d.Zsepa(k,:)))
      geo_input.Rsepa = profil0d.Rsepa(k,:);
      geo_input.Zsepa = profil0d.Zsepa(k,:);
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
	      geo_input.Rsepa = interp1(KH,Rsepa,index_full,'linear',NaN);
	      geo_input.Zsepa = interp1(KH,Zsepa,index_full,'linear',NaN);
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
      geo_input.z0    = (max(profil0d.Zsepa(k,:)) + min((profil0d.Zsepa(k,:)))) ./ 2;
      geo_input.Zsepa = geo_input.Zsepa - geo_input.z0; 
      box.a    = max(z0dstruct.z0dinput.geo.a);
      box.z0   = sum(z0dstruct.z0dinput.geo.z0 .* z0dstruct.z0dinput.cons.ip) ./ sum(z0dstruct.z0dinput.cons.ip);
      box.rmin = min(profil0d.Rsepa(:)) - box.a/2;
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
      box.rmin = min(z0dstruct.z0dinput.geo.R - z0dstruct.z0dinput.geo.a) - box.a/2;
      box.rmax = max(z0dstruct.z0dinput.geo.R + z0dstruct.z0dinput.geo.a) + box.a/2;
      box.zmin = min(- z0dstruct.z0dinput.geo.a .* z0dstruct.z0dinput.geo.K) - box.a/2*rzp + box.z0;
      box.zmax = max(z0dstruct.z0dinput.geo.a .* z0dstruct.z0dinput.geo.K) + box.a/2*rzp + box.z0;
    end
    %
    [profil,deuxd,moment,scalaire,facteur,Ip_lcfs] = metis2equi1t(prof,geo_input,phys,nbth,1, ...
                                             z0dstruct.z0dinput.option.morphing,0,factor_two_pi);
    % there is a problem with z0
    deuxd.Z     = deuxd.Z + geo_input.z0;
    moment.zaxe = moment.zaxe  +  geo_input.z0;
    Iprof = cumtrapz(profil0d.xli,profil0d.spr(k,:) .* profil0d.jli(k,:),2);
    Iprof_bis = profil0d.rmx(k,end) .* cumtrapz(profil0d.xli,profil0d.vpr_tor(k,:) .* profil0d.jli(k,:) ./ 2 ./ pi .* profil0d.ri(k,:),2);
%      figure(118);clf
%      plot(prof.x,Ip_lcfs,'b',profil0d.xli,Iprof,'r',1,equilibrium.global_param.i_plasma(k),'ok',profil0d.xli,Iprof_bis,'g');
%      drawnow
    error_2d = abs(equilibrium.global_param.i_plasma(k) - Ip_lcfs(end)) ./ max(1,abs(equilibrium.global_param.i_plasma(k)));
    fprintf(' error_2D = %g | ',error_2d);
    % rescale to preserve plasma current even with high pedestal
    Iprof  = abs(Iprof ./ Iprof(end) .* equilibrium.global_param.i_plasma(k));
    factor = abs(Iprof ./ max(1,abs(Ip_lcfs)));
    factor(1) = 1;
    %figure(119);clf;plot(factor);drawnow
    vth            = ones(1,size(deuxd.dPSIdx,2));
    deuxd.dPSIdx   = (factor' * vth) .* deuxd.dPSIdx;
    deuxd.BR       = (factor' * vth) .* deuxd.BR;
    deuxd.BZ       = (factor' * vth) .* deuxd.BZ;
    deuxd.dPSIdR   = (factor' * vth) .* deuxd.dPSIdR; 
    deuxd.dPSIdZ   = (factor' * vth) .* deuxd.dPSIdZ;
    %
    deuxd.Z     = deuxd.Z + z0_offset;
    moment.zaxe = moment.zaxe   + z0_offset;

    % remplissage de la structure pour le precalcul de la grille rectangulaire
    % remplissage de la structure equivide
    % passage en coodornnees de flux
    xout = profil.rho ./ max(profil.rho);
	
    % donnees scalaire
    equicronos.errit          = 0;               % pas de convergence
    equicronos.fail           = 0;               % fonctionne toujours
    equicronos.amix           = NaN;                 % final value of convergence mixing parameter
    equicronos.bnorme         = NaN;
    equicronos.psi0           = profil.psi(1) ./ factor_two_pi;
    equicronos.ip             = scalaire.ipout;              % plasma current from equilibrium (A)
    equicronos.li             = scalaire.liout;              % internal inductance
    equicronos.betap          = scalaire.betap;              % poloidal normalised pressure
    equicronos.betat          = scalaire.betat;              % toroidal normalised pressure
    equicronos.betan          = scalaire.betan;              % current normalised pressure
	
     % calcul de rhomax
     equicronos.rhomax          = profil.rho(end);
     equicronos.q(1,:)          = profil.q;
     equicronos.q0              = equicronos.q(1,1);          
     % calcul du shear
     
     equicronos.shear(1,:)           = pdederive(xout,profil.q,0,2,2,1) ./ profil.q .* xout;
     equicronos.raxe(1,:)            = moment.raxe';  
     equicronos.zaxe(1,:)            = moment.zaxe' - moment.zaxe(1);;  
     equicronos.d(1,:)               = equicronos.raxe(1,:) - equicronos.raxe(1,end);
     equicronos.a(1,:)               = moment.a';
     equicronos.rhog(1,:)            = equicronos.a(1,:) ./ equicronos.a(1,end);
     equicronos.e(1,:)               = moment.e';
     equicronos.trh(1,:)             = moment.trh';
     equicronos.trl(1,:)             = moment.trl';
     equicronos.indh(1,:)            = NaN .* xout;
     equicronos.indl(1,:)            = NaN .* xout;
     equicronos.b2i(1,:)             = profil.b2i;  
     equicronos.b2(1,:)              = profil.b2;
     equicronos.vpr(1,:)             = profil.vpr;  
     equicronos.volume(1,:)          = profil.volume; 
     equicronos.spr(1,:)             = profil.spr; 
     equicronos.surface(1,:)         = profil.surface;
     equicronos.grho2r2(1,:)         = profil.grho2r2;
     equicronos.grhor(1,:)           = profil.grhor;        
     equicronos.ri(1,:)              = profil.ri;
     equicronos.grho2(1,:)           = profil.grho2;
     equicronos.c2c(1,:)             = profil.C2;
     equicronos.grho(1,:)            = profil.grho;
     equicronos.rmoy(1,:)            = profil.rmoy;
     equicronos.jmoy(1,:)            = profil.jmoy;
     equicronos.ptot(1,:)            = profil.ptot;
     equicronos.r2(1,:)              = profil.r2;
     equicronos.r2i(1,:)             = profil.r2i;
     equicronos.grho2b2(1,:)         = profil.grho2b2;    
     equicronos.errcur    	       = 0;                 % final relative error on current
	
     % utilisation des donnnees du mapping pour les indice de temps problematique
     equicronos.ftrap(1,:)           = profil.ftrap;
     equicronos.phi(1,:)             = profil.phi; 
     equicronos.vpr(1,:)             = profil.vpr; 
     equicronos.spr(1,:)             = profil.spr;
     equicronos.c2c(1,:)             = profil.C2;  
 
     % plus utiliser dans cronos
     % equi.r3tau3         = NaN .* equi.r2tau2;
     % equi.r3tau          = NaN .* equi.r3tau; 
     % equi.r2tau2         = NaN .* equi.r2tau2;
		
    % le R et le Z
    equicronos.R(1,:,:)          = shiftdim(deuxd.R,-1);
    equicronos.Z(1,:,:)          = shiftdim(deuxd.Z,-1) - moment.zaxe(1);;
    equicronos.rhoRZ(1,:)        = profil.rho';
    equicronos.psiRZ(1,:)        = profil.psi' ./ factor_two_pi;
    equicronos.df2RZ(1,:)        = (2 .* profil.fdia .* pdederive(profil.psi,profil.fdia,2,2,2,1))';
    equicronos.dprRZ(1,:)        = pdederive(profil.psi,profil.ptot,2,2,2,1)';
    nb_frmode         = NaN;
    equicronos.frmode(1,:)       = NaN;
    % La carte de champ
    equicronos.BR(1,:,:)     = shiftdim(deuxd.BR,-1) .* sign(factor_two_pi);
    equicronos.BZ(1,:,:)     = shiftdim(deuxd.BZ,-1) .* sign(factor_two_pi);
    equicronos.BPHI(1,:,:)   = shiftdim(deuxd.BPHI,-1);


    % remplissage des profils;
    equilibrium.profiles_1d.psi(k,:) 		= profil0d.psi(k,:);
    equilibrium.profiles_1d.phi(k,:) 		= profil.phi;
    equilibrium.profiles_1d.pressure(k,:) 	= profil0d.ptot(k,:);
    equilibrium.profiles_1d.F_dia(k,:) 		= profil0d.fdia(k,:);
    equilibrium.profiles_1d.pprime(k,:) 	= pdederive(profil0d.psi(k,:),profil0d.ptot(k,:),0,2,2,1);
    equilibrium.profiles_1d.ffprime(k,:) 	= pdederive(profil0d.psi(k,:),profil0d.fdia(k,:) .^ 2,0,2,2,1) ./ 2 ;
    equilibrium.profiles_1d.jphi(k,:) 		= profil0d.jli(k,:);
    equilibrium.profiles_1d.jparallel(k,:) 	= profil0d.jeff(k,:);
    equilibrium.profiles_1d.q(k,:) 		= profil0d.qjli(k,:);

    equilibrium.profiles_1d.r_inboard(k,:) 	= min(deuxd.R,[],2);
    equilibrium.profiles_1d.r_outboard(k,:) = max(deuxd.R,[],2);

    equilibrium.profiles_1d.li(k,:)      	= 2 .* cumtrapz(profil.rho,profil.bpol .^ 2 .* profil.vpr,2) ./  ...
                                                  ((phys.mu0 .* max(1,cumtrapz(profil.rho,profil.jmoy .* profil.spr,2)) ) .^ 2  .* profil.Raxe);
    equilibrium.profiles_1d.beta_pol(k,:)	= cumtrapz(profil.rho,profil.ptot .* profil.spr,2) ./  ...
						  max(eps,cumtrapz(profil.rho,profil.spr,2)) ./ max(eps,profil.bpol .^2 ./ 2 ./ phys.mu0);   
    equilibrium.profiles_1d.rho_tor(k,:) 	= profil0d.rmx(k,:);
    equilibrium.profiles_1d.rho_vol(k,:) 	= sqrt(profil.volume ./max(profil.volume));
    equilibrium.profiles_1d.dpsidrho_tor(k,:) 	= pdederive(profil0d.rmx(k,:),profil0d.psi(k,:),0,2,2,1);
    equilibrium.profiles_1d.elongation(k,:) 	= moment.e;
    equilibrium.profiles_1d.tria_upper(k,:) 	= moment.trh;
    equilibrium.profiles_1d.tria_lower(k,:) 	= moment.trl;
    equilibrium.profiles_1d.volume(k,:) 	= profil.volume;
    % equilibrium.profiles_1d.vprime(k,:) 	= profil.vpr;
    equilibrium.profiles_1d.vprime(k,:) 	= pdederive(profil.psi,profil.volume,0,2,2,1);
    equilibrium.profiles_1d.vprime(k,:) 	= equilibrium.profiles_1d.vprime(k,:) .* profil.vpr ./ ...
                                              (equilibrium.profiles_1d.vprime(k,end) .* equilibrium.profiles_1d.dpsidrho_tor(k,end));
    equilibrium.profiles_1d.area(k,:) 		= profil.surface;
    %equilibrium.profiles_1d.aprime(k,:) 	= profil.spr;
    equilibrium.profiles_1d.aprime(k,:) 	= pdederive(profil.psi,profil.surface,0,2,2,1);
    equilibrium.profiles_1d.aprime(k,:) 	= equilibrium.profiles_1d.aprime(k,:) .* profil.spr ./ ...
                                              (equilibrium.profiles_1d.aprime(k,end) .* equilibrium.profiles_1d.dpsidrho_tor(k,end));
    equilibrium.profiles_1d.area(k,:) 		= profil.surface;
    sext =  interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.zerod.sext ,profil0d.temps(k),'pchip','extrap');
    equilibrium.profiles_1d.surface(k,:) 	= profil.vpr .* profil.grho; 
    equilibrium.profiles_1d.surface(k,:)    = (sext ./ equilibrium.profiles_1d.surface(k,end)) .*  equilibrium.profiles_1d.surface(k,:);  
    equilibrium.profiles_1d.ftrap(k,:) 		= profil0d.ftrap(k,:);
    equilibrium.profiles_1d.gm1(k,:) 		= profil.r2i;
    equilibrium.profiles_1d.gm2(k,:) 		= profil.grho2r2;
    equilibrium.profiles_1d.gm3(k,:) 		= profil.grho2;
    equilibrium.profiles_1d.gm4(k,:) 		= profil.b2i;
    equilibrium.profiles_1d.gm5(k,:) 		= profil.b2;
    equilibrium.profiles_1d.gm6(k,:) 		= profil.grho2b2;
    equilibrium.profiles_1d.gm7(k,:)		= profil.grho;
    equilibrium.profiles_1d.gm8(k,:) 		= profil.rmoy;
    equilibrium.profiles_1d.gm9(k,:) 		= profil.ri;
    equilibrium.profiles_1d.b_av(k,:) 		= sqrt(profil.b2);
    
    bmap = sqrt(deuxd.BR .^ 2 +deuxd.BZ .^ 2 + ((profil.fdia' * ones(1,size(deuxd.BR,2))) ./ deuxd.R) .^ 2);
    b_min = min(bmap,[],2);
    b_max = max(bmap,[],2);
    equilibrium.profiles_1d.b_min(k,:) 		= b_min;
    equilibrium.profiles_1d.b_max(k,:) 		= b_max;
    
    equilibrium.profiles_1d.omega(k,:)      = profil0d.omega(k,:);
    equilibrium.profiles_1d.omegaprime(k,:) = pdederive(profil0d.psi(k,:),profil0d.omega(k,:),0,2,2,1);
    
    % calcul du nombre de mach
    rhoi = 4.57e-3 .* sqrt(meff(k)) .* sqrt(profil0d.tip(k,:) ./ 1e3) ./ sqrt(profil.b2);
    valf = sqrt(profil.b2) ./ sqrt(phys.mu0 .* rhoi);
    equilibrium.profiles_1d.mach_a(k,:) = sqrt(profil0d.vtor(k,:) .^ 2  + profil0d.vtheta(k,:) .^ 2) ./ valf;
    
    equilibrium.profiles_1d.phi_flow(k,:) =  profil0d.rmx(k,:) .*  profil0d.vtheta(k,:) .* profil0d.bpol(k,:);
        
    % moment pour coherence interne
    equilibrium.eqgeometry.geom_axis.r(k)	= profil0d.Raxe(k,1);
    equilibrium.eqgeometry.geom_axis.z(k) 	= moment.zaxe(1);
    equilibrium.eqgeometry.a_minor(k) 	 	= moment.a(end);
    equilibrium.eqgeometry.elongation(k) 	= moment.e(end);
    equilibrium.eqgeometry.tria_upper(k) 	= moment.trh(end);
    equilibrium.eqgeometry.tria_lower(k) 	= moment.trl(end); 
  
    % calcul de Jphi et jpar
    pprim  = equilibrium.profiles_1d.pprime(k,:)'  * ones(1,size(deuxd.R,2));
    ffprim = equilibrium.profiles_1d.ffprime(k,:)' * ones(1,size(deuxd.R,2));
    pressure = equilibrium.profiles_1d.pressure(k,:)' * ones(1,size(deuxd.R,2)); 
    phi      = equilibrium.profiles_1d.phi(k,:)' * ones(1,size(deuxd.R,2)); 
    jphi  = pprim .* deuxd.R + ffprim ./ deuxd.R ./ phys.mu0; 
    jpar  = jphi .* deuxd.BPHI + ffprim ./ phys.mu0 ./ (equilibrium.profiles_1d.F_dia(k,:)' * ones(1,size(deuxd.R,2))) .*  ...
            (deuxd.BR .^ 2 +deuxd.BZ .^ 2);
    vphi   = (profil0d.omega(k,:)' * ones(1,size(deuxd.R,2))) .* deuxd.R;
    bpol2d = abs(deuxd.BR + sqrt(-1) .*  deuxd.BZ);
    vtheta =  (profil0d.utheta(k,:)' * ones(1,size(bpol2d,2))) .* bpol2d;
    b0    = equilibrium.profiles_1d.F_dia(k,end) ./ profil0d.Raxe(k,end);
    jpar  = jpar ./ b0;
    RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
    ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
    theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
    % selon le choxi de grille
    if opt_grille == 1
        % preparatiopn des donnees
        RZ = deuxd.R(:) + sqrt(-1) .*  deuxd.Z(:);
        [RZ,ind] = unique(RZ);
        BP           = deuxd.BR(:) + sqrt(-1) .*  deuxd.BZ(:);
        BP           = BP(ind);
        BPHI         = deuxd.BPHI(:);
        BPHI         = BPHI(ind);
        PSI          = deuxd.PSI(:);
        PSI          = PSI(ind);
        %PHI          = profil.phi' * ones(1,size(deuxd.R,2));
        %PHI          = PHI(:);
        %PHI          = PHI(ind);
        jphi 	     = jphi(:);
        jphi         = jphi(ind);
        jpar 	     = jpar(:);
        jpar         = jpar(ind);
        vphi 	     = vphi(:);
        vphi         = vphi(ind);
        vtheta 	     = vtheta(:);
        vtheta       = vtheta(ind);
        phi 	     = phi(:);
        phi          = phi(ind);
        pressure     = pressure(:);
        pressure     = pressure(ind);
        theta        = theta(:);
        theta        = theta(ind);

        R            = real(RZ);
        Z	         = imag(RZ);
        BR           = real(BP);
        BZ	         = imag(BP);
        % multi matlab version
        try
        	tri          = delaunay(R,Z,{'Qt','Qbb','Qc','Qz'});
        catch
        	tri          = delaunay(R,Z);
	end
        if isempty(equilibrium.profiles_2d.grid.connect)
              equilibrium.profiles_2d.grid.connect    = 0 .* ones(size(profil0d.psi,1),size(tri,1),3);
              equilibrium.profiles_2d.grid.dim1       = NaN .* ones(size(profil0d.psi,1),length(R));
              equilibrium.profiles_2d.grid.dim2       = NaN .* ones(size(profil0d.psi,1),length(R));
              equilibrium.profiles_2d.psi 	  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.jphi 	  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.jpar 	  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.br 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.bz 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.bphi 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.vphi 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.vtheta 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.phi 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.pressure 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
              equilibrium.profiles_2d.theta 		  = NaN .* ones(size(profil0d.psi,1),length(R),1);
        end
        equilibrium.profiles_2d.grid.dim1(k,:,1)       = R; 	
        equilibrium.profiles_2d.grid.dim2(k,:,1)       = Z;
        if size(tri,1) <= size(equilibrium.profiles_2d.grid.connect,2)
            equilibrium.profiles_2d.grid.connect(k,1:size(tri,1),1:3)  = tri;
        else
            connect_mem = equilibrium.profiles_2d.grid.connect(1:max(1,k-1),:,:);
            equilibrium.profiles_2d.grid.connect    = 0 .* ones(size(profil0d.psi,1),size(tri,1),size(connect_mem,3));
            for lk =1:max(1,k-1)
                equilibrium.profiles_2d.grid.connect(lk,1:size(connect_mem,2),:)  = connect_mem(lk,:,:);
            end
            equilibrium.profiles_2d.grid.connect(k,1:size(tri,1),1:3)  = tri;
        end
        equilibrium.profiles_2d.psi(k,:,1) 	   = PSI;
        equilibrium.profiles_2d.jphi(k,:,1)	   = jphi;
        equilibrium.profiles_2d.jpar(k,:,1) 	   = jpar;
        equilibrium.profiles_2d.br(k,:,1) 		   = BR;
        equilibrium.profiles_2d.bz(k,:,1) 		   = BZ;
        equilibrium.profiles_2d.bphi(k,:,1) 	           = BPHI;
        equilibrium.profiles_2d.vtheta(k,:,1) 		   = vtheta;
        equilibrium.profiles_2d.vphi(k,:,1) 	         = vphi;
        equilibrium.profiles_2d.theta(k,:,1) 		 = theta;
        equilibrium.profiles_2d.phi(k,:,1) 	         = phi;
        equilibrium.profiles_2d.pressure(k,:,1) 	 = pressure;
    else 
        equilibrium.profiles_2d.grid.dim1(k,:)           = deuxd.PSI(:,1);
        equilibrium.profiles_2d.grid.dim2(k,:)           = deuxd.th(1,:);
        equilibrium.profiles_2d.r(k,:,:)        = deuxd.R;
        equilibrium.profiles_2d.z(k,:,:)        = deuxd.Z;
        equilibrium.profiles_2d.psi(k,:,:) 	 = deuxd.PSI;
        equilibrium.profiles_2d.jphi(k,:,:)	 = jphi;
        equilibrium.profiles_2d.jpar(k,:,:) 	 = jpar;
        equilibrium.profiles_2d.br(k,:,:) 		 = deuxd.BR;
        equilibrium.profiles_2d.bz(k,:,:) 		 = deuxd.BZ;
        equilibrium.profiles_2d.bphi(k,:,:) 	         = deuxd.BPHI;
        equilibrium.profiles_2d.vtheta(k,:,:) 		 = vtheta;
        equilibrium.profiles_2d.vphi(k,:,:) 	         = vphi;
        equilibrium.profiles_2d.theta(k,:,:) 		 = theta;
        equilibrium.profiles_2d.phi(k,:,:) 	         = phi;
        equilibrium.profiles_2d.pressure(k,:,:) 	 = pressure;
   end

    % remplissage de coord_sys
    grid2.dim1(k,:)       = deuxd.PSI(:,1);
    grid2.dim2(k,:)       = deuxd.th(1,:);
    equilibrium.coord_sys.grid_type       	   = 'polar';
    jacobian                                       = deuxd.R .* (deuxd.dRdPSI .* deuxd.dZdth - deuxd.dZdPSI .* deuxd.dRdth);
    jacobian(1,:)                                  = 0;
    equilibrium.coord_sys.jacobian(k,:,:)   	   = abs(jacobian);
    equilibrium.coord_sys.g_11(k,:,:)  	           = (deuxd.dPSIdR .^ 2 + deuxd.dPSIdZ .^ 2);
    equilibrium.coord_sys.g_12(k,:,:)   	   = deuxd.dthdR .* deuxd.dPSIdR +  deuxd.dthdZ .* deuxd.dPSIdZ;
    equilibrium.coord_sys.g_13(k,:,:)   	   = 0;
    equilibrium.coord_sys.g_22(k,:,:)   	   =  deuxd.dthdR .^2 + deuxd.dthdZ .^ 2;
    equilibrium.coord_sys.g_23(k,:,:)   	   =  0;
    equilibrium.coord_sys.g_33(k,:,:)   	   =  1./ deuxd.R .^2;
    equilibrium.coord_sys.position.r(k,:,:)        = deuxd.R;
    equilibrium.coord_sys.position.z(k,:,:)        = deuxd.Z;


    % equilibre etendu hors du plasma
    % EVALUATION IN THE VACCUUM
    method_extrap = NaN;
    if equi_extrap == 0
      equivide = zfitvide_interp(equicronos,0,[]);
      method_extrap = 0;
    else
      method_extrap = 1;
      equivide = zfitvide_jorek(equicronos,0,[]); 
      if (equivide.dpsi_error > 0.1) || (equivide.dB_error > 0.2)   
	  fprintf('changing to direct methode |');
	  method_extrap = 0;
	  dpsi_error_mem = equivide.dpsi_error;
	  equivide = zfitvide_interp(equicronos,0,[]);
	  equivide.dpsi_error = equivide.dpsi_error + dpsi_error_mem;
      end 
    end

    % (R,Z) GRID
    a   = (max(equivide.R(:)) - min(equivide.R(:))) ./ 2;
    if isfield(z0dstruct.z0dinput.option,'fixed_grid') && (z0dstruct.z0dinput.option.fixed_grid == 1)
	grid3.dim1(k,:) = linspace(box.rmin,box.rmax,nbeqdsk);
	grid3.dim2(k,:) = linspace(box.zmin,box.zmax,ceil(nbeqdsk .* rzp)) - moment.zaxe(1); 
    else
	grid3.dim1(k,:) = linspace(min(equivide.R(:)) - a ./ 4 ,max(equivide.R(:)) + a ./ 4,nbeqdsk);
	grid3.dim2(k,:) = linspace(min(equivide.Z(:)) - ( a ./ 4 .* rzp),max(equivide.Z(:)) + ( a ./ 4 .* rzp),ceil(nbeqdsk .* rzp));
    end
    [r2d,z2d] = meshgrid(grid3.dim1(k,:),grid3.dim2(k,:));

    % multi matlab version
    try
   		tri          = delaunay(r2d',z2d',{'Qt','Qbb','Qc','Qz'});
    catch
   		tri          = delaunay(r2d',z2d');
    end
    if ~isfield(grid3,'connect')
            grid3.connect(k,1:size(tri,1),1:3)  = tri;
    elseif size(tri,1) <= size(grid3.connect,2)
            grid3.connect(k,1:size(tri,1),1:3)  = tri;
    else
            connect_mem = grid3.connect(1:max(1,k-1),:,:);
            grid3.connect    = 0 .* ones(size(profil0d.psi,1),size(tri,1),size(connect_mem,3));
            for lk =1:max(1,k-1)
                grid3.connect(lk,1:size(connect_mem,2),:)  = connect_mem(lk,:,:);
            end
            grid3.connect(k,1:size(tri,1),1:3)  = tri;
    end
    if method_extrap == 0
      [psi2d,BR2d,BZ2d] = zpsivide_interp(equivide,r2d,z2d); 
    else 
      [psi2d,BR2d,BZ2d] = zpsivide_updown_16_jorek(equivide,r2d,z2d); 
    end
%      if equi_extrap == 0
%  	[psi2d BR2d,BZ2d] = zpsivide_updown_16(equivide,r2d,z2d); 
%  	psi2d = zsmooth_psi(psi2d);
%      else
%   	[psi2d BR2d,BZ2d] = zpsivide_updown_16_jorek(equivide,r2d,z2d);    
%      end
    BR2d  = BR2d .* sign(factor_two_pi);
    BZ2d  = BZ2d .* sign(factor_two_pi);
    psi2d = psi2d .* factor_two_pi;
    mask = zinout(equivide.R(end,:),equivide.Z(end,:),r2d,z2d);
    BPHI2d   = profil0d.fdia(k,end) ./ r2d;
    warning off 
    if all(diff(equilibrium.profiles_1d.psi(k,:)) > 0) || all(diff(equilibrium.profiles_1d.psi(k,:)) < 0)
      BPHI2d(mask)  = interp1(equilibrium.profiles_1d.psi(k,:),equilibrium.profiles_1d.F_dia(k,:),psi2d(mask),'pchip','extrap') ./ r2d(mask);
    else
      warning('Problem for B_field_tor interolation in 2D');
      BPHI2d(mask)  = profil0d.fdia(k,1)./ r2d(mask);
    end
    %BPHI2d(mask)  = griddata(equivide.R,equivide.Z,equivide.BZ,r2d(mask),z2d(mask),'cubic');
    %bpol_test = griddata(r2d,z2d,abs(BR2d + sqrt(-1) .* BZ2d),equivide.R,equivide.Z);
    %bpol_test = mean(bpol_test,2);
    %bpol_test_alt = mean(abs(deuxd.BR + sqrt(-1) .* deuxd.BZ),2);
    %figure(32);plot(profil0d.xli,profil0d.bpol(k,:),'r',profil0d.xli,bpol_test,'b',profil0d.xli,bpol_test_alt,'k');drawnow

    if all(diff(equilibrium.profiles_1d.psi(k,:)) > 0) || all(diff(equilibrium.profiles_1d.psi(k,:)) < 0)
	pprim  = interp1(equilibrium.profiles_1d.psi(k,:),equilibrium.profiles_1d.pprime(k,:), ...
			psi2d,'pchip','extrap');
	pprim(~mask) = 0;            
	ffprim  = interp1(equilibrium.profiles_1d.psi(k,:),equilibrium.profiles_1d.ffprime(k,:), ...
			psi2d,'pchip','extrap');
	ffprim(~mask) = 0;            
	jphi  = pprim .* r2d + ffprim ./ r2d ./ phys.mu0; 
	jphi(~mask) = 0;
	jpar  = jphi .* BPHI2d + ffprim ./ phys.mu0 ./ (BPHI2d .* r2d) .*  ...
		(BR2d .^ 2 + BZ2d .^ 2);
	jpar(~mask) = 0;	
	b0    = equilibrium.profiles_1d.F_dia(k,end) ./ profil0d.Raxe(k,end);
	jpar  = jpar ./ b0;

	psin    =  abs(profil0d.psi(k,:) - profil0d.psi(k,end)) ./ abs(profil0d.psi(k,1) - profil0d.psi(k,end));
	psin_loc = abs(equilibrium.profiles_1d.psi(k,:) - equilibrium.profiles_1d.psi(k,end)) ./ ...
		  abs(equilibrium.profiles_1d.psi(k,:) - equilibrium.profiles_1d.psi(k,end));
	omega    = interp1(psin,profil0d.omega(k,:),psin_loc,'pchip','extrap');
	vphi    = interp1(equilibrium.profiles_1d.psi(k,:),omega,psi2d,'pchip','extrap') .* r2d;
	vphi(~mask) = 0;
	utheta    = interp1(psin,profil0d.utheta(k,:),psin_loc,'pchip','extrap');
	vtheta    = interp1(equilibrium.profiles_1d.psi(k,:),utheta,psi2d,'pchip','extrap') .* abs(BR2d + sqrt(-1) .* BZ2d);
	vtheta(~mask) = 0;

	phi  = interp1(equilibrium.profiles_1d.psi(k,:),equilibrium.profiles_1d.phi(k,:), ...
			psi2d,'pchip','extrap');  
			
	phi(~mask) = 0;
	RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
	ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
	theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
	theta  = griddata(deuxd.R,deuxd.Z,theta,r2d,z2d,'cubic');
	theta(~mask) = 0;   
  
    else
	pprim  = equilibrium.profiles_1d.pprime(k,:)'  * ones(1,size(equivide.R,2));
	pprim  = griddata(equivide.R,equivide.Z,pprim,r2d,z2d,'cubic');
	ffprim = equilibrium.profiles_1d.ffprime(k,:)'  * ones(1,size(equivide.R,2));
	ffprim  = griddata(equivide.R,equivide.Z,ffprim,r2d,z2d,'cubic');
	jphi  = pprim .* r2d + ffprim ./ r2d ./ phys.mu0; 
	jphi(~mask) = 0;
	jpar  = jphi .* BPHI2d + ffprim ./ phys.mu0 ./ (BPHI2d .* r2d) .*  ...
		(BR2d .^ 2 + BZ2d .^ 2);
	jpar(~mask) = 0;
	vphi   = (profil0d.omega(k,:)' * ones(1,size(deuxd.R,2))) .* deuxd.R;
	bpol2d = abs(deuxd.BR + sqrt(-1) .*  deuxd.BZ);
	vtheta =  (profil0d.utheta(k,:)' * ones(1,size(bpol2d,2))) .* bpol2d;
	b0    = equilibrium.profiles_1d.F_dia(k,end) ./ profil0d.Raxe(k,end);
	jpar  = jpar ./ b0;
	vphi    = griddata(deuxd.R,deuxd.Z,vphi,r2d,z2d,'cubic');
	vphi(~mask) = 0;
	vtheta  = griddata(deuxd.R,deuxd.Z,vtheta,r2d,z2d,'cubic');
	vtheta(~mask) = 0;
	pressure = equilibrium.profiles_1d.pressure(k,:)' * ones(1,size(deuxd.R,2)); 
	pressure      = griddata(deuxd.R,deuxd.Z,pressure,r2d,z2d,'cubic');
	pressure(~mask) = 0;
	phi      = equilibrium.profiles_1d.phi(k,:)' * ones(1,size(deuxd.R,2)); 
	phi      = griddata(deuxd.R,deuxd.Z,phi,r2d,z2d,'cubic');
	phi(~mask) = 0;
	RC    = (min(deuxd.R,[],2) + max(deuxd.R,[],2)) ./ 2;
	ZC    = (min(deuxd.Z,[],2) + max(deuxd.Z,[],2)) ./ 2;
	theta = unwrap(angle((deuxd.R - RC * (ones(1,size(deuxd.R,2)))) + sqrt(-1) .*  (deuxd.Z - ZC * (ones(1,size(deuxd.Z,2))))));
	theta  = griddata(deuxd.R,deuxd.Z,theta,r2d,z2d,'cubic');
	theta(~mask) = 0;
    end
    warning on
    profiles_2d_2.psi(k,:,:) 	         = psi2d';
    profiles_2d_2.jphi(k,:,:)            = jphi';
    profiles_2d_2.jpar(k,:,:)            = jpar';
    profiles_2d_2.br(k,:,:)	         = BR2d';
    profiles_2d_2.bz(k,:,:) 	         = BZ2d';
    profiles_2d_2.bphi(k,:,:) 	         = BPHI2d';
    profiles_2d_2.r(k,:,:)               = r2d';
    profiles_2d_2.z(k,:,:)               = z2d' + moment.zaxe(1);
    profiles_2d_2.vtheta(k,:,:) 	 = vtheta';
    profiles_2d_2.vphi(k,:,:) 	         = vphi';
    profiles_2d_2.theta(k,:,:) 		 = theta';
    profiles_2d_2.phi(k,:,:) 	         = phi';
    profiles_2d_2.pressure(k,:,:) 	 = pressure';

    disp(' ')
    
    
%      %% graphe for control
%      profiles_2d_1 = equilibrium.profiles_2d;
%      figure(31)
%      clf
%      plot(squeeze(profiles_2d_1.r(k,:,:))',squeeze(profiles_2d_1.z(k,:,:))','r');
%      hold on
%      contour(squeeze(profiles_2d_2.r(k,:,:)),squeeze(profiles_2d_2.z(k,:,:)),squeeze(profiles_2d_2.psi(k,:,:)),equilibrium.profiles_1d.psi(k,:),'color','b','linestyle','--')
%      contour(squeeze(profiles_2d_2.r(k,:,:)),squeeze(profiles_2d_2.z(k,:,:)),squeeze(profiles_2d_2.psi(k,:,:)),101,'color','g','linestyle','-')
%      quiver(squeeze(profiles_2d_2.r(k,:,:)),squeeze(profiles_2d_2.z(k,:,:)),squeeze(profiles_2d_2.br(k,:,:)),squeeze(profiles_2d_2.bz(k,:,:)),'color','c');
%      quiver(squeeze(profiles_2d_1.r(k,:,:)),squeeze(profiles_2d_1.z(k,:,:)),squeeze(profiles_2d_1.br(k,:,:)),squeeze(profiles_2d_1.bz(k,:,:)),'color','m');
%      if isfield(profil0d,'Rsepa') &&isfield(profil0d,'Zsepa')
%        plot(profil0d.Rsepa(k,:),profil0d.Zsepa(k,:) + z0_offset ,'k');
%      end	
%      %wall = load('jt60sawall');
%      %rwall = wall.wall.data(1:end,1);
%      %zwall = wall.wall.data(1:end,2);
%      %plot(rwall,zwall,'k');
%        
%      xlabel('R (m)');
%      ylabel('Z (m)');
%      drawnow
%      %keyboard
%      %pause(1);

end
equilibrium.coord_sys.grid = grid2;
profiles_2d_2.grid = grid3;
fprintf('\n');

% convention ordre ITM
profiles_2d_1 = equilibrium.profiles_2d;
equilibrium   = rmfield(equilibrium,'profiles_2d');     
equilibrium.profiles_2d(1) = profiles_2d_2;
equilibrium.profiles_2d(2) = profiles_2d_1;




% COCOS
s_ip  = sign(mean(sign(equilibrium.global_param.i_plasma)));
s_b0  = sign(mean(sign(equilibrium.global_param.toroid_field.b0)));
s_psi = sign(mean(sign(equilibrium.profiles_1d.psi(:,end) - equilibrium.profiles_1d.psi(:,1))));
s_bp  = s_psi .* s_ip;
s_rho_phi_theta = sign(mean(sign(equilibrium.profiles_1d.q(:)))) .* s_ip .* s_b0;
if sign(mean(sign(equilibrium.profiles_1d.F_dia(:)))) ~= s_b0
    warning('mapequilibrium : sign missmatch between B0 and Fdia');
end
if sign(mean(sign(equilibrium.profiles_1d.phi(:)))) ~= s_b0
    warning('mapequilibrium : sign missmatch between B0 and Phi');
end

if  sign(mean(sign(equilibrium.profiles_1d.jphi(:)))) ~= s_ip
    warning('mapequilibrium : sign missmatch between Ip and jmoy');
end
if  sign(mean(sign(equilibrium.profiles_1d.jparallel(:)))) ~= s_ip
    warning('mapequilibrium : sign missmatch between Ip and jeff');
end


% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_itm(t,y,tt,varargin)


if length(t) <= 1
	rep = y(end,:);
else
	rep = interp1(t,y,tt,varargin{:});
end

% surcouche pour traiter le cas a 1 temps en meme temps
function rep = griddata_itm(t,x,y,tt,xx,methode)


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

%  % function to test alternative to griddata
%  function out = griddata_alt(xin,yin,var,xout,yout,methode)
%  ss = size(var);
%  
%  F = scatteredInterpolant(xin(:),yin(:),var(:),'natural','linear');
%  out = F(xout,yout);
%  out2 = griddata(xin,yin,var,xout,yout,methode);
%  delta = out(:) - out2(:);
%  sqrt(sum((delta(isfinite(delta)) .^ 2))./ sum(out2(isfinite(delta)) .^ 2))

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
  equivide.F_PSI = scatteredInterpolant(equivide.Rext,equivide.Zext,equivide.PSIext,'natural','linear');
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



