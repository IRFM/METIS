% cette fonction extrait les informations du CPO pulse_schedule pour en faire une structure z0dinput de metis
function z0dinput = ids2metis_input(pulse_schedule)

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
pulse_schedule = make_reference_data(pulse_schedule);

% parametre et info des parametres 
info = metis4imas(1);
z0dinput.info          = info.info;
z0dinput.option        = pulse_schedule2option_imas(pulse_schedule);


% langue
z0dinput.langue        =  'anglais';
% variable de sorties
z0dinput.zsinfo        = zero1t;
z0dinput.profinfo      = z0dprofinfo;
z0dinput.mode_exp      = 118;


% generate pulse schedule with homogeneous time
if isempty(pulse_schedule.time)
    pulse_schedule = make_time_pulse_schedule(pulse_schedule);    
end

% recherche des infos
z0dinput.cons.temps    = pulse_schedule.time;
% choix de la composition
switch z0dinput.option.gaz
    case 3
        if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
            z0dinput.cons.iso      = pulse_schedule.density_control.n_t_over_n_d.reference.data;
        else
            z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
        end
    case 5
        if isfield(pulse_schedule.density_control,'n_he3_over_n_d')
            if ~isempty(pulse_schedule.density_control.n_he3_over_n_d.reference.data)
                z0dinput.cons.iso      = pulse_schedule.density_control.n_he3_over_n_d.reference.data;
            else
                z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
            end
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                z0dinput.cons.iso      = z0dinput.cons.iso + sqrt(-1) .* pulse_schedule.density_control.n_t_over_n_d.reference.data;
            end
        else
            disp('data "pulse_schedule.density_control.n_he3_over_n_d" is not defined in this dataset,')
            disp('using instead "pulse_schedule.density_control.n_t_over_n_d"');
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                z0dinput.cons.iso      = pulse_schedule.density_control.n_t_over_n_d.reference.data;
            else
                z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
            end
        end
    case 11
        if isfield(pulse_schedule.density_control,'n_b_over_n_h')
            if ~isempty(pulse_schedule.density_control.n_b_over_n_h.reference.data)
                z0dinput.cons.iso      = pulse_schedule.density_control.n_b_over_n_h.reference.data;
            else
                z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
            end
        else
            disp('data "pulse_schedule.density_control.n_b_over_n_h" is not defined in this dataset,')
            disp('using instead "pulse_schedule.density_control.n_t_over_n_d"');
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                z0dinput.cons.iso      = pulse_schedule.density_control.n_t_over_n_d.reference.data;
            else
                z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
            end
        end
    otherwise
        z0dinput.cons.iso     = zeros(size(z0dinput.cons.temps));
end
% cette donnees doit etre non vide
z0dinput.cons.ip       = max(eps,abs(pulse_schedule.flux_control.i_plasma.reference.data)); % any COCOS
if ~isempty(pulse_schedule.flux_control.loop_voltage.reference.data)
    dt = diff(pulse_schedule.time);
	flux_rebuilt = -cumsum(sign(mean(pulse_schedule.flux_control.loop_voltage.reference.data)) .* pulse_schedule.flux_control.loop_voltage.reference.data(2:end) .* dt) ./ 2 ./ pi; % any COCOS ?
	flux_rebuilt = cat(1,0,flux_rebuilt);
	z0dinput.cons.flux     = flux_rebuilt;
else
	z0dinput.cons.flux     = zeros(size(z0dinput.cons.temps));
end
% cette donnees doit etre non vide
% security
if isfield(pulse_schedule.density_control.n_e_line,'n_e_line_method') && ...
   ~isempty(pulse_schedule.density_control.n_e_line.n_e_line_method.index) && ...
   (pulse_schedule.density_control.n_e_line.n_e_line_method.index ~= 3)
      error('Only supported method for pulse.density_control.n_e_line is 3: integral of a 1D core profile over rho_tor_norm up to the LCFS');
end
z0dinput.cons.nbar     = pulse_schedule.density_control.n_e_line.reference.data;
% gaspuff -> we must implement a loop over valve{k}
gaspuff = zeros(size(z0dinput.cons.temps));
if length(pulse_schedule.density_control.valve) >=1
  for k = 1:length(pulse_schedule.density_control.valve)
      if ~isempty(pulse_schedule.density_control.valve{k}.flow_rate.reference.data)
	    %loop over elements
	    zsum = 0;
	    if isfield(pulse_schedule.density_control.valve{k},'species') && (length(pulse_schedule.density_control.valve{k}.species) >=1)
	      for l=1:length(pulse_schedule.density_control.valve{k}.species)
		  if length(pulse_schedule.density_control.valve{k}.species{l}.element) >=1
		    for m=1:length(pulse_schedule.density_control.valve{k}.species{l}.element);
		      zsum = zsum + double(pulse_schedule.density_control.valve{k}.species{l}.element{m}.z_n) .* ...
		                    double(pulse_schedule.density_control.valve{k}.species{l}.fraction);
		    end
		  elseif isfinite(pulse_schedule.density_control.valve{k}.species{l}.fraction) && (pulse_schedule.density_control.valve{k}.species{l}.fraction > 0)
		      zsum = zsum + 1 .* double(pulse_schedule.density_control.valve{k}.species{l}.fraction);
		  else
		      zsum = zsum + 1		      
		  end
		  gaspuff = gaspuff + zsum .* pulse_schedule.density_control.valve{k}.flow_rate.reference.data;	  
	      end
	    else
		% assume H2,D2,T2 or He4
		gaspuff = gaspuff + pulse_schedule.density_control.valve{k}.flow_rate.reference.data .* 2;
	    end
      end
  end
end
z0dinput.cons.nbar     = max(1,z0dinput.cons.nbar) + sqrt(-1) .* gaspuff;

%% IC
if isfield(pulse_schedule.ic,'antenna')
    if (length(pulse_schedule.ic.antenna) >= 1) && ~isempty(pulse_schedule.ic.antenna{1}.power.reference.data)
	  z0dinput.cons.picrh    = pulse_schedule.ic.antenna{1}.power.reference.data;
    else
	  z0dinput.cons.picrh    = zeros(size(z0dinput.cons.temps));
    end
else
    if (length(pulse_schedule.ic.launcher) >= 1) && ~isempty(pulse_schedule.ic.launcher{1}.power.reference.data)
	  z0dinput.cons.picrh    = pulse_schedule.ic.launcher{1}.power.reference.data;	
    else
	  z0dinput.cons.picrh    = zeros(size(z0dinput.cons.temps));
    end
end

%% LH
if isfield(pulse_schedule.lh,'antenna')
    if (length(pulse_schedule.lh.antenna)>=1) && ~isempty(pulse_schedule.lh.antenna{1}.power.reference.data)
	    z0dinput.cons.plh      = pulse_schedule.lh.antenna{1}.power.reference.data; 
    else
	    z0dinput.cons.plh      = zeros(size(z0dinput.cons.temps));
    end
else
    if (length(pulse_schedule.lh.launcher)>=1) && ~isempty(pulse_schedule.lh.launcher{1}.power.reference.data)
	    z0dinput.cons.plh      = pulse_schedule.lh.launcher{1}.power.reference.data; 	
    else
	    z0dinput.cons.plh      = zeros(size(z0dinput.cons.temps));
    end

end

%% NBI
pnbi = zeros(length(z0dinput.cons.temps),8);
% loop on NBI unit
if length(pulse_schedule.nbi.unit) >= 1
  for k = 1:length(pulse_schedule.nbi.unit)
      if ~isempty(pulse_schedule.nbi.unit{k}.power.reference.data) && (k<=8)
          pnbi(:,k) = pnbi(:,k) + pulse_schedule.nbi.unit{k}.power.reference.data;
      end
  end
end
switch z0dinput.option.gaz
    case 3
        z0dinput.cons.pnbi     =  (pnbi(:,2) + pnbi(:,3)) + sqrt(-1) .* (pnbi(:,5) + pnbi(:,6));
        z0dinput.cons.ftnbi    =  pnbi(:,3) ./ max(1,real(z0dinput.cons.pnbi)) + sqrt(-1) .* pnbi(:,6) ./ max(1,imag(z0dinput.cons.pnbi));
    case 5
        % No he3 in NBI only D
        z0dinput.cons.pnbi     =  pnbi(:,2) + sqrt(-1) .* pnbi(:,5);
        z0dinput.cons.ftnbi    =  zeros(size(z0dinput.cons.pnbi));
    case 11
        % boron on 7 and 8
        z0dinput.cons.pnbi     =  (pnbi(:,1) + pnbi(:,7)) + sqrt(-1) .* (pnbi(:,4) + pnbi(:,8));
        z0dinput.cons.ftnbi    =  pnbi(:,7) ./ max(1,real(z0dinput.cons.pnbi)) + sqrt(-1) .* pnbi(:,8) ./ max(1,imag(z0dinput.cons.pnbi));
    otherwise
        z0dinput.cons.pnbi     =  (pnbi(:,2) + pnbi(:,1)) + sqrt(-1) .* (pnbi(:,5) + pnbi(:,4));
        z0dinput.cons.ftnbi    =  pnbi(:,1) ./ max(1,real(z0dinput.cons.pnbi)); + sqrt(-1) .* pnbi(:,4) ./ max(1,imag(z0dinput.cons.pnbi));
end

%%EC
if isfield(pulse_schedule.ec,'beam')
 	if (length(pulse_schedule.ec.beam)>=1) && ~isempty(pulse_schedule.ec.beam{1}.power_launched.reference.data)
		z0dinput.cons.pecrh     = pulse_schedule.ec.beam{1}.power_launched.reference.data;
	else
		z0dinput.cons.pecrh     = zeros(size(z0dinput.cons.temps));
	end
   
elseif isfield(pulse_schedule.ec,'antenna')
	if (length(pulse_schedule.ec.antenna)>=1) && ~isempty(pulse_schedule.ec.antenna{1}.power.reference.data)
		z0dinput.cons.pecrh     = pulse_schedule.ec.antenna{1}.power.reference.data;
	else
		z0dinput.cons.pecrh     = zeros(size(z0dinput.cons.temps));
	end
else
	if (length(pulse_schedule.ec.launcher)>=1) && ~isempty(pulse_schedule.ec.launcher{1}.power.reference.data)	
		z0dinput.cons.pecrh     = pulse_schedule.ec.launcher{1}.power.reference.data;	
	else
		z0dinput.cons.pecrh     = zeros(size(z0dinput.cons.temps));
	end
end
if isfield(pulse_schedule.ec,'beam')
      if ~isempty(pulse_schedule.ec.beam) && (length(pulse_schedule.ec.beam) > 1) && ...
	~isempty(pulse_schedule.ec.beam{2}.power_launched.reference.data)
	      z0dinput.cons.plh     = pulse_schedule.ec.beam{2}.power_launched.reference.data;
	      z0dinput.option.lhmode  = 5;
      end
elseif isfield(pulse_schedule.ec,'antenna')
      if ~isempty(pulse_schedule.ec.antenna) && (length(pulse_schedule.ec.antenna) > 1) && ...
	~isempty(pulse_schedule.ec.antenna{2}.power.reference.data)
	      z0dinput.cons.plh     = pulse_schedule.ec.antenna{2}.power.reference.data;
	      z0dinput.option.lhmode  = 5;
      end
else
      if ~isempty(pulse_schedule.ec.launcher) && (length(pulse_schedule.ec.launcher) > 1) && ...
	~isempty(pulse_schedule.ec.launcher{2}.power.reference.data)
	      z0dinput.cons.plh     = pulse_schedule.ec.launcher{2}.power.reference.data;	
	      z0dinput.option.lhmode  = 5;
      end
end
if isfield(pulse_schedule.ec,'beam')
    if (length(pulse_schedule.ec.beam)>=1) && ~isempty(pulse_schedule.ec.beam{1}.deposition_rho_tor_norm.reference.data)
	   z0dinput.cons.xece      = pulse_schedule.ec.beam{1}.deposition_rho_tor_norm.reference.data;
    else
	   z0dinput.cons.xece     = zeros(size(z0dinput.cons.temps));
    end
elseif isfield(pulse_schedule.ec,'antenna')
    if (length(pulse_schedule.ec.antenna)>=1) && ~isempty(pulse_schedule.ec.antenna{1}.deposition_rho_tor_norm.reference.data)
	   z0dinput.cons.xece      = pulse_schedule.ec.antenna{1}.deposition_rho_tor_norm.reference.data;
    else
	    z0dinput.cons.xece     = zeros(size(z0dinput.cons.temps));
    end
else
    if (length(pulse_schedule.ec.launcher)>=1) && ~isempty(pulse_schedule.ec.launcher{1}.deposition_rho_tor_norm.reference.data)
	    z0dinput.cons.xece      = pulse_schedule.ec.launcher{1}.deposition_rho_tor_norm.reference.data;	

    else
	    z0dinput.cons.xece     = zeros(size(z0dinput.cons.temps));
    end

end
%% enhancement factor
if ~isempty(pulse_schedule.flux_control.li_3.reference.data) && ~isempty(pulse_schedule.flux_control.beta_normal.reference.data)
    z0dinput.cons.hmore = pulse_schedule.flux_control.beta_normal.reference.data ./ pulse_schedule.flux_control.li_3.reference.data ./ 4;
else
    z0dinput.cons.hmore     = ones(size(z0dinput.cons.temps));
end

% cette donnees doit etre non vide
% security
if isfield(pulse_schedule.density_control.zeff,'zeff_method') && ...
   ~isempty(pulse_schedule.density_control.zeff.zeff_method.index) && ...
   (pulse_schedule.density_control.zeff.zeff_method.index ~= 3)
      error('Only supported method for pulse.density_control.zeff is 3: average of a 1D core profile over rho_tor_norm up to the LCFS');
end
z0dinput.cons.zeff     = max(1,pulse_schedule.density_control.zeff.reference.data);
z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = 3;
% Geometry
% selon les donnees disponibles

if ~isempty(pulse_schedule.position_control.boundary_outline{1}.r.reference.data) &&  ~isempty(pulse_schedule.position_control.boundary_outline{1}.z.reference.data)
	% cas separatrice donnees par points
	% rebuilt of LCFS from references
	nbpts  = length(pulse_schedule.position_control.boundary_outline);
	sepa.R = NaN .* ones(length(z0dinput.cons.temps),nbpts);
	sepa.Z = sepa.R;
	for k=1:nbpts
	      sepa.R(:,k) = pulse_schedule.position_control.boundary_outline{k}.r.reference.data;
	      sepa.Z(:,k) = pulse_schedule.position_control.boundary_outline{k}.z.reference.data;	      
	end
	% calcul des moments
	% la courbe doit etre fermee
	if (sepa.R(1,1) ~= sepa.R(1,end)) || (sepa.Z(1,1) ~= sepa.Z(1,end))
		sepa.R(:,end+1) = sepa.R(:,1);
		sepa.Z(:,end+1) = sepa.Z(:,1);
	end

	% calcul des moments
	% centre pour angle d'integration
	rc = mean(sepa.R,2);
	zc = mean(sepa.Z,2);
	vc = ones(1,size(sepa.R,2));
	uc = unwrap(angle((sepa.R-rc*vc) + sqrt(-1) .* (sepa.Z  -zc*vc)));
	uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
	uc(:,1)   = uc(:,end) + 2 .* pi;
	xu    = linspace(0,1,length(vc));
	%dudx  = pdederive(xu,uc,2,2,2,1);
	%dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
	%dudx(:,end) = dudx(:,1);
	dRdx  = pdederive(xu,sepa.R,2,2,2,1);
	%dZdx  = pdederive(xu,sepa.Z,2,2,2,1);
	% calcul de R0 et Z0
	maskrmax  = (sepa.R == (max(sepa.R,[],2) * vc));
	geo.z0        = sum(sepa.Z .* maskrmax,2) ./ sum(maskrmax,2);
	% recalcul des parametres sur le vecteur final
	rmin  = min(sepa.R,[],2);
	rmax  = max(sepa.R,[],2);
	geo.a = 0.5 .* (rmax - rmin);
	geo.R = 0.5 .* (rmax + rmin);
	zmin  = min(sepa.Z,[],2);
	zmax  = max(sepa.Z,[],2);
	%geo.z0      =(zmax + zmin) ./ 2;
	geo.K    = abs(trapz(xu,sepa.Z .*  dRdx,2) ./ pi ./ geo.a .^ 2);
	%geo.K  = (zmax -zmin) ./ 2 ./ geo.a;
	rzmax = geo.R;
	rzmin = geo.R;
	for k = 1:size(sepa.Z,1)
		rzmax(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmax(k))));
		rzmin(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmin(k))));
	end
	uu   =  angle(rzmax - geo.R + sqrt(-1) .* (zmax - geo.z0));
	ul   =  angle(rzmin - geo.R + sqrt(-1) .* (zmin - geo.z0));
	tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
	tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
	tm   =  (tl + tu) ./ 2;
	geo.d = sin(tm);


	z0dinput.geo.R       = max(0.1,geo.R);      % grand rayon du plasma (m)
	z0dinput.geo.z0      = geo.z0;     % centre geometricque du plasma en Z (m)
	z0dinput.geo.a       = max(eps,geo.a);      % petit rayon du plasma (m)
	z0dinput.geo.K       = max(0.1,geo.K);     % elongation (b/a)
	z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
	z0dinput.exp0d.Rsepa = max(0.1 - eps,sepa.R);       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)
else
	% cas separatrice donnees par moments	
	if isfield(z0dinput,'exp0d') && isfield(z0dinput.exp0d,'Rsepa')&& isfield(z0dinput.exp0d,'Zsepa')
	    z0dinput.exp0d = rmfield(z0dinput.exp0d,'Rsepa');
	    z0dinput.exp0d = rmfield(z0dinput.exp0d,'Zsepa');
	end
   	z0dinput.geo.a     	= pulse_schedule.position_control.minor_radius.reference.data;
   	if isfield(pulse_schedule.position_control,'geometric_axis')
	      z0dinput.geo.R     	= max(0.1,pulse_schedule.position_control.geometric_axis.r.reference.data);
	      z0dinput.geo.z0    	= pulse_schedule.position_control.geometric_axis.z.reference.data;
   	else
	      disp('test mode for pulse_schedule: you are connected to a obsolete version of IMAS');
	      z0dinput.geo.R     	= max(0.1,pulse_schedule.position_control.magnetic_axis.r.reference.data);
	      z0dinput.geo.z0    	= pulse_schedule.position_control.magnetic_axis.z.reference.data;
	end  
	if isempty(pulse_schedule.position_control.elongation)
	      z0dinput.geo.K     	= max(0.1,(pulse_schedule.position_control.elongation_upper.reference.data + ...
					  pulse_schedule.position_control.elongation_lower.reference.data) ./ 2);
	else
	      z0dinput.geo.K     	= max(0.1,pulse_schedule.position_control.elongation.reference.data);
	end
	if isempty(pulse_schedule.position_control.triangularity)
	      z0dinput.geo.d     	= (pulse_schedule.position_control.triangularity_upper.reference.data + ...
					  pulse_schedule.position_control.triangularity_lower.reference.data) ./ 2;
	else
	      z0dinput.geo.d     	= pulse_schedule.position_control.triangularity.reference.data;
	end
	%[vps,sps,sexts]         = zgeo0(z0dinput.geo);
	%z0dinput.geo.vp    	= vps;
  	%z0dinput.geo.sp    	= sps;
   	%z0dinput.geo.sext    	= sexts;
end

% vaccum magnetic field (temporary waiting for final version of pulse schedule)
z0dinput.geo.b0 = max(50e-6,abs(pulse_schedule.tf.b_field_tor_vacuum_r.reference.data) ./ z0dinput.geo.R);

% informations generales   
z0dinput.machine   = z0dinput.option.machine;
z0dinput.shot      = z0dinput.option.shot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% securite 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mise a jour de la structure experimentale vide
noms = fieldnames(z0dinput.zsinfo);
if ~isfield(z0dinput,'exp0d')
    z0dinput.exp0d=[];
end
exp0d  = z0dinput.exp0d;
if isfield(exp0d,'temps')
	texp = exp0d.temps;
else
	texp = z0dinput.cons.temps;
	exp0d.temps = texp;
end
nbt  = length(texp);
%
vtnan = NaN .* ones(nbt,1);
for k = 1:length(noms)
	nomc = noms{k};
	if isfield(exp0d,nomc)
		var = getfield(exp0d,nomc);
		if length(var) ~= nbt
			disp('dimension mismatch')
			var = mean(var(isfinite(var))) .* ones(nbt,1);
	  		exp0d = setfield(exp0d,nomc,var);
		else
			% si donnnees non valides
			fnan = imag(var);
			var  = real(var);
			ind  = find(fnan~=0 & var == 0);
			if  ~isempty(ind)
				var(ind) = NaN;
			end 
	  		exp0d  = setfield(exp0d,nomc,var);
		end
	else
	  	exp0d = setfield(exp0d,nomc,vtnan);
	end
end

% donnees experimentale
z0dinput.exp0d = exp0d;


if ~isfield(z0dinput.cons,'xece')
      z0dinput.cons.xece = zeros(size(z0dinput.cons.temps));
end

% mise a jour de cons.iso
if ~isfield(z0dinput.cons,'iso')
   z0dinput.cons.iso = zeros(size(z0dinput.cons.temps)); 
elseif length(z0dinput.cons.iso) == 1
   z0dinput.cons.iso = z0dinput.cons.iso .* ones(size(z0dinput.cons.temps)); 
end
% consigne d'injection de tritium par nbi (fraction de la puissance)
if ~isfield(z0dinput.cons,'ftnbi')
   z0dinput.cons.ftnbi = min(1,z0dinput.cons.iso .* 0.5);
elseif length(z0dinput.cons.iso) == 1
   z0dinput.cons.ftnbi = z0dinput.cons.ftnbi .* ones(size(z0dinput.cons.temps)); 
end

% securite mise en forme et NaN
noms = fieldnames(z0dinput.cons);
for k=1:length(noms)
   nomc = noms{k};
   val = z0dinput.cons.(nomc);
   val(~isfinite(val)) = 0;
   z0dinput.cons = setfield(z0dinput.cons,nomc,val(:));
end
noms = fieldnames(z0dinput.geo);
for k=1:length(noms)
   nomc = noms{k};
   val = z0dinput.geo.(nomc);
   val(~isfinite(val)) = 0;
   z0dinput.geo = setfield(z0dinput.geo,nomc,val(:));
end

% securite sur le zeff
if z0dinput.option.gaz == 4
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(2.2,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
else
   z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = z0dinput.option.zmax - 0.1;
   z0dinput.cons.zeff = max(1.1,min(z0dinput.cons.zeff,z0dinput.option.zmax - 0.1));
end


% gestion auto du ripple
if strcmp(z0dinput.machine,'TS')
   z0dinput.option.rip = 1;
else
   z0dinput.option.rip = 0;

end

%securite largeur LH
if ~isfinite(z0dinput.option.dlh)
	z0dinput.option.dlh = 0.2;
	z0dinput.option.xlh = 0.2;
end



if ~isfinite(z0dinput.option.npar0)
	z0dinput.option.npar = 2;
end

% securite geo
z0dinput.geo.a = max(z0dinput.geo.a,1e-2);
z0dinput.geo.R = max(z0dinput.geo.R,3e-2);
z0dinput.geo.K = max(z0dinput.geo.K,0.1);
z0dinput.geo.b0 = max(abs(z0dinput.geo.b0),1e-4);
