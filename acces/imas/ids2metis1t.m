%
% lecture des consignes pour la version evolution
%
function [option,time,cons1t,geo1t,sepa1t] = ids2metis1t(pulse_schedule)

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

% extraction des options METIS de codeparam.
option = pulse_schedule2option_imas(pulse_schedule);
if isempty(option)
	iof    = metis4imas(1);
	option = iof.valeur;
end
% valeur du temps
if ~isempty(pulse_schedule.time)
	time = pulse_schedule.time(end);
else
	time = [];
end
cons1t.temps = time;
 
% les consignes a 1 temps
% recherche des infos
% choix de la composition
switch option.gaz
    case 3
        if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
            cons1t.iso       = pulse_schedule.density_control.n_t_over_n_d.reference.data(end);
        else
            cons1t.iso     = 0;
        end
    case 5
        if isfield(pulse_schedule.density_control,'n_he3_over_n_d')
            if ~isempty(pulse_schedule.density_control.n_he3_over_n_d.reference.data)
                cons1t.iso       = pulse_schedule.density_control.n_he3_over_n_d.reference.data(end);
            else
                cons1t.iso     = 0;
            end
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                cons1t.iso       = cons.iso + sqrt(-1) .* pulse_schedule.density_control.n_t_over_n_d.reference.data(end);
            end
        else
            disp('data "pulse_schedule.density_control.n_he3_over_n_d" is not defined in this dataset,')
            disp('using instead "pulse_schedule.density_control.n_t_over_n_d"');
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                cons1t.iso       = pulse_schedule.density_control.n_t_over_n_d.reference.data(end);
            else
                cons1t.iso     = 0;
            end           
        end
    case 11
        if isfield(pulse_schedule.density_control,'n_b_over_n_h')
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                cons1t.iso       = pulse_schedule.density_control.n_b_over_n_h.reference.data(end);
            else
                cons1t.iso     = 0;
            end
        else
            disp('data "pulse_schedule.density_control.n_b_over_n_h" is not defined in this dataset,')
            disp('using instead "pulse_schedule.density_control.n_t_over_n_d"');
            if ~isempty(pulse_schedule.density_control.n_t_over_n_d.reference.data)
                cons1t.iso       = pulse_schedule.density_control.n_t_over_n_d.reference.data(end);
            else
                cons1t.iso     = 0;
            end           
        end
    otherwise
        cons1t.iso     = 0;
end
% cette donnees doit etre non vide
cons1t.ip       = abs(pulse_schedule.flux_control.i_plasma.reference.data(end)); % any COCOS
if ~isempty(pulse_schedule.flux_control.loop_voltage.reference.data)
        % temprorary field, will be translated later in poloidal edge flux reference
	cons1t.vloop = abs(pulse_schedule.flux_control.loop_voltage.reference.data(end)); % It is COCOS dependent and orientation dependent; it is assumed to be >=0
else
	cons1t.vloop = 0;
end
% cette donnees doit etre non vide
% security
if isfield(pulse_schedule.density_control.n_e_line,'n_e_line_method') && ...
   ~isempty(pulse_schedule.density_control.n_e_line.n_e_line_method.index) && ...
   (pulse_schedule.density_control.n_e_line.n_e_line_method.index ~= 3)
      error('Only supported method for pulse.density_control.n_e_line is 3: integral of a 1D core profile over rho_tor_norm up to the LCFS');
end
cons1t.nbar     = pulse_schedule.density_control.n_e_line.reference.data(end);
% gaspuff
% gaspuff -> we must implement a loop over valve{k}
gaspuff = 0;
if length(pulse_schedule.density_control.valve) >=1
    for k = 1:length(pulse_schedule.density_control.valve)
        if ~isempty(pulse_schedule.density_control.valve{k}.flow_rate.reference.data)
            %loop over elements
            zsum = 0;
            if length(pulse_schedule.density_control.valve{k}.species) >=1
                for l=1:length(pulse_schedule.density_control.valve{k}.species)
                    if length(pulse_schedule.density_control.valve{k}.species{l}.element) >=1
                        for m=1:length(pulse_schedule.density_control.valve{k}.species{l}.element);
                            zsum = zsum + double(pulse_schedule.density_control.valve{k}.species{l}.element{m}.z_n) .* ...
                                double(pulse_schedule.density_control.valve{k}.species{l}.fraction);
                        end
                    elseif isfinite(pulse_schedule.density_control.valve{k}.species{l}.fraction) && (pulse_schedule.density_control.valve{k}.species{l}.fraction > 0)
                        zsum = zsum + 1 .* double(pulse_schedule.density_control.valve{k}.species{l}.fraction);
                    else
                        zsum = zsum + 1;
                    end
                    gaspuff = gaspuff + zsum .* pulse_schedule.density_control.valve{k}.flow_rate.reference.data(end);
                end
            else
                % assume H2,D2,T2 or He4
                gaspuff = gaspuff + pulse_schedule.density_control.valve{k}.flow_rate.reference.data(end) .* 2;
            end
        end
    end
end
cons1t.nbar = cons1t.nbar + sqrt(-1) .* gaspuff;

%% IC
if isfield(pulse_schedule.ic,'antenna')
    if ~isempty(pulse_schedule.ic.antenna) && ~isempty(pulse_schedule.ic.antenna{1}.power.reference.data)
	    cons1t.picrh    = pulse_schedule.ic.antenna{1}.power.reference.data(end);
    else
	    cons1t.picrh    = 0;
    end
else
    if ~isempty(pulse_schedule.ic.launcher) && ~isempty(pulse_schedule.ic.launcher{1}.power.reference.data)
            cons1t.picrh    = pulse_schedule.ic.launcher{1}.power.reference.data(end);	  
    else
	    cons1t.picrh    = 0;
    end
end

%% LH
if isfield(pulse_schedule.lh,'antenna')
    if ~isempty(pulse_schedule.lh.antenna) && ~isempty(pulse_schedule.lh.antenna{1}.power.reference.data)
            cons1t.plh      = pulse_schedule.lh.antenna{1}.power.reference.data(end); 
    else
	    cons1t.plh      = 0;
    end
else
    if ~isempty(pulse_schedule.lh.launcher) && ~isempty(pulse_schedule.lh.launcher{1}.power.reference.data)
	    cons1t.plh      = pulse_schedule.lh.launcher{1}.power.reference.data(end); 	
    else
	    cons1t.plh      = 0;
    end
end

%% NBI
pnbi = zeros(1,8);
% loop on NBI unit
if length(pulse_schedule.nbi.unit) >= 1
    for k = 1:length(pulse_schedule.nbi.unit)
        if ~isempty(pulse_schedule.nbi.unit{k}.power.reference.data) && (k<=8)
            pnbi(1,k) = pnbi(1,k) + pulse_schedule.nbi.unit{k}.power.reference.data(end);
        end
    end
end
switch option.gaz
    case 3
        cons1t.pnbi     =  (pnbi(1,2) + pnbi(1,3)) + sqrt(-1) .* (pnbi(1,5) + pnbi(1,6));
        cons1t.ftnbi    =  pnbi(1,3) ./ max(1,real(cons1t.pnbi)) + sqrt(-1) .* pnbi(1,6) ./  max(1,imag(cons1t.pnbi));
    case 5
        % No he3 in NBI only D
        cons1t.pnbi     =  pnbi(1,2)  + sqrt(-1) .* pnbi(1,5);
        cons1t.ftnbi    =  zeros(size(cons1t.pnbi));
    case 11   
        % boron on 7 and 8
        cons1t.pnbi     =  (pnbi(1,1) + pnbi(1,7)) + sqrt(-1) .* (pnbi(1,4) + pnbi(1,8));
        cons1t.ftnbi    =  pnbi(1,7) ./ max(1,real(cons1t.pnbi)) + sqrt(-1) .* pnbi(1,8) ./ max(1,imag(cons1t.pnbi));        
    otherwise
        cons1t.pnbi     =  (pnbi(1,2) + pnbi(1,1)) + sqrt(-1) .* (pnbi(1,5) + pnbi(1,4));
        cons1t.ftnbi    =  pnbi(1,1) ./ max(1,real(cons1t.pnbi)) + sqrt(-1) .* pnbi(1,4) ./ max(1,imag(cons1t.pnbi));
end


%% EC
if isfield(pulse_schedule.ec,'beam')
    if ~isempty(pulse_schedule.ec.beam) && ~isempty(pulse_schedule.ec.beam{1}.power_launched.reference.data)
	    cons1t.pecrh     = pulse_schedule.ec.beam{1}.power_launched.reference.data(end);
    else
	    cons1t.pecrh     = 0;
    end
elseif isfield(pulse_schedule.ec,'antenna')
    if ~isempty(pulse_schedule.ec.antenna) && ~isempty(pulse_schedule.ec.antenna{1}.power.reference.data)
	    cons1t.pecrh     = pulse_schedule.ec.antenna{1}.power.reference.data(end);
    else
	    cons1t.pecrh     = 0;
    end
else 
    if ~isempty(pulse_schedule.ec.launcher) && ~isempty(pulse_schedule.ec.launcher{1}.power.reference.data)
            cons1t.pecrh     = pulse_schedule.ec.launcher{1}.power.reference.data(end);	
    else
	    cons1t.pecrh     = 0;
    end
end
if isfield(pulse_schedule.ec,'beam')
    if ~isempty(pulse_schedule.ec.beam) && (length(pulse_schedule.ec.beam) > 1) && ...
      ~isempty(pulse_schedule.ec.beam{2}.power_launched.reference.data)
	    cons1t.plh     = pulse_schedule.ec.beam{2}.power_launched.reference.data(end);
	    option.lhmode  = 5;
    end
elseif isfield(pulse_schedule.ec,'antenna')
    if ~isempty(pulse_schedule.ec.antenna) && (length(pulse_schedule.ec.antenna) > 1) && ...
      ~isempty(pulse_schedule.ec.antenna{2}.power.reference.data)
	    cons1t.plh     = pulse_schedule.ec.antenna{2}.power.reference.data(end);
	    option.lhmode  = 5;
    end
else
    if ~isempty(pulse_schedule.ec.launcher) && (length(pulse_schedule.ec.launcher) > 1) && ...
      ~isempty(pulse_schedule.ec.launcher{2}.power.reference.data)
	    cons1t.plh     = pulse_schedule.ec.launcher{2}.power.reference.data(end);	
	    option.lhmode  = 5;
    end

end
if isfield(pulse_schedule.ec,'beam')
    if (length(pulse_schedule.ec.beam)>=1) && ~isempty(pulse_schedule.ec.beam{1}.deposition_rho_tor_norm.reference.data)
	    cons1t.xece      = pulse_schedule.ec.beam{1}.deposition_rho_tor_norm.reference.data(end);
    else
	    cons1t.xece     = 0;
    end
elseif isfield(pulse_schedule.ec,'antenna')
    if (length(pulse_schedule.ec.antenna)>=1) && ~isempty(pulse_schedule.ec.antenna{1}.deposition_rho_tor_norm.reference.data)
	    cons1t.xece      = pulse_schedule.ec.antenna{1}.deposition_rho_tor_norm.reference.data(end);
    else
	    cons1t.xece     = 0;
    end
else 
    if (length(pulse_schedule.ec.launcher)>=1) && ~isempty(pulse_schedule.ec.launcher{1}.deposition_rho_tor_norm.reference.data)
	    cons1t.xece      = pulse_schedule.ec.launcher{1}.deposition_rho_tor_norm.reference.data(end);	
    else
	    cons1t.xece     = 0;
    end
end
%% Hmore
%% enhancement factor
if ~isempty(pulse_schedule.flux_control.li_3.reference.data) && ~isempty(pulse_schedule.flux_control.beta_normal.reference.data)
    cons1t.hmore = pulse_schedule.flux_control.beta_normal.reference.data(end) ./ pulse_schedule.flux_control.li_3.reference.data(end) ./ 4;
else
    cons1t.hmore     = 1;
end

% cette donnees doit etre non vide
% security
if isfield(pulse_schedule.density_control.zeff,'zeff_method') && ...
   ~isempty(pulse_schedule.density_control.zeff.zeff_method.index) && ...
   (pulse_schedule.density_control.zeff.zeff_method.index ~= 3)
      error('Only supported method for pulse.density_control.zeff is 3: average of a 1D core profile over rho_tor_norm up to the LCFS');
end
cons1t.zeff     = pulse_schedule.density_control.zeff.reference.data(end);
cons1t.zeff(~isfinite(cons1t.zeff)) = 3;

% selon les donnees disponibles
% Geometry
% selon les donnees disponibles
if ~isempty(pulse_schedule.position_control.boundary_outline{1}.r.reference.data) && ...
   ~isempty(pulse_schedule.position_control.boundary_outline{1}.r.reference.data) 
   
	% cas separatrice donnees par points
	% rebuilt of LCFS from references
	nbpts  = length(pulse_schedule.position_control.boundary_outline);
	sepa.R = NaN .* ones(1,nbpts);
	sepa.Z = sepa.R;
	for k=1:nbpts
	      sepa.R(1,k) = pulse_schedule.position_control.boundary_outline{k}.r.reference.data(end);
	      sepa.Z(1,k) = pulse_schedule.position_control.boundary_outline{k}.z.reference.data(end);	      
	end
	% calcul des moments
	% la courbe doit etre fermee
	if (sepa.R(1,1) ~= sepa.R(1,end)) | (sepa.Z(1,1) ~= sepa.Z(1,end))
		sepa.R(1,end+1) = sepa.R(1,1);
		sepa.Z(1,end+1) = sepa.Z(1,1);
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
	rzmax = sepa.R(1,min(find(sepa.Z(1,:) == zmax)));
	rzmin = sepa.R(1,min(find(sepa.Z(1,:) == zmin)));
	uu   =  angle(rzmax - geo.R + sqrt(-1) .* (zmax - geo.z0));
	ul   =  angle(rzmin - geo.R + sqrt(-1) .* (zmin - geo.z0));
	tu   =  abs((acos((rzmax - geo.R) ./ geo.a) - acos(cos(uu))) ./ sin(uu));
	tl   =  abs((acos((rzmin - geo.R) ./ geo.a) - acos(cos(ul))) ./ sin(ul));
	tm   =  (tl + tu) ./ 2;
	geo.d = sin(tm);


	geo1t.R       = geo.R(end);      % grand rayon du plasma (m)
	geo1t.z0      = geo.z0(end);   % centre geometricque du plasma en Z (m)
	geo1t.a       = geo.a(end)  ;      % petit rayon du plasma (m)
	geo1t.K       = geo.K(end)  ;     % elongation (b/a)
	geo1t.d       = geo.d(end)  ;    % triangularite haute (definition entree de helena)
	%
	sepa1t.Rsepa = sepa.R(end,:);       % vecteur R des points de la separatrice (m)
	sepa1t.Zsepa = sepa.Z(end,:)   - geo.z0(end)   * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)
else
	% cas separatrice donnees par moments	
	sepa1t = [];
   	geo1t.a     	= pulse_schedule.position_control.minor_radius.reference.data(end);
   	if isfield(pulse_schedule.position_control,'geometric_axis')
	      geo1t.R     	= pulse_schedule.position_control.geometric_axis.r.reference.data(end);
	      geo1t.z0    	= pulse_schedule.position_control.geometric_axis.z.reference.data(end);
   	else
	      disp('test mode for pulse_schedule: you are connected to a obsolete version of IMAS');
	      geo1t.R     	= pulse_schedule.position_control.magnetic_axis.r.reference.data(end);
	      geo1t.z0    	= pulse_schedule.position_control.magnetic_axis.z.reference.data(end);
	end  
	if isempty(pulse_schedule.position_control.elongation)
	      geo1t.K     	= (pulse_schedule.position_control.elongation_upper.reference.data(end) + ...
					  pulse_schedule.position_control.elongation_lower.reference.data(end)) ./ 2;
	else
	      geo1t.K     	= pulse_schedule.position_control.elongation.reference.data(end);
	end
	if isempty(pulse_schedule.position_control.triangularity)
	      geo1t.d     	= (pulse_schedule.position_control.triangularity_upper.reference.data(end) + ...
					  pulse_schedule.position_control.triangularity_lower.reference.data(end)) ./ 2;
	else
	      geo1t.d     	= pulse_schedule.position_control.triangularity.reference.data(end);
	end
end

% vaccum magnetic field (temporary waiting for final version of pulse schedule)
geo1t.b0 = pulse_schedule.tf.b_field_tor_vacuum_r.reference.data(end) ./ ...
                  pulse_schedule.position_control.magnetic_axis.r.reference.data(end);
% converting frequency in position
% freq = phys.e .* geo.b0 .* geo.R ./ (cons.xece .* geo.a  + geo.R)./ phys.me;
% xece = (phys.e .* geo.b0 ./ phys.me ./ freq - 1) .* geo.R ./ geo.a;
% waiting for proper definition !
%  if ~isempty(pulse_schedule.ec.antenna) && ~isempty(pulse_schedule.ec.antenna{1}.frequency.reference.data)
%          freq = pulse_schedule.ec.antenna{1}.frequency.reference.data(end);
%  	xece = (phys.e .* geo1t.b0 ./ phys.me ./ freq  - 1) .* geo1t.R  ./ geo1t.a;
%          cons1t.xece     = max(0,min(1,abs(xece))); 
%          cons1t.xece(cons1t.xece < sqrt(eps)) =  0; 
%  else
%  	cons1t.xece     = 0;
%  end
% vaccum magnetic field (temporary waiting for final version of pulse schedule)
geo1t.b0 = pulse_schedule.tf.b_field_tor_vacuum_r.reference.data ./ geo1t.R;
