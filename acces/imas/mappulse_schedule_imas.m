% test pulse 
%  metis_load ...
%  idx = imas_create('ids',1,1);
%  idx = imas_open('ids',1,1);
%  ps_caneva = ids_gen('pulse_schedule');
%  pulse = mappulse_schedule_imas(post,post.zerod,1);
%  ids_put(idx,'pulse_schedule',pulse);
%  test = ids_get(idx,'pulse_schedule');
%  zcompstruct(pulse,test);
%  imas_close(idx);
%
% convert data from metis to ids pulse_schedule (mapping)
%
function pulse = mappulse_schedule_imas(z0dstruct,data_zerod,plotonoff)

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

% test of input
if nargin < 3
	plotonoff = 0;
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


% Number of points of control for LCFS
nbpts = nbpts_max;  % number has been decreased for testing new UAL.

% creation of an empty structure to help to detect bug
pulse = ids_gen('pulse_schedule');

% backward compatible 
if isfield(pulse.flux_control.i_plasma.reference,'data')
    purge_reference_data = false;
else
    purge_reference_data = true;
    pulse = make_reference_data(pulse);
end    

% properties
pulse.ids_properties.homogeneous_time = 1;
switch z0dstruct.z0dinput.mode_exp
case -2
	pulse.ids_properties.comment = 'METIS evolution';
case -1
	pulse.ids_properties.comment = 'METIS simulation from scratch';	
case -3
	pulse.ids_properties.comment = 'METIS simulation from scenaio generator';	
case 0
	pulse.ids_properties.comment = 'METIS simulation from CRONOS data set';
case 1
	pulse.ids_properties.comment = 'METIS simulation for Tore Supra';	
case 2 
	pulse.ids_properties.comment = 'METIS simulation for JET';		
case 3 
	pulse.ids_properties.comment = 'METIS simulation for DIII-D';		
case 4
	pulse.ids_properties.comment = 'METIS simulation for ASDEX-U';		
case 5
	pulse.ids_properties.comment = 'METIS simulation for COMPASS';		
case 6
	pulse.ids_properties.comment = 'METIS simulation for EAST';		
case 7
	pulse.ids_properties.comment = 'METIS simulation for TCV';		
case 11
	pulse.ids_properties.comment = 'METIS simulation for Tore Supra (preparation)';	
case 12
	pulse.ids_properties.comment = 'METIS simulation for WEST';	
case 13
	pulse.ids_properties.comment = 'METIS simulation for WEST (preparation)';	
case 117 
	pulse.ids_properties.comment = 'METIS simulation from ITM data set';	
case 118 
	pulse.ids_properties.comment = 'METIS simulation from IMAS data set';	
case 119 
	pulse.ids_properties.comment = 'METIS simulation from code system data set';	
case 201
	pulse.ids_properties.comment = 'METIS simulation for ST-40';	
otherwise
	pulse.ids_properties.comment = sprintf('METIS simlation process by %s',getenv('USER'));	
end
pulse.ids_properties.source    = 'METIS';
if ispc
  pulse.ids_properties.provider  = getenv('USERNAME');
else
  pulse.ids_properties.provider  = getenv('USER');
end
pulse.ids_properties.creation_date = sprintf('%s (julian date in second = %f)',datestr(now,'dd-mmm-yyyy HH:MM:SS'),clock2julday);

% code info
[repository,commit,version] = metis_info_imas;
pulse.code.name        = 'METIS4IMAS';
pulse.code.version     = version;
pulse.code.commit      = commit;
pulse.code.repository  = repository;
% creation des donnees pour codeparam
data = z0dstruct.z0dinput.option;
tpn = tempname;
xml_write(tpn,data);
% lecture du fichier
fid = fopen(tpn,'r');
if fid > 0
	parameters = char(fread(fid,Inf,'char')');
	fclose(fid);
else
	parameters = 'unable to read parameters xml file';
end
delete(tpn);
pulse.code.parameters = parameters;
tol0d = z0dstruct.z0dinput.option.tol0d;
pulse.code.output_flag = fix(double((data_zerod.dw <= tol0d) & (data_zerod.dpfus <= tol0d) & ...
                                     (data_zerod.dini <= tol0d) & (data_zerod.diboot <= tol0d)));

% derived data
vloop = - 2 .* pi .* z0dxdt(z0dstruct.z0dinput.cons.flux,z0dstruct.z0dinput.cons.temps);
% extract time slices from z0dinput.cons, z0dinput.geo and sepa.
if (length(data_zerod.temps) == length(z0dstruct.z0dinput.cons.temps)) && all( data_zerod.temps == z0dstruct.z0dinput.cons.temps)
	% notthing to do
    if isfield(z0dstruct.z0dinput.exp0d,'Rsepa') && ~isempty(z0dstruct.z0dinput.exp0d.Rsepa)
        Rsepa = z0dstruct.z0dinput.exp0d.Rsepa;
        Zsepa = z0dstruct.z0dinput.exp0d.Zsepa;
    else
        Rsepa = [];
        Zsepa = [];
    end
    cons = z0dstruct.z0dinput.cons;
    geo = z0dstruct.z0dinput.geo;
else
	cons = z0dstruct.z0dinput.cons;
	temps = cons.temps;
	noms = fieldnames(cons);
	for k = 1:length(noms)
		var = cons.(noms{k});
		if length(var) == length(temps)
			cons.(noms{k}) = interp1(temps,var,data_zerod.temps,'nearest','extrap');
		end	
	end
        vloop = interp1(temps,vloop,data_zerod.temps,'nearest','extrap');
	
	geo = z0dstruct.z0dinput.geo;
	noms = fieldnames(geo);
	for k = 1:length(noms)
		var = geo.(noms{k});
		if length(var) == length(temps)
			geo.(noms{k}) = interp1(temps,var,data_zerod.temps,'nearest','extrap');
		end	
	end

	if isfield(z0dstruct.z0dinput.exp0d,'Rsepa') && ~isempty(z0dstruct.z0dinput.exp0d.Rsepa)
		temps = z0dstruct.z0dinput.exp0d.temps;
                if isempty(temps) || any(~isfinite(temps))
                       temps = z0dstruct.z0dinput.cons.temps;
		end
	        if all(isfinite(temps)) & (size(z0dstruct.z0dinput.exp0d.Rsepa,1) == length(temps))
			Rsepa = interp1(temps,z0dstruct.z0dinput.exp0d.Rsepa,data_zerod.temps,'nearest','extrap');
			Zsepa = interp1(temps,z0dstruct.z0dinput.exp0d.Zsepa,data_zerod.temps,'nearest','extrap');
                else
			Rsepa = [];
			Zsepa = [];
		end
 	else
		Rsepa = [];
		Zsepa = [];
	end
end
% homogeneous time
pulse.time = cons.temps;
vt = ones(size(pulse.time));

% for all METIS reference there is no enveloppe and all references are absolute
%reference_type = 1

% convert reference to ids data
% check internal config of METIS to know what go where
% z0dstruct.z0dinput.cons
%        iso: check internal metis config to split between n_H/n_D and n_T/n_D
%       xece: 
%         ip: 
%       flux: must be convertd in vloop
%       nbar: split between density (real) and gas puff (imag)
%      picrh: depending on heating scheme
%        plh: switch betwwen LHCD and second ECCD launcher 
%       pnbi: split real part (NBI1) and imag part (NBI2) and between H D T depending on ftnbi
%      ftnbi: see above
%      pecrh:  
%      hmore: 
%       zeff: 

% density and ISO
nbar    = real(cons.nbar);
gaspuff = imag(cons.nbar);

%
switch z0dstruct.z0dinput.option.gaz
    case 1
        n_t_over_n_d   = 0 .* vt;
        n_h_over_n_d   = 0 .* vt;
        n_b_over_n_h   = 0 .* vt;
        n_he3_over_n_d = 0 .* vt;
    case 2
        n_t_over_n_d   = 0 .* vt;
        n_b_over_n_h   = 0 .* vt;
        n_he3_over_n_d = 0 .* vt;
        switch z0dstruct.z0dinput.option.mino
            case 'H'
                n_h_over_n_d =  z0dstruct.z0dinput.option.cmin * vt;
            otherwise
                n_h_over_n_d =  cons.iso;
        end
    case 3
        n_t_over_n_d   =  cons.iso;
        n_b_over_n_h   = 0 .* vt;
        n_he3_over_n_d = 0 .* vt;
        switch z0dstruct.z0dinput.option.mino
            case 'H'
                n_h_over_n_d =  z0dstruct.z0dinput.option.cmin * vt;
            otherwise
                n_h_over_n_d = 0 .* vt;
        end
    case 5
        n_he3_over_n_d = cons.iso;
        n_t_over_n_d   = 0 .* vt;
        n_h_over_n_d   = 0 .* vt;
        n_b_over_n_h   = 0 .* vt;
    case 11
        n_he3_over_n_d = 0 .* vt
        n_t_over_n_d   = 0 .* vt;
        n_h_over_n_d   = 0 .* vt;
        n_b_over_n_h   = cons.iso;        
    case 4
        switch z0dstruct.z0dinput.option.mino
            case 'H'
                n_h_over_n_d =  z0dstruct.z0dinput.option.cmin * vt;
            otherwise
                n_h_over_n_d = 0 .* vt;
        end
        n_t_over_n_d   =  cons.iso;
        n_b_over_n_h   = 0 .* vt;
        n_he3_over_n_d = 0 .* vt;
end
%
pulse.density_control.n_e_line.reference_type = 1;
pulse.density_control.n_e_line.reference_name = 'nbar: line averaged density (m^-3)'; 
pulse.density_control.n_e_line.reference.data = nbar;
pulse.density_control.n_e_line.n_e_line_method.name  = 'nbar';
pulse.density_control.n_e_line.n_e_line_method.index = 3;
pulse.density_control.n_e_line.n_e_line_method.description = 'integral of a 1D core profile over rho_tor_norm up to the LCFS';

% data structure initialisation
species = pulse.density_control.valve{1}.species{1};
element = species.element{1};
%
switch z0dstruct.z0dinput.option.gaz
    case 1
        element.a   = 1;
        element.z_n = 1;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'H2';
        species.name  = 'H2';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [H2/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2;  % we assume injected gas is H2, D2, T2 or He4.
    case 2
        pulse.density_control.valve{2} = pulse.density_control.valve{1};
        element.a   = 1;
        element.z_n = 2;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'D2';
        species.name = 'D2';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [D2/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_h_over_n_d);  % we assume injected gas is H2, D2, T2 or He4.
        element.a   = 1;
        element.z_n = 1;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'H2';
        species.name = 'H2';
        species.fraction = 1;
        pulse.density_control.valve{2}.species{1} = species;
        pulse.density_control.valve{2}.flow_rate.reference_type = 1;
        pulse.density_control.valve{2}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [H2/s]';
        pulse.density_control.valve{2}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_h_over_n_d) .* n_h_over_n_d;  % we assume injected gas is H2, D2, T2 or He4.
    case 3
        pulse.density_control.valve{2} = pulse.density_control.valve{1};
        element.a   = 1;
        element.z_n = 2;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'D2';
        species.name = 'D2';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [D2/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_t_over_n_d);  % we assume injected gas is H2, D2, T2 or He4.
        element.a   = 1;
        element.z_n = 3;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'T2';
        species.name = 'T2';
        species.fraction = 1;
        pulse.density_control.valve{2}.species{1} = species;
        pulse.density_control.valve{2}.flow_rate.reference_type = 1;
        pulse.density_control.valve{2}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [T2/s]';
        pulse.density_control.valve{2}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_t_over_n_d) .* n_t_over_n_d;  % we assume injected gas is H2, D2, T2 or He4.
    case 5
        pulse.density_control.valve{2} = pulse.density_control.valve{1};
        element.a   = 1;
        element.z_n = 2;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'D2';
        species.name = 'D2';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [D2/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_he3_over_n_d);  % we assume injected gas is D2.
        element.a   = 3;
        element.z_n = 2;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'He3';
        species.name  = 'He3';
        species.fraction = 1;
        pulse.density_control.valve{2}.species{1} = species;
        pulse.density_control.valve{2}.flow_rate.reference_type = 1;
        pulse.density_control.valve{2}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [He3/s]';
        pulse.density_control.valve{2}.flow_rate.reference.data = gaspuff  ./ (1 + n_he3_over_n_d) .* n_he3_over_n_d;  % we assume injected gas is He3.
   case 4
        element.a   = 4;
        element.z_n = 2;
        species.element{1} = element;
        species.label = 'He4';
        species.name = 'He4';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [He4/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2;  % we assume injected gas is H2, D2, T2 or He4.
    case 11
        pulse.density_control.valve{2} = pulse.density_control.valve{1};
        element.a   = 1;
        element.z_n = 2;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'H2';
        species.name = 'H2';
        species.fraction = 1;
        pulse.density_control.valve{1}.species{1} = species;
        pulse.density_control.valve{1}.flow_rate.reference_type = 1;
        pulse.density_control.valve{1}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [H2/s]';
        pulse.density_control.valve{1}.flow_rate.reference.data = gaspuff ./  2 ./ (1 + n_b_over_n_h);  % we assume injected gas is H2.
        element.a   = 11;
        element.z_n = 5;
        species.element{1} = element;
        species.element{2} = element;
        species.label = 'B';
        species.name = 'B';
        species.fraction = 1;
        pulse.density_control.valve{2}.species{1} = species;
        pulse.density_control.valve{2}.flow_rate.reference_type = 1;
        pulse.density_control.valve{2}.flow_rate.reference_name = 'gaspuff: sum of gas puff flow rate [B11/s]';
        pulse.density_control.valve{2}.flow_rate.reference.data = gaspuff  ./ (1 + n_h_over_n_h) .* n_b_over_n_h;  % we assume injected gas is B11 (eevn if it is solid)
end
% 
pulse.density_control.zeff.reference_type = 1;
pulse.density_control.zeff.reference_name = 'zeff: line averaged effective charge computed without He ashes content'; 
pulse.density_control.zeff.reference.data = cons.zeff;
pulse.density_control.zeff.zeff_method.name  = 'zeff_bar';
pulse.density_control.zeff.zeff_method.index = 3;
pulse.density_control.zeff.zeff_method.description = 'average of a 1D core profile over rho_tor_norm up to the LCFS';
%
pulse.density_control.n_t_over_n_d.reference_type = 1;
pulse.density_control.n_t_over_n_d.reference_name = 'ISO (n_T/n_D): isotopic composition of the plasma ; ratio between tritum and deuterium density'; 
switch z0dstruct.z0dinput.option.gaz
    case 5
        if ~isfield(pulse.density_control,'n_he3_over_n_d')
            pulse.density_control.n_t_over_n_d.reference_name = 'ISO (n_He3/n_D): isotopic composition of the plasma ; ratio between helium 3 and deuterium density';
            pulse.density_control.n_t_over_n_d.reference.data = n_he3_over_n_d;           
        else
            pulse.density_control.n_t_over_n_d.reference.data = n_t_over_n_d;          
        end
    case 11
        if ~isfield(pulse.density_control,'n_b_over_n_h')
            pulse.density_control.n_t_over_n_d.reference_name = 'ISO (n_B/n_H): isotopic composition of the plasma ; ratio between bron 11 and hydrogen density';
            pulse.density_control.n_t_over_n_d.reference.data = n_b_over_n_h;
        else
            pulse.density_control.n_t_over_n_d.reference.data = n_t_over_n_d;          
        end
    otherwise
        pulse.density_control.n_t_over_n_d.reference.data = n_t_over_n_d;
end
%
pulse.density_control.n_h_over_n_d.reference_type = 1;
pulse.density_control.n_h_over_n_d.reference_name = 'ISO (n_H/n_D): isotopic composition of the plasma ; ratio between hydrogen and deuterium density'; 
pulse.density_control.n_h_over_n_d.reference.data = n_h_over_n_d;
%
pulse.density_control.n_he3_over_n_d.reference_type = 1;
pulse.density_control.n_he3_over_n_d.reference_name = 'ISO (n_He3/n_D): isotopic composition of the plasma ; ratio between helium 3 and deuterium density'; 
pulse.density_control.n_he3_over_n_d.reference.data = n_he3_over_n_d;
%
pulse.density_control.n_b_over_n_h.reference_type = 1;
pulse.density_control.n_b_over_n_h.reference_name = 'ISO (n_B/n_H): isotopic composition of the plasma ; ratio between bron 11 and hydrogen density'; 
pulse.density_control.n_b_over_n_h.reference.data = n_b_over_n_h;
% mode control
if isfield(pulse.density_control.mode,'data')
    pulse.density_control.mode.data = z0dstruct.z0dinput.option.neasser .* vt;
else
    pulse.density_control.mode = z0dstruct.z0dinput.option.neasser .* vt;
end

% EM boundary condition 
pulse.flux_control.i_plasma.reference_type = 1;
pulse.flux_control.i_plasma.reference_name = 'Ip: plasma current (A)';
pulse.flux_control.i_plasma.reference.data = cons.ip;
%
pulse.flux_control.loop_voltage.reference_type = 1;
pulse.flux_control.loop_voltage.reference_name = 'Vloop: plasma loop voltage (V, at plasma LCFS)';
pulse.flux_control.loop_voltage.reference.data = vloop;
% verification
if plotonoff 
	dt = diff(pulse.time);
	flux_rebuilt = -cumsum(pulse.flux_control.loop_voltage.reference.data(2:end) .* dt) ./ 2 ./ pi;
	flux_rebuilt = cat(1,0,flux_rebuilt);
	figure;plot(cons.temps,cons.flux - cons.flux(1),':r', pulse.time,flux_rebuilt,'b');drawnow
end
% mode control
if isfield(pulse.flux_control.mode,'data')
    pulse.flux_control.mode.data = data_zerod.asser;
else
    pulse.flux_control.mode = data_zerod.asser;    
end

% initalisation sub structure
inter   = pulse.nbi.unit{1};
species = inter.species{1};
element = species.element{1};
% NBI power is splitted in 2 injector and 3 species = 6 refrences
% nbi_1 = injector 1 hydrogen
% nbi_2 = injector 1 deuterium
% nbi_3 = injector 1 tritium
% nbi_4 = injector 2 hydrogen
% nbi_5 = injector 2 deuterium
% nbi_6 = injector 2 tritium
% nbi_7 = injector 1 boron
% nbi_8 = injector 2 boron
% slpit of references
pnbi1   = real(cons.pnbi);
pnbi2   = imag(cons.pnbi);
ftnbi1 = real(cons.ftnbi);
ftnbi2 = imag(cons.ftnbi);
% case of Hydrogen NBI in DT plasma
if isfield(z0dstruct.z0dinput.option,'forced_H_NBI') && (z0dstruct.z0dinput.option.forced_H_NBI ~= 0)
   gas_nbi = -1; 
else
   gas_nbi = z0dstruct.z0dinput.option.gaz; 
end

switch gas_nbi
    case -1
        pnbi{3} = 0 .* vt;
        pnbi{2} = 0 .* vt;
        pnbi{1} = pnbi1;
        pnbi{6} = 0 .* vt;
        pnbi{5} = 0 .* vt;
        pnbi{4} = pnbi2;
        pnbi{7} = 0 .* vt;
        pnbi{8} = 0 .* vt;
    case 3
        pnbi{1} = 0 .* vt;
        pnbi{2} = pnbi1 .* (1-ftnbi1);
        pnbi{3} = pnbi1 .* ftnbi1;
        pnbi{4} = 0 .* vt;
        pnbi{5} = pnbi2 .* (1-ftnbi2);
        pnbi{6} = pnbi2 .* ftnbi2;
        pnbi{7} = 0 .* vt;
        pnbi{8} = 0 .* vt;
    case 5
        pnbi{1} = 0 .* vt;
        pnbi{2} = pnbi1;
        pnbi{3} = 0 .* vt;
        pnbi{4} = 0 .* vt;
        pnbi{5} = pnbi2;
        pnbi{6} = 0 .* vt;
        pnbi{7} = 0 .* vt;
        pnbi{8} = 0 .* vt;
    case 11
        pnbi{2} = 0 .* vt;
        pnbi{1} = pnbi1 .* (1-ftnbi1);
        pnbi{7} = pnbi1 .* ftnbi1;
        pnbi{5} = 0 .* vt;
        pnbi{4} = pnbi2 .* (1-ftnbi2);
        pnbi{8} = pnbi2 .* ftnbi2;
        pnbi{3} = 0 .* vt;
        pnbi{6} = 0 .* vt;
    otherwise
        pnbi{3} = 0 .* vt;
        pnbi{2} = pnbi1 .* (1-ftnbi1);
        pnbi{1} = pnbi1 .* ftnbi1;
        pnbi{6} = 0 .* vt;
        pnbi{5} = pnbi2 .* (1-ftnbi2);
        pnbi{4} = pnbi2 .* ftnbi2;
        pnbi{7} = 0 .* vt;
        pnbi{8} = 0 .* vt;
end
for k=1:8
	switch	k
	case {1,4}
		A    = 1;
		elem = 'hydrogen';
		element.a   = 1;
		element.z_n = 1;
		species.element{1} = element;
		species.label = 'H';
		species.name = 'H';
		species.fraction = 1;
	case {2,5}
		A    = 2;
		elem = 'deuterium';
		element.a   = 1;
		element.z_n = 2;
		species.element{1} = element;
		species.label = 'D';
		species.name = 'D';
		species.fraction = 1;
	case {7,8}
		A    = 11;
		elem = 'boron';
		element.a   = 11;
		element.z_n = 5;
		species.element{1} = element;
		species.label = 'B';
		species.name = 'B';
		species.fraction = 1;
	otherwise
		A    = 3;
		elem = 'tritium';
		element.a   = 3;
		element.z_n = 1;
		species.element{1} = element;
		species.label = 'T';
		species.name = 'T';
		species.fraction = 1;
	end
	switch k
	    case {1,2,3,7}
		num_nbi = 1;
		einj    = z0dstruct.z0dinput.option.einj * vt;
	    otherwise
		num_nbi = 2;
		einj    = z0dstruct.z0dinput.option.einj2 * vt;
        end
        inter = pulse.nbi.unit{1};
        inter.species{1}  = species;
	inter.power.reference_name = sprintf('NBI power for injector %d : injection of %s [Z=%d, A=%d] (W)',num_nbi,elem,1,A);
	inter.power.reference_type = 1;
	inter.power.reference.data = pnbi{k};
	inter.energy.reference_name = sprintf('NBI fast neutral energy for injector %d : injection of %s [Z=%d, A=%d] (Hz)',num_nbi,elem,1,A);
	inter.energy.reference_type = 1;
	inter.energy.reference.data = einj;
	pulse.nbi.unit{k} = inter;
end


%  	% altenative implementation	
%  	pnbi{1}  = real(cons.pnbi);
%  	pnbi{2}   = imag(cons.pnbi);
%  	ftnbi{1} = real(cons.ftnbi);
%  	ftnbi{2}  = imag(cons.ftnbi);
%          for k = 1:2
%  		if k ==1
%  			num_nbi = 1;
%  			einj    = z0dstruct.z0dinput.option.einj * vt ./ phys.h;
%  		else
%  			num_nbi = 2;
%  			einj    = z0dstruct.z0dinput.option.einj2 * vt ./ phys.h;
%  		end
%  		inter.power.reference_name = sprintf('NBI power for injector %d (W)',num_nbi);
%  		inter.power.reference_type = 1;
%  		inter.power.reference.data = pnbi{k};
%  		inter.frequency.reference_name = sprintf('NBI fast neutral energy for injector %d (Hz)',num_nbi);
%  		inter.frequency.reference_type = 1;
%  		inter.frequency.reference.data = einj ;
%  		inter.phase.reference_name = sprintf('NBI fast neutral composition for injector %d (power fraction to tritium, either hydrogen depending on plasma composition)',num_nbi);
%  		inter.phase.reference_type = 1;
%  		inter.phase.reference.data = ftnbi{k} ;
%  		pulse.nbi.antenna{k} = inter;
%  	end


% ICRH: only one antenna - antenna phasing is not used in METIS
% initialisation sub structure
if isfield(pulse.ic,'antenna')
    inter = pulse.ic.antenna{1};
else
    inter = pulse.ic.launcher{1};
end
inter.power.reference_name = 'ICRH power(W)';
inter.power.reference_type = 1;
inter.power.reference.data = cons.picrh;
inter.frequency.reference_name = 'ICRH frequency (Hz)';
inter.frequency.reference_type = 1;
inter.frequency.reference.data = z0dstruct.z0dinput.option.freq .* 1e6 .* vt;
if isfield(pulse.ic,'antenna')
    pulse.ic.antenna{1} = inter;
else
    pulse.ic.launcher{1} = inter;
end
% ECRH: one or two antenna depending on METIS internal configuration
% initialisation sub structure
if isfield(pulse.ec,'antenna')
    inter = pulse.ec.antenna{1};
elseif isfield(pulse.ec,'launcher')
    inter = pulse.ec.launcher{1};
else
    inter = pulse.ec.beam{1};
end
% the deposition position is encoded in frequency using resonnance relation in midplane without doppler effect in perpendicular injection
switch z0dstruct.z0dinput.option.lhmode
    case 5
        
        % convert radial position in frequency
        freq1 = phys.e .* geo.b0 .* geo.R ./ (cons.xece .* geo.a  + geo.R)./ phys.me;
        freq2 = phys.e .* geo.b0 .* geo.R ./ (z0dstruct.z0dinput.option.xlh .* geo.a  + geo.R)./ phys.me;
        % first EC antenna
        inter.power.reference_name = sprintf('EC power for antenna %d (W)',1);
        inter.power.reference_type = 1;
        inter.power.reference.data = cons.pecrh;
        % for new dictionnary version ( <= 3.42)
        inter.power_launched.reference_name = sprintf('EC power for antenna %d (W)',1);
        inter.power_launched.reference_type = 1;
        inter.power_launched.reference.data = cons.pecrh;
        %inter.frequency.reference_name = sprintf('EC frequency (encoding deposition position) for antenna %d (Hz)',1);
        %inter.frequency.reference_type = 1;
        %inter.frequency.reference.data = freq1;
        inter.deposition_rho_tor_norm.reference_name = 'EC maximum power deposition position in Lao coordinate\n(make a small error compare to real rho_tor_norm in METIS;\nthe use of rho_tor_norm require equilibrium knowledge before the simulation).';
        inter.deposition_rho_tor_norm.reference_type = 1;
        inter.deposition_rho_tor_norm.reference.data = cons.xece;
        %inter.phase.reference_name = sprintf('EC poloidal angle deposition position for antenna %d (rad)',1);
        %inter.phase.reference_type = 1;
        %inter.phase.reference.data = z0dstruct.z0dinput.option.angle_ece ./ 180 .* pi .* vt;
        if isfield(pulse.ec,'antenna')
            pulse.ec.antenna{1} = inter;
            pulse.ec.antenna{2} = pulse.ec.antenna{1};
        elseif isfield(pulse.ec,'launcher')
            pulse.ec.launcher{1} = inter;
            pulse.ec.launcher{2} = pulse.ec.launcher{1};
        else
            pulse.ec.launcher{1} = inter;
            pulse.ec.launcher{2} = pulse.ec.beam{1};
        end
        % second EC antenna
        inter.power.reference_name = sprintf('EC power for antenna %d (W)',2);
        inter.power.reference_type = 1;
        inter.power.reference.data = cons.plh;
        % for new dictionnary version ( <= 3.42)
        inter.power_launched.reference_name = sprintf('EC power for antenna %d (W)',2);
        inter.power_launched.reference_type = 1;
        inter.power_launched.reference.data = cons.plh;
        %inter.frequency.reference_name = sprintf('EC frequency (encoding deposition position) for antenna %d (Hz)',2);
        %inter.frequency.reference_type = 1;
        %inter.frequency.reference.data = freq2;
        inter.deposition_rho_tor_norm.reference_name = 'EC maximum power deposition position in Lao coordinate\n(make a small error compare to real rho_tor_norm in METIS;\nthe use of rho_tor_norm require equilibrium knowledge before the simulation).';
        inter.deposition_rho_tor_norm.reference_type = 1;
        inter.deposition_rho_tor_norm.reference.data = z0dstruct.z0dinput.option.xlh * vt;
        %inter.phase.reference_name = sprintf('EC poloidal angle deposition position for antenna %d (rad)',2);
        %inter.phase.reference_type = 1;
        %inter.phase.reference.data = z0dstruct.z0dinput.option.angle_ece2 ./ 180 .* pi .* vt;
        if isfield(pulse.ec,'antenna')
            pulse.ec.antenna{2} = inter;
        elseif isfield(pulse.ec,'launcher')
            pulse.ec.launcher{2} = inter;
        else
            pulse.ec.beam{2} = inter;
        end
    otherwise
        
        % Launcher are not descibed in METIS, so launching_angle_pol and launching_angle_tor is not required.
        % convert radial position in frequency
        freq1 = phys.e .* geo.b0 .* geo.R ./ (cons.xece .* geo.a  + geo.R)./ phys.me;
        % ECCD antenna
        inter.power.reference_name = 'EC power (W)';
        inter.power.reference_type = 1;
        inter.power.reference.data = cons.pecrh;
        % for new dictionnary version ( <= 3.42)
        inter.power_launched.reference_name = 'EC power (W)';
        inter.power_launched.reference_type = 1;
        inter.power_launched.reference.data = cons.pecrh;
        %inter.frequency.reference_name = 'EC frequency for first hamonic @ maximum power deposition position';
        %inter.frequency.reference_type = 1;
        %inter.frequency.reference.data = freq1;
        inter.deposition_rho_tor_norm.reference_name = 'EC maximum power deposition position in Lao coordinate\n(make a small error compare to real rho_tor_norm in METIS;\nthe use of rho_tor_norm require equilibrium knowledge before the simulation).';
        inter.deposition_rho_tor_norm.reference_type = 1;
        inter.deposition_rho_tor_norm.reference.data = cons.xece;
        % this is not the right way to store this information
        %inter.phase.reference_name = 'EC poloidal angle deposition position';
        %inter.phase.reference_type = 1;
        %inter.phase.reference.data = z0dstruct.z0dinput.option.angle_ece ./ 180 .* pi .* vt;
        if isfield(pulse.ec,'antenna')
            pulse.ec.antenna{1} = inter;
            %pulse.ec.antenna{2} = pulse.ec.antenna{1};
        elseif isfield(pulse.ec,'launcher')
            pulse.ec.launcher{1} = inter;
            %pulse.ec.launcher{2} = pulse.ec.launcher{1};
        else
            pulse.ec.launcher{1} = inter;
            %pulse.ec.launcher{2} = pulse.ec.beam{1};
        end
%         % second EC antenna
%         inter.power.reference_name = sprintf('EC power for antenna %d (W)',2);
%         inter.power.reference_type = 1;
%         inter.power.reference.data = cons.plh;
%         %inter.frequency.reference_name = sprintf('EC frequency (encoding deposition position) for antenna %d (Hz)',2);
%         %inter.frequency.reference_type = 1;
%         %inter.frequency.reference.data = freq2;
%         inter.deposition_rho_tor_norm.reference_name = 'EC maximum power deposition position in Lao coordinate\n(make a small error compare to real rho_tor_norm in METIS;\nthe use of rho_tor_norm require equilibrium knowledge before the simulation).';
%         inter.deposition_rho_tor_norm.reference_type = 1;
%         inter.deposition_rho_tor_norm.reference.data = z0dstruct.z0dinput.option.xlh * vt;
%         %inter.phase.reference_name = sprintf('EC poloidal angle deposition position for antenna %d (rad)',2);
%         %inter.phase.reference_type = 1;
%         %inter.phase.reference.data = z0dstruct.z0dinput.option.angle_ece2 ./ 180 .* pi .* vt;
%         if isfield(pulse.ec,'antenna')
%             pulse.ec.antenna{2} = inter;
%         elseif isfield(pulse.ec,'launcher')
%             pulse.ec.launcher{2} = inter;
%         else
%             pulse.ec.beam{2} = inter;
%         end
%     otherwise
%         
%         % Launcher are not descibed in METIS, so launching_angle_pol and launching_angle_tor is not required.
%         % convert radial position in frequency
%         freq1 = phys.e .* geo.b0 .* geo.R ./ (cons.xece .* geo.a  + geo.R)./ phys.me;
%         % ECCD antenna
%         inter.power.reference_name = 'EC power (W)';
%         inter.power.reference_type = 1;
%         inter.power.reference.data = cons.pecrh;
%         %inter.frequency.reference_name = 'EC frequency for first hamonic @ maximum power deposition position';
%         %inter.frequency.reference_type = 1;
%         %inter.frequency.reference.data = freq1;
%         inter.deposition_rho_tor_norm.reference_name = 'EC maximum power deposition position in Lao coordinate\n(make a small error compare to real rho_tor_norm in METIS;\nthe use of rho_tor_norm require equilibrium knowledge before the simulation).';
%         inter.deposition_rho_tor_norm.reference_type = 1;
%         inter.deposition_rho_tor_norm.reference.data = cons.xece;
%         % this is not the right way to store this information
%         %inter.phase.reference_name = 'EC poloidal angle deposition position';
%         %inter.phase.reference_type = 1;
%         %inter.phase.reference.data = z0dstruct.z0dinput.option.angle_ece ./ 180 .* pi .* vt;
%         if isfield(pulse.ec,'antenna')
%             pulse.ec.antenna{1} = inter;
%         elseif isfield(pulse.ec,'launcher')
%             pulse.ec.launcher{1} = inter;
%         else
%             pulse.ec.beam{1} = inter;
%         end
end

% LHCD: zero or one antenna depending on METIS internal configuration
% phase will be used to encode n_//
% initialisation sub structure
if isfield(pulse.lh,'antenna')
    inter = pulse.lh.antenna{1};
else
    inter = pulse.lh.launcher{1};
end
switch z0dstruct.z0dinput.option.lhmode
    case 5
        % no LHCD in this case
    otherwise
        % LH antenna
        inter.power.reference_name = 'LH power (W)';
        inter.power.reference_type = 1;
        inter.power.reference.data = cons.plh;
        inter.frequency.reference_name = 'LH wave frequency (Hz)';
        inter.frequency.reference_type = 1;
        inter.frequency.reference.data = z0dstruct.z0dinput.option.freqlh .* 1e9 .* vt ;
        inter.n_parallel.reference_name = 'LH wave parralel index at launcher (main spectrum lob center)';
        inter.n_parallel.reference_type = 1;
        inter.n_parallel.reference.data = z0dstruct.z0dinput.option.npar0 .* vt ;
        if isfield(pulse.lh,'antenna')
            pulse.lh.antenna{1} = inter;
        else
            pulse.lh.launcher{1} = inter;
        end
end

%% compatibility with new version of IMAS
% for k=1:length(pulse.ec.antenna)
%     pulse.ec.beam{k} = pulse.ec.antenna{k};
% end
% for k=1:length(pulse.ec.antenna)
%     pulse.ec.launcher{k} = pulse.ec.antenna{k};
% end
% for k=1:length(pulse.lh.antenna)
%     pulse.lh.launcher{k} = pulse.lh.antenna{k};
% end


%% HMORE in betan and li ; hmore = betan ./ 4 ./ li
pulse.flux_control.li_3.reference_type = 1;
pulse.flux_control.li_3.reference_name = 'li_3 targeted; will be used in the control of betan (confinement factor = beta_normal / 4 / li_3)';
pulse.flux_control.li_3.reference.data = (2/3) .* ones(size(cons.temps));
pulse.flux_control.beta_normal.reference_type = 1;
pulse.flux_control.beta_normal.reference_name = 'Beta_N targeted; will be used in the control of betan (confinement factor = beta_normal / 4 / li_3)';
pulse.flux_control.beta_normal.reference.data =  (2/3) .* 4 .* cons.hmore;

%  initialisation sub structure
boundary_outline = pulse.position_control.boundary_outline{1};
% z0dstruct.z0dinput.geo
%         a: minor radius 
%         R: major raidus
%         K: averaged up/down elongation
%         d: averaged up/down triangularities
%        b0: vacuum toroidal magnetic field @ R
%        z0: vertical position (of magnetic axis and LCFS shift)
%   
if ~isempty(Rsepa) && ~isempty(Zsepa)
    control =sepa_moments(Rsepa,Zsepa + geo.z0 * ones(1,size(Zsepa,2)),nbpts);
    pulse.position_control.magnetic_axis.r.reference_name = 'Magnetic axis major radius (m)';
    pulse.position_control.magnetic_axis.r.reference_type = 1;
    pulse.position_control.magnetic_axis.r.reference.data = control.R0 + data_zerod.d0;
    pulse.position_control.magnetic_axis.z.reference_name = 'Magnetic axis vertical shift (m)';
    pulse.position_control.magnetic_axis.z.reference_type = 1;
    pulse.position_control.magnetic_axis.z.reference.data = control.z0;
    pulse.position_control.geometric_axis.r.reference_name = 'LCFS major radius (m)';
    pulse.position_control.geometric_axis.r.reference_type = 1;
    pulse.position_control.geometric_axis.r.reference.data = control.R0;
    pulse.position_control.geometric_axis.z.reference_name = 'LCFS vertical shift (m)';
    pulse.position_control.geometric_axis.z.reference_type = 1;
    pulse.position_control.geometric_axis.z.reference.data = control.z0_geo;
    pulse.position_control.minor_radius.reference_name    = 'Minor radius of the plasma (m)';
    pulse.position_control.minor_radius.reference_type    = 1;
    pulse.position_control.minor_radius.reference.data    = control.a;
    pulse.position_control.elongation.reference_name            = 'Vertical elongation of the plasma';
    pulse.position_control.elongation.reference_type            = 1;
    pulse.position_control.elongation.reference.data            = control.K;
    pulse.position_control.elongation_upper.reference_name      = 'Vertical elongation of upper part of the plasma (Z >=  vertical_shift)';
    pulse.position_control.elongation_upper.reference_type      = 1;
    pulse.position_control.elongation_upper.reference.data      = control.Ku;
    pulse.position_control.elongation_lower.reference_name      = 'Vertical elongation of the lower part of the plasma (Z <=  vertical_shift)';
    pulse.position_control.elongation_lower.reference_type      = 1;
    pulse.position_control.elongation_lower.reference.data      = control.Kl;
    pulse.position_control.triangularity.reference_name            = 'Averaged triangularity of the plasma';
    pulse.position_control.triangularity.reference_type            = 1;
    pulse.position_control.triangularity.reference.data            = control.d;
    pulse.position_control.triangularity_upper.reference_name      = 'Triangularity of upper part of the plasma (Z >=  vertical_shift)';
    pulse.position_control.triangularity_upper.reference_type      = 1;
    pulse.position_control.triangularity_upper.reference.data      = control.du;
    pulse.position_control.triangularity_lower.reference_name      = 'Triangularity of the lower part of the plasma (Z <=  vertical_shift)';
    pulse.position_control.triangularity_lower.reference_type      = 1;
    pulse.position_control.triangularity_lower.reference.data      = control.dl;
    % x_point
    if length(control.x_point) > 0
        for lz=1:length(control.x_point)
            x_point = pulse.position_control.x_point{1};
            x_point.r.reference_name  = sprintf('Radial position of the X-point %d (m)',lz);
            x_point.r.reference_type  = 1;
            x_point.r.reference.data = control.x_point{lz}.r;
            x_point.z.reference_name  = sprintf('Vertical position of the X-point %d (m)',lz);
            x_point.z.reference_type  = 1;
            x_point.z.reference.data = control.x_point{lz}.z;
            pulse.position_control.x_point{lz} = x_point;
        end
    end
    for lz=1:size(control.Rsepa,2)
      boundary_outline  = pulse.position_control.boundary_outline{1};
      boundary_outline.r.reference_name  = 'LCFS list of points: radial position (m)';
      boundary_outline.r.reference_type  = 1;
      boundary_outline.r.reference.data  = control.Rsepa(:,lz);
      boundary_outline.z.reference_name  = 'LCFS list of points: vertical position (m)';
      boundary_outline.z.reference_type  = 1;
      boundary_outline.z.reference.data  = control.Zsepa(:,lz);	  
      pulse.position_control.boundary_outline{lz} = boundary_outline;
    end
    disp('Pulse_schedule: position defined by LCFS');
else
    pulse.position_control.magnetic_axis.r.reference_name = 'Magnetic axis major radius (m)';
    pulse.position_control.magnetic_axis.r.reference_type = 1;
    pulse.position_control.magnetic_axis.r.reference.data = geo.R + data_zerod.d0;
    pulse.position_control.magnetic_axis.z.reference_name = 'Magnetic axis vertical shift (m)';
    pulse.position_control.magnetic_axis.z.reference_type = 1;
    pulse.position_control.magnetic_axis.z.reference.data = geo.z0;
    pulse.position_control.geometric_axis.r.reference_name = 'LCFS major radius (m)';
    pulse.position_control.geometric_axis.r.reference_type = 1;
    pulse.position_control.geometric_axis.r.reference.data = geo.R;
    pulse.position_control.geometric_axis.z.reference_name = 'LCFS vertical shift (m)';
    pulse.position_control.geometric_axis.z.reference_type = 1;
    pulse.position_control.geometric_axis.z.reference.data = geo.z0;
    pulse.position_control.minor_radius.reference_name    = 'Minor radius of the plasma (m)';
    pulse.position_control.minor_radius.reference_type    = 1;
    pulse.position_control.minor_radius.reference.data    = geo.a;
    pulse.position_control.elongation.reference_name      = 'Vertical elongation of the plasma';
    pulse.position_control.elongation.reference_type      = 1;
    pulse.position_control.elongation.reference.data      = geo.K;
    pulse.position_control.triangularity.reference_name   = 'Averaged triangularity of the plasma';
    pulse.position_control.triangularity.reference_type   = 1;
    pulse.position_control.triangularity.reference.data   = geo.d;
    disp('Pulse_schedule: position defined by moments');  
end
%
pulse.tf.b_field_tor_vacuum_r.reference_name = 'Vacuum magnetic rigidity (T.m)';
pulse.tf.b_field_tor_vacuum_r.reference_type = 1;
pulse.tf.b_field_tor_vacuum_r.reference.data = geo.R .* geo.b0;

% backward compatible 
if  purge_reference_data
    pulse = make_reference_without_data(pulse);
end    



function control =sepa_moments(Rsepa,Zsepa,nbpts)

if nargin < 3
  nbpts = nbpts_max;
end

% memory allocation
vnan        = NaN .* ones(size(Rsepa,1),1);
control.R0  = vnan;
control.a   = vnan;
control.z0  = vnan; 
control.Ku  = vnan; 
control.Kl  = vnan; 
control.d   = vnan; 
control.du  = vnan; 
control.dl  = vnan; 
control.x_point{1}.r = vnan;
control.x_point{1}.z = vnan;   
control.x_point{2}.r = vnan;
control.x_point{2}.z = vnan;   
if size(Rsepa,2) <= nbpts
  control.Rsepa = Rsepa;
  control.Zsepa = Zsepa;
else
  control.Rsepa = vnan * ones(1,nbpts);
  control.Zsepa = vnan * ones(1,nbpts);
end

% calcul des moments
for l=1:size(Rsepa,1)
    % calcul de R0 et Z0
    % recalcul des parametres sur le vecteur final
    rmin  = min(Rsepa(l,:),[],2);
    rmax  = max(Rsepa(l,:),[],2);
    a = 0.5 .* (rmax - rmin);
    R0 = 0.5 .* (rmax + rmin);
    control.R0(l)  = R0;
    control.a(l)   = a;
    zmin  = min(Zsepa(l,:),[],2);
    zmax  = max(Zsepa(l,:),[],2);
    control.K(l)    = (zmax - zmin) ./ 2 ./ a;
    rzmax = Rsepa(l,min(find(Zsepa(l,:) == zmax)));
    rzmin = Rsepa(l,min(find(Zsepa(l,:) == zmin)));
    z0    = Zsepa(l,min(find(Rsepa(l,:) == rmax)));
    % IMAS definition
    control.z0_geo(l)    = (zmax + zmin) ./ 2;
    control.z0(l)  = z0;
    control.Ku(l)    = (zmax - z0)  ./ a;
    control.Kl(l)    = (z0 - zmin)  ./ a;
    control.d(l)     = abs(rzmax + rzmin -  2 .* R0) ./ 2 ./ a;
    control.du(l)    = abs(rzmax - R0) ./ a;
    control.dl(l)    = abs(rzmin - R0) ./ a;
    % Xpoint detection
    indh     = find(Zsepa(l,:) == zmax,1);
    indh     = cat(2,indh - 3, indh - 2,indh - 1,indh,indh + 1,indh + 2,indh + 3); 
    indh     = mod(indh-1,size(Zsepa,2))+1;
    rh       = Rsepa(l,indh);
    zh       = Zsepa(l,indh);
    ph       = polyfit(rh,zh,2);
    eh       = sqrt(mean((zh - polyval(ph,rh)).^ 2)) ./ (max(rh) - min(rh));
    indl     = find(Zsepa(l,:) == zmin,1);
    indl     = cat(2,indl - 3, indl - 2,indl - 1,indl,indl + 1,indl + 2,indl + 3);
    indl     = mod(indl-1,size(Zsepa,2))+1;
    rl       = Rsepa(l,indl);
    zl       = Zsepa(l,indl);
    pl       = polyfit(rl,zl,2);
    el       = sqrt(mean((zl - polyval(pl,rl)).^ 2)) ./ (max(rl) - min(rl));
    if el > 2e-2
      indlz  = find(zl == min(zl),1);
      control.x_point{1}.r(l) = rl(indlz);
      control.x_point{1}.z(l) = zl(indlz);
    end
    if eh > 2e-2
      indhz  = find(zh == max(zh),1);
      control.x_point{2}.r(l) = rh(indhz);
      control.x_point{2}.z(l) = zh(indhz);
    end
    % loop on time slices to extract nbpts points for LCFS
    if size(Rsepa,2) >  nbpts
	  indice = union(find(Rsepa(l,:) == min(Rsepa(l,:))),find(Rsepa(l,:) == max(Rsepa(l,:))));
	  indice = union(indice,find(Zsepa(l,:) == min(Zsepa(l,:))));
	  indice = union(indice,find(Zsepa(l,:) == max(Zsepa(l,:))));
	  indice = union(union(indice,indice - 1), indice + 1);
	  if isfinite(control.x_point{1}.r(l))
		indice_x = find((control.x_point{1}.r(l) == Rsepa(l,:)) & (control.x_point{1}.z(l) == Zsepa(l,:)),1);
		indice   = union(indice,indice_x);
	  end
	  if isfinite(control.x_point{2}.r(l))
		indice_x = find((control.x_point{2}.r(l) == Rsepa(l,:)) & (control.x_point{2}.z(l) == Zsepa(l,:)),1);
		indice   = union(indice,indice_x);
	  end
	  indice(indice < 1) = [];
	  indice(indice > size(Rsepa,2)) = [];
	  while length(indice) < nbpts
	    reste  = setxor(1:size(Rsepa,2),indice);
	    lcomp  = nbpts - length(indice);
	    indice = union(indice,reste(ceil(linspace(1,length(reste),lcomp))));
	  end
	  indice = sort(indice);
	  control.Rsepa(l,:) = Rsepa(l,indice);
	  control.Zsepa(l,:) = Zsepa(l,indice);
	  %figure(31);plot(Rsepa(l,:),Zsepa(l,:),'r',control.Rsepa(l,:),control.Zsepa(l,:),'.b');
	  %drawnow
    end
end

% set nbpts max depending on IMAS version
function nbpts = nbpts_max

% default value
nbpts = 49;

% test IMAS version
if isappdata(0,'UAL_DATAVERSION')
   rep = getappdata(0,'UAL_DATAVERSION'); 
   [v,r] = strtok( rep,'.');
   if str2double(v) > 3
      nbpts = 301;
   else
       [v,r] = strtok( r,'.');
       if str2double(v) > 20
           nbpts = 301;
       end
   end
 
end


