function [option_feeqs_jt60sa,fbe_lcfs] = compute_jt60sa_fbe_inverse_evol(option_feeqs_jt60sa)

% test if FEEQS is in the path
set_feeqs_path_jt60sa;
fbe_lcfs = [];

% GUI generation
if nargin < 1

    valeur.metisfilename	= '';
    type.metisfilename          = 'string';                    
    borne.metisfilename         = '';  
    defaut.metisfilename        = '';                
    info.metisfilename          = 'path and name of METIS simulation file used for the computation (if empty, a file selector will be open at the next step)';

    valeur.number_of_time_slices    = 101;   
    type.number_of_time_slices      = 'float';
    borne.number_of_time_slicese    = [3,Inf];  
    defaut.number_of_time_slices    = 101;
    info.number_of_time_slices      = 'Number of time slices extracted from METIS simulation at which inverse equilibrium would be computed';

    valeur.Vol_poly_degree    = 2;   
    type.Vol_poly_degree      = 'float';
    borne.Vol_poly_degree     = [2,Inf];  
    defaut.Vol_poly_degree    = 2;
    info.Vol_poly_degree      = 'Degree of polynomial representation of coil voltage waveforms.\nUsed to smooth the feedforward and speed up the computation.\nif equal to number of time step, switch to local times value (no more smoothing).';

    valeur.first_time     = 0;   
    type.first_time       = 'float';
    borne.first_time      = [0,100];  
    defaut.first_time     = 0;
    info.first_time      = 'First time of the METIS simulation at which inverse equilibrium would be computed';

    valeur.last_time     = 100;   
    type.last_time       = 'float';
    borne.last_time      = [1,Inf];  
    defaut.last_time     = 100;
    info.last_time      = 'Last time of the METIS simulation for which inverse equilibrium would be computed';

    valeur.reinforce     = 30;   
    type.reinforce       = 'float';
    borne.reinforce      = [0,100];  
    defaut.reinforce     = 30;
    info.reinforce      = 'weigths of extrema in LCFS (to have a better match on X-points and RIG and ROG)';
    
    valeur.mesh_type = 'standard';     
    type.mesh_type   = 'string';
    borne.mesh_type  = {'test','standard','fine'};  
    defaut.mesh_type = 'standard';
    info.mesh_type   = 'number of elements in the 2D mesh (test is reserved for testing and must be not used for production)';    
    
    valeur.run_name 	    = '';
    type.run_name           = 'string';                    
    borne.run_name          = '';  
    defaut.run_name         = '';                
    info.run_name           = 'run_name is added to name of save files (can be leave empty)';
    
    valeur.current_margins     = 0;   
    type.current_margins       = 'float';
    borne.current_margins      = [-0.2,0.2];  
    defaut.current_margins     = 0;
    info.current_margins       = 'margin to be keep to coil current limits (on absolute value, if < 0, allows to go above the limit);\n Pay attention to be compatible with plasma initialisation value;\n Applied only to EF coils';
    
    valeur.voltage_margins     = 0;   
    type.voltage_margins       = 'float';
    borne.voltage_margins      = [-0.2,0.2];  
    defaut.voltage_margins     = 0;
    info.voltage_margins       = 'margin to be keep to coil voltage limits (on absolute value, if < 0, allows to go above the limit);\n Pay attention to be compatible with plasma initialisation value;\n Applied on all coils';
    
    valeur.mode_init_ps_vv = 'METIS Vloop';     
    type.mode_init_ps_vv   = 'string';
    borne.mode_init_ps_vv  = {'METIS Vloop', 'Precribed','zero','undetermined'};  
    defaut.mode_init_ps_vv = 'METIS Vloop';
    info.mode_init_ps_vv   = 'Methode to initialise currents in passive strcutures at first time of computation (first_time):\nif = Metis Vloop, used vloop from METIS simulation;\nif = Prescribed, used parameter value_init_ps_vv;\nif = zero, set initial current in passive structures to zero;\nif = undetermined, no information are provided on current in passive structures';    
    
    valeur.value_init_ps_vv     = 0;   
    type.value_init_ps_vv       = 'float';
    borne.value_init_ps_vv      = [-100,100];  
    defaut.value_init_ps_vv     = 0;
    info.value_init_ps_vv       = 'Value of loop voltage in passive structure elements at first time of computation (V)';
    
    valeur.init_mode = 'full';     
    type.init_mode   = 'string';
    borne.init_mode  = {'current','full','simple'};  
    defaut.init_mode = 'full';
    info.init_mode   = 'Methode used to initialise first equilibrium coil current and poloidal flux map:\nif = simple, start from pre-magnetisation current correct of poloidal flux consumption computed by METIS and simple analytical poloidal flux map;\nif = current, used fast fixed boundary identification mode to compute first coil current values and simple analytical poloidal flux map;\nif = full, used fast fixed boundary identification mode to compute first coil current values and interpolated on mesh poloidal flux map from  fast fixed boundary identification mode.';    
 
    valeur.plotonoff    = 2;      
    type.plotonoffs     = 'integer';
    borne.plotonoff     = {0,1,2};  
    defaut.plotonoff    = 2;
    info.plotonoff      = 'Set level of displayed graphs:  0 = no graphs; 1 = results only; 2 = all graphs (debug)';
    
    valeur.movie_onoff  = 0;      
    type.movie_onoff    = 'integer';
    borne.movie_onoff   = {0,1,2};  
    defaut.movie_onoff  = 0;
    info.movie_onoff    = 'if = 1, record movie made of 2D equilibrium time slice graphs; if  = 2, record movie made of 2D equilibrium time slice graphs and save a .fig file for each time slice';
    
    valeur.premag_file  = '';
    type.premag_file    = 'string';                    
    borne.premag_file   = '';  
    defaut.premag_file  = '';                
    info.premag_file    = 'if available, name of file containing premagnetisation data (can be leave empty; if = workspace, reads data structure premagnetisation_jt60sa in Matlab workspace)';
    
    valeur.flux_constant_mode  = 'manual';      
    type.flux_constant_mode    = 'string';
    borne.flux_constant_mode   = {'manual','natural','start of flux','balanced','end of flux','xpoint','xpoint with margin','auto'};  
    defaut.flux_constant_mode  = 'manual';
    info.flux_constant_mode   = ['This parameter allows to choose the way poloidal flux offset used to match METIS current diffusion LCFS poloidal flux and FEEQS.M LCFS poloidal flux is computed:\n', ...
         'if = manual,offset is prescribed by the user\n', ...
         'in other offset is automatically computed taking into account flux consumption estimation for breakdown computed in METIS and flux leakage computed with the help of current diffusion equation between end of breakdown and  first_time compare to flux computed by premagnetisation tool:', ...
         '    case natural: compute CS2 current, for first computed time proportionnal to METIS poloidal flux consumption\n', ...
         '    case start_of_flux: offset is computed to have maximum available flux for the first plasma, taking into account breakdonw flux\n', ...
         '    case balanced: offset is comuted to have flux margin for ramp-up and ramp down if possible\n', ...
         '    case end_of_flux:offset is computed to have minimum flux value (negative) at the end of rampdwon\n', ...
         '    case xpoint: set offset to be on the maximum allowed flux given by in Urano paper at X-point formattion\n', ...
         '    case xpoint with margin: set offset to be near the maximum allowed flux ( - 4 Wb) given by in Urano paper at X-point formattion\n', ...
         '    case auto: old automatic computation method of the offset to match prescribed CS2 current. Must be not any more be used.'];
    
    valeur.voltage_limits = 'flattop';     
    type.voltage_limits   = 'string';
    borne.voltage_limits  = {'ramp-up','flattop'};  
    defaut.voltage_limits = 'flattop';
    info.voltage_limits   = 'Voltage limits range:\nif = ramp-up, use extended value allowed by booster and discharge resistor;\nif = flattop, used limited range allowed by steady state power supply';    

    valeur.flux_offset    = 17.4;   
    type.flux_offset      = 'float';
    borne.flux_offset     = [-50,50];  
    defaut.flux_offset    = 17.4;
    info.flux_offset      = 'poloidal flux offset used to match METIS current diffusion LCFS poloidal flux and FEEQS.M LCFS poloidal flux (Wb)';
   
    valeur.breakdown_flux   = 0.3;   
    type.breakdown_flux     = 'float';
    borne.breakdown_flux    = [0, 20];  
    defaut.breakdown_flux   = 0.3;
    info.breakdown_flux   = 'poloidal flux consumed for plasma break-down and burn-through (Wb);\nif METIS breakdown model is turn on, take the maximum of this value and METIS value.';
   
    valeur.init_voltages  = 'METIS Vloop';      
    type.init_voltages    = 'string';
    borne.init_voltages   = {'zero','METIS Vloop','resitive'};  
    defaut.init_voltages  = 'METIS Vloop';
    info.init_voltages    = 'Method used to compute initial guess of voltages provided by power supplies for optimizer:\nif = zero, start with null values;\nif = resistive, start with voltages equal to resistive losses;\nif = METIS Vloop, use METIS loop voltage to estimate first guess for voltages provided by power supplies';

    valeur.sigma_PS    = 0.72e-6;   
    type.sigma_PS      = 'float';
    borne.sigma_PS     = [1e-16,1];  
    defaut.sigma_PS    = 0.72e-6;
    info.sigma_PS      = 'Resistivity of passive structure material (Ohm.m)';
 
    valeur.sigma_VV    = 0.72e-6;   
    type.sigma_VV      = 'float';
    borne.sigma_VV     = [1e-16,1];  
    defaut.sigma_VV    = 0.72e-6;
    info.sigma_VV      = 'Resistivity of vacuum vessel material (Ohm.m)';
 
    valeur.strike_point  = 'off';      
    type.strike_point    = 'string';
    borne.strike_point   = {'on','off'};  
    defaut.strike_point  = 'off';
    info.strike_point    = 'For diverted plasma: if = off, no constraint on strike point position;if = on, add constraints on strike point positions';

    valeur.R_strike_point    = '1.92, 2.64';  
    type.R_strike_point      = 'string';
    borne.R_strike_point     = '';  
    defaut.R_strike_point    = '1.92, 2.64';
    info.R_strike_point      = 'list of strike point radial positions (m, quote separated)';

    valeur.Z_strike_point    = '-2.48, -2.83';  
    type.Z_strike_point      = 'string';
    borne.Z_strike_point     = '';  
    defaut.Z_strike_point    = '-2.48, -2.83';
    info.Z_strike_point      = 'list of strike point vertical positions (m, quote separated)';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to be put in weight_struct %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    valeur.weight_I2     = 1e-14;   
    type.weight_I2       = 'float';
    borne.weight_I2      = [1e-21,1e-4];  
    defaut.weight_I2     = 1e-14;
    info.weight_I2       = 'Weight or regularisation term son sum(Icoils^2) in first inverse mode';

    valeur.factor_w_pfcur     = 1;   
    type.factor_w_pfcur       = 'float';
    borne.factor_w_pfcur      = [0,10];  
    defaut.factor_w_pfcur     = 1;
    info.factor_w_pfcur       = 'relative weight on penalisation term for current limits in coils; used only for first inverse equilibrium computation';
    
    valeur.ef5_w    = 3;   
    type.ef5_w      = 'float';
    borne.ef5_w     = [0,10];  
    defaut.ef5_w    = 3;
    info.ef5_w      = 'penalisation of current in EF5 coil to prevent it be used as main coil to make X-point (power of 10); used only for first inverse equilibrium computation';
    
    valeur.weight_V2     = 1e-14;   
    type.weight_V2       = 'float';
    borne.weight_V2      = [1e-21,1e-4];  
    defaut.weight_V2     = 1e-14;
    info.weight_V2       = 'Weight or regularisation term son sum(Voltages^2) in inverse evolutive mode';

    valeur.weight_LCFS     = 1;   
    type.weight_LCFS       = 'float';
    borne.weight_LCFS      = [1e-15,10];  
    defaut.weight_LCFS     = 1;
    info.weight_LCFS       = 'relative weight on penalisation term on LCFS desired points in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';
    
    valeur.weight_psibd     = 1;   
    type.weight_psibd       = 'float';
    borne.weight_psibd      = [0,10];  
    defaut.weight_psibd     = 1;
    info.weight_psibd       = 'relative weight on penalisation term on prescribed poloidal flux at LCFS in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';    
    
    valeur.weight_PS     = 1;   
    type.weight_PS       = 'float';
    borne.weight_PS      = [0,10];  
    defaut.weight_PS     = 1;
    info.weight_PS       = 'relative weight on penalisation term on integral(J_PS^2*dS) on passive structure elements in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';    
    
    valeur.weight_VV     = 1;   
    type.weight_VV       = 'float';
    borne.weight_VV      = [0,10];  
    defaut.weight_VV     = 1;
    info.weight_VV       = 'relative weight on penalisation term on integral(J_VV^2*dS) on vacuum vessel in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';    
    
    valeur.weight_current     = 1;   
    type.weight_current       = 'float';
    borne.weight_current      = [0,10];  
    defaut.weight_current     = 1;
    info.weight_current       = 'relative weight on penalisation term for coil current limits compiliance in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';    

    valeur.weight_Voltage     = 1;   
    type.weight_Voltage       = 'float';
    borne.weight_Voltage      = [0,10];  
    defaut.weighweight_Voltaget_Voltage     = 1;
    info.weight_current       = 'relative weight on penalisation term for coil voltage limits compiliance in inverse evolutive mode\n(must be set to zero to switch off this term in cost function)';    

    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

    option_feeqs_jt60sa.valeur     = valeur;
    option_feeqs_jt60sa.type       = type;
    option_feeqs_jt60sa.borne      = borne;
    option_feeqs_jt60sa.defaut     = defaut;
    option_feeqs_jt60sa.info       = info;
    option_feeqs_jt60sa.interface  = interface;
    
    option_feeqs_jt60sa.description = 'Computation of pre magnetisation before breakdown for JT-60SA';   % description (une ligne) de la fonction

    option_feeqs_jt60sa.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    option_feeqs_jt60sa.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    option_feeqs_jt60sa.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

    % end of GUI form declaration
    return

end

% strike points informations
switch option_feeqs_jt60sa.strike_point
case 'on'
    Rsp = str2num(option_feeqs_jt60sa.R_strike_point);
    Zsp = str2num(option_feeqs_jt60sa.Z_strike_point);
    strike_point = cat(2,Rsp(:),Zsp(:));
otherwise
    strike_point = [];
end

% select file is needed
[pathstr,name,ext] = fileparts(option_feeqs_jt60sa.metisfilename);
if isempty(ext)
  option_feeqs_jt60sa.metisfilename = sprintf('%s.mat',option_feeqs_jt60sa.metisfilename);
end
if isempty(option_feeqs_jt60sa.metisfilename) || ~exist(option_feeqs_jt60sa.metisfilename)
    [FileName,PathName] = uigetfile('*.mat','Select a METIS file');
    if isempty(FileName) || isnumeric(FileName)
	disp('Call to FEEQS.M has been canceled');
	return	
    end
    option_feeqs_jt60sa.metisfilename = fullfile(PathName,FileName);
end
option_feeqs_jt60sa.premag_file = strtrim(option_feeqs_jt60sa.premag_file);
switch option_feeqs_jt60sa.premag_file
  case 'workspace'
      try 
	  premagnetisation_jt60sa = evalin('base','premagnetisation_jt60sa');
      catch
	  disp('data premagnetisation_jt60sa do not exist in workspace');
	  option_feeqs_jt60sa.premag_file = '?';
      end
end
switch option_feeqs_jt60sa.premag_file
  case 'workspace'
      % nothing else
  otherwise
      [pathstr,name,ext] = fileparts(option_feeqs_jt60sa.premag_file);
      if isempty(ext) && ~isempty(option_feeqs_jt60sa.premag_file)
	option_feeqs_jt60sa.premag_file = sprintf('%s.mat',option_feeqs_jt60sa.premag_file);
      end
      if isempty(option_feeqs_jt60sa.premag_file) 
	  % no premagnetisation
	  %fprintf('Using default configuration for premagnetisation\n');
      elseif ~exist(option_feeqs_jt60sa.premag_file)
	  [FileName,PathName] = uigetfile('*.mat','Select a premagnetisation file');
	  if isempty(FileName) || isnumeric(FileName)
	      % no file selected
	      option_feeqs_jt60sa.premag_file = '';
	  else
	      option_feeqs_jt60sa.premag_file = fullfile(PathName,FileName);
	  end
      end
      % load premagnetisation if needed
      if ~isempty(option_feeqs_jt60sa.premag_file)
	  premagnetisation_jt60sa = load(option_feeqs_jt60sa.premag_file);
	  fprintf('Premagnetisation file %s has been loaded\n',option_feeqs_jt60sa.premag_file);
      else
	  premagnetisation_jt60sa = [];
	  fprintf('Using default configuration for premagnetisation\n');
      end
end
% update parameters in workspace
zassignin('base','premagnetisation_jt60sa',premagnetisation_jt60sa);

% weight structure
weight_struct.weight_I2         = option_feeqs_jt60sa.weight_I2;
weight_struct.weight_V2         = option_feeqs_jt60sa.weight_V2;
weight_struct.f_ef5w            = option_feeqs_jt60sa.ef5_w;
weight_struct.factor_w_pfcur    = option_feeqs_jt60sa.factor_w_pfcur;
weight_struct.weight_LCFS       = option_feeqs_jt60sa.weight_LCFS;
weight_struct.weight_PS         = option_feeqs_jt60sa.weight_PS;
weight_struct.weight_VV         = option_feeqs_jt60sa.weight_VV;
weight_struct.weight_current    = option_feeqs_jt60sa.weight_current;
weight_struct.weight_Voltage    = option_feeqs_jt60sa.weight_Voltage;
weight_struct.weight_psibd      = option_feeqs_jt60sa.weight_psibd;

%%%%%%%%% security on polynomial order
option_feeqs_jt60sa.Vol_poly_degree = max(2,min(option_feeqs_jt60sa.number_of_time_slices,option_feeqs_jt60sa.Vol_poly_degree));

% initialisation of current in PS and VV
switch option_feeqs_jt60sa.mode_init_ps_vv
    case 'Precribed'     
        mode_init_ps_vv =option_feeqs_jt60sa.value_init_ps_vv;
    otherwise
        mode_init_ps_vv = option_feeqs_jt60sa.mode_init_ps_vv;
end

disp(' ');
% call FEESQ.M
dirmem = pwd;
fprintf('Initial directory is %s\n',dirmem);
switch option_feeqs_jt60sa.flux_constant_mode
    case 'manual'
        offset_flux = option_feeqs_jt60sa.flux_offset / 2 / pi;
    case 'auto'
        offset_flux = NaN;
    otherwise
        offset_flux = option_feeqs_jt60sa.flux_constant_mode;
end
% make a local copy of the project
projectpath = fileparts(which('compute_inverse4metis'));
[~,projectname]  = fileparts(projectpath);
% temprary file 
if isdir('/dev/shm')
  templocaldir = fullfile('/','dev','shm',sprintf('%s_%s',fileparts(tempname),getenv('USER')));
else
    templocaldir = tempname;
end
if ~isdir(templocaldir)
    mkdir(templocaldir);
end
unix(sprintf('cp -rpf %s %s',projectpath,templocaldir));
try
    cd(fullfile(templocaldir,projectname));
    
    [output,fbe_lcfs] = compute_inverse4metis_evol(option_feeqs_jt60sa.metisfilename,option_feeqs_jt60sa.number_of_time_slices,option_feeqs_jt60sa.Vol_poly_degree,option_feeqs_jt60sa.first_time,option_feeqs_jt60sa.last_time,...
                                                   option_feeqs_jt60sa.reinforce,option_feeqs_jt60sa.mesh_type,option_feeqs_jt60sa.plotonoff,option_feeqs_jt60sa.movie_onoff,premagnetisation_jt60sa, ...
                                                   option_feeqs_jt60sa.run_name,offset_flux,dirmem,strike_point,option_feeqs_jt60sa.breakdown_flux, ...
                                                   option_feeqs_jt60sa.current_margins,option_feeqs_jt60sa.voltage_margins,mode_init_ps_vv, ...
                                                   option_feeqs_jt60sa.voltage_limits,weight_struct,option_feeqs_jt60sa.init_mode, ...
                                                   option_feeqs_jt60sa.init_voltages,option_feeqs_jt60sa.sigma_PS,option_feeqs_jt60sa.sigma_VV);
    cd(dirmem);
catch
  cd(dirmem);
  error(lasterror);
end
% remove local temporary project directory
unix(sprintf('rm  -rf %s',fullfile(templocaldir,projectname)));
% set ouput
option_mem = option_feeqs_jt60sa;
option_feeqs_jt60sa = output;
option_feeqs_jt60sa.flux_offset = 2 .* pi .* output.psi_off;
option_feeqs_jt60sa.strike_point = option_mem.strike_point;
noms = fieldnames(option_mem);
for k=1:length(noms)
    if ~isfield(option_feeqs_jt60sa,noms{k})
        option_feeqs_jt60sa.(noms{k}) = option_mem.(noms{k});
    end
end
