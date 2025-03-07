function [option_feeqs_jt60sa,fbe_lcfs] = compute_jt60sa_fbe_inverse(option_feeqs_jt60sa)

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
    
    valeur.priority  = 'Coils';      
    type.priority    = 'string';
    borne.priority   = {'Coils','LCFS'};  
    defaut.priority  = 'Coils';
    info.priority    = 'Choice of main objectif: if = Coils, gives priority to the respect of Coils limits and if = LCFS, gives priority to respect plasma shape';

    valeur.mode_init_ps_vv = 'METIS Vloop';     
    type.mode_init_ps_vv   = 'string';
    borne.mode_init_ps_vv  = {'METIS Vloop', 'Precribed','zero','undetermined'};  
    defaut.mode_init_ps_vv = 'METIS Vloop';
    info.mode_init_ps_vv   = 'Methode used to prescrived currents in passive strcutures at each time of computation :\nif = Metis Vloop, used vloop from METIS simulation;\nif = Prescribed, used parameter value_init_ps_vv;\nif = zero, set null current in passive structures\nif = undetermined, no information are provided on current in passive structures';    
    
    valeur.value_init_ps_vv     = 0;   
    type.value_init_ps_vv       = 'float';
    borne.value_init_ps_vv      = [-100,100];  
    defaut.value_init_ps_vv     = 0;
    info.value_init_ps_vv       = 'Value of loop voltage in passive structure elements at first time of computation (V)';
    
    valeur.current_margins     = 0;   
    type.current_margins       = 'float';
    borne.current_margins      = [-0.2,0.2];  
    defaut.current_margins     = 0;
    info.current_margins       = 'margin to be keep to coil current limits (on absolute value, if < 0, allows to go above the limit);\n Pay attention tobe compatible with plasma initialisation value;\n Applied only to EF coils';
    
    valeur.factor_w_pfcur     = 1;   
    type.factor_w_pfcur       = 'float';
    borne.factor_w_pfcur      = [0,10];  
    defaut.factor_w_pfcur     = 1;
    info.factor_w_pfcur       = 'relative weight on penalisation term for current limits in coils';
    
    valeur.weight_force  = 1;      
    type.weight_force    = 'integer';
    borne.weight_force   = {0,1};  
    defaut.weight_force  = 1;
    info.weight_force    = 'if = 1, constraint forces on coil to be in engineering limits';

    valeur.weight_force_diff  = 1;      
    type.weight_force_diff    = 'integer';
    borne.weight_force_diff   = {0,1};  
    defaut.weight_force_diff  = 1;
    info.weight_force_diff    = 'if = 1, constraint difference of force between coils to be in engineering limits';

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
   
    valeur.wmin    = 17;   
    type.wmin      = 'float';
    borne.wmin     = [8,21];  
    defaut.wmin    = 14;
    info.wmin      = 'minimal weight of regularisation term in initialisation loop (power of 10)';

    valeur.wmax    = 17;   
    type.wmax     = 'float';
    borne.wmax     = [14,27];  
    defaut.wmax    = 21;
    info.wmax     = 'maximum weight of regularisation term in initialisation loop (power of 10)';
    
    valeur.ef5_w    = 3;   
    type.ef5_w     = 'float';
    borne.ef5_w     = [0,10];  
    defaut.ef5_w    = 3;
    info.ef5_w     = 'penalisation of current in EF5 coil to prevent it be used as main coil to make X-point (power of 10)';
    
    valeur.delta_i  = 'on';      
    type.delta_i    = 'string';
    borne.delta_i   = {'on','off'};  
    defaut.delta_i  = 'on';
    info.delta_i    = 'Set regularisation on sum(Icoil^2) or on variation of coil currents between two time slices';
    
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

    valeur.init_file  = '';
    type.init_file    = 'string';                    
    borne.init_file   = '';  
    defaut.init_file  = '';                
    info.init_file    = 'if available, for very similar simulation and parameters, name of file containing FEEQS initilalisation data (mesh, initial psi map, ...): this allows to save time;\notherwise initilialisation phase will be rerun. Skiping initialisation phase shotcut make mesh selection, weight optimisation are not working';

    valeur.init_time    = 0;   
    type.init_time      = 'float';
    borne.init_time     = [0,100];  
    defaut.init_time    = 0;
    info.init_time      = 'time selected for equilibrium initialisation (s);\nif = 0, the selection is automatic; must be between first_time and last_time';
    
 
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

% init time
if option_feeqs_jt60sa.init_time > 0
    option_feeqs_jt60sa.init_time = max(option_feeqs_jt60sa.first_time,min(option_feeqs_jt60sa.last_time,option_feeqs_jt60sa.init_time));
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

% initialisation of current in PS and VV
switch option_feeqs_jt60sa.mode_init_ps_vv
    case 'Precribed'     
        mode_init_ps_vv =option_feeqs_jt60sa.value_init_ps_vv;
    otherwise
        mode_init_ps_vv = option_feeqs_jt60sa.mode_init_ps_vv;
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
  [output,fbe_lcfs] =  compute_inverse4metis(option_feeqs_jt60sa.metisfilename,option_feeqs_jt60sa.number_of_time_slices,option_feeqs_jt60sa.first_time,option_feeqs_jt60sa.last_time, ...
                                             option_feeqs_jt60sa.reinforce,option_feeqs_jt60sa.weight_force,option_feeqs_jt60sa.weight_force_diff,option_feeqs_jt60sa.mesh_type, ...
                                             option_feeqs_jt60sa.plotonoff, option_feeqs_jt60sa.movie_onoff, premagnetisation_jt60sa, option_feeqs_jt60sa.run_name, ...
                                             option_feeqs_jt60sa.wmin,option_feeqs_jt60sa.wmax,offset_flux,dirmem,strike_point,option_feeqs_jt60sa.delta_i, ...
                                             option_feeqs_jt60sa.breakdown_flux,option_feeqs_jt60sa.init_file,option_feeqs_jt60sa.ef5_w, ...
                                             option_feeqs_jt60sa.priority,option_feeqs_jt60sa.current_margins,option_feeqs_jt60sa.factor_w_pfcur, ...
                                             option_feeqs_jt60sa.init_time,mode_init_ps_vv,option_feeqs_jt60sa.sigma_PS,option_feeqs_jt60sa.sigma_VV);
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
