function option_feeqs_jt60sa = compute_jt60sa_fbe_inverse_fast(option_feeqs_jt60sa)

% test if FEEQS is in the path
set_feeqs_path_jt60sa;


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
    
    valeur.run_name 	    = '';
    type.run_name           = 'string';                    
    borne.run_name          = '';  
    defaut.run_name         = '';                
    info.run_name           = 'run_name is added to name of save files (can be leave empty)';
    
    valeur.plotonoff    = 2;      
    type.plotonoff     = 'integer';
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
   
    valeur.weight    = 5e-4;   
    type.weight      = 'float';
    borne.weight     = [1e-6,1];  
    defaut.weight    = 5e-4;
    info.weight      = 'regularisation term wieght in cost function';

    valeur.mode_profile    = 'on';   
    type.mode_profile     = 'string';
    borne.mode_profile     = {'on','off'};  
    defaut.mode_profile    = 'on';
    info.mode_profile     = 'if = on, use P'' and FF'', otherwise use beta-p and li_3';
    
    valeur.V_Force    = 'off';      
    type.V_Force      = 'integer';
    borne.V_Force     = {'on','off'};  
    defaut.V_Force    = 'off';
    info.V_Force      = 'Compute (on) or not (off) volatages and forces';



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
projectpath = fileparts(which('compute_inverse4metis_fast'));
[~,projectname]  = fileparts(projectpath);
% temprary file 
if isdir('/dev/shm')
  templocaldir = fullfile('/','dev','shm',sprintf('%s_%s',fileparts(tempname),getenv('USER')));
else
  templocaldir = tempname;
end
%  if ~isdir(templocaldir)
%      mkdir(templocaldir);
%  end
%  unix(sprintf('cp -rpf %s %s',projectpath,templocaldir));
try
%  cd(fullfile(templocaldir,projectname));
  
  output =  compute_inverse4metis_fast(option_feeqs_jt60sa.metisfilename,option_feeqs_jt60sa.number_of_time_slices,option_feeqs_jt60sa.first_time,option_feeqs_jt60sa.last_time, ...
                                             option_feeqs_jt60sa.reinforce,option_feeqs_jt60sa.plotonoff, option_feeqs_jt60sa.movie_onoff, premagnetisation_jt60sa, option_feeqs_jt60sa.run_name, ...
                                             option_feeqs_jt60sa.weight,offset_flux,dirmem,option_feeqs_jt60sa.breakdown_flux,option_feeqs_jt60sa.mode_profile,option_feeqs_jt60sa.V_Force);
 
  cd(dirmem);
catch
    keyboard
  cd(dirmem);
  error(lasterror);
end
% remove local temporary project directory
%unix(sprintf('rm  -rf %s',fullfile(templocaldir,projectname)));
% set ouput
option_mem = option_feeqs_jt60sa;
option_feeqs_jt60sa = output;
option_feeqs_jt60sa.flux_offset = 2 .* pi .* output.psi_off;
noms = fieldnames(option_mem);
for k=1:length(noms)
    if ~isfield(option_feeqs_jt60sa,noms{k})
        option_feeqs_jt60sa.(noms{k}) = option_mem.(noms{k});
    end
end
