%% METIS4IMAS : interface between IMAS data and METIS data
%-------------------------------------------------------------
% fonction Matlab 7; metis4imas.m ->  metis4imas, mapsummary_imas, test_metis_imas, load_metis_imas, summary2option_imas
%					ids2metis_input, ids2metis1t, mapcore_profiles_ids, interp1_imas, griddata_imas, z0rot_imas,
%					mappulse_schedule, 
%
%
% This function is the interface between IMAS/IDS and METIS data structure
%
% syntax : 
%
%   * create xml an xsd files (without output):
%	metis4imas;
%
%   * parameters declaration stucture (with one output):
%	info = metis4imas;
%
%   * returning xml an xsd text (with 2 outputs):
%	[xml,xsd] = metis4imas;
%
%   * computing :
%	[error_flag,output_data] = metis4imas(shot,run,occurrence,method,time,codeparam_filename,interpolation_methode);
%
%   * post processing (with UAl data writing) :
%	error_flag = metis4imas(shot,run,occurrence,metis_data_structure);
%
%   * post processing of zerodevolution (time to time simulation, write only one time slice)
%   error_flag = metis4imas(shot,run,occurrence,zerodevolution_data_structure);
%
%   * post processing (without UAl data writing = return only matlab data structure counterpart of IDSs structures) :
%   [error_flag,output_data] = metis4imas(shot,run,occurrence,<metis_data_structure>,[],optionnal_codeparam_filename);
%
%   * post processing of zerodevolution (time to time simulation, write only one time slice)
%   error_flag = metis4imas(shot,run,occurrence,zerodevolution_data_structure,[],optionnal_codeparam_filename);
%
%   * wrinting a METIS file in database:
%	error_flag = metis4imas(shot,run,occurrence,'metis_file_name',[],optionnal_codeparam_filename);
%
%   * reading a UAL reccord:
%	[error_flag,output_data] = metis4imas(shot,run,occurrence,'read',[],optionnal_codeparam_filename);
%
%   Data are kept in COCOS format and orientations store in IMAS database.
%   You should use "ids_generic_cocos_nodes_transformation_symbolic" to
%   change COCOS and orientations (see https://gitlab.epfl.ch/spc/cocos for
%   details).
%
%
% input :
%
%     shot                    = shot number (integer >0)
%				if is complex, imag(shot) is used as reference shot
%     run                     = run number for this shot (integer >0)
%				if is complex, imag(shot) is used as reference run
%
%     occurrence              = input ids scenario occurrence in Kepler (default = [])
%                               or occurrence can be a structure to use an alternative database :
%					               occurrence.tokamak = new UAL tokamak database name
%					               occurrence.user    = new UAL user database name
%                                  occurrence.dataversion = selected a differente version of data (not the last one)
%		                           occurrence.occurrence = input ids scenario occurrence in Kepler (default = [])
%		                           occurrence.backend = selection of IMAS backend (default = [])
% 
%     method                  = metis command for Kepler actor or METIS file name (default = '') or metis data structure
%     time                    = time scalar or vector for which MESTIS must be run (scalar). 
%                               If time is scalar, METIS is run in evolution mode.
%                               If time is vector, METIS is run for whole shot simulation.
%                               If empty, reads all time slice of the IDS pulse_schedule and uses pulse_schedule IDS vector time as input.
%
%     codeparam_filename      = is name of the file containing the METIS parameters (XML),
%                               or XML string coding for codeparam,
%                               or matlab structure coding for METIS (all or some parameters,codeparam_filename.key = value).
%
%                               this input is optionnal. Defaults METIS parameters values are used for missing parameters. 
%                               this input is not usedin test mode.
%                               if method is 'read' ,'writing' or 'post processing' , the codeparam input allows to sellected the IDSs
%                               that are read or write (only the IDS pulse_schedule is mandatory)
%
%     interpolation_methode   = interpolation method used for data extraction (see UAL user manual).
%                               1 = closest sample, 2 = previous sample & 3 = linear interpolation.
%                               optionnal input.
%
% output :
%
%      err_flag		      = error code, must be = 0 if no error (as for unix program)
%      output_data        = matlab structure containing the couterpart of IDSs structures in matlab compatible format.
%                           Data are re-read in the UAL juste after wrinting.
%
%  method list: 
% 
%       'test'                 = test the function (call METIS with test data) and create IDSs.
%       'auto_test'            = complete test of metis4imas (all methods)
%       'full'                 = complete shot computation (all time slices given in pulse_schedule), full METIS computation mode
%       'fast'                 = complete shot computation (selected time slices given in pulse_schedule), fast METIS computation mode
%       'init'                 = initialisation of the evolution mode of METIS, pulse_schedule must contain the fisrt time slice for references.
%       'one_time'             = compute plasma evolution for one time step (used evolution mode of METIS), 
%                                pulse_schedule must contain the next time slice for references.
%        'read'                = read UAL data and return in matlab structure
%
%       <filename>             = load a METIS simulation and create associated IDS matlab data structures.
%                                (this is also the method for make a restart)
%
%  if method is empty, then the choice is automatic:
%      - if pulse_schedule is empty, method = 'test'
%      - if pulse_schedule contains only one time slice, a the fisrt call, method = 'init', 
%        for next calls, method = 'one_time' 
%      - if pulse_schedule contains more of one time slice, method = 'full'
%
%  stucture zerodevolution_data_structure:
%  Concatenate zerodevolution outputs in a stucture
%    [zs,profil,z0dstruct] = zerodevolution( ...
%    zerodevolution_data_structure.zs = zs;
%    zerodevolution_data_structure.profil = profil;
%    zerodevolution_data_structure.z0dstruct = z0dstruct;
%    zerodevolution_data_structure.z0dstruct = z0dstruct;
%
%  Information about METIS is avalaible in METIS technical documentation that is included in the METIS distribution or with the button 
%  <HELP> of METIS GUI (PDF reader must be available).
%
%  Script to call METIS from a workflow :
%
%	addpath <path to cronos project>
%	zineb_path
%
%       % mapping of workflow matlab actor input variables
%
%       shot        = workflow_input_num_shot_variable;
%       run         = workflow_input_num_shot_run;
%
%       % optionnal :
%	    % occurrence = workflow_input_metis_ids_occurrence;
%       occurrence  = [];
%
%       % optionnal :
%	    % method      = workflow_input_metis_method;
%       method      = '';
%
%	    % optionnal :
%       % time        = workflow_vetor_or_scalar_input_time;
%	    time        = [];
%
%	    % optionnal
%	    % codeparam_filename = workflow_file_name_or_xml_string;
%	    codeparam_filename = '';
%
%	    % optionnal : 
%       % interpolation_method = 2;
%
%       % error_flag  = metis4imas(shot,run,occurrence,method,time,codeparam_filename,interpolation_method);
%        error_flag  = metis4imas(shot,run,occurrence,method,time);
%
%
%	    % return to workflow
%	    workflow_return_value = error_flag;
%
%
%  access to metis internal data :
%
%	z0dstruct = getappdata(0,'IMAS_Z0DSTRUCT');
%
%  simple test of the function :
%
%       [error_lfag,output] = metis4imas(1,1);
% 
%  full test of the function :
%
%  	metis4imas(1,1,'','auto_test')  
%
% fonction ecrite par J-F Artaud
% version Git (last change 2024/06/18)
%-----------------------------------------------------------------------
%
function [error_flag,xsd] = metis4imas(shot,run,occurrence,methode,time,codeparam_filename,interpolation_methode)

% test availability of tools for cocos
if isempty(which('ids_check_cocos'))
    init_cocos_test_in_metis([],true);
end


%% MANAGE NUMBER OF INPUT VARIABLEScore_sources
if (nargin == 0) && ((nargout == 0)  || (nargout == 2))
  help metis4imas
  %% Create codeparam files
  module_xsd_make('metis4imas',1);
  %% Ouput codeparam
  [xsd,error_flag]=module_xsd_make('metis4imas',1);
  return
end
if nargin <= 1
  error_flag = zerod;
  %% Variables specific to metis4imas
  names_before = fieldnames(error_flag.valeur);
  

  %% IDS initialization at first call
  error_flag.valeur.init_output_ids  = 0;
  error_flag.type.init_output_ids    = 'integer';
  error_flag.borne.init_output_ids   = {0,1}; 
  error_flag.defaut.init_output_ids  = 0;
  error_flag.info.init_output_ids    = 'UAL control:\nif = 0, continue to write at the end of existing record;\nif=1 on init call, initialise output IDS (reset and write the 1st time slice)';
  error_flag.section.init_output_ids = 'UAL';
        
  %% Restart file name
  error_flag.valeur.restart  = '';
  error_flag.type.restart    = 'string';
  error_flag.borne.restart   = ''; 
  error_flag.defaut.restart  = '';
  error_flag.info.restart    = 'UAL control:name of the restart file;\nif empty, no restart file is saved; otherwise, after each call, the restart file is saved';
  error_flag.section.restart = 'UAL';

  %% ---------------------------
  %% DECLARATION OF OUTPUT IDSS
  %% ---------------------------
  error_flag.valeur.summary   = 1;
  error_flag.type.summary     = 'integer';
  error_flag.borne.summary    = {0,1}; 
  error_flag.defaut.summary   = 1;
  error_flag.info.summary     = 'IDS selection:\nif = 0, do not write summary IDS;\nif=1, write summary IDS';
  error_flag.section.summary  = 'UAL';

  error_flag.valeur.core_profiles   = 1;
  error_flag.type.core_profiles     = 'integer';
  error_flag.borne.core_profiles    = {0,1}; 
  error_flag.defaut.core_profiles   = 1;
  error_flag.info.core_profiles     = 'IDS selection:\nif = 0, do not write core_profiles IDS;\nif=1, write core_profiles IDS';
  error_flag.section.core_profiles  = 'UAL';

  error_flag.valeur.core_transport  = 1;
  error_flag.type.core_transport    = 'integer';
  error_flag.borne.core_transport   = {0,1}; 
  error_flag.defaut.core_transport  = 1;
  error_flag.info.core_transport    = 'IDS selection:\nif = 0, do not write core_transport IDS;\nif=1, write core_transport IDS';
  error_flag.section.core_transport = 'UAL';

  error_flag.valeur.core_sources  = 1;
  error_flag.type.core_sources    = 'integer';
  error_flag.borne.core_sources   = {0,1}; 
  error_flag.defaut.core_sources  = 1;
  error_flag.info.core_sources    = 'IDS selection:\nif = 0, do not write core_sources IDS;\nif=1, write core_sources IDS';
  error_flag.section.core_sources = 'UAL';

  error_flag.valeur.radiation  = 1;
  error_flag.type.radiation    = 'integer';
  error_flag.borne.radiation   = {0,1}; 
  error_flag.defaut.radiation  = 1;
  error_flag.info.radiation    = 'IDS selection:\nif = 0, do not write radiation IDS;\nif=1, write core_sources IDS';
  error_flag.section.radiation = 'UAL';

  error_flag.valeur.edge  = 1;
  error_flag.type.edge    = 'integer';
  error_flag.borne.edge   = {0,1}; 
  error_flag.defaut.edge  = 1;
  error_flag.info.edge    = 'IDS selection:\nif = 0, do not write edge IDSs;\nif=1, write edge IDSs (edge_profiles and edge_transport)';
  error_flag.section.edge = 'UAL';

  error_flag.valeur.numerics  = 1;
  error_flag.type.numerics    = 'integer';
  error_flag.borne.numerics   = {0,1}; 
  error_flag.defaut.numerics  = 1;
  error_flag.info.numerics    = 'IDS selection:\nif = 0, do not write numerics IDSs;\nif=1, write numerics IDSs (transport_solver_numerics)';
  error_flag.section.numerics = 'UAL';

  error_flag.valeur.equilibrium = 1;
  error_flag.type.equilibrium	= 'integer';
  error_flag.borne.equilibrium  = {0,1}; 
  error_flag.defaut.equilibrium	= 1;
  error_flag.info.equilibrium	= 'IDS selection:\nif = 0, do not write equilibrium IDS;\nif=1, write equilibrium IDS';
  error_flag.section.equilibrium= 'UAL';
  
%    error_flag.valeur.grid_equi 	= 0;
%    error_flag.type.grid_equi	= 'integer';
%    error_flag.borne.grid_equi  	= {0,1}; 
%    error_flag.defaut.grid_equi	= 0;
%    error_flag.info.grid_equi  	= 'mesh grid used for the equilibrum 2 D profile :\n  if = 0, 2D grid is in (psi,theta);\n  if = 1, 2D grid is a Delaunay grid';
%    error_flag.section.grid_equi  = 'UAL';  

  error_flag.valeur.equi_extrap  = 0;
  error_flag.type.equi_extrap    = 'integer';
  error_flag.borne.equi_extrap   = {0,1,2,3}; 
  error_flag.defaut.equi_extrap  = 0;
  error_flag.info.equi_extrap    = '2D equilibrium - method for extrapolation of Psi oustside the LCFS:\nif = 0, interpolation of Psi using a polynomial G-S solution on each LCFS point;\nif = 1, hybrid method: Analitical solution of GS, or simple extrapolation if non converged;\nif = 2, recompute equilibrium with fixed boundary equilibrium solver of FEEQS.M code (if available, FEEEQS.M should have been installed separately);\nif = 3, as option 2 but used polynomial solution of G-S constained with flux and magnetic field on each point of LCFS for extrapolation ouside the LCFS instead of a simple interpolation';
  error_flag.section.equi_extrap = 'UAL';
  
  error_flag.valeur.Convex_LCFS  = 1;
  error_flag.type.Convex_LCFS    = 'integer';
  error_flag.borne.Convex_LCFS   = {0,1}; 
  error_flag.defaut.Convex_LCFS  = 1;
  error_flag.info.Convex_LCFS    = '2D equilibrium: force LCFS used for 2D extrapolation to be convex:\nif = 0, keep LCFS as it is provided;\nif = 1, force LCFS to be convex';
  error_flag.section.Convex_LCFS = 'UAL';

  error_flag.valeur.fixed_grid  = 1;
  error_flag.type.fixed_grid    = 'integer';
  error_flag.borne.fixed_grid   = {0,1}; 
  error_flag.defaut.fixed_grid  = 1;
  error_flag.info.fixed_grid    = '2D equilibrium, for rectangular grid:\nif = 0, uses floating grid following plasma displacement;\nif = 1 uses same grid for all time slices';
  error_flag.section.fixed_grid = 'UAL';
  
  error_flag.valeur.nb_points_pol  = 65;
  error_flag.type.nb_points_pol    = 'integer';
  error_flag.borne.nb_points_pol   = [35,255]; 
  error_flag.defaut.nb_points_pol  = 65;
  error_flag.info.nb_points_pol    = '2D equilibrium: number of points in poloidal direction for inverse (rho,theta) grid of equilibrium';
  error_flag.section.nb_points_pol = 'UAL';
  
  error_flag.valeur.nb_points_radial  = 51;
  error_flag.type.nb_points_radial    = 'integer';
  error_flag.borne.nb_points_radial   = [33,301]; 
  error_flag.defaut.nb_points_radial  = 51;
  error_flag.info.nb_points_radial    = '2D equilibrium: number of points in radial direction for rectangular (R,Z) grid of equilibrium';
  error_flag.section.nb_points_radial = 'UAL';
  
  %% ---------------------------------
  %% OCCURRENCE STRING OF OUTPUT IDSS
  %% ---------------------------------
  error_flag.valeur.pulse_schedule_occurrence  = '';
  error_flag.type.pulse_schedule_occurrence    = 'string';
  error_flag.borne.pulse_schedule_occurrence   = ''; 
  error_flag.defaut.pulse_schedule_occurrence  = '';
  error_flag.info.pulse_schedule_occurrence    = 'UAL control: output ids occurrence for pulse_schedule IDS;\nif empty, use default occurence (0)';
  error_flag.section.pulse_schedule_occurrence = 'Occurrence UAL';

  error_flag.valeur.summary_occurrence  = '';
  error_flag.type.summary_occurrence    = 'string';
  error_flag.borne.summary_occurrence   = ''; 
  error_flag.defaut.summary_occurrence  = '';
  error_flag.info.summary_occurrence    = 'UAL control: output ids occurrence for summary IDS;\nif empty, use default occurence (0)';
  error_flag.section.summary_occurrence = 'Occurrence UAL';

  error_flag.valeur.core_profiles_occurrence  = '';
  error_flag.type.core_profiles_occurrence    = 'string';
  error_flag.borne.core_profiles_occurrence   = ''; 
  error_flag.defaut.core_profiles_occurrence  = '';
  error_flag.info.core_profiles_occurrence    = 'UAL control: output ids occurrence for core_profiles IDS;\nif empty, use default occurence (0)';
  error_flag.section.core_profiles_occurrence = 'Occurrence UAL';

  error_flag.valeur.core_transport_occurrence  = '';
  error_flag.type.core_transport_occurrence    = 'string';
  error_flag.borne.core_transport_occurrence   = ''; 
  error_flag.defaut.core_transport_occurrence  = '';
  error_flag.info.core_transport_occurrence    = 'UAL control: output ids occurrence for core_transport IDS;\nif empty, use default occurence (0)';
  error_flag.section.core_transport_occurrence = 'Occurrence UAL';
  
  error_flag.valeur.core_sources_occurrence  = '';
  error_flag.type.core_sources_occurrence    = 'string';
  error_flag.borne.core_sources_occurrence   = ''; 
  error_flag.defaut.core_sources_occurrence  = '';
  error_flag.info.core_sources_occurrence    = 'UAL control: output ids occurrence for core_sources IDS;\nif empty, use default occurence (0)';
  error_flag.section.core_sources_occurrence = 'Occurrence UAL';

  error_flag.valeur.radiation_occurrence  = '';
  error_flag.type.radiation_occurrence    = 'string';
  error_flag.borne.radiation_occurrence   = ''; 
  error_flag.defaut.radiation_occurrence  = '';
  error_flag.info.radiation_occurrence    = 'UAL control: output ids occurrence for radiation IDS;\nif empty, use default occurence (0)';
  error_flag.section.radiation_occurrence = 'Occurrence UAL';

  error_flag.valeur.edge_occurrence  = '';
  error_flag.type.edge_occurrence    = 'string';
  error_flag.borne.edge_occurrence   = ''; 
  error_flag.defaut.edge_occurrence  = '';
  error_flag.info.edge_occurrence    = 'UAL control: output ids occurrence for edge IDSs (egde_profiles and edge_transport);\nif empty, use default occurence (0)';
  error_flag.section.edge_occurrence = 'Occurrence UAL';

  error_flag.valeur.numerics_occurrence  = '';
  error_flag.type.numerics_occurrence    = 'string';
  error_flag.borne.numerics_occurrence   = ''; 
  error_flag.defaut.numerics_occurrence  = '';
  error_flag.info.numerics_occurrence    = 'UAL control: output ids occurrence for transport_solver_numerics IDSs;\nif empty, use default occurence (0)';
  error_flag.section.numerics_occurrence = 'Occurrence UAL';

  error_flag.valeur.equilibrium_occurrence  = '';
  error_flag.type.equilibrium_occurrence    = 'string';
  error_flag.borne.equilibrium_occurrence   = ''; 
  error_flag.defaut.equilibrium_occurrence  = '';
  error_flag.info.equilibrium_occurrence    = 'UAL control: output ids occurrence for equilibrium IDS;\nif empty, use default occurence (0)';
  error_flag.section.equilibrium_occurrence = 'Occurrence UAL';

  error_flag.valeur.COCOS  = 11;
  error_flag.type.COCOS    = 'integer';
  error_flag.borne.COCOS   = [1,18]; 
  error_flag.defaut.COCOS  = 11;
  error_flag.info.COCOS    = 'UAL control: choice for the output COCOS';
  error_flag.section.COCOS = 'UAL';
  
  error_flag.valeur.COCOS_method   = 'Sauter';
  error_flag.type.COCOS_method     = 'string';
  error_flag.borne.COCOS_method    = {'Sauter','Native'}; 
  error_flag.defaut.COCOS_method   = 'Sauter';
  error_flag.info.COCOS_method     = 'Method use to make the COCOS mapping:\nNative = native METIS method.\nSauter = use O.Sauter tools (see https://gitlab.epfl.ch/spc/cocos)';
  error_flag.section.COCOS_method  = 'UAL';
  
  error_flag.valeur.COCOS_verbose   = 1;
  error_flag.type.COCOS_verbose     = 'integer';
  error_flag.borne.COCOS_verbose    = {0,1,2,3}; 
  error_flag.defaut.COCOS_verbose   = 1;
  error_flag.info.COCOS_verbose     = 'if COCOS_method = Sauter, level of verbosity of COCOS transformation (0 = no text output).';
  error_flag.section.COCOS_verbose  = 'UAL';
  
  error_flag.valeur.COCOS_check   = 'on';
  error_flag.type.COCOS_check     = 'string';
  error_flag.borne.COCOS_check    = {'off','on'}; 
  error_flag.defaut.COCOS_check   = 'on';
  error_flag.info.COCOS_check     = 'if COCOS_method = Sauter, siwtch on the control of COCOS in obtained IDS for core_profiles and equilbrium.';
  error_flag.section.COCOS_check  = 'UAL';
  
  noms = fieldnames(error_flag.valeur);
  for k=1:length(noms)
      if isempty(strmatch(noms{k},names_before))
	  error_flag.mode.(noms{k}) = 'advanced';
      end
  end

  return
end

if nargin < 2
  error_flag = - 9999;
  return
end

%% garbage collection to prevent java heap overflow
% not useful with current IMAS Matlab API
% try
%   java.lang.System.gc()
% catch
%   try
%     java.lang.Runtime.getRuntime().gc
%   end
% end

%% HLI mex configuration
if exist('imas_set_mex_params')
    imas_set_mex_params('get_int_as_double',true);
    imas_set_mex_params('put_int_from_double',true);
    imas_set_mex_params('get_empty_as_nan',false);
    imas_set_mex_params('use_cell_array_for_array_of_structures',true);
    imas_set_mex_params('error_on_missing_field',true);
    %imas_set_mex_params('convert_whole_ids',false);
    imas_set_mex_params('verbosity',0);        
end

% no more handle with environment variables
%% CHANGE DATABASE
%setappdata(0,'UAL_TOKAMAK','');
%setappdata(0,'UAL_USER','');
%setappdata(0,'UAL_DATAVERSION','');

%% THE OCCURRENCE PARAMETER CAN BE A STRUCTURE
if nargin < 3
  occurrence = '';
elseif ischar(occurrence)
  occurrence = fix(str2num(occurrence));
  if occurrence == 0
    occurrence = '';
  else
    occurrence = num2str(occurrence);		
  end
elseif isnumeric(occurrence)
  occurrence = fix(occurrence);
  if occurrence == 0
    occurrence = '';
  else
    occurrence = num2str(occurrence);		
  end
elseif isstruct(occurrence)
  setappdata(0,'UAL_TOKAMAK',strtrim(occurrence.tokamak));
  setappdata(0,'UAL_USER',strtrim(occurrence.user));
  setappdata(0,'UAL_DATAVERSION',strtrim(occurrence.dataversion));
  if isfield(occurrence,'backend')
      setappdata(0,'UAL_BACKEND',occurrence.backend);
  else
      setappdata(0,'UAL_BACKEND','');      
  end
  if isfield(occurrence,'occurrence')
    occurrence = occurrence.occurrence;
  else
    occurrence = '';
  end
  if ischar(occurrence)
    occurrence = fix(str2num(deblank(occurrence)));
    if occurrence == 0
      occurrence = '';
    else
      occurrence = num2str(occurrence);		
    end
  elseif isnumeric(occurrence)
    occurrence = fix(occurrence);
    if occurrence == 0
      occurrence = '';
    else
      occurrence = num2str(occurrence);		
    end
  end
end
if nargin < 4
  methode = '';
end
if nargin < 5
  time =[];
end
if nargin < 6
  codeparam_filename ='';
end
if nargin < 7
  interpolation_methode = 2;
  % 1 = closest sample, 2 = previous sample & 3 = linear interpolation
end

%% DATA DECLARATION
option 	= [];
cons1t 	= [];
geo1t 	= [];
sepa1t 	= [];
optin   = metis4imas(1);
z0dstruct.z0dinput.option = optin.valeur;

%% EMPTY IDS TO START WITH
pulse_schedule         = [];
dataset_description    = [];
dataset_fair           = [];
summary                = [];
core_profiles          = [];
core_transport         = [];
core_sources           = [];
%  core_sources_lhcd      = [];
%  core_sources_eccd      = [];
%  core_sources_icrh      = [];
%  core_sources_nbicd     = [];
radiation = [];
%  core_sources_cyclotron = [];
%  core_sources_neutral   = [];
%  core_sources_fusion    = [];
%  core_sources_ohm_boot  = [];
%  core_sources_full      = [];
equilibrium            = [];

%% POINTER TOWARD INPUT IDS
pulse_schedule_ids_in = [];
dataset_description_in = [];	
dataset_fair_in = [];	

%% ERROR HANDLING: NO ERROR A PRIORI
error_flag = 0;

%% INITIALIZATION FLAG FOR OUTPUT IDSS
init_output_ids = 0;

%% 
%% LOGFILE
%%
switch get(0,'diary')
case 'on'
	diary off
	diary_file = get(0,'DiaryFile');
	diary on
	diary_already_on = 1;
otherwise
	diary off
	diary_file = tempname;
	diary(diary_file);
	diary_already_on = 0;
end
disp('Starting metis4imas')
%% MANAGE INTERNAL MODELS
%% WRITE_MODE : IF = 1, UAL WRITE ACTIVE
write_mode = 1;
%% READ_INPUT : IF = 1, UAL READ OF PUSLE_SCHEDULE IDS ACTIVE
read_input = 0;
%% EXIT AFTER READING METIS INPUT DATA
exit_after_read_input = 0;
%% READ_OUTPUT: IF = 1, UAL READ OF IDS OUTPUT ACTIVE
read_output = 0;
post = [];
if isstruct(methode)
  read_input = 0;
  write_mode = 0;
  post = methode;
  if isfield(post,'z0dstruct')
        methode = 'post_one_time';
        post.z0dinput = post.z0dstruct.z0dinput;
  else
        methode = 'post';
  end
  if nargout < 2
    write_mode = 1;
  end
	
  %% COMPLETE OPTIONS IF NECESSARY
  %% OPTION
  model = metis4imas(1);
  noms = fieldnames(model.valeur);
  for k = 1:length(noms)
    if ~isfield(post.z0dinput.option,noms{k})
      post.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
    end
  end
  post.z0dinput.option = secure_option_imas(post.z0dinput.option);
  % with updated options
  post.z0dstruct.z0dinput = post.z0dinput;
else
  switch methode
   case 'read'
    write_mode  = 0;
    read_input  = 0;	
    read_output = 1;
   case 'metis_from_ual'
    write_mode  = 0;
    read_input  = 1;	
    read_output = 0;
    exit_after_read_input = 1;    
   case {'full','fast','init','one_time'}
    read_input = 1;	
  end
  if nargout == 2 
    read_output = 1;
  end
end

%% DATABASE IDENTIFICATION
if (write_mode == 1) || (read_input == 1) || (read_output == 1)
  [treename,user,ver] = which_MDSdatabase_metis_imas;
  if  isempty(ver)
    disp('No database available')
    err_flag = -1001;
    return
  end
end

% handling backend
id_backend = get_imas_backend_id;


%% OPEN DATABASE IF NEEDED
if read_input == 1
    try
        disp('Opening tree');
        if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
            if ~isempty(id_backend)
                expIdx = imas_open_env_backend(real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'),id_backend);
            else
                expIdx = imas_open_env('ids',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'));
            end
        else
            expIdx = imas_open('ids',real(shot),real(run));
        end
        expIdx = mem_cache_imas(expIdx);
    catch
        %% NO INPUT DATA
        warning('no initial IDS available');
        expIdx = NaN;
    end
    mode_full = 0;
    switch methode
        case {'fast','full'}
            mode_full = 1;
    end
    if isfinite(expIdx) & (expIdx >=0)
        %% OCCURENCE MANAGEMENT
        if isempty(occurrence)
            nids = 'pulse_schedule';
        else
            nids = sprintf('pulse_schedule/%s',occurrence);
        end
        %% GET PULSE_SCHEDULE SCENARIO
        disp('Reading pulse_schedule IDS');
        if isempty(time)  | (mode_full == 1)
            try
                pulse_schedule_ids_in =ids_get(expIdx,nids);
            catch
                pulse_schedule_ids_in = [];
            end
        else
            try
                pulse_schedule_ids_in = ids_get_slice(expIdx,nids,time,interpolation_methode);
            catch
                pulse_schedule_ids_in = [];
            end
        end
        %% GET DATASET_DESCRIPTION
        %% OCCURENCE MANAGEMENT
        if isempty(occurrence)
            nids = 'dataset_description';
        else
            nids = sprintf('dataset_description/%s',occurrence);
        end
        disp('Reading dataset_descriptionp IDS');
        if isempty(time)  | (mode_full == 1)
            try
                dataset_description_in =ids_get(expIdx,nids);
            catch
                dataset_description_in = [];
            end
        else
            try
                dataset_description_in = ids_get_slice(expIdx,nids,time,interpolation_methode);
            catch
                dataset_description_in = [];
            end
        end
        %% GET DATASET_FAIR
        %% OCCURENCE MANAGEMENT
        if isempty(occurrence)
            nids = 'dataset_fair';
        else
            nids = sprintf('dataset_fair/%s',occurrence);
        end
        disp('Reading dataset_fair IDS');
        if isempty(time)  | (mode_full == 1)
            try
                dataset_fair_in =ids_get(expIdx,nids);
            catch
                dataset_fair_in = [];
            end
        else
            try
                dataset_fair_in = ids_get_slice(expIdx,nids,time,interpolation_methode);
            catch
                dataset_fair_in = [];
            end
        end
        
        
        %% CLOSE THE CURRENTLY OPEN DATABASE
        imas_close(expIdx);
    else
        pulse_schedule_ids_in = [];
        dataset_description_in = [];
        dataset_fair_in = [];
    end
    
end

%% MANAGE INPUT DATA
if ~isempty(pulse_schedule_ids_in)
    %% ACCORDING TO NUMBER OF TIME STEPS
    
    if length(pulse_schedule_ids_in) == 1
        if isempty(methode)
            methode = 'one_time';
        end
    end
    switch methode
        case {'init','one_time'}
            
            %% TRANSFORM THE IDS INTO METIS STRUCTURE
            disp('Parsing input data from IDS pulse_schedule to METIS for one time');
            %% EVOLUTION MODE
            [option,time_ids,cons1t,geo1t,sepa1t] = ids2metis1t(pulse_schedule_ids_in);
            if isempty(time)
                time = time_ids;
            elseif isempty(time_ids)
                cons1t.temps = time;
            end
            option = secure_option_imas(option);
            % recompute flux reference from vloop one and remove vloop from data structure
            switch methode
                case 'init'
                    cons1t.flux = 0;
                    cons1t.vloop = NaN; % this is no a real field of cons1t
                otherwise
                    if isappdata(0,'IMAS_Z0DSTRUCT')
                        temp_strcut0d = getappdata(0,'IMAS_Z0DSTRUCT');
                        if time <= temp_strcut0d.z0dinput.cons.temps(end)
                            cons1t.flux   = interp1(temp_strcut0d.z0dinput.cons.temps,temp_strcut0d.z0dinput.cons.flux,time,'linear','extrap');
                        else
                            cons1t.flux   = temp_strcut0d.z0dinput.cons.flux(end) - ...
                                cons1t.vloop .* (time - temp_strcut0d.z0dinput.cons.temps(end)) ./ 2 ./ pi ;
                        end
                    else
                        error('unknown internal state for METIS4IMAS');
                    end
                    cons1t.vloop = NaN; % this is no a real field of cons1t
            end
        otherwise
            
            %% TRANSFORM IDS INTO CRONOS DATA
            disp('Parsing input data from IDS pulse_schedule to METIS');
            % mode simulation complete
            z0dinput = ids2metis_input(pulse_schedule_ids_in);
            z0dinput.option = secure_option_imas(z0dinput.option);
            if ~isempty(dataset_description_in)
                z0dinput.dataset_description = dataset_description_in;
            end
            z0dstruct.z0dinput = z0dinput;
            if isempty(methode)
                methode = 'full';
            end
    end
elseif isempty(methode)
    methode = 'test';
elseif read_input == 1
    switch methode
        case {'full','fast','init','one_time','read'}
            disp('No data available in database of UAL memory')
            err_flag = -3007;
            return
    end
end

% lecture des codeparam
if  ~isempty(codeparam_filename) && isstruct(codeparam_filename)
  if isempty(option)
      info = metis4imas(1);
      option = info.valeur;
  end
  % codeparam_filename est une structure matlab
  info = codeparam_filename;
  noms = fieldnames(info);
  for k=1:length(noms)
    if isfield(option,noms{k})
      option.(noms{k}) = info.(noms{k});
    end
  end
  option = secure_option_imas(option);
  z0dstruct.z0dinput.option = option;
  z0dinput.option = option;
  fprintf('METIS4IMAS using parameters from input :\n')
  option
  
elseif ~isempty(codeparam_filename)
  % codeparam_filename designe un fichier xml
  % lecture du fichier
  fid = fopen(codeparam_filename,'r');
  if fid > 0
    codeparam = char(fread(fid,Inf,'char')');
    fclose(fid);
    % add codeparam to options
    info = metis4imas(1);
    option = info.valeur;
    if ~isempty(codeparam)
      info  = xml_read(codeparam_filename);
      noms = fieldnames(info);
      for k=1:length(noms)
	if isfield(option,noms{k})
	  option.(noms{k}) = info.(noms{k});
	end
      end
    end 
    option = secure_option_imas(option);
    z0dstruct.z0dinput.option = option;
    z0dinput.option = option;
    fprintf('METIS4IMAS  using parameters from %s :\n',codeparam_filename)
    option
  elseif ~isempty(codeparam_filename)    
    
    % codeparam_filename est une chaine xml
    
    tnp = tempname;
    fid = fopen(tnp,'w');
    fprintf(fid,'%s\n',codeparam_filename);
    fclose(fid);
    % add codeparam to options
    info = metis4imas(1);
    option = info.valeur;
    info  = xml_read(tnp);
    noms = fieldnames(info);
    for k=1:length(noms)
        if isfield(option,noms{k})
            option.(noms{k}) = info.(noms{k});
        end
    end
    option = secure_option_imas(option);
    z0dstruct.z0dinput.option = option;
    z0dinput.option = option;
    fprintf('METIS4IMAS using parameters from input :\n')
    option
    delete(tnp);
  else
      disp(sprintf('Unable to read codeparam file %s',codeparam_filename));
      err_flag = -7001;
      return
  end
  
end
if exit_after_read_input == 1
  % retourne les donnees pour metis init (zerod_init)
  error_flag = 0;
  xsd = z0dstruct.z0dinput;
  return
end    

% securite sur les option

% cas auto_test
fprintf('METIS4IMAS called with method %s\n',methode);
switch methode
 case 'auto_test'
  error_flag = 1;
  try
    metis4imasautotest(real(shot),real(run));
    error_flag = 0;
    if ~diary_already_on
    	diary off
    	delete(diary_file);
    end
  catch
    if ~diary_already_on
    	diary off
    end
    lasterror 
    keyboard
  end
  error_flag = 0;
  return
end

% restart file
save_restart_file = 0;

% error handling 
% par defaut si le programme s'arrete
error_flag = 1;
try
    % selon le methode
    switch methode
        case 'test'
            % cette methode appel metis en mode test pour genere des donnees de test;
            disp('Test launched')
            if exist('option','var') && isstruct(option) && ~isempty(option)
                z0dstruct  = test_metis_imas(option.signe,option.orientation,option.COCOS,option.COCOS_method);
            else
                z0dstruct  = test_metis_imas;              
            end
            data_zerod = z0dstruct.zerod;
            profil0d   = z0dstruct.profil0d;
            z0dinput   = z0dstruct.z0dinput;
            % compatibilite ascendante
            z0dstruct.zs = z0dstruct.zerod;
            
        case 'full'
            % cette methode calcul la simulation complete en mode fast : le IDS d'entree doit contenir au moins 3 temps
            disp('Full computation')
            [z0dstruct.zerod,void,z0dstruct.profil0d] = zerod(z0dstruct.z0dinput.option,z0dstruct.z0dinput.cons, ...
                z0dstruct.z0dinput.geo,z0dstruct.z0dinput.exp0d);
            data_zerod = z0dstruct.zerod;
            profil0d   = z0dstruct.profil0d;
            z0dinput   = z0dstruct.z0dinput;
            % compatibilite ascendante
            z0dstruct.zs = z0dstruct.zerod;
            
        case 'fast'
            % cette methode calcul la simulation complete en mode fast : le IDS d'entree doit contenir au moins 5 temps
            disp('Fast computation')
            [z0dstruct.zerod,void,z0dstruct.profil0d] = zerodfast(z0dstruct.z0dinput.option,z0dstruct.z0dinput.cons,...
                z0dstruct.z0dinput.geo,z0dstruct.z0dinput.exp0d);
            data_zerod = z0dstruct.zerod;
            profil0d   = z0dstruct.profil0d;
            z0dinput   = z0dstruct.z0dinput;
            % compatibilite ascendante
            z0dstruct.zs = z0dstruct.zerod;
            
        case 'init'
            % cette methode precede one-time, elle doit etre appelee juste avant le calcul pour le premier intervalle de temps
            disp('Init METIS for evolution computation')
            if isempty(time)
                disp('The time must be given with method init')
                error_flag = 4001;
                return
            end
            if ~isempty(sepa1t)
                [data_zerod,profil0d,z0dstruct] = zerodevolution([],option,time,cons1t,geo1t,[],sepa1t);
            else
                [data_zerod,profil0d,z0dstruct] = zerodevolution([],option,time,cons1t,geo1t,[]);
            end
            % UAL option is not always preserve (security rule on parameters)
            noms = fieldnames(option);
            for kl=1:length(noms)
                if ~isfield(z0dstruct.z0dinput.option,noms{kl})
                    z0dstruct.z0dinput.option.(noms{kl}) = option.(noms{kl});
                end
            end
            % compatibilite ascendante
            z0dstruct.zerod = z0dstruct.zs;
            
            % add possible missign field in option
            model = metis4imas(1);
            noms = fieldnames(model.valeur);
            for k = 1:length(noms)
                if ~isfield(z0dstruct.z0dinput.option,noms{k})
                    z0dstruct.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
                end
            end
            z0dstruct.z0dinput.option = secure_option_imas(z0dstruct.z0dinput.option);

            
            % initialisation des IDS de sortie
            if option.init_output_ids
                init_output_ids = 1;
            end
        case 'one_time'
            if isempty(time)
                disp('The time must be given with method one_time')
                error_flag = 4002;
                return
            end
            if isappdata(0,'IMAS_Z0DSTRUCT')
                disp('One time step METIS computation (using internal memorized data)')
                z0dstruct = getappdata(0,'IMAS_Z0DSTRUCT');
                % restart ?
                if ~isempty(option.restart)
                    save_restart_file = 1;
                end
                
            else
                % equivalent a init
                disp('Init METIS for evolution computation')
                z0dstruct = [];
                % initialisation des IDS de sortie
                if option.init_output_ids
                    init_output_ids = 1;
                end
            end
            % cette methode calcul l'evolution du plasma pour un intervalle de temps
            if ~isempty(sepa1t)
                [data_zerod,profil0d,z0dstruct] = zerodevolution(z0dstruct,option,time,cons1t,geo1t,[],sepa1t);
            else
                [data_zerod,profil0d,z0dstruct] = zerodevolution(z0dstruct,option,time,cons1t,geo1t,[]);
            end
            % compatibilite ascendante
            z0dstruct.zerod = z0dstruct.zs;
            % UAL option is not always preserve (security rule on parameters)
            noms = fieldnames(option);
            for kl=1:length(noms)
                if ~isfield(z0dstruct.z0dinput.option,noms{kl})
                    z0dstruct.z0dinput.option.(noms{kl}) = option.(noms{kl});
                end
            end

            
        case 'post'
            fprintf('Post processing of METIS simulation data\n');
            z0dstruct = post;
            % compatibilite ascendante
            z0dstruct.zs = z0dstruct.zerod;
            z0dstruct.profil = z0dstruct.profil0d;
            %
            data_zerod = z0dstruct.zerod;
            profil0d   = z0dstruct.profil0d;
            z0dinput   = z0dstruct.z0dinput;
            
        case 'post_one_time'
            fprintf('Post processing of METIS evolution simulation data\n');
            z0dstruct = post.z0dstruct;
            % compatibilite ascendante
            if ~isfield(z0dstruct,'zs')
                z0dstruct.zs = z0dstruct.zerod;
            elseif ~isfield(z0dstruct,'zerod')
                 z0dstruct.zerod = z0dstruct.zs;
            end
            if ~isfield(z0dstruct,'profil')
                 z0dstruct.profil = z0dstruct.profil0d;
            elseif ~isfield(z0dstruct,'profil0d')
                  z0dstruct.profil0d = z0dstruct.profil;
            end
            %
            if isfield(post,'data_zerod')
                    data_zerod = post.data_zerod;
            else
                    data_zerod = post.zs;
            end  
            if isfield(post,'profil0d')
                 profil0d   = post.profil0d;
            else
                 profil0d   = post.profil;              
            end
            z0dinput   = post.z0dinput;
            
        case 'read'
            % rien dans ce cas
        otherwise
            % we assume that the string method is a filename with complete path
            % this is the load method
            fprintf('loading METIS simulation %s\n',methode)
            z0dstruct = load_metis_imas(methode);
            % add codeparam to options
            info = metis4imas(1);
            optf = info.valeur;
            noms = fieldnames(optf);
            for k=1:length(noms)
                if ~isfield(z0dstruct.z0dinput.option,noms{k})
                    z0dstruct.z0dinput.option.(noms{k}) = optf.(noms{k});
                end
            end
            z0dstruct.z0dinput.option = secure_option_imas(z0dstruct.z0dinput.option);
            % compatibilite ascendante
            z0dstruct.zs = z0dstruct.zerod;
            z0dstruct.profil = z0dstruct.profil0d;
            %
            data_zerod = z0dstruct.zerod;
            profil0d   = z0dstruct.profil0d;
            z0dinput   = z0dstruct.z0dinput;
            
    end
    error_flag = 0;
    setappdata(0,'IMAS_Z0DSTRUCT',z0dstruct);
  
catch	
  ers = lasterror;
  if isstruct(ers)	
    fprintf('error message : %s\n',ers.message);	
    fprintf('error identifier : %s\n',ers.identifier);	
    for ke=1:length(ers.stack)
      fprintf('error stack %d in %s (%s) at %d\n',ke,ers.stack(ke).file,ers.stack(ke).name,ers.stack(ke).line);	
      
    end
  else
    disp(ers);
  end
  error_flag = sum(abs(ers.identifier));
  return
end

% debut traitement COCOS
% before COCOS transformation if any
sigma_B0_eff = 1;
COCOS_Sauter = false;
z0dstruct_loc = z0dstruct;
if ~isempty(which('ids_check_cocos'))
    switch  z0dstruct.z0dinput.option.COCOS_method
        case 'Sauter'
            COCOS_Sauter = true;
    end
end
if COCOS_Sauter
    % use Olivier Sauter tools for COCOS in this case
    sigma_B0_in = 1;
    sigma_Ip_in = 1;
    factor_two_pi = 1;
    sigma_bvac_r  = 1;
    %COCOS_Sauter = true;
    
    % change orientation and signe in z0dstruct if COCOS_Sauter to 1:
    z0dstruct_loc.z0dinput.option.signe = 1;
    z0dstruct_loc.z0dinput.option.orientation = 1;
   
     
elseif exist('profil0d','var') && exist('data_zerod','var')

  % applied sign that is not completely defined in METIS
  if isfield(z0dstruct.z0dinput.option,'orientation')
      sigma_B0_in = z0dstruct.z0dinput.option.orientation;
  else
      sigma_B0_in = 1;
  end
  disp('Signs are:')
  sigma_B0_in
  sigma_Ip_in = sign(z0dstruct.z0dinput.option.signe .* sigma_B0_in)
  % test neutral if sign of METIS
%    data_zerod_test = data_zerod;
%    profil0d_test   = profil0d;
%    [data_zerod,profil0d] = makesign_imas(data_zerod,profil0d,z0dstruct,1,1);
%    zcompstruct(data_zerod,data_zerod_test);
%    zcompstruct(profil0d,profil0d_test);
%    keyboard
  %
  [data_zerod,profil0d] = makesign_imas(data_zerod,profil0d,z0dstruct,sigma_Ip_in,sigma_B0_in);
  
  % changement normalisation psi
  factor_two_pi   = 2 .* pi;
  profil0d.psi    = + profil0d.psi    .* 2 .* pi;
  profil0d.dpsidt = + profil0d.dpsidt .* 2 .* pi;
  profil0d.qjli   = - profil0d.qjli;  
  profil0d.utheta = + profil0d.utheta ./ 2 ./ pi;
  data_zerod.q0   = - data_zerod.q0; 
  data_zerod.q95  = - data_zerod.q95;
  data_zerod.qa   = - data_zerod.qa;
  data_zerod.qmin = - data_zerod.qmin;
  data_zerod.qeff = - data_zerod.qeff;
  
  % pour test
  %  z0dstruct.z0dinput.option.signe = -1
  %  z0dstruct.z0dinput.option.COCOS = 11;
  
%    % mise en place signe Bphi
%    if z0dstruct.z0dinput.option.signe < 0
%      % on change le signe de b0
%      data_zerod.qa           = - data_zerod.qa;
%      data_zerod.q95          = - data_zerod.q95;
%      data_zerod.qmin         = - data_zerod.qmin;
%      data_zerod.qeff         = - data_zerod.qeff;
%      data_zerod.q0           = - data_zerod.q0;
%      data_zerod.phiplasma    = - data_zerod.phiplasma;
%      profil0d.qjli           = - profil0d.qjli;
%      profil0d.fdia           = - profil0d.fdia;
%      profil0d.phi            = - profil0d.phi;
%      profil0d.dphidx         = - profil0d.dphidx;
%    end
    
  % appel de la fonction qui regle le bon cocos
  [data_zerod,profil0d,sigma_B0_eff,sigma_bvac_r,factor_two_pi] = ...
      makecocos_imas(data_zerod,profil0d,z0dstruct,z0dstruct.z0dinput.option.COCOS,factor_two_pi);
      
end 

% Structure for O. Sauter tools
cocos_struct.COCOS_Sauter = COCOS_Sauter;
% METIS native COCOS is 7:
% q > 0, Psi decreasing, dP/dPsi positive and phi is counted counter clock-wise.
% So (R,phi,Z) for computing BR  and BZ with sigma_Bp = -1 and  sigma_RZphi = 1.
% Theta from front is clockwise (sigma_rho_theta_phi = 1 and coordinate are
% (rho,theta,phi)
% By default with signe = 1 and orientation = 1, Bphi and the
% plasma current are in the same direction and Bphi is positive and counter
% clockwise the tokamak se from above.
% Psi is in Wb/rad so e_Bp = 0
cocos_struct.cocos_in     = 7; % METIS native
cocos_struct.cocos_out    = z0dstruct.z0dinput.option.COCOS;
cocos_struct.ipsign_out   = sign(z0dstruct.z0dinput.option.signe .* z0dstruct.z0dinput.option.orientation);
cocos_struct.b0sign_out   = z0dstruct.z0dinput.option.orientation;
cocos_struct.ipsign_in    = 1; % METIS default
cocos_struct.b0sign_in    = 1; % METIS default
cocos_struct.error_bar    = 'none'; % no error bar in METIS output
cocos_struct.verbose      = z0dstruct.z0dinput.option.COCOS_verbose; % can be switch for debug
cocos_struct.check        = z0dstruct.z0dinput.option.COCOS_check; % can be switch for debug

disp('COCOS transform information:')
cocos_struct

% reading logfile    
diary off
% lecture du fichier
fid = fopen(diary_file,'r');
if fid > 0
  texte_diary = char(fread(fid,Inf,'char')');
  fclose(fid);
else
  texte_diary = 'unable to read diary file';
end
if diary_already_on
	diary on
else
	delete(diary_file);
end
% restart file
if save_restart_file
  metis4imas_save_restart(option.restart,z0dstruct,texte_diary);
end

% access to UAL
if write_mode == 1
    disp('METIS4IMAS: Opening database');
    
    %  creation of imasdb
    if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK')) && ~strcmp(treename,getappdata(0,'UAL_TOKAMAK'))

        al_version = getenv('AL_VERSION'); 
        if ~isempty(al_version)
            al_version = str2num(al_version(1));
        else
            al_version = 4;
        end
        if al_version < 5
            [s,t] = unix(sprintf('%s %s',fullfile(fileparts(fileparts( ...
                which('imas_open_env'))),'bin','imasdb'),getappdata(0,'UAL_TOKAMAK')));
        else
            s = 0;
        end
        if s == 0
            fprintf('METIS4IMAS: creating a new trename entry: %s\n',getappdata(0,'UAL_TOKAMAK'));
        else
            fprintf('METIS4IMAS: error creating a new treename %s:\n%s\n',getappdata(0,'UAL_TOKAMAK'),t);
            return
        end
        if isappdata(0,'UAL_USER') && isempty(getappdata(0,'UAL_USER'))
            setappdata(0,'UAL_USER',deblank(user));
        end
        if isappdata(0,'UAL_DATAVERSION') && isempty(getappdata(0,'UAL_DATAVERSION'))
            setappdata(0,'UAL_DATAVERSION',deblank(ver));
        end
    end
    try
        fprintf('METIS4IMAS: testing if entry exist in database -> ');
        if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK')) &&  ...
                isappdata(0,'UAL_USER') && ~isempty(getappdata(0,'UAL_USER')) &&  ...
                isappdata(0,'UAL_DATAVERSION') && ~isempty(getappdata(0,'UAL_DATAVERSION'));
            if ~isempty(id_backend)
                expIdx = mem_cache_imas( imas_open_env_backend(real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'),id_backend));
            else
                expIdx = mem_cache_imas( imas_open_env('ids',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
            end
        else
            expIdx = mem_cache_imas(imas_open('ids',real(shot),real(run)));
        end
        if expIdx < 0
            error('unexisting database');
        end
        fprintf('reusing previously created database entry\n');
    catch
        if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK')) &&  ...
                isappdata(0,'UAL_USER') && ~isempty(getappdata(0,'UAL_USER')) &&  ...
                isappdata(0,'UAL_DATAVERSION') && ~isempty(getappdata(0,'UAL_DATAVERSION'));
            if ~isempty(id_backend)
                expIdx = imas_create_env_backend(real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'),id_backend);
            elseif imag(shot) ~= 0
                expIdx = imas_create_env('ids',real(shot),real(run),imag(shot),imag(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'));
            else
                expIdx =imas_create_env('ids',real(shot),real(run),0,0,getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'));
            end
            disp('METIS4IMAS: creating a new entry in database');
        else
            if imag(shot) ~= 0
                expIdx = imas_create('ids',real(shot),real(run),imag(shot),imag(run));
            else
                expIdx = imas_create('ids',real(shot),real(run),0,0);
            end
            disp('METIS4IMAS: creating a new entry in database');
        end
    end
    
    if ~(isfinite(expIdx) && (expIdx >=0))
        disp('Unable to open data base or to connect to UAL')
        error_flag = -9001;
        return
    else
        expIdx = mem_cache_imas(expIdx);
    end
    
    % mapping IMAS structures
    ntime_summary = length(data_zerod.temps);
    ntime_core_profiles = length(profil0d.temps);
    %
    % IDS scenario rule is now play by pulse_chedule
    disp('Mapping IDS data structure pulse_schedule')
    pulse_schedule = mappulse_schedule_imas(z0dstruct_loc,data_zerod);
    pulse_schedule = change_cocos_OS(pulse_schedule,'pulse_schedule',cocos_struct);
    % recopie des informations
    summary             = ids_gen('summary');
    summary.code        = pulse_schedule.code;
    core_profiles       = ids_gen('core_profiles');
    core_profiles.code  = pulse_schedule.code;
    core_profiles.ids_properties  = pulse_schedule.ids_properties;
    core_transport      = ids_gen('core_transport');
    core_transport.code = pulse_schedule.code;
    core_transport.ids_properties  = pulse_schedule.ids_properties;
    core_sources       = ids_gen('core_sources');
    core_sources.code  = pulse_schedule.code;
    core_sources.ids_properties  = pulse_schedule.ids_properties;
    radiation       = ids_gen('radiation');
    radiation.code   = pulse_schedule.code;
    radiation.ids_properties  = pulse_schedule.ids_properties;
    transport_solver_numerics       = ids_gen('transport_solver_numerics');
    transport_solver_numerics.code  = pulse_schedule.code;
    transport_solver_numerics.ids_properties  = pulse_schedule.ids_properties;
    equilibrium         = ids_gen('equilibrium');
    equilibrium.code    = pulse_schedule.code;
    equilibrium.ids_properties  = pulse_schedule.ids_properties;
    %
    disp('Writing IDS data structure pulse_schedule')
    if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
        nids = 'pulse_schedule';
    else
        nids = sprintf('pulse_schedule/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
    end
    if init_output_ids
        ids_put_non_timed(expIdx,nids,pulse_schedule);
    elseif  ntime_summary > 1
        ids_put(expIdx,nids,pulse_schedule);
    else
        ids_put_slice(expIdx,nids,pulse_schedule);
    end
    % information are also in  dataset_description
    disp('Mapping IDS data structure dataset_description')
    dataset_description = mapdataset_description_imas(z0dstruct_loc,data_zerod,texte_diary,error_flag,run,shot,pulse_schedule.code);
    dataset_description = change_cocos_OS(dataset_description,'dataset_description',cocos_struct);
    %
    disp('Writing IDS data structure dataset_description')
    if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
        nids = 'dataset_description';
    else
        nids = sprintf('dataset_description/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
    end
    if init_output_ids
        ids_put_non_timed(expIdx,nids,dataset_description);
    elseif  ntime_summary > 1
        ids_put(expIdx,nids,dataset_description);
    else
        ids_put_slice(expIdx,nids,dataset_description);
    end
    
    % dataset_fair IDS
    disp('Mapping IDS data structure dataset_fair')
    dataset_fair = map_dataset_fair(dataset_description,dataset_fair_in,real(shot),real(run));
    %%%dataset_fair = change_cocos_OS(dataset_fair,'dataset_fair',cocos_struct);
    if ~isempty(dataset_fair)
        %
        disp('Writing IDS data structure dataset_fair')
        if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
            nids = 'dataset_fair';
        else
            nids = sprintf('dataset_fair/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,dataset_fair);
        elseif  ntime_summary > 1
            ids_put(expIdx,nids,dataset_fair);
        else
            ids_put_slice(expIdx,nids,dataset_fair);
        end
    end
    
    % now summary must be seen as other IDSs
    if z0dstruct.z0dinput.option.summary == 1
        disp('Mapping IDS data structure summary')
        summary = mapsummary_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,summary,run,occurrence,sigma_B0_eff,sigma_bvac_r);
        summary = change_cocos_OS(summary,'summary',cocos_struct);
        disp('Writing IDS data structure summary')
        if isempty(strtrim(z0dstruct.z0dinput.option.summary_occurrence))
            nids = 'summary';
        else
            nids = sprintf('summary/%s',strtrim(z0dstruct.z0dinput.option.summary_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,summary);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,summary);
        else
            ids_put_slice(expIdx,nids,summary);
        end
    end
    %
    if z0dstruct.z0dinput.option.core_profiles == 1
        disp('Mapping IDS data structure core_profiles')
        core_profiles = mapcore_profiles_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,core_profiles,sigma_B0_eff);
        core_profiles = change_cocos_OS(core_profiles,'core_profiles',cocos_struct);
        disp('Writing IDS data structure core_profiles')
        if isempty(strtrim(z0dstruct.z0dinput.option.core_profiles_occurrence))
            nids = 'core_profiles';
        else
            nids = sprintf('core_profiles/%s',strtrim(z0dstruct.z0dinput.option.core_profiles_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,core_profiles);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,core_profiles);
        else
            ids_put_slice(expIdx,nids,core_profiles);
        end
    end
    if z0dstruct.z0dinput.option.core_transport == 1
        disp('Mapping IDS data structure core_transport')
        core_transport = mapcore_transport_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,core_transport,sigma_B0_eff);
        core_transport = change_cocos_OS(core_transport,'core_transport',cocos_struct);
        disp('Writing IDS data structure core_transport')
        if isempty(strtrim(z0dstruct.z0dinput.option.core_transport_occurrence))
            nids = 'core_transport';
        else
            nids = sprintf('core_transport/%s',strtrim(z0dstruct.z0dinput.option.core_transport_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,core_transport);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,core_transport);
        else
            ids_put_slice(expIdx,nids,core_transport);
        end
    end
    count_sources = 0; % Need to be initialized when we map no source at all
    if z0dstruct.z0dinput.option.core_sources == 1
        model_source = core_sources.source{1};
        disp('Mapping IDS data structure core_sources (generic part)')
        core_sources = mapcore_sources_generic_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,core_sources,sigma_B0_eff);
        count_sources = 1;
        disp('Mapping IDS data structure core_sources for lhcd')
        core_sources.source{count_sources} = mapcore_sources_lhcd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources for eccd')
        core_sources.source{count_sources} = mapcore_sources_eccd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources for radiation')
        core_sources.source{count_sources} = mapcore_sources_radiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        core_sources.source{count_sources} = mapcore_sources_brem_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources for cyclotron radiation')
        core_sources.source{count_sources} = mapcore_sources_cyclotron_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources ICRH')
        core_sources.source{count_sources} = mapcore_sources_icrh_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources NBI')
        core_sources.source{count_sources} = mapcore_sources_nbicd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources for fusion')
        core_sources.source{count_sources} = mapcore_sources_fusion_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources ohm & boot')
        [core_sources.source{count_sources},core_sources.source{count_sources + 1}] = mapcore_sources_ohm_boot_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 2;
        disp('Mapping IDS data structure core_sources neutrals')
        core_sources.source{count_sources} = mapcore_sources_neutral_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources equipartition')
        core_sources.source{count_sources} = mapcore_sources_equipartition_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        disp('Mapping IDS data structure core_sources_full')
        core_sources.source{count_sources} = mapcore_sources_full_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
        count_sources = count_sources + 1;
        % extra source for runaways
        if z0dstruct.z0dinput.option.runaway > 0
            disp('Mapping IDS data structure core_sources for runaway electron current')
            core_sources.source{count_sources} = mapcore_sources_runaways_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
            count_sources = count_sources + 1;
        end
        core_sources = change_cocos_OS(core_sources,'core_sources',cocos_struct);
    end
    
    % there i no IDS neoclassic ? core transport ?
    %neoclassic = mapneoclassic_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,neoclassic,vtor,vpol,sigma_B0_eff);
    
    % it is core_profiles.profil_1d{kt}.neutral{m}
    %coreneutrals = mapcore_neutrals_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreneutrals,sigma_B0_eff);
    %coreneutrals = mapcore_neutrals_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreneutrals,sigma_B0_eff);
    
    % for super thermal energy the data must be stored in distribution(:)/global_quantities(:)
    
    % ou sont le neutrons ?
    
    % ecriture ids core_sources
    if count_sources > 1
        disp('Writing IDS data structure core_sources')
        if isempty(strtrim(z0dstruct.z0dinput.option.core_sources_occurrence))
            nids = 'core_sources';
        else
            nids = sprintf('core_sources/%s',strtrim(z0dstruct.z0dinput.option.core_sources_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,core_sources);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,core_sources);
        else
            ids_put_slice(expIdx,nids,core_sources);
        end
    end
    
    if z0dstruct.z0dinput.option.radiation == 1
        disp('Mapping IDS data structure radiation')
        radiation = mapradiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,radiation,sigma_B0_eff);
        radiation = change_cocos_OS(radiation,'radiation',cocos_struct);
        disp('Writing IDS data structure radiation')
        if isempty(strtrim(z0dstruct.z0dinput.option.radiation_occurrence))
            nids = 'radiation';
        else
            nids = sprintf('radiation/%s',strtrim(z0dstruct.z0dinput.option.radiation_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,radiation);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,radiation);
        else
            ids_put_slice(expIdx,nids,radiation);
        end
    end
    
    switch z0dstruct.z0dinput.option.sol_model
        case '2_points'
            if z0dstruct.z0dinput.option.edge == 1
                disp('Mapping IDSs data structure edge_profiles and edge_transport')
                [edge_profiles,edge_transport] = mapedge_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,pulse_schedule,sigma_B0_eff);
                edge_profiles = change_cocos_OS(edge_profiles,'edge_profiles',cocos_struct);
                edge_transport = change_cocos_OS(edge_transport,'edge_profiles',cocos_struct);
                disp('Writing IDS data structure edge_profiles')
                if isempty(strtrim(z0dstruct.z0dinput.option.edge_occurrence))
                    nids = 'edge_profiles';
                else
                    nids = sprintf('edge_profiles/%s',strtrim(z0dstruct.z0dinput.option.edge_occurrence));
                end
                if init_output_ids
                    ids_put_non_timed(expIdx,nids,edge_profiles);
                elseif ntime_core_profiles  > 1
                    ids_put(expIdx,nids,edge_profiles);
                else
                    ids_put_slice(expIdx,nids,edge_profiles);
                end
                disp('Writing IDS data structure edge_transport')
                if isempty(strtrim(z0dstruct.z0dinput.option.edge_occurrence))
                    nids = 'edge_transport';
                else
                    nids = sprintf('edge_transport/%s',strtrim(z0dstruct.z0dinput.option.edge_occurrence));
                end
                if init_output_ids
                    ids_put_non_timed(expIdx,nids,edge_transport);
                elseif ntime_core_profiles  > 1
                    ids_put(expIdx,nids,edge_transport);
                else
                    ids_put_slice(expIdx,nids,edge_transport);
                end
            end
    end
    if z0dstruct.z0dinput.option.numerics == 1
        disp('Mapping IDS data structure transport_solver_numerics')
        transport_solver_numerics = maptransport_solver_numerics(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,transport_solver_numerics,run,occurrence,sigma_B0_eff,sigma_bvac_r);
        transport_solver_numerics = change_cocos_OS(transport_solver_numerics,'transport_solver_numerics',cocos_struct);
        disp('Writing IDS data structure transport_solver_numerics')
        if isempty(strtrim(z0dstruct.z0dinput.option.numerics_occurrence))
            nids = 'transport_solver_numerics';
        else
            nids = sprintf('transport_solver_numerics/%s',strtrim(z0dstruct.z0dinput.option.numerics_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,transport_solver_numerics);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,transport_solver_numerics);
        else
            ids_put_slice(expIdx,nids,transport_solver_numerics);
        end
    end
    
    if z0dstruct.z0dinput.option.equilibrium == 1
        disp('Mapping IDS data structure equilibrium')
        equilibrium = mapequilibrium_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,equilibrium,summary, ...
            1,sigma_B0_eff,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi);
        equilibrium = change_cocos_OS(equilibrium,'equilibrium',cocos_struct);
        disp('Writing IDS data structure equilibrium')
        if isempty(strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence))
            nids = 'equilibrium';
        else
            nids = sprintf('equilibrium/%s',strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence));
        end
        if init_output_ids
            ids_put_non_timed(expIdx,nids,equilibrium);
        elseif ntime_core_profiles  > 1
            ids_put(expIdx,nids,equilibrium);
        else
            ids_put_slice(expIdx,nids,equilibrium);
        end
    end
    
    if getappdata(0,'imas_enable_mem_cache_imas') == 1
        imas_flush_all(expIdx)
    end
    imas_close(expIdx);
    disp('End of metis4imas writing');
    
else
  switch methode
   case {'post','post_one_time'}
    %% CREATE MATLAB OUTPUT STRUCTURE    
    disp('Mapping IDS data structure pulse_schedule')
    xsd.pulse_schedule = mappulse_schedule_imas(z0dstruct_loc,data_zerod);
    disp('Mapping IDS data structure dataset_description')
    xsd.dataset_description = mapdataset_description_imas(z0dstruct_loc,data_zerod,texte_diary,error_flag,run,shot,xsd.pulse_schedule.code);
    xsd.dataset_fair        = map_dataset_fair(xsd.dataset_description,dataset_fair_in);
    xsd.summary             = ids_gen('summary');
    xsd.summary.code        = xsd.pulse_schedule.code;
    xsd.core_profiles       = ids_gen('core_profiles');
    xsd.core_profiles.code  = xsd.pulse_schedule.code;
    xsd.core_profiles.ids_properties  =  xsd.pulse_schedule.ids_properties;
    xsd.core_transport      = ids_gen('core_transport');
    xsd.core_transport.code = xsd.pulse_schedule.code;
    xsd.core_transport.ids_properties  = xsd.pulse_schedule.ids_properties;
    xsd.core_sources       = ids_gen('core_sources');
    xsd.core_sources.code  = xsd.pulse_schedule.code;
    xsd.core_sources.ids_properties  = xsd.pulse_schedule.ids_properties;
    xsd.radiation          = ids_gen('radiation');
    xsd.radiation.code     = xsd.pulse_schedule.code;
    xsd.radiation.ids_properties  = xsd.pulse_schedule.ids_properties;
    xsd.transport_solver_numerics                 = ids_gen('transport_solver_numerics');
    xsd.transport_solver_numerics.code            = xsd.pulse_schedule.code;
    xsd.transport_solver_numerics.ids_properties  =  xsd.pulse_schedule.ids_properties;
    xsd.equilibrium         = ids_gen('equilibrium');
    xsd.equilibrium.code    = xsd.pulse_schedule.code;
    xsd.equilibrium.ids_properties  = xsd.pulse_schedule.ids_properties;

     
    if z0dstruct.z0dinput.option.summary == 1
      disp('Mapping IDS data structure summary')
      xsd.summary = mapsummary_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.summary,run,occurrence,sigma_B0_eff,sigma_bvac_r);
    end
    if z0dstruct.z0dinput.option.core_profiles == 1
      disp('Mapping IDS data structure core_profiles')
      xsd.core_profiles = mapcore_profiles_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.core_profiles,sigma_B0_eff);
    end
    if z0dstruct.z0dinput.option.core_transport == 1
      disp('Mapping IDS data structure core_transport')
      xsd.core_transport = mapcore_transport_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.core_transport,sigma_B0_eff);
    end
    if z0dstruct.z0dinput.option.core_sources == 1
      model_source = xsd.core_sources.source{1};
      disp('Mapping IDS data structure core_sources (generic part)')
      xsd.core_sources = mapcore_sources_generic_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.core_sources,sigma_B0_eff);
      count_sources = 1;
      disp('Mapping IDS data structure core_sources for lhcd')
      xsd.core_sources.source{count_sources} = mapcore_sources_lhcd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for eccd')
      xsd.core_sources.source{count_sources} = mapcore_sources_eccd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for radiation')
      xsd.core_sources.source{count_sources} = mapcore_sources_radiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      xsd.core_sources.source{count_sources} = mapcore_sources_brem_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for cyclotron')
      xsd.core_sources.source{count_sources} = mapcore_sources_cyclotron_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for icrh')
      xsd.core_sources.source{count_sources} = mapcore_sources_icrh_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for nbicd')
      xsd.core_sources.source{count_sources} = mapcore_sources_nbicd_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for fusion')
      xsd.core_sources.source{count_sources} = mapcore_sources_fusion_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources ohm & boot')
      [xsd.core_sources.source{count_sources},xsd.core_sources.source{count_sources + 1}] = mapcore_sources_ohm_boot_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 2;
      disp('Mapping IDS data structure core_sources neutrals')
      xsd.core_sources.source{count_sources} = mapcore_sources_neutral_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources equipartition')
      xsd.core_sources.source{count_sources} = mapcore_sources_equipartition_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      disp('Mapping IDS data structure core_sources for sum of sources')
      xsd.core_sources.source{count_sources} = mapcore_sources_full_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
      count_sources = count_sources + 1;
      % extra source for runaways
      if z0dstruct.z0dinput.option.runaway > 0
		disp('Mapping IDS data structure core_sources for runaway electron current')
		core_sources.source{count_sources} = mapcore_sources_runaways_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff,model_source);
		count_sources = count_sources + 1;
      end
    end
    if z0dstruct.z0dinput.option.radiation == 1
      disp('Mapping IDS data structure radiation')
      xsd.radiation = mapradiation_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.radiation,sigma_B0_eff);
    end
    switch z0dstruct.z0dinput.option.sol_model
    case '2_points'
      if z0dstruct.z0dinput.option.edge == 1
	disp('Mapping IDS data structure edge_profiles and edge_transport')
	[xsd.edge_profiles,xsd.edge_transport] = mapedge_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.pulse_schedule,sigma_B0_eff);
      end
    end
    if z0dstruct.z0dinput.option.numerics == 1 
	disp('Mapping IDS data structure transport_solver_numerics')
        xsd.transport_solver_numerics = maptransport_solver_numerics(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.transport_solver_numerics,run,occurrence,sigma_B0_eff,sigma_bvac_r);
    end
    if z0dstruct.z0dinput.option.equilibrium == 1
      disp('Mapping IDS data structure equilibrium')
      if isfield(xsd,'summary')
        xsd.equilibrium = mapequilibrium_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.equilibrium,xsd.summary, ...
					      1,sigma_B0_eff,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi);
      else
        xsd.equilibrium = mapequilibrium_imas(z0dstruct_loc,data_zerod,profil0d,texte_diary,error_flag,xsd.equilibrium,[], ...
					       1,sigma_B0_eff,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi);     
      end
    end 
    
    % change all COCOS is xsd
    noms = fieldnames(xsd);
    for k=1:length(noms)
        xsd.(noms{k}) = change_cocos_OS(xsd.(noms{k}),noms{k},cocos_struct);
    end
  end
end

if read_output == 1
  %% OUTPUT VARIABLE
  xsd = [];
  try
      if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
          if ~isempty(id_backend)
              expIdx = mem_cache_imas( imas_open_env_backend(real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION'),id_backend));
          else
              expIdx = mem_cache_imas( imas_open_env('ids',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
          end
      else
          expIdx = mem_cache_imas(imas_open('ids',real(shot),real(run)));
      end
  end
  
  if ~(isfinite(expIdx) & (expIdx >=0))
    disp('Unable to open data base or to connect to UAL')
    error_flag = -9001;
    return
  end
  if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
    nids = 'pulse_schedule';
  else
    nids = sprintf('pulse_schedule/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
  end
  %% GET IDS 
  disp('Reading IDS data structure pulse_schedule')
  try
    pulse_schedule_ids_read = ids_get(expIdx,nids);
  catch	
    pulse_schedule_ids_read = [];	    
  end
  if ~isempty(pulse_schedule_ids_read)
    xsd.pulse_schedule = pulse_schedule_ids_read;
  else
    disp('Unable to read pulse_schedule IDS');	
    xsd.pulse_schedule = [];
  end  
  % dataset_description
  if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
    nids = 'dataset_description';
  else
    nids = sprintf('dataset_description/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
  end
  % GET  IDS
  disp('Reading IDS data structure dataset_description')
  try
    dataset_description_ids_read = ids_get(expIdx,nids);
  catch	
    dataset_description_ids_read = [];	    
  end
  if ~isempty(pulse_schedule_ids_read)
    xsd.dataset_description = dataset_description_ids_read;
  else
    disp('Unable to read dataset_description IDS');	
    xsd.dataset_description = [];
  end  
  % dataset_fair
  if isempty(strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence))
    nids = 'dataset_fair';
  else
    nids = sprintf('dataset_fair/%s',strtrim(z0dstruct.z0dinput.option.pulse_schedule_occurrence));
  end
  % GET  IDS
  disp('Reading IDS data structure dataset_fair')
  try
    dataset_fair_ids_read = ids_get(expIdx,nids);
  catch	
    dataset_fair_ids_read = [];	    
  end
  if ~isempty(pulse_schedule_ids_read)
    xsd.dataset_fair = dataset_fair_ids_read;
  else
    disp('Unable to read dataset_fair IDS');	
    xsd.dataset_fair = [];
  end  
  % summary
  if z0dstruct.z0dinput.option.summary == 1
    if isempty(strtrim(z0dstruct.z0dinput.option.summary_occurrence))
	nids = 'summary';
    else
	nids = sprintf('summary/%s',strtrim(z0dstruct.z0dinput.option.summary_occurrence));
    end
    %% GET IDS 
    disp('Reading IDS data structure summary')
    try
	summary_ids_read = ids_get(expIdx,nids);
    catch	
	summary_ids_read = [];		
    end
    if ~isempty(summary_ids_read)
	xsd.summary = summary_ids_read;
    else
        disp('Unable to read summary IDS');	
	xsd.summary = [];
     end  
  end
  if z0dstruct.z0dinput.option.core_profiles == 1
    %% GET IDS 
    disp('Reading IDS data structure core_profiles')
    if isempty(strtrim(z0dstruct.z0dinput.option.core_profiles_occurrence))
      nids = 'core_profiles';
    else
      nids = sprintf('core_profiles/%s',strtrim(z0dstruct.z0dinput.option.core_profiles_occurrence));
    end
    try
      core_profiles_ids_read =ids_get(expIdx,nids);
    catch	
      core_profiles_ids_read = [];	
		  
    end
    if ~isempty(core_profiles_ids_read)
      xsd.core_profiles = core_profiles_ids_read;
    else
      disp('Unable to read core_profiles IDS');	
      xsd.core_profiles = [];
    end
    
  end
  if z0dstruct.z0dinput.option.core_transport == 1
    %% GET IDS 
    disp('Reading IDS data structure core_transport')
    if isempty(strtrim(z0dstruct.z0dinput.option.core_transport_occurrence))
      nids = 'core_transport';
    else
      nids = sprintf('core_transport/%s',strtrim(z0dstruct.z0dinput.option.core_transport_occurrence));
    end
    try
      core_transport_ids_read =ids_get(expIdx,nids);
    catch	
      core_transport_ids_read = [];	     
    end
    if ~isempty(core_transport_ids_read)
      xsd.core_transport = core_transport_ids_read;
    else
      disp('Unable to read core_transport IDS');	
      xsd.core_transport = [];
    end
  end
  if isempty(strtrim(z0dstruct.z0dinput.option.core_sources_occurrence))
    nids = 'core_sources';
  else
    nids = sprintf('core_sources/%s',strtrim(z0dstruct.z0dinput.option.core_sources_occurrence));
  end
  %% GET IDS 
  disp('Reading IDS data structure core_sources')
  try
    core_sources_ids_read =ids_get(expIdx,nids);
  catch	
    core_sources_ids_read = [];	
    
  end
  if ~isempty(core_sources_ids_read)
    xsd.core_sources = core_sources_ids_read;
  else
    disp('Unable to read core_sources IDS');	
    xsd.core_sources = [];
  end
  if isempty(strtrim(z0dstruct.z0dinput.option.radiation_occurrence))
    nids = 'radiation';
  else
    nids = sprintf('radiation/%s',strtrim(z0dstruct.z0dinput.option.radiation_occurrence));
  end
  %% GET IDS 
  disp('Reading IDS data structure radiation')
  try
    radiation_ids_read =ids_get(expIdx,nids);
  catch	
    radiation_ids_read = [];	    
  end
  if ~isempty(radiation_ids_read)
    xsd.radiation = radiation_ids_read;
  else
    disp('Unable to read radiation IDS');	
    xsd.radiation = [];
  end
  if z0dstruct.z0dinput.option.edge == 1
    %% GET IDS 
    disp('Reading IDS data structure edge profiles and edge_transport')
    if isempty(strtrim(z0dstruct.z0dinput.option.edge_occurrence))
      nids = 'edge_profiles';
    else
      nids = sprintf('edge_profiles/%s',strtrim(z0dstruct.z0dinput.option.edge_occurrence));
    end
    try
      edge_profiles_ids_read =ids_get(expIdx,nids);
    catch	
      edge_profiles_ids_read = [];	     
    end
    if ~isempty(edge_profiles_ids_read)
      xsd.edge_profiles = edge_profiles_ids_read;
    else
      disp('Unable to read ege_profiles IDS');	
      xsd.edge_transport = [];
    end
    if isempty(strtrim(z0dstruct.z0dinput.option.edge_occurrence))
      nids = 'edge_transport';
    else
      nids = sprintf('edge_transport/%s',strtrim(z0dstruct.z0dinput.option.edge_occurrence));
    end
    try
      edge_transport_ids_read =ids_get(expIdx,nids);
    catch	
      edge_transport_ids_read = [];	     
    end
    if ~isempty(edge_transport_ids_read)
      xsd.edge_transport = edge_transport_ids_read;
    else
      disp('Unable to read edge_transport IDS');	
      xsd.edge_transport = [];
    end
  end
  
  if z0dstruct.z0dinput.option.numerics == 1
    %% GET IDS 
    disp('Reading IDS data structure transport_solver_numerics')
    if isempty(strtrim(z0dstruct.z0dinput.option.numerics_occurrence))
      nids = 'transport_solver_numerics';
    else
      nids = sprintf('transport_solver_numerics/%s',strtrim(z0dstruct.z0dinput.option.numerics_occurrence));
    end
    try
      transport_solver_numerics_ids_read =ids_get(expIdx,nids);
    catch	
      transport_solver_numerics_ids_read = [];			  
    end
    if ~isempty(transport_solver_numerics_ids_read)
      xsd.transport_solver_numerics = transport_solver_numerics_ids_read;
    else
      disp('Unable to read transport_solver_numerics IDS');	
      xsd.transport_solver_numerics = [];
    end
    
  end

  
  if z0dstruct.z0dinput.option.equilibrium == 1
    disp('Reading IDS data structure equilibrium')
    if isempty(strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence))
      nids = 'equilibrium';
    else
      nids = sprintf('equilibrium/%s',strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence));
    end
    try
      equilibrium_ids_read = ids_get(expIdx,nids);
    catch	
      equilibrium_ids_read = [];	
      
    end
    if ~isempty(equilibrium_ids_read)
      xsd.equilibrium = equilibrium_ids_read;
    else
      disp('Unable to read core_sources IDS');	
      xsd.equilibrium = [];
    end
  end
  
  imas_close(expIdx);
  disp('End of metis4imas reading');
  
end

