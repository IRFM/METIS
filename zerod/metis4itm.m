% METIS4ITM : interface between ITM data and METIS data
%-------------------------------------------------------------
% fonction Matlab 7; metis4itm.m -> metis4itm, mapscenario, test_metis_itm, load_metis_itm, scenario2option
%					cpo2metis_input, cpo2metis1t, mapcoreprof, interp1_itm, griddata_itm, z0rot_itm
%
%
% This function is the interface between ITM/CPO and METIS data structure
%
% syntax : 
%
%   * create xml an xsd files (without output):
%	metis4itm;
%
%   * parameters declaration stucture (with one output):
%	info = metis4itm;
%
%   * returning xml an xsd text (with 2 outputs):
%	[xml,xsd] = metis4itm;
%
%   * computing :
%	[error_flag,output_data] = metis4itm(shot,run,occurrence,method,time,codeparam_filename,interpolation_methode);
%
%   * post processing (with UAl data writing) :
%	error_flag = metis4itm(shot,run,occurrence,metis_data_structure);
%
%   * post processing (without UAl data writing = return only matlab data structure counterpart of CPOs structures) :
%	[error_flag,output_data] = metis4itm(shot,run,occurrence,<metis_data_structure>,[],optionnal_codeparam_filename);
%
%   * wrinting a METIS file in database:
%	error_flag = metis4itm(shot,run,occurrence,'metis_file_name',[],optionnal_codeparam_filename);
%
%   * reading a UAL reccord:
%	[error_flag,output_data] = metis4itm(shot,run,occurrence,'read',[],optionnal_codeparam_filename);
%
%
% input :
%
%     shot                    = shot number (integer >0)
%				if is complex, imag(shot) is used as reference shot
%     run                     = run number for this shot (integer >0)
%				if is complex, imag(shot) is used as reference run
%
%     occurrence              = input cpo scenario occurrence in Kepler (default = [])
%                               or occurrence can be a structure to use an alternative database :
%					occurrence.tokamak = new UAL tokamak database name
%					occurrence.user    = new UAL user database name
%                                       occurrence.dataversion = selected a differente version of data (not the last one)
%		                        occurrence.occurrence = input cpo scenario occurrence in Kepler (default = [])
%
%     method                  = metis command for Kepler actor or METIS file name (default = '') or metis data structure
%     time                    = time scalar or vector for which MESTIS must be run (scalar). 
%                               If time is scalar, METIS is run in evolution mode.
%                               If time is vector, METIS is run for whole shot simulation.
%                               If empty, reads all time slice of the CPO scenario and uses scenario CPO vector time as input.
%
%     codeparam_filename      = is name of the file containing the METIS parameters (XML),
%                               or XML string coding for codeparam,
%                               or matlab structure coding for METIS (all or some parameters,codeparam_filename.key = value).
%
%                               this input is optionnal. Defaults METIS parameters values are used for missing parameters. 
%                               this input is not usedin test mode.
%                               if method is 'read' ,'wrinting' or 'post processing' , the codeparam input allows to sellected the CPOs
%                               that are read or write (only the CPO scenario is mandatory)
%
%     interpolation_methode   = interpolation method used for data extraction (see UAL user manual).
%                               1 = closest sample, 2 = previous sample & 3 = linear interpolation.
%                               optionnal input.
%
% output :
%
%      err_flag		      = error code, must be = 0 if no error (as for unix program)
%      output_data            = matlab structure containing the couterpart of CPOs structures in matlab compatible format.
%                               Data are re-read in the UAL juste after wrinting.
%
%  method list: 
% 
%	'test'                 = test the function (call METIS with test data) and create CPOs.
%       'auto_test'            = complete test of metis4itm (all methods)
%       'full'                 = complete shot computation (all time slices given in scenario), full METIS computation mode
%       'fast'                 = complete shot computation (selected time slices given in scenario), fast METIS computation mode
%       'init'                 = initialisation of the evolution mode of METIS, scenario must contain the fisrt time slice for references.
%       'one_time'             = compute plasma evolution for one time step (used evolution mode of METIS), 
%                                scenario must contain the next time slice for references.
%        'read'                = read UAL data and return in matlab structure
%
%       <filename>             = load a METIS simulation and create associated CPO matlab data structures.
%                                (this is also the method for make a restart)
%
%  if method is empty, then the choice is automatic:
%      - if scenario is empty, method = 'test'
%      - if scenario contains only one time slice, a the fisrt call, method = 'init', 
%        for next calls, method = 'one_time' 
%      - if scenario contains more of one time slice, method = 'full'
%
%
%  Information about METIS is avalaible in METIS technical documentation that is included in the CRONOS distribution or with the button 
%  <HELP> of METIS GUI (PDF reader must be available).
%
%  Script to call METIS from Kepler :
%
%	addpath <path to cronos project>
%	zineb_path
%
%       % mapping of Kepler matlab actor input variables
%
%       shot        = Kepler_input_num_shot_variable;
%       run         = Kepler_input_num_shot_run;
%
%       % optionnal :
%	%occurrence = Kepler_input_metis_cpo_occurrence;
%       occurrence  = [];
%
%       % optionnal :
%	%method      = Kepler_input_metis_method;
%       method      = '';
%
%	% optionnal :
%       %time        = Kepler_vetor_or_scalar_input_time;
%	time        = [];
%
%	%optionnal
%	%codeparam_filename = Kepler_file_name_or_xml_string;
%	codeparam_filename = '';
%
%	% optionnal : 
%       %interpolation_method = 2;
%
%       %error_flag  = metis4itm(shot,run,occurrence,method,time,codeparam_filename,interpolation_method);
%       error_flag  = metis4itm(shot,run,occurrence,method,time);
%
%
%	% return to kepler
%	Kepler_return_value = error_flag;
%
%
%  access to metis internal data :
%
%	z0dstruct = getappdata(0,'ITM_Z0DSTRUCT');
%
%  simple test of the function :
%
%       [error_lfag,output] = metis4itm(1,1);
% 
%  full test of the function :
%
%  	metis4itm(1,1,'','auto_test')  
%
% fonction ecrite par J-F Artaud
% version CVS (created the 10/20/08)
%-----------------------------------------------------------------------
%
function [error_flag,xsd] = metis4itm(shot,run,occurrence,methode,time,codeparam_filename,interpolation_methode)

%link with UAL
if ~isappdata(0,'UALVERSION')
	try
		rep = javaclasspath('-all');
		for l=1:length(rep)
			sl = rep{l}; 
			ind = findstr(sl,'ualjava.jar');
			if ~isempty(ind)
				[f,s,e] = fileparts(fileparts(fileparts(sl)));
				%warning off
				%addpath('/afs/efda-itm.eu/project/switm/');
				%addpath(sprintf('/afs/efda-itm.eu/project/switm/ual/%s%s/matlabinterface',s,e));
				%warning on
		
		setappdata(0,'UALVERSION',sprintf('%s%s',s,e));
		% use :  any(getappdata(0,'UALVERSION') >'4.08b') to test if
		% ual versionis greater
		
			end
		end
	end
end
%  addpath('/afs/efda-itm.eu/project/switm/');
%  addpath('/afs/efda-itm.eu/project/switm/ual/4.07b/matlabinterface');
if isempty(import) && ~isappdata(0,'JAVAINTERFACESET') && isempty(getappdata(0,'JAVAINTERFACESET'))
	import ualmemory.javainterface.*;
end
setappdata(0,'JAVAINTERFACESET','done')

% gestion du  nombre d'entrees
if (nargin == 0) && ((nargout == 0)  || (nargout == 2))
	help metis4itm
	% create codeparam files
	module_xsd_make('metis4itm',1);
	% ouput codeparam
	[xsd,error_flag]=module_xsd_make('metis4itm',1);
	return
end
if nargin <= 1
	error_flag = zerod;
	% variables propres a metis4itm
        % initialisation des CPO au premier appel
    	error_flag.valeur.init_output_cpo    	= 0;
    	error_flag.type.init_output_cpo      	= 'integer';
    	error_flag.borne.init_output_cpo     	= {0,1}; 
    	error_flag.defaut.init_output_cpo   	= 0;
    	error_flag.info.init_output_cpo     	= 'if=1, on init call, initialise output CPO (reset and write the first time slice)';
        error_flag.section.init_output_cpo      ='UAL';
        
        % restart file name
    	error_flag.valeur.restart    	= '';
    	error_flag.type.restart     	= 'string';
    	error_flag.borne.restart     	= ''; 
    	error_flag.defaut.restart   	= '';
    	error_flag.info.restart     	= 'if empty, nor restart file save, otherwise after each call, the restart file is saved';
        error_flag.section.restart      ='UAL';

        % declaration des cop a ecrire
    	error_flag.valeur.coreprof    	= 1;
    	error_flag.type.coreprof      	= 'integer';
    	error_flag.borne.coreprof     	= {0,1}; 
    	error_flag.defaut.coreprof   	= 1;
    	error_flag.info.coreprof      	= 'if=1, write coreprof cpo';
        error_flag.section.coreprof      ='UAL';

        error_flag.valeur.coretransp    	= 1;
    	error_flag.type.coretransp      	= 'integer';
    	error_flag.borne.coretransp     	= {0,1}; 
    	error_flag.defaut.coretransp   		= 1;
    	error_flag.info.coretransp      	= 'if=1, write coretransp cpo';
        error_flag.section.coretransp      ='UAL';

    	error_flag.valeur.coresource_lhcd    	= 1;
    	error_flag.type.coresource_lhcd      	= 'integer';
    	error_flag.borne.coresource_lhcd     	= {0,1}; 
    	error_flag.defaut.coresource_lhcd   	= 1;
    	error_flag.info.coresource_lhcd      	= 'if=1, write coresource_lhcd cpo (coresource occurrence 1)';
        error_flag.section.coresource_lhcd      ='UAL';

    	error_flag.valeur.coresource_eccd    	= 1;
    	error_flag.type.coresource_eccd      	= 'integer';
    	error_flag.borne.coresource_eccd     	= {0,1}; 
    	error_flag.defaut.coresource_eccd   	= 1;
    	error_flag.info.coresource_eccd      	= 'if=1, write coresource_eccd cpo (coresource occurrence 2)';
        error_flag.section.coresource_eccd      ='UAL';

    	error_flag.valeur.coresource_icrh    	= 1;
    	error_flag.type.coresource_icrh      	= 'integer';
    	error_flag.borne.coresource_icrh     	= {0,1}; 
    	error_flag.defaut.coresource_icrh   	= 1;
    	error_flag.info.coresource_icrh      	= 'if=1, write coresource_icrh cpo (coresource occurrence 5)';
        error_flag.section.coresource_icrh      ='UAL';

    	error_flag.valeur.coresource_nbicd    	= 1;
    	error_flag.type.coresource_nbicd      	= 'integer';
    	error_flag.borne.coresource_nbicd     	= {0,1}; 
    	error_flag.defaut.coresource_nbicd   	= 1;
    	error_flag.info.coresource_nbicd      	= 'if=1, write coresource_nbicd cpo (coresource occurrence 6)';
        error_flag.section.coresource_nbicd      ='UAL';

    	error_flag.valeur.coresource_fusion    	= 1;
    	error_flag.type.coresource_fusion      	= 'integer';
    	error_flag.borne.coresource_fusion     	= {0,1}; 
    	error_flag.defaut.coresource_fusion   	= 1;
    	error_flag.info.coresource_fusion      	= 'if=1, write coresource_fusion cpo (coresource occurrence 7)';
        error_flag.section.coresource_fusion    ='UAL';

    	error_flag.valeur.coreneutrals    = 1;
    	error_flag.type.coreneutrals      = 'integer';
    	error_flag.borne.coreneutrals     = {0,1}; 
    	error_flag.defaut.coreneutrals    = 1;
    	error_flag.info.coreneutrals      = 'if=1, write coresource_neutral cpo (coreneutrals + coresource occurrence 8)';
        error_flag.section.coreneutrals    ='UAL';

    	error_flag.valeur.coresource_radiation    = 1;
    	error_flag.type.coresource_radiation      = 'integer';
    	error_flag.borne.coresource_radiation     = {0,1}; 
    	error_flag.defaut.coresource_radiation    = 1;
    	error_flag.info.coresource_radiation      = 'if=1, write coresource_radiation cpo (coresource occurrence 3)';
        error_flag.section.coresource_radiation    ='UAL';

    	error_flag.valeur.coresource_cyclotron    = 1;
    	error_flag.type.coresource_cyclotron      = 'integer';
    	error_flag.borne.coresource_cyclotron     = {0,1}; 
    	error_flag.defaut.coresource_cyclotron    = 1;
    	error_flag.info.coresource_cyclotron      = 'if=1, write coresource_cyclotron cpo (coresource occurrence 4)';
        error_flag.section.coresource_cyclotron   ='UAL';

     	error_flag.valeur.neoclassic   		  = 1;
    	error_flag.type.neoclassic      	  = 'integer';
    	error_flag.borne.neoclassic     	  = {0,1}; 
    	error_flag.defaut.neoclassic    	  = 1;
    	error_flag.info.neoclassic            = 'if=1, write neoclassic cpo (neoclassic + coresource occurrence 9)';
        error_flag.section.neoclassic         ='UAL';

   	error_flag.valeur.coresource_full    	= 1;
    	error_flag.type.coresource_full      	= 'integer';
    	error_flag.borne.coresource_full     	= {0,1}; 
    	error_flag.defaut.coresource_full   	= 1;
    	error_flag.info.coresource_full      	= 'if=1, write coresource_full cpo (coresource occurrence 10)';
        error_flag.section.coresource_full      ='UAL';

    	error_flag.valeur.equilibrium    	= 1;
    	error_flag.type.equilibrium      	= 'integer';
    	error_flag.borne.equilibrium     	= {0,1}; 
    	error_flag.defaut.equilibrium   	= 1;
    	error_flag.info.equilibrium      	= 'if=1, write equilibrium cpo';
        error_flag.section.equilibrium      ='UAL';
        
    	error_flag.valeur.grid_equi    	= 0;
    	error_flag.type.grid_equi      	= 'integer';
    	error_flag.borne.grid_equi     	= {0,1}; 
    	error_flag.defaut.grid_equi   	= 0;
    	error_flag.info.grid_equi     	= 'mesh grid used for the equilibrum 2 D profile :\n  if = 0, 2D grid is in (psi,theta);\n  if = 1, 2D grid is a Delaunay grid';
        error_flag.section.grid_equi      ='UAL';        

   	error_flag.valeur.equi_extrap   = 1;
    	error_flag.type.equi_extrap    	= 'integer';
    	error_flag.borne.equi_extrap  	= {0,1}; 
    	error_flag.defaut.equi_extrap   = 1;
    	error_flag.info.equi_extrap     = 'method for extrapolation of Psi oustside the LCFS: if = 1, smooth fields but accuracy decrease inside LCFS; if = 0, surface current at LCFS but better accuracy inside LCFS';
        error_flag.section.equi_extrap      ='UAL';
        
	error_flag.valeur.Convex_LCFS  = 1;
	error_flag.type.Convex_LCFS    = 'integer';
	error_flag.borne.Convex_LCFS   = {0,1}; 
	error_flag.defaut.Convex_LCFS  = 1;
	error_flag.info.Convex_LCFS    = 'Force LCFS used for 2D extrapolation to be convex:\n if = 0, keep LCFS as it is provided;\nif = 1, force LCFS to be convex';
	error_flag.section.Convex_LCFS = 'UAL';

	error_flag.valeur.fixed_grid  = 0;
	error_flag.type.fixed_grid    = 'integer';
	error_flag.borne.fixed_grid   = {0,1}; 
	error_flag.defaut.fixed_grid  = 0;
	error_flag.info.fixed_grid    = 'For retangular grid, if = 1 uses same grid for all time slices, otherwise uses floating grid following plasma displacement';
	error_flag.section.fixed_grid = 'UAL';
	
	error_flag.valeur.nb_points_pol  = 65;
	error_flag.type.nb_points_pol    = 'integer';
	error_flag.borne.nb_points_pol   = [35,255]; 
	error_flag.defaut.nb_points_pol  = 65;
	error_flag.info.nb_points_pol    = 'number of points in poloidal direction for inverse (rho,theta) grid of equilibrium';
	error_flag.section.nb_points_pol = 'UAL';
	
	error_flag.valeur.nb_points_radial  = 51;
	error_flag.type.nb_points_radial    = 'integer';
	error_flag.borne.nb_points_radial   = [33,301]; 
	error_flag.defaut.nb_points_radial  = 51;
	error_flag.info.nb_points_radial    = 'number of points in radial direction for rectangular (R,Z) grid of equilibrium';
	error_flag.section.nb_points_radial = 'UAL';
	
        % occurence string of the  output cpos
    	error_flag.valeur.scenario_occurrence    	= '';
    	error_flag.type.scenario_occurrence      	= 'string';
    	error_flag.borne.scenario_occurrence     	= ''; 
    	error_flag.defaut.scenario_occurrence   	= '';
    	error_flag.info.scenario_occurrence      	= ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.scenario_occurrence      ='Occurrence UAL';

   	error_flag.valeur.coreprof_occurrence    	= '';
    	error_flag.type.coreprof_occurrence      	= 'string';
    	error_flag.borne.coreprof_occurrence     	= ''; 
    	error_flag.defaut.coreprof_occurrence   	= '';
    	error_flag.info.coreprof_occurrence      	= ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.coreprof_occurrence      ='Occurrence UAL';

    	error_flag.valeur.coretransp_occurrence    = '';
    	error_flag.type.coretransp_occurrence      = 'string';
    	error_flag.borne.coretransp_occurrence     = ''; 
    	error_flag.defaut.coretransp_occurrence    = '';
    	error_flag.info.coretransp_occurrence      = ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.coretransp_occurrence   ='Occurrence UAL';
 
    	error_flag.valeur.coreneutrals_occurrence    = '';
    	error_flag.type.coreneutrals_occurrence      = 'string';
    	error_flag.borne.coreneutrals_occurrence     = ''; 
    	error_flag.defaut.coreneutrals_occurrence    = '';
    	error_flag.info.coreneutrals_occurrence      = ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.coreneutrals_occurrence   ='Occurrence UAL';

    	error_flag.valeur.neoclassic_occurrence    = '';
    	error_flag.type.neoclassic_occurrence      = 'string';
    	error_flag.borne.neoclassic_occurrence     = ''; 
    	error_flag.defaut.neoclassic_occurrence    = '';
    	error_flag.info.neoclassic_occurrence      = ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.neoclassic_occurrence   = 'Occurrence UAL';

    	error_flag.valeur.equilibrium_occurrence    = '';
    	error_flag.type.equilibrium_occurrence      = 'string';
    	error_flag.borne.equilibrium_occurrence     = ''; 
    	error_flag.defaut.equilibrium_occurrence  	= '';
    	error_flag.info.equilibrium_occurrence      = ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.equilibrium_occurrence   ='Occurrence UAL';

     	error_flag.valeur.coresources_occurrence    = '';
    	error_flag.type.coresources_occurrence     = 'string';
    	error_flag.borne.coresources_occurrence     = ''; 
    	error_flag.defaut.coresources_occurrence    = '';
    	error_flag.info.coresources_occurrence      = ' output cpo occurrence; if empty, use default occurence (0)';
        error_flag.section.coresources_occurrence   ='Occurrence UAL';

    	error_flag.valeur.COCOS    = 13;
    	error_flag.type.COCOS      = 'integer';
    	error_flag.borne.COCOS     = [1,18]; 
    	error_flag.defaut.COCOS    = 13;
    	error_flag.info.COCOS      = 'choice for the output COCOS';
        error_flag.section.COCOS    ='UAL';

%      	error_flag.valeur.TOKAMAK    = '';
%      	error_flag.type.TOKAMAK      = 'string';
%      	error_flag.borne.TOKAMAK     = ''; 
%      	error_flag.defaut.TOKAMAK    = '';
%      	error_flag.info.TOKAMAK      = 'Tokamak name to open another database than the default one';
%          error_flag.section.TOKAMAK   ='UAL';
%  
%      	error_flag.valeur.USER    = '';
%      	error_flag.type.USER      = 'string';
%      	error_flag.borne.USER     = ''; 
%      	error_flag.defaut.USER    = '';
%      	error_flag.info.USER      = 'User name to open another database than the default one';
%       error_flag.section.USER   = 'UAL';

	return
end 
if nargin < 2
	error_flag = - 9999;
	return
end

% changement de base de donnee
setappdata(0,'UAL_TOKAMAK','');
setappdata(0,'UAL_USER','');
setappdata(0,'UAL_DATAVERSION','');
% le parametre occurence peut etre une structure
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
	setappdata(0,'UAL_TOKAMAK',occurrence.tokamak);
	setappdata(0,'UAL_USER',occurrence.user);
	setappdata(0,'UAL_DATAVERSION',occurrence.dataversion);
	if isfield(occurrence,'occurrence')
		 occurrence = occurrence.occurrence;
	else
		 occurrence = '';
	end
	if ischar(occurrence)
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

% declaration des donnees
option 	= [];
cons1t 	= [];
geo1t 	= [];
sepa1t 	= [];
optin   = metis4itm(1);
z0dstruct.z0dinput.option = optin.valeur;

% cpo vide pour commencer
scenario = [];
coreprof = [];
coretransp = [];
coresources  = [];
coresource_lhcd  = [];
coresource_eccd  = [];
coresource_icrh  = [];
coresource_nbicd = [];
coresource_radiation = [];
coresource_cyclotron = [];
coreneutrals   = [];
coresource_neutral    = [];
coresource_fusion    = [];
coresource_ohm_boot    = [];
neoclassic  = [];
coresource_full  = [];
equilibrium = [];

% pointeur vers le cpo d'entree
scenario_cpo_in = [];

% error handling
% par d'erreur pour le debut
error_flag = 0;

% flag d'initialisation des CPOs de sortie
init_output_cpo = 0;

% 
% logfile
%
diary off
diary_file = tempname;
diary(diary_file);
disp('starting metis4itm')

% gestion des modes internes
% write_mode :si = 1, l'ecriture dans l'UAL est activee
write_mode = 1;
% read_input : si = 1, lecture du cpo scenario dans l'UAL
read_input = 0;       
% ressort apres la lecture des entrees pour METIS init
exit_after_read_input = 0;    
% read_output: si = 1, lecture des cpos en sortie depuis l'UAL
read_output = 0;
post = [];
if isstruct(methode)
	read_input = 0;
	write_mode = 0;
	post = methode;
	methode = 'post';
	if nargout < 2
		write_mode = 1;
	end
	
	% complete les options si necessaires
    	% option
	model = metis4itm(1);
    	noms = fieldnames(model.valeur);
    	for k = 1:length(noms)
        	if ~isfield(post.z0dinput.option,noms{k})
               		post.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
        	end
   	 end
         post.z0dinput.option = secure_option(post.z0dinput.option);
	
else
	switch methode
	case 'read'
		write_mode = 0;
		read_input = 0;	
		read_output = 1;
        case 'metis_from_ual'
		write_mode = 0;
		read_input = 1;	
		read_output = 0;
        exit_after_read_input = 1;    
	case {'full','fast','init','one_time'}
		read_input = 1;	
	end
	if nargout == 2 
		read_output = 1;
	end
end

% identification de la base
if (write_mode == 1) | (read_input == 1) | (read_output == 1)
	[treename,user,ver] = which_MDSdatabase_metis;
	if isempty(treename) | isempty(user) |isempty(ver)
		disp('No data base available')
		err_flag = -1001;
		return
	end
end


% open data base if need
if read_input == 1
	try
		disp('opening tree');
                if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
			expIdx =mem_cache( euitm_open_env('euitm',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
		else	
			expIdx = mem_cache(euitm_open('euitm',real(shot),real(run)));
		end
	catch
		% pas de donnees d'entree
		warning('no initial CPO available');
		expIdx = NaN;
	end
	mode_full = 0; 
	switch methode
	case {'fast','full'}
		mode_full = 1; 		
	end
	if isfinite(expIdx) & (expIdx >=0)
		% gestion occurrence
		if isempty(occurrence)
			ncpo = 'scenario';
		else
			ncpo = sprintf('scenario/%s',occurrence);    
		end
		% get scenario cpo 
		disp('reading CPO scenario');	
		if isempty(time)  | (mode_full == 1)
			try
				scenario_cpo_in =euitm_get(expIdx,ncpo);
			catch
				scenario_cpo_in = [];	
			end
		else
			try
				scenario_cpo_in = euitm_get_slice(expIdx,ncpo,time,interpolation_methode);
			catch
				scenario_cpo_in = [];	
			end
		end
		
		% Close the currently open database
		euitm_close(expIdx);
	else
		scenario_cpo_in = [];	
	end

end

% gestion des donnees d'entree
if ~isempty(scenario_cpo_in)
	% selon le nombre de temps
	
	if length(scenario_cpo_in) == 1
		if isempty(methode)
			methode = 'one_time';
		end	
	end
	switch methode
	case {'init','one_time'}

		% transforme le CPO en structure cronos	
		disp('parsing input data from CPO scenario to METIS for one time');
		scenario_cpo_in(end).codeparam = scenario_cpo_in(1).codeparam;
		scenario = swaptime2cronos(scenario_cpo_in(end));
		% mode evolution
		[option,time_cpo,cons1t,geo1t,sepa1t] = cpo2metis1t(scenario);
		if isempty(option)
			iof    = metis4itm(1);
			option = iof.valeur;
		end
		if isempty(time)
			time = 	time_cpo;
		end
        option = secure_option(option);
	otherwise

		% transforme le CPO en structure cronos	
		disp('parsing input data from CPO scenario to METIS');	
		scenario = swaptime2cronos(scenario_cpo_in);

		% mode simulation complete
		z0dinput = cpo2metis_input(scenario);
                z0dinput.option = secure_option(z0dinput.option);
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
    % codeparam_filename est une structure matlab
    info = codeparam_filename;
    noms = fieldnames(info);
    for k=1:length(noms)
        if isfield(option,noms{k})
            option.(noms{k}) = info.(noms{k});
        end
    end
    option = secure_option(option);
    z0dstruct.z0dinput.option = option;
    z0dinput.option = option;
    fprintf('METIS4ITM using parameters from input :\n')
    option
    
elseif ~isempty(codeparam_filename)
    % codeparam_filename designe un fichier xml
	% lecture du fichier
	fid = fopen(codeparam_filename,'r');
	if fid > 0
		codeparam = char(fread(fid,Inf,'char')');
		fclose(fid);
		% add codeparam to options
		info = metis4itm(1);
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
                option = secure_option(option);
		z0dstruct.z0dinput.option = option;
		z0dinput.option = option;
		fprintf('METIS4ITM  using parameters from %s :\n',codeparam_filename)
		option
	elseif ~isempty(codeparam_filename)    
        
        % codeparam_filename est une chaine xml

		tnp = tempname;
		fid = fopen(tnp,'w');
		fprintf(fid,'%s\n',codeparam_filename);
		fclose(fid);
		% add codeparam to options
		info = metis4itm(1);
		option = info.valeur;
		info  = xml_read(tnp);
		noms = fieldnames(info);
		for k=1:length(noms)
			if isfield(option,noms{k})
				option.(noms{k}) = info.(noms{k});
			end
		end
                option = secure_option(option);
		z0dstruct.z0dinput.option = option;
		z0dinput.option = option;
		fprintf('METIS4ITM using parameters from input :\n')
		option
		delete(tnp);
	else
		disp(sprintf('unable to read codeparam file %s',codeparam_filename));
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
fprintf('METIS4ITM called with method %s\n',methode);
switch methode
case 'auto_test'
        disp('METIS4ITM: starting suite of tests')
	error_flag = 1;
	try
		metis4itmautotest(real(shot),real(run));
		error_flag = 0;
		diary off
		delete(diary_file);
	catch
		diary off
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
		disp('test launched')
		z0dstruct  = test_metis_itm;
		data_zerod = z0dstruct.zerod;
		profil0d   = z0dstruct.profil0d;
		z0dinput   = z0dstruct.z0dinput;
		% compatibilite ascendante
		z0dstruct.zs = z0dstruct.zerod;

	case 'full'
		% cette methode calcul la simulation complete en mode fast : le CPO d'entree doit contenir au moins 3 temps
		disp('full computation')
		[z0dstruct.zerod,void,z0dstruct.profil0d] = zerod(z0dstruct.z0dinput.option,z0dstruct.z0dinput.cons, ...		
									z0dstruct.z0dinput.geo,z0dstruct.z0dinput.exp0d);
		data_zerod = z0dstruct.zerod;
		profil0d   = z0dstruct.profil0d;
		z0dinput   = z0dstruct.z0dinput;
		% compatibilite ascendante
		z0dstruct.zs = z0dstruct.zerod;

	case 'fast'
		% cette methode calcul la simulation complete en mode fast : le CPO d'entree doit contenir au moins 5 temps
		disp('fast computation')
		[z0dstruct.zerod,void,z0dstruct.profil0d] = zerodfast(z0dstruct.z0dinput.option,z0dstruct.z0dinput.cons,...
									z0dstruct.z0dinput.geo,z0dstruct.z0dinput.exp0d);
		data_zerod = z0dstruct.zerod;
		profil0d   = z0dstruct.profil0d;
		z0dinput   = z0dstruct.z0dinput;
		% compatibilite ascendante
		z0dstruct.zs = z0dstruct.zerod;

	case 'init'
		% cette methode precede one-time, elle doit etre appelee juste avant le calcul pour le premier intervalle de temps
		disp('init METIS for evolution computation')
		if isempty(time)
			disp('the time must be given with method init')
			error_flag = 4001;
			return
		end	
		if ~isempty(sepa1t)
			[data_zerod,profil0d,z0dstruct] = zerodevolution([],option,time,cons1t,geo1t,[],sepa1t);
		else
			[data_zerod,profil0d,z0dstruct] = zerodevolution([],option,time,cons1t,geo1t,[]);
		end
		% compatibilite ascendante
		z0dstruct.zerod = z0dstruct.zs;
        
        % initialisation des CPO de sortie
        if option.init_output_cpo
            init_output_cpo = 1;
        end
	case 'one_time'
		if isempty(time)
			disp('the time must be given with method one_time')
			error_flag = 4002;
			return
		end	
		if isappdata(0,'ITM_Z0DSTRUCT')
			disp('one time step METIS computation (using internal memorized data)')
			z0dstruct = getappdata(0,'ITM_Z0DSTRUCT');
            % restart ?
            if ~isempty(option.restart)
                save_restart_file = 1;
            end    

		else
			% equivalent a init
			disp('init METIS for evolution computation')
			z0dstruct = [];
            % initialisation des CPO de sortie
            if option.init_output_cpo
                init_output_cpo = 1;
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
                
	case 'post'
		fprintf('post processing of METIS simulation  data\n');
		z0dstruct = post;
		% compatibilite ascendante
		z0dstruct.zs = z0dstruct.zerod;
		z0dstruct.profil = z0dstruct.profil0d;
		%
		data_zerod = z0dstruct.zerod;
		profil0d   = z0dstruct.profil0d;
		z0dinput   = z0dstruct.z0dinput;
	case 'read'
		% rien dans ce cas	
	otherwise
		% we assume that the string method is a filename with complete path
		% this is the load method 
		fprintf('loading METIS simulation %s\n',methode)
		z0dstruct = load_metis_itm(methode);
		% add codeparam to options
		info = metis4itm(1);
		optf = info.valeur;
		noms = fieldnames(optf);
		for k=1:length(noms)
			if ~isfield(z0dstruct.z0dinput.option,noms{k})
				z0dstruct.z0dinput.option.(noms{k}) = optf.(noms{k});
			end
		end
                z0dstruct.z0dinput.option = secure_option(z0dstruct.z0dinput.option);
		% compatibilite ascendante
		z0dstruct.zs = z0dstruct.zerod;
		z0dstruct.profil = z0dstruct.profil0d;
		%
		data_zerod = z0dstruct.zerod;
		profil0d   = z0dstruct.profil0d;
		z0dinput   = z0dstruct.z0dinput;

	end
	error_flag = 0;
	setappdata(0,'ITM_Z0DSTRUCT',z0dstruct);

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
sigma_B0_eff = 1;
if exist('profil0d','var') && exist('data_zerod','var')

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
	
	% mise en place signe Bphi
	if z0dstruct.z0dinput.option.signe < 0
		% on change le signe de b0
		data_zerod.qa           = - data_zerod.qa;
		data_zerod.q95          = - data_zerod.q95;
		data_zerod.qmin         = - data_zerod.qmin;
		data_zerod.qeff         = - data_zerod.qeff;
		data_zerod.q0           = - data_zerod.q0;
		data_zerod.phiplasma    = - data_zerod.phiplasma;
		profil0d.qjli           = - profil0d.qjli;
		profil0d.fdia           = - profil0d.fdia;
		profil0d.phi            = - profil0d.phi;
		profil0d.dphidx         = - profil0d.dphidx;
	end


	% appel de la fonction qui regle le bon cocos
	[data_zerod,profil0d,sigma_B0_eff,sigma_bvac_r,factor_two_pi] = ...
 			makecocos(data_zerod,profil0d,z0dstruct,z0dstruct.z0dinput.option.COCOS,factor_two_pi);
        
        
    if z0dstruct.z0dinput.option.Sn_fraction > 0
        error('METIS ITM CPO interface is not compatible whith tin (Sn) in plasma composition (option.Sn_fraction should be 0)');
    end
end

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
delete(diary_file);


% restart file
if save_restart_file
    metis4itm_save_restart(option.restart,z0dstruct,texte_diary);
end

% access to UAL
if write_mode == 1
	try
               	if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
			expIdx =mem_cache( euitm_open_env('euitm',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
		else	
			expIdx = mem_cache(euitm_open('euitm',real(shot),real(run)));
		end
	catch
               	if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
			if imag(shot) ~= 0
				expIdx =mem_cache( euitm_create_env('euitm',real(shot),real(run),imag(shot),imag(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
			else
				expIdx =mem_cache( euitm_create_env('euitm',real(shot),real(run),real(shot),0,getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
			end
		else	
			if imag(shot) ~= 0
				expIdx = mem_cache(euitm_create('euitm',real(shot),real(run),imag(shot),imag(run)));
			else
				expIdx = mem_cache(euitm_create('euitm',real(shot),real(run),real(shot),0));
			end
		end
	end
	
	if ~(isfinite(expIdx) & (expIdx >=0))
		disp('Unable to open data base or to connect to UAL')
		error_flag = -9001;
		return
	end

	% mapping ITM structures
	ntime_scenario = length(data_zerod.temps);
	ntime_coreprof = length(profil0d.temps);
	%
	disp('allocating data structure scenario ')
	if ntime_scenario > 1
		scenario_cpo_out   = CPO_gen('scenario',ntime_scenario);
	else
		scenario_cpo_out   = CPO_gen('scenario');
	end
	%
	disp('mapping CPO data structure scenario')
	scenario = mapscenario(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,[],run,occurrence,sigma_B0_eff,sigma_bvac_r);
	% recopie des informations
	coreprof.datainfo = scenario.datainfo;
	coretransp.datainfo = scenario.datainfo;
  	coresources.datainfo = scenario.datainfo;
%  	coresource_lhcd.datainfo = scenario.datainfo;
%  	coresource_eccd.datainfo = scenario.datainfo;
%  	coresource_icrh.datainfo = scenario.datainfo;
%  	coresource_nbicd.datainfo = scenario.datainfo;
%  	coresource_full.datainfo = scenario.datainfo;
%  	coresource_radiation.datainfo = scenario.datainfo;
%  	coresource_cyclotron.datainfo = scenario.datainfo;
	coreneutrals.datainfo   = scenario.datainfo;
%  	coresource_neutral.datainfo    = scenario.datainfo;
%  	coresource_fusion.datainfo    = scenario.datainfo;
%  	coresource_ohm_boot.datainfo    = scenario.datainfo;
	neoclassic.datainfo  = scenario.datainfo;
	equilibrium.datainfo = scenario.datainfo;
	%
	disp('casting CPO data structure scenario from matlab to java')
	scenario_cpo_out = swaptime2ual(scenario_cpo_out,scenario,'scenario');
	disp('writing CPO data structure scenario')
	if isempty(strtrim(z0dstruct.z0dinput.option.scenario_occurrence))
		ncpo = 'scenario';
	else
		ncpo = sprintf('scenario/%s',strtrim(z0dstruct.z0dinput.option.scenario_occurrence));
    end
    if init_output_cpo
        euitm_put_non_timed(expIdx,ncpo,scenario_cpo_out(1));
    elseif  ntime_scenario > 1
		euitm_put(expIdx,ncpo,scenario_cpo_out);
	else
		euitm_put_slice(expIdx,ncpo,scenario_cpo_out(1));
	end
	%
	if z0dstruct.z0dinput.option.coreprof == 1 
		disp('allocating data structure coreprof')
		if isempty(strtrim(z0dstruct.z0dinput.option.coreprof_occurrence))
			ncpo = 'coreprof';
		else
			ncpo = sprintf('coreprof/%s',strtrim(z0dstruct.z0dinput.option.coreprof_occurrence));
		end
		if ntime_coreprof  > 1
			coreprof_cpo_out   = CPO_gen('coreprof',ntime_coreprof);
		else
			coreprof_cpo_out   = CPO_gen('coreprof');
		end
		disp('mapping CPO data structure coreprof')
		[coreprof,vtor,vpol] = mapcoreprof(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreprof,sigma_B0_eff);
		disp('casting CPO data structure coreprof from matlab to java')
		coreprof_cpo_out = swaptime2ual(coreprof_cpo_out,coreprof,'coreprof');
		disp('writing CPO data structure coreprof')		
		if init_output_cpo
            euitm_put_non_timed(expIdx,ncpo,coreprof_cpo_out(1));
        elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,coreprof_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,coreprof_cpo_out(1));
		end
	end
	if z0dstruct.z0dinput.option.coretransp == 1
		disp('allocating data structure coretransp')
		if isempty(strtrim(z0dstruct.z0dinput.option.coretransp_occurrence))
			ncpo = 'coretransp';
		else
			ncpo = sprintf('coretransp/%s',strtrim(z0dstruct.z0dinput.option.coretransp_occurrence));
		end
		if ntime_coreprof  > 1
			coretransp_cpo_out = CPO_gen('coretransp',ntime_coreprof);
		else
			coretransp_cpo_out = CPO_gen('coretransp');
		end
		disp('mapping CPO data structure coretransp')
		coretransp = mapcoretransp(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coretransp,sigma_B0_eff);
		disp('casting CPO data structure coretransp from matlab to java')
		coretransp_cpo_out = swaptime2ual(coretransp_cpo_out,coretransp,'coretransp');
		disp('writing CPO data structure coretransp')
		if init_output_cpo
            euitm_put_non_timed(expIdx,ncpo,coretransp_cpo_out(1));
        elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,coretransp_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,coretransp_cpo_out(1));
		end
	end
	disp('allocating data structure coresources')
	if ntime_coreprof  > 1
			coresources_cpo_out = CPO_gen('coresource',ntime_coreprof);
	else
			coresources_cpo_out = CPO_gen('coresource');
	end
	disp('mapping CPO data structure coresource (generic part)')
	coresources = mapcoresource_generic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coresources,sigma_B0_eff);
	count_sources = 1;

	if z0dstruct.z0dinput.option.coresource_lhcd == 1
		disp('mapping CPO data structure coresource for lhcd')
		coresources.values(count_sources) = mapcoresource_lhcd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_eccd == 1
		disp('mapping CPO data structure coresource for eccd')
		coresources.values(count_sources) = mapcoresource_eccd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_radiation == 1
		disp('mapping CPO data structure coresource for radiation')
		coresources.values(count_sources) = mapcoresource_radiation(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
 	        coresources.values(count_sources) = mapcoresource_brem(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_cyclotron == 1
		disp('mapping CPO data structure coresource for cyclotron radiation')
		coresources.values(count_sources) = mapcoresource_cyclotron(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_icrh == 1
		disp('mapping CPO data structure coresource ICRH')
		coresources.values(count_sources) = mapcoresource_icrh(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_nbicd == 1
		disp('mapping CPO data structure coresource NBI')
		coresources.values(count_sources) = mapcoresource_nbicd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coresource_fusion == 1
		disp('mapping CPO data structure coresource for fusion')
		coresources.values(count_sources) = mapcoresource_fusion(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.coreneutrals == 1

		disp('allocating data structure coreneutrals')
		if isempty(strtrim(z0dstruct.z0dinput.option.coreneutrals_occurrence))
			ncpo = 'coreneutrals';
		else
			ncpo = sprintf('coreneutrals/%s',strtrim(z0dstruct.z0dinput.option.coreneutrals_occurrence));
		end

		if ntime_coreprof  > 1
			coreneutrals_cpo_out = CPO_gen('coreneutrals',ntime_coreprof);
		else
			coreneutrals_cpo_out = CPO_gen('coreneutrals');
		end
		disp('mapping CPO data structure coreneutrals')
		coreneutrals = mapcoreneutrals(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreneutrals,sigma_B0_eff);
		disp('casting CPO data structure coreneutrals from matlab to java')
		coreneutrals_cpo_out = swaptime2ual(coreneutrals_cpo_out,coreneutrals,'coreneutrals');
		disp('writing CPO data structure coreneutrals')
		if init_output_cpo
            		euitm_put_non_timed(expIdx,ncpo,coreneutrals_cpo_out(1));
        	elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,coreneutrals_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,coreneutrals_cpo_out(1));
		end
		%
		disp('mapping CPO data structure coresource for neutrals')
		coresources.values(count_sources) = mapcoresource_neutral(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end
	if z0dstruct.z0dinput.option.neoclassic == 1
		disp('mapping CPO data structure coresource_ohm_boot')
		coresources.values(count_sources) = mapcoresource_ohm_boot(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
		%
		
	        disp('allocating data structure neoclassic ')
		if isempty(strtrim(z0dstruct.z0dinput.option.neoclassic_occurrence))
			ncpo = 'neoclassic';
		else
			ncpo = sprintf('neoclassic/%s',strtrim(z0dstruct.z0dinput.option.neoclassic_occurrence));
		end
		if ntime_coreprof  > 1
			neoclassic_cpo_out = CPO_gen('neoclassic',ntime_coreprof);
		else
			neoclassic_cpo_out = CPO_gen('neoclassic');
		end
		disp('mapping CPO data structure neoclassic')
		neoclassic = mapneoclassic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,neoclassic,vtor,vpol,sigma_B0_eff);
		disp('casting CPO data structure neoclassic from matlab to java')
		neoclassic_cpo_out = swaptime2ual(neoclassic_cpo_out,neoclassic,'neoclassic');
		disp('writing CPO data structure neoclassic ')
		if init_output_cpo
            		euitm_put_non_timed(expIdx,ncpo,neoclassic_cpo_out(1));
        	elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,neoclassic_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,neoclassic_cpo_out(1));
		end
	end
	
	if z0dstruct.z0dinput.option.coresource_full == 1
		disp('mapping CPO data structure coresource_full')
		coresources.values(count_sources) = mapcoresource_full(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end

        % extra source for runaways
	if z0dstruct.z0dinput.option.runaway > 0
		disp('mapping CPO data structure coresource for runaway electron current')
		coresources.values(count_sources) = mapcoresource_runaways(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 	count_sources = count_sources + 1;
	end

	% ecriture cpo coresource
	if count_sources > 1
		disp('casting CPO data structure coresources from matlab to java')
		coresources_cpo_out = swaptime2ual(coresources_cpo_out,coresources,'coresource');
		disp('writing CPO data structure coresources')
		if isempty(strtrim(z0dstruct.z0dinput.option.coresources_occurrence))
			ncpo = 'coresource';
		else
			ncpo = sprintf('coresource/%s',strtrim(z0dstruct.z0dinput.option.coresources_occurrence));
		end
		if init_output_cpo
            		euitm_put_non_timed(expIdx,ncpo,coresources_cpo_out(1));
        	elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,coresources_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,coresources_cpo_out(1));
		end
	end

	if z0dstruct.z0dinput.option.equilibrium == 1
		disp('allocating data structure equilibrium')
		if isempty(strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence))
			ncpo = 'equilibrium';
		else
			ncpo = sprintf('equilibrium/%s',strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence));
		end
		if ntime_coreprof  > 1
			equilibrium_cpo_out = CPO_gen('equilibrium',ntime_coreprof);
		else
			equilibrium_cpo_out = CPO_gen('equilibrium');
		end
		disp('mapping CPO data structure equilibrium')
		equilibrium = mapequilibrium(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,equilibrium,scenario, ...
                      z0dstruct.z0dinput.option.grid_equi,sigma_B0_eff,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi);
		disp('casting CPO data structure equilibrium from matlab to java')
		equilibrium_cpo_out = swaptime2ual(equilibrium_cpo_out,equilibrium,'equilibrium');
		disp('writing CPO data structure equilibrium')
		if init_output_cpo
            		euitm_put_non_timed(expIdx,ncpo,equilibrium_cpo_out(1));
        	elseif ntime_coreprof  > 1
			euitm_put(expIdx,ncpo,equilibrium_cpo_out);
		else
			euitm_put_slice(expIdx,ncpo,equilibrium_cpo_out(1));
		end
	end
	
	if getappdata(0,'euitm_enable_mem_cache') == 1
		euitm_flush_all(expIdx)
	end
	euitm_close(expIdx);
	disp('end of metis4itm writing');
    
else
	switch methode
	case 'post'
		% creation de la structure de sortie vers matlab
		disp('mapping CPO data structure scenario')
		xsd.scenario = mapscenario(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,[],run,occurrence,sigma_B0_eff,sigma_bvac_r);

		xsd.coreprof.datainfo = xsd.scenario.datainfo;
		xsd.coretransp.datainfo = xsd.scenario.datainfo;
 		xsd.coresources.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_lhcd.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_eccd.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_icrh.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_nbicd.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_full.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_radiation.datainfo = xsd.scenario.datainfo;
%  		xsd.coresource_cyclotron.datainfo = xsd.scenario.datainfo;
		xsd.coreneutrals.datainfo   = xsd.scenario.datainfo;
%  		xsd.coresource_neutral.datainfo    = xsd.scenario.datainfo;
%  		xsd.coresource_fusion.datainfo    = xsd.scenario.datainfo;
%  		xsd.coresource_ohm_boot.datainfo    = xsd.scenario.datainfo;
		xsd.neoclassic.datainfo  = xsd.scenario.datainfo;
		xsd.equilibrium.datainfo =xsd.scenario.datainfo;

		disp('mapping CPO data structure coresource (generic part)')
		coresources = mapcoresource_generic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coresources,sigma_B0_eff);
		count_sources = 1;

		if z0dstruct.z0dinput.option.coreprof == 1
			disp('mapping CPO data structure coreprof')
			[xsd.coreprof,vtor,vpol] = mapcoreprof(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,xsd.coreprof,sigma_B0_eff);
		end
		if z0dstruct.z0dinput.option.coretransp == 1
			xsd.coretransp = mapcoretransp(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,xsd.coretransp,sigma_B0_eff);
			disp('casting CPO data structure coretransp from matlab to java')
		end
		if z0dstruct.z0dinput.option.coresource_lhcd == 1
			disp('mapping CPO data structure coresource for lhcd')
			xsd.coresource.values(count_sources) = mapcoresource_lhcd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_eccd == 1
			disp('mapping CPO data structure coresource for eccd')
			xsd.coresource.values(count_sources) = mapcoresource_eccd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_radiation == 1
			disp('mapping CPO data structure coresource for radiation')
			xsd.coresource.values(count_sources) = mapcoresource_radiation(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
 	        	xsd.coresources.values(count_sources) = mapcoresource_brem(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_cyclotron == 1
			disp('mapping CPO data structure coresource for cyclotron')
			xsd.coresource.values(count_sources) = mapcoresource_cyclotron(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_icrh == 1
			disp('mapping CPO data structure coresource for icrh')
			xsd.coresource.values(count_sources) = mapcoresource_icrh(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_nbicd == 1
			disp('mapping CPO data structure coresource for nbicd')
			xsd.coresource.values(count_sources) = mapcoresource_nbicd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coresource_fusion == 1
			disp('mapping CPO data structure coresource for fusion')
			xsd.coresource.values(count_sources) = mapcoresource_fusion(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.coreneutrals == 1
			disp('mapping CPO data structure coreneutrals')
			xsd.coreneutrals = mapcoreneutrals(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,xsd.coreneutrals,sigma_B0_eff);
			%
			disp('mapping CPO data structure coresource for neutrals')
			xsd.coresource.values(count_sources) = mapcoresource_neutral(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.neoclassic == 1
			disp('mapping CPO data structure coresource for ohm & boot')
			xsd.coresource.values(count_sources) = mapcoresource_ohm_boot(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
			%
			disp('mapping CPO data structure neoclassic')
			xsd.neoclassic = mapneoclassic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,xsd.neoclassic,vtor,vpol,sigma_B0_eff);
		end
		if z0dstruct.z0dinput.option.coresource_full == 1
			disp('mapping CPO data structure coresource for sum of sources')
			xsd.coresource.values(count_sources) = mapcoresource_full(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff);
	 		count_sources = count_sources + 1;
		end
		if z0dstruct.z0dinput.option.equilibrium == 1
			disp('mapping CPO data structure equilibrium')
			xsd.equilibrium = mapequilibrium(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,xsd.equilibrium,xsd.scenario, ...
                                             z0dstruct.z0dinput.option.grid_equi,sigma_B0_eff,z0dstruct.z0dinput.option.equi_extrap,factor_two_pi);
		end

		
	end
end

if read_output == 1
	% varaible de sortie
	xsd = [];
	%
	try
                if isappdata(0,'UAL_TOKAMAK') && ~isempty(getappdata(0,'UAL_TOKAMAK'))
			expIdx =mem_cache( euitm_open_env('euitm',real(shot),real(run),getappdata(0,'UAL_USER'),getappdata(0,'UAL_TOKAMAK'),getappdata(0,'UAL_DATAVERSION')));
		else	
			expIdx = mem_cache(euitm_open('euitm',real(shot),real(run)));
		end
	end
	
	if ~(isfinite(expIdx) & (expIdx >=0))
		disp('Unable to open data base or to connect to UAL')
		error_flag = -9001;
		return
	end
%  	if isempty(occurrence)
%  		ncpo = 'scenario';
%  	else
%  		ncpo = sprintf('scenario/%s',occurrence);
%  	end
	if isempty(strtrim(z0dstruct.z0dinput.option.scenario_occurrence))
		ncpo = 'scenario';
	else
		ncpo = sprintf('scenario/%s',strtrim(z0dstruct.z0dinput.option.scenario_occurrence));
	end
	% get  cpo 
	disp('reading CPO data structure scenario')
	try
		scenario_cpo_read = euitm_get(expIdx,ncpo);
	catch	
		scenario_cpo_read = [];	
	
	end
	if ~isempty(scenario_cpo_read)
		disp('parsing CPO to matlab structure');	
		xsd.scenario = swaptime2cronos(scenario_cpo_read);
	else
		disp('unable to read scenario CPO');	
		xsd.scenario = [];
	end

	
	%
	if z0dstruct.z0dinput.option.coreprof == 1
		% get  cpo 
		disp('reading CPO data structure coreprof')
		if isempty(strtrim(z0dstruct.z0dinput.option.coreprof_occurrence))
			ncpo = 'coreprof';
		else
			ncpo = sprintf('coreprof/%s',strtrim(z0dstruct.z0dinput.option.coreprof_occurrence));
		end
		try
			coreprof_cpo_read =euitm_get(expIdx,ncpo);
		catch	
			coreprof_cpo_read = [];	
	
		end
		if ~isempty(coreprof_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.coreprof = swaptime2cronos(coreprof_cpo_read);
		else
			disp('unable to read coreprof CPO');	
			xsd.coreprof = [];
		end

	end
	if z0dstruct.z0dinput.option.coretransp == 1
		% get  cpo 
		disp('reading CPO data structure coretransp')
		if isempty(strtrim(z0dstruct.z0dinput.option.coretransp_occurrence))
			ncpo = 'coretransp';
		else
			ncpo = sprintf('coretransp/%s',strtrim(z0dstruct.z0dinput.option.coretransp_occurrence));
		end
		try
			coretransp_cpo_read =euitm_get(expIdx,ncpo);
		catch	
			coretransp_cpo_read = [];	
	
		end
		if ~isempty(coretransp_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.coretransp = swaptime2cronos(coretransp_cpo_read);
		else
			disp('unable to read coretransp CPO');	
			xsd.coretransp = [];
		end
	end
	if isempty(strtrim(z0dstruct.z0dinput.option.coresources_occurrence))
			ncpo = 'coresource';
	else
			ncpo = sprintf('coresource/%s',strtrim(z0dstruct.z0dinput.option.coresources_occurrence));
	end
	% get  cpo 
	disp('reading CPO data structure coresource')
	try
		coresource_cpo_read =euitm_get(expIdx,ncpo);
	catch	
		coresource_cpo_read = [];	
	
	end
	if ~isempty(coresource_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.coresource = swaptime2cronos(coresource_cpo_read);
	else
			disp('unable to read coresource CPO');	
			xsd.coresource = [];
	end
	if z0dstruct.z0dinput.option.coreneutrals == 1
		disp('reading CPO data structure coreneutrals')
		if isempty(strtrim(z0dstruct.z0dinput.option.coreneutrals_occurrence))
			ncpo = 'coreneutrals';
		else
			ncpo = sprintf('coreneutrals/%s',strtrim(z0dstruct.z0dinput.option.coreneutrals_occurrence));
		end
		try
			coreneutrals_cpo_read = euitm_get(expIdx,ncpo);
		catch	
			coreneutrals_cpo_read = [];	
	
		end
		if ~isempty(coreneutrals_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.coreneutrals = swaptime2cronos(coreneutrals_cpo_read);
		else
			disp('unable to read coresource CPO');	
			xsd.coreneutrals = [];
		end
	end
	if z0dstruct.z0dinput.option.neoclassic == 1
		%
		disp('reading CPO data structure neoclassic')
		if isempty(strtrim(z0dstruct.z0dinput.option.neoclassic_occurrence))
			ncpo = 'neoclassic';
		else
			ncpo = sprintf('neoclassic/%s',strtrim(z0dstruct.z0dinput.option.neoclassic_occurrence));
		end
		try
			neoclassic_cpo_read = euitm_get(expIdx,ncpo);
		catch	
			neoclassic_cpo_read = [];	
	
		end
		if ~isempty(neoclassic_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.neoclassic = swaptime2cronos(neoclassic_cpo_read);
		else
			disp('unable to read coresource CPO');	
			xsd.neoclassic = [];
		end
	end
	if z0dstruct.z0dinput.option.equilibrium == 1
		disp('reading CPO data structure equilibrium')
		if isempty(strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence))
			ncpo = 'equilibrium';
		else
			ncpo = sprintf('equilibrium/%s',strtrim(z0dstruct.z0dinput.option.equilibrium_occurrence));
		end
		try
			equilibrium_cpo_read = euitm_get(expIdx,ncpo);
		catch	
			equilibrium_cpo_read = [];	
	
		end
		if ~isempty(equilibrium_cpo_read)
			disp('parsing CPO to matlab structure');	
			xsd.equilibrium = swaptime2cronos(equilibrium_cpo_read);
		else
			disp('unable to read coresource CPO');	
			xsd.equilibrium = [];
		end
	end
	
	euitm_close(expIdx);
	disp('end of metis4itm reading');

end

% datainfo empty
function datainfo = datainfo_empty

datainfo.dataprovider 		= '';
datainfo.putdate      		= '';
datainfo.source       		= '';
datainfo.comment      		= '';
datainfo.isref        		= [];
datainfo.putinfo.putmethod 	= ''; 
datainfo.putinfo.putaccess 	= ''; 
datainfo.putinfo.putlocation 	= ''; 
datainfo.putinfo.rigths 	= ''; 
datainfo.whatref.user		= '';
datainfo.whatref.machine 	= ''; 
datainfo.whatref.shot 		= [];
datainfo.whatref.run 		= []; 
datainfo.whatref.occurrence    = [];

%
% convert data from metis to cpo scenario (mapping)
%
function scenario = mapscenario(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,scenario,run,occurrence,sigma_B0_eff,sigma_bvac_r)

% pour les donnees independantes du temps qui sont dependantes du temps dans le cpo
vt = ones(size(data_zerod.temps));


% calcul du fuelling
fuelling_data = metis_gaz(z0dstruct,data_zerod,profil0d);
if ~isfield(scenario,'datainfo')
	scenario.datainfo = datainfo_empty;      
end
if isempty(scenario.datainfo.dataprovider) 
	switch z0dstruct.z0dinput.mode_exp
	case -2
		scenario.datainfo.dataprovider = 'METIS evolution';
	case -1
		scenario.datainfo.dataprovider = 'METIS simulation from scratch';	
	case 0
		scenario.datainfo.dataprovider = 'CRONOS';
	case 1
		scenario.datainfo.dataprovider = 'Tore Supra';	
	case 2 
		scenario.datainfo.dataprovider = 'JET';		
	case 11
		scenario.datainfo.dataprovider = 'Tore Supra (preparation)';	
	otherwise
		scenario.datainfo.dataprovider = getenv('USER');	
	end
end
if isempty(scenario.datainfo.putdate)
	scenario.datainfo.putdate      = sprintf('%f',clock2julday);
end
if isempty(scenario.datainfo.source)
	scenario.datainfo.source       = 'METIS';
end
if isempty(scenario.datainfo.comment)
	scenario.datainfo.comment      = 'ITM implementation of METIS';
end
if isempty(scenario.datainfo.isref)
	scenario.datainfo.isref        = 0;
end
%
if isempty(scenario.datainfo.whatref.user)
	scenario.datainfo.whatref.user = scenario.datainfo.dataprovider;
end
if isempty(scenario.datainfo.whatref.machine)
	scenario.datainfo.whatref.machine =  z0dstruct.z0dinput.machine;
end
if isempty(scenario.datainfo.whatref.shot)
	scenario.datainfo.whatref.shot    =  real(z0dstruct.z0dinput.shot(1));
end
if isempty(scenario.datainfo.whatref.run)
	scenario.datainfo.whatref.run    =  real(run);
end
if isempty(scenario.datainfo.whatref.occurrence)
	scenario.datainfo.whatref.occurrence    =  fix(str2num(occurrence));
end
scenario.datainfo.cocos = z0dstruct.z0dinput.option.COCOS;
%
scenario.centre.te0.value    = data_zerod.te0;
scenario.centre.ti0.value    = interp1_itm(profil0d.temps,profil0d.tip(:,1),data_zerod.temps,'pchip','extrap');
scenario.centre.ne0.value    = data_zerod.ne0;
scenario.centre.ni0.value    = data_zerod.ni0;
scenario.centre.shift0.value = data_zerod.d0;
scenario.centre.psi0.value   =  interp1_itm(profil0d.temps,profil0d.psi(:,1),data_zerod.temps,'pchip','extrap');
scenario.centre.phi0.value   = zeros(size(data_zerod.temps));
scenario.centre.q0.value     = data_zerod.q0;
scenario.centre.Rmag.value   = interp1_itm(profil0d.temps,profil0d.Raxe(:,1),data_zerod.temps,'pchip','extrap');
if length(scenario.centre.Rmag.value) == 1
	scenario.centre.Zmag.value   = z0dstruct.z0dinput.geo.z0(end);
else
	scenario.centre.Zmag.value   = z0dstruct.z0dinput.geo.z0;
end
scenario.centre.vtor_0.value = interp1_itm(profil0d.temps,profil0d.vtor(:,1),data_zerod.temps,'pchip','extrap');

%
scenario.configs.config.value = 0 * vt;
switch z0dstruct.z0dinput.option.configuration
case 0
	scenario.configs.config.value = 2  * vt;
case 1
	scenario.configs.config.value = 4  * vt;
case 2
	scenario.configs.config.value = 2 + 5 .* data_zerod.modeh;
case 3
	scenario.configs.config.value =  4 + 3 .* data_zerod.modeh;
case 4
	scenario.configs.config.value = 7  * vt;
end
scenario.configs.ecrh_freq.value 		= [];
scenario.configs.ecrh_loc.value 		= mean(z0dstruct.z0dinput.cons.xece) ;
scenario.configs.ecrh_mode.value 		= []; 
scenario.configs.ecrh_tor_ang.value 	= (z0dstruct.z0dinput.option.sens .* (pi / 8) ) * vt; 
scenario.configs.ecrh_pol_ang.value 	= (z0dstruct.z0dinput.option.angle_ece ./180 .* pi) * vt;
scenario.configs.ecrh_harm.value 		= [];
scenario.configs.enbi.value 		= z0dstruct.z0dinput.option.einj * vt ;
scenario.configs.r_nbi.value 		= z0dstruct.z0dinput.option.rtang * vt ;
scenario.configs.grad_b_drift.value 	= 1;
scenario.configs.icrh_freq.value 		= z0dstruct.z0dinput.option.freq *vt ;
scenario.configs.icrh_phase.value         = [];

if z0dstruct.z0dinput.option.fwcd == 2	
	scenario.configs.icrh_scheme = 'FW';	
elseif z0dstruct.z0dinput.option.fwcd == -1
	scenario.configs.icrh_scheme = 'FW_CCD';	
elseif z0dstruct.z0dinput.option.fwcd == 1
	scenario.configs.icrh_scheme = 'FW_CD';	
else
	scenario.configs.icrh_scheme = sprintf('%s_min_%d',z0dstruct.z0dinput.option.mino,ceil(mean(data_zerod.harm)));	
end

scenario.configs.LH_freq.value 		= (z0dstruct.z0dinput.option.freqlh .* 1e9) * vt ;
scenario.configs.LH_npar.value 		= z0dstruct.z0dinput.option.npar0 *vt; 
scenario.configs.pellet_ang.value 	= [];
scenario.configs.pellet_v.value	        = [];
scenario.configs.pellet_nba.value 	= [];

% champs manquant    
switch z0dstruct.z0dinput.option.scaling
case 0
	scenario.configs.lmode_sc 	= 'ITERH-96P(th)';
	scenario.configs.hmode_sc 	= 'ITERH-98P(y,2)';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= '(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
case 1
	scenario.configs.lmode_sc 	= 'OH scaling law (G. Bracco and K. Thomsen, NF, 37 ,1997)';
	scenario.configs.hmode_sc 	= 'ITERH-98P(y,2)';
	scenario.configs.core_sc  	=  'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	='(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
case 2
	scenario.configs.lmode_sc 	= 'W from core scaling from McDonald NF 47 (2007) p147';
	scenario.configs.hmode_sc 	= 'W from core scaling  +W from  pedestal scaling from McDonald NF 47 (2007) p147';
	scenario.configs.core_sc  	= 'core scaling from McDonald NF 47 (2007) p147'; 
	scenario.configs.pedestal_sc  	= 'pedestal scaling from McDonald NF 47 (2007) p147'; 
case 3
	scenario.configs.lmode_sc 	= 'W form GS03 - W from pedestal scaling';
	scenario.configs.hmode_sc 	= 'GS03 (ITPA @AIEA 2004 from McDonald  PPCF 46 ,2004,A215- , #11';
	scenario.configs.core_sc  	= 'W form GS03 - W from pedestal scaling';
	scenario.configs.pedestal_sc  	= 'T. Hoang for ITPA 2004';
case 4
	scenario.configs.lmode_sc 	= 'fit of experimental measurement of Wdia';
	scenario.configs.hmode_sc 	= 'fit of experimental measurement of Wdia';
	scenario.configs.core_sc  	= 'fit of experimental measurement of Wdia - Wpedestal';
	scenario.configs.pedestal_sc  	= '(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2'
case 5
	scenario.configs.lmode_sc 	= 'ITERH-96P(th)';
	scenario.configs.hmode_sc 	= 'ITERH-EIV(y,2)  (Mc Donald, NF 47 (2007) 147-174)';
	scenario.configs.core_sc  	= 'difference between Wtot and Wpedestal';
	scenario.configs.pedestal_sc  	= 'W from ITERH-EIV(y,2) - W from ITERH-96P(th)';
case 6
	scenario.configs.lmode_sc 	= 'OH scaling law from Wesson Tokamak with transition to ITERH-96P(th) when additional heating is applied';
	scenario.configs.hmode_sc 	= 'ITERH-98P(y,2)';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	='(W from ITERH-98P(y,2) - W from ITERH-96P(th)) / 2';
case 7
	scenario.configs.lmode_sc 	= 'ITERH-98P(y,2)/2';
	scenario.configs.hmode_sc 	= 'ITERH-98P(y,2)';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= 'ITERH-98P(y,2)/2';
case 8
	scenario.configs.lmode_sc 	= 'user defined scaling';
	scenario.configs.hmode_sc 	= 'user defined scaling';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= 'user defined scaling';
case 9
	scenario.configs.lmode_sc 	= 'ITERH-96P(th)';
	scenario.configs.hmode_sc 	= 'J_pol = 0 (Hybrid or advance scenario)';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= 'difference between H mode and L mode';
case 10
	scenario.configs.lmode_sc 	= 'Sauter & Martin';
	scenario.configs.hmode_sc 	= 'quasi analytical scaling law';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= 'difference between H mode and L mode';
case 11
	scenario.configs.lmode_sc 	= 'Sauter & Martin';
	scenario.configs.hmode_sc 	= 'Elbeze EIV 2005';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal';
	scenario.configs.pedestal_sc  	= 'difference between H mode and L mode';
case 12
	scenario.configs.lmode_sc 	= 'Sauter & Martin';
	scenario.configs.hmode_sc 	= 'Elbeze EIV 2005';
	scenario.configs.core_sc  	= 'W from L-mode or H-mode - Wpedestal with limitation on beta_N';
	scenario.configs.pedestal_sc  	= 'difference between H mode and L mode';
end

if z0dstruct.z0dinput.option.tauhemul  == 0
	scenario.configs.helium_sc 	= 'ITER,physics basis, Nuclear Fusion 39, 1999, p 2226';
else
	scenario.configs.helium_sc 	= sprintf('%g * tau_E',z0dstruct.z0dinput.option.tauhemul);
end

switch z0dstruct.z0dinput.option.zeff
case 0
	scenario.configs.impurity_sc = 'None, line averaged Zeff given';
case 1
	scenario.configs.impurity_sc = 'None, line averaged Zeff given + impurities neoclassical accumulation';
case 2
	scenario.configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
case 3
	scenario.configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
case 4
	scenario.configs.impurity_sc = 'undefined scaling ; must be not used (reserved)';
case 5
	scenario.configs.impurity_sc = 'Tore Supra scaling law';
case 6
	scenario.configs.impurity_sc = 'Matthews scaling law';
otherwise
	scenario.configs.impurity_sc = 'other scaling law';
end

switch z0dstruct.z0dinput.option.l2hscaling   
case 0
	scenario.configs.l2h_sc		= 'LH99(1)';	
case 1
	scenario.configs.l2h_sc		= 'LH2002 without Zeff effect';	
case 2
	scenario.configs.l2h_sc		= 'LH2002 with Zeff effect';
case 3
	scenario.configs.l2h_sc		= 'YR Martin 2008';	
case 4
	scenario.configs.l2h_sc		= 'NLM-7 Murari 2012';	
case 5
	scenario.configs.l2h_sc		= 'NLM-11 Murari 2012';	
case 6
	scenario.configs.l2h_sc		= 'Jpol change of signe in edge region (E. R. Solano rule)';	
case 10
	scenario.configs.l2h_sc		= 'multimachine scaling law from Murami 2013 (BUEMS)';	
otherwise
	scenario.configs.l2h_sc		= 'criterion base on plasma rotation ( Gamma_ExB / Gamma_ITG)';	
end
if z0dstruct.z0dinput.option.taurotmul == 0
  scenario.configs.tor_rot_sc	= 'ion heat confinement time';
else
  scenario.configs.tor_rot_sc	= sprintf(' somptaneous rotation given by Rice NF  47, (2007) p 1618-1624 + NBI torque with %g * tau_E confinement time', ...
                                          		z0dstruct.z0dinput.option.taurotmul);
end
%
scenario.configs.coordinate	= 'sqrt(PSI_TOR/(pi * B0))';
% matieres 
%  scenario.configs.wall_mat 	= '';
%  scenario.configs.evap_mat 	= '';
%  scenario.configs.lim_mat 	= '';
%  scenario.configs.div_mat 	= '';

%
scenario.confinement.tau_e.value         = data_zerod.taue;
scenario.confinement.tau_l_sc.value      = data_zerod.tauthl;
scenario.confinement.tau_h_sc.value      = data_zerod.tauh;
scenario.confinement.tau_he.value        = data_zerod.tauh;
scenario.confinement.tau_e_ee.value      = data_zerod.tauee;
scenario.confinement.tau_e_ii.value      = data_zerod.tauii; 
scenario.confinement.tau_e_ei.value      = data_zerod.tauei;  
scenario.confinement.tau_cur_diff.value  = data_zerod.tauj; 
scenario.confinement.tau_i_rol.value     = data_zerod.tauip;
%
scenario.currents.RR.value  		= data_zerod.RR;
scenario.currents.i_align.value  	= data_zerod.ialign;
scenario.currents.i_boot.value  	= data_zerod.iboot; 
scenario.currents.i_cd_tot.value  	= data_zerod.icd; 
scenario.currents.i_eccd.value 		= data_zerod.ieccd;  
scenario.currents.i_fast_ion.value 	= data_zerod.ifus;
scenario.currents.i_fwcd.value 		= data_zerod.ifwcd; 
scenario.currents.i_lhcd.value 		= data_zerod.ilh;  
scenario.currents.i_nbicd.value 	= data_zerod.inbicd; 
scenario.currents.i_ni_tot.value 	= data_zerod.ini;  
scenario.currents.i_ohm.value 		= data_zerod.iohm;
scenario.currents.i_par.value 		= data_zerod.ipar;
scenario.currents.i_runaway.value 	= data_zerod.irun;
scenario.currents.v_loop.value 		= data_zerod.vloop; 
scenario.currents.v_meas.value 		= data_zerod.vmes; 
%
scenario.edge.te_edge.value 		= data_zerod.tebord; 
scenario.edge.ti_edge.value 		= data_zerod.tibord;  
scenario.edge.ne_edge.value 		= data_zerod.nebord; 
scenario.edge.ni_edge.value 		= data_zerod.nibord;
scenario.edge.psi_edge.value 		=  interp1_itm(profil0d.temps,profil0d.psi(:,end),data_zerod.temps,'pchip','extrap') ;
scenario.edge.phi_edge.value 		= interp1_itm(profil0d.temps,profil0d.phi(:,end),data_zerod.temps,'pchip','extrap'); 
scenario.edge.rho_edge.value 		= data_zerod.rm; 
scenario.edge.drho_edge_dt.value 	= data_zerod.drmdt;
scenario.edge.q_edge.value 		= data_zerod.qa;
scenario.edge.neutral_flux.value 	= data_zerod.n0a;
scenario.edge.phi_plasma.value 		= data_zerod.phiplasma;
scenario.edge.vtor_edge.value 		= interp1_itm(profil0d.temps,profil0d.vtor(:,end),data_zerod.temps,'pchip','extrap');
%
scenario.energy.w_tot.value 		= data_zerod.w;
scenario.energy.w_b_pol.value 		= data_zerod.wbp;
scenario.energy.w_dia.value 		= data_zerod.wdia;

if length(data_zerod.temps) == 1
	wdia  = z0dstruct.zerod.wdia;
        t     = z0dstruct.zerod.temps;
	scenario.energy.dwdia_dt.value 	= (wdia(end) - wdia(end-1)) ./ (t(end) - t(end -1));
else
	scenario.energy.dwdia_dt.value 	= z0dxdt(data_zerod.wdia,data_zerod.temps);
end
scenario.energy.w_b_tor_pla.value 	= data_zerod.wmagtor;
scenario.energy.w_th.value 		= data_zerod.wth;
scenario.energy.dwtot_dt.value 		= data_zerod.dwdt;
scenario.energy.dwbpol_dt.value 	= data_zerod.dwbpdt;
scenario.energy.dwbtorpla_dt.value 	= data_zerod.dwmagtordt;
scenario.energy.dwth_dt.value 		= data_zerod.dwthdt;
scenario.energy.esup_icrhtot.value 	= data_zerod.esup_icrh;
scenario.energy.esup_icrhper.value 	= ((1/3) + (1/3) .*  max(0,tanh(data_zerod.einj_icrh ./ data_zerod.ecrit_icrh))) .* ...
                                  		data_zerod.esup_icrh;
scenario.energy.esup_nbitot.value 	= data_zerod.esup_nbi;
scenario.energy.esup_nbiperp.value 	= ((1/3) + (2/3) .*  sin(z0dstruct.z0dinput.option.angle_nbi./180*pi) .*  ...
                                  		max(0,tanh(z0dstruct.z0dinput.option.einj ./ data_zerod.ecrit_nbi))) .* ...
                                  		data_zerod.esup_nbi;
scenario.energy.esup_lhcd.value 	= data_zerod.esup_lh;
scenario.energy.esup_alpha.value 	= data_zerod.esup_fus;
%
scenario.eqgeometry.source              = 'METIS';
%   for k=1:length(data_zerod.temps)
%   	if data_zerod.xpoint(k) == 1
%   		scenario.eqgeometry.boundarytype{k}  = 'separatrix';
%   	else 
%   		scenario.eqgeometry.boundarytype{k}  = 'limiter';
%   	end
%   end
%scenario.eqgeometry.boundarytype = sprintf('%f',data_zerod.xpoint);
scenario.eqgeometry.boundarytype = data_zerod.xpoint;

if isfield(profil0d,'Rsepa') && isfield(profil0d,'Zsepa')
	scenario.eqgeometry.boundary(1).r 		= interp1_itm(profil0d.temps,profil0d.Rsepa,data_zerod.temps,'nearest','extrap');
	scenario.eqgeometry.boundary(1).z 		= interp1_itm(profil0d.temps,profil0d.Zsepa,data_zerod.temps,'nearest','extrap');
else
	scenario.eqgeometry.boundary.r 		= [];
	scenario.eqgeometry.boundary.z 		= [];

end
if length(data_zerod.temps) == 1
	scenario.eqgeometry.geom_axis.r 		= z0dstruct.z0dinput.geo.R(end);
	scenario.eqgeometry.geom_axis.z 		= z0dstruct.z0dinput.geo.z0(end);
	scenario.eqgeometry.a_minor 	 		= z0dstruct.z0dinput.geo.a(end);
	scenario.eqgeometry.elongation 	 		= z0dstruct.z0dinput.geo.K(end);
	scenario.eqgeometry.tria_upper 	 		= z0dstruct.z0dinput.geo.d(end);
	scenario.eqgeometry.tria_lower 	 		= z0dstruct.z0dinput.geo.d(end);
else
	scenario.eqgeometry.geom_axis.r 		= z0dstruct.z0dinput.geo.R;
	scenario.eqgeometry.geom_axis.z 		= z0dstruct.z0dinput.geo.z0;
	scenario.eqgeometry.a_minor 	 		= z0dstruct.z0dinput.geo.a;
	scenario.eqgeometry.elongation 	 		= z0dstruct.z0dinput.geo.K;
	scenario.eqgeometry.tria_upper 	 		= z0dstruct.z0dinput.geo.d;
	scenario.eqgeometry.tria_lower 	 		= z0dstruct.z0dinput.geo.d;
end
scenario.eqgeometry.xpts.r 			= [];
scenario.eqgeometry.xpts.z 			= [];
scenario.eqgeometry.left_low_st.r 		= [];
scenario.eqgeometry.left_low_st.z 		= [];
scenario.eqgeometry.right_low_st.r 		= [];
scenario.eqgeometry.right_low_st.z 		= [];
scenario.eqgeometry.left_up_st.r 		= []; 
scenario.eqgeometry.left_up_st.z 		= []; 
scenario.eqgeometry.right_up_st.r 		= [];
scenario.eqgeometry.right_up_st.z 		= [];
scenario.eqgeometry.active_limit.r 		= [];  
scenario.eqgeometry.active_limit.z 		= [];
%
scenario.global_param.ip.value		= data_zerod.ip;
if length(data_zerod.temps) == 1
 	ip    = z0dstruct.zerod.ip;
    	t     = z0dstruct.zerod.temps;
 	scenario.global_param.dip_dt.value 	= (ip(end) - ip(end-1)) ./ (t(end) - t(end -1));
else
 	scenario.global_param.dip_dt.value 		= z0dxdt(data_zerod.ip,data_zerod.temps); 
end
scenario.global_param.beta_pol.value 		= data_zerod.betaptot;
mu0 =  4 .* pi .* 1e-7;
if length(data_zerod.temps) == 1
	scenario.global_param.beta_tor.value 		= (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0(end) .^ 2 ./2 ./ mu0); 
else
	scenario.global_param.beta_tor.value 		= (2/3) .* (data_zerod.w ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0 .^ 2 ./2 ./ mu0); 
end	
scenario.global_param.beta_normal.value 	= data_zerod.betan;
scenario.global_param.li.value 		= data_zerod.li; 
scenario.global_param.volume.value 		= data_zerod.vp;
scenario.global_param.area_pol.value 	= data_zerod.sp;
scenario.global_param.area_ext.value 	= data_zerod.sext; 
scenario.global_param.len_sepa.value 	= data_zerod.peri; 
scenario.global_param.beta_pol_th.value 	= data_zerod.betap;
if length(data_zerod.temps) == 1
	scenario.global_param.beta_tor_th.value	=  (2/3) .* (data_zerod.wth ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0(end) .^ 2 ./2 ./ mu0); 
else
	scenario.global_param.beta_tor_th.value	=  (2/3) .* (data_zerod.wth ./ data_zerod.vp) ./ (z0dstruct.z0dinput.geo.b0 .^ 2 ./2 ./ mu0); 
end
scenario.global_param.beta_n_th.value 	= scenario.energy.w_th.value ./ max(eps,scenario.energy.w_tot.value ) .* data_zerod.betan;
scenario.global_param.disruption.value 	= data_zerod.disrup;
scenario.global_param.mode_h.value 		= data_zerod.modeh;
scenario.global_param.s_alpha.value 		= data_zerod.salpha;

%
scenario.heat_power.plh.value 			= data_zerod.plh;
scenario.heat_power.pohmic.value 		= data_zerod.pohm;
scenario.heat_power.picrh.value 		= data_zerod.picrh;
scenario.heat_power.pecrh.value 		= data_zerod.pecrh;
scenario.heat_power.pnbi.value 			= data_zerod.pnbi;
snbi  = sin(z0dstruct.z0dinput.option.angle_nbi ./ 180 * pi);
scenario.heat_power.pnbi_co_cur.value 	= snbi .* (snbi > 0) .* data_zerod.pnbi;
scenario.heat_power.pnbi_counter.value 	= abs(snbi) .* (snbi < 0) .* data_zerod.pnbi;
scenario.heat_power.plh_th.value 		= data_zerod.plh_th;
scenario.heat_power.picrh_th.value 		= data_zerod.picrh_th;
scenario.heat_power.pecrh_th.value 		= data_zerod.pecrh;    % this is not a mistake
scenario.heat_power.pnbi_th.value 		= data_zerod.pnbi_th;
scenario.heat_power.ploss_icrh.value 		= data_zerod.frloss_icrh .* data_zerod.picrh;
scenario.heat_power.ploss_nbi.value 		= data_zerod.frnbi .* data_zerod.pnbi;
scenario.heat_power.pbrem.value		= data_zerod.pbrem;
scenario.heat_power.pcyclo.value 		= data_zerod.pcyclo; 
scenario.heat_power.prad.value 		= data_zerod.prad; 
scenario.heat_power.pdd_fus.value 		= data_zerod.pddfus;
scenario.heat_power.pei .value		= data_zerod.pei;
scenario.heat_power.pel_tot.value 		= data_zerod.pel;
scenario.heat_power.pel_fus.value 		= data_zerod.pel_fus;
scenario.heat_power.pel_icrh.value 		= data_zerod.pel_icrh;
scenario.heat_power.pel_nbi.value 		= data_zerod.pel_nbi;    % activer a partir de la version 4.07 de l'UAL
scenario.heat_power.pfus_dt.value 		= data_zerod.pfus;
scenario.heat_power.ploss_fus.value 		= data_zerod.pfus_loss;
scenario.heat_power.pfus_nbi.value 		= data_zerod.pfus_nbi;
scenario.heat_power.pfus_th.value 		= data_zerod.pfus_th; 
scenario.heat_power.padd_tot.value 		= data_zerod.pin; 
scenario.heat_power.pion_tot.value 		= data_zerod.pion; 
scenario.heat_power.pion_fus.value 		= data_zerod.pion_fus; 
scenario.heat_power.pion_icrh.value 		= data_zerod.pion_icrh; 
scenario.heat_power.pion_nbi.value 		= data_zerod.pion_nbi;
scenario.heat_power.pioniz.value 		= data_zerod.pioniz; 
scenario.heat_power.ploss.value 		= data_zerod.ploss; 
scenario.heat_power.p_wth.value 		= data_zerod.pth; 
scenario.heat_power.p_w.value		= data_zerod.pw; 
scenario.heat_power.p_l2h_thr.value 		= data_zerod.plhthr;
scenario.heat_power.p_l2h_sc.value 		= data_zerod.plossl2h;
scenario.heat_power.p_nbi_icrh.value 		= data_zerod.pnbi_icrh;
%
scenario.itb.q_min.value 		= data_zerod.qmin;
[xx,tt] = meshgrid(profil0d.xli,profil0d.temps);
scenario.itb.te_itb.value		= griddata_itm(tt,xx,profil0d.tep,data_zerod.temps,data_zerod.xitb,'cubic');
scenario.itb.ti_itb.value 		= griddata_itm(tt,xx,profil0d.tip,data_zerod.temps,data_zerod.xitb,'cubic'); 
scenario.itb.ne_itb.value 		= griddata_itm(tt,xx,profil0d.nep,data_zerod.temps,data_zerod.xitb,'cubic'); 
scenario.itb.ni_itb.value 		= griddata_itm(tt,xx,profil0d.nip,data_zerod.temps,data_zerod.xitb,'cubic'); 
scenario.itb.psi_itb.value 		= griddata_itm(tt,xx,profil0d.psi,data_zerod.temps,data_zerod.xitb,'cubic');  
scenario.itb.phi_itb.value 		= griddata_itm(tt,xx,profil0d.phi,data_zerod.temps,data_zerod.xitb,'cubic');
scenario.itb.rho_itb .value		= griddata_itm(tt,xx,profil0d.rmx,data_zerod.temps,data_zerod.xitb,'cubic');
scenario.itb.h_itb .value		= data_zerod.hitb;
width = interp1_itm(profil0d.temps,mean(double(profil0d.xieshape_itb < profil0d.xieshape),2),data_zerod.temps,'pchip','extrap') .* data_zerod.rm;
scenario.itb.width_itb.value 		= width;
scenario.itb.vtor_itb.value 		= griddata_itm(tt,xx,profil0d.vtor,data_zerod.temps,data_zerod.xitb,'cubic');
scenario.itb.itb_type.value 		= ((data_zerod.hitb .* data_zerod.hmhd) > 0) .* 8;  % a completer ulterieument
%
scenario.lim_div_wall.te_lim_div.value 	= data_zerod.telim;
scenario.lim_div_wall.ti_lim_div.value 	= [];
scenario.lim_div_wall.ne_lim_div.value 	= data_zerod.nelim;
scenario.lim_div_wall.ni_lim_div.value 	= []; 
scenario.lim_div_wall.p_peak_div.value 	= data_zerod.peakdiv;
scenario.lim_div_wall.surf_temp.value 	= []; 
scenario.lim_div_wall.p_lim_div.value 	= data_zerod.plim;
scenario.lim_div_wall.p_rad_div.value 	= [];
scenario.lim_div_wall.wall_temp.value 	= []; 
scenario.lim_div_wall.wall_state.value 	= [];
scenario.lim_div_wall.detach_state.value 	= double(((data_zerod.prad + data_zerod.pbrem + data_zerod.pcyclo + data_zerod.pioniz + data_zerod.pradsol) >= data_zerod.ploss) | (data_zerod.telim < 1));

% compatibilite de version
%  if any(getappdata(0,'UALVERSION')  > '4.08b')
scenario.lim_div_wall.pump_flux.value 	= fuelling_data.pumping;
%  else
%     scenario.lim_div_wall.pump_flux 	= fuelling_data.pumping;
%end

%
scenario.line_ave.ne_line.value 		= data_zerod.nbar;
scenario.line_ave.zeff_line.value 		= data_zerod.zmszl ./ max(eps,data_zerod.zeff);
zeffne = trapz(profil0d.xli,profil0d.nep .* profil0d.zeff,2) ./ max(eps,trapz(profil0d.xli,profil0d.nep,2));
scenario.line_ave.ne_zeff_line.value 		= interp1_itm(profil0d.temps,zeffne,data_zerod.temps,'pchip','extrap'); 
if length(data_zerod.temps) == 1
	nbar = z0dstruct.zerod.nbar;
        t     = z0dstruct.zerod.temps;
	scenario.line_ave.dne_line_dt.value 	= (nbar(end) - nbar(end-1)) ./ (t(end) - t(end -1));
else
	scenario.line_ave.dne_line_dt.value 	= z0dxdt(data_zerod.nbar,data_zerod.temps);
end
%
scenario.neutron.ndd_tot.value 		= data_zerod.ndd;
scenario.neutron.ndd_th.value		= data_zerod.ndd_th; 
scenario.neutron.ndd_nbi_th.value	= data_zerod.ndd_nbi_th; 
scenario.neutron.ndd_nbi_nbi.value 	= data_zerod.ndd_nbi_nbi;
scenario.neutron.ndt_tot.value 		= data_zerod.pfus ./(3.56e6 .* 1.602176462e-19);
scenario.neutron.ndt_th.value 		= max(0,data_zerod.pfus - data_zerod.pfus_nbi) ./ (3.56e6 .* 1.602176462e-19); 
%
tt = profil0d.temps * ones(size(profil0d.xli));
psin = (profil0d.psi - profil0d.psi(:,end) *  ones(size(profil0d.xli))) ./ ((profil0d.psi(:,1) - profil0d.psi(:,end)) *  ones(size(profil0d.xli)));
v95  = 0.95 .* ones(size(data_zerod.temps));
scenario.ninety_five.q_95.value 		= griddata_itm(tt,psin,profil0d.qjli,data_zerod.temps,v95,'cubic');
scenario.ninety_five.elong_95.value 		= griddata_itm(tt,psin,profil0d.kx,data_zerod.temps,v95,'cubic');
scenario.ninety_five.tria_95.value 		= griddata_itm(tt,psin,profil0d.dx,data_zerod.temps,v95,'cubic');
scenario.ninety_five.tria_up_95.value 	= [];
scenario.ninety_five.tria_lo_95.value 	= []; 
scenario.ninety_five.te_95.value 		= griddata_itm(tt,psin,profil0d.tep,data_zerod.temps,v95,'cubic'); 
scenario.ninety_five.ti_95.value 		= griddata_itm(tt,psin,profil0d.tip,data_zerod.temps,v95,'cubic'); 
scenario.ninety_five.ne_95.value 		= griddata_itm(tt,psin,profil0d.nep,data_zerod.temps,v95,'cubic'); 
scenario.ninety_five.ni_95.value 		= griddata_itm(tt,psin,profil0d.nip,data_zerod.temps,v95,'cubic'); 
scenario.ninety_five.phi_95.value 		= griddata_itm(tt,psin,profil0d.phi,data_zerod.temps,v95,'cubic');   
scenario.ninety_five.rho_95.value 		= griddata_itm(tt,psin,profil0d.rmx,data_zerod.temps,v95,'cubic');  
scenario.ninety_five.vtor_95.value 		= griddata_itm(tt,psin,profil0d.vtor,data_zerod.temps,v95,'cubic');
%
scenario.pedestal.te_ped.value 		= data_zerod.teped;
scenario.pedestal.ti_ped.value 		= data_zerod.tiped;
scenario.pedestal.ne_ped.value 		= data_zerod.neped; 
scenario.pedestal.ni_ped.value 		= data_zerod.niped; 
scenario.pedestal.psi_ped.value 	= data_zerod.modeh .* interp1_itm(profil0d.temps,profil0d.psi(:,end-1),data_zerod.temps,'pchip','extrap') + ...
					  (~data_zerod.modeh) .* scenario.edge.psi_edge.value; 
scenario.pedestal.phi_ped.value 		= data_zerod.modeh .* interp1_itm(profil0d.temps,profil0d.phi(:,end-1),data_zerod.temps,'pchip','extrap') + ...
					  (~data_zerod.modeh) .* scenario.edge.phi_edge.value; 
scenario.pedestal.rho_ped.value 		= data_zerod.modeh .* interp1_itm(profil0d.temps,profil0d.rmx(:,end-1),data_zerod.temps,'pchip','extrap') + ...
					  (~data_zerod.modeh) .* scenario.edge.rho_edge.value; 
scenario.pedestal.q_ped.value 		= data_zerod.modeh .* interp1_itm(profil0d.temps,profil0d.qjli(:,end-1),data_zerod.temps,'pchip','extrap') + ...
					  (~data_zerod.modeh) .* scenario.edge.q_edge.value;  
scenario.pedestal.pressure_ped.value 		= data_zerod.modeh .* data_zerod.pped;
scenario.pedestal.vtor_ped.value 		=  data_zerod.modeh .* interp1_itm(profil0d.temps,profil0d.vtor(:,end-1),data_zerod.temps,'pchip','extrap') + ...
					  (~data_zerod.modeh) .* scenario.edge.vtor_edge.value;
%
if length(data_zerod.temps) == 1
	scenario.references.plh.value 		= z0dstruct.z0dinput.cons.plh(end);
	scenario.references.picrh.value 	= z0dstruct.z0dinput.cons.picrh(end);
	scenario.references.pecrh.value 	= z0dstruct.z0dinput.cons.pecrh(end); 
	scenario.references.pnbi.value 		= z0dstruct.z0dinput.cons.pnbi(end);  
	scenario.references.ip.value		= z0dstruct.z0dinput.cons.ip(end); 
else
	scenario.references.plh.value 		= z0dstruct.z0dinput.cons.plh;
	scenario.references.picrh.value 	= z0dstruct.z0dinput.cons.picrh;
	scenario.references.pecrh.value 	= z0dstruct.z0dinput.cons.pecrh; 
	scenario.references.pnbi.value 		= z0dstruct.z0dinput.cons.pnbi;  
	scenario.references.ip.value		= z0dstruct.z0dinput.cons.ip; 
end
if length(data_zerod.temps) == 1
	scenario.references.bvac_r.value 	= sigma_bvac_r .* z0dstruct.z0dinput.geo.R(end) .* z0dstruct.z0dinput.geo.b0(end); 
else
	scenario.references.bvac_r.value 	= sigma_bvac_r .* z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0; 
end 
if length(data_zerod.temps) == 1
	scenario.references.zeffl.value 	= z0dstruct.z0dinput.cons.zeff(end);
	scenario.references.nbar.value		= real(z0dstruct.z0dinput.cons.nbar(end));
        if imag(z0dstruct.z0dinput.cons.nbar(end))
	    scenario.references.gas_puff.value      = imag(z0dstruct.z0dinput.cons.nbar(end));
	end
	scenario.references.xecrh.value 	= z0dstruct.z0dinput.cons.xece(end);
	scenario.references.pol_flux.value 	= z0dstruct.z0dinput.cons.flux(end) .* 2 .* pi; 
	scenario.references.enhancement.value 	= z0dstruct.z0dinput.cons.hmore(end); 
	scenario.references.isotopic.value 	= real(z0dstruct.z0dinput.cons.iso(end)); 
	scenario.references.nbi_td_ratio.value	= z0dstruct.z0dinput.cons.ftnbi(end);
else
	scenario.references.zeffl.value 	= z0dstruct.z0dinput.cons.zeff;
	scenario.references.nbar.value		= real(z0dstruct.z0dinput.cons.nbar);
        if any(imag(z0dstruct.z0dinput.cons.nbar(end)))
	    scenario.references.gas_puff.value      = imag(z0dstruct.z0dinput.cons.nbar);
	end
	scenario.references.xecrh.value 	= z0dstruct.z0dinput.cons.xece;
	scenario.references.pol_flux.value 	= z0dstruct.z0dinput.cons.flux .* 2 .* pi; 
	scenario.references.enhancement.value 	= z0dstruct.z0dinput.cons.hmore; 
	scenario.references.isotopic.value 	= real(z0dstruct.z0dinput.cons.iso); 
	scenario.references.nbi_td_ratio.value	= z0dstruct.z0dinput.cons.ftnbi;
end
%
if isfield(data_zerod,'dsol')
	scenario.sol.l_ne_sol.value 			= data_zerod.dsol;
elseif length(data_zerod.temps) == 1
	scenario.sol.l_ne_sol.value 			= z0dstruct.z0dinput.geo.a(end) ./ 100;
else
	scenario.sol.l_ne_sol.value 			= z0dstruct.z0dinput.geo.a ./ 100;
end
scenario.sol.l_ni_sol.value 			= scenario.sol.l_ne_sol.value;
scenario.sol.l_qe_sol.value 			= scenario.sol.l_ne_sol.value ./ 0.62;
scenario.sol.l_qi_sol.value 			= scenario.sol.l_qe_sol.value;
scenario.sol.l_te_sol.value 			= 3 .* scenario.sol.l_qe_sol.value;  
scenario.sol.l_ti_sol.value 			= 3 .* scenario.sol.l_qe_sol.value;
scenario.sol.p_rad_sol.value			= data_zerod.pradsol;
scenario.sol.gaz_puff       			= fuelling_data.gas_puff;

%
scenario.vol_ave.te_ave.value 		= data_zerod.tem;
scenario.vol_ave.ti_ave.value 		= data_zerod.tem .* data_zerod.tite;
scenario.vol_ave.ne_ave.value 		= data_zerod.nem;
if length(data_zerod.temps) == 1
	nem   = z0dstruct.zerod.nem;
        t     = z0dstruct.zerod.temps;
	scenario.vol_ave.dne_ave_dt.value 	= (nem(end) - nem(end-1)) ./ (t(end) - t(end -1));
else
	scenario.vol_ave.dne_ave_dt.value 	= zdxdt(data_zerod.nem,data_zerod.temps);
end
scenario.vol_ave.ni_ave.value 		= data_zerod.nim;
scenario.vol_ave.zeff_ave.value 	= data_zerod.zeff; 
scenario.vol_ave.ti_o_te_ave.value	= data_zerod.tite;
scenario.vol_ave.meff_ave.value 	= data_zerod.meff;
scenario.vol_ave.pellet_flux.value 	= data_zerod.frac_pellet .* data_zerod.n0a;
scenario.vol_ave.omega_ave.value 	= data_zerod.wrad;
switch z0dstruct.z0dinput.option.mino
case 'He3'
	switch z0dstruct.z0dinput.option.gaz
	case 4
		nHe3m = z0dstruct.z0dinput.option.cmino .* data_zerod.nhem;
		nHem  = max(0,data_zerod.nhem - nHe3m);
	otherwise
		nHe3m = z0dstruct.z0dinput.option.cmino .* data_zerod.n1m;
		nHem  = max(0,data_zerod.nhem - nHe3m);
	end
otherwise
	nHem  = data_zerod.nhem;
	nHe3m = 0 .* nHem;
end
%scenario.vol_ave.nions_ave.value 		= cat(2,max(0,data_zerod.n1m - data_zerod.nDm - data_zerod.nTm),data_zerod.nDm,data_zerod.nTm,nHe3m,nHem, ...
%						data_zerod.nimpm,z0dstruct.z0dinput.option.rimp .* data_zerod.nimpm);
scenario.vol_ave.nions_ave		= cat(2,max(0,data_zerod.n1m - data_zerod.nDm - data_zerod.nTm),data_zerod.nDm,data_zerod.nTm,nHe3m,nHem, ...
						data_zerod.nimpm,z0dstruct.z0dinput.option.rimp .* data_zerod.nimpm,data_zerod.nwm);

% 
% composition
% 
scenario.composition.amn 		= [1,2,3,3,4,ceil(7/3 .* z0dstruct.z0dinput.option.zimp),ceil(7/3 .* z0dstruct.z0dinput.option.zmax),183.84];
scenario.composition.zn	 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
scenario.composition.zion 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
scenario.composition.imp_flag 		= zeros(size(scenario.composition.zion));
scenario.composition.rot_imp_flag 	= [0,0,0,0,0,1,0,0];
switch z0dstruct.z0dinput.option.gaz
case 1
	scenario.composition.pellet_amn = 1;
	scenario.composition.pellet_zn 	= 1;
	scenario.composition.nbi_amn 	= [1,2];
	scenario.composition.nbi_zn 	= [1,1];
case 2
	scenario.composition.pellet_amn = 2;
	scenario.composition.pellet_zn 	= 1;
	scenario.composition.nbi_amn 	= [1,2];
	scenario.composition.nbi_zn 	= [1,1];
case 3
	scenario.composition.pellet_amn = [2,3];
	scenario.composition.pellet_zn 	= [1,1];
	scenario.composition.nbi_amn 	= [2,3];
	scenario.composition.nbi_zn 	= [1,1];
case 4
	scenario.composition.pellet_amn = [1,2];
	scenario.composition.pellet_zn 	= [1,1];
	scenario.composition.nbi_amn 	= [1,2];
	scenario.composition.nbi_zn 	= [1,1];
otherwise
	error(sprintf('unknown gaz option %d',z0dstruct.z0dinput.option.gaz));
end 
scenario.compositions = copy_composition_to_compositionstype(scenario.composition);
%
% temps
%
scenario.time = data_zerod.temps;

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
%
scenario.codeparam.codename    = 'METIS4ITM';
scenario.codeparam.codeversion = num2str(zinebversion);
scenario.codeparam.parameters  = parameters;
%
data             = [];
data.diboot      = data_zerod.diboot;
data.dw          = data_zerod.dw;
data.dpfus       = data_zerod.dpfus;
data.dini        = data_zerod.dini;
data.difcurconv  = data_zerod.difcurconv;
data.stf         = data_zerod.stf;
data.nb          = data_zerod.nb;
data.texte_diary = num2str(abs(texte_diary));
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
scenario.codeparam.output_diag  = parameters;
%
scenario.codeparam.output_flag = error_flag * ones(size(data_zerod.temps));

%
% end wrapping CPO scenario
% 


% appel du mode test
function post = test_metis_itm


b0    = 5.3;
R     = 6.2;
a     = 2;
K95     = 1.7;
d95     = 0.33;
ip    = 15e6;
nbar  = 0.85 .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
plh   = 0;
picrh = 17e6;
pecrh = 0e6;
pnbi  = 33e6;
zeff  = 1.5;
li    = 0.7;
hmore = 1.0;
iso   = 1;
ftnbi = 0;
ane   = 1.01;
einj  = 1e6;
frad  = 0.6;
rw    = 0.7;
fpped = 1;
fprad  = 1/3;
sepa_option = z0dsepanew2;
sepa_option = sepa_option.valeur;
sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
sepa_option.zxup      = 1.687;    % upper altitude X point (minor radius unit)
sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
sepa_option.za        = 0.65;       % altitude of the magnetic axis (m) [0.9]
sepa_option.a         = 2;         % minor radius (m) [2]
sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
sepa_option.zxdo      = 2.001;       % lower altitude X point (minor radius unit)
sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
sepa_option.b0        = 5.3 ./ (1 - 2 ./ 6.2 - 1 / 6.2);      % magnetic field at R0
sepa_option.delta     = 1;      % magnetic field at R0
sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]


% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K  = 2 .* (K95 - 0.5 .* (1 - 0.95 .^4)) ./ (1 + 0.95 .^4);

d  = d95 ./ (0.95 .^ 2);

temps    = linspace(1,301,301)';

xece = 0.3;

if rand(1) > 0.33
	ga = 2;
elseif rand(1) > 0.5
	ga = 3;
else
	ga = 4;
end
fprintf('testing for gas option %d\n',ga);

z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plh,picrh,pecrh,pnbi,zeff,xece,hmore,iso,ftnbi,0, ...
                          'gaz',ga,'frhe0',0,'tauhemul',5,'ane',4,'vane',ane, ...
			  'scaling',0,'l2hscalin',0,'modeh',1,'l2hmul',0,'fpped',fpped, ...
			  'xiioxie',2,'kishape',3,'qdds',0.95,'kidds',3,'vloop',0,'runaway',0,'modeboot',1,'vref',0, ...
			  'zeff',0,'zmax',18,'zimp',4,'rimp',0.06,'matthews',0, ...
			  'frad',frad,'rw',rw,'angle_ece',90,'synergie',1, ...
			  'sens',1,'angle_nbi',50,'einj',einj,'rtang',R+a/4,'lhmode',3,'etalh',0.8, ...
			  'npar0',2,'freqlh',5,'wlh',a/2,'xlh',0,'dlh',0,'fwcd',0, ...
			  'mino','T','cmin',1,'nphi',25,'freq',72,'sitb',0,'tae',0, ...
			  'smhd',0,'tmhd',inf,'rip',0,'fprad',fprad,'li',1,'configuration',2);

%z0dinput.cons.ip = z0dinput.cons.ip .* min(1, 0.1 + temps ./ 100);
z0dinput.cons.nbar = z0dinput.cons.nbar .* min(1, 0.1 + temps ./ 300);
z0dinput = z0separatrix(z0dinput,sepa_option);

[zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);

info = metis4itm(1);
fv = info.valeur;
noms = fieldnames(fv);
for k = 1:length(noms)
	if ~isfield(z0dinput.option,noms{k})
		z0dinput.option.(noms{k})= fv.(noms{k});	
	end
end


post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;



% appel du mode  load
function data = load_metis_itm(filename)

data = load(filename,'post');
if isempty(data)
	error('This is not a Metis file');

else
	data =data.post;
end
if ~isfield(data,'zerod')
	error('This is not a Metis file');
end
if ~isfield(data,'z0dinput')
	error('This is not a Metis file');
end
         
% compatibilite entre version
if isfield( data.z0dinput,'exp') & ~isfield( data.z0dinput,'exp0d')
	data.z0dinput.exp0d = data.z0dinput.exp;
end   

% mise a niveau des anciennes versions
model = zerod_init(-2,data.z0dinput.shot,data.z0dinput.option.gaz,data.z0dinput.cons.temps);

% option
noms = fieldnames(model.option);
for k = 1:length(noms)
	if ~isfield(data.z0dinput.option,noms{k})
		data.z0dinput.option.(noms{k}) = model.option.(noms{k});
	end
end
% info
noms = fieldnames(model.info);
for k = 1:length(noms)
	if ~isfield(data.z0dinput.info,noms{k})
		data.z0dinput.info.(noms{k}) = model.info.(noms{k});
	end
end
% zsinfo
noms = fieldnames(model.zsinfo);
for k = 1:length(noms)
	if ~isfield(data.z0dinput.zsinfo,noms{k})
		data.z0dinput.zsinfo.(noms{k}) = model.zsinfo.(noms{k});
	end
end
% profinfo
if isfield(data.z0dinput,'profinfo')
	noms = fieldnames(model.profinfo);
	for k = 1:length(noms)
		if ~isfield(data.z0dinput.profinfo,noms{k})
			data.z0dinput.profinfo.(noms{k}) = model.profinfo.(noms{k});
		end
	end
else
	noms = fieldnames(model.profinfo);
	for k = 1:length(noms)
		data.z0dinput.profinfo.(noms{k}) = model.profinfo.(noms{k});
	end        
end
% consignes
noms = fieldnames(model.cons);
for k = 1:length(noms)
	if ~isfield(data.z0dinput.cons,noms{k})
		data.z0dinput.cons.(noms{k}) = model.cons.(noms{k});
	end
end
% geo
noms = fieldnames(model.geo);
for k = 1:length(noms)
	if ~isfield(data.z0dinput.geo,noms{k})
		data.z0dinput.geo.(noms{k}) = model.geo.(noms{k});
	end
end

%zerod data
noms = fieldnames(model.zsinfo);
for k = 1:length(noms)
	if ~isfield(data.zerod,noms{k})
		data.zerod.(noms{k}) = NaN .* data.z0dinput.cons.temps;
	end
end
    
% profil 0d
if isfield(data,'profil0d')   
	if isfield(data.profil0d,'temps')
		prnan = NaN .* ones(length(data.profil0d.temps),21);
	else
		prnan = NaN .* ones(length(data.z0dinput.cons.temps),21);
	end
	if ~isfield(data.profil0d,'xli')
		data.profil0d.xli = linspace(0,1,21);
	end

	noms = fieldnames(model.profinfo);
	for k = 1:length(noms)
		if ~isfield(data.profil0d,noms{k})
			data.profil0d.(noms{k}) = prnan;
		end
	end        
else
	prnan = NaN .* ones(length(data.z0dinput.cons.temps),21);
	noms = fieldnames(model.profinfo);
	for k = 1:length(noms)
		data.profil0d.(noms{k}) = prnan;
	end       
	data.profil0d.xli   = linspace(0,1,21);
	data.profil0d.temps = data.z0dinput.cons.temps;
end
   
% mise a jour de la structure experimentale vide
noms = fieldnames(data.z0dinput.zsinfo);
if ~isfield(data.z0dinput,'exp0d')
	data.z0dinput.exp0d=[];
end
exp0d  = data.z0dinput.exp0d;
if isfield(exp0d,'temps')
	texp = exp0d.temps;
else
	texp = data.z0dinput.cons.temps;
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
data.z0dinput.exp0d = exp0d;

% cas de deux injecteur de neutre
if isfield(data.z0dinput.option,'nb_nbi')
    if data.z0dinput.option.nb_nbi == 2 
	error('METIS4ITM: this version is incompatible with simulation made with 2 neutral beam injectors');
    end
end


% cas du couplage avec simulink
if isfield(data.z0dinput,'system')

    	zassignin('base','post.zerod',data.zerod);
    	zassignin('base','post.z0dinput',data.z0dinput);
    	zassignin('base','post.profil0d',data.profil0d);
    	zassignin('base','z0dinput',data.z0dinput);

    	if isfield(data,'simout')
   	 	zassignin('base','post.simout',data.simout);
   	 	zassignin('base','simout',data.simout);
    	end	

    	vv = ver('simulink');
	if length(vv) > 0
    		% restauration du model simulink
		chem = which(data.z0dinput.system.name);
		if isempty(chem)
			fprintf('Simulink system %s does not exist in Matlab path\n',data.z0dinput.system.name);
			if isdir(fileparts(data.z0dinput.system.fullname))
				fprintf('Simulink system %s will be created in %s \n',data.z0dinput.system.name, ...
				            fileparts(data.z0dinput.system.fullname));
				newname = data.z0dinput.system.fullname;
				try
					addpath(fileparts(data.z0dinput.system.fullname));	
				end		
			else
				fprintf('Simulink system %s will be created in %s \n',data.z0dinput.system.name,pwd);						
				[void,newname,ext] = fileparts(data.z0dinput.system.fullname);
				newname = strcat(newname,ext);
				try
					addpath(pwd)
				end			
			end
			[fid,mess] = fopen(newname,'w');
			if fid >= 3
				fprintf(fid,'%s',data.z0dinput.system.mdl);
				fclose(fid);
				drawnow
				void=which(data.z0dinput.system.name);
				evalin('base',data.z0dinput.system.name);
			else 
				error(mess)
			end
		else
			%[s,mdl_loc]  = unix(sprintf('cat %s',chem));
			% lecture du fichier
			fid = fopen(chem,'r');
			if fid > 0
				mdl_loc = char(fread(fid,Inf,'char')');
				fclose(fid);
			else
				mdl_loc ='';
			end
			if strcmp(mdl_loc,data.z0dinput.system.mdl)
				% model identique ouverture
				evalin('base',data.z0dinput.system.name);
			else
				fprintf('Simulink system %s as some differences with system saved in Metis data\n',data.z0dinput.system.name);
				tage = num2str(fix(datenum(clock)*1e5));
                		reste = 62 - length(tage);
                		reste = min(length(data.z0dinput.system.name),reste);
				newname =strcat(data.z0dinput.system.name(1:reste),tage);
				fprintf('Simulink system %s will be created in same directory\n',newname);
				[fid,mess] = fopen(fullfile(fileparts(chem),strcat(newname,'.mdl')),'w');
				if fid >= 3
					fprintf(fid,'%s',data.z0dinput.system.mdl);
					fclose(fid);
					drawnow
					void=which(newname);
					evalin('base',newname);
				else 
					error(mess)
				end
			end				
		end
	       zassignin('base','z0dinput',data.z0dinput.system.z0dinput);

	end
end
    

% cette fonction extrait les options de metis du cpo scenario
function options = scenario2option(scenario)

% add codeparam to options
info = metis4itm(1);
options = info.valeur;
codeparam = scenario.codeparam.parameters;
if ~isempty(codeparam)
	tpn = tempname;
	fid = fopen(tpn,'w');
	if iscell(codeparam)
		fprintf(fid,'%s\n',codeparam{end});
	else
		fprintf(fid,'%s\n',codeparam);
	end
	fclose(fid);
	info  = xml_read(tpn);
	delete(tpn);
	noms = fieldnames(info);
	for k=1:length(noms)
		if isfield(options,noms{k})
			options.(noms{k}) = info.(noms{k});
		end
	end
end 

% cette fonction extrait les informations du CPO scenario pour en faire une structure z0dinput de metis
function z0dinput = cpo2metis_input(scenario)

% parametre et info des parametres 
info = metis4itm(1);
z0dinput.info          = info.info;
z0dinput.option        = scenario2option(scenario);


% langue
z0dinput.langue        =  'anglais';
% variable de sorties
z0dinput.zsinfo        = zero1t;
z0dinput.profinfo      = z0dprofinfo;
z0dinput.mode_exp      = -sum(abs('ITM'));

% recherche des infos
z0dinput.cons.temps    = scenario.time;
% choix de la composition
if ~isempty(scenario.references.isotopic.value)
	z0dinput.cons.iso      = real(scenario.references.isotopic.value);		
else
	z0dinput.cons.iso      = zeros(size(z0dinput.cons.temps));
end 
if ~isempty(scenario.references.nbi_td_ratio.value)
	z0dinput.cons.ftnbi    = scenario.references.nbi_td_ratio.value;
else
	z0dinput.cons.ftnbi    = zeros(size(z0dinput.cons.temps));
end
% cette donnees doit etre non vide
z0dinput.cons.ip       = scenario.references.ip.value;
if ~isempty(scenario.references.pol_flux.value)
	z0dinput.cons.flux     = scenario.references.pol_flux.value ./ 2 ./ pi;
else
	z0dinput.cons.flux     = zeros(size(z0dinput.cons.temps));
end
% cette donnees doit etre non vide
z0dinput.cons.nbar     = scenario.references.nbar.value;
% gaspuff
if isfield(scenario.references,'gaspuff')
    if ~isempty(scenario.references.gaspuff.value)
	  z0dinput.cons.nbar     = z0dinput.cons.nbar + sqrt(-1) .* scenario.references.gaspuff.value;
    end
end

if ~isempty(scenario.references.picrh.value)
	z0dinput.cons.picrh    = scenario.references.picrh.value;
else
	z0dinput.cons.picrh    = zeros(size(z0dinput.cons.temps));
end
if ~isempty(scenario.references.plh.value)
	z0dinput.cons.plh      = scenario.references.plh.value; 
else
	z0dinput.cons.plh      = zeros(size(z0dinput.cons.temps));
end
if ~isempty(scenario.references.pnbi.value)
	z0dinput.cons.pnbi     = scenario.references.pnbi.value;
else
	z0dinput.cons.pnbi     = zeros(size(z0dinput.cons.temps));
end
if ~isempty(scenario.references.pecrh.value)
	z0dinput.cons.pecrh    = scenario.references.pecrh.value;
else
	z0dinput.cons.pecrh    = zeros(size(z0dinput.cons.temps));
end
if ~isempty(scenario.references.pecrh.value)
	z0dinput.cons.hmore    = scenario.references.enhancement.value;
else
	z0dinput.cons.hmore    = ones(size(z0dinput.cons.temps));
end
% cette donnees doit etre non vide
z0dinput.cons.zeff     = scenario.references.zeffl.value;
z0dinput.cons.zeff(~isfinite(z0dinput.cons.zeff)) = 3;
if ~isempty(scenario.references.xecrh.value)
	z0dinput.cons.xece     = scenario.references.xecrh.value; 
else
	z0dinput.cons.xece     = zeros(size(z0dinput.cons.temps));
end

% selon les donnees disponibles
if ~isempty(scenario.eqgeometry.boundary.r) & ~isempty(scenario.eqgeometry.boundary.z)
	% cas separatrice donnees par points
	sepa.R = squeeze(double(scenario.eqgeometry.boundary.r));
	sepa.Z = squeeze(double(scenario.eqgeometry.boundary.z));
	% calcul des moments
	% la courbe doit etre fermee
	if (sepa.R(1,1) ~= sepa.R(1,end)) | (sepa.Z(1,1) ~= sepa.Z(1,end))
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
	dZdx  = pdederive(xu,sepa.Z,2,2,2,1);
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


	z0dinput.geo.R       = geo.R;      % grand rayon du plasma (m)
	z0dinput.geo.z0      = geo.z0;     % centre geometricque du plasma en Z (m)
	z0dinput.geo.a       = geo.a;      % petit rayon du plasma (m)
	z0dinput.geo.K       = geo.K;     % elongation (b/a)
	z0dinput.geo.d       = geo.d;    % triangularite haute (definition entree de helena)
	z0dinput.geo.b0      = abs(scenario.references.bvac_r.value ./ z0dinput.geo.R); % champ toroidal a R = 6.2 m (T)


	z0dinput.exp0d.Rsepa = sepa.R;       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa = sepa.Z - geo.z0 * ones(1,size(sepa.Z,2));       % vecteur Z des points de la separtrice (m)
else
	% cas separatrice donnees par moments	
	z0dinput.exp0d.Rsepa 	= [];       % vecteur R des points de la separatrice (m)
	z0dinput.exp0d.Zsepa 	= [];       % vecteur Z des points de la separtrice (m)
   	z0dinput.geo.a     	= scenario.eqgeometry.a_minor;
   	z0dinput.geo.R     	= scenario.eqgeometry.geom_axis.r ;
   	z0dinput.geo.K     	= scenario.eqgeometry.elongation;
   	z0dinput.geo.d     	= (scenario.eqgeometry.tria_upper + scenario.eqgeometry.tria_lower) ./ 2;
   	z0dinput.geo.b0    	= abs(scenario.references.bvac_r ./ z0dinput.geo.R);
  	z0dinput.geo.z0    	= scenario.eqgeometry.geom_axis.z;
	if ~isempty(scenario.global.volume)
   		z0dinput.geo.vp    	= scenario.global.volume;
	else
   		z0dinput.geo.vp    	= NaN .* ones(size(z0dinput.cons.temps));
	end
 	if ~isempty(scenario.global.area_pol)
  		z0dinput.geo.sp    	= scenario.global.area_pol;
	else
   		z0dinput.geo.sp    	= NaN .* ones(size(z0dinput.cons.temps));
	end
  	if ~isempty(scenario.global.area_ext)
  		z0dinput.geo.sext  	= scenario.global.area_ext;
	else
   		z0dinput.geo.sext   	= NaN .* ones(size(z0dinput.cons.temps));
	end
end

% informations generales   
if iscell(scenario.datainfo.whatref.machine)
	z0dinput.machine   = scenario.datainfo.whatref.machine{1};
else
	z0dinput.machine   = scenario.datainfo.whatref.machine;
end
z0dinput.shot      = real(scenario.datainfo.whatref.shot);

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
   z0dinput.cons.iso = real(z0dinput.cons.iso) .* ones(size(z0dinput.cons.temps)); 
end
% consigne d'injection de tritium par nbi (fraction de la puissance)
if ~isfield(z0dinput.cons,'ftnbi')
   z0dinput.cons.ftnbi = min(1,real(z0dinput.cons.iso) .* 0.5);
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


%
% lecture des consignes pour la version evolution
%
function [option,time,cons1t,geo1t,sepa1t] = cpo2metis1t(scenario)

% extraction des options METIS de codeparam.
option = scenario2option(scenario);


% valeur du temps
time = scenario.time(end);
cons1t.temps = time;
 
% les consignes a 1 temps
% recherche des infos
% choix de la composition
if ~isempty(scenario.references.isotopic.value)
	cons1t.iso      = real(scenario.references.isotopic.value(end));		
else
	cons1t.iso      = 0;
end 
if ~isempty(scenario.references.nbi_td_ratio.value)
	cons1t.ftnbi    = scenario.references.nbi_td_ratio.value(end);
else
	cons1t.ftnbi    = 0;
end
% cette donnees doit etre non vide
cons1t.ip       = scenario.references.ip.value(end);
if ~isempty(scenario.references.pol_flux.value)
	cons1t.flux     = scenario.references.pol_flux.value(end);
else
	cons1t.flux     = 0;
end
% cette donnees doit etre non vide
cons1t.nbar     = scenario.references.nbar.value(end);
% gaspuff
if isfield(scenario.references,'gaspuff')
    if ~isempty(scenario.references.gaspuff.value)
	  cons1t.nbar     = cons1t.nbar + sqrt(-1) .* scenario.references.gaspuff.value(end);
    end
end

if ~isempty(scenario.references.picrh.value)
	cons1t.picrh    = scenario.references.picrh.value(end);
else
	cons1t.picrh    = 0;
end
if ~isempty(scenario.references.plh.value)
	cons1t.plh      = scenario.references.plh.value(end); 
else
	cons1t.plh      = 0;
end
if ~isempty(scenario.references.pnbi.value)
	cons1t.pnbi     = scenario.references.pnbi.value(end);
else
	cons1t.pnbi     = 0;
end
if ~isempty(scenario.references.pecrh.value)
	cons1t.pecrh    = scenario.references.pecrh.value(end);
else
	cons1t.pecrh    = 0;
end
if ~isempty(scenario.references.pecrh.value)
	cons1t.hmore    = scenario.references.enhancement.value(end);
else
	cons1t.hmore    = 1;
end
% cette donnees doit etre non vide
cons1t.zeff     = scenario.references.zeffl.value(end);
if ~isempty(scenario.references.xecrh.value)
	cons1t.xece     = scenario.references.xecrh.value(end); 
else
	cons1t.xece     = 0;
end


% selon les donnees disponibles
if ~isempty(scenario.eqgeometry.boundary.r) & ~isempty(scenario.eqgeometry.boundary.z)
	% cas separatrice donnees par points
	sepa.Rsepa = squeeze(double(scenario.eqgeometry.boundary.r(end,:)));
        sepa.Rsepa = sepa.Rsepa(:)'; 
	sepa.Zsepa = squeeze(double(scenario.eqgeometry.boundary.z(end,:)));
        sepa.Zsepa = sepa.Zsepa(:)'; 
	% calcul des moments
	% la courbe doit etre fermee
	if (sepa.Rsepa(1,1) ~= sepa.Rsepa(1,end)) | (sepa.Zsepa(1,1) ~= sepa.Zsepa(1,end))
		sepa.Rsepa(:,end+1) = sepa.Rsepa(:,1);
		sepa.Zsepa(:,end+1) = sepa.Zsepa(:,1);
	end

	% calcul des moments
	% centre pour angle d'integration
	rc = mean(sepa.Rsepa,2);
	zc = mean(sepa.Zsepa,2);
	vc = ones(1,size(sepa.Rsepa,2));
	uc = unwrap(angle((sepa.Rsepa-rc*vc) + sqrt(-1) .* (sepa.Zsepa  -zc*vc)));
	uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
	uc(:,1)   = uc(:,end) + 2 .* pi;
	xu    = linspace(0,1,length(vc));
	%dudx  = pdederive(xu,uc,2,2,2,1);
	%dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
	%dudx(:,end) = dudx(:,1);
	dRdx  = pdederive(xu,sepa.Rsepa,2,2,2,1);
	%dZdx  = pdederive(xu,sepa.Zsepa,2,2,2,1);
	% calcul de R0 et Z0
	maskrmax  = (sepa.Rsepa == (max(sepa.Rsepa,[],2) * vc));
	geo.z0        = sum(sepa.Zsepa .* maskrmax,2) ./ sum(maskrmax,2);
	% recalcul des parametres sur le vecteur final
	rmin  = min(sepa.Rsepa,[],2);
	rmax  = max(sepa.Rsepa,[],2);
	geo.a = 0.5 .* (rmax - rmin);
	geo.R = 0.5 .* (rmax + rmin);
	zmin  = min(sepa.Zsepa,[],2);
	zmax  = max(sepa.Zsepa,[],2);
	%geo.z0      =(zmax + zmin) ./ 2;
	geo.K    = abs(trapz(xu,sepa.Zsepa .*  dRdx,2) ./ pi ./ geo.a .^ 2);
	%geo.K  = (zmax -zmin) ./ 2 ./ geo.a;
	rzmax = geo.R;
	rzmin = geo.R;
	for k = 1:size(sepa.Zsepa,1)
		rzmax(k) = sepa.Rsepa(k,min(find(sepa.Zsepa(k,:) == zmax(k))));
		rzmin(k) = sepa.Rsepa(k,min(find(sepa.Zsepa(k,:) == zmin(k))));
	end
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
	geo1t.b0      = abs(scenario.references.bvac_r.value(end) ./ geo.R(end))  ; % champ toroidal a R = 6.2 m (T)
	%
	sepa1t.Rsepa = sepa.Rsepa(end,:);       % vecteur R des points de la separatrice (m)
	sepa1t.Zsepa = sepa.Zsepa(end,:)   - geo.z0(end)   * ones(1,size(sepa.Zsepa,2));       % vecteur Z des points de la separtrice (m)
else
	% cas separatrice donnees par moments	
	sepa1t = [];       % pas de separatrice
   	geo1t.a     	= scenario.eqgeometry.a_minor(end);
   	geo1t.R     	= scenario.eqgeometry.geom_axis.r(end) ;
   	geo1t.K     	= scenario.eqgeometry.elongation(end);
   	geo1t.d     	= (scenario.eqgeometry.tria_upper(end) + scenario.eqgeometry.tria_lower(end)) ./ 2;
   	geo1t.b0    	= abs(scenario.references.bvac_r(end) ./ z0dinput.geo.R);
  	geo1t.z0    	= scenario.eqgeometry.geom_axis.z(end);
end


% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function [coreprof,vtor,vpol] = mapcoreprof(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreprof,sigma_B0_eff)

if z0dstruct.z0dinput.option.Sn_fraction > 0
    error('METIS ITM CPO interface is not compatible whith tin (Sn) in plasma composition (option.Sn_fraction should be 0)');
end

% coreprof.datainfo est le meme que celui de scenario
if ~isfield(coreprof,'datainfo')
	coreprof.datainfo = datainfo_empty;
end

coreprof.time			= profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coreprof.rho_tor_norm 		= xli;
coreprof.rho_tor       		= profil0d.rmx;
if length(profil0d.temps) == 1
	rb0 =   z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0;
	coreprof.toroid_field.time        = profil0d.temps;
	coreprof.toroid_field.r0          = mean(z0dstruct.z0dinput.geo.R);
	coreprof.toroid_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0(end) ./ coreprof.toroid_field.r0;
	coreprof.toroid_field.b0prime     = (rb0(end) - rb0(1)) ./ (z0dstruct.z0dinput.cons.temps(end) -z0dstruct.z0dinput.cons.temps(1));
	coreprof.drho_dt       		  = data_zerod.drmdt ./ profil0d.rmx(end) .*profil0d.rmx;
else 
	rb0 =   interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
	coreprof.toroid_field.time        = profil0d.temps;
	coreprof.toroid_field.r0          = mean(z0dstruct.z0dinput.geo.R);
	coreprof.toroid_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ coreprof.toroid_field.r0;
	coreprof.toroid_field.b0prime     = z0dxdt(coreprof.toroid_field.b0,coreprof.toroid_field.time);
	coreprof.drho_dt       		  = zdxdt(profil0d.rmx,profil0d.temps);
end
%
coreprof.composition.amn 		= [1,2,3,3,4,ceil(7/3 .* z0dstruct.z0dinput.option.zimp),ceil(7/3 .* z0dstruct.z0dinput.option.zmax),183.84];
coreprof.composition.zn	 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
coreprof.composition.zion 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax,74];
coreprof.composition.imp_flag 		= zeros(size(coreprof.composition.zion));
coreprof.compositions = copy_composition_to_compositionstype(coreprof.composition);

%
coreprof.psi.value 		= profil0d.psi;
coreprof.psi.derivative =  pdederive(profil0d.rmx,profil0d.psi,0,2,2,1);
coreprof.psi.source             = 'METIS current diffusion';
coreprof.psi.flag               = 2 .* ones(size(coreprof.time));
asser = fix(interp1_itm(data_zerod.temps,data_zerod.asser,profil0d.temps,'pchip','extrap'));
ip    = interp1_itm(data_zerod.temps,data_zerod.asser,profil0d.temps,'pchip','extrap');
value     =  ip .* (asser == 0) + profil0d.psi(:,end) .* asser;
coreprof.psi.boundary.value     =  value * cat(2,1,0,0);
coreprof.psi.boundary.source    =  'METIS';
coreprof.psi.boundary.type      =  1 + (asser == 0);
% en attendant la correction
%coreprof.psi.boundary.rho_tor   = profil0d.rmx(:,end);
coreprof.psi.boundary.rho  = profil0d.rmx(:,end);
coreprof.psi.jni.value 		= profil0d.jni;
coreprof.psi.jni.integral 	= cumtrapz(profil0d.xli,profil0d.jni .* profil0d.spr,2);
coreprof.psi.jni.source   	=  'METIS';
coreprof.psi.sigma_par.value 	= 1 ./ max(eps,profil0d.eta);
coreprof.psi.sigma_par.source   =  'Sauter law';
%
coreprof.te.value 		= profil0d.tep;
coreprof.te.derivative =  pdederive(profil0d.rmx,profil0d.tep,0,2,2,1);
coreprof.te.source 		= 'METIS';
coreprof.te.flag 		= 2 .* ones(size(profil0d.temps));
coreprof.te.boundary.value     	= profil0d.tep(:,end) * cat(2,1,0,0);
coreprof.te.boundary.source    	= 'METIS';
coreprof.te.boundary.type      	= ones(size(profil0d.temps));
coreprof.te.boundary.rho_tor   	= profil0d.rmx(:,end);
coreprof.te.source_term.value 	= profil0d.source_el;
coreprof.te.source_term.integral = cumtrapz(profil0d.xli,profil0d.source_el .* profil0d.vpr,2);
coreprof.te.source_term.source 	= 'METIS';
coreprof.te.transp_coef.diff 	= profil0d.xie;
coreprof.te.transp_coef.vconv 	= zeros(length(profil0d.temps),length(profil0d.xli));
coreprof.te.transp_coef.source 	= 'METIS';
%
vmati    = ones(1,length(coreprof.composition.zn));
smati    = cat(2,size(profil0d.tip,1),size(profil0d.tip,2),length(vmati));
coreprof.ti.value 		= reshape(profil0d.tip(:) * vmati,smati);
dtidrho                         = pdederive(profil0d.rmx,profil0d.tip,0,2,2,1);
coreprof.ti.derivative          = reshape(dtidrho(:) * vmati,smati);
coreprof.ti.source 		= 'METIS';
coreprof.ti.flag 		= 2 .* ones(size(profil0d.temps));     %bug dans l'UAL
coreprof.ti.boundary.value     	= profil0d.tip(:,end)  * cat(2,1,0,0);
coreprof.ti.boundary.source    	= 'METIS';
coreprof.ti.boundary.type      	= ones(size(profil0d.temps));
coreprof.ti.boundary.rho_tor   	= profil0d.rmx(:,end);
coreprof.ti.source_term.value 	= reshape(profil0d.source_ion(:) * vmati,smati);
intsti                          = cumtrapz(profil0d.xli,profil0d.source_ion .* profil0d.vpr,2);
coreprof.ti.source_term.integral = reshape(intsti(:) * vmati,smati);
coreprof.ti.source_term.source 	= 'METIS';
coreprof.ti.transp_coef.diff 	= profil0d.xii;
coreprof.ti.transp_coef.vconv 	= zeros(length(profil0d.temps),length(profil0d.xli),length(vmati));
coreprof.ti.transp_coef.source 	= 'METIS';
%
coreprof.ne.value 		= profil0d.nep;
coreprof.ne.derivative =  pdederive(profil0d.rmx,profil0d.nep,0,2,2,1);
coreprof.ne.source 		= 'METIS';
coreprof.ne.flag 		= 2 .* ones(size(profil0d.temps));
coreprof.ne.boundary.value     	= profil0d.nep(:,end) * cat(2,1,0,0);
coreprof.ne.boundary.source    	= 'METIS';
coreprof.ne.boundary.type      	= ones(size(profil0d.temps));
coreprof.ne.boundary.rho_tor    = profil0d.rmx(:,end);
ee   = 1.602176462e-19;
snbi = profil0d.pnbi ./ (z0dstruct.z0dinput.option.einj .* ee);
stot = profil0d.s0m + profil0d.s0 + snbi + profil0d.spellet;
coreprof.ne.source_term.value 	= stot;
coreprof.ne.source_term.integral = cumtrapz(profil0d.xli,stot .* profil0d.vpr,2);
coreprof.ne.source_term.source 	= 'METIS';
profil0d.dn(~isfinite(profil0d.dn)) = 0;
profil0d.vn(~isfinite(profil0d.vn)) = 0;
coreprof.ne.transp_coef.diff 	= profil0d.dn;
coreprof.ne.transp_coef.vconv 	= -profil0d.vn;
coreprof.ne.transp_coef.source 	= 'METIS';

%
% densite ioniques
%
ve    = ones(size(profil0d.xli));
nDm   = interp1_itm(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap');
nTm   = interp1_itm(data_zerod.temps,data_zerod.nTm,profil0d.temps,'pchip','extrap');
nDp   = max(1,profil0d.n1p .* ((nDm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1,profil0d.n1p .* ((nTm./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nHp   = max(1,profil0d.n1p - nTp - nDp);
nhep  = max(1,profil0d.nhep);
nz1p  = max(1,profil0d.nzp);
nz2p  = max(1,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);   
nwp   = max(1,profil0d.nwp);
%
switch z0dstruct.z0dinput.option.mino
case 'He3'
	switch z0dstruct.z0dinput.option.gaz
	case 4
		nHe3m = z0dstruct.z0dinput.option.cmino .* data_zerod.nhem;
		nHem  = max(0,data_zerod.nhem - nHe3m);
	otherwise
		nHe3m = z0dstruct.z0dinput.option.cmino .* data_zerod.n1m;
		nHem  = max(0,data_zerod.nhem - nHe3m);
	end
otherwise
	nHem  = data_zerod.nhem;
	nHe3m = 0 .* nHem;
end
frhe3  = nHe3m ./ max(1e11,nHe3m + nHem);
frhe3  = interp1_itm(data_zerod.temps,frhe3,profil0d.temps,'pchip','extrap') * ve;
nions  = NaN * ones(length(profil0d.temps),length(profil0d.xli),8);
nions(:,:,1) = nHp;
nions(:,:,2) = nDp;
nions(:,:,3) = nTp;
nions(:,:,4) = nhep .* frhe3;
nions(:,:,5) = nhep .* (1 - frhe3);
nions(:,:,6) = nz1p;
nions(:,:,7) = nz2p;
nions(:,:,8) = nwp;
% retabli l'electroneutralite (probleme arrondis numeriques)
% enforced Zeff et ne -> improved precision 
% due to resampling, small error on electroneutrality is introduce and tiny difference in Zeff
nizi2 = zeros(size(profil0d.nep));
nizi  = zeros(size(profil0d.nep));
for k= 1:length(coreprof.composition.zn)
  if k < length(coreprof.composition.zn)
        nizi	= nizi  + coreprof.composition.zn(k)     .* squeeze(nions(:,:,k));
        nizi2	= nizi2 + coreprof.composition.zn(k) .^2 .* squeeze(nions(:,:,k));
  else
        nizi	= nizi  + z0wavez(profil0d.tep)      .* squeeze(nions(:,:,k));
        nizi2	= nizi2 + z0wavez(profil0d.tep)  .^2 .* squeeze(nions(:,:,k));      
  end  
end
fprintf('Relative difference on electoneutrality = %g\n',sqrt(mean((profil0d.nep(:) - nizi(:)) .^ 2 ./ profil0d.nep(:) .^ 2 ./ size(profil0d.nep,1))));
fprintf('Relative difference on Zeff = %g\n',sqrt(mean((profil0d.zeff(:) - nizi2(:) ./ nizi(:)) .^ 2 ./ size(profil0d.nep,1))));
profil0d.nep   = nizi;
profil0d.zeff  = nizi2 ./ nizi;
%  figure;
%  subplot(2,1,1)
%  zplotprof(gca,profil0d.temps,profil0d.xli,profil0d.nep,'color','r');
%  zplotprof(gca,profil0d.temps,profil0d.xli,nizi,'color','b','marker','o','linestyle',':');
%  subplot(2,1,2)
%  zplotprof(gca,profil0d.temps,profil0d.xli,profil0d.zeff,'color','r');
%  zplotprof(gca,profil0d.temps,profil0d.xli,nizi2 ./ nizi,'color','b','marker','o','linestyle',':');
%  keyboard

% nefromnions = zeros(size(profil0d.nep));
% for klm = 1:size(nions,3)
%     nefromnions  =  nefromnions + squeeze(nions(:,:,klm)) .* coreprof.composition.zn(klm);
% end
% rapnenions      = profil0d.nep ./ max(1,nefromnions);
% rapnenions(nefromnions == 0)  = 1;
% rapnenions = reshape(rapnenions(:) * ones(1,length(coreprof.composition.zn(:))), ...
%              size(rapnenions,1),size(rapnenions,2),length(coreprof.composition.zn));
% nions      = nions .* rapnenions;
% keyboard

coreprof.ni.value = nions;
for klm = 1:size(nions,3)
    coreprof.ni.derivative(:,:,klm) =  pdederive(profil0d.rmx,squeeze(nions(:,:,klm)),0,2,2,1);
end
coreprof.ni.source 		= 'METIS';
coreprof.ni.flag 		= 2 .* ones(size(profil0d.temps,1),size(nions,3)); % bug dan l'UAL
%
[vtor,vpol,omega,mtor] = z0rot_itm(data_zerod,profil0d,z0dstruct.z0dinput.option,frhe3,z0dstruct.z0dinput.geo,z0dstruct.z0dinput.cons);
coreprof.vtor.value 		= vtor;
for klm = 1:size(vtor,3)
    coreprof.vtor.derivative(:,:,klm) =  pdederive(profil0d.rmx,squeeze(vtor(:,:,klm)),0,2,2,1);
end
coreprof.vtor.source 		= 'METIS';
coreprof.vtor.flag 		= 2 .* ones(size(profil0d.temps,1),size(vtor,3)); % bug dan l'UAL

coreprof.profiles1d.vpol.value 		= vpol;
coreprof.profiles1d.vpol.source 	= 'METIS';

%
% profiles1d
%
coreprof.profiles1d.pe.value     	=  ee .* profil0d.tep .* profil0d.nep;
coreprof.profiles1d.pe.source 	 	= 'METIS'; 
if length(profil0d.temps) > 1
    coreprof.profiles1d.dpedt.value     =  zdxdt(coreprof.profiles1d.pe.value,profil0d.temps);
    coreprof.profiles1d.dpedt.source 	= 'METIS'; 
end
coreprof.profiles1d.pi.value  	 	=  ee .* profil0d.tip .* profil0d.nip;
coreprof.profiles1d.pi.source 	 	= 'METIS'; 
coreprof.profiles1d.pi_tot.value  	=  ee .* profil0d.tip .* profil0d.nip;
coreprof.profiles1d.pi_tot.source 	 = 'METIS'; 
if length(profil0d.temps) > 1
    coreprof.profiles1d.dpi_totdt.value     =  zdxdt(coreprof.profiles1d.pi_tot.value,profil0d.temps);
    coreprof.profiles1d.dpi_totdt.source    = 'METIS'; 
end

coreprof.profiles1d.pr_th.value  	=  coreprof.profiles1d.pe.value + coreprof.profiles1d.pi.value;
coreprof.profiles1d.pr_th.source 	= 'METIS'; 
coreprof.profiles1d.pr_perp.value  	= profil0d.ptot;
coreprof.profiles1d.pr_perp.source 	= 'METIS'; 
coreprof.profiles1d.jtot.value     	= profil0d.jli;
coreprof.profiles1d.jtot.source    	= 'METIS'; 
coreprof.profiles1d.jni.value      	= profil0d.jni;
coreprof.profiles1d.jni.source     	= 'METIS'; 
coreprof.profiles1d.jphi.value      	= profil0d.jeff;
coreprof.profiles1d.jphi.source     	= 'METIS'; 
coreprof.profiles1d.joh.value      	= profil0d.jeff - profil0d.jni;
coreprof.profiles1d.joh.source     	= 'METIS'; 
coreprof.profiles1d.vloop.value      	= - 2 .* pi .* profil0d.dpsidt;
coreprof.profiles1d.vloop.source     	= 'METIS'; 
coreprof.profiles1d.eparallel.value     = profil0d.epar;
coreprof.profiles1d.eparallel.source    = 'METIS'; 
coreprof.profiles1d.e_b.value     = profil0d.epar .* (coreprof.toroid_field.b0 * ones(1,size(profil0d.epar,2)));
coreprof.profiles1d.e_b.source    = 'METIS'; 
coreprof.profiles1d.q.value     	= profil0d.qjli;
coreprof.profiles1d.q.source    	= 'METIS'; 
coreprof.profiles1d.shear.value      	= (ones(size(profil0d.temps)) * profil0d.xli) ./ profil0d.qjli .* pdederive(profil0d.xli,profil0d.qjli,2,2,2,1);
coreprof.profiles1d.shear.source     	= 'METIS'; 
coreprof.profiles1d.mtor.value     	= mtor;
coreprof.profiles1d.mtor.source    	= 'METIS'; 

coreprof.profiles1d.wtor.value     	= omega;
coreprof.profiles1d.wtor.source    	= 'METIS'; 
coreprof.profiles1d.zeff.value      	= profil0d.zeff;
coreprof.profiles1d.zeff.source     	= 'METIS'; 
coreprof.profiles1d.bpol.value      	= profil0d.bpol;
coreprof.profiles1d.bpol.source     	= 'METIS'; 
coreprof.profiles1d.qoh.value      	= profil0d.pohm;
coreprof.profiles1d.qoh.source     	= 'METIS'; 
warning off
coreprof.profiles1d.qei.value      	= pdederive(profil0d.xli,profil0d.qei,0,2,2,1) ./ pdederive(profil0d.xli,profil0d.vpr_tor,0,2,2,1);
warning on
coreprof.profiles1d.qei.value(:,1)      = 2 .* coreprof.profiles1d.qei.value(:,2) - coreprof.profiles1d.qei.value(:,3);
coreprof.profiles1d.qei.source     	= 'METIS'; 
coreprof.profiles1d.sigmapar.value     	= 1./max(1e-308,profil0d.eta);
coreprof.profiles1d.sigmapar.source    	= 'METIS'; 

psid1    = pdederive(profil0d.xli,profil0d.psi,0,2,2,1);
coreprof.profiles1d.dpsidt.value        = profil0d.dpsidt;
coreprof.profiles1d.dpsidt.source       = 'METIS current diffusion';
if length(profil0d.temps) == 1
    coreprof.profiles1d.dvprimedt.value 	    =  profil0d.vpr_tor .* data_zerod.drmdt ./ profil0d.rmx(end);
    coreprof.profiles1d.dvprimedt.source       = 'METIS';
    coreprof.profiles1d.dpsidt_phi.value    = coreprof.profiles1d.dpsidt.value -  (profil0d.xli .* data_zerod.drmdt(end) ./ profil0d.rmx(end,end))  .* psid1;
    coreprof.profiles1d.dpsidt_phi.source   = 'METIS current diffusion';
else
    coreprof.profiles1d.dvprimedt.value 	    = zdxdt(profil0d.vpr_tor,profil0d.temps);
    coreprof.profiles1d.dvprimedt.source       = 'METIS';
    drmdt = interp1_itm(data_zerod.temps,data_zerod.drmdt,profil0d.temps,'pchip','extrap');
    coreprof.profiles1d.dpsidt_phi.value    = coreprof.profiles1d.dpsidt.value - ((drmdt ./ profil0d.rmx(:,end)) * profil0d.xli)  .* psid1;
    coreprof.profiles1d.dpsidt_phi.source   = 'METIS current diffusion';
end

% fast particle density
nfast  = compute_nfast(data_zerod,profil0d,z0dstruct.z0dinput.option,z0dstruct.z0dinput.cons,nions);
coreprof.profiles1d.ns.value = nfast;
coreprof.profiles1d.ns.source 	= 'METIS';

% pression //
rap_p = (max(0,data_zerod.w - data_zerod.wdia) + data_zerod.wrot) ./ max(eps,data_zerod.w);
if length(data_zerod.temps) > 1
    rap_p = interp1(data_zerod.temps,rap_p,profil0d.temps,'linear','extrap');
end
coreprof.profiles1d.pr_parallel.value  	= profil0d.ptot .* (rap_p * ones(size(profil0d.xli)));
coreprof.profiles1d.pr_parallel.source 	= 'METIS'; 




% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_itm(t,y,tt,varargin)


if length(t) <= 1
	rep = y;
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

rep(~isfinite(rep)) = -9.99999e99;


% dans metis la direction toroidal est dans le sens du courant plasma, 
% de telle sorte que Btheta est toujours positif
% le moment injecte par l'IDN est posistif si l'injection est co courant
% le champs toroidal est positif s'il est dans le sens du courant
% calcul de la vitesse de rotation moyenne
% cette fonction calcul la rotation toroidal pour chaque especes 
% ainsi que le rotation poloidal
function [vtor,vpol,omega,mtor] = z0rot_itm(zs,profli,option,frhe3,geo,cons)

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

% calcul du profil de rotation toroidal
% palsma de fond
% switch option.gaz
% case 1
%    zj = 1;
%    aj = 1;
% case 2
%    zj = 1;
%    aj = 2;
% case 3
%    zj = 1;
%    aj = mean(2 .* (1 - cons.iso) +  3 .* cons.iso);
% case 4
%    zj = 2;
%    aj = 4;
% end

% impurete principale
zimp = option.zimp;
aimp = ceil(zimp .* (7/3));

% 2ieme impurete
zmax = option.zmax;
amax = ceil(zmax .* (7/3));


% pour chaque espece d'ions
x     = profli.xli;
ve    = ones(size(x));
vt    = ones(size(profli.n1p,1),1);
nDm   = interp1_itm(zs.temps,zs.nDm,profli.temps,'pchip','extrap');
nTm   = interp1_itm(zs.temps,zs.nTm,profli.temps,'pchip','extrap');
nDp   = max(1e13,profli.n1p .* ((nDm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nTp   = max(1e13,profli.n1p .* ((nTm./ max(1,trapz(x,profli.vpr .* abs(profli.n1p),2)) .* trapz(x,profli.vpr,2)) * ve));
nHp   = max(1e13,profli.n1p - nTp - nDp);
nhep  = max(1e13,profli.nhep .* (1 - frhe3));
nhep3 = max(1e13,profli.nhep .* frhe3);
nz1p  = max(1e13,profli.nzp);
nz2p  = max(1e11,profli.nzp .* option.rimp);  
nwp   = max(1, profli.nwp);
% masse
Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + aimp .* nz1p + amax .* nz2p + 183.84 .* nwp);
% omega homothetic a Ti (cas Jet avec NBI)
% calcul de la rottaion poloidal 
warning off


% formulaire ORNL
lnii    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(1); 
lnhe    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2); 
lnhe3   = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(2); 
lnz1    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zimp); 
lnz2    = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(zmax); 
lnw     = 23 - 0.5 .* log(profli.nep ./ 1e20) + 3/2 .* log(profli.tip ./ 1e3) - log(z0wavez(profli.tep)); 

% pour l'espece principale
% Tokamaks, Wesson p 663
taui_h     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 1)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nHp ./ lnii ./ 1 .^ 4 ;
taui_d     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 2)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nDp ./ lnii ./ 1 .^ 4 ;
taui_t     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nTp ./ lnii ./ 1 .^ 4 ;
taui_he     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 4)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe ./ 2 .^ 4 ;
taui_he3     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 3)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nhep ./ lnhe ./ 2 .^ 4 ;
taui_z1     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* aimp)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz1p ./ lnz1 ./ zimp .^ 4 ;
taui_z2     =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* amax)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nz2p ./ lnz2 ./ zmax .^ 4 ;
taui_w      =  12 .* pi .^ (3/2) .* phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.mp .* 183.84)  .* ...
	                   (phys.e .* profli.tip) .^ (3/2) ./ nwp ./ lnw ./ max(1,z0wavez(profli.tep)) .^ 4 ;

% Plasma rotation in Tokamaks, V. Rozhansky and M. Tendler,Review of plasma physics, tome 19, p 163
vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 1 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_h;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkh     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 2 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi./ taui_d;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkd     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_t;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkt     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 4 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkhe     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 3 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_he;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkhe3     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ aimp ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz1     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ amax ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_z1;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkz2     = max(-2.1,min(1.7,fk));

vthi   = sqrt(2 .* phys.e .* profli.tip  ./ 183.84 ./ phys.mp);
nuis   = profli.Raxe .* profli.qjli .* profli.epsi .^ (-3/2) ./ vthi ./ taui_w;
nuis(:,1) =  2 .* nuis(:,2) - nuis(:,3);
fk     = - (1.17 - 0.35 .* sqrt(nuis) - 2.1 .* nuis .^ 2 .* profli.epsi .^ 2) ./ ...
	            (1 + 0.7 .* sqrt(nuis) + nuis .^ 2 .* profli.epsi .^ 3);
fkw     = max(-2.1,min(1.7,fk));

warning on



% changement de repere
dpsidrho = abs(pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve));
dpsidrho(:,1) = NaN;

% vitesse poloidal 
b2 = (profli.fdia .* profli.ri) .^ 2; 
% standard  + orbit squeezing
gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .*nHp,0,2,2,1) ./ nHp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_h  = - fkh .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_h(:,1) = 2 .* utheta_h(:,2) - utheta_h(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nDp,0,2,2,1) ./ nDp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_d  = - fkd .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_d(:,1) = 2 .* utheta_d(:,2) - utheta_d(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nTp,0,2,2,1) ./ nTp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_t  = - fkt .* option.signe .* profli.fdia ./ b2 ./ phys.e .* gtheta ./ dpsidrho; 
utheta_t(:,1) = 2 .* utheta_t(:,2) - utheta_t(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep,0,2,2,1) ./ nhep;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_he  = - fkhe .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho; 
utheta_he(:,1) = 2 .* utheta_he(:,2) - utheta_he(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nhep3,0,2,2,1) ./ nhep;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_he3  = - fkhe3 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ 2 .* gtheta ./ dpsidrho; 
utheta_he3(:,1) = 2 .* utheta_he3(:,2) - utheta_he3(:,3);

gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
gtheta2 = pdederive(x,phys.e .* profli.tip .* nwp,0,2,2,1) ./ nwp;
stheta  = min(sign(gtheta1),sign(gtheta2));
gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
utheta_w  = - fkw .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ max(1,z0wavez(profli.tep)) .* gtheta ./ dpsidrho; 
utheta_w(:,1) = 2 .* utheta_w(:,2) - utheta_w(:,3);

%  gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
%  gtheta2 = pdederive(x,phys.e .* profli.tip .* nz1p,0,2,2,1) ./ nz1p;
%  stheta  = min(sign(gtheta1),sign(gtheta2));
%  gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
%  utheta_z1  = - fkz1 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zimp .* gtheta ./ dpsidrho; 
%  utheta_z1(:,1) = 2 .* utheta_z1(:,2) - utheta_z1(:,3);
%  
%  gtheta1 = pdederive(x,phys.e .* profli.tip,0,2,2,1);
%  gtheta2 = pdederive(x,phys.e .* profli.tip .* nz2p,0,2,2,1) ./ nz2p;
%  stheta  = min(sign(gtheta1),sign(gtheta2));
%  gtheta  = stheta .* max(abs(gtheta1),abs(gtheta2)) ./ (profli.rmx(:,end) * ve);
%  utheta_z2  = - fkz2 .* option.signe .* profli.fdia ./ b2 ./ phys.e ./ zmax .* gtheta ./ dpsidrho; 
%  utheta_z2(:,1) = 2 .* utheta_z2(:,2) - utheta_z2(:,3);

% les champ en Rmax
a = interp1_itm(cons.temps,geo.a,profli.temps,'pchip','extrap');
rmax         = profli.Raxe + a * x;
btor         = option.signe .* (profli.fdia ./rmax);
grho         = abs((profli.rmx(:,end) * ve) ./ max(eps,pdederive(x,rmax,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,profli.psi,0,2,2,1)./ rmax .* grho ./ (profli.rmx(:,end) * ve);
btot         = sqrt(btor .^ 2 + bpol .^ 2);

switch option.mino
case 'He3'
	ag = 3;
	zg = 2;
case 'T'
	ag = 3;
	zg = 1;
case 'He4'
	ag = 4;   
	zg = 2;
case 'D'
	ag = 2;   
	zg = 1;
otherwise
	ag = 1;   
	zg = 1;
end

switch option.gaz
case 1
   zj = 1;
   aj = 1;
case 2
   zj = 1;
   aj = 2;
case 3
   zj = 1;
   aj = mean((2  + 3 .* real(cons.iso)) ./ (1 +  real(cons.iso)));
case 4
   zj = 2;
   aj = 4;
end


% calcul de la rotation poloidal pour l'impurete principale
% Y. B. Kim et all Phys. Fluids. B 3  (8) 1991  p 2050-
switch option.gaz
case 4
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.nhep) .* 4);
	nii   = max(1e13,profli.nhep);
otherwise
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* 1); 
	nii   = max(1e13,profli.n1p);
end
% ontraite les imurepte a l'etat de trace ...
alpha  = min(zimp,alpha);
beta   = (27/4) .^ 2 .* (aj ./ aimp) .^ 2 ./ (15/2 + sqrt(2 .* alpha) .* sqrt(aimp ./ aj));
g      = profli.ftrap ./ max(0.01,1 - profli.ftrap);
mui_00 = g .* (          alpha + sqrt(2)             -           log(1 + sqrt(2)));
mui_10 = g .* ((3/2)  .* alpha + 4  ./ sqrt(2)       - (5/2)  .* log(1 + sqrt(2)));
mui_01 = mui_10;
mui_11 = g .* ((13/4) .* alpha + 39 ./ (4 * sqrt(2)) - (25/4) .* log(1 + sqrt(2)));
D      = mui_00 .* (mui_11 + sqrt(2) + alpha - alpha .* beta) - mui_01 .^ 2;
D(D==0) = 1e38;
K1     = mui_01 ./ D .* (sqrt(2) + alpha - alpha .* beta);
K2     = (mui_00 .* mui_11 - mui_01 .* mui_10) ./ D;
vth    = sqrt(2 .* profli.tip .* phys.e ./ (phys.mp .* aj));
b2     = sqrt(profli.bpol .^ 2 + (profli.fdia .* profli.ri) .^ 2);
rhoi   = 4.57e-3 .* sqrt(aj .* profli.tip ./ 1e3 ./ b2);
ltim1  = pdederive(x,profli.tip,0,2,2,1) ./ profli.tip ./ (profli.rmx(:,end) * ve); 
lpiim1 = pdederive(x,profli.tip .* nii,0,2,2,1)  ./ (profli.tip .* nii)  ./ (profli.rmx(:,end) * ve); 
lpiIm1 = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.tip .* nz1p) ./ (profli.rmx(:,end) * ve); 
%figure(21);plot(x,1./ltim1,'b',x,1./lpiim1,'r',x,1./lpiIm1,'g');drawnow
%
vtehta_z1 = option.signe .*  0.5 .* vth .* rhoi .* ((K1 + (3/2) .* K2) .* ltim1 - lpiim1 + (zj ./ zimp) .* 1 .* lpiIm1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_z1 = vtehta_z1 ./ max(eps,profli.bpol);
utheta_z1(:,1) = 0;

% calcul de la rotation poloidal pour l'impurete principale
% Y. B. Kim et all Phys. Fluids. B 3  (8) 1991  p 2050-
switch option.gaz
case 4
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.nhep) .* 4);
	nii   = max(1e13,profli.nhep);
otherwise
	alpha = (nz1p .* zimp .^ 2) ./ (max(1e13,profli.n1p) .* 1); 
	nii   = max(1e13,profli.n1p);
end
% ontraite les imurepte a l'etat de trace ...
alpha  = min(zmax,alpha);
beta   = (27/4) .^ 2 .* (aj ./ aimp) .^ 2 ./ (15/2 + sqrt(2 .* alpha) .* sqrt(amax ./ aj));
g      = profli.ftrap ./ max(0.01,1 - profli.ftrap);
mui_00 = g .* (          alpha + sqrt(2)             -           log(1 + sqrt(2)));
mui_10 = g .* ((3/2)  .* alpha + 4  ./ sqrt(2)       - (5/2)  .* log(1 + sqrt(2)));
mui_01 = mui_10;
mui_11 = g .* ((13/4) .* alpha + 39 ./ (4 * sqrt(2)) - (25/4) .* log(1 + sqrt(2)));
D      = mui_00 .* (mui_11 + sqrt(2) + alpha - alpha .* beta) - mui_01 .^ 2;
D(D==0) = 1e38;
K1     = mui_01 ./ D .* (sqrt(2) + alpha - alpha .* beta);
K2     = (mui_00 .* mui_11 - mui_01 .* mui_10) ./ D;
vth    = sqrt(2 .* profli.tip .* phys.e ./ (phys.mp .* aj));
b2     = sqrt(profli.bpol .^ 2 + (profli.fdia .* profli.ri) .^ 2);
rhoi   = 4.57e-3 .* sqrt(aj .* profli.tip ./ 1e3 ./ b2);
ltim1  = pdederive(x,profli.tip,0,2,2,1) ./ profli.tip ./ (profli.rmx(:,end) * ve); 
lpiim1 = pdederive(x,profli.tip .* nii,0,2,2,1)  ./ (profli.tip .* nii)  ./ (profli.rmx(:,end) * ve); 
lpiIm1 = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.tip .* nz1p) ./ (profli.rmx(:,end) * ve); 
%figure(21);plot(x,1./ltim1,'b',x,1./lpiim1,'r',x,1./lpiIm1,'g');drawnow
%
vtehta_z2 = option.signe .*  0.5 .* vth .* rhoi .* ((K1 + (3/2) .* K2) .* ltim1 - lpiim1 + (zj ./ zimp) .* 1 .* lpiIm1) .* (profli.fdia .* profli.ri) ./ sqrt(b2);
utheta_z2 = vtehta_z2 ./ max(eps,profli.bpol);
utheta_z2(:,1) = 0;



% changement de repere
dpsidrho =  pdederive(x,profli.psi,0,2,2,1) ./ (profli.rmx(:,end) * ve);
dpsidrho(:,1) = 0;
% calcul du champ electrique radial (Er gradient(rho))
Ptor    = phys.mp .* pdederive(x,profli.tip .* nHp,0,2,2,1) + ...
          2 .* phys.mp .* pdederive(x,profli.tip .* nDp,0,2,2,1) + ...
          3 .* phys.mp .* pdederive(x,profli.tip .* nTp,0,2,2,1) + ...
          2 .* phys.mp .* pdederive(x,profli.tip .* nhep,0,2,2,1) + ...
          aimp ./ zimp .* phys.mp .* pdederive(x,profli.tip .* nz1p,0,2,2,1) + ...
          amax ./ zmax .* phys.mp .* pdederive(x,profli.tip .* nz2p,0,2,2,1) + ...
          183.84  ./ max(1,z0wavez(profli.tep)) .* phys.mp .* pdederive(x,profli.tip .* nwp,0,2,2,1);
Ptor    = Ptor ./ (profli.rmx(:,end) * ve);			



% sorties
vtor         = ones(length(profli.temps),length(profli.xli),8);
vpol         = ones(length(profli.temps),length(profli.xli),8);

% calul de la vitessse toroidal en Rmax pour l'impurete principale
vtheta_z1   = utheta_z1 .* bpol;
vtheta_z2   = utheta_z2 .* bpol;
vtheta_he   = utheta_he  .* bpol;
vtheta_he3  = utheta_he3 .* bpol;
vtheta_h    = utheta_h .* bpol;
vtheta_d    = utheta_d .* bpol;
vtheta_t    = utheta_t .* bpol;
vtheta_w    = utheta_w .* bpol;


% calculde la roration toroidale
warning off
inter          = pdederive(x,profli.tip .* nz1p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zimp .* nz1p);
omega_z1       = (profli.er - inter) ./ dpsidrho; 
vtor_z1        = omega_z1 .* rmax + utheta_z1 .* option.signe .* profli.fdia ./ rmax;  
vtor_z1(:,1) = 2 .* vtor_z1(:,2) - vtor_z1(:,3);

inter          = pdederive(x,profli.tip .* nz2p,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (zmax .* nz2p);
omega_z2       = (profli.er - inter) ./ dpsidrho; 
vtor_z2        = omega_z2 .* rmax + utheta_z2 .* option.signe .* profli.fdia ./ rmax;  
vtor_z2(:,1) = 2 .* vtor_z2(:,2) - vtor_z2(:,3);

inter          = pdederive(x,profli.tip .* nwp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (max(1,z0wavez(profli.tep)) .* nwp);
omega_w       = (profli.er - inter) ./ dpsidrho; 
vtor_w        = omega_w .* rmax + utheta_w .* option.signe .* profli.fdia ./ rmax;  
vtor_w(:,1) = 2 .* vtor_w(:,2) - vtor_w(:,3);

inter          = pdederive(x,profli.tip .* nhep,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep);
omega_he       = (profli.er - inter) ./ dpsidrho; 
vtor_he        = omega_he .* rmax + utheta_he .* option.signe .* profli.fdia ./ rmax;  
vtor_he(:,1) = 2 .* vtor_he(:,2) - vtor_he(:,3);

inter          = pdederive(x,profli.tip .* nhep3,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (2 .* nhep3);
omega_he3       = (profli.er - inter) ./ dpsidrho; 
vtor_he3        = omega_he3 .* rmax + utheta_he3 .* option.signe .* profli.fdia ./ rmax;  
vtor_he3(:,1) = 2 .* vtor_he3(:,2) - vtor_he3(:,3);

inter          = pdederive(x,profli.tip .* nHp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nHp);
omega_H       = (profli.er - inter) ./ dpsidrho; 
vtor_H        = omega_H .* rmax + utheta_h .* option.signe .* profli.fdia ./ rmax;  
vtor_H(:,1) = 2 .* vtor_H(:,2) - vtor_H(:,3);

inter          = pdederive(x,profli.tip .* nDp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nDp);
omega_D       = (profli.er - inter) ./ dpsidrho; 
vtor_D        = omega_D .* rmax + utheta_d .* option.signe .* profli.fdia ./ rmax;  
vtor_D(:,1) = 2 .* vtor_D(:,2) - vtor_D(:,3);


inter          = pdederive(x,profli.tip .* nTp,0,2,2,1) ./ (profli.rmx(:,end) * ve) ./ (1 .* nTp);
omega_T       = (profli.er - inter) ./ dpsidrho; 
vtor_T        = omega_T .* rmax + utheta_t .* option.signe .* profli.fdia ./ rmax;  
vtor_T(:,1) = 2 .* vtor_T(:,2) - vtor_T(:,3);
		 
		 
		 
warning on
		 
%
% mise en forme
%
vtor(:,:,1) = vtor_H;
vtor(:,:,2) = vtor_D;
vtor(:,:,3) = vtor_T;
vtor(:,:,4) = vtor_he3;
vtor(:,:,5) = vtor_he;
vtor(:,:,6) = vtor_z1;
vtor(:,:,7) = vtor_z2;
vtor(:,:,8) = vtor_w;

vpol(:,:,1) = vtheta_h;
vpol(:,:,2) = vtheta_d;
vpol(:,:,3) = vtheta_t;
vpol(:,:,4) = vtheta_he3;
vpol(:,:,5) = vtheta_he;
vpol(:,:,6) = vtheta_z1;
vpol(:,:,7) = vtheta_z2;
vpol(:,:,8) = vtheta_w;

omega(:,:,1) = omega_H;
omega(:,:,2) = omega_D;
omega(:,:,3) = omega_T;
omega(:,:,4) = omega_he3;
omega(:,:,5) = omega_he;
omega(:,:,6) = omega_z1;
omega(:,:,7) = omega_z2;
omega(:,:,8) = omega_w;

iso   = interp1_itm(cons.temps,real(cons.iso),profli.temps,'pchip','extrap');
ftnbi = interp1_itm(cons.temps,cons.ftnbi,profli.temps,'pchip','extrap');
mtor_H = 0 .* omega_H;
mtor_D = 0 .* omega_D;
mtor_T = 0 .* omega_T;
mtor_he = 0 .* omega_he;
switch option.gaz
case 1
   mtor_H = mtor_H + profli.rot_n0;   
case 2
   mtor_D = mtor_D + profli.rot_n0;
case 3
   mtor_D = mtor_D + profli.rot_n0 ./ (1 + iso * ones(1,size(mtor_D,2)));
   mtor_T = mtor_T + profli.rot_n0 .* (iso * ones(1,size(mtor_T,2))) ./ (1 + iso * ones(1,size(mtor_D,2)));
case 4
   mtor_he = mtor_he + profli.rot_n0;  
end

switch option.gaz
case 3
    mtor_D = mtor_D + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_D,2)));
    mtor_T = mtor_T + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_T,2)));
otherwise
    mtor_D = mtor_D + profli.rot_nbi .* (1 - ftnbi * ones(1,size(mtor_D,2)));
    mtor_H = mtor_H + profli.rot_nbi .* (ftnbi * ones(1,size(mtor_H,2)));    
end

mtor(:,:,1) = mtor_H;
mtor(:,:,2) = mtor_D;
mtor(:,:,3) = mtor_T;
mtor(:,:,4) = 0;
mtor(:,:,5) = mtor_he;
mtor(:,:,6) = 0;
mtor(:,:,7) = 0;
mtor(:,:,8) = 0;


% fonction d'auto test de metis4itm
function metis4itmautotest(shot,run)

%[error_flag,xsd] = metis4itm(shot,run,occurrence,methode,time,codeparam_filename,interpolation_methode)
% simple test
disp('simple test, creation of xml and xsd')
metis4itm;
disp('simple test, run of test case')
[error_flag,output] = metis4itm(shot,run,'','test');
% re lecture
[error_flag,output] = metis4itm(shot,run,'','read');
% test calcul complet avec entree
disp('run in fast mode')
metis4itm(shot,run,'','fast')
disp('run in full mode')
metis4itm(shot,run,'','full')
disp('initialisation evolution mode')
metis4itm(shot,run,'','init',1);
for k=2:10
	disp('one time evolution mode')
	k
	% k is time
	metis4itm(shot,run,'','one_time',k);
end
% save with restart
disp('test restart');
option.restart='test_retsart'
metis4itm(shot,run,'','one_time',k+1,option);
metis4itm(shot,run,'','test_retsart',k-2);
metis4itm(shot,run,'','one_time',k-1);
metis4itm(shot,run,'','one_time',k);
% same with reset cpo
disp('test init output CPOs');
option.init_output_cpo =1;
metis4itm(shot,run,'','init',1,option);
for k=2:10
	disp('one time evolution mode')
	k
	% k is time
	metis4itm(shot,run,'','one_time',k,option);
end




function [machine,user,ver] = which_MDSdatabase_metis
%from LUKE - Function that check which ITM MDS+ database is connected
%
%Function that check which ITM MDS+ database is connected
%
%by Y.Peysson CEA-IRFM <yves.peysson@cea.fr> and Joan Decker CEA-IRFM (joan.decker@cea.fr)
%
%
%
hdf5base = getenv('HDF5_BASE');
mdstreebase = getenv('MDSPLUS_TREE_BASE_0');
%
if isempty(hdf5base) | isempty(mdstreebase),
	warning('Environment variables for ITM MDS+/HDF5 database are not set-up.');
	machine = '';
	user = '';
	ver = '';
	return
end
%
is = 1;
str = {};
remain = hdf5base;
while true
	[strtemp,remain] = strtok(remain,'/');	
	if isempty(strtemp)
		break
	end
	str{is} = strtemp;
	is = is + 1;
end
%
if strcmp(str{1},'afs'),%Private MDS+ database (may be own by any ITM user)
	ver = str{length(str)-1};
	machine = str{length(str)-2};
	user = str{length(str)-6};
else,%Public MDS+ database
	ver = str{length(str)-1};
	machine = str{length(str)-2};
	user = str{length(str)-3};
end	



% other cpos
% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coretransp = mapcoretransp(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coretransp,sigma_B0_eff)

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



% coretransp.datainfo est le meme que celui de scenario
if ~isfield(coretransp,'datainfo')
	coretransp.datainfo = datainfo_empty;
end

coretransp.time			= profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coretransp.rho_tor_norm 		= xli;
coretransp.rho_tor       		= profil0d.rmx;
%  if length(profil0d.temps) == 1
%  	coretransp.drho_dt       		  = data_zerod.drmdt ./ profil0d.rmx(end) .*profil0d.rmx;
%  else 
%  	coretransp.drho_dt       		  = zdxdt(profil0d.rmx,profil0d.temps);
%  end
%  rb0 =   interp1_itm(data_zerod.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
%  coretransp.toroid_field.time        = profil0d.temps;
%  coretransp.toroid_field.r0          = mean(z0dstruct.z0dinput.geo.R);
%  coretransp.toroid_field.b0          = rb0 ./ coretransp.toroid_field.r0;
%  coretransp.toroid_field.b0prime     = z0dxdt(coretransp.toroid_field.b0,coretransp.toroid_field.time);
%
coretransp.composition.amn 		= [1,2,3,3,4,ceil(7/3 .* z0dstruct.z0dinput.option.zimp),ceil(7/3 .* z0dstruct.z0dinput.option.zmax)];
coretransp.composition.zn	 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax];
coretransp.composition.zion 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax];
coretransp.composition.imp_flag 		= zeros(size(coretransp.composition.zion));
coretransp.compositions = copy_composition_to_compositionstype(coretransp.composition);
%
coretransp.values(1).sigma = 1./max(1e-307,profil0d.eta);
%
coretransp.values(1).te_transp.diff_eff = profil0d.xie;
coretransp.values(1).te_transp.flux     = profil0d.qe;
coretransp.values(1).te_transp.flag     = 1;
% un seul Ti pour ni = sum(nj)
coretransp.values(1).ti_transp.diff_eff = profil0d.xie;
coretransp.values(1).ti_transp.flux     = profil0d.qi;
coretransp.values(1).values(1).ti_transp.exchange = - profil0d.qei;
coretransp.values(1).ti_transp.flag     = 1;
%
coretransp.values(1).ne_transp.diff_eff  = cat(3,profil0d.dn,1.5 .* profil0d.dn,2.5 .* profil0d.dn);
coretransp.values(1).ne_transp.vconv_eff = -cat(3,profil0d.vn,1.5 .* profil0d.vn,2.5 .* profil0d.vn);
coretransp.values(1).ne_transp.flux      = profil0d.ge;
coretransp.values(1).ne_transp.flag      = 2;
% un seul Rtor 

ve = ones(size(profil0d.xli));
nDp   = max(1e13,profil0d.n1p .* ((interp1_itm(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap')./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nTp   = max(1e13,profil0d.n1p .* ((interp1_itm(data_zerod.temps,data_zerod.nDm,profil0d.temps,'pchip','extrap')./ max(1,trapz(profil0d.xli,profil0d.vpr .* abs(profil0d.n1p),2)) .* trapz(profil0d.xli,profil0d.vpr,2)) * ve));
nHp   = max(1e13,profil0d.n1p - nTp - nDp);
nhep  = max(1e13,profil0d.nhep);
nz1p  = max(1e13,profil0d.nzp);
nz2p  = max(1e11,profil0d.nzp .* z0dstruct.z0dinput.option.rimp);   
nwp  = max(1,profil0d.nwp);
% une seule vitesse d'ensemble
Mtor    = phys.mp .*  max(1e13,nHp +  2 .* nDp + 3 .* nTp + 4 .* nhep + ...
          ceil(7/3 .* z0dstruct.z0dinput.option.zimp) .* nz1p + ceil(7/3 .* z0dstruct.z0dinput.option.zmax) .* nz2p + 183.84 .* nwp);

coretransp.values(1).vtor_transp.diff_eff  = profil0d.drot;
coretransp.values(1).vtor_transp.vconv_eff = -profil0d.vrot;
coretransp.values(1).vtor_transp.flux      = profil0d.frot ./ Mtor .* profil0d.r2i .* profil0d.Raxe;
coretransp.values(1).vtor_transp.flag        = 2;



% mapping du CPO coreprof (seule les variables d'interet sont remplies, un model vide doit etre fourni)
function coresource = mapcoresource_generic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coresource,sigma_B0_eff)

% coresource.datainfo est le meme que celui de scenario
if ~isfield(coresource,'datainfo')
	coresource.datainfo = datainfo_empty;
end
coresource.datainfo.comment = 'METIS sources';
coresource.time			= profil0d.temps;
%  if length(profil0d.temps) == 1
%  	coresource.drho_dt       		  = data_zerod.drmdt ./ profil0d.rmx(end) .*profil0d.rmx;
%  else 
%  	coresource.drho_dt       		  = zdxdt(profil0d.rmx,profil0d.temps);
%  end
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
coresource.toroid_field.r0          = mean(z0dstruct.z0dinput.geo.R);
coresource.toroid_field.b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ coresource.toroid_field.r0;
%
coresource.composition.amn 		= [1,2,3,3,4,ceil(7/3 .* z0dstruct.z0dinput.option.zimp),ceil(7/3 .* z0dstruct.z0dinput.option.zmax)];
coresource.composition.zn	 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax];
coresource.composition.zion 		= [1,1,1,2,2,z0dstruct.z0dinput.option.zimp,z0dstruct.z0dinput.option.zmax];
coresource.composition.imp_flag 		= zeros(size(coresource.composition.zion));
coresource.compositions = copy_composition_to_compositionstype(coresource.composition);

%

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_lhcd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'LHCD';
coresource.sourceid.id = 'lh';
coresource.sourceid.flag = 3;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  		= profil0d.jlh  .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));

coresource.qe.exp 	= profil0d.plh;
coresource.se.exp 	= 0.* profil0d.plh;
coresource.qi.exp 	= 0.* profil0d.plh;
coresource.ui.exp	= profil0d.rot_lh ; % contient aussi la contribution de eccd 

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_eccd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'ECCD';
coresource.sourceid.id = 'ec';
coresource.sourceid.flag = 2;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  		= profil0d.jeccd .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= profil0d.pecrh;
coresource.se.exp 	= 0.* profil0d.pecrh;
coresource.qi.exp 	= 0.* profil0d.pecrh;
coresource.ui.exp 	= 0 .* profil0d.rot_lh ; % contenu dans LH

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_radiation(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'Line radiation';
coresource.sourceid.id = 'lineradiation';
coresource.sourceid.flag = 19;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

coresource.qe.exp 	= profil0d.prad;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_brem(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'Brehmstrahlung';
coresource.sourceid.id = 'brehmstrahlung';
coresource.sourceid.flag = 15;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

coresource.qe.exp 	= profil0d.pbrem;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_cyclotron(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'Cyclotron radiation';
coresource.sourceid.id = 'cyclotronradiation';
coresource.sourceid.flag = 16;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

coresource.qe.exp 	= profil0d.pcyclo;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_neutral(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'Cold neutrals';
coresource.sourceid.id = 'coldneutralcooling';
coresource.sourceid.flag = 24;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
coresource.j  	= 0 .* profil0d.pioniz;
coresource.qe.exp 	= profil0d.pioniz;
coresource.se.exp 	= profil0d.s0 + profil0d.s0m + profil0d.spellet;
coresource.qi.exp 	= 0.* profil0d.pioniz;
coresource.ui.exp 	= profil0d.rot_n0;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_nbicd(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'NBICD';
coresource.sourceid.id = 'nbi';
coresource.sourceid.flag = 1;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  	= profil0d.jnbicd .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= max(0,profil0d.pnbi - profil0d.pnbi_ion);
coresource.se.exp 	= real(profil0d.nbinesource) + imag(profil0d.nbinesource);
coresource.qi.exp 	= profil0d.pnbi_ion;
coresource.ui.exp 	= profil0d.rot_nbi;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_icrh(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'ICRH';
coresource.sourceid.id = 'ic';
coresource.sourceid.flag = 4;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  	= profil0d.jfwcd .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= profil0d.pfweh + max(0,profil0d.picrh - profil0d.picrh_ion);
coresource.se.exp 	= 0 .* profil0d.nbinesource;
coresource.qi.exp 	= profil0d.picrh_ion;
coresource.ui.exp 	= 0 .* profil0d.rot_nbi;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_fusion(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'Fusion';
coresource.sourceid.id = 'fusion';
coresource.sourceid.flag = 5;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  	= profil0d.jfus .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= max(0,profil0d.pfus - profil0d.pfus_ion);
coresource.se.exp 	= 0 .* profil0d.nbinesource;
coresource.qi.exp 	= profil0d.pfus_ion;
coresource.si.exp 	= zeros(size(profil0d.pfus,1),size(profil0d.pfus,2),8);
coresource.si.exp(:,:,5) = profil0d.salf;
coresource.ui.exp 	= 0 .* profil0d.rot_nbi;

% continuer ici ---->>>>>>>

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_ohm_boot(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'ohmic & bootstrap';
coresource.sourceid.id = 'ohmic';
coresource.sourceid.flag = 14;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
%coresource.j  	= (profil0d.jrun + profil0d.jboot) .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.j  	= profil0d.jboot .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= profil0d.pohm;
coresource.se.exp 	= 0 .* profil0d.pohm;
coresource.qi.exp 	= 0 .* profil0d.pohm;
coresource.ui.exp 	= 0 .* profil0d.pohm;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_full(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'sum of all source terms';
coresource.sourceid.id = 'unspecified';
coresource.sourceid.flag = 0;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  	= profil0d.jni .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= profil0d.source_el;
coresource.se.exp 	= real(profil0d.nbinesource) + imag(profil0d.nbinesource) + profil0d.s0 + profil0d.s0m + profil0d.spellet;
coresource.qi.exp 	= profil0d.source_ion;
coresource.si.exp 	= zeros(size(profil0d.pfus,1),size(profil0d.pfus,2),8);
coresource.si.exp(:,:,5) = profil0d.salf;
coresource.ui.exp 	= profil0d.rot_nbi + profil0d.rot_n0 + profil0d.rot_lh;

% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coresource = mapcoresource_runaways(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,sigma_B0_eff)

coresource.sourceid.description = 'runaway electron source';
coresource.sourceid.id = 'runaways';
coresource.sourceid.flag = 34;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coresource.rho_tor_norm 		= xli;
coresource.rho_tor       		= profil0d.rmx;

zz              = zeros(size(profil0d.rmx));
coresource.j  	= zz;
coresource.qe.exp 	= zz;
coresource.se.exp 	= zz;
coresource.qi.exp 	= zz;
coresource.si.exp 	= zeros(size(profil0d.rmx,1),size(profil0d.rmx,2),8);
coresource.ui.exp 	= zz;
coresource.ujxb.exp 	= zz;
coresource.sigma 	= zz;

%
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
toroid_field.r0  = mean(z0dstruct.z0dinput.geo.R);
%toroid_field.b0  = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ toroid_field.r0;
% les courants ont déja le bon signe
toroid_field.b0  = rb0 ./ toroid_field.r0;
coresource.j  		= profil0d.jrun .* ((toroid_field.b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
coresource.qe.exp 	= 0.* profil0d.jrun;
coresource.se.exp 	= 0.* profil0d.jrun;
coresource.qi.exp 	= 0.* profil0d.jrun;
coresource.ui.exp 	= 0 .* profil0d.rot_lh ; % contenu dans LH



% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function coreneutrals = mapcoreneutrals(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,coreneutrals,sigma_B0_eff)

% coreneutrals.datainfo est le meme que celui de scenario
if ~isfield(coreneutrals,'datainfo')
	coreneutrals.datainfo = datainfo_empty;
end
coreneutrals.time			= profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
coreneutrals.rho_tor_norm 		= xli;
coreneutrals.rho_tor       		= profil0d.rmx;

switch z0dstruct.z0dinput.option.gaz
case 1
	coreneutrals.composition.atomlist.amn = 1;
	coreneutrals.composition.atomlist.zn  = 1;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
case 2
	coreneutrals.composition.atomlist.amn = 2;
	coreneutrals.composition.atomlist.zn  = 1;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
case 3
	coreneutrals.composition.atomlist.amn = [2 3];
	coreneutrals.composition.atomlist.zn  = [1 1];
	coreneutrals.composition.neutallist.ncomp = [1 1]; 
	coreneutrals.composition.neutallist.tatm = [1 2];
	coreneutrals.composition.neutallist.multatm = [1 1];
	coreneutrals.composition.typelist.ntype  = [2 2]; 
	coreneutrals.composition.typelist.type = [0 1;0 1]; 
case 4
	coreneutrals.composition.atomlist.amn = 2;
	coreneutrals.composition.atomlist.zn  = 4;
	coreneutrals.composition.neutallist.ncomp = 1; 
	coreneutrals.composition.neutallist.tatm = 1;
	coreneutrals.composition.neutallist.multatm = 1; 
	coreneutrals.composition.typelist.ntype  = 2; 
	coreneutrals.composition.typelist.type = [0 1]; 
otherwise 
	error('plasma compostion not yet implemented in METIITM');
end
compo_void = coreneutrals.composition.atomlist;
compo_void.zion = compo_void.zn;
coreneutrals.compositions = copy_composition_to_compositionstype(compo_void);

%
n0a = interp1_itm(data_zerod.temps,data_zerod.n0a,profil0d.temps,'pchip','extrap');
%
switch z0dstruct.z0dinput.option.gaz
case 3
	iso = interp1_itm(z0dstruct.z0dinput.cons.temps,real(z0dstruct.z0dinput.cons.iso),profil0d.temps,'pchip','extrap');
	iso = iso * ones(size(profil0d.xli));
	coreneutrals.profiles.n0.value            = NaN .* ones(size(profil0d.n0,1),size(profil0d.n0,2),2,2);
	coreneutrals.profiles.n0.value(:,:,1,1)   = reshape(profil0d.n0m ./ (1+iso),size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,1,2)   = reshape(profil0d.n0  ./ (1+iso),size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,2,1)   = reshape(profil0d.n0m ./ (1+iso) .* iso,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,2,2)   = reshape(profil0d.n0  ./ (1+iso) .* iso,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.flux             = NaN .* coreneutrals.profiles.n0.value;
	coreneutrals.profiles.n0.flux(:,:,1,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,1,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,2,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso) .* iso, ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,2,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho) ./ (1+iso) .* iso, ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.boundary.value   = zeros(size(profil0d.n0,1),3,2,2);
	coreneutrals.profiles.n0.boundary.value(:,1,1,1)   = n0a ./ (1+iso(:,1));
	coreneutrals.profiles.n0.boundary.value(:,1,1,2)   =0;
	coreneutrals.profiles.n0.boundary.value(:,1,2,1)   = n0a ./ (1+iso(:,1)) .* iso(:,1);
	coreneutrals.profiles.n0.boundary.value(:,1,2,2)   =0;
	coreneutrals.profiles.n0.boundary.type    = NaN .* ones(size(profil0d.n0,1),2,2);
	coreneutrals.profiles.n0.boundary.type(:)    = 4;  
	coreneutrals.profiles.n0.boundary.rho_tor = NaN .* ones(size(profil0d.n0,1),2,2);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,2) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,2,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,2,2) = coreneutrals.rho_tor(:,end);
otherwise
	coreneutrals.profiles.n0.value            = NaN .* ones(size(profil0d.n0,1),size(profil0d.n0,2),1,2);
	coreneutrals.profiles.n0.value(:,:,1,1)   = reshape(profil0d.n0m,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.value(:,:,1,2)   = reshape(profil0d.n0,size(coreneutrals.profiles.n0.value(:,:,1,1)));
	coreneutrals.profiles.n0.flux             = NaN .* coreneutrals.profiles.n0.value;
	coreneutrals.profiles.n0.flux(:,:,1,1)    = reshape(cumtrapz(profil0d.xli,profil0d.s0m .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.flux(:,:,1,2)    = reshape(cumtrapz(profil0d.xli,profil0d.s0 .* profil0d.vpr,2) ./ ...
                                            		max(eps,profil0d.vpr_tor .* profil0d.grho), ...
							size(coreneutrals.profiles.n0.flux(:,:,1,1)));
	coreneutrals.profiles.n0.boundary.value   = zeros(size(profil0d.n0,1),3,1,2);
	coreneutrals.profiles.n0.boundary.value(:,1,1,1)   =n0a;
	coreneutrals.profiles.n0.boundary.value(:,1,1,2)   =0;
	coreneutrals.profiles.n0.boundary.type    = NaN .* ones(size(profil0d.n0,1),1,2);
	coreneutrals.profiles.n0.boundary.type(:)    = 4;  
	coreneutrals.profiles.n0.boundary.rho_tor = NaN .* ones(size(profil0d.n0,1),1,2);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,1) = coreneutrals.rho_tor(:,end);
	coreneutrals.profiles.n0.boundary.rho_tor(:,1,2) = coreneutrals.rho_tor(:,end);
end


% mapping du CPO coreprof (seule les variables d'interret sont remplie, un model vide doit etre fourni)
function neoclassic = mapneoclassic(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,neoclassic,vtor,vpol,sigma_B0_eff)

% neoclassic.datainfo est le meme que celui de scenario
if ~isfield(neoclassic,'datainfo')
	neoclassic.datainfo = datainfo_empty;
end
neoclassic.time			        = profil0d.temps;
xli    = profil0d.rmx ./ (profil0d.rmx(:,end) * ones(1,size(profil0d.rmx,2)));
neoclassic.rho_tor_norm 		= xli;
neoclassic.rho_tor       		= profil0d.rmx;
rb0 =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R .* z0dstruct.z0dinput.geo.b0,profil0d.temps,'pchip','extrap');
r0  =   interp1_itm(z0dstruct.z0dinput.cons.temps,z0dstruct.z0dinput.geo.R ,profil0d.temps,'pchip','extrap');
r0  = mean(z0dstruct.z0dinput.geo.R);
%b0          = sigma_B0_eff .* z0dstruct.z0dinput.option.signe .* rb0 ./ r0;
% les courants ont déja le bon signe
b0  = rb0 ./ r0;
%
neoclassic.jboot = (profil0d.jrun + profil0d.jboot) .* ((b0 ./ (rb0 ./ r0)) * ones(size( profil0d.xli)));
neoclassic.sigma = 1./max(1e-307,profil0d.eta);
neoclassic.vpol  = vpol;
neoclassic.er    = profil0d.er;
neoclassic.ne_neo.vconv_eff = profil0d.ware;


% securite sur les occurrences de cpo
function option = secure_option(option)

option.scenario_occurrence = secure_occurrence(option.scenario_occurrence);
option.coreprof_occurrence = secure_occurrence(option.coreprof_occurrence);
option.coretransp_occurrence = secure_occurrence(option.coretransp_occurrence);
option.coreneutrals_occurrence = secure_occurrence(option.coreneutrals_occurrence);
option.neoclassic_occurrence = secure_occurrence(option.neoclassic_occurrence);
option.equilibrium_occurrence = secure_occurrence(option.equilibrium_occurrence);
option.coresources_occurrence = secure_occurrence(option.coresources_occurrence);

function occ_out = secure_occurrence(occ)

if isempty(occ)
	occ_out = '';
elseif isnumeric(occ)
	occ = fix(occ);
	if le(occ,0)
            occ_out ='';
        else
            occ_out = num2str(occ);
        end
else
	occ = fix(str2num(occ));
	if le(occ,0)
            occ_out ='';
        else
            occ_out = num2str(occ);
        end
end

function mat = metis_gaz(z0dstruct,data_zerod,profil0d)


% compatibilite
if ~isfield(z0dstruct,'profil')
   z0dstruct.profil =  z0dstruct.profil0d;
end
% script d'estimation des flux de matiere
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3
 
% compatibilite ascendante
if ~isfield(z0dstruct.zerod,'nbar_nat') || all(~isfinite(z0dstruct.zerod.nbar_nat))
     ulh = 0.25;
     fh  = z0dstruct.zerod.modeh;
     z0dstruct.zerod.nbar_nat   = min(z0dstruct.zerod.negr,max(1e13, 1e20 .* (z0dstruct.z0dinput.geo.b0 ./ z0dstruct.zerod.q95 ./ z0dstruct.z0dinput.geo.R) .^ 0.6  .*  (ulh + (1 - ulh) .* fh)));
     data_zerod.nbar_nat = inpterp1(z0dstruct.zerod.temps,z0dstruct.zerod.nbar_nat,data_zerod.temps);
     z0dstruct.zerod.frac_pellet = z0dstruct.zerod.frac_pellet .* (z0dstruct.zerod.frac_pellet > 1e-2);
     data_zerod.frac_pellet = inpterp1(z0dstruct.zerod.temps,z0dstruct.zerod.frac_pellet,data_zerod.temps);
end
pellet_fraction = data_zerod.frac_pellet;
flag_compute = 'new';

switch flag_compute

case 'new'

    % choix du coefficiet de recyclage
    eta_p = 0.5;
    filter_width = 11;
    Recycling = z0dstruct.z0dinput.option.Recycling;
    eta_g_default = 0.1;

    % temps
    mat.temps        = data_zerod.temps;
    % cette formule n'est pas precise numeriquement
    % flux de matiere sortant
    mat.output      = profil0d.ge(:,end) .* profil0d.grho2(:,end) .* profil0d.vpr_tor(:,end) ./ (4.41e-4 .* phys.avo);
    if length(data_zerod.temps) > 1
	  mat.output      = interp1(profil0d.temps,mat.output,data_zerod.temps,'linear','extrap');
    end
    % bilan contenu du plasma
    ntot              = trapz(z0dstruct.profil.xli,z0dstruct.profil.vpr .* z0dstruct.profil.nep,2) ./ (4.41e-4 .* phys.avo);
    mat.contents      = interp1(z0dstruct.profil.temps,ntot,mat.temps,'linear','extrap');
    mat.plasmanetflux = interp1(z0dstruct.profil.temps,z0dxdt(ntot,z0dstruct.profil.temps),data_zerod.temps,'linear','extrap'); % flux net qui sort du plasam

    % source de recylcage dans le plasma (a partir de la source de neutre qui entre dans le plasam)
    mat.s0_in       = data_zerod.n0a ./ (4.41e-4 .* phys.avo);  % recyclage  + gas puff

    % effect a recycling fraction (at LCFS not in the divertor)
    % sorte de minimum pour le flux de gaz qui circule dans la SOL
    switch z0dstruct.z0dinput.option.configuration
    case {0,1}
	mat.s0       = mat.s0_in  ./ max(eps,z0dstruct.z0dinput.option.fn0a); 
    case {2,3}
	mat.s0       = mat.s0_in  ./ max(eps,data_zerod.xpoint .* z0dstruct.z0dinput.option.fn0a_div + (~data_zerod.xpoint) .* z0dstruct.z0dinput.option.fn0a);
    otherwise
	mat.s0       = mat.s0_in  ./ max(eps,z0dstruct.z0dinput.option.fn0a_div); 
    end
    % s0_star est la source qui sort du plasma
    % computation of suggested fn0a
    fn0a_div = mat.output(data_zerod.xpoint ~= 0) ./ mat.s0(data_zerod.xpoint ~= 0);
    fn0a_div = mean(fn0a_div(isfinite(fn0a_div) & (fn0a_div>0)));
    fn0a     = mat.output(data_zerod.xpoint == 0) ./ mat.s0(data_zerod.xpoint == 0);
    fn0a     = mean(fn0a(isfinite(fn0a) & (fn0a>0)));
    try
	if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
		mat.s0          = medfilt1(mat.s0,filter_width);
		mat.s0_in       = medfilt1(mat.s0_in,filter_width);
		mat.output      = medfilt1(mat.output ,filter_width);
	end
    catch
	disp('No signal toolbox licence available : filtering removed in metis_gaz');
    end
    if length(data_zerod.temps) > 1
	mat.bilan_s0    = cumtrapz(data_zerod.temps,mat.s0);
	mat.bilan_s0_in = cumtrapz(data_zerod.temps,mat.s0_in);
    else
	mat.bilan_s0    = 0;
	mat.bilan_s0_in = 0;    
    end
    % security
    output_mem       = mat.output;
    mat.output       = max(mat.output,1.1 .* mat.s0);
    if length(data_zerod.temps) > 1
	  mat.bilan_output = cumtrapz(data_zerod.temps,mat.output);
    else
	  mat.bilan_output = 0;
    end
    % source due au glacon (dans le plasma)
    mat.pellet_star   = trapz(profil0d.xli,profil0d.spellet .* profil0d.vpr,2)  ./ (4.41e-4 .* phys.avo);
    mat.pellet_star(~isfinite(mat.pellet_star)) = 0;
    if (filter_width > 1) &&  (z0dstruct.z0dinput.option.pif < 1) && (length(data_zerod.temps) > filter_width)
        try
	    mat.pellet_star          = medfilt1(mat.pellet_star,filter_width);
	end 
    end
    if length(data_zerod.temps) > 1
	    mat.pellet_star          = interp1(profil0d.temps,mat.pellet_star,data_zerod.temps,'linear','extrap');
    end
    if length(data_zerod.temps) > 1
	    mat.bilan_pellet_star  = max(eps,cumtrapz(data_zerod.temps,mat.pellet_star));
    else
 	    mat.bilan_pellet_star  = 0;  
    end
    % flux de  matiere du a l'injection de gla�on
    mat.pellet = mat.pellet_star ./ eta_p;
    mat.bilan_pellet = mat.bilan_pellet_star ./ eta_p;

    % source due a l'injection de neutres
    mat.snbi_star        = real(data_zerod.pnbi) ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo) + ...
			  imag(data_zerod.pnbi) ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo);
    if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
        try
	    mat.snbi_star          = medfilt1(mat.snbi_star,filter_width);
	end
    end
    if length(data_zerod.temps) > 1
	    mat.bilan_nbi_star   = max(eps,cumtrapz(data_zerod.temps,mat.snbi_star));
    else
	    mat.bilan_nbi_star   = eps;  
    end
    % flux de matiere du a l'injection de neutres
    mat.snbi             = real(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(data_zerod.frnbi)) + ...
			  imag(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(data_zerod.frnbi));
    if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
         try
	    mat.snbi          = medfilt1(mat.snbi,filter_width);
	 end
    end
    mat.perte_snbi       = real(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo)  ./ max(eps,real(data_zerod.frnbi)) .* (1 - real(data_zerod.frnbi)) + ...
			  imag(data_zerod.pnbi)  ./ z0dstruct.z0dinput.option.einj2 ./ phys.e ./ (4.41e-4 .* phys.avo) ./ max(eps,imag(data_zerod.frnbi)) .* (1 - imag(data_zerod.frnbi));                      
    if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
         try
	    mat.perte_snbi          = medfilt1(mat.perte_snbi,filter_width);
	 end
    end
     if length(data_zerod.temps) > 1
	    mat.bilan_nbi   = max(eps,cumtrapz(data_zerod.temps,mat.snbi));
    else
	    mat.bilan_nbi   = eps;  
    end

    % calcul gaz puff (a partir du scaling pour la densite naturelle)
    fact_nbar = data_zerod.nem ./ max(1,data_zerod.nbar); 
    if (z0dstruct.z0dinput.option.tauhemul < 0) && (z0dstruct.z0dinput.option.Recycling < 1)
	    % modele qui prend en compte le confinement reel et le recyclage dans le divertor
	    % ref Stangeby section 6.7 
	    % ref originale : D. Reiter et al, PPCF vol 33 (1991) p 1579-1600
	    tau_ref    = data_zerod.tauhe - z0dstruct.z0dinput.option.Recycling ./ (1 - z0dstruct.z0dinput.option.Recycling) .* data_zerod.taup;
    else
	    tau_ref    = data_zerod.tauhe;
    end
    % gaz fuelling effciency
    mat.eta_g_low  =  data_zerod.taup ./ tau_ref;
    % source correspondant a la difference entre la densit� et la densit� naturelle du plasma
    mat.dfuelling_bilan_dt  = fact_nbar .* (data_zerod.nbar - data_zerod.nbar_nat) ./ tau_ref .* data_zerod.vp ./ (4.41e-4 .* phys.avo) - ...
			      (mat.snbi_star + mat.pellet_star);
    if (filter_width > 1) && (length(data_zerod.temps) > filter_width)
        try
	    mat.dfuelling_bilan_dt     = medfilt1(mat.dfuelling_bilan_dt,filter_width);
	end
    end
    if length(data_zerod.temps) > 1
	    mat.fuelling_bilan         = cumtrapz(data_zerod.temps,mat.dfuelling_bilan_dt);
    else
	     mat.fuelling_bilan        = 0;  
    end
    % source de gaz puff dans le plasma si efficacite de 1
    mat.gas_puff          =  mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt > 0);
    mat.wall_pumping      = -mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt < 0);
    if ~isfinite(fn0a)
      fn0a = z0dstruct.z0dinput.option.fn0a;
    end
    if ~isfinite(fn0a_div)    
      fn0a_div = z0dstruct.z0dinput.option.fn0a_div;
    end
    switch z0dstruct.z0dinput.option.configuration
    case {0,1}
	mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a); 
    case {2,3}
	mat.gas_puff      = mat.gas_puff ./ max(eps,data_zerod.xpoint .* fn0a_div + (~data_zerod.xpoint) .* fn0a);
    otherwise
	mat.gas_puff      = mat.gas_puff ./ max(eps,fn0a_div); 
    end

    % source totale dans le plasma
    mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
    % bilan = input - losses = 0 + changement de densite du plasma
    mat.bilan = mat.input - mat.plasmanetflux;
    % correction de gaz puff
    mat.gas_puff = mat.gas_puff - 1.1 .* mat.bilan .* (mat.bilan < 0);
    %mise a jour  bilan = input - losses = 0 + changement de densite du plasma
    mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
    mat.bilan = mat.input - mat.plasmanetflux;
    % pumping
    mat.pumping = max(mat.bilan,mat.wall_pumping);
    if length(data_zerod.temps) > 1
	    mat.bilan_pumping = cumtrapz(data_zerod.temps,mat.pumping);
    else
	    mat.bilan_pumping =  0;  
    end
    % correction to gas puff
    puff_correc  = (mat.pumping - mat.bilan);
    mat.gas_puff = mat.gas_puff + puff_correc .* (puff_correc> 0);
    if length(data_zerod.temps) > 1
	    mat.bilan_gas_puff = cumtrapz(data_zerod.temps,mat.gas_puff);
    else
	    mat.bilan_gas_puff = 0;  
    end
    mat.input         = mat.pellet + mat.snbi + mat.gas_puff;
    % bilan final
    mat.bilan = mat.input - mat.plasmanetflux  -mat.pumping;
    if length(data_zerod.temps) > 1
	    mat.bilan_input   = cumtrapz(data_zerod.temps,mat.input);
    else
	    mat.bilan_input   = 0;  
    end


    % flux to divertor/limiter
    mat.flux_divlim_total = mat.output + mat.gas_puff  + (1 - eta_p) .* mat.pellet + mat.perte_snbi - mat.s0_in;
    mat.recycle           = mat.flux_divlim_total - mat.pumping;
    mat.Reff              = mat.recycle ./ mat.flux_divlim_total;
    if length(data_zerod.temps) > 1
	    mat.bilan_recycle     = cumtrapz(data_zerod.temps,mat.recycle);
    else
	    mat.bilan_recycle     = 0;  
    end
    mat.eta_g_high        = mat.s0_in ./ max(eps,mat.output + mat.gas_puff + mat.recycle);




otherwise
    % flux de matiere sortant
    mat.temps        = z0dstruct.zerod.temps;
    mat.output       = interp1(z0dstruct.profil.temps,z0dstruct.profil.ge(:,end) .* z0dstruct.profil.grho2(:,end) .* z0dstruct.profil.vpr_tor(:,end), ...
			      z0dstruct.zerod.temps,'pchip','extrap') ./ (4.41e-4 .* phys.avo);
    mat.recycle      = z0dstruct.zerod.n0a ./ (4.41e-4 .* phys.avo);  % recyclage  + gas puff
    mat.pellet       = interp1(z0dstruct.profil.temps,trapz(z0dstruct.profil.xli,z0dstruct.profil.spellet .* z0dstruct.profil.vpr,2), ...
			      z0dstruct.zerod.temps,'pchip','extrap')  ./ (4.41e-4 .* phys.avo);
    mat.pellet(~isfinite(mat.pellet)) = 0;
    mat.snbi         = z0dstruct.zerod.pnbi  ./ z0dstruct.z0dinput.option.einj ./ phys.e ./ (4.41e-4 .* phys.avo);
    mat.input        = mat.recycle + mat.pellet + mat.snbi;

    mat.bilan_input   = cumtrapz(z0dstruct.zerod.temps,mat.input);
    mat.bilan_recycle = cumtrapz(z0dstruct.zerod.temps,mat.recycle);
    mat.bilan_pellet  = max(eps,cumtrapz(z0dstruct.zerod.temps,mat.pellet));
    mat.bilan_nbi     = max(eps,cumtrapz(z0dstruct.zerod.temps,mat.snbi));
    mat.bilan_output  = cumtrapz(z0dstruct.zerod.temps,mat.output);     

    ntot                 = trapz(z0dstruct.profil.xli,z0dstruct.profil.vpr .* z0dstruct.profil.nep,2) ./ (4.41e-4 .* phys.avo);
    mat.contents         = interp1(z0dstruct.profil.temps,ntot,z0dstruct.zerod.temps,'pchip','extrap');
    mat.bilan_io         = mat.bilan_input - mat.bilan_output;
    % la precision numerique du calcul est limite (la source est mieux connu que le flux sortant)
    mat.bilan_error      = mat.contents - mat.bilan_io -mat.contents(1); 

    %mat.fuelling_pumping = z0dxdt(mat.bilan_io,z0dstruct.zerod.temps);
    mat.plasmanetflux    = z0dxdt(mat.contents,z0dstruct.zerod.temps); % flux net
    mat.flux_error       = z0dxdt(mat.bilan_error,z0dstruct.zerod.temps); % flux convecte et/ou diffuse
    % separation
    % correction flux sortant 
    mat.output        = mat.output - mat.flux_error ;
    mat.bilan_output  = mat.bilan_output - mat.bilan_error;     
    mat.bilan_io      = mat.bilan_input - mat.bilan_output;


    % calcul gaz puff (a partir du scaling pour la densite naturelle)
    % scaling
    fact_nbar = z0dstruct.zerod.nem ./ max(1,z0dstruct.zerod.nbar);   
    tauref = min(1e3,max(min(z0dstruct.zerod.tauhe,z0dstruct.zerod.taup) ./ max(eps,1 - z0dstruct.z0dinput.option.Recycling),1e-6));
    snem   = fact_nbar .* z0dstruct.zerod.nbar_nat./ tauref .* z0dstruct.zerod.vp + z0dstruct.zerod.pnbi_th ./ z0dstruct.z0dinput.option.einj;
    [tntot,mat.ntot_nat] = z0ode(z0dstruct.zerod.temps,snem,tauref,fact_nbar(1) .* z0dstruct.zerod.nbar_nat(1) .* z0dstruct.zerod.vp(1));
    %else
    %  mat.ntot_nat = fact_nbar .* z0dstruct.zerod.nbar_nat .* z0dstruct.zerod.vp;
    %end
    % facteur de conversion
    mat.fuelling_bilan         = z0dstruct.zerod.vp .* z0dstruct.zerod.nem  - mat.ntot_nat - (mat.bilan_pellet + mat.bilan_nbi) .* (4.41e-4 .* phys.avo);
    mat.dfuelling_bilan_dt     = z0dxdt(mat.fuelling_bilan,z0dstruct.zerod.temps);
    mat.gas_puff               = mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt> 0) ./ (4.41e-4 .* phys.avo);
    mat.pumping                = -mat.dfuelling_bilan_dt .* (mat.dfuelling_bilan_dt < 0) ./ (4.41e-4 .* phys.avo);
    mat.recycle                = max(0,mat.recycle - mat.gas_puff);
    mat.coef_recycle           = mat.recycle  ./ max(mat.input,mat.output); % securite a cause du bruit sur les donnees + convection
    mat.bilan_recycle          = cumtrapz(z0dstruct.zerod.temps,mat.recycle);
    mat.bilan_gas_puff         = cumtrapz(z0dstruct.zerod.temps,mat.gas_puff);

    % extraction du temps d'interet
    if length(data_zerod.temps) == 1
	noms = fieldnames(mat);
	for k=1:length(noms)
	    mat.(noms{k}) = interp1(z0dstruct.zerod.temps,mat.(noms{k}),data_zerod.temps,'nearest','extrap');
	end
    end
end


function metis4itm_save_restart(filename,post,diary_text)

% variable global importante pour les test
% utilisation du mexfile si disponible	 
if isappdata(0,'MEXSOLVER_IN_METIS')
	mexsolver = getappdata(0,'MEXSOLVER_IN_METIS');
else	
	mexsolver =  [];
end
if isempty(mexsolver)
	repmex=which(strcat('mexpde1dsolver.',mexext));
	if ~isempty(repmex)
		mexsolver = 1;
	else
		mexsolver = 0;	
	end
	setappdata(0,'MEXSOLVER_IN_METIS',mexsolver);
end  

langue      =  lower(getappdata(0,'langue_cronos'));
if isempty(post)
	warning('No data to be saved');
	return
end
if ~isfield(post,'zerod')
	warning('No data to be saved');
	return
end
if ~isfield(post,'z0dinput')
	warning('No data to be saved');
	return
end
post.profil0d = post.profil;
data.zerod    = post.zerod;
data.z0dinput = post.z0dinput;
data.profil0d = post.profil0d;
% sauvegarde des donnees externes
noms = fieldnames(getappdata(0));
for k = 1:length(noms)
	if findstr(noms{k},'_EXP')
		data.appdata.(noms{k}) = getappdata(0,noms{k});
	end
end

% sauvegarde de la sortie standard du model simulink
try
	simout = evalin('base','simout');
catch
	simout = [];
end
zassignin('base','post.simout',simout);
data.simout = simout;

% donnees LUKE  dans METIS
if isfield(post,'lukeinmetis')
	data.lukeinmetis = post.lukeinmetis;
end   


% ajout des donnees de certication
root = getappdata(0,'root');
post = data;
clear data
save(filename,'post','diary_text');

% appel de la fonction qui regle le bon cocos
function [data_zerod,profil0d,sigma_B0_eff,sigma_bvac_r,factor_two_pi] = makecocos(data_zerod,profil0d,z0dstruct,COCOS_out,factor_two_pi)

% par defaut
sigma_B0_eff = 1;
%  Btor		Ip		Psi		Phi		safety factor
%  positive 	positive 	decreasing 	increasing 	negative
%  positive 	negative 	increasing 	increasing 	positive
%  negative 	positive 	decreasing 	decreasing 	positive
%  negative 	negative 	increasing 	decreasing 	negative
% consitency tests
% by choice in METIS ip > 0 and sign of Btor is given by option.signe
s_psi = sign(sum(sign(profil0d.psi(:,1) - profil0d.psi(:,end))));
s_phi = sign(sum(sign(profil0d.phi(:,1) - profil0d.phi(:,end))));
s_q   = sign(sum(sign(profil0d.qjli(:))));
switch z0dstruct.z0dinput.option.signe
case -1
	if s_psi < 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi < 0
	  	error('Unconsistent Phi variation : Phi must be decreasing from magnetic axis to LCFS');
	end
	if s_q < 0
	  	error('Unconsistent safety factor sign :  safety factor must be positive');
	end
otherwise
	if s_psi < 0
	  	error('Unconsistent Psi variation : Psi must be decreasing from magnetic axis to LCFS');
	end
	if s_phi > 0
	  	error('Unconsistent Phi variation : Phi must be increasing from magnetic axis to LCFS');
	end
	if s_q > 0
	  	error('Unconsistent safety factor sign :  safety factor must be negative');
	end
end
% CoCos de METIS
sigma_Ip_in = 1; % par definition dans metis
sigma_B0_in = z0dstruct.z0dinput.option.signe;
sigma_Bp_in = sigma_Ip_in .* sign(mean(sign(profil0d.psi(:,end) - profil0d.psi(:,1))));
s_dpdspi    = sign(sum(sign((profil0d.ptot(:,end) - profil0d.ptot(:,1)) ./ (profil0d.psi(:,end) - profil0d.psi(:,1)))));
if s_dpdspi ~= (-sigma_Bp_in * sigma_Ip_in)
	error('consisency error :sign(dP/dPsi) is not - sigma_Bp_in * sigma_Ip_in');
end
sigma_rhothetaphi_in = sigma_Ip_in .* sigma_B0_in .* s_q;
% we assume coordinate (R,phi,Z) right-handed
if  sigma_Bp_in > 0
    if sigma_rhothetaphi_in > 0
        COCOS_in = 1;
    else
        COCOS_in = 5;
    end
else
    if sigma_rhothetaphi_in > 0
        COCOS_in = 7;
    else
        COCOS_in = 3;
    end
end
COCOS_in = COCOS_in + 10;

% COCOS METIS
[Kexp_Bp_in,Ksigma_Bp_in,Ksigma_RphiZ_in,Ksigma_rhothetaphi_in,Ksign_q_pos_in,Ksign_pprime_pos_in] = cocos(COCOS_in);
% COCOS cible
[Kexp_Bp_out,Ksigma_Bp_out,Ksigma_RphiZ_out,Ksigma_rhothetaphi_out,Ksign_q_pos_out,Ksign_pprime_pos_out] = cocos(COCOS_out);

% verifications
if any(sign(profil0d.qjli(:).* Ksigma_rhothetaphi_in .* sigma_Ip_in .* sigma_B0_in)<= 0)
    fprintf('WARNING: sign(q) is not consistent with COCOS_in= %d\n',COCOS_in)
    fprintf('qedge = %g\n',mean(profil0d.qjli(:,end)))
    fprintf('sig_rhothetaphi*sign(Ip)*sign(B0) = %d * %d *%d = %d\n', ...
            Ksigma_rhothetaphi_in, sigma_Ip_in,sigma_B0_in, ...
            Ksigma_rhothetaphi_in*sigma_Ip_in*sigma_B0_in);
end 
if any(sign(profil0d.fdia(:) .* sigma_B0_in)<= 0)
    fprintf('WARNING: Signs of F and B0 are not consistent\n');
end
if any(sign((profil0d.psi(:,end) - profil0d.psi(:,1)).*Ksigma_Bp_in.*sigma_Ip_in) <= 0)
    if  any(sign(profil0d.psi(:,end) - profil0d.psi(:,1)) <= 0)
      fprintf('WARNING: psi should be decreasing with : sign(Ip)= %d and %d  for COCOS= %d\n',sigma_Ip_in,sigma_Ip_in,COCOS_in);
    else
      fprintf('WARNING: psi should be increasing with : sign(Ip)= %d and %d  for COCOS=%d\n',sigma_Ip_in,sigma_Ip_in,COCOS_in);
    end
end

% for cpo scenario
sigma_bvac_r = sigma_B0_in;

if COCOS_in == COCOS_out 
	return
end

%  
%  Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
%  sign(Ip) in output:
%     
KIPsign_out = sigma_Ip_in;
KB0sign_out = sigma_B0_in;
sigma_RphiZ_eff  = Ksigma_RphiZ_out * Ksigma_RphiZ_in;
%sigma_IP_eff = sigma_RphiZ_eff
sigma_Ip_eff = sigma_Ip_in * KIPsign_out;
sigma_Ip_out = sigma_Ip_in * sigma_Ip_eff;
%sigma_B0_eff = Ksigma_RphiZ_in * Ksigma_RphiZ_out
sigma_B0_eff = sigma_B0_in * KB0sign_out;
sigma_B0_out = sigma_B0_in * sigma_B0_eff;
sigma_Bp_eff = Ksigma_Bp_out * Ksigma_Bp_in;
exp_Bp_eff = Kexp_Bp_out - Kexp_Bp_in;
sigma_rhothetaphi_eff  = Ksigma_rhothetaphi_out * Ksigma_rhothetaphi_in;
fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2.*pi) ^ exp_Bp_eff;
fact_q = sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff;
factor_two_pi = factor_two_pi .* abs(fact_psi);
% transformation
% changement normalisation psi
profil0d.psi            = profil0d.psi    .* fact_psi;
profil0d.dpsidt         = profil0d.dpsidt .* fact_psi;
profil0d.qjli           = profil0d.qjli   .* fact_q;  
profil0d.fdia           = profil0d.fdia .* sigma_B0_eff;
profil0d.phi            = profil0d.phi .* sigma_B0_eff;
profil0d.dphidx         = profil0d.dphidx .* sigma_B0_eff;
profil0d.jeff           = profil0d.jeff .* sigma_Ip_eff;
profil0d.jboot          = profil0d.jboot .* sigma_Ip_eff;
profil0d.jnbicd         = profil0d.jnbicd .* sigma_Ip_eff;
profil0d.jlh            = profil0d.jlh   .* sigma_Ip_eff;
profil0d.jeccd          = profil0d.jeccd .* sigma_Ip_eff;
profil0d.jfwcd          = profil0d.jfwcd .* sigma_Ip_eff;
profil0d.jfus           = profil0d.jfus .* sigma_Ip_eff;
profil0d.jrun           = profil0d.jrun  .* sigma_Ip_eff; 
profil0d.jni            = profil0d.jni .* sigma_Ip_eff; 
profil0d.jfusshape      = profil0d.jfusshape .* sigma_Ip_eff; 
profil0d.epar           = profil0d.epar .* sigma_Ip_eff;
profil0d.omega          = profil0d.omega *  sigma_RphiZ_eff;
profil0d.utheta         = profil0d.utheta * sigma_rhothetaphi_eff  *  sigma_RphiZ_eff ./ fact_psi;
profil0d.vtheta         = profil0d.vtheta * sigma_rhothetaphi_eff  *  sigma_RphiZ_eff;
profil0d.vtor           = profil0d.vtor    *  sigma_RphiZ_eff;
profil0d.rot_nbi        = profil0d.rot_nbi *  sigma_RphiZ_eff;
profil0d.rot_n0         = profil0d.rot_n0  *  sigma_RphiZ_eff;
profil0d.rot_lh         = profil0d.rot_lh  *  sigma_RphiZ_eff;

data_zerod.q0           = data_zerod.q0   .* fact_q;  
data_zerod.q95          = data_zerod.q95  .* fact_q; 
data_zerod.qa           = data_zerod.qa   .* fact_q; 
data_zerod.qmin         = data_zerod.qmin .* fact_q; 
data_zerod.phiplasma    = data_zerod.phiplasma .* sigma_B0_eff;
data_zerod.ifus         = data_zerod.ifus  .* sigma_Ip_eff;
data_zerod.jxfus        = data_zerod.jxfus .* sigma_Ip_eff;
data_zerod.j0fus        = data_zerod.j0fus .* sigma_Ip_eff;
data_zerod.iohm         = data_zerod.iohm  .* sigma_Ip_eff;
data_zerod.vloop        = data_zerod.vloop .* sigma_Ip_eff;
data_zerod.ip           = data_zerod.ip    .* sigma_Ip_eff;
data_zerod.ifwcd        = data_zerod.ifwcd .* sigma_Ip_eff;
data_zerod.ieccd        = data_zerod.ieccd .* sigma_Ip_eff;
data_zerod.inbicd       = data_zerod.inbicd.* sigma_Ip_eff;
data_zerod.iboot        = data_zerod.iboot .* sigma_Ip_eff;
data_zerod.wrad         = data_zerod.wrad  .* sigma_RphiZ_eff;
data_zerod.irun         = data_zerod.irun  .* sigma_Ip_eff;
data_zerod.sn0fr        = data_zerod.sn0fr .* sigma_RphiZ_eff;
data_zerod.ilh          = data_zerod.ilh   .* sigma_Ip_eff;
data_zerod.icd          = data_zerod.icd   .* sigma_Ip_eff;
data_zerod.qeff         = data_zerod.qeff  .* fact_q; 
data_zerod.vmes         = data_zerod.vmes  .* sigma_Ip_eff;
data_zerod.ipar         = data_zerod.ipar  .* sigma_Ip_eff;
data_zerod.ini          = data_zerod.ini   .* sigma_Ip_eff;
data_zerod.edgeflux     = data_zerod.edgeflux .* sigma_Ip_eff * sigma_Bp_eff;

% for Bphi . grad(phi)
sigma_bvac_r = sigma_B0_eff .* sigma_RphiZ_eff .* sigma_bvac_r;


% fonction pour la gestion du cache de l'UAL
function [idx,state] = mem_cache(idx)

if isempty(import) && ~isappdata(0,'JAVAINTERFACESET') && isempty(getappdata(0,'JAVAINTERFACESET'))
	import ualmemory.javainterface.*;
end
setappdata(0,'JAVAINTERFACESET','done')

if nargin  == 0
	idx = [];
	state = 0;
        setappdata(0,'euitm_enable_mem_cache',0);
	return
elseif isempty(idx)
	state = 0;
        setappdata(0,'euitm_enable_mem_cache',0);
	return
end

% en attendant la mise en place d'une synchronisation avec Kepler, le mecanisme n'est pas utilisé.
setappdata(0,'euitm_enable_mem_cache',0);
state = 0;
return

try
	UALAccess.enableMemCache(idx)
	setappdata(0,'euitm_enable_mem_cache',1);
	state = 1;
catch
	fprintf('Unable to turn on UAL memory cache:\n%s\n',lasterr);
        setappdata(0,'euitm_enable_mem_cache',0);
	state = 0;
end


function nfast  = compute_nfast(data_zerod,profil0d,option,cons,model)

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

% 0 template
nfast = zeros(size(model));

% resample data zerod
temps = data_zerod.temps;
if length(temps ) > 1
    noms = fieldnames(data_zerod);
    for k=1:length(noms)
        var = data_zerod.(noms{k});
	if length(temps) == length(var)
		data_zerod.(noms{k}) = interp1(temps,var,profil0d.temps,'linear','extrap');
        end 
    end
    data_zerod.temps = profil0d.temps;
end
% resample cons
temps = cons.temps;
if length(temps ) > 1
    noms = fieldnames(cons);
    for k=1:length(noms)
        var = cons.(noms{k});
	if length(temps) == length(var)
		cons.(noms{k}) = interp1(temps,var,profil0d.temps,'linear','extrap');
        end 
    end
    cons.temps = profil0d.temps;
end

% useful vector
ve  = ones(size(profil0d.xli));
vt  = ones(size(profil0d.temps));

% ATTENTION : normalisation des mass sur le gaz principal 
% density of fast ions
psup_alpha  = profil0d.pfus;
psup_alpha  = psup_alpha ./ (max(1, trapz(profil0d.xli,psup_alpha .* profil0d.vpr,2)) * ve) .*  (data_zerod.esup_fus * ve)  ./ (3/2);
nfast_alpha = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	      data_zerod.meff * ve,4,2,3.56e6,psup_alpha);
talpha      = psup_alpha ./ max(1,nfast_alpha) ./  phys.e;

switch option.gaz
case 3
	  minj = 2 .* (1-cons.ftnbi ) + 3 .* cons.ftnbi ;
otherwise
	  minj = 2 .* (1-cons.ftnbi ) + 1 .* cons.ftnbi ;
end
psup_nbi1  = real(profil0d.nbinesource );
psup_nbi1  = psup_nbi1 ./ (max(1, trapz(profil0d.xli,psup_nbi1 .* profil0d.vpr ,2)) * ve) .*  (real(data_zerod.esup_nbi )* ve) ./ (3/2);
nfast_nbi1 = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	      data_zerod.meff * ve,minj*ve,1,option.einj,psup_nbi1);

psup_nbi2  = imag(profil0d.nbinesource );
psup_nbi2  = psup_nbi2 ./ (max(1, trapz(profil0d.xli,psup_nbi2 .* profil0d.vpr ,2)) * ve) .*  (imag(data_zerod.esup_nbi ) * ve) ./ (3/2);
nfast_nbi2 = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	     data_zerod.meff * ve ,minj*ve,1,option.einj2,psup_nbi2);

tnbi      = (psup_nbi1 + psup_nbi2) ./ max(1,nfast_nbi1 +nfast_nbi2) ./  phys.e;

% choix du minoritaire
switch option.mino
case 'He3'
   ag = 3;
   zg = 2;
   lg = 7.92e-3;
case 'T'
   ag = 3;
   zg = 1;
   lg = 7.92e-3;
case 'He4'
   ag = 4;   
   zg = 2;
   lg = 4.55e-3;
case 'D'
   ag = 2;   
   zg = 1;
   lg = 6.46e-3;
otherwise
   ag = 1;   
   zg = 1;
   lg = 4.576e-3;
end
psup_icrh   = profil0d.picrh ;
psup_icrh   = psup_icrh ./ (max(1, trapz(profil0d.xli,psup_icrh .* profil0d.vpr ,2)) * ve) .*  (data_zerod.esup_icrh  *ve) ./ (3/2);
nfast_icrh = zerod_fast_ions_density(profil0d.nep ,profil0d.tep ,profil0d.zeff , ...
	     data_zerod.meff * ve,ag,zg,data_zerod.einj_icrh*ve ,psup_icrh);
ticrh       = psup_icrh ./ max(1,nfast_icrh) ./  phys.e;
warning off
% pas de meilleurs formulation pour le moment
nfast_lh  = max(0,data_zerod.esup_lh ) ./ (data_zerod.einj_lh  .* phys.e) .* 2 ./ data_zerod.vp ;
warning on
psup_lh   = profil0d.plh ;
nfast_lh  = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr ,2)) * ve) .*  (nfast_lh * ve);
psup_lh   = psup_lh ./ (max(1, trapz(profil0d.xli,psup_lh .* profil0d.vpr ,2)) * ve) .*  (data_zerod.esup_lh * ve)  ./ (3/2);
tlh       = psup_lh ./ max(1,nfast_lh) ./  phys.e;


% hydrogenoid density
nDp   = zeros(size(profil0d.n1p));
nTp   = zeros(size(profil0d.n1p));
nHp   = zeros(size(profil0d.n1p));

% helium
nHep4  = zeros(size(profil0d.nhep));
nHep3  = zeros(size(profil0d.nhep));

% compute corrected density depending on minority scheme
switch option.mino
case 'H'
  nHp = nfast_icrh;  
case 'T'
  nTp = nfast_icrh;
case 'He3'
  nHep3 = nfast_icrh;  
case 'He4'
  nHep4 = nfast_icrh;
otherwise
  error(sprintf('minority species %s not yet implemanted',option.mino));
end

% correction of fast NBI ions
switch option.gaz
case 3
      nTp = nTp + real(cons.ftnbi * ve)  .* nfast_nbi1;
      nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi1;
      nTp = nTp + imag(cons.ftnbi * ve)  .* nfast_nbi2;
      nDp = nDp + (1 - imag(cons.ftnbi * ve)) .* nfast_nbi2;
otherwise
      nHp = nHp + real(cons.ftnbi * ve)  .* nfast_nbi1;
      nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi1;
      nHp = nHp + real(cons.ftnbi * ve)  .* nfast_nbi2;
      nDp = nDp + (1 - real(cons.ftnbi * ve)) .* nfast_nbi2;
end

% contribution from fast alpha
nHep4 = nHep4 + nfast_alpha;


% fill template 
nfast(:,:,1) = nHp;
nfast(:,:,2) = nDp;
nfast(:,:,3) = nTp;
nfast(:,:,4) = nHep3;
nfast(:,:,5) = nHep4;
