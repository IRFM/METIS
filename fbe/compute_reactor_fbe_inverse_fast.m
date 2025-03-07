function option_feeqs_reactor = compute_reactor_fbe_inverse_fast(option_feeqs_reactor)



% GUI generation
if nargin < 1
    
    valeur.coils_description	      = '';
    type.coils_description          = 'string';
    borne.coils_description         = '';
    defaut.coils_description        = '';
    info.coils_description          = 'path and name of coils description file used for the computation (if empty, a file selector will be open at the next step);\nsee example in directory of project Reactor in FEEQS sources for the format of this file.';
    
    valeur.metisfilename	      = '';
    type.metisfilename          = 'string';
    borne.metisfilename         = '';
    defaut.metisfilename        = '';
    info.metisfilename          = 'path and name of METIS simulation file used for the computation (if empty, a file selector will be open at the next step)';
    
    valeur.first_wall_name	      = '';
    type.first_wall_name          = 'string';
    borne.first_wall_name         = '';
    defaut.first_wall_name        = '';
    info.first_wall_name          = 'path and name of first wall description file used for the computation (can be left empty, optional parameter);\nsee example in directory of project Reactor in FEEQS sources for the format of this file.';
    
    
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
    
    
    valeur.flux_constant_mode  = 'manual';
    type.flux_constant_mode    = 'string';
    borne.flux_constant_mode   = {'manual','start of flux','balanced','end of flux'};
    defaut.flux_constant_mode  = 'manual';
    info.flux_constant_mode   = ['This parameter allows to choose the way poloidal flux offset used to match METIS current diffusion LCFS poloidal flux and FEEQS.M LCFS poloidal flux is computed:\n', ...
        'if = manual,offset is prescribed by the user\n', ...
        'in other offset is automatically computed taking into account flux consumption estimation for breakdown computed in METIS and flux leakage computed with the help of current diffusion equation between end of breakdown and  first_time compare to flux computed by premagnetisation tool:', ...
        '    case start_of_flux: offset is computed to have maximum available flux for the first plasma, taking into account breakdonw flux\n', ...
        '    case balanced: offset is comuted to have flux margin for ramp-up and ramp down if possible\n', ...
        '    case end_of_flux:offset is computed to have minimum flux value (negative) at the end of rampdwon.' ...
        ];
    
    valeur.flux_offset    = 0;
    type.flux_offset      = 'float';
    borne.flux_offset     = [-1000,1000];
    defaut.flux_offset    = 0;
    info.flux_offset      = 'poloidal flux offset used to match METIS current diffusion LCFS poloidal flux and FEEQS.M LCFS poloidal flux (Wb)';
    
    valeur.breakdown_flux   = 0.3;
    type.breakdown_flux     = 'float';
    borne.breakdown_flux    = [0, 20];
    defaut.breakdown_flux   = 0.3;
    info.breakdown_flux   = 'poloidal flux consumed for plasma break-down and burn-through (Wb);\nif METIS breakdown model is turn on, take the maximum of this value and METIS value.';
    
    valeur.weight    = 1;
    type.weight      = 'float';
    borne.weight     = [1e-16,1e6];
    defaut.weight    = 1;
    info.weight      = 'regularisation term weight in cost function (sum(I^2) and dipole)';
    
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
    
    option_feeqs_reactor.valeur     = valeur;
    option_feeqs_reactor.type       = type;
    option_feeqs_reactor.borne      = borne;
    option_feeqs_reactor.defaut     = defaut;
    option_feeqs_reactor.info       = info;
    option_feeqs_reactor.interface  = interface;
    
    option_feeqs_reactor.description = 'Computation of coil currents for reactors using FEEQS (Fast coil currents identification)';   % description (une ligne) de la fonction
    
    option_feeqs_reactor.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    option_feeqs_reactor.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    option_feeqs_reactor.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
    
    % end of GUI form declaration
    return
    
end
% select file is needed
[pathstr,name,ext] = fileparts(option_feeqs_reactor.metisfilename);
if isempty(ext)
    option_feeqs_reactor.metisfilename = sprintf('%s.mat',option_feeqs_reactor.metisfilename);
end
if isempty(option_feeqs_reactor.metisfilename) || ~exist(option_feeqs_reactor.metisfilename)
    [FileName,PathName] = uigetfile('*.mat','Select a METIS file');
    if isempty(FileName) || isnumeric(FileName)
        disp('Call to FEEQS.M has been canceled');
        return
    end
    option_feeqs_reactor.metisfilename = fullfile(PathName,FileName);
end
%
% select file is needed
[pathstr,name,ext] = fileparts(option_feeqs_reactor.coils_description);
if isempty(ext)
    option_feeqs_reactor.coils_description = sprintf('%s.mat',option_feeqs_reactor.coils_description);
end
if isempty(option_feeqs_reactor.coils_description) || ~exist(option_feeqs_reactor.coils_description)
    [FileName,PathName] = uigetfile('*.txt','Select a file for coils desctiption');
    if isempty(FileName) || isnumeric(FileName)
        disp('Call to FEEQS.M has been canceled');
        return
    end
    option_feeqs_reactor.coils_description = fullfile(PathName,FileName);
end

disp(' ');
% call FEESQ.M
dirmem = pwd;
fprintf('Initial directory is %s\n',dirmem);
switch option_feeqs_reactor.flux_constant_mode
    case 'manual'
        offset_flux = option_feeqs_reactor.flux_offset / 2 / pi;
    case 'auto'
        offset_flux = NaN;
    otherwise
        offset_flux = option_feeqs_reactor.flux_constant_mode;
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
try
    output =  reactor_inverse4metis_fast(option_feeqs_reactor.coils_description,option_feeqs_reactor.metisfilename, ...
        option_feeqs_reactor.number_of_time_slices,option_feeqs_reactor.first_time,option_feeqs_reactor.last_time, ...
        option_feeqs_reactor.reinforce,option_feeqs_reactor.plotonoff, option_feeqs_reactor.movie_onoff,option_feeqs_reactor.run_name, ...
        option_feeqs_reactor.weight,offset_flux,dirmem,option_feeqs_reactor.breakdown_flux, ...
        option_feeqs_reactor.mode_profile,option_feeqs_reactor.V_Force,option_feeqs_reactor.first_wall_name);
    
    cd(dirmem);
catch
    keyboard
    cd(dirmem);
    error(lasterror);
end
option_mem = option_feeqs_reactor;
option_feeqs_reactor = output;
option_feeqs_reactor.flux_offset = 2 .* pi .* output.psi_off;
noms = fieldnames(option_mem);
for k=1:length(noms)
    if ~isfield(option_feeqs_reactor,noms{k})
        option_feeqs_reactor.(noms{k}) = option_mem.(noms{k});
    end
end
