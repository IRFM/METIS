% METIS4MUSCLE3: this fonction is the bridge between MUSCLE3 and METIS
%
%  This function connects METIS to a MUSCLE3 workflow
%  (https://muscle3.readthedocs.io/en/latest/). Data are exchanged in IDS
%  format (https://github.com/iterorganization). Various comptutation modes
%  are implemented in this actor. 
%
% % syntax:
%
%    metis4muscle3(debug_level,checkpoint_on_off,port_name_legacy)
%
% inputs:
%
%   debug_level = level of debug and Matlab profiler state. If degug_level
%                 is negative, the Matlab profiler is switch on and
%                 a repport is created at the end of the simulation in the 
%                 workdir of the METIS instance. 
%                 The possible values of  abs(debug_level) are:
%                      0 = no debug mecanism
%                      1 = Matlab logfile is written in workdir of the 
%                          METIS instance (debug_metis_in_muscle3.txt)
%                      2 = as 1 + debuger stop on error
%                      3 = as 2 + Matlab stop at the start of METIS4MUSCLE3
%                          (keyboard insertion, work only in interactive
%                          mode of Matlab).
%                  Additionnaly, if debug_level is non zero, internal
%                  variables of the function metis4muscle3 are saves in the
%                  workdir of the METIS instance 
%                  (meti4muscle3_caller_*_context.mat) and the Matlab 
%                   workspace is saved in the workdir of the 
%                  METIS instance (meti4muscle3_workspace_context.mat). 
%
%   checkpoint_on_off = if true, turns on the checkpointing features of the
%                       actor and allows creation of checkpoint dataset
%                       and cold and hot restart of the worflow from a  
%                       checkpoint dataset.
%
%   port_name_legacy = if true, use original port name before
%                      standardization by ITER (default = false)
%
% remark:
%  The actor metis4muscle3 can be launched in a interactive Matlab (option
%  -r of matlab command) or in batch mode (option -batch of matlab
%  command). The interactive mode is recommended for debugging and the
%  batch mode is recommended for standard execution as the start of Matlab
%  is faster.
%
%  Following execution modes are available:
%
%     * resuming = restart of a previously stopped simulation.
%
%     * replay   = replay a METIS simulation file and broacast IDSs
%                  messages preserving time_stamp from original simulation.
%
%     * interpolate = Broadcast IDSs messages from a METIS simulation file
%                     using time_stamp provided by MUSCLE3 workflow.
% 
%     * working_point = compute the working point solution (the stationnary
%                       time independant spolution). The purpose of the
%                       mode is the coupling to system code.
% 
%     * metis_init = initialisation of the METIS actor. There two way to do
%                    it: cold initialisation where intitiale data are read
%                    from a METIS simulation file or hot initialisation
%                    where the intiale data are computed by running METIS.
%
%     * metis_computation = METIS will compute new data at the time
%                           prescribed by the actor using previously know
%                           internal state. This time can be in future for
%                           compute evolution of plasma state, be the same
%                           than the previous one to perform some
%                           convergence in the workflow or in past to
%                           restart from previous converged state of the
%                           workflow in case of non convergence in some
%                           actor.
%
% The datation in this actor can be provided either by:
%    1/ another actor on pulse_schedule port, equilibrium port, *_profiles
%       port, *_sources port or clock port using time_stamp or next_time_stamp
%    2/ METIS simulation file
%    3/ interal actor clock.
%
% The actor has 3 sets of input ports:
%   1/ one set for initialisation (F_INIT)
%   2/ two sets for time evolution (F_INIT and S with precedance of S ports).
%
% Input ports are:
%   * clock -> use only time_stamp or next_time_stamp information from
%              message, can be any IDS.
%
%   * pulse_schedule -> should contain main input data for METIS
%
%   * equilibrium -> data for extenal equilibrium and/or external current
%                    diffusion.
% 
%   * plasma_profiles or core_profiles -> kinetic data for external data
%
%   * plasma_source or core_sources -> data for external sources.
% 
% The actor has 2 set of output ports that can be both use (O_F and O_I)
% Output port IDSs are:
%
%   * pulse_schedule
%   * equilibrium
%   * plasma_profiles and core_profiles
%   * plasma_sources and core_sources
%   * dataset_description (for DD version before 4.x)
%   * dataset_fair
%   * summary
%   * plasma_transport and core_transport
%   * radiation
%   * edge (2 IDSs: edge_profiles and edge_transport, produced only if 
%           two point model is switch on in METIS)
%   * transport_solver_mumerics.
%
% Corresponding IDS are mapped only if the port is connected.
%
% Examples of configuration files are available in the directory:
%      <METIS_ROOT_DIR>/worflow/muscle3/YMMSL
%
% These examples provides some information on the setting of the actor.
%
% Examples of scripts to run a worklow are available in the directory:
%    <METIS_ROOT_DIR>/worflow/muscle3/run_script
%
% The configuration of the environment, set in this examples, should be 
% updated depending of the local installation of Matlab, IMAS and MUSCLE3.
% 
% More comprehensive informations are available here:
%    < insert ref here>.
%
% J-F. Artaud 2024, 2025, CEA
% R. Nouailletas 2024, CEA

function metis4muscle3(debug_level,checkpoint_on_off,port_name_legacy)
%function metis4muscle3(debug_level,checkpoint_on_off,overlay_on_off)
%function metis4muscle3(debug_level,checkpoint_on_off,overlay_on_off,operator_init,operator_input,operator_output)



% all data neede for restart (checkpoint) should be in internal_state
% persistant variable
persistent internal_state
% initialisation on main substructure
internal_state.metis    = [];
internal_state.workflow = [];
internal_state.actor    = [];


% test of inputs
if (nargin == 0) || isempty(debug_level)
    debug_level = 3;
    % later it will be set to 0 in production
end
if debug_level < 0
    internal_state.profiler_status = true; 
    debug_level = abs(debug_level);
    profile off
    profile clear
    profile on
else
    internal_state.workflow.profiler_status = false;     
end
if debug_level > 0
    diary debug_metis_in_muscle3.txt
    diary on
end
if debug_level > 2
    dbstop if all error
end
if debug_level > 3
    keyboard;
end
internal_state.workflow.debug_level = debug_level;
%
if nargin < 2 || isempty(checkpoint_on_off)
    internal_state.workflow.checkpoint_on_off = false;
else
    internal_state.workflow.checkpoint_on_off = checkpoint_on_off;
end  
% the flag is now imposed in the actor as the other option is not used
internal_state.workflow.overlay_on_off = true;
% if nargin < 3 || isempty(overlay_on_off)
%     internal_state.workflow.overlay_on_off = true;
% else
%     internal_state.workflow.overlay_on_off = overlay_on_off;
% end   
% if nargin < 4 || isempty(operator_init)
%     operator_init = '';
% end
% if nargin < 5 || isempty(operator_input)
%     operator_input = '';
% end
% if nargin < 6 || isempty(operator_output)
%      operator_output = '';
% end
if nargin < 3 || isempty(port_name_legacy)
    port_name_legacy = false;
elseif port_name_legacy
    disp('Using legacy port naming "_input_" and "_output_" instead of "_in_" and "_out_"');
end

% intitial comment in diary file
fprintf('METIS bridge with Muscle3 initialization @ %s \n',datestr(now,0))

% set path and environement
% set environement for METIS if needed
if isempty(which('metis4imas'))
    addpath(fileparts(fileparts(fileparts(fileparts(which(mfilename))))));
    zineb_path;
    addpath(fileparts(which(mfilename)));
    disp('Used path is:')
    path
end
if debug_level > 1
    addpath(fullfile(fileparts(which(mfilename)),'finish'));
end 
% test python version from module
test_matlab_python_engine_version;
% declare muscle3 as non-codegen compatible
coder.extrinsic('py.libmuscle.Instance');
% set python environement
setdefault(py.os.environ,'MUSCLE_INSTANCE','macro');



% actor description
internal_state.actor.settings.list = metis4muscle3_setting_list;
%internal_state.actor.io.list       = metis4muscle3_port_list(operator_init,operator_input,operator_output);
if port_name_legacy
    internal_state.actor.io.list       = metis4muscle3_port_list_legacy;
else
    internal_state.actor.io.list       = metis4muscle3_port_list;
end

% open information logger in Python environement 
internal_state.actor.logger = open_python_logger_for_muscle3("Starting Metis actor");
% define instance flags
flags = instance_flags;
% Create libmuscle instance, specify all Ports
% Ports need to match the MUSCLE3 yaml description on the other side !
internal_state.actor.io.ports = metis4muscle3_make_ports(internal_state.actor.io);
% write instance object to matlab workspace. Currently, there is no other
% way than to retrieve it from there within simulink with nasty 'evalin'
% statements.
% only populate when not yet there
internal_state.actor.logger.info("Creation of muscle3 instance")
if ~isfield(internal_state.actor,'instance') || isempty(internal_state.actor.instance)
    if internal_state.workflow.checkpoint_on_off
        internal_state.actor.instance = py.libmuscle.Instance(internal_state.actor.io.ports,flags.USES_CHECKPOINT_API);
    else
        internal_state.actor.instance = py.libmuscle.Instance(internal_state.actor.io.ports);
    end
end

% test if instance exist
if ~isempty(internal_state.actor.instance)
    switch class(internal_state.actor.instance)
        case 'py.NoneType'
            internal_state.actor.logger.info("Unable to create Muscle3 instance")
            error('Communication with MUSCLE3 failed')            
        otherwise
            internal_state.actor.logger.info("Muscle3 instance created")
    end
else
    internal_state.actor.logger.info("Unable to create Muscle3 instance")
    error('Communication with MUSCLE3 failed')
end

% if ~internal_state.actor.instance.reuse_instance()
%     internal_state.actor.logger.info("Unable to create Muscle3 instance")
%     error('Communication with MUSCLE3 failed')
% else
%     internal_state.actor.logger.info("Muscle3 instance created")  
% end

% external data 
internal_state.external_data.data  = [];
internal_state.external_data.LCFS_nbpts  = NaN;
internal_state.external_data.list  = metis4muscle3_listofexternal_data;

% retrieve data from Muscle3
% get setting
internal_state.actor.logger.info("get parameters:")  
internal_state.actor.settings.value = get_muscle3_settings(internal_state.actor.instance,internal_state.actor.settings.list,'metis_option','metis_external_data',internal_state.external_data.list);
% get connected ports
internal_state.actor.logger.info("get list of connected ports:")  
internal_state.actor.io.status = metis4muscle3_port_status(internal_state.actor.instance,internal_state.actor.io);

% state of external data 
internal_state.external_data.state = get_external_data_state(internal_state.actor.io,internal_state.actor.settings.value);

% default workflow data and state
internal_state.workflow.clock.dt     =  NaN;
internal_state.workflow.clock.t_cur  =  internal_state.actor.settings.value.metis.clock.time.start;
internal_state.workflow.clock.t_next =  NaN;
internal_state.workflow.used_ports   = 'input';
internal_state.workflow.LCFS_nbpts   =  NaN;
% to keep memory of the r0 normalisation for magnetic rigidity during all the life of the actor
% as this factor is not stored in time dependant vector
internal_state.workflow.r0 = [];

% For mechanism mimicking sawtooth (used with NICE equilibrium + current
% diffusion).
internal_state.workflow.sawtooth.previous.delta_eta = [];
internal_state.workflow.sawtooth.previous.delta_jni = [];
internal_state.workflow.sawtooth.previous.x_ext     = [];
internal_state.workflow.sawtooth.next.delta_eta     = [];
internal_state.workflow.sawtooth.next.delta_jni     = [];
internal_state.workflow.sawtooth.next.x_ext         = [];
internal_state.workflow.sawtooth.LAST_ST_TIME       = -Inf;
internal_state.workflow.sawtooth.INDICE_INV_ST      = 0;
internal_state.workflow.sawtooth.MINIMAL_TAU_ST     = Inf;
internal_state.workflow.sawtooth.on_off             = 'off';

% default METIS data
internal_state.metis.loaded_data.post      = [];
internal_state.metis.loaded_data.ind_tcur  = NaN;
internal_state.metis.input                 = [];
internal_state.metis.LCFS_nbpts            = NaN;
internal_state.metis.output.zs             = [];
internal_state.metis.output.profil         = [];
internal_state.metis.output.z0dstruct      = [];   % METIS internal state
internal_state.metis.idses.error_flag      = 1;

% IMAS output data
internal_state.imas.ids_data = [];

% workflow state (initialisation)
internal_state.workflow.state = metis4muscle3_state_initialisation;
% reuse method on initialisation if needed
if internal_state.actor.settings.value.metis.synchronization.reuse
    %
    if ~internal_state.actor.instance.reuse_instance()
        internal_state.workflow.state.metis_exit = true;
    elseif internal_state.workflow.checkpoint_on_off 
        % get instance state for resuming and other
        internal_state.actor.logger.info("test if resuming mode:")
        internal_state.workflow.state.resuming       = internal_state.actor.instance.resuming();
        internal_state.workflow.state.should_init    = internal_state.actor.instance.should_init();
    end
elseif isany_f_init_used(internal_state.actor,{'f_init','in_f','input_f'})
    internal_state.actor.logger.info('Initial call of reuse_instance method');
    fprintf('Initial call of reuse_instance method\n');    
    %
    if ~internal_state.actor.instance.reuse_instance()
        internal_state.workflow.state.metis_exit = true;
    elseif internal_state.workflow.checkpoint_on_off 
        % get instance state for resuming and other
        internal_state.actor.logger.info("test if resuming mode:")
        internal_state.workflow.state.resuming       = internal_state.actor.instance.resuming();
        internal_state.workflow.state.should_init    = internal_state.actor.instance.should_init();
    end
else
    internal_state.actor.logger.info('No initial call of reuse_instance method at METIS actor initialization');
    fprintf('No initial call of reuse_instance method at METIS actor initialization\n');
end

if ~internal_state.actor.settings.value.metis.computation
    if isempty(internal_state.actor.settings.value.metis.filename)
        internal_state.actor.logger.info("metis_filename should define if setting metis_computation is 0")
        error('Empty METIS file name');
    elseif exist(internal_state.actor.settings.value.metis.filename,'file')
        internal_state.workflow.state.metis_computation = 0;
    else
        internal_state.actor.logger.info("metis_filename should point to a existing  METIS simulation file if setting metis_computation is 0")
        error('Invalid METIS file name');
    end
else
        internal_state.workflow.state.metis_computation = 1;    
end
% working point
if internal_state.actor.settings.value.metis.workingpoint
   internal_state.workflow.state.working_point = true;
else
   internal_state.workflow.state.working_point = false;  
end
                                      
% plot actor layout
metis4muscle3_layout(internal_state)

% start loop on instance existance
while ~internal_state.workflow.state.metis_exit
    
    % itialisation of METIS
    if internal_state.workflow.state.should_init
        % force reinitialisation of METIS
        internal_state.actor.logger.info("Force METIS reinitialisation from scratch")
        internal_state.workflow.state.metis_init = false;
    end
    
    % use either init or input ports depending of state
    if  internal_state.workflow.state.metis_init
        if internal_state.actor.settings.value.metis.alwaysreadinitport
            if isany_f_init_used(internal_state.actor,{'f_init'})
                if port_name_legacy
                    internal_state.workflow.used_ports   = {'f_init','input_f','input_s'};
                else
                    internal_state.workflow.used_ports   = {'f_init','in_f','in_s'};
                end
                
            elseif port_name_legacy
                internal_state.workflow.used_ports   = {'input_f','input_s'};
            else
                internal_state.workflow.used_ports   = {'in_f','in_s'};                
            end
        elseif port_name_legacy
            internal_state.workflow.used_ports   = {'input_f','input_s'};   %  now any input
        else
            internal_state.workflow.used_ports   = {'in_f','in_s'};   %  now any input            
        end
    elseif isany_f_init_used(internal_state.actor,{'f_init'})
        internal_state.workflow.used_ports   = {'f_init'};
    elseif isany_f_init_used(internal_state.actor,{'input_f','in_f'})
        if port_name_legacy
            internal_state.workflow.used_ports   = {'input_f'};
        else
            internal_state.workflow.used_ports   = {'in_f'};
        end
    elseif port_name_legacy
        internal_state.workflow.used_ports   = {'input_s'}; 
    else
        internal_state.workflow.used_ports   = {'in_s'};         
    end
    
    % get sawtooth mechanism state
    internal_state.workflow.sawtooth.on_off = get_st_state(internal_state.actor.io.status, ...
                                                           internal_state.actor.settings.value.metis, ...
                                                           internal_state.external_data.state, ...
                                                            internal_state.workflow.used_ports);

    
    % read messages and update the new one
    if ~internal_state.workflow.state.resuming
        if  internal_state.workflow.state.metis_init
            % read input ports
            internal_state.actor.logger.info("test if input port should be read:")
            if internal_state.actor.settings.value.metis.alwaysreadinitport
                internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'init', internal_state.workflow.overlay_on_off,internal_state.actor.io.data);             
            end
            internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'input', internal_state.workflow.overlay_on_off,internal_state.actor.io.data);
        else
            % most of the case init port should be read, but if no init
            % port is connected, read input ports in backup
            switch internal_state.workflow.used_ports{1}
                case 'f_init'
                    % read init ports
                    internal_state.actor.logger.info("test if init port should be read:")
                    internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'init', internal_state.workflow.overlay_on_off);
                otherwise
                    % read input ports in backup
                    internal_state.actor.logger.info("No init port connected: in backup,test if input port should be read:")
                    internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'init', internal_state.workflow.overlay_on_off);
                    
            end
        end
        if debug_level > 1
            %display message contents
            noms = fieldnames(internal_state.actor.io.data);
            for l=1:length(noms)
                if ~isempty(internal_state.actor.io.data.(noms{l}).data_raw)
                    disp('--------------------------------------');
                    fprintf('Message content for %s:\n',noms{l});
                    internal_state.actor.io.data.(noms{l})
                    disp('with data structure:')
                    internal_state.actor.io.data.(noms{l}).data
                end
            end
            
        end
    else
        if internal_state.workflow.state.should_init
            if  internal_state.workflow.state.metis_init
                % read input ports
                internal_state.actor.logger.info("read input ports on resuming")
                if internal_state.actor.settings.value.metis.alwaysreadinitport
                    internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'init', internal_state.workflow.overlay_on_off,internal_state.actor.io.data);
                end
                internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'input', internal_state.workflow.overlay_on_off,internal_state.actor.io.data);
            else
                % read init ports
                internal_state.actor.logger.info("read init ports on resuming")
                internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'init', internal_state.workflow.overlay_on_off);
            end
            if debug_level > 1
                %display message contents
                noms = fieldnames(internal_state.actor.io.data);
                for l=1:length(noms)
                    if ~isempty(internal_state.actor.io.data.(noms{l}).data_raw)
                        disp('--------------------------------------');
                        fprintf('Message content for %s:\n',noms{l});
                        internal_state.actor.io.data.(noms{l})
                        disp('with data structure:')
                        internal_state.actor.io.data.(noms{l}).data
                    end
                end
                
            end
        end
        
    end
    
    % get Mode
    internal_state.actor.logger.info("Detecting actor running mode:")
    internal_state.workflow.state = get_workflow_mode(internal_state.workflow.state,internal_state.actor,internal_state.workflow.used_ports);
    
    % for other mode than replay
    if ~internal_state.workflow.state.metis_replay
        % get current time and next time depending on actor configuration
        [internal_state.workflow.clock.t_cur, ...
            internal_state.workflow.clock.t_next, ...
            internal_state.workflow.clock.dt] = ...
            metis4muscle3_time_operator(internal_state);
    end
    
    %
    % section for METIS call in different mode depending of connected ports:
    %    * no port connected -> read a METIS simulation file and stream it
    %      on outputed ports.
    %
    %    * only port clock connected -> read a METIS simulation file,
    %      compute metis results in evolution mode using METIS setting at
    %      the requested time and stream results on outputed ports.
    %
    %    * if pusle_schedule port is connected -> get data from
    %      pulse_schedule, compute metis results in evolution mode using
    %      and stream results on outputed ports.
    %
    %    If other ports are connected, data are used to create external data
    %    for METIS depening of the settings (that included external
    %    equilibrium).
    %
    %    If METIS computation mode is switch off, the data from the METIS
    %    simulation will be interpolated of request time whatever is the
    %    input conneted ports.
    %
    if internal_state.workflow.state.resuming
        internal_state.actor.logger.info("Resuming simulation from previous run")
        % load previous internal_state
% for testing        
%         rep =load('snapshot_message.mat');
%         rep.timestamp = 32;
%         rep.next_timestamp = 34;
%         resume_message = rep;
        resume_message = internal_state.actor.instance.load_snapshot();
        internal_state = metis4muscle3_resume_simulation(internal_state, ...
                         resume_message,internal_state.actor.settings.value.metis.updatetimeonrestart);
        % restored external data for METIS
        metis4muscle3_restore_external_data(internal_state);
        
        % change states
        internal_state.workflow.state.metis_init = true;
        internal_state.workflow.state.resuming   = false;
        internal_state.workflow.state.should_init = false;
        
        % this is the last action  for this loop as no computation will be
        % done on resume
        
    elseif internal_state.workflow.state.metis_replay
        % load METIS file
        if isempty(internal_state.metis.loaded_data.post)
            if exist(internal_state.actor.settings.value.metis.filename,'file')
                internal_state.metis.loaded_data.post = load_metis_imas(internal_state.actor.settings.value.metis.filename);
                internal_state.actor.logger.info(sprintf("METIS file %s has been loaded and METIS data will be stream in replay mode",...
                    internal_state.actor.settings.value.metis.filename));
                fprintf("METIS file %s has been loaded and METIS data will be stream in replay mode\n",...
                    internal_state.actor.settings.value.metis.filename);
                % get inital time
                internal_state.metis.loaded_data.ind_tcur = find(internal_state.metis.loaded_data.post.z0dinput.cons.temps >= ...
                    max(internal_state.actor.settings.value.metis.clock.time.start, ...
                    internal_state.metis.loaded_data.post.z0dinput.cons.temps(1)),1);
            else
                internal_state.actor.logger.info("In replay mode, METIS simulation file should be exist and ")
                error("replay mode: METIS data file is not defined or do not exist");
            end
        else
            % update time
            internal_state.metis.loaded_data.ind_tcur = internal_state.metis.loaded_data.ind_tcur + 1;
            % + exit condition
            if internal_state.metis.loaded_data.ind_tcur >= length( internal_state.metis.loaded_data.post.z0dinput.cons.temps)
                % exit at the next while loop
                internal_state.workflow.state.metis_exit = true;
            end
        end
        internal_state.workflow.clock.t_cur = internal_state.metis.loaded_data.post.z0dinput.cons.temps(internal_state.metis.loaded_data.ind_tcur);
        if  internal_state.metis.loaded_data.ind_tcur < length( internal_state.metis.loaded_data.post.z0dinput.cons.temps)
            internal_state.workflow.clock.t_next =  internal_state.metis.loaded_data.post.z0dinput.cons.temps( ...
                internal_state.metis.loaded_data.ind_tcur + 1);
        else
            internal_state.workflow.clock.t_next = Inf;
        end
        internal_state.workflow.clock.dt = internal_state.workflow.clock.t_next - internal_state.workflow.clock.t_cur;
        %
        % extract data at this time
        internal_state.actor.logger.info(sprintf("METIS replay data at %g s",internal_state.workflow.clock.t_cur));
        fprintf("METIS replay data at %g s\n",internal_state.workflow.clock.t_cur);
        [internal_state.metis.output.zs,internal_state.metis.output.profil,internal_state.metis.output.z0dstruct] =  ...
            make_z0dstruct_from_post(internal_state.metis.loaded_data.post,internal_state.metis.loaded_data.ind_tcur);
        
        % change state
        internal_state.workflow.state.metis_init  = true;
        internal_state.workflow.state.resuming    = false;
        internal_state.workflow.state.should_init = false;
        
    elseif ~internal_state.workflow.state.metis_computation
        % not computing METIS, just stream interpolated data
        % load METIS file
        if isempty(internal_state.metis.loaded_data.post)
            if exist(internal_state.actor.settings.value.metis.filename,'file')
                internal_state.metis.loaded_data.post = load_metis_imas(internal_state.actor.settings.value.metis.filename);
            else
                internal_state.actor.logger.info("In replay mode, METIS simulation file should be exist and ")
                error("replay mode: METIS data file is not defined or do not exist");
            end
        end
        % interpolate data on requested time
        internal_state.actor.logger.info(sprintf("Streaming METIS interpolated data at %g s",internal_state.workflow.clock.t_cur));
        [internal_state.metis.output.zs,internal_state.metis.output.profil,internal_state.metis.output.z0dstruct] =  ...
            metis4muscle3_interpol(internal_state.metis.loaded_data.post,internal_state.workflow.clock.t_cur);
        
        % change state
        internal_state.workflow.state.metis_init  = true;
        internal_state.workflow.state.resuming    = false;
        internal_state.workflow.state.should_init = false;
   
    elseif internal_state.workflow.state.working_point
        % METIS is not yet initialized
        internal_state.actor.logger.info("Launch METIS in working point mode")
        disp("Launch METIS in working point mode");
        % load METIS file
        if isempty(internal_state.metis.loaded_data.post)
            if exist(internal_state.actor.settings.value.metis.filename,'file')
                internal_state.metis.loaded_data.post = load_metis_imas(internal_state.actor.settings.value.metis.filename);
            end
        end
        
        % decode pulse_schedule if available, otherwise load METIS file
        [internal_state.metis.input,internal_state.metis.LCFS_nbpts] = get_metis_input_from_muscle3(internal_state.metis,internal_state.actor,internal_state.workflow.clock,internal_state.workflow.used_ports);
        % make initial external data for METIS depending on port connected to the actor
        [internal_state.external_data.data,internal_state.external_data.LCFS_nbpts]  = metis4muscle3_init_external_data(internal_state);
        % here the change in number of points in LCFS should be handled !
        if isfinite(internal_state.metis.LCFS_nbpts) && isfinite(internal_state.external_data.LCFS_nbpts) && ...
                (internal_state.metis.LCFS_nbpts ~= internal_state.external_data.LCFS_nbpts)
                    [internal_state.metis.input.sepa1t,internal_state.metis.output.z0dstruct,internal_state.workflow.LCFS_nbpts] =  ...
                        manage_LCFS_nbpoints(internal_state.metis.input.sepa1t, ...
                                             internal_state.metis.output.z0dstruct, ...
                                             internal_state.metis.LCFS_nbpts, ...
                                             internal_state.external_data.LCFS_nbpts);
        end
        % for external equilibrium
        internal_state.external_data.data = manage_LCFS_external_equilibrium(internal_state.external_data.data,internal_state.workflow.LCFS_nbpts);
        
        % update cons and geo if external equilibrium (
        if isfield(internal_state.external_data.data,'equilibrium') && ~isempty(internal_state.external_data.data.equilibrium)
           %  update geo and cons
           internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts);
        elseif isfield( internal_state.external_data.data,'current_diffusion') && ~isempty(internal_state.external_data.data.current_diffusion)
           % update ip only 
           internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts,'only_ip');
        end

        % call METIS in working point mode 
        [internal_state.metis.output.zs,internal_state.metis.output.profil, ...
                internal_state.metis.output.z0dstruct] = metis4muscle3_working_point( ...
                                        internal_state.metis.input.option, ...
                                        internal_state.metis.input.time, ...
                                        internal_state.metis.input.cons1t, ...
                                        internal_state.metis.input.geo1t, ...
                                        internal_state.metis.input.sepa1t);

        
        % UAL option is not always preserve (security rule on parameters)
        noms = fieldnames(internal_state.metis.input.option);
        for kl=1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{kl})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{kl}) = internal_state.metis.input.option.(noms{kl});
            end
        end
        % compatibilite ascendante
        if ~isfield(internal_state.metis.output.z0dstruct,'zs') && isfield(internal_state.metis.output.z0dstruct,'zerod')
            internal_state.metis.output.z0dstruct.zs = internal_state.metis.output.z0dstruct.zerod;
        elseif isfield(internal_state.metis.output.z0dstruct,'zs')
            internal_state.metis.output.z0dstruct.zerod = internal_state.metis.output.z0dstruct.zs;
        end
        if ~isfield(internal_state.metis.output.z0dstruct,'profil') && isfield(internal_state.metis.output.z0dstruct,'profil0d')
            internal_state.metis.output.z0dstruct.profil = internal_state.metis.output.z0dstruct.profil0d;
        elseif isfield(internal_state.metis.output.z0dstruct,'profil')
            internal_state.metis.output.z0dstruct.profil0d = internal_state.metis.output.z0dstruct.profil;
        end
        
        % add possible missign field in option
        model = metis4imas(1);
        noms = fieldnames(model.valeur);
        for k = 1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{k})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
            end
        end
        internal_state.metis.output.z0dstruct.z0dinput.option = secure_option_imas(internal_state.metis.output.z0dstruct.z0dinput.option);
        
        % change state (should be updated)
        internal_state.workflow.state.metis_init  = true;
        internal_state.workflow.state.resuming    = false;
        internal_state.workflow.state.should_init = false;

    elseif ~internal_state.workflow.state.metis_init
        % METIS is not yet initialized
        internal_state.actor.logger.info("Starting METIS initialisation")
        disp("Starting METIS initialisation");
        % load METIS file if available
        if isempty(internal_state.metis.loaded_data.post)
            if exist(internal_state.actor.settings.value.metis.filename,'file')
                internal_state.metis.loaded_data.post = load_metis_imas(internal_state.actor.settings.value.metis.filename);
            end
        end
        
        % decode pulse_schedule if available, otherwise load METIS file
        [internal_state.metis.input,internal_state.metis.LCFS_nbpts] = get_metis_input_from_muscle3(internal_state.metis,internal_state.actor,internal_state.workflow.clock,internal_state.workflow.used_ports);
        % make initial external data for METIS depending on port connected to the actor
        [internal_state.external_data.data,internal_state.external_data.LCFS_nbpts]  = metis4muscle3_init_external_data(internal_state);
        % here the change in number of points in LCFS should be handled !
        if isfinite(internal_state.metis.LCFS_nbpts) && isfinite(internal_state.external_data.LCFS_nbpts) && ...
                (internal_state.metis.LCFS_nbpts ~= internal_state.external_data.LCFS_nbpts)
                    [internal_state.metis.input.sepa1t,internal_state.metis.output.z0dstruct,internal_state.workflow.LCFS_nbpts] =  ...
                        manage_LCFS_nbpoints(internal_state.metis.input.sepa1t, ...
                                             internal_state.metis.output.z0dstruct, ...
                                             internal_state.metis.LCFS_nbpts, ...
                                             internal_state.external_data.LCFS_nbpts);
        end
        % for external equilibrium
        internal_state.external_data.data = manage_LCFS_external_equilibrium(internal_state.external_data.data,internal_state.workflow.LCFS_nbpts);
        
        % update cons and geo if external equilibrium (
        if isfield(internal_state.external_data.data,'equilibrium') && ~isempty(internal_state.external_data.data.equilibrium)
           %  update geo and cons
           internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts);
        elseif isfield( internal_state.external_data.data,'current_diffusion') && ~isempty(internal_state.external_data.data.current_diffusion)
           % update ip only 
           internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts,'only_ip');
        end
       

        
        % should add call start depending on settings
        if internal_state.actor.settings.value.metis.coldinit
            internal_state.actor.logger.info('Cold initialisation of METIS');
            disp('Cold initialisation of METIS');
            % interpolate data on requested time
            internal_state.actor.logger.info(sprintf("Streaming METIS interpolated data at %g s",internal_state.workflow.clock.t_cur));
            [internal_state.metis.output.zs,internal_state.metis.output.profil,internal_state.metis.output.z0dstruct] =  ...
                metis4muscle3_interpol(internal_state.metis.loaded_data.post,internal_state.workflow.clock.t_cur);
            
            % call METIS evolution in initialisation mode
        elseif ~isempty(internal_state.metis.input.sepa1t)
            metis4muscle3_external_data_info(internal_state.actor.logger);
            %
            [internal_state.metis.output.zs,internal_state.metis.output.profil, ...
                internal_state.metis.output.z0dstruct] = zerodevolution([], ...
                internal_state.metis.input.option, ...
                internal_state.metis.input.time, ...
                internal_state.metis.input.cons1t, ...
                internal_state.metis.input.geo1t,[], ...
                internal_state.metis.input.sepa1t);
        else
            metis4muscle3_external_data_info(internal_state.actor.logger);
            %
            [internal_state.metis.output.zs,internal_state.metis.output.profil, ...
                internal_state.metis.output.z0dstruct] = zerodevolution([], ...
                internal_state.metis.input.option, ...
                internal_state.metis.input.time, ...
                internal_state.metis.input.cons1t, ...
                internal_state.metis.input.geo1t,[]);
        end
        % UAL option is not always preserve (security rule on parameters)
        noms = fieldnames(internal_state.metis.input.option);
        for kl=1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{kl})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{kl}) = internal_state.metis.input.option.(noms{kl});
            end
        end
        % compatibilite ascendante
        if ~isfield(internal_state.metis.output.z0dstruct,'zs') && isfield(internal_state.metis.output.z0dstruct,'zerod')
            internal_state.metis.output.z0dstruct.zs = internal_state.metis.output.z0dstruct.zerod;
        elseif isfield(internal_state.metis.output.z0dstruct,'zs')
            internal_state.metis.output.z0dstruct.zerod = internal_state.metis.output.z0dstruct.zs;
        end
        if ~isfield(internal_state.metis.output.z0dstruct,'profil') && isfield(internal_state.metis.output.z0dstruct,'profil0d')
            internal_state.metis.output.z0dstruct.profil = internal_state.metis.output.z0dstruct.profil0d;
        elseif isfield(internal_state.metis.output.z0dstruct,'profil')
            internal_state.metis.output.z0dstruct.profil0d = internal_state.metis.output.z0dstruct.profil;
        end
        
        % add possible missign field in option
        model = metis4imas(1);
        noms = fieldnames(model.valeur);
        for k = 1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{k})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
            end
        end
        internal_state.metis.output.z0dstruct.z0dinput.option = secure_option_imas(internal_state.metis.output.z0dstruct.z0dinput.option);
        
        % change state
        internal_state.workflow.state.metis_init  = true;
        internal_state.workflow.state.resuming    = false;
        internal_state.workflow.state.should_init = false;
        
    else
        % compute METIS evolution or convergence
        % decode pulse_schedule if available, otherwise use METIS data
        % loaded at initialisation to get current reference and configuration
        [internal_state.metis.input,internal_state.metis.LCFS_nbpts] = get_metis_input_from_muscle3(internal_state.metis,internal_state.actor,internal_state.workflow.clock,internal_state.workflow.used_ports);
        
        
        % time evolve ?
        % if internal_state last time is older than new t_cur, this is a
        % evolution phase otherwise it is a convergence loop
        if internal_state.metis.output.zs.temps > internal_state.workflow.clock.t_cur
            internal_state.workflow.state.metis_rewind      = true;
            internal_state.workflow.state.metis_convergence = false;
            internal_state.workflow.state.metis_evolution   = false;
            fprintf('METIS_EVOLUTION will rewind from t= %g s to t = %g s .\n',  ...
                internal_state.metis.output.zs.temps, ...
                internal_state.workflow.clock.t_cur);
            internal_state.actor.logger.info( ...
                sprintf('METIS_EVOLUTION will rewind from t= %g s to t = %g s .\n',  ...
                internal_state.metis.output.zs.temps, ...
                internal_state.workflow.clock.t_cur));
            
        elseif internal_state.metis.output.zs.temps == internal_state.workflow.clock.t_cur
            internal_state.workflow.state.metis_rewind      = false;
            internal_state.workflow.state.metis_convergence = true;
            internal_state.workflow.state.metis_evolution   = false;
            fprintf('METIS_EVOLUTION will perform a convergence step at  t= %g s .\n',  ...
                internal_state.metis.output.zs.temps);
            internal_state.actor.logger.info( ...
                sprintf('METIS_EVOLUTION will perform a convergence step at  t= %g s .\n',  ...
                internal_state.metis.output.zs.temps));
        else
            internal_state.workflow.state.metis_rewind      = false;
            internal_state.workflow.state.metis_convergence = false;
            internal_state.workflow.state.metis_evolution   = true;
            fprintf('METIS_EVOLUTION will compute a evolution from t= %g s to t = %g s .\n',  ...
                internal_state.metis.output.zs.temps, ...
                internal_state.workflow.clock.t_cur);
            internal_state.actor.logger.info( ...
                sprintf('METIS_EVOLUTION will compute a evolution from t= %g s to t = %g s .\n',  ...
                internal_state.metis.output.zs.temps, ...
                internal_state.workflow.clock.t_cur));
        end
        
        % mechanism mimicking sawtooth when METIS is coupled with NICE
        % solving G-S and dFdia/dt
        switch internal_state.workflow.sawtooth.on_off
            case 'on'
                internal_state.workflow.sawtooth = metis4muscle3_mimic_sawtooth(internal_state);
                % METIS outputed dans are modified juste before IDSs creation
        end
        
        % update external data for METIS depending on port connected to the actor
        [internal_state.external_data.data,internal_state.external_data.LCFS_nbpts] =  metis4muscle3_update_external_data(internal_state);
        
        % here the change in number of points in LCFS should be handled !
        [internal_state.metis.input.sepa1t,internal_state.metis.output.z0dstruct,internal_state.workflow.LCFS_nbpts] =  ...
                      manage_LCFS_nbpoints(internal_state.metis.input.sepa1t, ...
                      internal_state.metis.output.z0dstruct, ...
                      internal_state.metis.LCFS_nbpts,internal_state.external_data.LCFS_nbpts);
         % for external equilibrium
         internal_state.external_data.data = manage_LCFS_external_equilibrium(internal_state.external_data.data,internal_state.workflow.LCFS_nbpts);
                  
         % update cons and geo if external equilibrium (
         if isfield(internal_state.external_data.data,'equilibrium') && ~isempty(internal_state.external_data.data.equilibrium)
             %  update geo and cons
             internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts);
         elseif isfield( internal_state.external_data.data,'current_diffusion') && ~isempty( internal_state.external_data.data.current_diffusion)
             % update ip only
             internal_state.metis.input = metis4muscle3_update_metis_input_on_equilibrium(internal_state.metis.input,internal_state.external_data.data.equilibrium,internal_state.workflow.LCFS_nbpts,'only_ip');
         end
        
        % call METIS evolution (evolution or convergence is already handle
        % internaly by zerodevolution).
        if ~isempty(internal_state.metis.input.sepa1t)
            metis4muscle3_external_data_info(internal_state.actor.logger);
            %
            [internal_state.metis.output.zs,internal_state.metis.output.profil, ...
                internal_state.metis.output.z0dstruct] = zerodevolution(internal_state.metis.output.z0dstruct, ...
                internal_state.metis.input.option, ...
                internal_state.metis.input.time, ...
                internal_state.metis.input.cons1t, ...
                internal_state.metis.input.geo1t,[], ...
                internal_state.metis.input.sepa1t);
        else
            metis4muscle3_external_data_info(internal_state.actor.logger);
            %
            [internal_state.metis.output.zs,internal_state.metis.output.profil, ...
                internal_state.metis.output.z0dstruct] = zerodevolution(internal_state.metis.output.z0dstruct, ...
                internal_state.metis.input.option, ...
                internal_state.metis.input.time, ...
                internal_state.metis.input.cons1t, ...
                internal_state.metis.input.geo1t,[]);
        end
        % UAL option is not always preserve (security rule on parameters)
        noms = fieldnames(internal_state.metis.input.option);
        for kl=1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{kl})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{kl}) = internal_state.metis.input.option.(noms{kl});
            end
        end
        % compatibilite ascendante
        if ~isfield(internal_state.metis.output.z0dstruct,'zs') && isfield(internal_state.metis.output.z0dstruct,'zerod')
            internal_state.metis.output.z0dstruct.zs = internal_state.metis.output.z0dstruct.zerod;
        elseif isfield(internal_state.metis.output.z0dstruct,'zs')
            internal_state.metis.output.z0dstruct.zerod = internal_state.metis.output.z0dstruct.zs;
        end
        if ~isfield(internal_state.metis.output.z0dstruct,'profil') && isfield(internal_state.metis.output.z0dstruct,'profil0d')
            internal_state.metis.output.z0dstruct.profil = internal_state.metis.output.z0dstruct.profil0d;
        elseif isfield(internal_state.metis.output.z0dstruct,'profil')
            internal_state.metis.output.z0dstruct.profil0d = internal_state.metis.output.z0dstruct.profil;
        end
        
        % add possible missign field in option
        model = metis4imas(1);
        noms = fieldnames(model.valeur);
        for k = 1:length(noms)
            if ~isfield(internal_state.metis.output.z0dstruct.z0dinput.option,noms{k})
                internal_state.metis.output.z0dstruct.z0dinput.option.(noms{k}) = model.valeur.(noms{k});
            end
        end
        internal_state.metis.output.z0dstruct.z0dinput.option = secure_option_imas(internal_state.metis.output.z0dstruct.z0dinput.option);
                
    end

    % get value of r0 at first pass
    if isempty(internal_state.workflow.r0)
        internal_state.workflow.r0  =  mean(internal_state.metis.output.z0dstruct.z0dinput.geo.R);
    end
    
    % here should be IDS creation depending on connected ports using metis4imas
    % secure option for AL if needed
    internal_state.metis.output.z0dstruct.z0dinput.option = set_missing_al_parameters(internal_state.metis.output.z0dstruct.z0dinput.option);
    % compatibility with older tools that can be involve in IMAS METIS
    if ~isfield(internal_state.metis.output.z0dstruct,'zs') && isfield(internal_state.metis.output.z0dstruct,'zerod')
        internal_state.metis.output.z0dstruct.zs = internal_state.metis.output.z0dstruct.zerod;
    elseif  isfield(internal_state.metis.output.z0dstruct,'zs')
        internal_state.metis.output.z0dstruct.zerod = internal_state.metis.output.z0dstruct.zs;
    end
    if ~isfield(internal_state.metis.output.z0dstruct,'profil') && isfield(internal_state.metis.output.z0dstruct,'profil0d')
        internal_state.metis.output.z0dstruct.profil = internal_state.metis.output.z0dstruct.profil0d;
    elseif  isfield(internal_state.metis.output.z0dstruct,'profil')
        internal_state.metis.output.z0dstruct.profil0d = internal_state.metis.output.z0dstruct.profil;
    end
    % set list opf generated IDSes by metis4imas
    internal_state.metis.output.z0dstruct.z0dinput.option = switch_al_parameters(internal_state.metis.output.z0dstruct.z0dinput.option,internal_state.actor.io.status,port_name_legacy);
    % generate IDSes data structures
    internal_state.actor.logger.info("Generating IDSes messages")
    disp("Generating IDSes messages")
    dbstate = dbstatus;
    dbclear all
    [internal_state.metis.idses.error_flag,internal_state.imas.ids_data] =  ...
        metis4muscle3_make_idses(internal_state.metis.output.zs, ...
        internal_state.metis.output.profil, ...
        internal_state.metis.output.z0dstruct, ...
        internal_state.actor.settings.value.metis.psioffset, ...
        internal_state.workflow.sawtooth,internal_state.workflow.r0);
    if ~isempty(dbstate)
        dbstop if all error
    end
   
   
    % endding on time condition
    if internal_state.workflow.clock.t_cur > internal_state.actor.settings.value.metis.clock.time.stop
        internal_state.workflow.state.metis_exit = true;
        % send information to muscle3 through message
        internal_state.workflow.clock.t_next = Inf;
    elseif ~isfinite(internal_state.workflow.clock.t_next)
        internal_state.workflow.state.metis_exit = true;
    end

    % here we send messages to the wild
    internal_state.actor.logger.info("Sending messages");
    disp("Sending messages");
    send_muscle3_messages(internal_state,internal_state.actor.settings.value.metis.synchronization.flush);
    
    % check and display reuse_instance state
    if isany_f_init_used(internal_state.actor,internal_state.workflow.used_ports)
        internal_state.actor.logger.info('at least one F_INIT port is connected: reuse_instance method will be call');
        disp('at least one F_INIT port is connected: reuse_instance method will be call');
    elseif internal_state.actor.settings.value.metis.synchronization.reuse
        internal_state.actor.logger.info('reuse_instance method call is forced');        
        disp('reuse_instance method call is forced');        
    end
    
    % update checkpoint commands (only in reuse loop)
    if isany_f_init_used(internal_state.actor,internal_state.workflow.used_ports) || internal_state.actor.settings.value.metis.synchronization.reuse 
        if internal_state.workflow.checkpoint_on_off && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.workflow.state.checkpoint   = internal_state.actor.instance.should_save_snapshot(internal_state.workflow.clock.t_cur);
            if ~internal_state.workflow.state.checkpoint
                internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
            end
        end
        % send checkpoint message if requested
        if internal_state.workflow.state.checkpoint && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.actor.logger.info(sprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur));
            fprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur);
            internal_state.actor.instance.save_snapshot(make_checkpoint_message(internal_state));
            internal_state.workflow.state.checkpoint = false;
            internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
        end
    else
        if internal_state.workflow.checkpoint_on_off && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.workflow.state.checkpoint   = internal_state.actor.instance.should_save_snapshot(internal_state.workflow.clock.t_cur);
        end
        % send checkpoint message if requested
        if internal_state.workflow.state.checkpoint && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.actor.logger.info(sprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur));
            fprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur);
            internal_state.actor.instance.save_snapshot(make_checkpoint_message(internal_state));
            internal_state.workflow.state.checkpoint = false;
        end
    end
    % it look like in MUSCLE3 spirit this commande should be here !!!!
    % finalisation by sending last message
    % test if final checkpoint is requested
    if internal_state.workflow.state.final_checkpoint && ~ internal_state.workflow.state.final_checkpoint_done
        internal_state.actor.logger.info(sprintf('saving final snapshot at t = %g s\n',internal_state.workflow.clock.t_cur));
        fprintf('saving final snapshot at t = %g s\n',internal_state.workflow.clock.t_cur);
        internal_state.actor.instance.save_final_snapshot(make_checkpoint_message(internal_state));
        internal_state.workflow.state.final_checkpoint = false;
        internal_state.workflow.state.final_checkpoint_done = true;
        if isany_f_init_used(internal_state.actor,internal_state.workflow.used_ports)
             internal_state.workflow.state.metis_exit = true;
        elseif ~isfinite(internal_state.workflow.clock.t_next)
              internal_state.workflow.state.metis_exit = true;
        end    
    % end of workflow by default if not another condition has been set
    elseif internal_state.workflow.state.metis_replay && (internal_state.actor.settings.value.metis.synchronization.flush == 1)
        % just do loop but should be not to fast to do not break muscle3
        % messaging: sleep and synchronise to prevent problem on overloaded machine
        drawnow;      % flush matlab buffer (including graphics)
        py.os.sync(); % flush systemp buffer
        pause(1);     % wait a bit
        if internal_state.workflow.clock.t_cur > internal_state.actor.settings.value.metis.clock.time.stop
            internal_state.workflow.state.metis_exit = true;
        elseif ~isfinite(internal_state.workflow.clock.t_next)
              internal_state.workflow.state.metis_exit = true;
        end
    end

    % reuse test if needed
    if internal_state.actor.settings.value.metis.synchronization.reuse 
        internal_state.actor.logger.info('Call of instance.reuse method');
        disp('Call of instance.reuse method');       
        %
        if ~internal_state.actor.instance.reuse_instance()
            internal_state.workflow.state.metis_exit = true;
            % check for resuming in reuse loop
        elseif internal_state.workflow.checkpoint_on_off
            % not use by METIS at this point but should be read as asked
            % by MUSCLE3            
            internal_state.actor.instance.resuming();
            internal_state.actor.instance.should_init();
        end
    elseif isany_f_init_used(internal_state.actor,internal_state.workflow.used_ports)
        internal_state.actor.logger.info('Call of reuse_instance method');
        disp('Call of reuse_instance method');       
        %
        if ~internal_state.actor.instance.reuse_instance()
            internal_state.workflow.state.metis_exit = true;
        % check for resuming in reuse loop
        elseif internal_state.workflow.checkpoint_on_off
            % not use by METIS at this point but should be read as asked
            % by MUSCLE3            
            internal_state.actor.instance.resuming();
            internal_state.actor.instance.should_init();
        end
    else
       internal_state.actor.logger.info('"reuse_instance" method was not call during this loop');
       disp('"reuse_instance" method was not call during this loop');       
    end    
    if internal_state.workflow.state.metis_exit 
        internal_state.actor.logger.info('End of METIS actor internal loop');
        disp('End of METIS actor interbal loop');       
    end
end


% update checkpoint commands (only in reuse loop): case with only one
% call of reuse.
if ~(isany_f_init_used(internal_state.actor,internal_state.workflow.used_ports) || internal_state.actor.settings.value.metis.synchronization.reuse)
    if internal_state.workflow.checkpoint_on_off && ~internal_state.workflow.state.final_checkpoint_done
        internal_state.workflow.state.checkpoint   = internal_state.actor.instance.should_save_snapshot(internal_state.workflow.clock.t_cur);
        if ~internal_state.workflow.state.checkpoint
            internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
        end
    end
    % send checkpoint message if requested
    if internal_state.workflow.state.checkpoint && ~internal_state.workflow.state.final_checkpoint_done
        internal_state.actor.logger.info(sprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur));
        fprintf('saving snapshot at t = %g s\n',internal_state.workflow.clock.t_cur);
        internal_state.actor.instance.save_snapshot(make_checkpoint_message(internal_state));
        internal_state.workflow.state.checkpoint = false;
        internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
    end
    if internal_state.workflow.state.final_checkpoint && ~ internal_state.workflow.state.final_checkpoint_done
        internal_state.actor.logger.info(sprintf('saving final snapshot at t = %g s\n',internal_state.workflow.clock.t_cur));
        fprintf('saving final snapshot at t = %g s\n',internal_state.workflow.clock.t_cur);
        internal_state.actor.instance.save_final_snapshot(make_checkpoint_message(internal_state));
        internal_state.workflow.state.final_checkpoint = false;
        internal_state.workflow.state.final_checkpoint_done = true;
    end   
end

% no more needed with new if condition at the end of the loop
% the previous one is not logical.
% % wait for end of message sending and reading by other actors
% if internal_state.workflow.state.metis_replay
%     internal_state.actor.instance.reuse_instance();
% end


% manage profiler results
profiler_status = profile('status');
switch profiler_status.ProfilerStatus
    case 'on'
        disp('Preparing profiler report');
        dbclear all
        profsave(profile('info'),'metis4muscle3_profile_results')
        profile report
        profile off
end

% then exit from Matlab
exit;









