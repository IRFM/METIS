%% METIS muscle bridge initialization script.
% R. Nouailletas 2024, CEA
% J-F. Artaud 2024, CEA
% based on work of T. Ravensbergen et M. Schneider 2023

% This script initializes MUSCLE3 (logging, ports),
% initializes the parameters, and then launches a simulation.
% The bridge receives a equilibrium IDS via MUSCLE3, and send a coreprofiles IDS
function metis4muscle3_data_provider(debug_level,operator_input,operator_output)



%% all data neede for restart (checkpoint) should be in internal_state
% persistant variable
persistent internal_state
% initialisation on main substructure
internal_state.metis    = [];
internal_state.workflow = [];
internal_state.actor    = [];


%% test of inputs
if (nargin == 0) || isempty(debug_level)
    debug_level = 3;
    % later it will be set to 0 in production
end
if debug_level > 0
    diary debug_dataprovider_in_muscle3.txt
    diary on
end
if debug_level > 1
    addpath(fullfile(fileparts(which(mfilename)),'finish'));
end 
if debug_level > 2
    dbstop if all error
end
if debug_level > 3
    keyboard;
end
internal_state.workflow.debug_level = debug_level;
%
internal_state.workflow.checkpoint_on_off = true;
internal_state.workflow.overlay_on_off = true;

if nargin < 2 || isempty(operator_input)
    operator_input = 'F_INIT';
end
if nargin < 3 || isempty(operator_output)
    switch upper(operator_input)
        case 'NONE'
            operator_output = 'O_I';
        otherwise
            operator_output = 'O_F';
    end
end

% intitial comment in diary file
fprintf('Data provider for metis4muscle3 initialization @ %s \n',datestr(now,0))

%% set path and environement
% set environement for METIS if needed
if isempty(which('metis4imas'))
    addpath(fileparts(fileparts(fileparts(fileparts(which(mfilename))))));
    zineb_path;
    addpath(fileparts(which(mfilename)));
    disp('Used path is:')
    path
end
% test python version from module
test_matlab_python_engine_version;
%% declare muscle3 as non-codegen compatible
coder.extrinsic('py.libmuscle.Instance');
% set python environement
setdefault(py.os.environ,'MUSCLE_INSTANCE','macro');



%% actor description
internal_state.actor.settings.list = metis4muscle3_dataprovider_setting_list;
internal_state.actor.io.list       = metis4muscle3_dataprovider_port_list(operator_input,operator_output);
% open information logger in Python environement 
internal_state.actor.logger = open_python_logger_for_muscle3("Starting DATAPROVIDER actor");
% define instance flags
flags = instance_flags;
%% Create libmuscle instance, specify all Ports
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
if ~internal_state.actor.instance.reuse_instance()
    internal_state.actor.logger.info("Unable to create Muscle3 instance")
    error('Communication with MUSCLE3 failed')
else
    internal_state.actor.logger.info("Muscle3 instance created")  
end
%% retrieve data from Muscle3
% workflow state (initialisation)
internal_state.workflow.state = metis4muscle3_state_initialisation;
% get setting
internal_state.actor.settings.value = get_muscle3_settings(internal_state.actor.instance,internal_state.actor.settings.list);
% get connected ports
internal_state.actor.io.status = metis4muscle3_port_status(internal_state.actor.instance,internal_state.actor.io);

% read metis file
internal_state.metis.loaded_data.post = load_metis_imas(internal_state.actor.settings.value.dataprovider.metis.filename);

% create IMAS dataset if required:
if internal_state.actor.settings.value.dataprovider.imas.writedata ~= 0
    % write both core_ and plasma_ IDSes for matter of compatibility
    internal_state.metis.loaded_data.post.z0dinput.option.plasma_ids  = 'both';
    % tune orientation
    internal_state.metis.loaded_data.post.z0dinput.option.sign = ...
                   internal_state.actor.settings.value.dataprovider.metis.sign;
    internal_state.metis.loaded_data.post.z0dinput.option.orientation = ...
                   internal_state.actor.settings.value.dataprovider.metis.orientation;
    % generate reference for IMAS
    occurrence = [];
    occurrence.user    = internal_state.actor.settings.value.dataprovider.imas.user;
    occurrence.dataversion = sprintf('%d',internal_state.actor.settings.value.dataprovider.imas.version);
    occurrence.tokamak =  internal_state.actor.settings.value.dataprovider.imas.machine;
    occurrence.occurrence = [];
    occurrence.backend = internal_state.actor.settings.value.dataprovider.imas.backend;
    error_flag = metis4imas(internal_state.actor.settings.value.dataprovider.imas.pulse, ...
                            internal_state.actor.settings.value.dataprovider.imas.run, ...
                            occurrence, ...
                            internal_state.metis.loaded_data.post);
end

% open IMAS database
internal_state.imas.idx = imas_open_env_backend(internal_state.actor.settings.value.dataprovider.imas.pulse, ...
                                                internal_state.actor.settings.value.dataprovider.imas.run, ...
                                                internal_state.actor.settings.value.dataprovider.imas.user, ...
                                                internal_state.actor.settings.value.dataprovider.imas.machine, ...
                                                sprintf('%d',internal_state.actor.settings.value.dataprovider.imas.version), ...
                                                internal_state.actor.settings.value.dataprovider.imas.backend);

% here we have a simple loop on time
% find initial time 
t_start = max(min(internal_state.metis.loaded_data.post.profil0d.temps), ...
              internal_state.actor.settings.value.dataprovider.starttime);
index_start = find(internal_state.metis.loaded_data.post.profil0d.temps  >=  t_start,1);  


% get instance state for resuming and other
if internal_state.workflow.checkpoint_on_off 
    internal_state.actor.logger.info("test if resuming mode:")  
    internal_state.workflow.state.resuming       = internal_state.actor.instance.resuming();
    internal_state.workflow.state.should_init    = internal_state.actor.instance.should_init();
end


for index = index_start:length(internal_state.metis.loaded_data.post.profil0d.temps)

    %% default workflow data and state
    internal_state.workflow.clock.dt     =  mean(diff(internal_state.metis.loaded_data.post.profil0d.temps));
    internal_state.workflow.clock.t_cur  =  internal_state.metis.loaded_data.post.profil0d.temps(index);
    internal_state.workflow.clock.t_next =  internal_state.workflow.clock.t_cur + internal_state.workflow.clock.dt;
    
    % end for METIS actor
    if index > fix(length(internal_state.metis.loaded_data.post.profil0d.temps)/2)
        internal_state.workflow.clock.t_next = py.None;
    end
        
    % get imas data for this time
    internal_state.imas.ids_data = get_ids_time_slice(internal_state);

    
    % send message to connected ports
    send_muscle3_messages(internal_state);
    
    
    
    switch upper(operator_input)
        case 'NONE'
            % just sleep and synchronise to prevent problem on overloaded machine
            drawnow;
            py.os.sync();
            pause(1);
        otherwise
            
            % test for checkpointing
            % update checkpoint commands
            if internal_state.workflow.checkpoint_on_off && ~internal_state.workflow.state.final_checkpoint_done
                internal_state.workflow.state.checkpoint   = internal_state.actor.instance.should_save_snapshot(internal_state.workflow.clock.t_cur);
                if ~internal_state.workflow.state.checkpoint
                    internal_state.workflow.state.final_checkpoint = internal_state.actor.instance.should_save_final_snapshot();
                end
            end
            % send checkpoint message if requested
            if internal_state.workflow.state.checkpoint && ~internal_state.workflow.state.final_checkpoint_done
                internal_state.actor.instance.save_snapshot(make_checkpoint_message(internal_state));
                internal_state.workflow.state.checkpoint = false;
                internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
            end
            % it look like this command should be here !!!
            % test if final checkpoint is requested
            if internal_state.workflow.state.final_checkpoint && ~internal_state.workflow.state.final_checkpoint_done
                % finalisation by sending last message
                internal_state.actor.instance.save_final_snapshot(make_checkpoint_message(internal_state));
                internal_state.workflow.state.final_checkpoint = false;
                internal_state.workflow.state.final_checkpoint_done  = true;
            end


            
            
            % end of workflow by default if not another condition has been set
            % check for resuming in reuse loop
            if internal_state.workflow.checkpoint_on_off
                internal_state.actor.logger.info("test if resuming mode:")
                internal_state.workflow.state.resuming       = internal_state.actor.instance.resuming();
                internal_state.workflow.state.should_init    = internal_state.actor.instance.should_init();
            end
            %
            % test if actor is alive
            if ~internal_state.actor.instance.reuse_instance()
                break
            end
            % read message(s)
            internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,operator_input);
            % display message status
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
    
    switch class(internal_state.workflow.clock.t_next)
        case 'py.NoneType'
            break
    end
end  

% one per reuse loop
switch upper(operator_input)
    case 'NONE'
        % test for checkpointing
        % update checkpoint commands
        if internal_state.workflow.checkpoint_on_off && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.workflow.state.checkpoint   = internal_state.actor.instance.should_save_snapshot(internal_state.workflow.clock.t_cur);
            if ~internal_state.workflow.state.checkpoint
                internal_state.workflow.state.final_checkpoint = internal_state.actor.instance.should_save_final_snapshot();
            end
        end
        % send checkpoint message if requested
        if internal_state.workflow.state.checkpoint && ~internal_state.workflow.state.final_checkpoint_done
            internal_state.actor.instance.save_snapshot(make_checkpoint_message(internal_state));
            internal_state.workflow.state.checkpoint = false;
            internal_state.workflow.state.final_checkpoint     = internal_state.actor.instance.should_save_final_snapshot();
        end
        % it look like this command should be here !!!
        % test if final checkpoint is requested
        if internal_state.workflow.state.final_checkpoint && ~internal_state.workflow.state.final_checkpoint_done
            % finalisation by sending last message
            internal_state.actor.instance.save_final_snapshot(make_checkpoint_message(internal_state));
            internal_state.workflow.state.final_checkpoint = false;
            internal_state.workflow.state.final_checkpoint_done  = true;
        end
        
end

% close IMAS database
imas_close(internal_state.imas.idx);

% wait for end of message sending and reading by other actors
switch upper(operator_input)
    case 'NONE'
        internal_state.actor.instance.reuse_instance();
end


% ending
% the exit from matlab
exit;


 

