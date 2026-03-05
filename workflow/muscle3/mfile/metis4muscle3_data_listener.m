% METIS muscle bridge initialization script.
% R. Nouailletas 2024, CEA
% J-F. Artaud 2024, CEA
% based on work of T. Ravensbergen et M. Schneider 2023

% This script initializes MUSCLE3 (logging, ports),
% initializes the parameters, and then launches a simulation.
% The bridge receives a equilibrium IDS via MUSCLE3, and send a coreprofiles IDS
%function metis4muscle3_data_listener(debug_level,checkpoint_on_off,overlay_on_off,operator_input)
function metis4muscle3_data_listener(debug_level,operator_input)



% all data neede for restart (checkpoint) should be in internal_state
% persistant variable
persistent internal_state
% initialisation on main substructure
internal_state.workflow = [];
internal_state.actor    = [];


% test of inputs
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
%if nargin < 2 || isempty(checkpoint_on_off)
    internal_state.workflow.checkpoint_on_off = false;
%else
%    internal_state.workflow.checkpoint_on_off = checkpoint_on_off;
%end   
%if nargin < 3 || isempty(overlay_on_off)
    internal_state.workflow.overlay_on_off = true;
%else
%    internal_state.workflow.overlay_on_off = overlay_on_off;
%end   
if nargin < 2 || isempty(operator_input)
    operator_input = 'F_INIT';
end

% intitial comment in diary file
fprintf('Data listener for metis4muscle3 initialization @ %s \n',datestr(now,0))

% set path and environement
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
% declare muscle3 as non-codegen compatible
coder.extrinsic('py.libmuscle.Instance');
% set python environement
setdefault(py.os.environ,'MUSCLE_INSTANCE','macro');



% actor description
internal_state.actor.settings.list = metis4muscle3_listener_setting_list;
internal_state.actor.io.list       = metis4muscle3_listener_port_list(operator_input,'');
% open information logger in Python environement 
internal_state.actor.logger = open_python_logger_for_muscle3("Starting LISTENER actor");
% defins instance flags
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
%     if internal_state.workflow.checkpoint_on_off ||  internal_state.workflow.overlay_on_off
%         sum_flags = double(internal_state.workflow.checkpoint_on_off .* flags.USES_CHECKPOINT_API) + ...
%                     double(~internal_state.workflow.overlay_on_off .* flags.DONT_APPLY_OVERLAY) +  ...
%                     double(true .* flags.KEEPS_NO_STATE_FOR_NEXT_USE);
%         internal_state.actor.instance = py.libmuscle.Instance(internal_state.actor.io.ports,sum_flags);
%     else
        internal_state.actor.instance = py.libmuscle.Instance(internal_state.actor.io.ports,flags.KEEPS_NO_STATE_FOR_NEXT_USE);
%    end
end
if ~internal_state.actor.instance.reuse_instance()
    internal_state.actor.logger.info("Unable to create Muscle3 instance")
    error('Communication with MUSCLE3 failed')
elseif internal_state.workflow.checkpoint_on_off
    % not use by METIS at this point but should be read as asked
    % by MUSCLE3
    internal_state.actor.instance.resuming();
    internal_state.actor.instance.should_init();
else
    internal_state.actor.logger.info("Muscle3 instance created")
end

% retrieve data from Muscle3
% workflow state (initialisation)
internal_state.workflow.state = metis4muscle3_state_initialisation;
% get setting
internal_state.actor.settings.value = get_muscle3_settings(internal_state.actor.instance,internal_state.actor.settings.list);
% get connected ports
internal_state.actor.io.status = metis4muscle3_port_status(internal_state.actor.instance,internal_state.actor.io);

% here we have a simple loop on time
for k=1:internal_state.actor.settings.value.listener.maxloop
    % wait for status from METIS
    internal_state.actor.io.data = metis4muscle3_read(internal_state.actor,'input',internal_state.workflow.overlay_on_off);
    % display message contents
    noms = fieldnames(internal_state.actor.io.data);
    for l=1:length(noms)
        if ~isempty(internal_state.actor.io.data.(noms{l}).data_raw)
            disp('--------------------------------------');
            fprintf('Message content for %s:\n',noms{l});
            internal_state.actor.io.data.(noms{l})
            disp('with data structure:')
            internal_state.actor.io.data.(noms{l}).data
            switch class(internal_state.actor.io.data.(noms{l}).next_timestamp)
                case 'py.NoneType'
                    % End of time message
                    internal_state.workflow.clock.t_next = Inf;
                    %dt     = internal_state.workflow.clock.dt;
                    internal_state.workflow.clock.t_cur  =  internal_state.actor.io.data.(noms{l}).timestamp;                    
                    disp('Received from MUSCLE 3 request to end the simulation in next_timestamp: py.NoneType');
                otherwise                    
                    internal_state.workflow.clock.dt     =  internal_state.actor.io.data.(noms{l}).next_timestamp - ...
                        internal_state.actor.io.data.(noms{l}).timestamp;
                    internal_state.workflow.clock.t_cur  =  internal_state.actor.io.data.(noms{l}).timestamp;
                    internal_state.workflow.clock.t_next =  internal_state.actor.io.data.(noms{l}).next_timestamp;
            end
        end
    end
    
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
    % end of workflow by dfault if not another condition has been set
    % test if actor is alive
    if ~internal_state.actor.instance.reuse_instance()
        break
        % check for resuming in reuse loop
    elseif internal_state.workflow.checkpoint_on_off
        % not use by METIS at this point but should be read as asked
        % by MUSCLE3
        internal_state.actor.instance.resuming();
        internal_state.actor.instance.should_init();
    end
    
end          

% ending
% update checkpoint commands
if internal_state.workflow.checkpoint_on_off
    internal_state.workflow.state.final_checkpoint = internal_state.actor.instance.should_save_final_snapshot();
end
% test if final checkpoint is requested
if internal_state.workflow.state.final_checkpoint
    % finalisation by sending last message
    internal_state.actor.instance.save_final_snapshot(make_checkpoint_message(internal_state));
end

% the exit from matlab
exit;


 

