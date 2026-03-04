% new version with standardized port names: 
% Regarding port names, after iterations with Artaud Jean Francois EXT , Coster David EXT  and discussions at the TSVV11 workshop with Fox Peter EXT  and Casson Francis EXT , we ended up with the following convention:
% 	Port names will be fully lower case
% 	They will contain the name of the IDS, the direction (input/output), and optionally the operator (F_init, O_i, S, O_f) and occurrence number, i.e.:
% 
% <ids name>_<direction>[_operator][_occurrence] 
% 	
% 1/	ids_name (mandatory) = any IDS of the Data Dictionary
% 2/	_direction (mandatory) and _operator (optional): 
% 		_f_init (initialisation)   --> No operator
% 		_in (input)                --> Operator _f for F_init and _s for S
% 		_out (output)              --> Operator _i for O_i and _f for O_f
% 3/	_occurrence (optional): none when it is 0 (like for IDSs), optional as last position when non-zero: {{_x}}{ { (with x = integer).}}
%
% Such convention is compatible with the name of ports as already defined with iWrap (because the operator and occurrence are optional), 
% which means nothing needs to be changed in the ETS, it fits the need of complex workflows like the METIS one, and the anticipated needs for using occurrences.
% Examples:
% 
% metis.equilibrium_f_init
% metis.equilibrium_in_f
% metis.equilibrium_in_s 
% metis.equilibrium_out_i
% metis.equilibrium_out_f 
% 
% ets.interferometer_in
% ets.interferometer_in_1
% ets.core_profiles_out
% ets.core_profiles_out_1
% 
% torbeam.equilibrium_in
% torbeam.core_profiles_in
% torbeam.ec_launchers_in
% torbeam.waves_out
% etc.
%  
 


function port_list = metis4muscle3_port_list
%function port_list = metis4muscle3_port_list(operator_init,operator_input,operator_output) 

% port operator:
% F_INIT = intialisation msg
% O_F    = finalisation
% O_I    = output -> input (input of the actor during internal convergence  loop)
% S      = new state (updated data at the end of the  internal convergence loop)

% ports is only used if they are connected !
% INIT port are read by default only during actor initialisation
%
% clock on initialisation
%
port_list.clock_f_init.operator       = 'F_INIT';   % default operator can be change in workflow
port_list.clock_f_init.direction      = 'init';
port_list.clock_f_init.list_data_name = {'time'};
port_list.clock_f_init.list_type      = {'empty'}; % can be also "ids" for IDSes
port_list.clock_f_init.comment        = 'external clock time for the METIS actor if time is not provided on anoyher port or internally (can be any port as only time_stamp data is used)';
%
% IDS input ports on initialisation
%
port_list.pulse_schedule_f_init.operator       = 'F_INIT';
port_list.pulse_schedule_f_init.direction      = 'init';
port_list.pulse_schedule_f_init.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_f_init.comment        = 'pulse_schedule IDS for simulation control';
%
port_list.equilibrium_f_init.operator       = 'F_INIT';
port_list.equilibrium_f_init.direction      = 'init';
port_list.equilibrium_f_init.list_data_name = {'equilibrium'};
port_list.equilibrium_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_f_init.comment        = 'equilibirum IDS used as an external data in METIS';
%
port_list.plasma_profiles_f_init.operator       = 'F_INIT';
port_list.plasma_profiles_f_init.direction      = 'init';
port_list.plasma_profiles_f_init.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_f_init.comment        = 'plasma_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.plasma_sources_f_init.operator       = 'F_INIT';
port_list.plasma_sources_f_init.direction      = 'init';
port_list.plasma_sources_f_init.list_data_name = {'plasma_sources'};
port_list.plasma_sources_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_f_init.comment        = 'plasma_sources IDS used for externally provided source term in METIS';
%
port_list.core_profiles_f_init.operator       = 'F_INIT';
port_list.core_profiles_f_init.direction      = 'init';
port_list.core_profiles_f_init.list_data_name = {'core_profiles'};
port_list.core_profiles_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_profiles_f_init.comment        = 'core_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.core_sources_f_init.operator       = 'F_INIT';
port_list.core_sources_f_init.direction      = 'init';
port_list.core_sources_f_init.list_data_name = {'core_sources'};
port_list.core_sources_f_init.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_sources_f_init.comment        = 'core_sources IDS used for externally provided source term in METIS';

% Input port with operator F_INIT: these ports a read at each time step after
% initialisation. The data from port F_INIT are erase by data from ports S
% if both are connected.
%
% clock update
% This port is used to synchronise METIS by default on another actor is it
% is used. This can be change globaly by input parameter of actor
% metis4muscle3
port_list.clock_in_f.operator       = 'F_INIT';   % default operator can be change in workflow
port_list.clock_in_f.direction      = 'input';
port_list.clock_in_f.list_data_name = {'time'};
port_list.clock_in_f.list_type      = {'empty'}; % can be also "ids" for IDSes
port_list.clock_in_f.comment        = 'external clock time for the METIS actor if time is not provided on anoyher port or internally (can be any port as only time_in_stamp data is used)';
%
% IDS input ports on state update
%
port_list.pulse_schedule_in_f.operator       = 'F_INIT';
port_list.pulse_schedule_in_f.direction      = 'input';
port_list.pulse_schedule_in_f.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_in_f.comment        = 'pulse_schedule IDS for simulation control';
%
port_list.equilibrium_in_f.operator       = 'F_INIT';
port_list.equilibrium_in_f.direction      = 'input';
port_list.equilibrium_in_f.list_data_name = {'equilibrium'};
port_list.equilibrium_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_in_f.comment        = 'equilibirum IDS used as an external data in METIS';
%
port_list.plasma_profiles_in_f.operator       = 'F_INIT';
port_list.plasma_profiles_in_f.direction      = 'input';
port_list.plasma_profiles_in_f.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_in_f.comment        = 'plasma_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.plasma_sources_in_f.operator       = 'F_INIT';
port_list.plasma_sources_in_f.direction      = 'input';
port_list.plasma_sources_in_f.list_data_name = {'plasma_sources'};
port_list.plasma_sources_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_in_f.comment        = 'plasma_sources IDS used for externally provided source term in METIS';
%
port_list.core_profiles_in_f.operator       = 'F_INIT';
port_list.core_profiles_in_f.direction      = 'input';
port_list.core_profiles_in_f.list_data_name = {'core_profiles'};
port_list.core_profiles_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_profiles_in_f.comment        = 'core_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.core_sources_in_f.operator       = 'F_INIT';
port_list.core_sources_in_f.direction      = 'input';
port_list.core_sources_in_f.list_data_name = {'core_sources'};
port_list.core_sources_in_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_sources_in_f.comment        = 'core_sources IDS used for externally provided source term in METIS';


% Input port with operator S: these ports a read at each time step after
% initialisation. The data from port S overwrite data from  ports F_INIT if
% both are connected
%
% clock update
% This port is used to synchronise METIS by default on another actor is it
% is used. This can be change globaly by input parameter of actor
% metis4muscle3
port_list.clock_in_s.operator       = 'S';   % default operator can be change in workflow
port_list.clock_in_s.direction      = 'input';
port_list.clock_in_s.list_data_name = {'time'};
port_list.clock_in_s.list_type      = {'empty'}; % can be also "ids" for IDSes
port_list.clock_in_s.comment        = 'external clock time for the METIS actor if time is not provided on anoyher port or internally (can be any port as only time_in_stamp data is used)';
%
% IDS input ports on state update
%
port_list.pulse_schedule_in_s.operator       = 'S';
port_list.pulse_schedule_in_s.direction      = 'input';
port_list.pulse_schedule_in_s.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_in_s.comment        = 'pulse_schedule IDS for simulation control';
%
port_list.equilibrium_in_s.operator       = 'S';
port_list.equilibrium_in_s.direction      = 'input';
port_list.equilibrium_in_s.list_data_name = {'equilibrium'};
port_list.equilibrium_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_in_s.comment        = 'equilibirum IDS used as an external data in METIS';
%
port_list.plasma_profiles_in_s.operator       = 'S';
port_list.plasma_profiles_in_s.direction      = 'input';
port_list.plasma_profiles_in_s.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_in_s.comment        = 'plasma_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.plasma_sources_in_s.operator       = 'S';
port_list.plasma_sources_in_s.direction      = 'input';
port_list.plasma_sources_in_s.list_data_name = {'plasma_sources'};
port_list.plasma_sources_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_in_s.comment        = 'plasma_sources IDS used for externally provided source term in METIS';
%
port_list.core_profiles_in_s.operator       = 'S';
port_list.core_profiles_in_s.direction      = 'input';
port_list.core_profiles_in_s.list_data_name = {'core_profiles'};
port_list.core_profiles_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_profiles_in_s.comment        = 'core_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.core_sources_in_s.operator       = 'S';
port_list.core_sources_in_s.direction      = 'input';
port_list.core_sources_in_s.list_data_name = {'core_sources'};
port_list.core_sources_in_s.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_sources_in_s.comment        = 'core_sources IDS used for externally provided source term in METIS';

%
% IDS output ports (O_F)
%
port_list.pulse_schedule_out_f.operator       = 'O_F';   % default operator, can be change by the workflow
port_list.pulse_schedule_out_f.direction      = 'output';
port_list.pulse_schedule_out_f.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_out_f.comment        = 'pulse_in_schedule IDS computed in METIS';
%
port_list.equilibrium_out_f.operator       = 'O_F';
port_list.equilibrium_out_f.direction      = 'output';
port_list.equilibrium_out_f.list_data_name = {'equilibrium'};
port_list.equilibrium_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_out_f.comment        = 'equilibirum IDS computed by METIS';
%
port_list.plasma_profiles_out_f.operator       = 'O_F';
port_list.plasma_profiles_out_f.direction      = 'output';
port_list.plasma_profiles_out_f.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_out_f.comment        = 'plasma_profiles IDS computed by METIS';
%
port_list.plasma_sources_out_f.operator       = 'O_F';
port_list.plasma_sources_out_f.direction      = 'output';
port_list.plasma_sources_out_f.list_data_name = {'plasma_sources'};
port_list.plasma_sources_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_out_f.comment        = 'plasma_sources IDS computed by METIS';
%
port_list.core_profiles_out_f.operator       = 'O_F';
port_list.core_profiles_out_f.direction      = 'output';
port_list.core_profiles_out_f.list_data_name = {'core_profiles'};
port_list.core_profiles_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_profiles_out_f.comment        = 'core_profiles IDS computed by METIS';
%
port_list.core_sources_out_f.operator       = 'O_F';
port_list.core_sources_out_f.direction      = 'output';
port_list.core_sources_out_f.list_data_name = {'core_sources'};
port_list.core_sources_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_sources_out_f.comment        = 'core_sources IDS computed by METIS';
%
% other possible IDSs
%
if get_imas_version < 4
    port_list.dataset_description_out_f.operator       = 'O_F';
    port_list.dataset_description_out_f.direction      = 'output';
    port_list.dataset_description_out_f.list_data_name = {'dataset_description'};
    port_list.dataset_description_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
    port_list.dataset_description_out_f.comment        = 'dataset_description IDS computed by METIS';
end
%
port_list.dataset_fair_out_f.operator       = 'O_F';
port_list.dataset_fair_out_f.direction      = 'output';
port_list.dataset_fair_out_f.list_data_name = {'dataset_fair'};
port_list.dataset_fair_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.dataset_fair_out_f.comment        = 'dataset_fair IDS computed by METIS';
%
port_list.summary_out_f.operator       = 'O_F';
port_list.summary_out_f.direction      = 'output';
port_list.summary_out_f.list_data_name = {'summary'};
port_list.summary_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.summary_out_f.comment        = 'summary IDS computed by METIS';
%
port_list.plasma_transport_out_f.operator       = 'O_F';
port_list.plasma_transport_out_f.direction      = 'output';
port_list.plasma_transport_out_f.list_data_name = {'plasma_transport'};
port_list.plasma_transport_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_transport_out_f.comment        = 'plasma_transport IDS computed by METIS';
%
port_list.core_transport_out_f.operator       = 'O_F';
port_list.core_transport_out_f.direction      = 'output';
port_list.core_transport_out_f.list_data_name = {'core_transport'};
port_list.core_transport_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_transport_out_f.comment        = 'core_transport IDS computed by METIS';
%
port_list.radiation_out_f.operator       = 'O_F';
port_list.radiation_out_f.direction      = 'output';
port_list.radiation_out_f.list_data_name = {'radiation'};
port_list.radiation_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.radiation_out_f.comment        = 'radiation IDS computed by METIS';
%
port_list.edge_out_f.operator       = 'O_F';
port_list.edge_out_f.direction      = 'output';
port_list.edge_out_f.list_data_name = {'edge_profiles','edge_transport'};
port_list.edge_out_f.list_type      = {'ids','ids'}; % can be also "ids" for IDSes
port_list.edge_out_f.comment        = 'radiation IDS computed by METIS';
%
port_list.transport_solver_numerics_out_f.operator       = 'O_F';
port_list.transport_solver_numerics_out_f.direction      = 'output';
port_list.transport_solver_numerics_out_f.list_data_name = {'transport_solver_numerics'};
port_list.transport_solver_numerics_out_f.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.transport_solver_numerics_out_f.comment        = 'transport_in_solver_numerics IDS computed by METIS';
%
% IDS output ports (O_I)
%
port_list.pulse_schedule_out_i.operator       = 'O_I';   % default operator, can be change by the workflow
port_list.pulse_schedule_out_i.direction      = 'output';
port_list.pulse_schedule_out_i.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_out_i.comment        = 'pulse_in_schedule IDS computed in METIS';
%
port_list.equilibrium_out_i.operator       = 'O_I';
port_list.equilibrium_out_i.direction      = 'output';
port_list.equilibrium_out_i.list_data_name = {'equilibrium'};
port_list.equilibrium_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_out_i.comment        = 'equilibirum IDS computed by METIS';
%
port_list.plasma_profiles_out_i.operator       = 'O_I';
port_list.plasma_profiles_out_i.direction      = 'output';
port_list.plasma_profiles_out_i.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_out_i.comment        = 'plasma_profiles IDS computed by METIS';
%
port_list.plasma_sources_out_i.operator       = 'O_I';
port_list.plasma_sources_out_i.direction      = 'output';
port_list.plasma_sources_out_i.list_data_name = {'plasma_sources'};
port_list.plasma_sources_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_out_i.comment        = 'plasma_sources IDS computed by METIS';
%
port_list.core_profiles_out_i.operator       = 'O_I';
port_list.core_profiles_out_i.direction      = 'output';
port_list.core_profiles_out_i.list_data_name = {'core_profiles'};
port_list.core_profiles_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_profiles_out_i.comment        = 'plasma_profiles IDS computed by METIS';
%
port_list.core_sources_out_i.operator       = 'O_I';
port_list.core_sources_out_i.direction      = 'output';
port_list.core_sources_out_i.list_data_name = {'core_sources'};
port_list.core_sources_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_sources_out_i.comment        = 'core_sources IDS computed by METIS';
%
% other possible IDSs
%
if get_imas_version < 4    
    port_list.dataset_description_out_i.operator       = 'O_I';
    port_list.dataset_description_out_i.direction      = 'output';
    port_list.dataset_description_out_i.list_data_name = {'dataset_description'};
    port_list.dataset_description_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
    port_list.dataset_description_out_i.comment        = 'dataset_description IDS computed by METIS';
end
%
port_list.dataset_fair_out_i.operator       = 'O_I';
port_list.dataset_fair_out_i.direction      = 'output';
port_list.dataset_fair_out_i.list_data_name = {'dataset_fair'};
port_list.dataset_fair_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.dataset_fair_out_i.comment        = 'dataset_fair IDS computed by METIS';
%
port_list.summary_out_i.operator       = 'O_I';
port_list.summary_out_i.direction      = 'output';
port_list.summary_out_i.list_data_name = {'summary'};
port_list.summary_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.summary_out_i.comment        = 'summary IDS computed by METIS';
%
port_list.plasma_transport_out_i.operator       = 'O_I';
port_list.plasma_transport_out_i.direction      = 'output';
port_list.plasma_transport_out_i.list_data_name = {'plasma_transport'};
port_list.plasma_transport_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_transport_out_i.comment        = 'plasma_transport IDS computed by METIS';
%
port_list.core_transport_out_i.operator       = 'O_I';
port_list.core_transport_out_i.direction      = 'output';
port_list.core_transport_out_i.list_data_name = {'core_transport'};
port_list.core_transport_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.core_transport_out_i.comment        = 'core_transport IDS computed by METIS';
%
port_list.radiation_out_i.operator       = 'O_I';
port_list.radiation_out_i.direction      = 'output';
port_list.radiation_out_i.list_data_name = {'radiation'};
port_list.radiation_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.radiation_out_i.comment        = 'radiation IDS computed by METIS';
%
port_list.edge_out_i.operator       = 'O_I';
port_list.edge_out_i.direction      = 'output';
port_list.edge_out_i.list_data_name = {'edge_profiles','edge_transport'};
port_list.edge_out_i.list_type      = {'ids','ids'}; % can be also "ids" for IDSes
port_list.edge_out_i.comment        = 'radiation IDS computed by METIS';
%
port_list.transport_solver_numerics_out_i.operator       = 'O_I';
port_list.transport_solver_numerics_out_i.direction      = 'output';
port_list.transport_solver_numerics_out_i.list_data_name = {'transport_solver_numerics'};
port_list.transport_solver_numerics_out_i.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.transport_solver_numerics_out_i.comment        = 'transport_in_solver_numerics IDS computed by METIS';


%
% update operator if needed
% no more used as all type of ports is alerady defined
%port_list = change_port_operator(port_list,operator_init,operator_input,operator_output);
   

    

