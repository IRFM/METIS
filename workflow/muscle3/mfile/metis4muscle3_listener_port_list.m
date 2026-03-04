function port_list = metis4muscle3_listener_port_list(operator_input,operator_output) 

% port operator:
% F_INIT = intialisation msg
% O_F    = finalisation
% O_I    = output -> input (input of the actor during internal convergence  loop)
% S      = new state (updated data at the end of the  internal convergence loop)

% ports is only used if they are connected !
port_list.clock_output.operator       = 'F_INIT';   % default operator can be change in workflow
port_list.clock_output.direction      = 'input';
port_list.clock_output.list_data_name = {'time'};
port_list.clock_output.list_type      = {'empty'}; % can be also "ids" for IDSes
port_list.clock_output.comment        = 'external clock time for the METIS actor if time is not provided on anoyher port or internally (can be any port as only time_stamp data is used)';
%
% IDS input ports
port_list.pulse_schedule_output.operator       = 'S';   % default operator, can be change by the workflow
port_list.pulse_schedule_output.direction      = 'input';
port_list.pulse_schedule_output.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule_output.comment        = 'pulse_schedule IDS computed in METIS';
%
port_list.equilibrium_output.operator       = 'S';
port_list.equilibrium_output.direction      = 'input';
port_list.equilibrium_output.list_data_name = {'equilibrium'};
port_list.equilibrium_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium_output.comment        = 'equilibirum IDS computed by METIS';
%
port_list.plasma_profiles_output.operator       = 'S';
port_list.plasma_profiles_output.direction      = 'input';
port_list.plasma_profiles_output.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles_output.comment        = 'core_profiles IDS computed by METIS';
%
port_list.plasma_sources_output.operator       = 'S';
port_list.plasma_sources_output.direction      = 'input';
port_list.plasma_sources_output.list_data_name = {'plasma_sources'};
port_list.plasma_sources_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources_output.comment        = 'core_sources IDS computed by METIS';
%
% other possible IDSs
%
if get_imas_version < 4
    port_list.dataset_description_output.operator       = 'S';
    port_list.dataset_description_output.direction      = 'input';
    port_list.dataset_description_output.list_data_name = {'dataset_description'};
    port_list.dataset_description_output.list_type      = {'ids'}; % can be also "ids" for IDSes
    port_list.dataset_description_output.comment        = 'dataset_description IDS computed by METIS';
end
%
port_list.dataset_fair_output.operator       = 'S';
port_list.dataset_fair_output.direction      = 'input';
port_list.dataset_fair_output.list_data_name = {'dataset_fair'};
port_list.dataset_fair_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.dataset_fair_output.comment        = 'dataset_fair IDS computed by METIS';
%
port_list.summary_output.operator       = 'S';
port_list.summary_output.direction      = 'input';
port_list.summary_output.list_data_name = {'summary'};
port_list.summary_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.summary_output.comment        = 'summary IDS computed by METIS';
%
port_list.plasma_transport_output.operator       = 'S';
port_list.plasma_transport_output.direction      = 'input';
port_list.plasma_transport_output.list_data_name = {'plasma_transport'};
port_list.plasma_transport_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_transport_output.comment        = 'core_transport IDS computed by METIS';
%
port_list.radiation_output.operator       = 'S';
port_list.radiation_output.direction      = 'input';
port_list.radiation_output.list_data_name = {'radiation'};
port_list.radiation_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.radiation_output.comment        = 'radiation IDS computed by METIS';
%
port_list.edge_output.operator       = 'S';
port_list.edge_output.direction      = 'input';
port_list.edge_output.list_data_name = {'edge_profiles','edge_transport'};
port_list.edge_output.list_type      = {'ids','ids'}; % can be also "ids" for IDSes
port_list.edge_output.comment        = 'radiation IDS computed by METIS';
%
port_list.transport_solver_numerics_output.operator       = 'S';
port_list.transport_solver_numerics_output.direction      = 'input';
port_list.transport_solver_numerics_output.list_data_name = {'transport_solver_numerics'};
port_list.transport_solver_numerics_output.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.transport_solver_numerics_output.comment        = 'transport_solver_numerics IDS computed by METIS';

% update operator if needed
port_list = change_port_operator(port_list,operator_input,operator_output);
   

    

