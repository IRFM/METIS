function port_list = metis4muscle3_dataprovider_port_list(operator_input,operator_output)  

% port operator:
% F_INIT = intialisation msg
% O_F    = finalisation
% O_I    = output -> input (input of the actor during loop)
% S      = new state (updated data at the end of the loop)

% ports is only used if they are connected !
%
port_list.status.operator       = 'F_INIT';
port_list.status.direction      = 'input';
port_list.status.list_data_name = {'time'};
port_list.status.list_type      = {'int'}; % can be also "ids" for IDSes
port_list.status.comment        = 'synchornisation port with metis4muscle3';
%

%
%  provide data to metis4muscle3 actor
%
%
port_list.clock.operator       = 'O_I'; % default operator, can be change by the workflow
port_list.clock.direction      = 'output';
port_list.clock.list_data_name = {'time'};
port_list.clock.list_type      = {'float'}; % can be also "ids" for IDSes
port_list.clock.comment        = 'external clock time for the METIS actor if time is not provided on anoyher port or internally';

% will be not used in IMAS workflow (should be implemented later)
% port_list.references.operator       = 'O_I';
% port_list.references.direction      = 'output';
% port_list.references.list_data_name = {'cons1t','geo1t','exp0d','option','sepa'};
% port_list.references.list_type      = {'matlab_serialized_structure','matlab_serialized_structure','matlab_serialized_structure','matlab_serialized_structure','matlab_serialized_structure'}; % can be also "ids" for IDSes
% port_list.references.comment        = 'metis zerod evoltution direct control ports (used if the pulse_schedule IDS is not used). If connected overwrite data coming from pulse_schedule';
%
% IDS input ports
%
port_list.pulse_schedule.operator       = 'O_I';
port_list.pulse_schedule.direction      = 'output';
port_list.pulse_schedule.list_data_name = {'pulse_schedule'};
port_list.pulse_schedule.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.pulse_schedule.comment        = 'pulse_schedule IDS for simulation control';
%
port_list.equilibrium.operator       = 'O_I';
port_list.equilibrium.direction      = 'output';
port_list.equilibrium.list_data_name = {'equilibrium'};
port_list.equilibrium.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.equilibrium.comment        = 'equilibirum IDS used as an external data in METIS';
%
port_list.plasma_profiles.operator       = 'O_I';
port_list.plasma_profiles.direction      = 'output';
port_list.plasma_profiles.list_data_name = {'plasma_profiles'};
port_list.plasma_profiles.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_profiles.comment        = 'plasma_profiles IDS used for externally provided kinetic profile in METIS';
%
port_list.plasma_sources.operator       = 'O_I';
port_list.plasma_sources.direction      = 'output';
port_list.plasma_sources.list_data_name = {'plasma_sources'};
port_list.plasma_sources.list_type      = {'ids'}; % can be also "ids" for IDSes
port_list.plasma_sources.comment        = 'plasma_sources IDS used for externally provided source term in METIS';


% update operator if needed
port_list = change_port_operator(port_list,operator_input,operator_output);
   

