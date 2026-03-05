function [data,LCFS_nbpts]  = metis4muscle3_update_external_data(internal_state)

% update external data structure
[data,~,LCFS_nbpts] = metis4muscle3_make_external_data(internal_state);
% set external data structure
metis4muscle3_set_external_data(data, ...
                                internal_state.actor.settings.value.metis_external_data, ...
                                internal_state.external_data.state, ...
                                internal_state.workflow.sawtooth);

