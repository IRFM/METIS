function [data,LCFS_nbpts] = metis4muscle3_init_external_data(internal_state)

% clear existing external data
metis4muscle3_clear_external_data;
% generate new external data
internal_state.external_data.data = [];
% update external data structure
[data,~,LCFS_nbpts] = metis4muscle3_make_external_data(internal_state);
% set external data structure
metis4muscle3_set_external_data(data, ...
                                internal_state.actor.settings.value.metis_external_data, ...
                                internal_state.external_data.state, ...
                                internal_state.workflow.sawtooth);

