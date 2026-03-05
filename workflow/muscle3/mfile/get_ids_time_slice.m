function data = get_ids_time_slice(internal_state)

% initialisation
data =[];

% loop on connected port
noms = fieldnames(internal_state.actor.io.status);
for k=1:length(noms)
    if internal_state.actor.io.status.(noms{k})
        % only ouputed data
        switch internal_state.actor.io.list.(noms{k}).direction
            case 'output'
                list = internal_state.actor.io.list.(noms{k}).list_type;
                for l=1:length(list)
                    switch list{l}
                        case 'ids'
                            data.(noms{k}){l} = ids_get_slice(internal_state.imas.idx, ...
                                  internal_state.actor.io.list.(noms{k}).list_data_name{l}, ...
                                  internal_state.workflow.clock.t_cur,2);
                    end
                      
                end
        end
    end
end