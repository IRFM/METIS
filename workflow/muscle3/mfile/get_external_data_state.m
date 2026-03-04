function state = get_external_data_state(ports,settings)

% initialisation
state = [];
% loop on settings state
noms = fieldnames(settings.metis_external_data);
nomp = fieldnames(ports.list);
for k=1:length(noms)
    state.(noms{k}) = 0;
    if settings.metis_external_data.(noms{k}).state
        for l=1:length(nomp)
            switch ports.list.(nomp{l}).direction
                case {'input','init'}
                    if  strcmp(settings.metis_external_data.(noms{k}).ids_name, ...
                            ports.list.(nomp{l}).list_data_name)
                        if ports.status.(nomp{l})
                            state.(noms{k}) = 1;
                        end
                    end
                    if (state.(noms{k}) == 0) && startsWith(settings.metis_external_data.(noms{k}).ids_name,'plasma_')
                        if  strcmp(strrep(settings.metis_external_data.(noms{k}).ids_name,'plasma_','core_'), ...
                                ports.list.(nomp{l}).list_data_name)
                            if ports.status.(nomp{l})
                                state.(noms{k}) = 1;
                            end
                        end
                        
                    end
            end
        end
    end
end