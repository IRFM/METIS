function info = isany_f_init_used(actor,used_ports)

info = false;
if iscell(used_ports)
    for k=1:length(used_ports)
        switch  used_ports{k}
            case 'in_f'
                info = isany_f_init_used(actor,'in_f');
            case 'input_f'
                info = isany_f_init_used(actor,'input_f');
            case 'f_init'
                info = isany_f_init_used(actor,'f_init');
        end
        if info
            break;
        end
    end
else
    noms = fieldnames(actor.io.status);
    for k=1:length(noms)
        if endsWith(noms{k},used_ports)
           if  actor.io.status.(noms{k})
              info = true;
              break;
           end
        end
    end
end
if actor.settings.value.metis.alwaysreadinitport && ~info
    % check also f_init port
    info = isany_f_init_used(actor,'f_init');
end

