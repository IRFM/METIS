function on_off = get_st_state(io_status, settings_metis,external_data_state,used_ports)

% default
on_off = 'off';

switch settings_metis.sawtooth.onoff
    case 'off'
        return
end
if (external_data_state.equilibrium == 1) &&  ...
   (external_data_state.current_diffusion ==1)
    
    for k=1:length(used_ports)
        if io_status.(sprintf('equilibrium_%s',used_ports{k})) == 1
             on_off = 'on';
             break;
        end
    end
end

