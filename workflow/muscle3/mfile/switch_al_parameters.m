function option = switch_al_parameters(option,port_status,port_name_legacy)

% switch to DDv 4.x
option.plasma_ids = 'both';

if port_name_legacy
    % port name before standardization
    
    %dependaning on connected ports
    if port_status.equilibrium_output_f || port_status.equilibrium_output_i
        option.equilibrium = 1;
    else
        option.equilibrium = 0;
    end
    if port_status.pulse_schedule_output_f || port_status.pulse_schedule_output_i
        option.pulse_schedule = 1;
    else
        option.pulse_schedule = 0;
    end
    if port_status.plasma_profiles_output_f || port_status.plasma_profiles_output_i || ...
            port_status.core_profiles_output_f || port_status.core_profiles_output_i
        option.plasma_profiles = 1;
    else
        option.plasma_profiles = 0;
    end
    if port_status.plasma_sources_output_f || port_status.plasma_sources_output_i || ...
            port_status.core_sources_output_f || port_status.core_sources_output_i
        option.plasma_sources = 1;
    else
        option.plasma_sources = 0;
    end
    if get_imas_version >= 4
        option.dataset_description = 0;
    elseif port_status.dataset_description_output_f || port_status.dataset_description_output_i
        option.dataset_description = 1;
    else
        option.dataset_description = 0;
    end
    if port_status.dataset_fair_output_f || port_status.dataset_fair_output_i
        option.dataset_fair = 1;
    else
        option.dataset_fair = 0;
    end
    if port_status.summary_output_f || port_status.summary_output_i
        option.summary = 1;
    else
        option.summary = 0;
    end
    if port_status.plasma_transport_output_f || port_status.plasma_transport_output_i || ...
            port_status.core_transport_output_f || port_status.core_transport_output_i
        option.plasma_transport = 1;
    else
        option.plasma_transport = 0;
    end
    if port_status.radiation_output_f || port_status.radiation_output_i
        option.radiation = 1;
    else
        option.radiation = 0;
    end
    if port_status.edge_output_f || port_status.edge_output_i
        option.edge = 1;
    else
        option.edge = 0;
    end
    if port_status.transport_solver_numerics_output_f || port_status.transport_solver_numerics_output_i
        option.numerics = 1;
    else
        option.numerics = 0;
    end
    
else
    % port name after standardization
    %dependaning on connected ports
    if port_status.equilibrium_out_f || port_status.equilibrium_out_i
        option.equilibrium = 1;
    else
        option.equilibrium = 0;
    end
    if port_status.pulse_schedule_out_f || port_status.pulse_schedule_out_i
        option.pulse_schedule = 1;
    else
        option.pulse_schedule = 0;
    end
    if port_status.plasma_profiles_out_f || port_status.plasma_profiles_out_i || ...
            port_status.core_profiles_out_f || port_status.core_profiles_out_i
        option.plasma_profiles = 1;
    else
        option.plasma_profiles = 0;
    end
    if port_status.plasma_sources_out_f || port_status.plasma_sources_out_i || ...
            port_status.core_sources_out_f || port_status.core_sources_out_i
        option.plasma_sources = 1;
    else
        option.plasma_sources = 0;
    end
    if get_imas_version >= 4 
        option.dataset_description = 0;
    elseif port_status.dataset_description_out_f || port_status.dataset_description_out_i
        option.dataset_description = 1;
    else
        option.dataset_description = 0;
    end
    if port_status.dataset_fair_out_f || port_status.dataset_fair_out_i
        option.dataset_fair = 1;
    else
        option.dataset_fair = 0;
    end
    if port_status.summary_out_f || port_status.summary_out_i
        option.summary = 1;
    else
        option.summary = 0;
    end
    if port_status.plasma_transport_out_f || port_status.plasma_transport_out_i || ...
            port_status.core_transport_out_f || port_status.core_transport_out_i
        option.plasma_transport = 1;
    else
        option.plasma_transport = 0;
    end
    if port_status.radiation_out_f || port_status.radiation_out_i
        option.radiation = 1;
    else
        option.radiation = 0;
    end
    if port_status.edge_out_f || port_status.edge_out_i
        option.edge = 1;
    else
        option.edge = 0;
    end
    if port_status.transport_solver_numerics_out_f || port_status.transport_solver_numerics_out_i
        option.numerics = 1;
    else
        option.numerics = 0;
    end
end


