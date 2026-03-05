function state = get_workflow_mode(state,actor,used_ports)

for k=length(used_ports):-1:1
    if actor.io.status.(sprintf('pulse_schedule_%s',used_ports{k})) && ...
       ~isempty(actor.io.data.(sprintf('pulse_schedule_%s',used_ports{k})).timestamp)   
        state.time_source = sprintf('pulse_schedule_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from pulse_schedule port")
        return
    elseif actor.io.status.(sprintf('equilibrium_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('equilibrium_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('equilibrium_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from equilibrium port")
        return
    elseif actor.io.status.(sprintf('plasma_profiles_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('plasma_profiles_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('plasma_profiles_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from plasma_profiles port")
        return
    elseif actor.io.status.(sprintf('core_profiles_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('core_profiles_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('core_profiles_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from core_profiles port")
        return
    elseif actor.io.status.(sprintf('plasma_sources_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('plasma_sources_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('plasma_sources_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from plasma_sources port")
        return
    elseif actor.io.status.(sprintf('core_sources_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('core_sources_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('core_sources_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from core_sources port")
        return
    elseif actor.io.status.(sprintf('clock_%s',used_ports{k})) && ...
           ~isempty(actor.io.data.(sprintf('clock_%s',used_ports{k})).timestamp)
        state.time_source = sprintf('clock_%s',used_ports{k});
        state.one_time   = true;
        actor.logger.info("Muscle3 METIS actor is running in time loop with time stamping read from clock port")
        return
    elseif (k==1)
        state.time_source   = 'none';
        state.metis_replay  = true;
        state.metis_init    = true;
        actor.logger.info("Muscle3 METIS actor is running in replay mode")
        return
    end
end
   
