function [t_cur, t_next, dt] =  metis4muscle3_time_operator(internal_state)

switch internal_state.actor.settings.value.metis.clock.internal
    case 'on'
        disp('Using METIS actor internal clock');
        if ~isfinite(internal_state.workflow.clock.dt) || ~isfinite(internal_state.workflow.clock.t_next)
            disp('Initialisation of METIS actor internal clock')
            dt     = internal_state.actor.settings.value.metis.clock.time.dt;
            t_cur  = internal_state.actor.settings.value.metis.clock.time.start;
            t_next = t_cur + dt;
        else
            dt     = internal_state.actor.settings.value.metis.clock.time.dt;
            t_cur  = internal_state.workflow.clock.t_next;
            t_next = t_cur + dt;
        end
    otherwise
        fprintf('Using MUSCLE3 message datation based on port %s\n',internal_state.workflow.state.time_source);
        % get time from message
        t_cur  = internal_state.actor.io.data.(internal_state.workflow.state.time_source).timestamp;
        switch class(internal_state.actor.io.data.(internal_state.workflow.state.time_source).next_timestamp)
            case 'py.NoneType'
                % End of time message
                t_next = Inf;
                dt     = internal_state.actor.settings.value.metis.clock.time.dt;
                disp('Received from MUSCLE 3 request to end the simulation in next_timestamp: py.NoneType');
            otherwise             
                t_next = internal_state.actor.io.data.(internal_state.workflow.state.time_source).next_timestamp;
                dt     = t_next - t_cur;
        end
        
        % change synchornisation
        if internal_state.actor.settings.value.metis.usenexttimestamp && ...
           isfinite(t_next)
            disp('METIS will compute for next_timestamp');
            t_cur  = t_next;
            t_next = t_cur + dt;
        else
            disp('METIS will compute fo timestamp');
        end
        
end