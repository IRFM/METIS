function internal_state_out = metis4muscle3_resume_simulation(internal_state_in,resume_message,update_time)

% inputs handling
if (nargin < 3) || isempty(update_time)
    update_time = false;
end

% deserializee message
internal_state_out = get_matlab_serialized_structure(resume_message.data);

% update clock
if update_time
    internal_state_out.workflow.clock.t_cur  = resume_message.timestamp;
    internal_state_out.workflow.clock.t_next = resume_message.next_timestamp;
    internal_state_out.workflow.clock.dt     = resume_message.next_timestamp - resume_message.timestamp; 
end

% change pointer for instance and logger
internal_state_out.actor = internal_state_in.actor;

% update profiler status
if internal_state_out.workflow.profiler_status  &&  ~internal_state_in.workflow.profiler_status
    profile on
end


