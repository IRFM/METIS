function state = metis4muscle3_state_initialisation

% workflow state (initialisation)
state.should_init           = false;  % muscle3 requested
state.metis_init            = false;  % if metis initialisation has been done
state.metis_convergence     = false;  % convergence for zerodevolution
state.metis_rewind          = false;  % retart from previous time slice for zerodevolution
state.metis_evolution       = false;  % evolution for zerod_evolution
state.metis_replay          = false;  % stream a existing METIS simulation
state.metis_computation     = false;  % stream a existing METIS simulation
state.working_point         = false;  % force METIS to run in working point mode
state.one_time              = false;  % if METIS should compute one time after the other.
state.metis_exit            = false;  % final call of the actor
state.checkpoint            = false;  % muscle3 request for checkpointing 
state.final_checkpoint      = false;  % muscle3 request for final checkpointing on exit
state.resuming              = false;  % muscle3 request for simulation cold restarting
state.final_checkpoint_done = false; % if final checkpoint data ha been saved
state.time_source           = 'none'; % source for actor synchronisation on time.
