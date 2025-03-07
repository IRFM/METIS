%% SECURITY ON IDS OCCURENCES
function option = secure_option_imas(option)

option.pulse_schedule_occurrence = secure_occurrence_imas(option.pulse_schedule_occurrence);
option.summary_occurrence        = secure_occurrence_imas(option.summary_occurrence);
option.core_profiles_occurrence  = secure_occurrence_imas(option.core_profiles_occurrence);
option.core_transport_occurrence = secure_occurrence_imas(option.core_transport_occurrence);
option.equilibrium_occurrence    = secure_occurrence_imas(option.equilibrium_occurrence);
option.core_sources_occurrence   = secure_occurrence_imas(option.core_sources_occurrence);
option.radiation_occurrence      = secure_occurrence_imas(option.radiation_occurrence);
option.edge_occurrence           = secure_occurrence_imas(option.edge_occurrence);
option.numerics_occurrence       = secure_occurrence_imas(option.numerics_occurrence);
