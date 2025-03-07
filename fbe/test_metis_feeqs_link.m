% test metis_feeqs_link
addpath ./certification/metis
load parameters.mat
%parameters.restart = 0;
load feeqs_data4metis.mat
%% metis call
[metis_data4feeqs,internal_state,sepa_metis] = metis_feeqs_link([],feeqs_data4metis,[],parameters);
feeqs_data4metis.time	 = 605; %current time  [1] (s)
feeqs_data4metis.Ip      = metis_data4feeqs.ip;
feeqs_data4metis.Bp_error_vacuum = 0;

feeqs_data4metis.R_LCFS  = sepa_metis.Rsepa;
feeqs_data4metis.Z_LCFS  = sepa_metis.Zsepa;
feeqs_data4metis.Psi     = internal_state.profil0d.psi;
feeqs_data4metis.Raxe    = internal_state.profil0d.Raxe;
feeqs_data4metis.K       = internal_state.profil0d.kx;
[metis_data4feeqs2,internal_state2] = metis_feeqs_link(internal_state,feeqs_data4metis,[],parameters); 
zplotstruct(metis_data4feeqs,metis_data4feeqs2)
