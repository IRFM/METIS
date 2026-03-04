function [pped,self_dunne] = z0dunne_pedestal(zs,profli,geo,cons,option)
% allocate data 
pped              = NaN * ones(size(cons.temps));
self_dunne.te_sep = NaN * ones(size(cons.temps));
self_dunne.ne_sep = NaN * ones(size(cons.temps));
self_dunne.ne_ped = NaN * ones(size(cons.temps));
self_dunne.te_ped = NaN * ones(size(cons.temps));
% get default parameters
params = dunne_pedestal;
%
% loop on time slices
for k=1:length(cons.temps)
 % update parameters
 params.ip       = cons.ip(k);
 params.bt       = geo.b0(k);
 params.a        = geo.a(k);
 params.rmaj     = geo.R(k);
 params.kappa    = geo.K(k);
 params.delta    = geo.d(k);
 params.betapol  = zs.betaptot(k);
 % not defined             w_pre: 0.1000
 params.density  = zs.nbar(k);
 % not used        dens_shift: 0
 % not used        pres_shift: 0
 % can be exported later in option if needed       nesep_scale: 4
 % not if profile are not computed, some logical problem in the function tesep   
 params.zeff     = zs.zeff(k);
 params.Amain    = zs.meff(k);
 %          profiles: 'EPED'
 params.profiles_onoff =  0;
 % can be exported later in option if needed                 scaling: 'Faitsch'
 params.printinfo  = 0;
 % can be exported later in option if needed            self_consistent: 1
 % can be exported later in option if needed                       npts: 500
 % call the model   
 self = dunne_pedestal(params);
 
 % set data
 pped(k)              = self.p_ped;
 self_dunne.te_sep(k) = self.sep_params.te_sep;
 self_dunne.ne_sep(k) = self.sep_params.ne_sep;
 self_dunne.ne_ped(k) = self.ne_ped;
 self_dunne.te_ped(k) = self.te_ped;

end