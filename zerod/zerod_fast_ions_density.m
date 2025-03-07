% densite de particule rapides
function [n_fast,source_fast,emoy_fast,tau_nfast,ecrit] = zerod_fast_ions_density(nep,tep,zeff,meff,Afast,Zfast,Einj,p_supra)
% temps de thermalisation 
lnldei      = 15.2 - 0.5 .* log(nep./1e20) + log(tep ./1e3);
taus        = 6.27e8 .* Afast ./ Zfast .^ 2 .* tep  .^ (3/2) ./ (nep ./ 1e6) ./ lnldei;
ecrit       = 14.8 .* Afast .* tep  .* (zeff./ meff) .^ (2/3);
source_fast = 2 .* p_supra ./ (Einj .* 1.602176462e-19 .* taus);
tau_nfast   = taus ./ 3 .* log(1 + (Einj ./ ecrit) .^(3/2));           % temps d'arret des alpha
n_fast      = source_fast .* tau_nfast;                                          % nombre d'alpha dans la decharge nonthermalise
emoy_fast   = max(tep,min(Einj,p_supra ./ max(1,n_fast)./  1.602176462e-19));  % energy moyenne en eV pour F(x)
