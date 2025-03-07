function [rnew,znew] = change_sepa_points_number(rold,zold,nreduc)
%% ----------------------------------------------------------------------
%% REDUCE THE NUMBER OF POINTS OF THE SEPARATRIX:
%% [rnew,znew] = change_sepa_points_number(data.geo.R,data.geo.Z,nreduc)
%% ----------------------------------------------------------------------

ntime = size(rold,1);

for it=1:ntime
  rnew(it,:) = rold(it,1:nreduc:end);
  znew(it,:) = zold(it,1:nreduc:end);
end

