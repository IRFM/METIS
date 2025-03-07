%% ------------------------------------------------------------------
%% Written by M. Schneider - 05/11/07
%% ------------------------------------------------------------------
%%
%% PURPOSE: In the CRONOS general workspace, this routine corrects
%% possible pi variations in the separatrix coordinates (data.geo.R
%% and data.geo.Z) from the time "k" to the time "k+1"
%% 
%% This allows to avoid following error messages:
%% MONOTONIEFEHLER IN SPLINE: X(I-1) >= X(I) ??? FEHLER IN ROUTINE
%% ------------------------------------------------------------------

for it=1:size(data.geo.R,1)-1
  idxmaxref = find(data.geo.Z(it,:)==max(data.geo.Z(it,:)));
  if data.geo.Z(it+1,idxmaxref)<0
    
    disp(['Inversion for index ',num2str(it)])
    
    idxmaxnew = find(data.geo.Z(it+1,:)==max(data.geo.Z(it+1,:)));
    idxshift  = abs(idxmaxref-idxmaxnew);
    idxlim    = size(data.geo.R,2)-idxshift;

    rpnew = [data.geo.R(it+1,idxlim:end) data.geo.R(it+1,1:idxlim-1)];
    zpnew = [data.geo.Z(it+1,idxlim:end) data.geo.Z(it+1,1:idxlim-1)];

    data.geo.R(it+1,:) = rpnew;
    data.geo.Z(it+1,:) = zpnew;
  
  end
end
clear it idxmaxnew idxshift idxlim idxmaxref rpnew zpnew
