% ----------------------------------------------------------------
% execute fitsepa (to smooth separatrix) for all Cronos time steps
% ----------------------------------------------------------------

npts = 31;
rgeonew = zeros(size(data.geo.R,1),npts);
zgeonew = zeros(size(data.geo.Z,1),npts);
for i=1:size(data.geo.R,1)
  [rgeonew(i,:),zgeonew(i,:)] = fitsepa(data.geo.R(i,:),data.geo.Z(i,:),31);
end
data.geo.R = rgeonew;
data.geo.Z = zgeonew;
clear rgeonew zgeonew npts i




