function [n,knan] = ndims(g)

% PSITBXGRID/NDIMS	Grid dimension

n = length(g.x);
knan = zeros(1,n);
for k = 1:n
 knan(k) = any(isnan(g.x{k}(:)));
end
n = n - sum(knan);
knan = find(knan);
