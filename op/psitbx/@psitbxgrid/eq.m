function y = eq(a,b)

% PSITBXGRID/EQ

if isa(a,'psitbxgrid') & isa(b,'psitbxgrid')
 y = strcmp(a.type,b.type) & ...
     strcmp(a.storage,b.storage) & ...
     isequal(size(a.x),size(b.x)) & ...
     (isempty(a.t) | isempty(b.t) | isequal(a.t,b.t)) & ...
     (isempty(a.par.r0) | isempty(b.par.r0) | isequal(a.par.r0,b.par.r0)) & ...
     (isempty(a.par.z0) | isempty(b.par.z0) | isequal(a.par.z0,b.par.z0));
 for k = 1:length(a.x)*y % a trick to avoid this loop if y is already 0
  if any(~isnan(a.x{k}(:))) | any(~isnan(b.x{k}(:)))
   y = y & isequal(a.x{k},b.x{k});
  end 
 end
else
 y = logical(0);
end
