function [x,storage] = p2t(x,storage,nt)

if storage(1) == 'P' & nt > 1
 dim = [];
 for k = 1:length(x)
  if isempty(dim) & length(x{k}) > 1
   dim = size(x{k});
  end
 end
 for k = 1:length(x)
  if length(x{k}) > 1
   x{k} = repmat(x{k},[ones(1,length(dim)),nt]);
  elseif ~isnan(x{k})
   x{k} = repmat(x{k},[dim,nt]);
  end
 end
 storage = 'T';
end
