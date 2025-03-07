function varargout = size(g,dim)

% PSITBXGRID/SIZE	Grid size

if g.storage(1) == 'G'
 dimg = zeros(1,length(g.x));
 for k = 1:length(g.x)
  dimg(k) = length(g.x{k});
 end
else
 dimg = size(g.x{1});
end
if nargin > 1
 if strcmp(dim,'end')
  dimg = dimg(end);
 else
  dimg = dimg(dim);
 end
end
if nargout > 1
 dimg = [dimg,ones(1,nargout-length(dimg))];
 for k = 1:nargout-1
  varargout{k} = dimg(k);
 end
 varargout{nargout} = prod(dimg(nargout:end));
else
 varargout{1} = dimg;
end
