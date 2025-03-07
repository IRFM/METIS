function varargout = size(f,dim)

% PSITBXFUN/SIZE	PsiTbx-Function size

dimf = sizeck(size(f.x),size(f.grid));
if nargin > 1
 if strcmp(dim,'end')
  dimf = dimf(end);
 else
  dimf = dimf(dim);
 end
end
if nargout > 1
 dimf = [dimf,ones(1,nargout-length(dimf))];
 for k = 1:nargout-1
  varargout{k} = dimf(k);
 end
 varargout{nargout} = prod(dimf(nargout:end));
else
 varargout{1} = dimf;
end
