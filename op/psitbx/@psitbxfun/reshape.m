function [f,ks] = reshape(f,varargin)

% PSITBXFUN/RESHAPE	PsiTbx-Function reshape
% F = RESHAPE(F,DIM)
% [F,KS] = RESHAPE(F)

switch nargin
 case 1, [dim,ks] = sizeck(size(f.x),size(f.grid));
 case 2, dim = varargin{1};
 otherwise, dim = [varargin{:}];
end
if isempty(sizeck(dim,size(f.grid)))
 error('New size incompatible with function grid size')
end
f.x = reshape(f.x,dim);
