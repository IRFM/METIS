function f = psitbxfun(x,g,t)

% PSITBXFUN   PsiTbx Function object constructor
% PSITBXFUN(X,GRID[,T])

% type casting
if nargin == 1 & isa(x,'psitbxfun')
 f = x;
 return
end

if nargin < 1, x = []; end
if nargin < 2, g = psitbxgrid; end
if nargin < 3, t = []; end
if isa(g,'psitbxfun'), t = g.t; g = g.grid; end

dimx = size(x);
if ~isempty(t) & (ndims(t) > 2 | min(size(t)) > 1)
 error('T must be a vector')
end
if length(t) > 1 & dimx(end) ~= length(t)
 error('Last dimension of X incompatible with T')
end
if isempty(sizeck(dimx,size(g)))
 error('X size incompatible with GRID size')
end

f.x    = x;
f.grid = g;
f.t    = t(:)';
f = class(f,'psitbxfun');
