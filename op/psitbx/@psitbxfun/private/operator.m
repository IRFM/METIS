function y = operator(op,varargin)

gridfound = 0;
for k = 1:nargin-1
 if ~isa(varargin{k},'double')
  if gridfound
   if g ~= varargin{k}.grid, error(['Incompatible grids using ' op]), end
  else
   gridfound = 1;
   g = varargin{k}.grid;
   y = varargin{k};
  end
  varargin{k} = varargin{k}.x;
 end
 varargin{k} = squeeze(varargin{k});
end
y.x = feval(op,varargin{:});  
