% compatibilite matlab5/6/7
function varargout = zgetfield(varargin)

varargout{1} = [];

if length(varargin) > 2
   varargout = {getfield(varargin{:})};
elseif isstruct(varargin{1}) & ischar(varargin{2})
   if findstr(varargin{2},'.')
      s            = varargin{1};
      n            = varargin{2};
      varargout(1) = {eval(sprintf('s.%s',n))};
   elseif isfield(varargin{:})
      varargout = {getfield(varargin{:})};
   end
elseif isfield(varargin{:})
   varargout = {getfield(varargin{:})};
end

