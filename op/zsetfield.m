% compatibilite matlab5/6
function varargout = zsetfield(varargin)

if length(varargin) > 3
   varargout = {setfield(varargin{:})};
elseif isstruct(varargin{1}) & ischar(varargin{2})
   if findstr(varargin{2},'.')
      s            = varargin{1};
      n            = varargin{2};
      v            = varargin{3};
      eval(sprintf('s.%s=v;',n));
      varargout(1) = {s};
   else
      varargout = {setfield(varargin{:})};
   end
else
   varargout = {setfield(varargin{:})};
end
