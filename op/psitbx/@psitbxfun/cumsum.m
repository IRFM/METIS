function f = cumsum(f,dim,varargin)

% PSITBXFUN/CUMSUM

f = scd('c',f,dim,varargin{:});
