function f = diff(f,dim,pad)

% PSITBXFUN/DIFF

if nargin < 3, pad = 0; end
f = scd('d',f,dim,'no',pad);
