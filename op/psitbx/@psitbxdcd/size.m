function y = size(dcd,n)

% PSITBXDCD/SIZE	PsiTbx Diagnostic-Chord-Description size
%   SIZE(DCD) returns the size of the PSITBXDCD object DCD

y = size(dcd.rd);
if nargin > 1, y = y(n); end
