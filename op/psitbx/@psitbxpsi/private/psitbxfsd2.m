function psi = psitbxfsd2(psi,gs)

% PSITBXPSI/PSITBXFSD	Flux Surface Description
%   PSITBXFSD(PSI)

if nargin < 2, gs = psitbxgrid('Flux','Grid','Default'); end

% make sure psi is normalised
psi = psitbxmag(psi);

x = psi.psitbxfun.x;
g = psi.psitbxfun.grid;
[nr,nz,nt] = size(x);
[nrh,nth] = size(gs);
d = psitbxdcd(repmat(psi.rmag,nth,1),repmat(psi.zmag,nth,1),[],repmat(gs.x{2}(:),1,nt),[],[],nt);
