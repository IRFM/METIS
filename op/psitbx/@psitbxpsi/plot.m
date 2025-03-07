function plot(psi,varargin)

% PSITBXPSI/PLOT

if strcmp(psi.format,'FS')
 nrh = length(psi.psitbxfun.grid.x{1});
 rs = psi.rmag(1) - psi.psitbxfun.x(:,:,1) .* repmat(cos(psi.psitbxfun.grid.x{2}(:)'),nrh,1);
 zs = psi.zmag(1) + psi.psitbxfun.x(:,:,1) .* repmat(sin(psi.psitbxfun.grid.x{2}(:)'),nrh,1);
 plot(rs',zs','-r',varargin{:})
 axis equal
else
 plot(psi.psitbxfun,varargin{:})
end
title(['PsiTbx-Poloidal-Flux object "' inputname(1) '" (' psi.psitbxfun.grid.type ',' psi.psitbxfun.grid.storage ',' psi.format ')'])

