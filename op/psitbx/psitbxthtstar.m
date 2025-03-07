function f = psitbxthtstar(psi)

% PSITBXPSI/PSITBXTHTSTAR	straight field line angle

% THTSTAR = PSITBXTHTSTAR(PSI)
%
% MODIFICATION HISTORY:
%   22-03-01  HR  introduce eps to avoid division by zero.

g = psitbxgrid('Flux','Grid','Default');

% make sure psi contains rmag and zmag
%psi = psitbxmag(psi);
psi = psitbxp2p(psi,'01');

t   = psi.psitbxfun.t;
nt  = length(t);
nth = length(g.x{2});
k0  = repmat(reshape([1:nt],[1,1,nt]),size(g));

% cylindrical grid (with Time-Points storage)
gc  = psitbxg2g(g,'C',psitbxp2p(psi,'FS',g));

dpsidr = psitbxf2f(psi,gc,[1,0]);
dpsidz = psitbxf2f(psi,gc,[0,1]);

% integration of (R*J(theta,psi))^-1 over theta
F = cumsum( ...
  psitbxfun( ((gc.x{1}-psi.rmag(k0)+eps).^2 + (gc.x{2}-psi.zmag(k0)+eps).^2 )./ ...
    ( gc.x{1} .* ( (gc.x{1}-psi.rmag(k0)+eps).*dpsidr.x ...
	                 +(gc.x{2}-psi.zmag(k0)+eps).*dpsidz.x ) ), g,t), ...
  2);

% normalize to -pi on LFS
f = F.*(2*pi./repmat(F.x(:,nth,:),[1,nth,1]))-pi;
