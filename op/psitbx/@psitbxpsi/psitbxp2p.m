function psi = psitbxp2p(psi,mode,varargin)

% PSITBXPSI/PSITBXP2P

switch mode
 case {'-0' '+0'}
  switch psi.format
   case {'-0' '+0'}
    if ~strcmp(mode,psi.format), psi.psitbxfun = -psi.psitbxfun; end
   case '01', psi.psitbxfun = (1 - psi.psitbxfun) .* ...
    repmat(psi.psimag,[size(psi.psitbxfun.x(1)),size(psi.psitbxfun.x(2)),1]);
   case 'FS', error('"Flux-Surfaces" can not be converted')
  end
 case '01'
  psi = psitbxmag(psi);
 case 'FS'
  psi = psitbxfsd(psi,varargin{:});
 otherwise
  error('Invalid mode')
end
psi.format = mode;
