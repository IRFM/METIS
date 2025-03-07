function psi = psitbxfsd(psi,gs)

% PSITBXPSI/PSITBXFSD	Flux Surface Description
%   PSITBXFSD(PSI)

if nargin < 2, gs = psitbxgrid('Flux','Grid','Default'); end

% make sure psi is normalised
psi = psitbxmag(psi);

x = psi.psitbxfun.x;
g = psi.psitbxfun.grid;

% alloc
[nr,nz,nt] = size(x);
[nrh,nth] = size(gs);
as = repmat(NaN,[nrh,nth,nt]);

% fuq
sith = repmat(sin(gs.x{2}(:)')                    ,nrh,1);
coth = repmat(cos(gs.x{2}(:)')                    ,nrh,1);
s2th = repmat(sin(gs.x{2}(:)').^2                 ,nrh,1);
c2th = repmat(cos(gs.x{2}(:)').^2                 ,nrh,1);
scth = repmat(sin(gs.x{2}(:)').*cos(gs.x{2}(:)')*2,nrh,1);
psth = repmat(gs.x{1}(:).^2                       ,1,nth);
k0 = find(gs.x{1} == 0);

Mx = g.csp.Mx; xk = g.csp.xk;
ki = find(abs(diff(psi.rmag)) < (g.x{1}(2) - g.x{1}(1))/4 & ...
          abs(diff(psi.zmag)) < (g.x{2}(2) - g.x{2}(1))/4) + 1;
for kt = 1:nt

 % spline coefficients
 c = xmtimes(xmtimes(Mx{1},x(:,:,kt)),Mx{2}');

 % initial guess, previous a when two equilibria are similar, avoid 0
 if any(kt == ki), a = max(.9,1-1/nrh) * a;
 else,             a = repmat(psi.tol,nrh,nth); end

 % Gauss-Newton iteration
 for kiter = 1:psi.iter
  r = psi.rmag(kt) - a .* coth;
  z = psi.zmag(kt) + a .* sith;
  cr0 = permute(bspsum(xk{1},c,r,0,0,1),[3,1,2]);
  cr1 = permute(bspsum(xk{1},c,r,1,0,1),[3,1,2]);
  cr2 = permute(bspsum(xk{1},c,r,2,0,1),[3,1,2]); 
  p   = bspsum(xk{2},cr0,z,0,1,1);
  dp  = bspsum(xk{2},cr0,z,1,1,1) .* sith - ...
        bspsum(xk{2},cr1,z,0,1,1) .* coth;
  d2p = bspsum(xk{2},cr2,z,0,1,1) .* c2th - ...
        bspsum(xk{2},cr1,z,1,1,1) .* scth + ...
        bspsum(xk{2},cr0,z,2,1,1) .* s2th;
  da = (sqrt(max(dp.^2+2*d2p.*(psth-p),0)).*sign(dp) - dp) ./ d2p;
  a = a + da;
  if max(abs(da(:))) < psi.tol, break, end
 end

 % keep 0 on 0
 a(k0,:) = 0;
 
 % dispatch in outputs
 if kiter < psi.iter
  as(:,:,kt) = a;
 end
end
psi.format = 'FS';
psi.psitbxfun = psitbxfun(as,gs,psi.psitbxfun.t);
