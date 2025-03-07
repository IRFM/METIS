function psi = psitbxmag(psi)

% PSITBXPSI/PSITBXMAG

if ~isempty(psi.rmag) & strcmp(psi.format,'01'), return, end

x = psi.psitbxfun.x;
g = psi.psitbxfun.grid;
switch(psi.format)
 case '-0', x = -x;
 case '01', x = 1 - x;
 case 'FS', error('"Flux-Surfaces" can not be normalised')
end

% alloc
[nr,nz,nt] = size(x);
[psi.rmag,psi.zmag,psi.psimag] = deal(repmat(NaN,1,nt));

Mx = g.csp.Mx;
xk = g.csp.xk;
for kt = 1:nt

 % spline coefficients
 c = xmtimes(xmtimes(Mx{1},x(:,:,kt)),Mx{2}');
 
 % initial guess
 [kr,kz] = find(x(:,:,kt) == max(max(x(:,:,kt),[],1)));
 r = g.x{1}(kr(1)); z = g.x{2}(kz(1)); % (1) in case of two maxima

 % Gauss-Newton iteration
 for kiter = 1:psi.iter
 
  cr0 =  bspsum(xk{1},c,  r,0,0,1)';
  cr1 =  bspsum(xk{1},c,  r,1,0,1)';
  dp  = [bspsum(xk{2},cr1,z,0,0,1);
         bspsum(xk{2},cr0,z,1,0,1)];
  d2p = [bspsum(xk{2},bspsum(xk{1},c,r,2)',z,0,0,1), ...
         bspsum(xk{2},cr1,z,1,1); ...
         0, ...
	 bspsum(xk{2},cr0,z,2,0,1)];
  d2p(2,1) = d2p(1,2);
  drz = d2p \ dp;
  r = r - drz(1); z = z - drz(2);
  if max(abs(drz)) < psi.tol, break, end
 end
 
 if kiter < psi.iter
  psi.rmag(kt) = r; psi.zmag(kt) = z;
  psi.psimag(kt) = bspsum(xk{2},bspsum(xk{1},c,r,0,0,1)',z,0,1,1);
  x(:,:,kt) = 1 - x(:,:,kt) / psi.psimag(kt);
 else
  x(:,:,kt) = NaN;
 end
end

psi.format = '01';
psi.psitbxfun = psitbxfun(x,g,psi.psitbxfun.t);
