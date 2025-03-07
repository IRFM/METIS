function [sl,sh,pc] = psitbxdxp(dcd,psi,rh)

% PSITBXDCD/PSITBXDXP   Intersection between chords and flux surfaces

% check inputs
if any(dcd.tvd(:)) | length(dcd.t)
 error('Only for fixed poloidal diagnostic chords')
end

% default inputs, sizes and alloc
nd = size(dcd); Nd = prod(nd);
psi = psitbxp2p(psi,'01'); x = psi.x; [nr,nz,nt] = size(x);
g = psi.grid; Mx = g.csp.Mx; xk = g.csp.xk;
if nargin < 3, rh = psitbxgrid('Flux','Grid','Default'); rh = rh.x{1}; end
ps = rh.^2; nps = length(ps);
[sl,sh] = deal(repmat(NaN,[nps,Nd,nt]));

% fuq
sid  = sin(dcd.pvd(:));
cod  = cos(dcd.pvd(:));
s2d  = sid.^2;
c2d  = cod.^2;
scd  = 2*sid.*cod;
rdps = repmat(dcd.rd(:),1,nps);
zdps = repmat(dcd.zd(:),1,nps);
sips = repmat(sid      ,1,nps);
cops = repmat(cod      ,1,nps);
s2ps = repmat(s2d      ,1,nps);
c2ps = repmat(c2d      ,1,nps);
scps = repmat(scd      ,1,nps);
psps = repmat(ps(:)'   ,Nd,1);

% extrem values
s = intersect(dcd,g,2); smin = s(1,:)'; smax = s(2,:)';

% closest approach initial guess
s  = intersect(dcd,g);
pc = double(psitbxf2f(psi.psitbxfun,psitbxgrid('D','P',{s},dcd)));
sc = repmat(s,[ones(1,ndims(s)),nt]);
[pc,k] = min(pc,[],1); pc = reshape(pc,[Nd,nt]);
sc = reshape(sc,dcd.nd,Nd*nt);
sc = reshape(sc(sub2ind(size(sc),k(:)',1:Nd*nt)),[Nd,nt]);
 
for kt = 1:nt
 % spline coefficients
 c = xmtimes(xmtimes(Mx{1},x(:,:,kt)),Mx{2}');
 
 % closest approach
 sn = sc(:,kt); ds = Inf;
 % iterations
 for kiter = 1:psi.iter
  if max(abs(ds)) < psi.tol, break, end
  s = sn;
  r = dcd.rd(:) - s.*cod; z = dcd.zd(:) + s.*sid;
  cr0 = bspsum(xk{1},c,r,0,0,1)';
  cr1 = bspsum(xk{1},c,r,1,0,1)';
  cr2 = bspsum(xk{1},c,r,2,0,1)';
  dp  = bspsum(xk{2},cr0,z,1,1,1) .* sid - ...
        bspsum(xk{2},cr1,z,0,1,1) .* cod;
  d2p = bspsum(xk{2},cr2,z,0,1,1) .* c2d - ...
        bspsum(xk{2},cr1,z,1,1,1) .* scd + ...
	bspsum(xk{2},cr0,z,2,1,1) .* s2d;
  ds = - dp ./ d2p;
  sn = s + ds;
  ds(sn < smin(:) | sn > smax(:)) = 0;
  sn = min(max(sn,smin),smax);
 end
 p = bspsum(xk{2},cr0,z,0,1,1);
 % remove outside plasma
 k = find(d2p < 0 | abs(dp) > psi.tol | p > 1); s(k) = NaN; p(k) = NaN;
 % dispatch in output
 sc(:,kt) = s; pc(:,kt) = p;
 
 % intersections, initial guess
 s   = repmat(s  ,1,nps);
 p   = repmat(p  ,1,nps);
 dp  = repmat(dp ,1,nps);
 d2p = repmat(d2p,1,nps);
 ds = sqrt(max(dp.^2+2*d2p.*(psps-p),0));
 % remove flux surfaces which does not cross dcd
 k = find(p > psps); ds(k) = NaN;
 s0{1} = s - (ds + dp) ./ d2p;
 s0{2} = s + (ds - dp) ./ d2p;
 % Gauss-Newton iterations
 for k = 1:2
  s = s0{k};
  for kiter = 1:psi.iter
   r = rdps - s.*cops;
   z = zdps + s.*sips;
   cr0 = permute(bspsum(xk{1},c,r,0,0,1),[3,1,2]);
   cr1 = permute(bspsum(xk{1},c,r,1,0,1),[3,1,2]);
   cr2 = permute(bspsum(xk{1},c,r,2,0,1),[3,1,2]);
   p   = bspsum(xk{2},cr0,z,0,1,1);
   dp  = bspsum(xk{2},cr0,z,1,1,1) .* sips - ...
         bspsum(xk{2},cr1,z,0,1,1) .* cops;
   d2p = bspsum(xk{2},cr2,z,0,1,1) .* c2ps - ...
         bspsum(xk{2},cr1,z,1,1,1) .* scps + ...
	 bspsum(xk{2},cr0,z,2,1,1) .* s2ps;
   ds = (sqrt(max(dp.^2+2*d2p.*(psps-p),0)).*sign(dp) - dp) ./ d2p;
   s = s + ds;
   if max(abs(ds(:))) < psi.tol, break, end
  end
  s0{k} = s;
 end
 sl(:,:,kt) = s0{1}'; sh(:,:,kt) = s0{2}';

end

% outputs
sl = psitbxgrid('D','T',{reshape(sl,[nps,nd,nt])},dcd,psi.t);
sh = psitbxgrid('D','T',{reshape(sh,[nps,nd,nt])},dcd,psi.t);
pc = psitbxfun(reshape(pc,[nd,nt]),psitbxgrid('D','T',{reshape(sc,[nd,nt])},dcd,psi.t));
