function [lm,lp] = cxl(psi,RBt,gb,gl,w)

% PSITBXPSI/CXL   connection length of field lines
% CXL(PSI,RBT[,gb,gl])

% optional grid conversion
if nargin < 3, gb = []; end, gb = gf(gb,psi);
if nargin < 4, gl = []; end, gl = gf(gl,psi);
if nargin < 5, wc = 1; end

% interpolate on new grids
psi = psitbxp2p(psi,'+0');
if isa(gb,'psitbxgrid'), pb = psitbxf2f(psi.psitbxfun,gb);
else,                    pb = psi.psitbxfun; gb = pb.grid; end
if isa(gl,'psitbxgrid'), pl = psitbxf2f(psi,gl);
else,                    pl = psi.psitbxfun; gl = pl.grid; end
t = pl.t; nt = length(t); nl = size(pl); ll = prod(nl)/nt;

switch gl.storage(1)
 case 'G'
  [tmp,rzl] = meshgrid(gl.x{2},gl.x{1});
  rzl = repmat(rzl(:) + j*tmp(:),1,nt);
 case 'P'
  rzl = repmat(gl.x{1}(:) + j*gl.x{2}(:),1,nt);
 case 'T'
  rzl = reshape(gl.x{1}(:) + j*gl.x{2}(:),ll,nt);
end
plx = reshape(pl.x,ll,nt);

rb = gb.x{1}; zb = gb.x{2};
[tmp,rzb] = meshgrid(rb,zb); rzb = rzb(:) + j*tmp(:);

[lm,lp] = deal(repmat(Inf,size(plx)));
for kt = 1:nt
 pbx = pb.x(:,:,kt);
 
 cc = contourc(zb,rb,pbx,plx(:,kt));
 % get size of contours
 l = 1; nc = 0; mc = 0;
 while l <= size(cc,2)
  nc = max(nc,cc(2,l));
  l = l + cc(2,l) + 1;
  mc = mc + 1;
 end
 % extract contours
 xc = repmat(NaN+j*NaN,mc,nc);
 pc = repmat(NaN,mc,1);
 l = 1;
 for k = 1:mc
  pc(k) = cc(1,l);
  xc(k,1:cc(2,l)) = [j,1] * cc(:,l+1:l+cc(2,l));
  l = l + cc(2,l) + 1;
 end
 % select closest contours
 [kc,lc] = deal(zeros(ll,1));
 cc = repmat(NaN+j*NaN,ll,nc);
 for k = 1:ll;
  l = find(pc == plx(k,kt));
  if ~isempty(l)
   kc(k) = l(imin(min(abs(xc(l,:) - rzl(k,kt)),[],2)));
   lc(k) = imin(abs(xc(kc(k),:) - rzl(k,kt)));
   cc(k,:) = xc(kc(k),:);
  end
 end
  
 % find flux derivatives on contour points
 k = find(~isnan(cc));
 px = psitbxfun(pbx,gb);
 [pr,pz] = deal(repmat(NaN,size(cc)));
 tmp = psitbxgrid('c','p',{real(cc(k)),imag(cc(k))});
 pr(k) = double(psitbxf2f(px,tmp,[1,0]));
 pz(k) = double(psitbxf2f(px,tmp,[0,1]));
 % weighting
 if nargin == 5
  wc = repmat(NaN,size(cc));
  wc(k) = double(psitbxf2f(psitbxfun(w.x(:,:,kt),gb),tmp));
  wc = (wc(:,1:end-1) + wc(:,2:end))/2;
 end

 % integrate connexion length along each contour
 dr = diff(real(cc),1,2);
 dz = diff(imag(cc),1,2);
 pr = (pr(:,1:end-1) + pr(:,2:end))/2;
 pz = (pz(:,1:end-1) + pz(:,2:end))/2;
 dy = 2*pi*RBt(min(kt,length(RBt)))*(dz.*sign(pr) - dr.*sign(pz)) ./ (abs(pr) + abs(pz));
 dl = sqrt((dr.^2 + dz.^2 + dy.^2)) .* wc;
 dl(isnan(dl)) = 0;
 dy(isnan(dy)) = 0;

 for k = find(any(~isnan(cc),2))'
  lm(k,kt) = sum(dl(k,1:lc(k)-1));
  lp(k,kt) = sum(dl(k,lc(k):end));
 end
 k = find(mean(dy,2) < 0);
 [lm(k,kt),lp(k,kt)] = deal(lp(k,kt),lm(k,kt));;
  
end

if nargout == 1, lm = lm + lp; end
lm = psitbxfun(reshape(lm,nl),gl,t);
lp = psitbxfun(reshape(lp,nl),gl,t);

function g = gf(g,psi)
switch class(g)
 case 'double' % grid given as a rectangle
  if ~isempty(g)
   rg = psi.psitbxfun.grid.x{1};
   zg = psi.psitbxfun.grid.x{2};
   tmp = g;
   rg = [tmp(1) rg(rg > tmp(1) & rg < tmp(2)) tmp(2)];
   zg = [tmp(3) zg(zg > tmp(3) & zg < tmp(4)) tmp(4)];
   g = psitbxgrid('Cylindrical','Grid',{rg,zg});
  end
 case 'cell' % grid given as {R,Z}
  g = psitbxgrid('Cylindrical','Grid',{sort(g{1}),sort(g{2})});
 case 'psitbxgrid' % grid given as a grid
 otherwise, error('Invalid grid specification')
end

