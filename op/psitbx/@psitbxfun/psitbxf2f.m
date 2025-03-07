function f = psitbxf2f(f,g2,der)

% PSITBXFUN/PSITBXF2F
% PSITBXF2F(F[,G][,DER])

if nargin < 3
 switch class(g2)
  case 'psitbxgrid', der = [0,0,0];
  case 'double', der = g2; g2 = f.grid;
  otherwise, error('Invalid sencond argument class')
 end
end

% short cut if possible
if all(der == 0) & f.grid == g2, return, end

if f.grid.storage(1) ~= 'G'
 error('Only a function defined on a "Grid" can be mapped')
end

g1 = f.grid;
g21 = psitbxg2g(g2,g1.type,g1);
x21 = g21.x;

% check grid compatibility and push edge and outside points
for k = 1:min(length(x21),3)
 if length(g1.x{k}) == 1 & ~isnan(g1.x{k}) & length(x21{k}) == 1 & isnan(x21{k})
  error(sprintf('Function varies along coordinate #%d',k))
 else
  tmp = g1.x{k}(1) - x21{k};
  x21{k}(tmp > 0 & tmp < 10*eps) = g1.x{k}(1);
  x21{k}(tmp >= 10*eps) = NaN;
  tmp = x21{k} - g1.x{k}(end);
  x21{k}(tmp >= 10*eps) = NaN;
  x21{k}(tmp > 0 & tmp < 10*eps) = g1.x{k}(end);
 end
end

% map
[f,ks] = reshape(f);
z = f.x;
dimf = size(f); ndimf = length(dimf);	
if g21.storage(1) == 'G'			% z = [n1,n2,n3,Ng,Nf[,nt]]
 kper = [2,1,3;3,2,1;2,3,1];			% g = [n1',n2',n3',Ng[,nt]]
 for k = 1:3
  if length(g1.x{k}) == 1
   z = repmat(z,[length(x21{k}),1]);
  else
   z = xmtimes(g1.csp.Mx{k},z);
   z = bspsum(g1.csp.xk{k},z,x21{k},der(k),0,1);
  end
  z = permute(z,[kper(k,:),4:ndimf]);
 end
else
 dimg21 = size(g21);
 if length(dimg21) == 2
  switch sum(dimg21 == 1)
   case 1, dimg21(dimg21 == 1) = [];
   case 2, dimg21 = 1;
  end
 end
 ndimg21 = length(dimg21);			% g = [Ng[,nt]]
 % k = 1					% z = [n1,n2,n3,Nf[,nt]]
 if length(f.t) > 1 & length(g21.t) > 1
  pNg = prod(dimg21(1:end-1)); nt = dimg21(end);
  zz = zeros(pNg,nt,prod(dimf(2:end-1)));	% [prod(Ng),nt,n2*n3*prod(Nf)]
  z = permute(z,[1,ndimf,2:ndimf-1]);		% [n1,nt,n2,n3,Nf]
  if length(g1.x{1}) == 1
   for kt = 1:nt
    zz(:,kt,:) = repmat(z(1,kt,:),[pNg,1,1]);	% [prod(Ng),1,n2*n3*prod(Nf)]
   end
  else
   tmp = reshape(x21{1},[pNg,nt]);		% [prod(Ng),nt]
   for kt = 1:nt
    zz(:,kt,:) = ...				% [prod(Ng),1,n2*n3*prod(Nf)]
     bspsum(g1.csp.xk{1},...
      xmtimes(g1.csp.Mx{1},z(:,kt,:)),...	% [n1,1,n2,n3,Nf]
      tmp(:,kt),...				% [prod(Ng),1]
      der(1),0,1);
   end
  end
  z = reshape(zz,[dimg21,dimf(2:end-1)]);	% [Ng,nt,n2,n3,Nf]
 else
  if length(g1.x{1}) == 1
   z = repmat(z,[prod(dimg21),1]);		% [prod(Ng)[*nt],n2,n3,Nf[,nt]]
   z = reshape(z,[dimg21,dimf(2:end)]);		% [Ng[,nt],n2,n3,Nf[,nt]]
  else
   z = xmtimes(g1.csp.Mx{1},z);
   z = bspsum(g1.csp.xk{1},z,x21{1},der(1),0,1);	% [Ng[,nt],n2,n3,Nf[,nt]]
  end
 end
 for k = 2:3
  z = permute(z,...				% [n2,Ng[,nt],n3,Nf[,nt]] (k=2)
   [ndimg21+1,1:ndimg21,ndimg21+2:ndims(z)]);	% [n3,Ng[,nt],Nf[,nt]]	  (k=3) 
  if length(g1.x{k}) == 1
   z = shiftdim(z,1);				% [Ng[,nt][,n3]Nf[,nt]]
  else
   z = xmtimes(g1.csp.Mx{k},z);			% [Ng[,nt],n3,Nf[,nt]] (k=2)
   z = bspsum(g1.csp.xk{k},z,x21{k},der(k),1,1);	% [Ng[,nt],Nf[,nt]]    (k=3)
  end
 end
 if length(g21.t) > 1
  z = permute(z,[1:ndimg21-1,ndimg21+1:ndims(z),ndimg21]); % [Ng,Nf,nt]
 end
end
z = squeeze(z);

t = f.t;
if isempty(t), t = g21.t; end
f = psitbxfun(z,g2,t);
