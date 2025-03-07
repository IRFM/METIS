function [T,TT] = tmat(dcd,g,b,a)

% PSITBXDCD/TMAT   Transfer matrix between chords and cylindrical grid
%   T = TMAT(DCD,G) returns the transfer matrix such that if X is a
%   function defined on G, Y=T*X gives the signal seen by DCD assuming a
%   constant emissivity on each grid cell of G. TMAT(DCD,G,B) limits the view 
%   within the grid B or the box B=[R1,R2,Z1,Z2].
%   T is a sparse matrix of size [PROD(SIZE(DCD)),PROD(SIZE(G)-1)].  
%   [T,TT] = TMAT(DCD,G) computes also the TT matrix such that X=TT*(T'*Y) is
%   the l.s. solution of Y=T*X. TMAT(DCD,G,B,[A,B]) regularises the inversion
%   such that [Y;0;0]=[T*X;A*X;B*SUM(T).*X] is solved in a l.s.s.

if nargin < 3, b = [-Inf,Inf,-Inf,Inf]; end
if g.type(1) ~= 'C' | g.storage(1) ~= 'G'
 error('Only for "Cylindrical" coordinates with "Grid" storage')
end
rsub = g.x{1}; zsub = g.x{2};
nr = length(rsub) - 1; nz = length(zsub) - 1;
if isa(b,'psitbxgrid')
 if b.type(1) ~= 'C' | b.storage(1) ~= 'G'
  error('Only for "Cylindrical" coordinates with "Grid" storage')
 end
 r1 = b.x{1}(1); r2 = b.x{1}(end);
 z1 = b.x{2}(1); z2 = b.x{2}(end);
else
 r1 = b(1); r2 = b(2); z1 = b(3); z2 = b(4); 
end

% search intersection with boundary
dimd = size(dcd.rd);
s = intersect2dmex(prod(dimd),dcd.rd,dcd.zd,dcd.pvd,dcd.tvd,r1,r2,z1,z2);
smax = s(2,:);

T = sparse(prod(size(dcd)),nr*nz);
for kr = 1:nr
 for kz = 1:nz
  s = intersect2dmex(prod(dimd),dcd.rd,dcd.zd,dcd.pvd,dcd.tvd,...
      rsub(kr),rsub(kr+1),zsub(kz),zsub(kz+1));
  l = max(s(2,:) - s(1,:),0) .* (s(1,:) < smax) + ...
      (s(4,:) - s(3,:)) .* (s(3,:) < smax);
  T(:,sub2ind([nr,nz],kr,kz)) = l(:);
 end
end
if nargout > 1 % computes TT matrix
 if nargin < 4, a = []; end
 a = [a [0 0]];
 if all(a) == 0 & any(~any(T))
  warning('Grid cell(s) not seen by any chord: add non zero regularisation')
 end
 a1 = a(1); if a1 == 0, a1 = zeros(0,nr*nz); end
 a2 = a(2);
 if a2 == 0, a2 = zeros(0,nr*nz);
 else,       a2 = a2 / sqrt(sum(sum(T).^2)/nr/nz); end
 TT = qr([T; a1*speye(nr*nz); a2*spdiags(sum(T)',0,nr*nz,nr*nz)],0);
 TT = TT\(TT'\eye(size(TT)));
end
