function s = intersect(dcd,g,nd)

% PSITBXDCD/INTERSECT   Intersection between chords and cylindrical grid
%   S=INTERSECT(DCD,G) returns DCD.ND equispaced chord linear coordinates
%   covering the intersection of the chords and the cylindrical grid G.
%   S=INTERSECT(DCD,G,ND) computes ND coordinates. Size of S is [ND,SIZE(DCD)].
%   G can be replaced by the quadruplet [R1,R2,Z1,Z2].

if nargin < 3, nd = dcd.nd; end
if isa(g,'psitbxgrid')
 if g.type(1) ~= 'C' | g.storage(1) ~= 'G'
  error('Only for "Cylindrical" coordinates with "Grid" storage')
 end
  r1 = g.x{1}(1); r2 = g.x{1}(end);
  z1 = g.x{2}(1); z2 = g.x{2}(end);
else
 r1 = g(1); r2 = g(2); z1 = g(3); z2 = g(4); 
end

% search intersection with boundary
dimd = size(dcd.rd);
s = intersect2dmex(prod(dimd),dcd.rd,dcd.zd,dcd.pvd,dcd.tvd,r1,r2,z1,z2);

s1 = reshape(s(1,:),[1,dimd]);
s2 = reshape(s(2,:),[1,dimd]);
s1(s1 > s2) = NaN;
s = reshape(repmat(linspace(0,1,nd)',[prod(dimd),1]),[nd,dimd]) .* ...
    repmat(s2-s1,[nd,1]) + repmat(s1,[nd,1]);
