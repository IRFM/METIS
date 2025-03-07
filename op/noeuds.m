function [x,y,z]=noeuds(a,b,c)
%
% enleve les noueds en double de a et b
%
a = a(:);
b = b(:);
c = c(:);

co = a + i*b;

[nco,ind,jnd] = unique(co);

x = real(nco);
y = imag(nco);

z = c(ind);

% elimination du zero

ind=find(x==0 & y==0 & z ~= 0);
if ~isempty(ind)
  x(ind) = [];
  y(ind) = [];
  z(ind) = [];
end
