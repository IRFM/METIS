% derivation with polynomial
% x = linspace(0,2*pi,101);
% [d1,d2] = z0polyderive(x,sin(x),7,x)
% figure;plot(x,cos(x),'.r',x,-sin(x),'.b',x,d1,'r',x,d2,'b')
%
% x = linspace(0,2*pi - pi/4,11);
% xx = linspace(0,2*pi,101);
% [d1,d2] = z0polyderive(x,sin(x),7,xx)
% figure;plot(xx,cos(xx),'.r',xx,-sin(xx),'.b',xx,d1,'r',xx,d2,'b')

function [d1,d2] = z0polyderive(xin,v,order,xout)

if nargin < 4
  xout = xin
end

pp   = polyfit(xin,v,order);
ppd1 = polyder(pp);
if nargout  > 1
  ppd2 = polyder(ppd1);
end
  

d1   = polyval(ppd1,xout);
if nargout  > 1
  d2   = polyval(ppd2,xout);
end