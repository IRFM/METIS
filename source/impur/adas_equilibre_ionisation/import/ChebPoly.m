function [y]=ChebPoly(x,c,y1,y2)

if (nargin<3) ,y1 = [0];end
if (nargin<4) ,y2 = [0];end

nc = length(c);
if (nc == 1)
 y = 0.5*c(1) + x.*y1 - y2;
else
 y0= c(nc) + 2*y1.*x - y2;
 y = ChebPoly(x,c(1:nc-1),y0,y1);
end

