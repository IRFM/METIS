% wrapper matlab for interpos if not compiled
function [yout,youtp,youtpp,youtint]=interpos(xin,yin,xout,taus,nbc,ybc,sigma)

persistent warning_interpos
if isempty(warning_interpos)
        warning_interpos =fprintf('===============>>>>>>> INTERPOS have not be compiled !');
end
yout = interp1(xin,yin,xout,'pchip','extrap');
if size(xout,2) > 1
  d = 2;
else
  d = 1;
end
if nargout > 1
   youtp = pdederive(xout,yout,2,2,d,1);
end
if nargout > 2
   youtpp = pdederive(xout,yout,2,2,d,2);
end
if nargout > 3
   youtint = cumtrapz(xout,yout,d);
end

