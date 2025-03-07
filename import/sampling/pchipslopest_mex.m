function d = pchipslopest_mex(x,y,del)
%PCHIPSLOPES  Derivative values for shape-preserving Piecewise Cubic Hermite
% Interpolation.
% d = pchipslopes(x,y,del) computes the first derivatives, d(k) = P'(x(k)).

%  Special case n=2, use linear interpolation.

   n = length(x);
   if n==2  
      d = repmat(del(1),size(y));
      return
   end

%  Slopes at interior points.
%  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
%  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

   d = zeros(size(y));
  
   k = find(sign(del(1:n-2)).*sign(del(2:n-1)) > 0);
  
   h = diff(x);
   hs = h(k)+h(k+1);
   w1 = (h(k)+hs)./(3*hs);
   w2 = (hs+h(k+1))./(3*hs);
   dmax = max(abs(del(k)), abs(del(k+1)));
   dmin = min(abs(del(k)), abs(del(k+1)));
   d(k+1) = dmin./conj(w1.*(del(k)./dmax) + w2.*(del(k+1)./dmax));

%  Slopes at end points.
%  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

   d(1) = ((2*h(1)+h(2))*del(1) - h(1)*del(2))/(h(1)+h(2));
   if sign(d(1)) ~= sign(del(1))
      d(1) = 0;
      %disp('d(1) = 0');
   elseif (sign(del(1)) ~= sign(del(2))) && (abs(d(1)) > abs(3*del(1)))
      d(1) = 3*del(1);
      %disp('d(1) = 3*del(1)');
   end
   d(n) = ((2*h(n-1)+h(n-2))*del(n-1) - h(n-1)*del(n-2))/(h(n-1)+h(n-2));
   if sign(d(n)) ~= sign(del(n-1))
      d(n) = 0;
      %disp('d(n) = 0');
   elseif (sign(del(n-1)) ~= sign(del(n-2))) && (abs(d(n)) > abs(3*del(n-1)))
      d(n) = 3*del(n-1);
      %disp('d(n) = 3*del(n-1)');
   end
%   whos
%   d