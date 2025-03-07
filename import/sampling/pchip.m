function v = pchip(x,y,xx)
%PCHIP  Piecewise Cubic Hermite Interpolating Polynomial.
%   PP = PCHIP(X,Y) provides the piecewise polynomial form of a certain
%   shape-preserving piecewise cubic Hermite interpolant, to the values
%   Y at the sites X, for later use with PPVAL and the spline utility UNMKPP.
%   X must be a vector.
%   If Y is a vector, then Y(j) is taken as the value to be matched at X(j), 
%   hence Y must be of the same length as X.
%   If Y is a matrix or ND array, then Y(:,...,:,j) is taken as the value to
%   be matched at X(j),  hence the last dimension of Y must equal length(X).
%   
%   YY = PCHIP(X,Y,XX) is the same as YY = PPVAL(PCHIP(X,Y),XX), thus
%   providing, in YY, the values of the interpolant at XX.
%
%   The PCHIP interpolating function, p(x), satisfies:
%   On each subinterval,  X(k) <= x <= X(k+1),  p(x) is the cubic Hermite
%      interpolant to the given values and certain slopes at the two endpoints.
%   Therefore, p(x) interpolates Y, i.e., p(X(j)) = Y(:,j), and
%       the first derivative, Dp(x), is continuous, but
%       D^2p(x) is probably not continuous; there may be jumps at the X(j).
%   The slopes at the X(j) are chosen in such a way that
%      p(x) is "shape preserving" and "respects monotonicity". This means that,
%   on intervals where the data is monotonic, so is p(x);
%   at points where the data have a local extremum, so does p(x).
%
% Comparing PCHIP with SPLINE:
%   The function s(x) supplied by SPLINE is constructed in exactly the same way,
%   except that the slopes at the X(j) are chosen differently, namely to make 
%   even D^2s(x) continuous. This has the following effects.
%   SPLINE is smoother, i.e., D^2s(x) is continuous.
%   SPLINE is more accurate if the data are values of a smooth function.
%   PCHIP has no overshoots and less oscillation if the data are not smooth.
%   PCHIP is less expensive to set up.
%   The two are equally expensive to evaluate.
%
%   Example:
%
%     x = -3:3;
%     y = [-1 -1 -1 0 1 1 1];
%     t = -3:.01:3;
%     plot(x,y,'o',t,[pchip(x,y,t); spline(x,y,t)])
%     legend('data','pchip','spline',4)
%
%   Class support for inputs x, y, xx:
%      float: double, single
%
%   See also INTERP1, SPLINE, PPVAL, UNMKPP.

% References:
%   F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
%   Interpolation", SIAM J. Numerical Analysis 17, 1980, 238-246.
%   David Kahaner, Cleve Moler and Stephen Nash, Numerical Methods
%   and Software, Prentice Hall, 1988.
%
%   Copyright 1984-2004 The MathWorks, Inc.

% Check that data are acceptable and, if not, try to adjust them appropriately
[x,y,sizey] = chckxy(x,y);
h = diff(x); m = prod(sizey);

%ensure that the interpolation points are real

if nargin==3 && any(~isreal(reshape(xx,numel(xx),1))) 
  error(message('MATLAB:pchip:ComplexInterpPts')) 
end

% Compute slopes

del = diff(y,1,2)./repmat(h,m,1);
slopes = zeros(size(y));
for r = 1:m
     if isreal(del) && isreal(x) && isreal(y) 
      slopes(r,:) = pchipslopest_mex(x,y(r,:),del(r,:));  
%        if ~all(slopes(r,:) == pchipslopes(x,y(r,:),del(r,:)))
%  	keyboard
%        end
     elseif isreal(del)
      slopes(r,:) = pchipslopes(x,y(r,:),del(r,:));
     else
      realslopes = pchipslopes(x,y(r,:),real(del(r,:)));   
      imagslopes = pchipslopes(x,y(r,:),imag(del(r,:)));
      slopes(r,:) = complex(realslopes, imagslopes);
     end
end

% Compute piecewise cubic Hermite interpolant to those values and slopes

v = pwch(x,y,slopes,h,del); v.dim = sizey;

if nargin == 3   % if values are wanted instead, provide them
   v = ppval(v,xx);
end

% -------------------------------------------------------

function d = pchipslopes(x,y,del)
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
   
function [x,y,sizey,endslopes] = chckxy(x,y)
%CHCKXY check and adjust input for SPLINE and PCHIP
%   [X,Y,SIZEY] = CHCKXY(X,Y) checks the data sites X and corresponding data
%   values Y, making certain that there are exactly as many sites as values,
%   that no two data sites are the same, removing any data points that involve 
%   NaNs, reordering the sites if necessary to ensure that X is a strictly
%   increasing row vector and reordering the data values correspondingly,
%   and reshaping Y if necessary to make sure that it is a matrix, with Y(:,j)
%   the data value corresponding to the data site X(j), and with SIZEY the
%   actual dimensions of the given values. 
%   This call to CHCKXY is suitable for PCHIP.
%
%   [X,Y,SIZEY,ENDSLOPES] = CHCKXY(X,Y) also considers the possibility that
%   there are two more data values than there are data sites.
%   If there are, then the first and the last data value are removed from Y
%   and returned separately as ENDSLOPES. Otherwise, an empty ENDSLOPES is
%   returned.  This call to CHCKXY is suitable for SPLINE.
%
%   See also PCHIP, SPLINE.

%   Copyright 1984-2011 The MathWorks, Inc.

% make sure X is a vector:
if length(find(size(x)>1))>1 
  error(message('MATLAB:chckxy:XNotVector')) 
end

% ensure X is real
if any(~isreal(x)) 
  error(message('MATLAB:chckxy:XComplex')) 
end

% deal with NaN's among the sites:
nanx = find(isnan(x));
if ~isempty(nanx)
   x(nanx) = [];
   warning(message('MATLAB:chckxy:nan'))
end

n=length(x);
if n<2 
  error(message('MATLAB:chckxy:NotEnoughPts')) 
end

% re-sort, if needed, to ensure strictly increasing site sequence:
x=x(:).'; 
dx = diff(x);

if any(dx<0), [x,ind] = sort(x); dx = diff(x); else ind=1:n; end

if ~all(dx), error(message('MATLAB:chckxy:RepeatedSites')), end

% if Y is ND, reshape it to a matrix by combining all dimensions but the last:
sizey = size(y);


while length(sizey)>2&&sizey(end)==1, sizey(end) = []; end


yn = sizey(end); 
sizey(end)=[]; 
yd = prod(sizey);

if length(sizey)>1
   y = reshape(y,yd,yn);
else
   % if Y happens to be a column matrix, change it to the expected row matrix.
   if yn==1
       yn = yd;
       y = reshape(y,1,yn); 
       yd = 1; 
       sizey = yd;
   end
end

% determine whether not-a-knot or clamped end conditions are to be used:
nstart = n+length(nanx);
if yn==nstart
   endslopes = [];
elseif nargout==4&&yn==nstart+2
   endslopes = y(:,[1 n+2]); y(:,[1 n+2])=[];
   if any(isnan(endslopes))
      error(message('MATLAB:chckxy:EndslopeNaN'))
   end
   if any(isinf(endslopes))
       error(message('MATLAB:chckxy:EndslopeInf'))
   end
else
   error(message('MATLAB:chckxy:NumSitesMismatchValues',nstart, yn))
end

% deal with NaN's among the values:
if ~isempty(nanx)
    y(:,nanx) = [];
end

y=y(:,ind);
nany = find(sum(isnan(y),1));
if ~isempty(nany)
   y(:,nany) = []; x(nany) = [];
   warning(message('MATLAB:chckxy:IgnoreNaN'))
   n = length(x);
   if n<2 
     error(message('MATLAB:chckxy:NotEnoughPts')) 
   end
end

   