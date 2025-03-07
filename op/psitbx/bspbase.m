function b = bspbase(t,k,x)

% BSPBASE   B-spline base functions
%   BSPBASE(T,K,X) returns the B-spline base functions of order K (degree+1) for
%   the knot sequence T evaluated at X, one in each row. There are LENGTH(T)-K
%   base functions. X must be contained in the interval [T(K), T(END-K+1)].
%   BSPBASE uses a MEX-file much faster than SPCOL in the Spline Toolbox. See
%   also BSPSUM.

if length(t)-k <= 0, error('Not enough knots for specified order.'), end
x = x(:);
sortflag = any(diff(x) < 0);
if sortflag, [x,ksort] = sort(x); end
if x(1) < t(k) | x(end) > t(end-k+1)
 error('Values outside knot range.')
end
b = bspbasemex(t,k,x);
if sortflag, b(:,ksort) = b; end
