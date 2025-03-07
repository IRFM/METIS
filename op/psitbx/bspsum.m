function v = bspsum(t,c,x,d,p,q)

% BSPSUM   B-spline base function combination
%   BSPSUM(T,C,X) estimates at each X the combination of spline base functions
%   defined by the knot sequence T. C contains along its first dimension the
%   coefficients of the combination, and BSPSUM returns an array of size
%   [size(X),size(C)(2:end)] or [length(X),size(C)(2:end)] if X is a vector.
%   BSPSUM(T,C,X,D) allows to estimate the D-th derivative.
%   BSPSUM(T,C,X,D,1) expects a C of size [size(C,1),size(X),DIM] and
%   returns an array of size [size(X),DIM], having applied the corresponding
%   coefficient to the corresponding X value.
%   BSPSUM(T,C,X,D,*,1) assumes distinct knots and uses a faster coding.
%   Spline order K is length(T)-size(C,1).
%   BSPSUM uses a MEX-file much faster than SPVAL in the Spline Toolbox.

% test:  plot(bspsum(-2:12,eye(12,15),linspace(0,10,101)))

if nargin < 4, d = 0; end
if nargin < 5, p = 0; end
if nargin < 6, q = all(diff(t)); end

nt = length(t);
dimc = size(c);
dimx = size(x); if length(dimx) == 2 & min(dimx) < 2, dimx = length(x); end
if nt - dimc(1) < 1
 error('Too many coefficients for specified knot sequence')
end
if p & (length(dimc)+1 < length(dimx) | ~isequal(dimc(2:length(dimx)+1),dimx))
 error('Incompatible C and X size')
end

v = bspsummex(t,c,prod(dimc(2:end)),x,prod(dimx),d,p,q);
if p
 v = reshape(v,[dimx,dimc(length(dimx)+2:end),1]);
else
 v = reshape(v,[dimx,dimc(2:end)]);
end
