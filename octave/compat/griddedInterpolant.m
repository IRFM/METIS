function F = griddedInterpolant(X1, X2, V, method, extrap)
% Minimal Octave shim for the 2D ndgrid form used in METIS (z0icrh).
if nargin < 4, method = 'linear'; end
x1 = X1(:,1); x2 = X2(1,:);
F = @(q1, q2) interp2(x2, x1, V, q2, q1, method, NaN);
end
