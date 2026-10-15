function pp = pwch(x, y, s, dx, divdif)
% Octave shim for MATLAB's internal PWCH (piecewise cubic Hermite), as used by pchip.m.
n = numel(x); d = size(y, 1);
dxd = repmat(dx(:).', d, 1);
dzzdx  = (divdif - s(:, 1:n-1)) ./ dxd;
dzdxdx = (s(:, 2:n) - divdif) ./ dxd;
dnm1 = d * (n - 1);
coefs = [reshape((dzdxdx - dzzdx) ./ dxd, dnm1, 1), reshape(2*dzzdx - dzdxdx, dnm1, 1), ...
         reshape(s(:, 1:n-1), dnm1, 1), reshape(y(:, 1:n-1), dnm1, 1)];
pp = mkpp(x, coefs, d);
end
