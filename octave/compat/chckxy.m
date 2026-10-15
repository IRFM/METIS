function [x, y, sizey, endslopes] = chckxy(x, y)
% Octave shim for MATLAB's internal CHCKXY (data checks for pchip/spline).
n = numel(x); x = x(:).';
if isvector(y), y = y(:).'; end
sizey = size(y);
if sizey(end) ~= n
  if sizey(end) == n + 2
    endslopes = y(:, [1 end]); y = y(:, 2:end-1); sizey = size(y);
  else
    error('chckxy: the number of sites does not match the number of values');
  end
else
  endslopes = [];
end
sizey = sizey(1:end-1); if isempty(sizey), sizey = 1; end
y = reshape(y, prod(sizey), n);
if any(diff(x) < 0), [x, ind] = sort(x); y = y(:, ind); end
nany = any(isnan(y), 1);
if any(nany), x(nany) = []; y(:, nany) = []; end
if any(diff(x) == 0), error('chckxy: the data sites should be distinct'); end
end
