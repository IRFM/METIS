% fraction de l'energie deposee sur les ion
% Wesson, Tokamak p 227 2ieme edition
function fr = zfract0(ec,eb)

ec = max(eps,ec);
x   = max(1e-3,eb ./ ec);
fr  = ((1/3) .* log((1 - sqrt(x) +x) ./ (1 + sqrt(x)).^ 2) + ...
      2 ./ sqrt(3) .* (atan((2 .* sqrt(x) - 1) ./ sqrt(3) ) + pi/6)) ./ x; 

