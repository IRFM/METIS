% from O. Sauter original paper + erratum
function eta = z0eta_spitzer(ne,te,zion)

nz  = 0.58 + 0.74 ./ (0.76 + zion);   % OK
out = 1.9012e4 .* te .^ (3/2) ./ zion ./ nz ./ lambda_ln_e(ne,te); % OK
eta = 1 ./ out;

function out = lambda_ln_e(ne,te)

out = 31.3 - log(sqrt(ne)./te); % OK

