% calul de l'integrant pour sigmavnbidplasmahe3 and sigmavnbihe3plasmahed
function svel = svb2int_dhe3(v,vth,vb)


% constants
phys = cphys;

% for pure D and He3
mu = 3.016029322 * 2.0141017778 / (3.016029322 + 2.0141017778) * phys.ua;
ecm   = mu .* v .^ 2./ 2  ./ phys.e  /  1e3 ;% keV
sddhe3  =  1e-31 .* z0csdhe3(ecm); % m^2
svel  = sddhe3 .* v .^ 2 .* (exp(-((v - vb) ./ vth) .^ 2) - exp(-((v + vb) ./ vth) .^ 2));



