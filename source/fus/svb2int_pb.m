% calul de l'integrant pour sigmavnbitplasmad
function svel = svb2int_pb(v,vth,vb)


% constants
phys = cphys;

% cas  proton beam in boron
mu = 11.0093054 * 1.00782503207 / (11.0093054 + 1.00782503207) * phys.ua; % kg
ecm   = mu .* v .^ 2./ 2  ./ phys.e ; % eV
sdpb  =  pb11_cross_section_tentori(ecm); % m ^2
%%%% sdpb  =  pb11_cross_section_nevins(ecm); % m ^2
svel  = sdpb .* v .^ 2 .* (exp(-((v - vb) ./ vth) .^ 2) - exp(-((v + vb) ./ vth) .^ 2));



