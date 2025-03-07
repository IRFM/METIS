% calul de l'integrant pour sigmavnbitplasmad
function svel = svb2int(v,vth,vb,mode)

if nargin < 4
	mode = 0; %DT
end

% constantes utiles
mp    = 1.6726485e-27;            % masse au repos du proton (kg)
e     = 1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
switch mode
case 0
	% cas dt
	mu    = mp .* 3 .* 2 ./ (3 + 2);  
	ecm   = mu .* v .^ 2./ 2  ./ 1.602176462e-19 ./ 1e3; %(m/s ->keV)
	svel  = z0csdt(ecm) .* v .^ 2 .* (exp(-((v - vb) ./ vth) .^ 2) - exp(-((v + vb) ./ vth) .^ 2));
case 2
	% cas tt
	mu    = mp .* 3 .* 3 ./ (3 + 3);  
	ecm   = mu .* v .^ 2./ 2  ./ 1.602176462e-19 ./ 1e3; %(m/s ->keV)
	stt   = tt_cross_section(ecm .* 1e3) ./ 1e-31;
	svel  = stt .* v .^ 2 .* (exp(-((v - vb) ./ vth) .^ 2) - exp(-((v + vb) ./ vth) .^ 2));
	disp('here')

otherwise
	% cas dd
	mu    = mp .* 2 .* 2 ./ (2 + 2);  
	ecm   = mu .* v .^ 2./ 2  ./ 1.602176462e-19 ./ 1e3; %(m/s ->keV)
	svel  = z0csdd(ecm) .* v .^ 2 .* (exp(-((v - vb) ./ vth) .^ 2) - exp(-((v + vb) ./ vth) .^ 2));
end



svel = svel .* 1e-3 .* 1e-24 .* 1e-4; % mb *m / s -> m^3/s
