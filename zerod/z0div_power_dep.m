% ref PSI2012  T. Eich et al.
function [sb,qpdep,fnorm] = z0div_power_dep(angle,S,flux_exp,dsol,pref,Rdiv,sb_lim)

% s en m
fx    = flux_exp ./ sin(angle ./ 180 .* pi);
%  scale_hi = (4 + S .^ 2./ 4 ./ dsol .^ 2) .* dsol .* fx;
%  scale_hi = sum(pref .* scale_hi) ./ max(1,sum(pref))
%  scale_lo = (S ./ 2 ./ dsol - 3) .* S .* fx; 
%  scale_lo = sum(pref .* scale_lo) ./ max(1,sum(pref))
sb=linspace(-sb_lim / 2,sb_lim,301);
% sb_length=max(sb)-min(sb);
vt = ones(size(pref));
ve = ones(size(sb));
if length(pref) > 1
	qpl = pref * ve;
else
	qpl = pref;
end
if length(angle) > 1
	angle = angle * ve;
end
if length(S) > 1
	S  = S * ve;
end
if length(flux_exp) > 1
	flux_exp = flux_exp * ve;
end
if length(dsol) > 1
	dsol = dsol * ve;
end
if length(Rdiv) > 1
	Rdiv = Rdiv * ve;
end
Rdiv = abs(Rdiv);

% forme du depot
qpdep = qpl ./ 2  .* exp((S ./ (2 .* dsol)) .^ 2 - (vt * sb) ./ (dsol .* fx)) .* erfc((S ./ (2 .* dsol)) - (vt * sb) ./ (S .* fx));

% normalization (external leg of divertor)
% wetted_area_o = sol_width*flux_exp/sin(tilt_angle*pi/180)*2*pi*rso
pdep    = max(eps,trapz(sb,qpdep .* 2 .* pi .* (Rdiv + vt * sb .* sin(angle ./ 180 .* pi)),2));
qpdep   = ((pref ./ pdep) * ve)  .* qpdep;
fnorm   = pref ./ pdep; 