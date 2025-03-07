function [Rsepa_out,Zsepa_out] = z0reshape_sepa(geo,zerod,Rsepa,Zsepa,expo,x,nbp)


ve     = ones(size(x));
vt     = ones(size(geo.R,1),1);
if size(x,1) == 1
	x = x *vt;
end

th   = linspace(0,2.*pi,nbp);
vh   = ones(size(th));
sh   = sin(th);
ch   = cos(th);
sina = vt * sh;
cosa = vt * ch;
t      = asin(geo.d .* x .^ 3); 

% separatrice donnees par les moments
Rmom   = geo.R * vh + (geo.a * vh) .* cos(vt * th + t * sin(th));
Zmom   = (geo.a .* geo.K) * sin(th);	
%
Raxe   = geo.R + zerod.d0 .* (1 - x .^ 2);
kx     = max(1,geo.K  .* x - (1 - x));
Rp = Raxe * vh + (geo.a * vh) .* (x * vh) .*cos( vt * th + (t * vh) .* sin(vt * th));
Zp = (kx * vh) .* (geo.a * vh) .* (x * vh) .* sin(vt * th);
[Rsepa_out,Zsepa_out]   = z0morph(Rp,Zp,Rsepa,Zsepa,Rmom,Zmom,x*vh,expo,1);

