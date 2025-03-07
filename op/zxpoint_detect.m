% cette fonction sert a detecter si le plasma est en mode divertor pour le module de bord.
function xpoint = zxpoint_detect(geo,equi)


%recupere le frontiere
switch geo.mode
case 3
	Rsepa = double(squeeze(equi.R(1,end,:)));
	Zsepa = double(squeeze(equi.Z(1,end,:)));
case 2	
	Rsepa = double(geo.R);
	Zsepa = double(geo.Z);
case 1
	th  = asin(max(0,min(1,geo.trh1)));
	tl  = asin(max(0,min(1,geo.trb1)));
	u  = linspace(0,2 .* pi,35);
	su = sin(u);
	Rsepa  = geo.r0 + geo.a  .* cos(u + (th .* (su >= 0) + tl .* (su < 0)) .* su);
	Zsepa  = geo.z0 + (geo.a .* geo.e1) * su;
case 0
	t  = asin(max(0,min(1,(geo.trh1 + geo.trb1) ./ 2)));
	u  = linspace(0,2 .* pi,35);
	Rsepa  = geo.r0 + geo.a  .* cos(u + t .* sin(u));
	Zsepa  = geo.z0 + (geo.a .* geo.e1) * sin(u);

end

if all(isfinite(Rsepa)) & all(isfinite(Zsepa))

	zmax     = max(Zsepa);
	zmin     = min(Zsepa); 
	indh     = find(Zsepa == zmax,1);
	indh     = cat(2,indh - 3, indh - 2,indh - 1,indh,indh + 1,indh + 2,indh + 3); 
	indh     = mod(indh-1,size(Zsepa,2))+1;
	rh       = Rsepa(indh);
	zh       = Zsepa(indh);
	warning off
	ph       = polyfit(rh,zh,2);
	warning on
	if (max(rh) - min(rh)) > 0
		eh       = sqrt(mean((zh - polyval(ph,rh)).^ 2)) ./ (max(rh) - min(rh));
	else
		eh = 0;
	end
	indl     = find(Zsepa == zmin,1);
	indl     = cat(2,indl - 3, indl - 2,indl - 1,indl,indl + 1,indl + 2,indl + 3);
	indl     = mod(indl-1,size(Zsepa,2))+1;
	rl       = Rsepa(indl);
	zl       = Zsepa(indl);
	warning off
	pl       = polyfit(rl,zl,2);
	warning on
	if (max(rl) - min(rl)) > 0
		el       = sqrt(mean((zl - polyval(pl,rl)).^ 2)) ./ (max(rl) - min(rl));
	else	
		el = 0;
	end
	xpoint   = max(el > 2e-2 , eh > 2e-2);

else
	xpoint = NaN;
end