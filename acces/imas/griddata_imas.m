% surcouche pour traiter le cas a 1 temps en meme temps
function rep = griddata_imas(t,x,y,tt,xx,methode)


if size(t,1) <= 1
	if all(size(x) == size(xx)) && all(x == xx)
		rep = y;
	else
		switch methode
		case 'cubic'
			methode = 'pchip';
		end
		rep  = interp1(x,y,xx,methode,'extrap');
	end
else
    try
        rep = griddata(t,x,y,tt,xx,methode);
    catch
        rep = griddata(t,x,y,tt,xx,methode,{'Qt','Qbb','Qc','Qz'});        
    end
end

rep(~isfinite(rep)) = -9.99999e99;

