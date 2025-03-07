% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_imas(t,y,tt,varargin)


if length(t) <= 1
	rep = y;
elseif (length(t) == length(tt)) && all(t == tt)
	rep = y;
elseif any(~isfinite(y(:)))
	rep = interp1(t,y,tt,'nearest','extrap');
else
	rep = interp1(t,y,tt,varargin{:});
end
