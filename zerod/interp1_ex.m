% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_ex(t,y,tt,varargin)


if (length(t) <= 1) && (length(tt) <= 1)
	rep = y;
elseif (length(t) == length(tt)) && all(t(:) == tt(:));
	rep = y;
elseif  (length(t) <= 1)
    rep = ones(length(tt),1) * y;
    %whos
    %figure(21);clf;plot(y);drawnow
else
	rep = interp1(t,y,tt,varargin{:});
end
