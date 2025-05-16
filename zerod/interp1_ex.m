% surcouche pour traiter le cas a 1 temps en meme temps
function rep = interp1_ex(t,y,tt,varargin)

% modification
if isappdata(0,'EXTRAPOLATION_NEAREST_EX')
    extrapolation_ex = getappdata(0,'EXTRAPOLATION_NEAREST_EX');
end
if isappdata(0,'NONAN_EX')
    nonan_ex = getappdata(0,'NONAN_EX');
end
if nonan_ex
    mask = any(~isfinite(y),2);
    if any(mask)
        t(mask)   = [];
        y(mask,:) = [];
    end
end

if (length(t) <= 1) && (length(tt) <= 1)
	rep = y;
elseif (length(t) == length(tt)) && all(t(:) == tt(:));
	rep = y;
elseif  (length(t) <= 1)
    rep = ones(length(tt),1) * y;
else
    rep = interp1(t,y,tt,varargin{:});
    if extrapolation_ex
        mask = tt < min(t);
        if any(mask)
            rep(mask,:) =y(1,:);
        end
        mask = tt > max(t);
        if any(mask)
            rep(mask,:) =y(end,:);
        end
    end
end
