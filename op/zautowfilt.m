function y_out = zautowfilt(y,t)

if all(y == 0)
	y_out = y;
	return
end

if nargin == 2
	tt = (min(t):min(diff(t)):max(t))';
	y  = pchip(t,y,tt);
end

lev=min(11,2 .* ceil(log(length(y))/log(2)/2)+1)
y_out = sgolayfilt(y,1,lev);

if nargin == 2
	y_out  = pchip(tt,y_out,t);
end


