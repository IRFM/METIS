% inportalion au plus proche causale
function yy = zinterpnc(t,y,tt)

if size(t,1) == 1
	t = t';
	y = y';
end
if size(tt,1) == 1
	tt = tt';
	trans = 1;
else
	trans = 0;
end

tt = tt * ones(1,length(t));
t  = ones(length(tt),1) * t';
mask = (tt >= t);
ind_out = sum(mask,2); 
yy = y(ind_out);

if trans
	yy = yy';
end
