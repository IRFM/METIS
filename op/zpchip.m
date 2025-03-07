% comme spline mais supprin eles point dupliques
function yi = zpchip(x,y,xi)

[x,ind]  = sort(x);
y        = y(ind);
ind = find(diff(x)<= 0);
x(ind) = [];
y(ind) = [];

yi = pchip(x,y,xi);
