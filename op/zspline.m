% comme spline mais supprin eles point dupliques
function yi = zspline(x,y,xi)

[x,ind]  = sort(x);
y        = y(ind);
ind = find(diff(x)<= 0);
x(ind) = [];
y(ind) = [];

yi = spline(x,y,xi);
