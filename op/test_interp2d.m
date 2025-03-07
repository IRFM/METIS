% test de zgriddata
x = rand( 30,7);
y = rand(30,7);
f = sin(2*pi.*x) .* cos(2*pi.*y);
xx = cat(2,x,2.*rand(30,7));
yy = cat(2,y,2.*rand(30,7));
[ff,dffdx,dffdy,fail] = zinterp2d(x,y,f,xx,yy);
ind = find(ff >1e308);
ff(ind) = NaN .* ones(1,length(ind));
plot3(x,y,f,'o',xx,yy,ff,'+')