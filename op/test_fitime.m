% test de zgriddata
x = rand( 600,30);
y = rand(600,30);
f = sin(2*pi.*x) .* cos(2*pi.*y);
w = ones(size(x));
s = 1;
xx =cat(2,x,randn(600,30));
yy = cat(2,y,randn(600,30));

[NX,NY,LAMDA,MU,C,FP,DEFRANK,FAIL] = zfittime(x,y,f,w,s,400,20,[min(xx(:)),max(xx(:)),min(yy(:)),max(yy(:))]);
[F,FAIL] = zvalfit(NX,NY,LAMDA,MU,C,xx,yy)
%subplot(2,1,2);
plot3(x,y,f,'o',xx,yy,F,'+')
