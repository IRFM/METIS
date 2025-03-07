% test de zgriddata
x = rand( 30,7);
y = rand(30,7);
f = sin(2*pi.*x) .* cos(2*pi.*y);
w = ones(size(x));
s = 0.001;
[NX,NY,LAMDA,MU,C,FP,DEFRANK,FAIL] = zfitdata(x,y,f,w,s);
[F,FAIL] = zvalfit(NX,NY,LAMDA,MU,C,x,y)
%subplot(2,1,1);
plot3(x,y,f,'o',x,y,F,'+')

xx =cat(2,x,randn(30,7));
yy = cat(2,y,randn(30,7));

[NX,NY,LAMDA,MU,C,FP,DEFRANK,FAIL] = zfitdata(x,y,f,w,s,[min(xx(:)),max(xx(:)),min(yy(:)),max(yy(:))]);
[F,FAIL] = zvalfit(NX,NY,LAMDA,MU,C,xx,yy)
%subplot(2,1,2);
plot3(x,y,f,'o',xx,yy,F,'+')
