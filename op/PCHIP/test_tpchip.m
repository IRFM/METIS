% test mexfile
profile clear
profile on
for k =1:10000
	n = 2 + 100 .* rand(1);
	x = linspace(0,6.28.*rand(1),n);
	nn = fix(3 + 1000 .* rand(1));
	xx = linspace(randn(1),10.*rand(1),nn);
	y = randn(length(x),1);
%	clf
        yold = pchip(x,y,xx);
%	plot(x,y,'.r',xx,pchip(x,y,xx),'bo')
	[yy,d,fail] = tpchip(x,y,xx);
%  	fail
%  	hold on
%  	plot(xx,yy,'k');
%  	pause(1)
end
profile off
profile report

