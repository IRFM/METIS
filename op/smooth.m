function y=smooth(x,width,typ)

% SMOOTH	data smoothing gaussian or rectangular. 
%
% y=smooth(x,width)
% y=smooth(x,width,typ)
%
% y		smoothed data
%
% x		data to smooth. If x is a matrix, smooth each column.
% width		smoothing width
% typ		smoothing type : 'gauss'ian or 'rect'angular. default is 'gaus'
% EPFL-CRPP-TCA  J.-M. Moret (december 1988)

if nargin<3, typ='gaus'; end

if typ=='gaus'
	width=width*2;
	w=(-width:width)/width*2;
	w=exp(-w.^2);
else
	w=one(-width:width);
	end
w=w/sum(w);

n=size(x);
if n(1,1)==1, x=x'; end
for i=1:n(1,2)
	xf=conv(x(:,i),w);
	y(:,i)=xf(width+1:n(1,1)+width);
	end
