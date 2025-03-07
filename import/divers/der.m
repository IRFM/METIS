function y=der(x,n);

% y=DER(x,n)
%
% DER calculates the n'th difference of x :  y=dx/dt  (with dt=1)
% using a third order scheme
%	x	matrix to be derived (data sets are columns)
%	n	derivative order n=1,2 or 3 (default is 1)
%	y	derived matrix

%							ThDdW  5/89

[m,n2]=size(x);
if m==1, x=x.'; m=n2; n2=1; end
if nargin<2, n=1; end			% default is first derivative

if n==1,
	edge = 3*x([1,m],:) - 4*x([2,m-1],:) + x([3,m-2],:);    % x' at edge
	y = [-edge(1,:);  x(3:m,:)-x(1:m-2,:);  edge(2,:)]/2;
elseif n==2,
	edge = 2*x([1,m],:) - 5*x([2,m-1],:) + 4*x([3,m-2],:) ...
			    - x([4,m-3],:);   		   	% x' at edge
	y = [-edge(1,:);  x(3:m,:)-2*x(2:m-1,:)+x(1:m-2,:);  edge(2,:)];
elseif n==3,
	edge =   x([1:2,m-1:m],:) - 3*x([2:3,m-2:m-1],:) ...
	     + 3*x([3:4,m-3:m-2],:) - x([4:5,m-4:m-3],:);
	y = [x(5:m,:)-2*x(4:m-1,:)+2*x(2:m-3,:)-x(1:m-4,:)];
	y = [-edge(1:2,:); y/2; edge(3:4,:)]; 
end
	
