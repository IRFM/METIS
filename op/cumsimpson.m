function z=cumsimpson(x,y)
%CUMSIMPSON Cumulative Simpson Numerical Integration.
% Z=CUMSIMPSON(Y) computes an approximation of the cumulative integral of Y
% using a Simpson rule with unit spacing between the data points in Y. To
% compute the integral for spacing different from one, multiply Z by the
% spacing increment. When Y is a matrix, the cumulative integral is
% computed over each column of Y.
%
% Z=CUMSIMPSON(X,Y) computes the integral with respect to the data in
% vector X. X need NOT be EQUALLY spaced but must have the same number of
% elements as Y. When Y is a matrix, Y must have as many rows as X has
% elements.
% 
% Z has the same dimensions as Y.
%
% See also CUMSUM, CUMTRAPZ, QUAD, QUADV.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2005-09-30, 2006-02-21

if nargin<2
   y=x;
   [ry,cy]=size(y);
   if ry==1
      x=1:cy;
   else
      x=1:ry;
   end
elseif nargin==2
   [ry,cy]=size(y);
else
   error('One or Two Inputs Are Required.')
end
if ndims(y)~=2
   error('N-dimensional Data is Not Supported.')
end
if min(size(x))>1
   error('X Must be a Vector.')
end
x=x(:);   % make x a column
if ry==1 % y is a row vector, make it a column
   yisrow=true;
   y=y.';
   ry=cy;
   cy=1;
else
   yisrow=false;
end
if length(x)~=ry
   error('Length of X Must Match Length of Y or Rows of Matrix Y.')
end
if length(x)<3
   error('At Least 3 Data Points are Required.')
end
dx=repmat(diff(x),1,cy);
dy=diff(y);

dx1=dx(1:end-1,:);
dx2=dx(2:end,:);
dxs=dx1+dx2;
dy1=dy(1:end-1,:);
dy2=dy(2:end,:);

a=(dy2./(dx2.*dxs) - dy1./(dx1.*dxs))/3;
b=(dy2.*dx1./(dx2.*dxs) + dy1.*dx2./(dx1.*dxs))/2;
c=y(2:end-1,:);

i1=((a.*dx1-b).*dx1+c).*dx1; % left half integral
i2=((a.*dx2+b).*dx2+c).*dx2; % right half integral

z=zeros(ry,cy);              % pure cumlative Simpson
z(2:2:end-1,:)=i1(1:2:end,:);
z(3:2:end,:)=i2(1:2:end,:);
z(end,:)=i2(end,:);
z=cumsum(z);

% z=[zeros(1,cy);...           % original algorithm here has ~2X error
%    cumsum([i1(1,:); (i1(2:end,:)+i2(1:end-1,:))/2; i2(end,:)])];

if yisrow
   z=z.';
end