function z = xmtimes(x,y)

% XMTIMES(X,Y) computes inner product of x*y using fast libraries.

% for those who like slow computers
nx = size(x); ny = size(y);
if nx(end) ~= ny(1), error('Inner dimension mismatch'), end
z = reshape(reshape(x,prod(nx(1:end-1)),nx(end)) * reshape(y,ny(1),prod(ny(2:end))),[nx(1:end-1) ny(2:end)]);
