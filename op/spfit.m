function output = spfit(x,y,xb,n,xs,ys,js)
%SPFIT Fit a spline to noisy data.
%   PP = SPFIT(X,Y,XB) fits a piecewise cubic spline with break points
%   (or knots) XB to the noisy data (X,Y). Use PPVAL to evaluate PP.
%
%   PP = SPFIT(X,Y,XB,N) is a generalization to piecewise polynomial 
%   functions of order N (degree N-1), with continuous derivatives
%   up to order N-2. Default is a cubic spline with N = 4.
%
%       N     Function
%       ---------------------------------------------------------
%       2     Piecewise linear and continuous
%       3     Piecewise quadratic with continuous first derivative
%       4     Piecewice cubic with continuous second derivative (default)
%       5     Etc.
%
%   PP = SPFIT(X,Y,XB,N,XS,YS,JS) adds sharp data (exact conditions) to
%   the spline. The array YS specifies values or derivatives for each
%   location in XS. JS specifies the order of the derivative. If all the
%   data in YS is of order zero, JS can be omitted. Derivatives of order
%   N-1 or greater are rejected.
%
%   Example 1:
%       x = cumsum(rand(1,50));
%       y = sin(x/2) + cos(x/4) + 0.04*randn(size(x));
%       xb = x(1:3:end);
%       pp = spfit(x,y,xb);
%       xx = linspace(min(x),max(x),200);
%       yy = ppval(pp,xx);
%       yb = ppval(pp,xb);
%       plot(x,y,'bo',xx,yy,'r',xb,yb,'r.')
%
%   Example 2:
%       x = linspace(0,2*pi,75);
%       y = sin(x) + 0.07*randn(size(x));
%       xb = [0 2 4.2 6.2];
%       pp = spfit(x,y,xb,5);
%       xx = linspace(min(x),max(x),200);
%       yy = ppval(pp,xx);
%       yb = ppval(pp,xb);
%       plot(x,y,'bo',xx,yy,'r',xb,yb,'r.')
%
%   Example 3:
%	The exact conditions y(0) = 0, y'(0) = 1 and y'(pi) = -1 in Example 2
%   are set by
%       xs = [0, 0, pi];
%       ys = [0, 1, -1];
%       js = [0, 1, 1];
%       pp = spfit(x,y,xb,4,xs,ys,js);       
%
%   See also PPVAL, SPLINE

%   Author: jonas.lundgren@saabgroup.com, 2007.

%   2007-02-07  Code cleanup. 
%   2007-02-21  Generalization to piecewise polynomial splines of arbitrary order. 
%   2007-12-11  Exact conditions added. 
%   2008-09-25  Missing pre-allocation added.
%   2008-12-17  New polynomial base eliminating half the unknowns.
%               Short description of the numerical method added.
%   2009-02-24  Update of examples in help.

%--------------------------------------------------------------------------
%   THE NUMERICAL METHOD
%
%   PROBLEM: Minimize |A*u - a| subject to the constraint B*u = b.
%
%   Since A and B are sparse matrices the problem should be solved using
%   sparse Quadratic Programming. Until I know how to do I will use full
%   QR decomposition. What we need are any u satisfying B*u = b and a base
%   for the null space of B.
%
%   Compute [Q,R] = qr(B') and partition Q = [Q1 Q2] and R = [R1; R2],
%   where R1 is upper triangular, R2 = 0, B' = Q*R = Q1*R1 and B*Q2 = 0.
%   Q2 is a base for the null space of B. Use Q as a base for the solution
%   and set u = Q*v = Q1*v1 + Q2*v2 = u1 + u2, where ui = Qi*vi.
%   The constraint reads
%
%       b = B*u = B*Q1*v1 + B*Q2*v2 = R1'*v1 + 0*v2 = R1'*v1
%
%   We get v1 = R1'\b and u1 = Q1*(R1'\b). u1 is a particular solution of
%   the constraint equations. Now we must find v2 to minimize the norm of 
%
%       A*u - a = A*u1 + A*u2 - a = A*Q2*v2 - (a - A*u1) = G*v2 - g
%   
%   where G = A*Q2 and g = a - A*u1. We get v2 = G\g and u2 = Q2*(G\g).
%   The sparse matrix A is not assembled. The products A*Q2 and A*u1 are
%   computed sparsely.
%
%   NOTE: A straightforward but less efficient solution for u1 and Q2 is
%   u1 = B\b and Q2 = null(B).


if nargin < 1, help spfit, return, end
if nargin < 2, y = 1; end
if nargin < 3, xb = 0; end
if nargin < 4, n = 4; end
if nargin < 5, xs = []; end
if nargin < 6, ys = zeros(size(xs)); end
if nargin < 7, js = zeros(size(xs)); end


%   Check order

if n < 2
    msgid = 'SPFIT:BadOrder'; 
    message = 'The polynomial order N must be at least 2!';
    error(msgid,message)
end


%   Check noisy data

x = x(:);
y = y(:);

if length(y) ~= length(x)
    if length(y) == 1
        y = y*ones(size(x));
    else
        msgid = 'SPFIT:BadNoisyData';
        message = 'Vectors X and Y must have the same length!';
        error(msgid,message)
    end
end


%   Check sharp data

xs = xs(:);
ys = ys(:);
js = js(:);
ms = max([length(xs),length(ys),length(js)]);

msgid = 'SPFIT:BadSharpData';
message = 'Vectors XS, YS and JS must have the same length!';

if length(xs) ~= ms
    if length(xs) == 1
        xs = xs*ones(ms,1);
    else
        error(msgid,message)
    end
end
if length(ys) ~= ms
    if length(ys) == 1
        ys = ys*ones(ms,1);
    else
        error(msgid,message)
    end
end
if length(js) ~= ms
    if length(js) == 1
        js = js*ones(ms,1);
    else
        error(msgid,message)
    end
end

xs = xs(js < n-1);          % Reject derivatives of high order
ys = ys(js < n-1);
js = js(js < n-1);


%   Sort and check the break points

xb = unique(xb(:));
if length(xb) < 2
    xb = [min(x); max(x)];
    if xb(1) == xb(2)
        xb(2) = xb(1) + 1;
    end
end

hb = diff(xb);


%   Adjust limits

xlim = xb;
xlim(1) = -Inf;
xlim(end) = Inf;


%   Dimensions

n1 = floor(n/2);            % n = n1 + n2
n2 = ceil(n/2);

nb = length(xb);            % Number of break points

iu(1:2:nb) = n1;            % Number of unknowns per break point
iu(2:2:nb) = n2;
ju = [0, cumsum(iu)];       % Cumulative sum of unknowns
nu = ju(nb+1);              % Number of unknowns

mb = nu - nb - n + 2;       % Number of smoothness conditions
ms = length(xs);            % Number of sharp data points
m = mb + ms;                % Total number of exact conditions

mn = length(x);             % Number of noisy data points

if m == nu
    msgid = 'SPFIT:NoisyDataIgnored';
    message = 'Noisy data ignored. Sharp data determines the spline!';
    warning(msgid,message)
elseif m > nu
    msgid = 'SPFIT:TooMuchSharpData';
    message = 'Overdetermined system. Too much sharp data!';
    error(msgid,message)
end


%   Generate polynomial base

[P1,P2] = polybase(n);


%   Smoothness (exact conditions)

B = zeros(m,nu);
b = zeros(m,1);

jj = 1;
for k = 1:nb-2              % Loop over inner breaks
    ni = iu(k);
    kk = ju(k)+1;
    for j = n-ni:n-2
        p1 = mypolyval(P1,P2,hb(k),hb(k),j,ni);
        p2 = mypolyval(P1,P2,0,hb(k+1),j,n-ni);
        Brow = [p1, zeros(1,ni)] - [zeros(1,ni), p2];
        B(jj,kk:kk+n+ni-1) = Brow;
        jj = jj+1;
    end
end


%   Sharp data (exact conditions)

for i = 1:ms                % Loop over sharp data
    k = find(xs(i) >= xlim,1,'last');
    ni = iu(k);
    kk = ju(k)+1;
    Brow = mypolyval(P1,P2,xs(i)-xb(k),hb(k),js(i),ni);
    B(jj,kk:kk+n-1) = Brow;
    b(jj) = ys(i);
    jj = jj+1;
end


%   QR decomposition

if  m == 0
    u1 = zeros(nu,1);
    Q2 = eye(nu);
elseif m < nu
    [Q,R] = qr(B');
    Q1 = Q(:,1:m);
    R1 = R(1:m,:);
    u1 = Q1*(R1'\b);
    Q2 = Q(:,m+1:nu);
else % m == nu
    u1 = B\b;
    Q2 = zeros(nu,0);
end


%   Noisy data (weak conditions in the least square sense)

G = zeros(mn,nu-m);
g = zeros(mn,1);

if m < nu
    jj = 1;
    for k = 1:nb-1          % Loop over intervalls
        I = (x <= xlim(k+1)) & (x > xlim(k));
        xdata = x(I);
        ydata = y(I);
        d = length(xdata);
        if d > 0
            ni = iu(k);
            kk = ju(k)+1;
            Ablock = mypolyval(P1,P2,xdata-xb(k),hb(k),0,ni);
            G(jj:jj+d-1,:) = Ablock*Q2(kk:kk+n-1,:);
            g(jj:jj+d-1,1) = ydata - Ablock*u1(kk:kk+n-1);
            jj = jj+d;
        end
    end
end


%   Solve

u2 = Q2*(G\g);
u = u1 + u2;


%   Expand coefficients

coefs = zeros(nb-1,n);
for k = 1:nb-1              % Loop over intervalls
    ni = iu(k);
    kk = ju(k)+1;
    uk = u(kk:kk+n-1);
    for j = 0:n-1
        coefs(k,n-j) = mypolyval(P1,P2,0,hb(k),j,ni)*uk;
        uk = uk/(j+1);
    end
end


%   Make piecewise polynomial

output = mkpp(xb,coefs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P1,P2] = polybase(n)
%   POLYBASE computes the polynomial bases P1 and P2 of order N.
%   The rows of P1(:,:,1) are the polynoms of the base (see POLYVAL) and
%   P1(:,:,J+1) are the derivatives of order J. P1 is the base of odd
%   intervalls and P2 of even intervalls. If N is even the bases are
%   identical, so P2 is not computed.

n1 = floor(n/2);
n2 = ceil(n/2);
D = diag(n-1:-1:1,1);       % Derivation matrix

x0 = [zeros(n-1,1); 1];
x1 = ones(n,1);

Dx0 = zeros(n,n2);
Dx1 = zeros(n,n2);
for k = 1:n2
    Dx0(:,k) = x0;
    Dx1(:,k) = x1;
    x0 = D*x0;
    x1 = D*x1;
end

P1 = zeros(n,n,n);
P1(:,:,1) = inv([Dx0(:,1:n1), Dx1(:,1:n2)]);        % Polynomial base 1
for k = 2:n
    P1(:,:,k) = P1(:,:,k-1)*D;                      % Derivatives
end

if n2 > n1
    P2 = zeros(n,n,n);
    P2(:,:,1) = inv([Dx0(:,1:n2), Dx1(:,1:n1)]);    % Polinomial base 2
    for k = 2:n
        P2(:,:,k) = P2(:,:,k-1)*D;                  % Derivatives
    end
else
    P2 = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = mypolyval(P1,P2,x,h,jder,ni)
%   MYPOLYVAL evaluates the derivatives (of order JDER) of the polynoms
%   in the base P1 (or P2). X is a vector of relative locations in the
%   intervall (the left break point is zero), H is the length of the 
%   intervall and NI is the number of unknowns associated to the left
%   break point. NI determines whether P1 or P2 is the base.

n = size(P1,1);
m = length(x);
y = zeros(m,n);

n1 = floor(n/2);
if ni == n1
    P = P1(:,:,jder+1);
else
    P = P2(:,:,jder+1);
end

a0 = 1/h^jder;
a = a0;

if m == 1 && x == 0
    for k = 1:n
        y(k) = a*P(k,n);
        a = h*a;
        if k == ni, a = a0; end
    end
else
    t = x/h;
    for k = 1:n
        z = a*P(k,jder+1);
        for j = jder+2:n
            z = t.*z + a*P(k,j);
        end
        y(:,k) = z;
        a = h*a;
        if k == ni, a = a0; end
    end
end

