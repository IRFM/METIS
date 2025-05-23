function P = zinterx(L1,varargin)
%INTERX Intersection of curves
%   P = INTERX(L1,L2) returns the intersection points of two curves L1 
%   and L2. The curves L1,L2 can be either closed or open and are described
%   by two-row-matrices, where each row contains its x- and y- coordinates.
%   The intersection of groups of curves (e.g. contour lines, multiply 
%   connected regions etc) can also be computed by separating them with a
%   column of NaNs as for example
%
%         L  = [x11 x12 x13 ... NaN x21 x22 x23 ...;
%               y11 y12 y13 ... NaN y21 y22 y23 ...]
%
%   P has the same structure as L1 and L2, and its rows correspond to the
%   x- and y- coordinates of the intersection points of L1 and L2. If no
%   intersections are found, the returned P is empty.
%
%   P = INTERX(L1) returns the self-intersection points of L1. To keep
%   the code simple, the points at which the curve is tangent to itself are
%   not included. P = INTERX(L1,L1) returns all the points of the curve 
%   together with any self-intersection points.
%   
%   Example:
%       t = linspace(0,2*pi);
%       r1 = sin(4*t)+2;  x1 = r1.*cos(t); y1 = r1.*sin(t);
%       r2 = sin(8*t)+2;  x2 = r2.*cos(t); y2 = r2.*sin(t);
%       P = zinterx([x1;y1],[x2;y2]);
%       plot(x1,y1,x2,y2,P(1,:),P(2,:),'ro')

%   Author : NS
%   Version: 1.0, 12/12/08

%   Two words about the algorithm: Most of the code is self-explanatory.
%   The only trick lies in the calculation of C1 and C2. To be brief, this
%   is essentially the two-dimensional analog of the condition that needs
%   to be satisfied by a function F(x) that has a zero in the interval
%   [a,b], namely
%           F(a)*F(b) <= 0
%   C1 and C2 exactly do this for each segment of curves 1 and 2
%   respectively. If this condition is satisfied simultaneously for two
%   segments then we know that they will cross at some point. 
%   Each factor of the 'C' arrays is essentially a matrix containing 
%   the numerators of the signed distances between points of one curve
%   and line segments of the other.

    %...Argument checks and assignment of L2
    error(nargchk(1,2,nargin));
    if nargin == 1,
        L2 = L1; 
        %hF = @lt;   %...Avoid the inclusion of common points
        mode_int = 0;
    else
        L2 = varargin{1}; 
        %hF = @le;
        mode_int = 1;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'
    S1 = repmat(dx1.*y1(1:end-1) - dy1.*x1(1:end-1),1,size(L2,2)-1);
    S2 = repmat(dx2.*y2(1:end-1) - dy2.*x2(1:end-1),size(L1,2)-1,1);
   
%      C1 = feval(hF,(kron(dx1,y2(1:end-1)) - kron(dy1,x2(1:end-1)) - S1).*...
%           (kron(dx1,y2(2:end))- kron(dy1,x2(2:end)) - S1),0);       
%      
%      C2 = feval(hF,(kron(dx2,y1(1:end-1)) - kron(dy2,x1(1:end-1)) - S2).*...
%           (kron(dx2,y1(2:end)) - kron(dy2,x1(2:end)) - S2),0);

    switch mode_int
    case 0
	C1 = ((kron(dx1,y2(1:end-1)) - kron(dy1,x2(1:end-1)) - S1) .* ...
	      (kron(dx1,y2(2:end))- kron(dy1,x2(2:end)) - S1)) < 0;       
    
	C2 = ((kron(dx2,y1(1:end-1)) - kron(dy2,x1(1:end-1)) - S2) .* ...
	      (kron(dx2,y1(2:end)) - kron(dy2,x1(2:end)) - S2)) < 0;
    otherwise
	C1 = ((kron(dx1,y2(1:end-1)) - kron(dy1,x2(1:end-1)) - S1) .* ...
	      (kron(dx1,y2(2:end))- kron(dy1,x2(2:end)) - S1)) <= 0;       
    
	C2 = ((kron(dx2,y1(1:end-1)) - kron(dy2,x1(1:end-1)) - S2) .* ...
	      (kron(dx2,y1(2:end)) - kron(dy2,x1(2:end)) - S2)) <= 0;
    end

    %...Obtain the points where an intersection is expected
    [i,j] = find(C1&C2);
    x2 = x2';dx2=dx2';y2=y2';dy2=dy2';  
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    %...Solve system of eqs to get the common points
    P = unique([(dx2(j).*(dx1(i).*y1(i) - dy1(i).*x1(i)) ...
         + dx1(i).*(dy2(j).*x2(j) - dx2(j).*y2(j))),...
           dy1(i).*(dy2(j).*x2(j) - dx2(j).*y2(j))...
         + dy2(j).*(dx1(i).*y1(i) - dy1(i).*x1(i))]./[L L],'rows')';


