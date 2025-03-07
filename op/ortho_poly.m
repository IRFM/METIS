function pl=ortho_poly(kf,x,n)

% This is a code downloaded from the website of MIT.
% http://ceta.mit.edu/comp_spec_func/

%       ==========================================================
%       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
%                or Ln(x) or Hn(x), and their derivatives
%       Input :  KF --- Function code
%                       KF=1 for Chebyshev polynomial (First kind) Tn(x)
%                       KF=2 for Chebyshev polynomial (Second kind) Un(x)
%                n ---  Order of orthogonal polynomials
%                x ---  Argument of orthogonal polynomials
%       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
%                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
%       =========================================================

% The only improvement in this program is it accepts vector arguments for x

% make sure that x is a row or column vector and not a matrix.
[r,c]=size(x);
if r==1 | c==1
    rowvec = 0;
    if r==1
        x=x';
        rowvec = 1;
    end
else
    error('x must be a vector, and cannot be a matrix');
end
lenx = length(x);

if n==0
    if rowvec
        pl = ones(1,lenx);
    else
       pl = ones(lenx,1);
    end
else
    pl = zeros(lenx,n);
    
    a=2;
    b=0;
    c=1;
    y0=1;
    y1=2.*x;
    
    % the i'th position in pl corresponds to the i'th term
    % don't bother storing pl = 1;
    
    pl(:,1)=2.*x;
    
    if (kf == 1)
        y1=x;
        pl(:,1)=y1;
    end
    
    for  k=2:n
        yn=(a.*x+b).*y1-c*y0;
        pl(:,k)=yn;
        y0=y1;
        y1=yn;
    end
    if rowvec
        pl = pl(:,n)';
    else
       pl = pl(:,n);
    end
end