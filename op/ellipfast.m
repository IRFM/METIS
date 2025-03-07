function f = ellipfast(m,kind,scalar)

% ELLIPFAST	Fast computation of the complete elliptic integral
%   ELLIPFAST(M,KIND) computes the complete elliptic integral of the first kind
%   K(M) and second kind E(M) for 0 <= M <= 1 with an error of less than 2e-8.
%   KIND is either 1 or 2 or 'K' or 'E'. Don't confuse the modulus K with M : M
%   = K^2.

% Source: Abramowitz 17.3.34 (page 591)
% fonction ecrite par J-M Moret (EPFL/CRPP)
global flag_out
if isempty(flag_out) 
    flag_out =0;
    flag_out = validellipfast;
    if flag_out == 0
       error('unvalid return of function ellipfast; check the compilation ...');
    end
end


if ischar(kind)
       kind = 13.5 - double(kind)/6; 
end
f = ellipfastmex(double(m),kind);

% la // n'est pas efficace pour le moment !
%  if (numel(m) > 1000) && (nargin < 3) 
%    try
%        warning off
%        nbcore = maxNumCompThreads;
%        warning on
%    catch
%        nbcore = matlabpool('size');
%    end
%    m   = double(m);
%    ss  = size(m);
%    m   = m(:);
%    nb  = fix(length(m)/nbcore);
%    nl  = fix(length(m)/nb);
%    fi  = zeros(nl,nb);
%    parfor (k=1:nl)
%    %for k=1:nl
%          fi(k,:)  = par_ellipfast(m,nb,k,kind);
%    end
%    fi = fi';
%    in = max((nl-1) * nb + (1:nb));
%    if in < length(m)
%       ind = (in + 1):length(m);
%       mi  = m(ind);
%       fc  = ellipfastmex(mi,kind);
%       f   = cat(1,fi(:),fc(:));
%    end
%    f = reshape(f,ss);
%    %disp(max(abs(f(:) - ellipfastmex(double(m),kind))));
%  
%  else
%    f = ellipfastmex(double(m),kind);
%  end

function f = par_ellipfast(m,nb,k,kind)

f  = ellipfastmex(m((k-1) * nb + (1:nb)),kind);


