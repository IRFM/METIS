% existance of path to second stability = no real root to equation 10 from
% paper : An analytical model of how the negative triangularity cuts off the access to the
% second stable region in tokamak plasmas, Y. Zhang et al, Nucl. Fusion 64 (2024) 026020 (9pp) 
% https://doi.org/10.1088/1741-4326/ad1b94
% return true if second stability domain is accessible.
% syntax: flag = z0secondstab(epsi,q,K,dmetis)
% with: epsi = inverse aspect ratio, q = safety factor, K = elongation,
% dmetis = geometrical triangularity as define in METIS.
function flag = z0secondstab(epsi,q,K,dmetis)


% v0 & v1 for domain mapping 
if numel(epsi) > 1
    v0 = zeros(size(epsi));
    v1 = ones(size(epsi));
elseif numel(q) > 1
    v0 = zeros(size(q));
    v1 = ones(size(q));
elseif numel(K) > 1
    v0 = zeros(size(K));
    v1 = ones(size(K));
elseif numel(dmetis) > 1
    v0 = zeros(size(dmetis));
    v1 = ones(size(dmetis));
else
    v0 = 0;
    v1 = 1;
end
% from METIS to paper definition for triangularity
d = - (acos(dmetis) - pi/2);

% A
A = 2 ./ (K .^ 2 + 1);

% DM without alpha
dm0 =  - epsi .* (1 -1 ./ q.^ 2 - 3 .* (K .^ 2 -1) ./ (3 .* K .^ 2 +1) ./ (K .^2 + 1) .* ...
                 (K .^ 2 - 2 .* d ./ epsi));
             
% coefficients of polynome
c0 = v0;
c1 = -dm0;
c2 = A .^ 2 ./ 4 + 1 / 4 - A ./ 2  - (1 - A) .* (1 / 2 - A / 2); 
c3 = v0;
c4 = 5 / 64 .* A .^ 2 - 3 / 32 .* (1 - A) .* A .^ 2;

% size
if (numel(c2) == 1) && (numel(v1) > 1)
    c2 = v1 .* c2;
end
if (numel(c4) == 1) && (numel(v1) > 1)
    c4 = v1 .* c4;
end

% roots
[x1,x2,x3,x4]= zpolyroot(c0,c1,c2,c3,c4);    

% roots selection
xsol = max(real(x1) .* double(~iscomplex(x1)) .* (real(x1) >= 0), ...
       max(real(x2) .* double(~iscomplex(x2)) .* (real(x2) >= 0),  ...
       max(real(x3) .* double(~iscomplex(x3)) .* (real(x3) >= 0), ...
       real(x4) .* double(~iscomplex(x4)) .* (real(x4) >= 0))));
   
% flag with tolerance
flag = (xsol < sqrt(eps));
   

