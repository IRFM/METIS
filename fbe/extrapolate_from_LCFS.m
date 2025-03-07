% extrapolate on each raduius of center-LCFS using Psi as:
% psi = psi0 + psi_R * R + psi_Z *Z  + psi_RR * R ^2 + psi_ZZ * Z^2 + psi_RZ * R *Z
% assuming Psi, BR and BZ given on LCFS and
% G-S hold in vacuum and psi(R=0,Z) = 0 and psi_0 = 0.
%
function [FPSI,FBR,FBZ] = extrapolate_from_LCFS(R_LCFS,Z_LCFS,PSI_LCFS,BR_LCFS,BZ_LCFS,mode)


if nargin < 6
   mode = 2;
end

% formatting
R_LCFS      = R_LCFS(:);
Z_LCFS      = Z_LCFS(:);
PSI_LCFS    = PSI_LCFS(:);
BR_LCFS     = BR_LCFS(:);
BZ_LCFS     = BZ_LCFS(:);

% compute z0
mask = double(R_LCFS == max(R_LCFS));
z0   = sum(mask .* Z_LCFS) ./ max(1,sum(mask));
r0   = (max(R_LCFS) + min(R_LCFS)) / 2;


% make polar irregula grid
nbtheta = length(R_LCFS);
nbrho   = 301;
RR      = NaN * ones(nbrho,nbtheta); 
ZZ      = NaN * ones(nbrho,nbtheta); 
BR      = NaN * ones(nbrho,nbtheta); 
BZ      = NaN * ones(nbrho,nbtheta); 
PSI     = NaN * ones(nbrho,nbtheta); 

% loop on LCFS points
l = linspace(0,1,nbrho)';
% indic = [];
% figure;
% plot(R_LCFS,Z_LCFS,'k',R_LCFS,Z_LCFS,'.k');
% hold on
for k=1:nbtheta
    % make matrix for linear problem
    switch mode
        case 0
            matl = cat(1, ...
                cat(2,R_LCFS(k),   R_LCFS(k) .^ 2,     0,  R_LCFS(k) .* (Z_LCFS(k) - z0)), ...
                cat(2,0,   0,    - (Z_LCFS(k) - z0),  -R_LCFS(k)), ...
                cat(2,1,   2 .* R_LCFS(k),     0,  (Z_LCFS(k) - z0)), ...
                cat(2,-1,  -R_LCFS(k),  R_LCFS(k), (Z_LCFS(k)-z0)));
            data = cat(1,0,R_LCFS(k) * BR_LCFS(k),R_LCFS(k) * BZ_LCFS(k),0);
            rep  = lscov(matl,data);
            
            RR(:,k) = R_LCFS(k) + (R_LCFS(k) - r0) .* l;
            ZZ(:,k) = Z_LCFS(k) + (Z_LCFS(k) - z0) .* l;
            PSI(:,k) = PSI_LCFS(k) + rep(1) * RR(:,k) + rep(2) * RR(:,k) .^ 2 + ...
                rep(4) * RR(:,k)  .* (ZZ(:,k) - z0);
            BR(:,k) = -1 ./ R_LCFS(k) .* (rep(3) .* (ZZ(:,k) - z0) + rep(4) .* RR(:,k));
            BZ(:,k) = 1 ./ R_LCFS(k)  .* (rep(1) + 2 .* rep(2) .* RR(:,k) + rep(4) .* (ZZ(:,k) - z0));
        case 1
            % removed condition psi(R=0,Z) = 0
            matl = cat(1, ...
                   cat(2, R_LCFS(k),   R_LCFS(k) .^ 2, (Z_LCFS(k) - z0),  (Z_LCFS(k) - z0) .^ 2  ,  R_LCFS(k) .* (Z_LCFS(k) - z0)), ...
                   cat(2, 0,   0, -1,  - 2 .* (Z_LCFS(k) - z0)   ,  -R_LCFS(k) ), ...
                   cat(2, 1,   2 .* R_LCFS(k), 0,  0 ,  (Z_LCFS(k) - z0)), ...
                   cat(2, -1,   -R_LCFS(k) , 0,  R_LCFS(k)  ,  -R_LCFS(k) .* (Z_LCFS(k) - z0)));
            data = cat(1,0,R_LCFS(k) * BR_LCFS(k),R_LCFS(k) * BZ_LCFS(k));
            
            warning off
            rep  = lscov(matl,data);
            warning on
            
            RR(:,k) = R_LCFS(k) + (R_LCFS(k) - r0) .* l;
            ZZ(:,k) = Z_LCFS(k) + (Z_LCFS(k) - z0) .* l;
            PSI(:,k) = PSI_LCFS(k) +  ...
                       rep(1) * RR(:,k) + rep(2) * RR(:,k) .^ 2 + ...
                       rep(3) * (Z_LCFS(k) - z0) + rep(4) * (Z_LCFS(k) - z0) .^ 2 + ...
                       rep(5) * RR(:,k)  .* (ZZ(:,k) - z0);
            BR(:,k) =  -1 ./ RR(:,k) .* ( ...
                        rep(3) + 2.* rep(4) * (Z_LCFS(k) - z0)+ ...
                        rep(5) * RR(:,k));
            BZ(:,k) =  1 ./ RR(:,k) .* ( ...
                       rep(1)  +  2.* rep(2) * RR(:,k) + ...
                       rep(5) .* (ZZ(:,k) - z0));
        case 2
                        
            % start from true analytical polynomial solution
            [psi_coef,BR_coef,BZ_coef] = psivide_coef(R_LCFS(k),(Z_LCFS(k) - z0));
            matl = cat(1, psi_coef,BR_coef,BZ_coef);
            data = cat(1,  0,BR_LCFS(k),BZ_LCFS(k) - mean(BZ_LCFS)); 
          
            % the term for vertical magnetic field (2) is constant and not
            % useful is the fit. Keep least orders terms 
            ind_sel = [ 1   3   8 ];          
                       
            warning off
            rep  = lscov(matl(:,ind_sel),data);           
            warning on
            
            
            RR(:,k) = R_LCFS(k) + (R_LCFS(k) - r0) .* l;
            ZZ(:,k) = Z_LCFS(k) + (Z_LCFS(k) - z0) .* l;
            [psi_coef,BR_coef,BZ_coef] = psivide_coef(RR(:,k),(ZZ(:,k) - z0));
            PSI(:,k) = PSI_LCFS(k) + psi_coef(:,ind_sel) * rep + 0.5 .* mean(BZ_LCFS) .* (RR(:,k) .^ 2 - R_LCFS(k) .^ 2); % adding contribution of term 2 for vertical magentic field contribution.
            BR(:,k)  = BR_coef(:,ind_sel) * rep;
            BZ(:,k)  = mean(BZ_LCFS) + BZ_coef(:,ind_sel) * rep;
              
    end
end
% plot(RR,ZZ,'.');
% quiver(RR(1:10:end,:),ZZ(1:10:end,:),BR(1:10:end,:),BZ(1:10:end,:));
% contour(RR,ZZ,PSI,101)
%keyboard
% 
% create interpolant function
warning off
FPSI   = scatteredInterpolant(RR(:),ZZ(:),PSI(:),'natural','linear');
FBR    = scatteredInterpolant(RR(:),ZZ(:),BR(:),'natural','linear');
FBZ    = scatteredInterpolant(RR(:),ZZ(:),BZ(:),'natural','linear');
warning on




function [psi_coef,BR_coef,BZ_coef] = psivide_coef(R,Z)

% calcul (partie contribution du vide)
% ref :
% Analytical tokamak equilibrium for shaped plasmas
% S. B. Zheng,a) A. J. Wootton, and Emilia R. Solano
% Phys. Plasmas 3 (3), March 1996

% A.J. Cerfon and J.P. Freidberg, POP 17, 032502 (2010)
% choix R0 = 1

% keep only to order 2 in this case whith the constant on psi
% we have only 3 (or 4) constains 

% 1
psi_1 = 1.* ones(size(R));
BR_1  = 0.* ones(size(R));
BZ_1  = 0.* ones(size(R));

% 2
psi_2 = R .^ 2;
BR_2  = 0 .* ones(size(R));
BZ_2  = 2 .* ones(size(R));

% 3
psi_3 = Z .^ 2  - R .^ 2 .* log(R);
BR_3  = - 2 .* Z ./ R;
BZ_3  = - 2 .* log(R) - 1;

%4
psi_4 = R .^ 4 - 4 .* R .^ 2 .* Z .^ 2;
BR_4  = 8 .* R .* Z;
BZ_4  = 4 .* R .^ 2 - 8 .* Z .^ 2;

%5
psi_5 = 2 .* Z .^ 4  - 9 .* R .^ 2 .* Z .^ 2  + 3 .* R .^ 4  .* log(R)  - 12 .* R .^ 2 .* Z .^ 2 .* log(R);
BR_5  = - (8 .* Z .^ 3 + (-24 .* R .^ 2 .* log(R) - 18 .* R .^ 2) .* Z) ./ R;
BZ_5  =  -(24 .* log(R) + 30) .* Z .^ 2 + 12 .* R .^ 2 .* log(R) + 3 .* R .^ 2;

%6
psi_6 = R .^ 6  - 12 .* Z .^ 2 .* R .^ 4 + 8 .* Z .^ 4 .* R .^ 2;
BR_6  = -(32 .* R .* Z .^ 3 - 24 .* R .^ 3 .* Z);
BZ_6  = 16 .* Z .^ 4 - 48 .* R .^ 2 .* Z .^ 2 + 6 .* R .^ 4;


psi_7 = 8 .* Z .^ 6  - 140 .* Z .^ 4 .* R .^ 2 + 75 .* Z .^ 2 .* R .^ 4  - 15  .* R .^ 6 .* log(R) + 180 .* R .^ 4 .* Z .^ 2 .* log(R) - 120 .* R .^ 2 .* Z .^ 4 .* log(R);
BR_7  = -(48 * Z .^ 5 - 560 .* Z .^ 3  .* R .^ 2 + 150 .* Z  .* R .^ 4 + 360 .* R .^ 4 .* Z .* log(R) - 480 .* R .^ 2 .* Z .^ 3 .* log(R))./ R;
BZ_7  = - 280 .* Z .^ 4 + 300 .* Z .^ 2 .* R .^ 2 - 90  .* R .^ 4 .* log(R) - 15  .* R .^ 4 + 720 .* R .^ 2 .* Z .^ 2 .* log(R) + 180 .* R .^ 2 .* Z .^ 2 - 240  .* Z .^ 4 .* log(R) - 120  .* Z .^ 4;

% 8
psi_8 = Z;
BR_8  = - 1./ R;
BZ_8  = 0 .* ones(size(R));

% 9
psi_9 = Z .* R .^ 2;
BR_9  = - R;
BZ_9  = 2.* Z;

% 10
psi_10 = Z .^ 3 - 3 .* Z .* R .^ 2 .* log(R);
BR_10  = -(3 .* Z .^ 2 - 3 .* R .^ 2 .* log(R)) ./ R;
BZ_10  = (-6 .* log(R) - 3) .* Z;

% 11
psi_11 = 3 .* Z .* R .^ 4  - 4 .* Z .^ 3 .* R .^ 2;
BR_11  = - 3 .* R .^ 3 + 12 .* R .* Z .^ 2;
BZ_11 = 12 .* R .^ 2 .* Z - 8 .* Z .^ 3;

%12
psi_12 = 8 .* Z .^ 5 - 45 .* Z .* R .^ 4 - 80 .* Z .^ 3 .* R .^ 2 .* log(R) + 60 .* Z .* R .^ 4 .* log(R);
BR_12  = - (40 .* Z .^ 4 - 240 .* R .^ 2 .* log(R) .* Z .^ 2 + 60 .* R .^ 4 .* log(R) - 45 .* R .^ 4) ./ R;
BZ_12  = (-160 .* log(R) - 80) .* Z .^ 3 + (240 .* R .^ 2 .* log(R) - 120 .* R .^ 2) .* Z;

psi_coef  = cat(2,psi_1,psi_2,psi_3,psi_4,psi_5,psi_6,psi_7,psi_8,psi_9,psi_10,psi_11,psi_12);
BR_coef   = cat(2,BR_1,BR_2,BR_3,BR_4,BR_5,BR_6,BR_7,BR_8,BR_9,BR_10,BR_11,BR_12);
BZ_coef   = cat(2,BZ_1,BZ_2,BZ_3,BZ_4,BZ_5,BZ_6,BZ_7,BZ_8,BZ_9,BZ_10,BZ_11,BZ_12);
% 


