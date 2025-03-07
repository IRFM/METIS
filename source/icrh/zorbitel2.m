function [s0_broad,dr_tot,dr_banana,dr_potato] = zorbitel2(rho,ene0,amass,zcharge,s0_origin,equi,phys)
%% -----------------------------------------------------------------
%%  PURPOSE:
%%  TO BROADEN THE FAST ION SOURCE PROFILE ACCORDING TO ORBIT 
%%  WIDTH EFFECTS
%% -----------------------------------------------------------------
%%  INPUT:
%% -------
%%  RHO       = TOROIDAL FLUX COORDINATE                   (M)
%%  ENE0      = PARTICLE ENERGY                            (EV)
%%  AMASS     = PARTICLE MASS NUMBER                       (-)
%%  ZCHARGE   = PARTICLE CHARGE NUMBER                     (-)
%%  S0_ORIGIN = SOURCE PROFILE WITHOUT ORBIT WIDTH EFFECTS (ANY)
%%  EQUI      = EQUILIBRIUM STRUCTURE
%% -----------------------------------------------------------------
%% OUTPUT:
%% --------
%% S0_BROAD  = BROADENED SOURCE PROFILE (SAME UNIT AS S0_ORIGIN)
%% DR_TOT    = PROFILE OF CALCULATED/CHOSEN ORBIT WIDTH (M)
%% DR_BANANA = PROFILE OF BANANA ORBIT WIDTH            (M)
%% DR_POTATO = PROFILE OF POTATO ORBIT WIDTH            (M)
%% -----------------------------------------------------------------

%% EXTRACT VARIABLES FROM STRUCTURE "EQUI"
vpr   = equi.vpr;      % VOLUME ELEMENT PROFILE                        (M2)
q     = equi.q;        % SECURITY FACTOR PROFILE                       (-)
ftrap = equi.ftrap;    % PROFILE OF FRACTION OF TRAPPED PARTICLES      (-)
r     = equi.a;        % PROFILE OF MINOR RADIUS, (RMAX-RMIN)/2        (M)
R     = equi.raxe;     % PROFILE OF FLUX SURFACE CENTRE, (RMAX+RMIN)/2 (M)
bmod  = sqrt(equi.b2); % MAGNETIC FIELD PROFILE                        (T)

%% NEED THE "CLOSEST" FUNCTION TO BE FOUND
%  chartest = which('closest');
%  if length(chartest) < 1
%    addpath ~/Projet_Cronos/source/idn/nemo/
%  end

%% PHYSICS CONSTANTS
qel = phys.e;
wmp = phys.mp;

%% INVERSE ASPECT RATIO = r/R (PROFILE)
epsi = r./R;

%% PARTICLE MASS (KG) AND CHARGE (C) AND VELOCITY (M/S)
mass   = amass .* wmp;
charge = zcharge*qel;
vperp  = sqrt(2.*qel.*ene0 ./ mass);

%% INITIALIZATIONS FOR POTATO WIDTH CALCULATION
B0     = bmod(1);              % SCALAR
R0     = R(1)-R(end);          % SCALAR
q0     = q(1);                 % SCALAR
vperp0 = vperp.*sqrt(B0./bmod);% PROFILE

%% NORMALIZED TOROIDAL FLUX COORDINATE (-)
rhonorm = rho./rho(end);

%% ORBIT WIDTH (M)
dr_banana = (epsi).^(-0.5).*q.*vperp.*mass./(charge.*bmod);
dr_potato = (2.*q0.*vperp0./(charge.*B0./mass.*R0)).^(2./3.).*R0;
dr_tot    = min(dr_banana,dr_potato);

%% FORCE THE BANANA WIDTH TO BE USED AT THE EDGE
idxbreak             = find(dr_banana==min(dr_banana));
dr_tot(idxbreak:end) = dr_banana(idxbreak:end);

%% NORMALIZED ORBIT WIDTH (-)
dr_norm = dr_tot./rho(end);

%% VECTOR TO CREATE THE CONVOLUTION MATRIX
xv = ones(size(rho));

%% CONVOLUTION MATRIX: 1/SQRT(2PI*DR^2)*EXP(-X^2/DR^2)
dd = ((1./sqrt(2.*pi.*dr_norm'.^2))*xv)'.*exp(-((rhonorm'*xv-xv'*rhonorm).^2./(dr_norm'*xv)'.^2)); 

%% CONVOLUTION BETWEEN GAUSSIAN AND INITIAL SOURCE PROFILE
%% WITH A RENORMALIZATION IN ORDER TO CONSERVE THE NUMBER OF PARTICLES
norme    = xv' * trapz(rhonorm',dd,1);
s0_broad = trapz(rhonorm,dd./norme.*(xv'*s0_origin),2)';

%% DO NOT APPLY THE ORBIT WIDTH EFFECTS ON PASSING PARTICLES
s0_broad = ftrap.*s0_broad+(1-ftrap).*s0_origin; 

%% RENORMALIZATION SO THAT TO CONSERVE THE ENERGY
eout     = trapz(rhonorm,vpr.*s0_broad);
ein      = trapz(rhonorm,vpr.*s0_origin);
s0_broad = s0_broad .* ein ./ max(eps,eout);




