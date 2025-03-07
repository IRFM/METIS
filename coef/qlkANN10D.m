%%%%%%%%  ARTIFICIAL NEURAL NETWORK BASED ON QUALIKIZ REGRESSION
%%%%%%%   ITG ETG & TEM. 10D INPUTS
%%%%%%%   APPROX 6 ORDERS OF MAGNITUDE FASTER THAN QUALIKIZ
%%%%%%%%  JONATHAN CITRIN, JUAN REDONDO, SARAH BRETON, CLARISSE BOURDELLE, FREDERIC IMBEAUX, 
%%%%%%%%% Karel van de Plassche , Jean-Francois Artaud
%%%%%%%%% TRANSPORT MODEL DEVELOPED FOR USE IN RAPTOR (F. FELICI et al) and METIS (J.F. Artaud et al)
%
% This function is the main interface between METIS and QLK-NN 10D which transporm input in units of METIS in input in units of QUALIKIZ
% and output with QUALIKIZ normalisation to data in physical units for METIS 

function [Chie,Chii,De,Ve] = qlkANN10D(xli,te,ti,ne,w_ExB,zeff,q,factor_em,meff,R0,a0,B0,Amain,qlknn_params)

%INPUTS
%xli - Lao coordinate from METIS
%te - electron temperature in eV 
%ti - ion temperature in eV
%ne - electron density in m^-3
%w_ExB - shear of radial electric field (same definition in METIS and QUALIKIZ)
%zeff  - effective charge
%q - safety factor
%factor_em - factor of reduction of L_Ti see by  ITG due to fast particles  
%meff  - effective mass
%R0 - device major radius [m]
%a0 - device minor radius (Rmin-Rmax)/2 [m]
%B0 - vacuum B0 at R0 [T]
%Amain - Atomic number of main ion species
%qlknn_params - QLK-NN 10D parameters (see declare_QLKNN10D_parameters)

%OUTPUTS
%chi(e,i) - transport coefficient output. SI units: m^2/s
%De - electron diffusivity. SI units: m^2/s
%Ve - electron pinch. SI units: m/s

  
%%%%%%BEGINNING OF CODE BODY
% constante physique (phys)
phys.c           =   2.99792458e8;             % vitesse de la lumiere dans le vide (m/s)  (definition)
phys.h           =   6.62606876e-34;           % constante de Planck (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % charge de l'electron (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeabilite du vide (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivite du vide (F/m)  (definition)
phys.g           =   6.673e-11;                % constante de la gravitation (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % constante de Boltzmann (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % constante de structure fine (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % masse au repos de l'electron (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % masse au repos du proton (kg)
phys.ua          =   1.66053873e-27;           % 1 unite atomique (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % nombre d'avogadro (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % constante de stephan ( W*m^-2*K^-4) (+/- 0.000040e-8)

% path to NN data
path_to_NN_data = fullfile(fileparts(which('install_QLKNN10D')),'QLKNN10D','qlknn-hyper-namelists');
qlknn_network.namelist_dir = path_to_NN_data;
qlknn_network.name = 'QLKNN10D';

% set qlkparam if not defined
param = declare_QLKNN10D_parameters;
if (nargin < 14) || isempty(qlknn_params)
    qlknn_params  = param.valeur;
end
if isfield(qlknn_params,'Effect_TioTe')
   Effect_TioTe = qlknn_params.Effect_TioTe;
   qlknn_params = rmfield(qlknn_params,'Effect_TioTe');
else
   Effect_TioTe = 'real value';
end
% cast 0 and 1 in logical
noms = fieldnames(qlknn_params);
for k=1:length(noms)
    if iscell(param.borne.(noms{k})) && (length(param.borne.(noms{k})) == 2)
        switch qlknn_params.(noms{k})
            case {1,'on','yes'}
                qlknn_params.(noms{k}) = true;
            otherwise
                qlknn_params.(noms{k}) = false;
        end
    elseif ischar(qlknn_params.(noms{k}))
        try
            qlknn_params.(noms{k}) = str2num(qlknn_params.(noms{k}));
        end
    end
end
if isfield(qlknn_params, 'norms')
    norms = qlknn_params.norms;
    qlknn_params = rmfield(qlknn_params, 'norms');
else
    norms = [];
    if qlknn_params.apply_victor_rule
        error('qlknn_params.norms needs to be defined if applying victor rule')
    end
end


%Rescale temperatures to SI units
teSI  = te*phys.e; 
tiSI  = ti*phys.e;
% compute gradient length
r     = a0(:,1) * xli;
ve = ones(size(xli));
vt = ones(size(te,1),1);
% Local gradient (for maximum allows value)
ted1  = pdederive(xli,teSI,0,2,2,1) ./ a0;
tid1  = pdederive(xli,tiSI,0,2,2,1) ./ a0;
ned1  = pdederive(xli,ne,0,2,2,1) ./ a0;
rlti  = -tid1  .* R0 ./ tiSI;
rlte  = -ted1  .* R0 ./ teSI;
rlne  = -ned1  .* R0./ ne;
% 
% q is q
%
tite = tiSI ./ max(eps,teSI);
switch Effect_TioTe
    case 'forced to 2.5'
        tite = 2.5 .* ones(size(tite));
    case 'forced to 0.25'
        tite = 0.25 .* ones(size(tite));
end
%
% sQLK = r/q dq/dr
shearr = pdederive(xli,q,0,2,2,1) .* (vt * xli) ./ q;
% 
% nustar
ne_19  = ne / 10^19;
te_keV = te / 10^3;
Lam = 15.2 - 0.5 .* log(0.1 .* ne_19) + log(te_keV);
nu = 914.7 .* zeff .* ne_19 .* Lam .* te_keV.^-1.5;
bounce = q .* R0 ./ ((r ./ R0) .^ 1.5 .* sqrt(10^3 .* te_keV .* phys.e ./ phys.me));
nuestar = nu .* bounce .* qlknn_params.factor_nustar;
log10nuestar = log10(nuestar);
%
% GammaE_GB
% 
cref = sqrt(phys.e * 10^3 / phys.mp);
%vtpr = (Vtor * Bpol - Vtheta * Bphi) / Btot @ LFS midplane
%vtpr ~ Er_metis * grad_rho ./ Btot @ LFS midplane
% same definition of w_ExB in METIS and Qualikiz
gamE = w_ExB .* R0 / cref;

% Victor rule normalization parameters
if qlknn_params.apply_victor_rule
    qlknn_params.norms.A1 = Amain;
    qlknn_params.norms.a  = a0;
    qlknn_params.norms.R0 = R0;
end

% stabilisation by fast particles
if qlknn_params.em_stab == 1
  rlti = rlti .* factor_em;
end

%% Assemble inputs
% Order into QLKNN: [Zeff, Ati, Ate, An, q, smag, x, Ti_Te, ln(Nustar),gammaE_GB, Te]
if isfield(qlknn_params, 'force_q')
    q = qlknn_params.force_q * ones(size(q));
    dq_dx = zeros(size(dq_dx));
end
if isfield(qlknn_params, 'force_s')
    shearr = qlknn_params.force_s * ones(size(q));
end

% definition of xin
epsilonNN = 1/3; % epsilon used to train 10D NN (see documentation)
xin = ((a0(:,1)./ R0(:,1)) / epsilonNN) * xli; % normalized x

% concanetation inputs
% QLKNN inputs:
% R/LTi = -dTi/dr * R0/Ti
% R/LTe = -dTe/dr * R0/Te
% R/Lne = -dne/dr * R0/ne
% Ti/Te
% q = q-profile
% s = magnetic shear
% x = r / a   (where a is the minor radius and r is the local coordinate in the midplane averaged radius coordinate r = (Rmax-Rmin)/2 , where Rmax and Rmin correspond to the local flux surface and =r_LCFS)
% nustar (defined in PDF)
% Zeff (defined in PDF)
% gammaE_GB =abs( dvtor/dr) *epsilon/q / (c_sou/a)    where c_sou=sqrt(Te/mi) and epsilon = x * a/R0 
%
% Order into QLKNN: [Zeff, Ati, Ate, An, q, smag, x, Ti_Te, ln(Nustar), gammaE_GB or Machtor]
%
% RAPTOR NOTE: in RAPTOR it is: if merge_modes = true nout = 7, namely :  
% 1 = qe, 2 = qi, 3 = Gamma_{e},  4 = De, 5 = Ve, 6 = Di, 7 = Vi
% If merge_modes = false:
% 1 = 'efeETG_GB', 2 = 'efeITG_GB', 3 ='efeTEM_GB', 4 ='efiITG_GB', 5 ='efiTEM_GB'
% 6 = 'pfeITG_GB', 7 = 'pfeTEM_GB', 8 = 'dfeITG_GB', 9 = 'dfeTEM_GB',
% 10 = 'vteITG_GB', 11 = 'vteTEM_GB', 12 = 'vceITG_GB', 13 = 'vceTEM_GB',
% 14 = 'dfiITG_GB', 15 = 'dfiTEM_GB', 16 = 'vtiITG_GB', 17 = 'vtiTEM_GB',
% 18 = 'vciITG_GB', 19 = 'vciTEM_GB', 20 = 'gam_leq_GB',
% Input is not hard-coded, but currently:
% 1 = Zeff, 2 = RLTi, 3 = RLTe, 4 = RLn, 5 = q, 6 = s^\hat{s}, 7 = X (QLK-NN definition !), 8 = Ti/Te, 9 = log_10(nue_star), 10 = gamma_ExB, 11 = T (keV) is only needed if use_victor_rule = true

insize = size(q);
NNin = [zeff(:),rlti(:),rlte(:),rlne(:),q(:),shearr(:),xin(:),tite(:),log10nuestar(:),gamE(:),te(:)/1e3];

% Assign outputs
if qlknn_params.apply_victor_rule
    NNout = qlknn_mex(qlknn_network.namelist_dir, NNin, qlknn_params.verbosity, qlknn_params, norms);
else
    NNout = qlknn_mex(qlknn_network.namelist_dir, NNin, qlknn_params.verbosity, qlknn_params);
end

% un applied normalization to ouput
%% GB normalization for outputs
%% taken from RAPTOR interface

%Define gyroBohm normalization. NN outputs are in GB normalized units
%chiGB = sqrt(mi)/(a*e^2*B^2)*Te^1.5
sqrtmi = sqrt(Amain* phys.mp);
teSI12 = sqrt(teSI);
teSI32 = teSI12.*teSI;
chifac = (sqrtmi(:,1)./(phys.e.^2.*B0(:,1).^2)./a0(:,1)) * ve;
chiGB  = teSI32.*chifac; % for chie,chii,de

%% Specific conversion factors per output (from documentation table)
irlte = 1./rlte;
irlti = 1./rlti;
irlne = 1./rlne;
irlte(:,1) = 0;
irlti(:,1) = 0;
irlne(:,1) = 0;

kchie   = (R0./a0) .* chiGB .* irlte;
kchii   = (R0./a0) .* chiGB .* irlti;
kDeff   = (R0./a0) .* chiGB .* irlne;
kD      = chiGB;
kV      = chiGB ./ a0;

% set ouput
if qlknn_params.merge_modes

  Chie = reshape(NNout(:,1),insize) .* kchie;
  Chii = reshape(NNout(:,2),insize) .* kchii;
  De   = reshape(NNout(:,4),insize) .* kD;
  Ve   = reshape(NNout(:,5),insize) .* kV;

else

  Chie = (reshape(NNout(:,1),insize).*qlknn_params.factor_ETG+reshape(NNout(:,2),insize)+reshape(NNout(:,3),insize)) .* kchie;
  Chii = (reshape(NNout(:,4),insize)+reshape(NNout(:,5),insize)) .* kchii;
  De   = (reshape(NNout(:,8),insize)+reshape(NNout(:,9),insize)) .* kD;
  Ve   = reshape(sum(NNout(:,10:13),2),insize) .* kV;

end

