%%%%%%%%  ARTIFICIAL NEURAL NETWORK BASED ON QUALIKIZ REGRESSION
%%%%%%%   ITG ETG & TEM. 10D INPUTS
%%%%%%%   APPROX 6 ORDERS OF MAGNITUDE FASTER THAN QUALIKIZ
%%%%%%%%  JONATHAN CITRIN, JUAN REDONDO, SARAH BRETON, CLARISSE BOURDELLE, FREDERIC IMBEAUX, 
%%%%%%%%% Karel van de Plassche , Jean-Francois Artaud
%%%%%%%%% TRANSPORT MODEL DEVELOPED FOR USE IN RAPTOR (F. FELICI et al) and METIS (J.F. Artaud et al)
%
% This function is the main interface between METIS and QLK-NN 10D which transporm input in units of METIS in input in units of QUALIKIZ
% and output with QUALIKIZ normalisation to data in physical units for METIS 

% function [Chie,Chii,De,Ve] = qlkz_std(xli,te,ti,ne,w_ExB,zeff,q,factor_em,meff,R0,a0,B0,Amain,qlkparams,prof,option,zerod)
function varargout  = qlkz_std(xli,te,ti,ne,w_ExB,zeff,q,factor_em,meff,R0,a0,B0,Amain,qlkparams,prof,option,zerod)

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
%qlkparams - QLK-NN 10D parameters (see declare_QLKNN10D_parameters)

%OUTPUTS
%chi(e,i) - transport coefficient output. SI units: m^2/s
%De - electron diffusivity. SI units: m^2/s
%Ve - electron pinch. SI units: m/s

if ~isfield(option,'Sn_fraction')
    Sn_fraction = 0;
else
    Sn_fraction = option.Sn_fraction;
end



  
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
path_to_METIS_QLKZ = fullfile(fileparts(which('install_QLKNN10D')));
qlknn_network.namelist_dir = path_to_METIS_QLKZ;
qlknn_network.name = 'QLKNN10D';

% set qlkparam if not defined
param = declare_QLKZ_std_parameters;
if isempty(qlkparams)
    qlkparams  = param.valeur;
end
if isfield(qlkparams,'Effect_TioTe')
   Effect_TioTe = qlkparams.Effect_TioTe;
   qlkparams = rmfield(qlkparams,'Effect_TioTe');
else
   Effect_TioTe = 'real value';
end

qlpar = qlkparams;

% cast 0 and 1 in logical
noms = fieldnames(qlkparams);
for k=1:length(noms)
    if iscell(param.borne.(noms{k})) && (length(param.borne.(noms{k})) == 2)
        switch qlkparams.(noms{k})
            case {1,'on','yes'}
                qlkparams.(noms{k}) = true;
            otherwise
                qlkparams.(noms{k}) = false;
        end
    elseif ischar(qlkparams.(noms{k}))
        try
            qlkparams.(noms{k}) = str2num(qlkparams.(noms{k}));
        end
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
irho = min(abs(xli-qlkparams.domain_TiovTe));

switch Effect_TioTe
case 'forced to 0.25'
      disp(['Te/Ti modified inside ' num2str(xli(irho))])
      for hh=1:size(tite,1)
        tite(hh,1:irho) = 0.25;% .* ones(size(tite));
      end
case 'forced to 1'
      disp(['Te/Ti modified inside ' num2str(xli(irho))])
      for hh=1:size(tite,1)
        tite(hh,1:irho) = 1;% .* ones(size(tite));
      end
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
nuestar = nu .* bounce .* qlkparams.coll_mult;
log10nuestar = log10(nuestar);
%
% GammaE_GB
% 
cref = sqrt(phys.e * 10^3 / phys.mp);
%vtpr = (Vtor * Bpol - Vtheta * Bphi) / Btot @ LFS midplane
%vtpr ~ Er_metis * grad_rho ./ Btot @ LFS midplane
% same definition of w_ExB in METIS and Qualikiz
gamE = w_ExB .* R0 / cref;

%gamE = -w_ExB .* R0 / cref;

% stabilisation by fast particles
if qlkparams.em_stab == 1
  rlti = rlti .* factor_em;
end

insize = size(q);

NNin_str.zeff = zeff;
NNin_str.rlti = rlti;
NNin_str.rlte = rlte;
NNin_str.rlne = rlne;
NNin_str.q = q;
NNin_str.shearr = shearr;
NNin_str.tite = tite;
NNin_str.te = te./1e3;
NNin_str.ne = ne./1e19;

NNin_str.Raxe = prof.Raxe;
NNin_str.R0 = repmat(prof.Raxe(:,end),1,size(te,2));
NNin_str.epsi = prof.epsi;
NNin_str.rmin = prof.epsi.*prof.Raxe;
NNin_str.amin = repmat(prof.epsi(:,end).*prof.Raxe(:,end),1,size(te,2));

NNin_str.B0 = B0;

%IMPURITIES PART
% Three species probably stored in .mat
z1 = option.zimp;  %for this one we specify a Z and a profile
z2 = option.zmax;  %for this one we specify a Z and take the profile of z1
z3 = round((1 - Sn_fraction) .* z0wavez(NNin_str.te*1000) + Sn_fraction .* z0snavez(NNin_str.te*1000));  %Tungsten

%The specified Z
NNin_str.z1 = z1.*ones(size(te,1),size(te,2));  
NNin_str.z2 = z2.*ones(size(te,1),size(te,2));
NNin_str.z3 = z3;

%The specified profiles
NNin_str.nz1p = prof.nzp./1e19;
NNin_str.nz2p = prof.nzp.*option.rimp./1e19;  %rimp = ratio of impurities
NNin_str.nwp = prof.nwp./1e19;
NNin_str.nip = (NNin_str.ne - z1.*NNin_str.nz1p-z2.*NNin_str.nz2p-z3.*NNin_str.nwp)./(1+option.cmin);  %Ambipolarity of profils
NNin_str.nminp = NNin_str.nip.*option.cmin;  %Hydrogen profil
N = [NNin_str.nip,NNin_str.nz1p,NNin_str.nz2p,NNin_str.nwp,NNin_str.nminp,NNin_str.ne];
%Gradients
nz1d1  = pdederive(prof.xli,prof.nzp,0,2,2,1) ./ a0;
rlnz1  = -nz1d1  .* R0./ prof.nzp;
rlnz2  = rlnz1;

nz3d1  = pdederive(prof.xli,prof.nwp,0,2,2,1) ./ a0;
rlnz3  = -nz3d1  .* R0./ prof.nwp;

%4 impurities z1,z2,z3 (W) and cmin (H)
NNin_str.rlnz1 = rlnz1;
NNin_str.rlnz2 = rlnz2;
NNin_str.rlnz3 = rlnz3;
NNin_str.cmin = option.cmin.*ones(size(te,1),size(te,2));   %cmin = 1/(n_H/(n_H  + n_D)-1)
NNin_str.rova = NNin_str.rmin./NNin_str.amin;

% Put a warning when difference on one of QuaLiKiz input is too large between all 4 times considered
fftmp = fieldnames(NNin_str);
time_flag = 0;
for ii=1:length(fieldnames(NNin_str));
  eval(['dif_eva = diff(NNin_str.' fftmp{ii} './mean(NNin_str.' fftmp{ii} '));']);
  if any(dif_eva(:)>0.001)
    time_flag = 1;
  end
end

if time_flag
  for ii=1:length(fieldnames(NNin_str));
     eval(['NNin_str.' fftmp{ii} '=reshape(NNin_str.' fftmp{ii} ',1,size(te,1).*size(te,2));']);
  end
else
 for ii=1:length(fieldnames(NNin_str));
  eval(['NNin_str.' fftmp{ii} '=NNin_str.' fftmp{ii} '(3,:);']);
 end
end

% Uses to prepare QuaLiKiz standalone simulations
username = getenv('USER');
runname_qlk=['METIS_QLKZ_tmp_' username];
runpath_qlk = [path_to_METIS_QLKZ '/runs_qlkz/from_METIS/'];

tmp_dir = tempdir;
%runpath_qlk = [tmp_dir 'runs_qlkz/from_METIS/'];

switch qlpar.QLKZ_mode
    case 'NN-like'
      disp('=========== Runs QuaLiKiz wih only D, e and Carbon ============')
      build_qualikiz_input_from_METIS_QLKZ_NNlike(NNin_str,runname_qlk,path_to_METIS_QLKZ,qlkparams);
    case 'Full'
      disp('=========== Runs QuaLiKiz wih all METIS species ============')
      build_qualikiz_input_from_METIS_QLKZ(NNin_str,runname_qlk,path_to_METIS_QLKZ,qlkparams,Sn_fraction);
end

% Read back outputs
yout_qlkz = struct_qlkz(runname_qlk,tmp_dir,0);

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

%Computation of heat transport coefficients, used in z0qlkz_std
if ~time_flag
  NNout(:,1) = reshape(repmat(yout_qlkz.efe_GB,1,size(te,1)),1,size(te,1).*size(te,2));
  NNout(:,2) = reshape(repmat(sum(yout_qlkz.efi_GB,2),1,size(te,1)),1,size(te,1).*size(te,2));

  Chie = reshape(NNout(:,1),flip(insize))'.* kchie;
  Chii = reshape(NNout(:,2),flip(insize))'.* kchii;

else

  NNout(:,1) = yout_qlkz.efe_GB;
  NNout(:,2) = sum(yout_qlkz.efi_GB,2);

  Chie = reshape(NNout(:,1),insize).* kchie;
  Chii = reshape(NNout(:,2),insize).* kchii;
end

% Species transport coefficients extraction
if qlkparams.imp_trans ~= 0

  if ~time_flag

    for i = 1:size(yout_qlkz.dfi_GB,2)
      %For all ions species
      NNoutI(:,1,i) = reshape(repmat(yout_qlkz.dfi_SI(:,i),1,size(te,1)),1,size(te,1).*size(te,2));
      NNoutI(:,2,i) = reshape(repmat(yout_qlkz.vti_SI(:,i)+yout_qlkz.vri_SI(:,i)+yout_qlkz.vci_SI(:,i),1,size(te,1)),1,size(te,1).*size(te,2));
      NNoutI(:,3,i) = reshape(repmat(yout_qlkz.pfi_SI(:,i),1,size(te,1)),1,size(te,1).*size(te,2));
      
      D(:,:,i) = reshape(NNoutI(:,1,i),flip(insize))';
      V(:,:,i) = reshape(NNoutI(:,2,i),flip(insize))';
      G(:,:,i) = reshape(NNoutI(:,3,i),flip(insize))';

      %For electrons
      % NNout(:,3) = reshape(repmat(yout_qlkz.dfe_GB,1,size(te,1)),1,size(te,1).*size(te,2));
      % NNout(:,4) = reshape(repmat(yout_qlkz.vce_GB+yout_qlkz.vte_GB,1,size(te,1)),1,size(te,1).*size(te,2));
      % NNout(:,5) = reshape(repmat(yout_qlkz.pfe_GB,1,size(te,1)),1,size(te,1).*size(te,2));

      % D(:,:,size(yout_qlkz.dfi_GB,2)+1) = reshape(NNout(:,3),flip(insize))' .* kD;
      % V(:,:,size(yout_qlkz.dfi_GB,2)+1) = reshape(NNout(:,4),flip(insize))' .* kV;
      % G(:,:,size(yout_qlkz.dfi_GB,2)+1) = reshape(NNout(:,5),flip(insize))' .* kV .* N(:,size(yout_qlkz.dfi_GB,2)+1);

      %For electrons
      NNout(:,3) = reshape(repmat(yout_qlkz.dfe_SI,1,size(te,1)),1,size(te,1).*size(te,2));
      NNout(:,4) = reshape(repmat(yout_qlkz.vce_SI+yout_qlkz.vte_SI,1,size(te,1)),1,size(te,1).*size(te,2));
      NNout(:,5) = reshape(repmat(yout_qlkz.pfe_SI,1,size(te,1)),1,size(te,1).*size(te,2));

      D(:,:,size(yout_qlkz.dfi_SI,2)+1) = reshape(NNout(:,3),flip(insize))';
      V(:,:,size(yout_qlkz.dfi_SI,2)+1) = reshape(NNout(:,4),flip(insize))';
      G(:,:,size(yout_qlkz.dfi_SI,2)+1) = reshape(NNout(:,5),flip(insize))';
      
    end
 
  else

    for i = 1:size(yout_qlkz.dfi_GB,2)
      %For all ions species
      NNoutI(:,1,i) = yout_qlkz.dfi_SI(:,i);
      NNoutI(:,2,i) = yout_qlkz.vti_SI(:,i)+yout_qlkz.vri_SI(:,i)+yout_qlkz.vci_SI(:,i);
      NNoutI(:,3,i) = yout_qlkz.pfi_SI(:,i);

      D(:,:,i) = reshape(NNoutI(:,1,i),insize) ;
      V(:,:,i) = reshape(NNoutI(:,2,i),insize) ;
      G(:,:,i) = reshape(NNoutI(:,3,i),insize) ;

      %For electrons
      NNout(:,3) = yout_qlkz.dfe_GB;
      NNout(:,4) = yout_qlkz.vce_GB+yout_qlkz.vte_GB;
      NNout(:,5) = yout_qlkz.pfe_GB;

      D(:,:,size(yout_qlkz.dfi_GB,2)+1)   = reshape(NNout(:,3),insize) .* kD;
      V(:,:,size(yout_qlkz.dfi_GB,2)+1)   = reshape(NNout(:,4),insize) .* kV;
      G(:,:,size(yout_qlkz.dfi_GB,2)+1)   = reshape(NNout(:,5),insize) .* kV .* N(:,size(yout_qlkz.dfi_GB,2)+1);

      %For electrons
      NNout(:,3) = yout_qlkz.dfe_SI;
      NNout(:,4) = yout_qlkz.vce_SI+yout_qlkz.vte_SI;
      NNout(:,5) = yout_qlkz.pfe_SI;

      D(:,:,size(yout_qlkz.dfi_SI,2)+1)   = reshape(NNout(:,3),insize);
      V(:,:,size(yout_qlkz.dfi_SI,2)+1)   = reshape(NNout(:,4),insize);
      G(:,:,size(yout_qlkz.dfi_SI,2)+1)   = reshape(NNout(:,5),insize);
    end

  end 

  varargout{1} = Chie;
  varargout{2} = Chii;
  varargout{3} = D;
  varargout{4} = V;
  varargout{5} = G;

else %without transport coefficient extraction
  if ~time_flag

    NNout(:,3) = reshape(repmat(yout_qlkz.dfe_GB,1,size(te,1)),1,size(te,1).*size(te,2));
    NNout(:,4) = reshape(repmat(yout_qlkz.vce_GB+yout_qlkz.vte_GB,1,size(te,1)),1,size(te,1).*size(te,2));

    De   = reshape(NNout(:,3),flip(insize))' .* kD;
    Ve   = reshape(NNout(:,4),flip(insize))' .* kV;

  else

    NNout(:,3) = yout_qlkz.dfe_GB;
    NNout(:,4) = yout_qlkz.vce_GB+yout_qlkz.vte_GB;

    De   = reshape(NNout(:,3),insize) .* kD;
    Ve   = reshape(NNout(:,4),insize) .* kV;

    
  end
  varargout{1} = Chie;
  varargout{2} = Chii;
  varargout{3} = De;
  varargout{4} = Ve;
end
