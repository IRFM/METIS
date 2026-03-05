function self = dunne_pedestal(params,in)

args.ip              = 1e6;
args.bt              = 2.5;
args.a               = 0.5;
args.rmaj            = 1.65;
args.kappa           = 1.7;
args.delta           = 0.35;
args.betapol         = 1.5;
args.w_pre           = 0.1;
args.density         = 5e19;
args.dens_shift      = 0;
args.pres_shift      = 0;
args.nesep_scale     = 4;
args.tesep           = 100;
args.zeff            = 2;
args.Amain           = 2;
args.profiles        = 'EPED';
args.profiles_onoff  = false;
args.scaling         = 'Faitsch';
args.printinfo       = true;
args.self_consistent = true;
args.npts            = 500;
if  nargin == 0
    self = args;
    return
end
if nargin > 1
    self = in;
else
    self = [];
end
noms = fieldnames(args);
for k=1:length(noms)
    if ~isfield(params,noms{k}) || isempty(params.(noms{k}))
        params.(noms{k}) = args.(noms{k});
    end
end
self.params = params;
if ~isfield( self,'npts') || isempty(self.npts)
    self.npts = self.params.npts;
end
consts = cphys;
zeff = self.params.zeff;
%m_p = 9.382720813e8/consts.c .^ 2;
m_p =  consts.mp ./ consts.e;  % strange unit !!!
m_i = self.params.Amain * m_p;


self.ne_ped = self.params.density;
nGW = self.params.ip./ 1e6 ./ (pi * self.params.a .^ 2) .* 1e20;
nnorm = self.params.density / nGW;
        
self.bpol = consts.mu0 .* self.params.ip ./ (2. .* pi .* self.params.a .* (1. + self.params.kappa .^ 2) ./ 2.);
self.qhat = self.params.a .* self.params.bt ./ (self.bpol .* self.params.rmaj);

%Scaling from eq 105, Masters' thesis, J. Puchmayr
self.beta_polped = (0.686 * self.params.kappa .^ 0.5 * (1. + self.params.delta).^ 1.68 .* ...
                    self.qhat .^ 1.61 .* self.params.betapol .^ 0.33 .* ...
                    nnorm .^ 0.06 * self.params.w_pre .^1.29);
if self.params.self_consistent
    if self.params.printinfo
            fprintf('Doing self consistent runs: ')
    end
    bptemp = self.params.betapol - 0.1;
    count = 100;
    while (abs(bptemp - self.params.betapol) > 0.01) && (count> 0)
        if self.params.printinfo
            fprintf('.');
        end
        self.beta_polped = (0.686 .* self.params.kappa .^ 0.5 * (1. + self.params.delta) .^ 1.68 .* ...
                           self.qhat .^1.61 * self.params.betapol .^ 0.33 .* self.params.w_pre .^ 1.29);
        bptemp = self.params.betapol;
        self.params.betapol  = self.beta_polped * 3.;
        count = count -1;
    end
    if self.params.printinfo
        fprintf('\n');
    end
end

self.p_ped = self.beta_polped  .*  self.bpol .^ 2 / (2.  .*  consts.mu0);
self.eped_width = self.params.w_pre .* sqrt(self.beta_polped);

self.pav = self.params.betapol .*  self.bpol .^ 2 / (2.  .*  consts.mu0);
self.beta_tor_glob = self.pav ./ self.params.bt .^ 2  .*  (2.  .*  consts.mu0);
self.beta_n_glob = self.beta_tor_glob .* self.params.bt' .* self.params.a ./ self.params.ip  .*  1e6  .*  100;

self.te_ped = self.p_ped ./ 2. ./ consts.e ./ self.ne_ped;
self.ne_sep = self.ne_ped ./ 3.;

%recalculate q_cyl for more accurate q95 estimate
aspect = self.params.rmaj' ./ self.params.a;
kappa_eff = sqrt((1. + self.params.kappa  .^ 2  .*  (1. + 2.  .*  self.params.delta  .^ 2 - 1.2  .*  self.params.delta  .^ 3)) / 2.);
Bpol = consts.mu0  .*  self.params.ip  ./ (2.  .*  pi  .*  self.params.a   .*   kappa_eff);
self.q_cyl = self.params.bt  ./ Bpol  .*  kappa_eff ./ aspect;

zlnl_e = 31.3 - log(sqrt(self.ne_ped) / self.te_ped);
% '''
% nustar = (np.sqrt(2.)  .*  consts.e .^ 2 / (12.  .*  np.pi .^ 1.5  .*  consts.epsilon_0 .^ 2)
%  .*  self.qhat  .*  self.params.rmaj   .*  self.ne_ped  .*  zlnl_e  .*  zeff / (self.te_ped .^ 2
%  .*  (self.params.a  / self.params.rmaj ) .^ 1.5))
% '''
nustar = (6.921e-18  .*  self.q_cyl  .*  self.params.rmaj   .*  zeff  .*  zlnl_e ./ ...
         ((self.params.a  ./ self.params.rmaj ) .^ 1.5)  .*  self.ne_ped ./ self.te_ped .^ 2);
rhostar = (sqrt( 2.  .*  m_i  .*  self.te_ped ) ./ ...
         ( 1.  .*  self.params.a   .*  self.params.bt ) );

beta_tor = (self.p_ped  .*  2. ./ (self.params.bt  .^ 2 ./ (2. .* consts.mu0)))  .*  100.;
beta_n = (beta_tor  .*  self.params.bt   .*  self.params.a  ./ (self.params.ip  / 1e6));

if self.params.printinfo
    disp(" ")
    fprintf('beta_polped = %g, eped_width = %g ,p_ped = %g\n', ...
             self.beta_polped,self.eped_width,self.p_ped);
    fprintf('q_hat = %g, te_ped = %g, ne_ped = %g\n', ...
            self.qhat,self.te_ped,self.ne_ped/1e19);
    fprintf('rho^* = %g, nu^* = %g, beta_ped = %g, beta_n_ped = %g, beta_n_global = %g\n', ...
            rhostar,nustar,beta_tor,beta_n,self.beta_n_glob);
end
self.psin = linspace(0,1,self.npts);

%calculate average plasma pressure from beta_pol - scale up te to meet global pressure
self.pav = self.params.betapol   .*  self.bpol .^ 2 ./ (2  .*  consts.mu0);
if self.params.printinfo
	fprintf('pav = %g\n',self.pav);
end
if self.params.profiles_onoff
    [self.te,self.ne] = profiles_eped(self.te_ped, self.psin, self.eped_width, self.ne_ped, ...
        self.params.nesep_scale , self.params.tesep ,self.pav, ...
        self.params.dens_shift , self.params.pres_shift);
end
self.dimensionless.rhostar    = rhostar;
self.dimensionless.nustar     = nustar;
self.dimensionless.beta_tor   = beta_tor;
self.dimensionless.beta_n     = beta_n;
self.dimensionless.beta_pol   = self.beta_polped;
self.dimensionless.f_GW       = self.ne_ped / nGW;

%power scalings
%PLH (Martin):
surf = 4. .* pi .^ 2  .*  self.params.a   .*  self.params.rmaj; 
plh = 0.0488  .*  (self.ne_ped ./ 1e20) .^ 0.717  .*  self.params.bt  .^ 0.803  .*  surf .^ 0.941;
pmin = 1.2  .*  plh;
% q / qGB:
%psep_qqGB = 10.  .*  (self.ne_ped/1e19)  .*  (self.te_ped/1e3) .^ 2.5 / self.params.bt  .^ 2
psep_qqGB = pmin  .*  self.params.bt  .^ 2 ./ (self.ne_ped ./ 1e19  .*  (self.te_ped ./ 1e3) .^ 2.5);
%divertor challenge
%psep_div = 10.  .*  self.qhat  .*  self.params.rmaj  / self.params.a   .*  self.params.rmaj  / self.params.bt 
psep_div = pmin  .*  self.params.bt  ./ (self.qhat  .*  self.params.rmaj  ./ self.params.a   .*  self.params.rmaj );

self.powers.plh        = plh;
self.powers.psep_qqGB  = psep_qqGB;
self.powers.psep_div   = psep_div;
[tedge,l_t,l_p,l_n,alphat,alpha_mhd,beta] = sep_params(self.params.kappa ,self.params.delta ,self.params.ip , ....
                                            self.params.bt ,self.params.rmaj ,self.params.a , ...
                                            1.5,pmin .* 1e6,self.ne_sep,self.params.scaling);

zlnl_e = 31.3 - log(sqrt(self.ne_sep) ./ tedge);
nustar_sep = (6.921e-18  .*  self.q_cyl  .*  self.params.rmaj   .*  zeff  .*  zlnl_e ./ ...
             ((self.params.a  ./ self.params.rmaj ) .^ 1.5)  .*  self.ne_sep ./ tedge .^ 2);
%nustar_sep = (np.sqrt(2.)  .*  consts.e .^ 2 / (12.  .*  np.pi .^ 1.5  .*  consts.epsilon_0 .^ 2)
%           .*  self.qhat  .*  self.params.rmaj   .*  self.ne_sep  .*  zlnl_e  .*  zeff /
%              (tedge .^ 2  .*  (self.params.a  / self.params.rmaj ) .^ 1.5))
self.sep_params.te_sep     = tedge;
self.sep_params.ne_sep     = self.ne_sep;
self.sep_params.alpha_t    = alphat;
self.sep_params.nustar_sep = nustar_sep;




function rhos = rhos(Tesep,Bpol) 
    rhos = 1.44E-4  .*  sqrt(Tesep) ./ Bpol;

function [te,lambda_t] = tesep_det(psep,q_cyl,aspect,kappa_eff,zeff,rhos_pol,alpha_t,scaling)

if (nargin < 8) || isempty(scaling)
   scaling ='Eich';
end
switch scaling
    case 'Eich'
        lambda_t = 2.1  .*  (1. + 2.1  .*  alpha_t.^1.7)  .*  rhos_pol;
    case 'Faitsch'
        lambda_t = 2.57  .*  (1. + 0.38  .*  alpha_t.^3.3)  .*  rhos_pol; %actual scaling from M. Faitsch
    otherwise
        disp('Error, no such scaling (yet)!')
        lambda_t = 2.1  .*  (1. + 2.1  .*  alpha_t.^1.7) .*  rhos_pol;
end
        
lambda_q = 2. ./ 7.  .*  lambda_t;
f_zeff = 0.672 + 0.076  .*  zeff .^ 0.5 + 0.252 .*  zeff;
k_spitzer = 2600. ./ f_zeff;

te = ((7.  .*  psep  .*  q_cyl.^2  .*  aspect) ./ (16.  .*  k_spitzer  .*  kappa_eff  .*  lambda_q)).^(2./7.);


function [tedge,lambda_t,lambda_p,lambda_n,alpha_t,alpha_mhd,beta] = ...
          sep_params(kappa,delta,ip,bt,rgeo,aminor,zeff, ...
          psep,dens,scaling,tstart,tmin)
if (nargin < 8) || isempty(psep)
  psep = 7e6;
end
if (nargin < 9) || isempty(dens)
  dens = 4e19;
end
if (nargin < 10) || isempty(scaling)
    scaling='Faitsch';
end
if (nargin < 11) || isempty(tstart)
    tstart=false;
end
if (nargin < 12) || isempty(tmin)
    tmin=0.;
end

mu0 = pi  .*  4e-7;
if tstart
    tesep = tstart;
else
    tesep = 100.;
end

aspect = rgeo ./ aminor;
kappa_eff = sqrt((1. + kappa .^ 2  .*  (1. + 2. .*  delta .^ 2 - 1.2  .*  delta .^3)) ./ 2.);
Bpol = mu0  .*  ip ./ (2.  .*  pi  .*  aminor  .*   kappa_eff);
q_cyl = bt ./ Bpol  .*  kappa_eff ./ aspect;

%Modification of scaling for q_edge from data base
BpCirc = mu0  .*  ip ./ (2.  .*  pi  .*  aminor);
q_circ = bt ./ BpCirc ./ aspect;
q_edge = q_circ  .*  1.11  .*  kappa .^ 1.6  .*  (1.+delta) .^ 0.3;
tedge = tesep;
nedge = dens;
diff = 1.;
while diff > 0.5
    old = tesep;
    alpha_t = 3.13e-18  .*  rgeo  .*  q_cyl .^2  .*  nedge ./ tedge .^ 2  .*  zeff;
    rho_s = rhos(tesep,Bpol);
    [tesep,lambda_t] = tesep_det(psep,q_cyl,aspect,kappa_eff,zeff,rho_s,alpha_t,scaling);
    lambda_n = 2.9  .*  (1. + 10.4  .*  alpha_t .^ 2.5)  .*  rho_s;
    
    tedge = tesep  .*  exp(1e-3 ./ lambda_t); %tesep+tesep/lambda_t .* 1e-3
    nedge = dens  .*  exp(1e-3 ./ lambda_n); %dens+dens/lambda_n .* 1e-3
    diff = abs(tesep - old);
end

alpha_t = 3.13e-18 .*  rgeo  .*  q_cyl .^ 2  .*  nedge ./ tedge .^ 2  .*  zeff;
rho_s = rhos(tesep,Bpol);
[tesep,lambda_t] = tesep_det(psep,q_cyl,aspect,kappa_eff,zeff,rho_s,alpha_t,scaling);
lambda_n = 2.9  .*  (1. + 10.4  .*  alpha_t .^ 2.5) .*  rho_s;

tedge = tesep  .*  exp(1e-3 ./ lambda_t); %+tesep/lambda_t .* 1e-3
nedge = dens  .*  exp(1e-3 ./ lambda_n); %+dens/lambda_n .* 1e-3

alpha_t = 3.13e-18  .*  rgeo  .*  q_cyl .^ 2  .*  nedge ./ tedge .^ 2  .*  zeff;
%determine lambda_p and alphas
nue = pi  .*  q_cyl  .*  rgeo ./ (1.03e16)  .*  nedge ./ tedge .^ 2  .*  zeff;

switch scaling
    case 'Eich'
        lambda_p = 1.2  .*  (1. + 3.6  .*  alpha_t .^ 1.9)  .*  rho_s;
    case 'Faitsch'
        lambda_p = 1.55  .*  (1. + 0.61 .*  alpha_t .^ 2.35) .*  rho_s;
    otherwise
        print('Error, no such scaling (yet)!')
        lambda_p = 1.2  .*  (1. + 3.6 .*  alpha_t .^ 1.9) .*  rho_s;
        %lambda_p = 1.22  .*  (1. + 0.005  .*  nue .^ 1.8)  .*  rho_s
end
btot = sqrt(bt .^ 2 + Bpol .^ 2);
beta = 4. .*  mu0  .*  tedge  .*  nedge  .*  1.602e-19 ./ btot .^ 2;  
alpha_mhd = rgeo  .*  q_edge .^ 2  .*  beta ./ lambda_p;


function [te, dens]  = profiles_eped(tetop, psin, width, dens_top, nesep_scale, te_sep, pav,dens_shift, pres_shift)

if (nargin  < 8) || isempty(dens_shift)
    dens_shift=0;
end
if (nargin  < 9) || isempty(pres_shift)
    pres_shift=0;
end

mu0 = pi .* 4e-7;
npts = length(psin);

psi_mid = 1-width ./ 2;
psi_ped = 1-width+pres_shift;

alphan1=1.1;
alphan2=1.1;
an1=1e19  ;                                                                     %-> WHERE??
nesep = dens_top ./ nesep_scale; %included as input param for scans

dens = nesep + (dens_top-nesep) ./ 1.75  .*  (tanh(2  .*  (1-psi_mid) ./ width) - tanh(2  .*  (psin-psi_mid) ./ width)); %tanh's -> 2 .* tanh(1) vs tanh(1) - tanh(-inf) 1.52, 1.75
dens = dens  .*  double(dens >= 0);

%The latter elements of psin/psi_ped are larger than 1 -> negative powers are invalid
warning off
core_mod = an1  .*  double((1-psin ./ psi_ped) >= 0.5)   .*  (1-(psin ./ psi_ped).^ alphan1).^ alphan2;
warning on
core_mod(~isfinite(core_mod)) = 0;

dens=dens + core_mod;
%if have density shift, then move entire density profile
%want to do this, as experiment shows whole profile shifts, not
%just central pedestal location at constant separatrix value
if ~isempty(dens_shift) || (dens_shift == 0)
    dens =interp1(psin + dens_shift, dens,psin,'linear',NaN);
end

alphaT1=1.2;
alphaT2=1.4;

te = te_sep + (tetop-te_sep) ./ 1.75  .*  (tanh(2  .*  (1 - psi_mid) ./ width) - tanh(2  .*  (psin - psi_mid) ./ width));
te = te  .*  double(te >= 0);


%The latter elements of vol/vol_ped are larger than 1 -> negative powers are invalid -> set to 0 since irrelevant (heavyside)
warning off
core_heaviside = double((1 - psin ./ psi_ped) >= 0.5)   .* (1 - (psin ./ psi_ped) .^ alphaT1).^ alphaT2;
warning on
core_heaviside(~isfinite(core_heaviside)) = 0;  
%total_mod = sum(core_heaviside .* dvol) ./ vol_ped;
%disp(total_mod)



at1 = 3  .*  tetop;
%vol_sep = np.interp(1.0,psin,vol)
%at1 = (pav / (1.602e-19  .*  (1 + np.mean(scale))) - integral_to(dens  .*  te, x=vol, to=vol_sep) / vol_sep) / (integral_to(dens  .*  core_heaviside, x=vol, to=vol_sep) / vol_sep)

te = te + at1 .*  core_heaviside;

if isempty(pres_shift) ||  (pres_shift == 0)
    dens=interp1(psin+dens_shift, dens, psin, 'linear',NaN);
    te=interp1(psin+pres_shift, te,psin, 'linear',NaN);
end



