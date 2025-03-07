%[chie,chii,flux_qe,flux_qi]  = z0qlkANN(post.z0dinput.option,post.zerod,post.profil0d); 
function [chie,chii,chii_neo,flux_qe,flux_qi,flux_qi_neo,chie_err,chii_err,flux_qe_err,flux_qi_err] = z0qlkANN(option,zerod,profil,cons,test)

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

% hh_effect
if nargin < 4
  hh_effect = 0;
else
  hh_effect = interp1(cons.temps,cons.hmore,profil.temps,'linear','extrap');
end

% parameters for coefficient computation
if option.xiioxie > 0 
    qiqe_ratio      = option.xiioxie;
else
    qiqe_ratio      = 2.0;                       % Set the qi/qe ratio (since original network was for adiabatic electrons)
end
qi_multiplier   = 1.5;                       % Set the multiplier of qi above output (to translate from adiabatic to kinetic electron network)
chi_min         = 0.1;
% Amain
Amain = option.gaz;
ve     = ones(size(profil.xli));

% neural network data
net = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','adiabatic_e_5D_ITG.mat'));
net = net.net;

% loop on time slices
chii     = NaN .* ones(length(profil.temps),length(profil.xli));
chii_err = NaN .* ones(length(profil.temps),length(profil.xli));
for k=1:length(profil.temps)
  t     = profil.temps(k);
  ind0d = find(zerod.temps >= t,1);
  if isempty(ind0d)
      if zerod.temps(1) >= t
	  ind0d = 1;
      elseif zerod.temps(end) <= t
	  ind0d = length(zerod.temps);
      end
  end
  RLFS = profil.Raxe(k,:) .* (1 + profil.epsi(k,:));
  te   = profil.tep(k,:);
  ted  = pdederive(RLFS,te,0,2,2,1);
  ti   = profil.tip(k,:);
  tid  = pdederive(RLFS,ti,0,2,2,1);
  q    = profil.qjli(k,:);
  s    = pdederive(profil.rmx(k,:),q,0,2,2,1) .* profil.rmx(k,:) ./ q;
  geo.a  = profil.Raxe(k,end) .* profil.epsi(k,end);
  geo.r0 = profil.Raxe(k,end);
  geo.b0 = profil.fdia(k,end) ./ profil.Raxe(k,end);
  % NOTE!! The actual neural network has to be added by hand as net
  %timing: ~800 us for one call!
  %%%%%%%%%%%%%%%%%%%%%%%%% Call Qualikiz ANN %%%%%%%%%%%%%%%%%%%%%%%%%
  %chii(k,:) = qlkANN(te,ti,ted,tid,q,s,geo,Amain,net);
  %qlkANN(te,ti,tep,tip,q,shear,R0,a0,B0,Amain,pp)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  chi0 = NaN .* ones(101,length(profil.xli));
  for l = 1:101
	chi0(l,:) = qlkANN(te' .* (0.9 + 0.2 .* rand(1)),ti' .* (0.9 + 0.2 .* rand(1)),ted' .* (0.8 + 0.4 .* rand(1)),tid' .* (0.8 + 0.4 .* rand(1)),q'.* (0.9 + 0.2 .* rand(1)),s'.* (0.9 + 0.2 .* rand(1)), ...
                    geo.r0,geo.a,geo.b0,Amain,net);
  end
  for l = 1:21
      c = chi0(:,l);
      if any(c>0)
	  chii(k,l) = mean(c(c>0));
	  chii_err(k,l) = std(c(c>0));
      else
	  chii(k,l) = NaN;      
	  chii_err(k,l) = NaN;
      end
  end
end

% scaling
if any(hh_effect > 0)
  chii = chii * qi_multiplier ./ (hh_effect * ve);
  chii_err = chii_err * qi_multiplier ./ (hh_effect * ve);
else
  chii = chii * qi_multiplier;
  chii_err = chii_err * qi_multiplier;
end
chie = chii / qiqe_ratio;
chie_err = chii_err / qiqe_ratio;
% for graph
chie_raw = chie;
chii_raw = chii;
%filter outputs
chie(chie < chi_min) = 0;
chii(chii < chi_min) = 0;
chie(~isfinite(chie)) = 0;
chii(~isfinite(chii)) = 0;
chie_err(~isfinite(chie_err)) = 0;
chii_err(~isfinite(chii_err)) = 0;

% neoclassical part (use se same on ion and electron)
switch option.gaz
    case 5
        iso = zerod.nhem./ zerod.n1m;
    case 11
        iso = zerod.nTm ./ zerod.n1m;
    otherwise
        iso = zeros(size(zerod.nem));
end
chii_neo = neo_Hinton(profil,option,Amain,iso);
chii_neo_raw = chii_neo;

% for Temperature recomputation
chie_ref = max(0,profil.xie .* profil.grho2 ./ profil.grho - chii_neo / qiqe_ratio);
chii_ref = max(0,profil.xii .* profil.grho2 ./ profil.grho - chii_neo);
chie(chie == 0) = chie_ref(chie == 0);
chii(chii == 0) = chii_ref(chii == 0);
% to test numerical integration:
if nargin > 4
  chie = chie_ref;
  chii = chii_ref;
end

% Local gradient
rhomax = profil.rmx(:,end) * ve;
ted1  = pdederive(profil.xli,profil.tep,0,2,2,1) ./ pdederive(profil.xli,profil.rmx,2,2,2,1);
tid1  = pdederive(profil.xli,profil.tip,0,2,2,1) ./ pdederive(profil.xli,profil.rmx,2,2,2,1);

% output flux
flux_qe = - chie .* profil.nep .* ted1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax; 
flux_qe_err = - chie_err .* profil.nep .* ted1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax; 
flux_qi = - chii .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax; 
flux_qi_err = - chii_err .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax; 
flux_qineo = - chii_neo .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax; 

% recomputation of temperature 
ted1_out = - profil.qe ./ (max(chi_min,chie + chii_neo / qiqe_ratio) .* profil.nep .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
ted1_out_m = - profil.qe ./ (max(chi_min,chie + chie_err + chii_neo / qiqe_ratio) .* profil.nep .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
ted1_out_p = - profil.qe ./ (max(chi_min,chie - chie_err + chii_neo / qiqe_ratio) .* profil.nep .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
tid1_out = - profil.qi ./ (max(chi_min,chii + chii_neo) .* profil.nip .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
tid1_out_m = - profil.qi ./ (max(chi_min,chii + chii_err + chii_neo) .* profil.nip .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
tid1_out_p = - profil.qi ./ (max(chi_min,chii - chii_err + chii_neo) .* profil.nip .* phys.e .* profil.vpr .* profil.grho) .* rhomax .^ 2;
% integration + boudary condition
te_out = compute_temp(profil.xli,ted1_out,profil.tep);
te_out_m = compute_temp(profil.xli,ted1_out_m,profil.tep);
te_out_p = compute_temp(profil.xli,ted1_out_p,profil.tep);
ti_out = compute_temp(profil.xli,tid1_out,profil.tip);
ti_out_m = compute_temp(profil.xli,tid1_out_m,profil.tip);
ti_out_p = compute_temp(profil.xli,tid1_out_p,profil.tip);

% Conversion  --> Cronos transport coefficients
chie = chie_raw .* profil.grho ./ profil.grho2;     % Electron heat diffusivity
chii = chii_raw .* profil.grho ./ profil.grho2;     % Ion heat diffusivity
chii_neo = chii_neo_raw .* profil.grho ./ profil.grho2;     % Ion heat diffusivity
chie_err = chie_err .* profil.grho ./ profil.grho2;     % Electron heat diffusivity
chii_err = chii_err  .* profil.grho ./ profil.grho2;     % Ion heat diffusivity

% graph
if nargout > 0
  return
end

% flux sortant

hz =findobj(0,'type','figure','tag','Chi_net');
if isempty(hz)
  	  hz=figure('tag','Chi_net','name','gradient length');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)
zplotprof(gca,profil.temps,profil.xli,profil.xie,'color','b');
zplotprof(gca,profil.temps,profil.xli,chie + chii_neo ,'color','r');
zplotprof(gca,profil.temps,profil.xli,(chie-chie_err + chii_neo),'color','m');
zplotprof(gca,profil.temps,profil.xli,(chie+chie_err + chii_neo),'color','m');
set(gca,'YScale','linear');
ylabel('Chi_e (m^2/s)');
legend('METIS','Qlk ANN');
xlabel('r/a');
z0loglin(gca);

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.xii,'color','b');
zplotprof(gca,profil.temps,profil.xli,chii + chii_neo,'color','r');
zplotprof(gca,profil.temps,profil.xli,(chii-chii_err) + chii_neo,'color','m');
zplotprof(gca,profil.temps,profil.xli,(chii+chii_err) + chii_neo,'color','m');
zplotprof(gca,profil.temps,profil.xli,chii_neo,'color','g');
set(gca,'YScale','linear');
ylabel('Chi_i (m^2/s)');
legend('METIS','Qlk ANN');
xlabel('r/a');
z0loglin(gca);

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,profil.tep ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,te_out ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,te_out_m ./ 1e3,'color','m');
zplotprof(gca,profil.temps,profil.xli,te_out_p ./ 1e3,'color','m');
set(gca,'YScale','linear');
ylabel('T_e (keV)');
legend('METIS','Qlk ANN');
xlabel('r/a');
z0loglin(gca);

subplot(2,2,4)
zplotprof(gca,profil.temps,profil.xli,profil.tip ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,ti_out ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,ti_out_m ./ 1e3,'color','m');
zplotprof(gca,profil.temps,profil.xli,ti_out_p ./ 1e3,'color','m');
set(gca,'YScale','linear');
ylabel('T_i (keV)');
legend('METIS','Qlk ANN');
xlabel('r/a');
z0loglin(gca);

function te_out = compute_temp(xli,ted1_out,tep);

ve     = ones(size(xli));

te_out = fliplr(cumtrapz(fliplr(xli),fliplr(ted1_out),2));
te_out = te_out  + (tep(:,end - 1) - te_out(:,end-1)) * ve;
te_out(:,end) = tep(:,end);
te_out(:,1) = te_out(:,2);


%ref: Wesson
function chii_neo = neo_Hinton(profil,option,Amain,iso)

ve     = ones(size(profil.xli));
vt     = ones(size(profil.temps));
epsi   = profil.epsi;
delta  = profil.Raxe - profil.Raxe(:,end) * ve;
a      = (profil.Raxe(:,end) .* epsi(:,end)) * ve;
delta_prim = pdederive(profil.xli,delta,0,2,2,1) ./ a;

f1     = (1 + 3/2 .* (epsi .^ 2 +  epsi .* delta_prim) + 3/8 .* epsi .^ 3 .* delta_prim) ./ ( 1 + 1/2 .* epsi .*  delta_prim);
f2     = sqrt(1 - epsi .^ 2) .* (1 + epsi .* delta_prim ./ 2) ./ (1 + delta_prim ./ epsi .* (sqrt(1 - epsi .^ 2)  - 1));

nI_ZI2    = (profil.nep .* profil.zeff - profil.n1p);
switch option.gaz
    case 4
        alpha  = nI_ZI2 ./ profil.nhep ./ 4;
    case 5
        alpha  = nI_ZI2 ./ profil.n1p ./ (1 + 4 .* (iso * ones(1,size(profil.n1p,2))));
    case 11
        nI_ZI2 = nI_ZI2 + profil.n1p .* (1 + iso * ones(1,size(profil.n1p,2)));
        alpha  = nI_ZI2 ./ profil.n1p ./ (1 + 25 .* (iso * ones(1,size(profil.n1p,2))));        
    otherwise
        alpha  = nI_ZI2 ./ profil.n1p;
end
lnldi   = 17.3 - 0.5 .* log(profil.nep./1e20) + log(profil.tip ./1e3);
tau_i   = 6.60e17 .* sqrt(Amain) .* (profil.tip ./1e3).^(3/2) ./ profil.nip ./ lnldi;
vth_i   = 3.09e5  .* sqrt(profil.tip ./1e3 ./ Amain);
B       = sqrt(profil.bpol .^ 2 + (profil.fdia .* profil.ri) .^ 2);
rho_i   = 4.57e-3 .* sqrt(Amain .* profil.tip ./1e3) ./ B;
nustar  = profil.Raxe .* profil.qjli ./ (epsi .^ (3/2) .* vth_i .* tau_i);
mustar_i  = nustar .* (1 + 1.54 .* alpha);

part1   = 0.66 .* (1 + 1.54 .* alpha) + (1.88 .* sqrt(epsi) - 1.54 .* epsi) .* (1 + 3.75 .* alpha);
part1   = part1 ./ (1 + 1.03 .* sqrt(mustar_i) + 0.31 .* mustar_i);

part2   = 0.59 .* mustar_i .* epsi ./ (1 + 0.74 .* mustar_i .* epsi .^ (3/2));
part2   = part1 .* (1 + (1.33 .* alpha .* (1 + 0.60 .* alpha)) ./ (1 + 1.79 .* alpha));

chii_neo = epsi .^(-3/2) .* profil.qjli .^ 2 .* rho_i .^ 2 ./ tau_i .* ...
           (part1 .* f1 + part2 .* (f1 -f2));
           
