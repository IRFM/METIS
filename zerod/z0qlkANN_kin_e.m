% Attention sorties = definitions Qualikiz !
%[chie,chii,chii_neo,D,V,chie_err,chii_err,D_err,V_err] = z0qlkANN_kin_e(post.z0dinput.option,post.zerod,post.profil0d); 
function [chie,chii,chii_neo,D,V,chie_err,chii_err,D_err,V_err] = z0qlkANN_kin_e(option,zerod,profil,mode,test)

% recalcul des profils seulement
if nargin < 4
  mode = 'full';
elseif isnumeric(mode)
  % compatibility with older version
  void = mode;
  mode = 'full';
end

if ~isfield(option,'Sn_fraction')
    Sn_fraction = 0;
else
    Sn_fraction = option.Sn_fraction;
end
if ~isfield(option,'extended_qei')
    extended_qei = 'off';
else
    extended_qei = option.extended_qei;
end



first_call = 1;
% question on limit in which QLK is used if not already defined
if isfield(option,'prediction_domain')
  prediction_domain = option.prediction_domain;
  first_call = 0;
else
  prediction_domain = questdlg('Choose prediction domain ?','QLKNN-4Dkin prediction interval','Restricted','Clamp','Extrapolated','Restricted');
  first_call = 1;
end                         

% model used 
% question on model version
% valid values are (string):
% '2016', '2017'
if isfield(option,'model_qlk_version')
  model_qlk_version = option.model_qlk_version;
  first_call = 0;
else
  model_qlk_version = questdlg('Choose prediction domain ?','QLKNN-4Dkin prediction interval','2016','2017','2016');
  first_call = 1;
end                         

% warning for particles transport
if first_call == 1
  switch model_qlk_version
  case '2017'
    warndlg('Particles transport is not yet implemented in this version', 'Particles tranport');
  end
end

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
%rechantillonage en temps
meff  = interp1(zerod.temps,zerod.meff,profil.temps,'linear','extrap');
modeh = interp1(zerod.temps,zerod.modeh,profil.temps,'linear','extrap');
indice_inv  = interp1(zerod.temps,zerod.indice_inv,profil.temps,'linear','extrap');
wth  = interp1(zerod.temps,zerod.wth,profil.temps,'linear','extrap');
taue = interp1(zerod.temps,zerod.taue,profil.temps,'linear','extrap');
% parameters for coefficient computation
if option.xiioxie > 0 
    qiqe_ratio      = option.xiioxie;
else
    qiqe_ratio      = 2.0;                       % Set the qi/qe ratio (since original network was for adiabatic electrons)
end
qi_multiplier   = 1.5;                       % Set the multiplier of qi above output (to translate from adiabatic to kinetic electron network)
chi_min         = 0.1;
ve     = ones(size(profil.xli));


% profil de densites
ve = ones(size(profil.xli));
n1m  = interp1(zerod.temps,zerod.n1m,profil.temps,'linear','extrap');
nDm  = interp1(zerod.temps,zerod.nDm,profil.temps,'linear','extrap');
nTm  = interp1(zerod.temps,zerod.nTm,profil.temps,'linear','extrap');
nhep = profil.nhep;
nHp  = profil.n1p .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);
nDp  = profil.n1p .* ((nDm ./ n1m) * ve);
nTp  = profil.n1p .* ((nTm ./ n1m) * ve);
nzp  = profil.nzp;
if isfield(profil,'nwp')
	nwp  = profil.nwp;
else
	nwp  = 0 .* nzp;		
end

switch model_qlk_version
case '2017'
  if isdeployed
    net    = load('kinetic_e_5D_ITG_20170523.mat');
  else
    net    = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kinetic_e_5D_ITG_20170523'));  
  end
  data.net    = net.net;
  data.net.constrainrlti = true;
  data.net.constraintite = true;
  data.net.constrainq = true;
  data.net.constrains = true;
  data.net.rltimin = 2;
  data.net.rltimax = 12;
  data.net.titemin = 0.3;
  data.net.titemax = 3.0;
  data.net.qmin = 1;
  data.net.qmax = 5;
  data.net.shearmin = 0.1;
  data.net.shearmax = 3;
  data.net.shift_chi_e = 0.75;
  data.net.chiemin = 0.1;
  data.net.chiimin = 0.1;
  
  % display choosen parameters
  if first_call == 1
      disp('=============================================================')
      disp('Neural network is call with following parameters:')
      fprintf('constrainrlti = %g\n',data.net.constrainrlti);
      fprintf('constraintite = %g\n',data.net.constraintite);
      fprintf('constrainq =  %g\n',data.net.constrainq);
      fprintf('constrains =  %g\n',data.net.constrains);
      fprintf('rltimin =  %g\n',data.net.rltimin);
      fprintf('rltimax =  %g\n',data.net.rltimax);
      fprintf('titemin =  %g\n',data.net.titemin);
      fprintf('titemax =  %g\n',data.net.titemax);
      fprintf('qmin =  %g\n',data.net.qmin);
      fprintf('qmax =  %g\n',data.net.qmax);
      fprintf('shearmin =  %g\n',data.net.shearmin);
      fprintf('shearmax =  %g\n',data.net.shearmax);
      fprintf('shift_chi_e =  %g\n',data.net.shift_chi_e);
      fprintf('chiemin =  %g\n',data.net.chiemin);
      fprintf('chiimin =  %g\n',data.net.chiimin);
      disp('=============================================================')
      disp('');
  end
  
  % for density waiting for complesion of new version
  if isdeployed
      net         = load('kin_e_5D_ITG_dfe.mat');
      data.netdfe = net.net;
      net         = load('kin_e_5D_ITG_vce.mat');
      data.netvce = net.net;
      net         = load('kin_e_5D_ITG_vte.mat');
      data.netvte = net.net;  
  else
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_dfe.mat'));
      data.netdfe = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_vce.mat'));
      data.netvce = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_vte.mat'));
      data.netvte = net.net;
  end
otherwise
  % neural network data
  if isdeployed
      net         = load('adiabatic_e_5D_ITG.mat');
      data.net    = net.net;
      net         = load('kin_e_5D_ITG_ief.mat');
      data.netief = net.net;
      net         = load('kin_e_5D_ITG_eef.mat');
      data.neteef = net.net;
      net         = load('kin_e_5D_ITG_dfe.mat');
      data.netdfe = net.net;
      net         = load('kin_e_5D_ITG_vce.mat');
      data.netvce = net.net;
      net         = load('kin_e_5D_ITG_vte.mat');
      data.netvte = net.net;  
  else
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','adiabatic_e_5D_ITG.mat'));
      data.net    = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_ief.mat'));
      data.netief = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_eef.mat'));
      data.neteef = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_dfe.mat'));
      data.netdfe = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_vce.mat'));
      data.netvce = net.net;
      net         = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kin_e_5D_ITG_vte.mat'));
      data.netvte = net.net;
  end
end
% loop on time slices
chii     = NaN .* ones(length(profil.temps),length(profil.xli));
chii_err = NaN .* ones(length(profil.temps),length(profil.xli));
chie     = NaN .* ones(length(profil.temps),length(profil.xli));
chie_err = NaN .* ones(length(profil.temps),length(profil.xli));
D        = NaN .* ones(length(profil.temps),length(profil.xli));
D_err    = NaN .* ones(length(profil.temps),length(profil.xli));
Vt       = NaN .* ones(length(profil.temps),length(profil.xli));
Vt_err   = NaN .* ones(length(profil.temps),length(profil.xli));
Vc       = NaN .* ones(length(profil.temps),length(profil.xli));
Vc_err   = NaN .* ones(length(profil.temps),length(profil.xli));
V       = NaN .* ones(length(profil.temps),length(profil.xli));
V_err   = NaN .* ones(length(profil.temps),length(profil.xli));


switch mode
case 'full'
    if nargout == 0
      hwaitbar = waitbar(0,'Please wait...');
    end

    for k=1:length(profil.temps)
      if nargout == 0
	  waitbar(0.5 .* k ./ length(profil.temps),hwaitbar,'Please wait...');
      end
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
      %
      Amain = meff(k);
      % NOTE!! The actual neural network has to be added by hand as net
      %timing: ~800 us for one call!
      %%%%%%%%%%%%%%%%%%%%%%%%% Call Qualikiz ANN %%%%%%%%%%%%%%%%%%%%%%%%%
      %chii(k,:) = qlkANN(te,ti,ted,tid,q,s,geo,Amain,net);
      %qlkANN(te,ti,tep,tip,q,shear,R0,a0,B0,Amain,pp)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      chii0 = NaN .* ones(101,length(profil.xli));
      chie0 = NaN .* ones(101,length(profil.xli));
      D0    = NaN .* ones(101,length(profil.xli));
      Vt0    = NaN .* ones(101,length(profil.xli));
      Vc0    = NaN .* ones(101,length(profil.xli));
      for l = 1:101
	rlti = -geo.r0 .* tid./ti.* (0.9 + 0.2 .* rand(size(ti)));
	tite = ti./te .* (0.9 + 0.2 .* rand(size(ti)));
	qr = q .* (0.9 + 0.2 .* rand(size(q)));
	sr = s .* (0.9 + 0.2 .* rand(size(s)));
	%
	ter = te  .* (0.9 + 0.2 .* rand(size(te)));
	tir = ti  .* (0.9 + 0.2 .* rand(size(ti)));
        tedr  = pdederive(RLFS,ter,0,2,2,1);
        tidr  = pdederive(RLFS,tir,0,2,2,1);
	%
        switch model_qlk_version
        case '2017'
	    [chie0(l,:),chii0(l,:)] = qlkANN4Dkin(ter',tir',tedr',tidr',qr',sr',geo.r0,geo.a,geo.b0,Amain,data.net);
        otherwise
	    chii0(l,:) = qlkANNk(rlti',qr',tite',sr',data.netief.max_qlk,data.netief.min_qlk,data.netief.net,data.netief.prepros);
	    chie0(l,:) = qlkANNk(rlti',qr',tite',sr',data.neteef.max_qlk,data.neteef.min_qlk,data.neteef.net,data.neteef.prepros);
	end
	D0(l,:)    = qlkANNk(rlti',qr',tite',sr',data.netdfe.max_qlk,data.netdfe.min_qlk,data.netdfe.net,data.netdfe.prepros);
	Vt0(l,:)   = qlkANNk(rlti',qr',tite',sr',data.netvte.max_qlk,data.netvte.min_qlk,data.netvte.net,data.netvte.prepros);
	Vc0(l,:)   = qlkANNk(rlti',qr',tite',sr',data.netvce.max_qlk,data.netvce.min_qlk,data.netvce.net,data.netvce.prepros);
      end
      for l = 1:21
	  c = chii0(:,l);
	  if any(c>0)
	      chii(k,l) = mean(c(c>0));
	      chii_err(k,l) = std(c(c>0));
	  else
	      chii(k,l) = NaN;      
	      chii_err(k,l) = NaN;
	  end
	  c = chie0(:,l);
	  if any(c>0)
	      chie(k,l) = mean(c(c>0));
	      chie_err(k,l) = std(c(c>0));
	  else
	      chie(k,l) = NaN;      
	      chie_err(k,l) = NaN;
	  end
	  c = D0(:,l);
	  if any(c>0)
	      D(k,l) = mean(c(c>0));
	      D_err(k,l) = std(c(c>0));
	  else
	      D(k,l) = NaN;      
	      D_err(k,l) = NaN;
	  end
	  c = Vt0(:,l);
	  if any(c>0)
	      Vt(k,l) = mean(c(c>0));
	      Vt_err(k,l) = std(c(c>0));
	  else
	      Vt(k,l) = NaN;      
	      Vt_err(k,l) = NaN;
	  end
	  c = Vc0(:,l);
	  if any(c>0)
	      Vc(k,l) = mean(c(c>0));
	      Vc_err(k,l) = std(c(c>0));
	  else
	      Vc(k,l) = NaN;      
	      Vc_err(k,l) = NaN;
	  end   
      end
      chifac = (te.*phys.e).^1.5.*sqrt(meff(k).*phys.ua)./(phys.e.^2.*geo.b0.^2.*geo.a);

      % The NN was (mistakingly) fit to unnormalized D and V. We must normalize back to GB using the
      % Te, B0, and a values used when constructing the NN database, and then unnormalize again using
      % the simulation values

      % NN training set values: Te = 8 keV , a = 1 m, B0 = 3 T, Amain = 2, n = 5e19 [m^-3], R/Ln=2
      chifacNN = (8e3.*phys.e).^1.5.*sqrt(2.*phys.ua)./(phys.e.^2.*3.^2.*1);
      switch model_qlk_version
      case '2017'
        % already in SI
      otherwise
	chii(k,:)=(chifac .* chii(k,:)); %convert to SI
	chie(k,:)=(chifac .* chie(k,:)); %convert to SI
      end
      D(k,:)=D(k,:) .* chifac./chifacNN; %convert from SI to GB with NN units, then to SI with simulation units
      V(k,:)=(Vt(k,:)+ Vc(k,:)) .*chifac./chifacNN .* 3./geo.r0;
      switch model_qlk_version
      case '2017'
        % already in SI
      otherwise
	chii_err(k,:)=(chifac .* chii_err(k,:)); %convert to SI
	chie_err(k,:)=(chifac .* chie_err(k,:)); %convert to SI
      end
      D_err(k,:)=D_err(k,:) .* chifac./chifacNN; %convert from SI to GB with NN units, then to SI with simulation units
      V_err(k,:)=(Vt_err(k,:)+ Vc_err(k,:)) .*chifac./chifacNN .* 3./geo.r0;
      
    end

    %filter outputs
    %  chie(chie < chi_min) = chi_min;
    %  chii(chii < chi_min) = chi_min;
    %  D(D < chi_min) = chi_min;
    %  chie(~isfinite(chie)) = 0;
    %  chii(~isfinite(chii)) = 0;
    %  D(~isfinite(D)) = 0;
    %  V(~isfinite(V)) = 0;
    %  chie_err(~isfinite(chie_err)) = 0;
    %  chii_err(~isfinite(chii_err)) = 0;
    %  D_err(~isfinite(D_err)) = 0;
    %  V_err(~isfinite(V_err)) = 0;
end

% neoclassical part (use se same on ion and electron)
% Amain
Amain = meff;
chii_neo = neo_Hinton(profil,option,Amain * ones(size(profil.xli)));

% for Temperature recomputation
chie_ref = profil.xie .* profil.grho2 ./ profil.grho;
chii_ref = profil.xii .* profil.grho2 ./ profil.grho;

coef_min = 1e-3;
% from BgBs model
Dbgbs    = max(coef_min,profil.xie) .* max(coef_min,profil.xii) ./ max(2 .* coef_min,profil.xie + profil.xii) ./ profil.qjli;
Dbgbs(modeh ~= 0,end) = coef_min;
Dbgbs(modeh ~= 0,end - 1) = coef_min;
Dbgbs(modeh ~= 0,end - 2) = 0.5 .* Dbgbs(modeh ~= 0,end - 2) + 0.5 .* coef_min;

ned1  = pdederive(profil.xli,profil.nep,0,2,2,1);
Vware = profil.ware;
ve = ones(size(profil.xli));
warning off
D_ref = - (Vware .* profil.nep + profil.ge ) ./ ned1 .* (profil.rmx(:,end) * ve);
warning on
D_ref = max(Dbgbs,D_ref);
% calcul de V a partir du flux
V_ref = - (profil.ge + D_ref .* ned1 ./ (profil.rmx(:,end) * ve)) ./ profil.nep;
% separation du pinch neoclassique et anormale
V_ref  = V_ref - Vware;
% signe > 0
V_ref  = max(0,V_ref);
% calcul final de D
warning off
D_ref = - ((V_ref + Vware) .* profil.nep + profil.ge) ./ ned1 .* (profil.rmx(:,end) * ve);
D_ref(:,1) = D_ref(:,2);
warning on

V_ref = - V_ref .* profil.grho2 ./ profil.grho;
Vware_ref = - Vware  .* profil.grho2 ./ profil.grho;
% to test numerical integration:
if nargin > 4
  chie = chie_ref;
  chii = chii_ref;
  D    = D_ref;
  V    = V_ref;
end

% Local gradient
rhomax = profil.rmx(:,end) * ve;
% dervative ->rho_tor
ted1  = pdederive(profil.xli,profil.tep,0,2,2,1);
tid1  = pdederive(profil.xli,profil.tip,0,2,2,1);
ned1  = pdederive(profil.xli,profil.nep,0,2,2,1);


% graph
switch mode
case 'full'
  if nargout > 0
    return
  end
end
% recalcul des profils
% factor time dependant
if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') 
    factor_time = ones(size(profil.qjli));
else
    dwthdt      = z0dxdt(wth,profil.temps);
    taue_demi   = 0.5 .* (cat(1,taue(1),taue(1:end-1)) + cat(1,taue(2:end),taue(end)));
    wth_demi   = 0.5 .* (cat(1,wth(1),wth(1:end-1)) + cat(1,wth(2:end),wth(end)));
    rap         = taue_demi .* dwthdt ./ wth_demi;
    rap         = min(9,max(-0.9,rap));
    %figure(21);clf;subplot(2,1,1);plot(profil.temps,rap);
    factor_time = 1 ./ (1 + rap);
    %figure(21);subplot(2,1,2);plot(profil.temps,factor_time);drawnow
    factor_time = factor_time * ones(size(profil.xli));
end

% heat flux are kept as in METIS: not self consistent
Qe    = profil.qe + profil.qei;
Qi    = profil.qi - profil.qei;
Qei   = profil.qei;
% 
ae    = profil.nip ./ max(1e13,profil.nep);

% time independnat particles flux
geoft =   1 ./ max(eps,profil.vpr_tor .* profil.grho2) .*  ...
       cumtrapz(profil.xli, - z0dxdt(profil.nep .* profil.vpr,profil.temps) +   ...
      ((z0dxdt(profil.rmx(:,end),profil.temps) ./ profil.rmx(:,end)) * profil.xli) .* ...
       pdederive(profil.xli,profil.nep .* profil.vpr,0,2,2,1),2);
geoft(:,1) = 0;
ge_sts = max(0,profil.ge - geoft);
%figure(22);plot(profil.xli,ge_sts,'r',profil.xli,profil.ge - geoft,'b.');drawnow

% boucle de convergence
% filter
filter = find((chie_ref(:) <=0 ) | (chii_ref(:) <=0 ) | (D_ref(:) <=0) |  ...
	~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));

chii_ref(filter)=0;
chie_ref(filter)=0;
D_ref(filter)=0;
V_ref(filter)=0;

% limits
chiD_max        = 30;                       % max. value of chi (for convergence)
chi_min         = 0.1;                      % min. value of chi (for convergence)
che_min         = 0.1;                      % min. value of chi (for convergence)
D_min           = 0.3;                      % min. value of D (for convergence)
V_min           = 0.1;                      % min. value of V (for convergence). Positive for inward pinch
V_max           = 300;
if isfield(data.net,'chiemin');
 che_min = data.net.chiemin;
end
if isfield(data.net,'chiemin');
 chi_min = data.net.chiimin;
end
chie_ref(chie_ref <  che_min) = che_min;
chii_ref(chii_ref <  chi_min) = chi_min;
V_ref(D_ref < D_min) = -V_min;
D_ref(D_ref < D_min) = D_min;
chie_ref(chie_ref > chiD_max) = chiD_max;
chii_ref(chii_ref > chiD_max) = chiD_max;
V_ref(V_ref > V_max)          = V_max;
V_ref(V_ref < -V_max)         = -V_max;
D_ref(D_ref > chiD_max)       = chiD_max;

% initialisation
chie_loop = ones(size(profil.tep));
chii_loop = 3 .* ones(size(profil.tep));
D_loop    = ones(size(profil.tep));
V_loop    = Vware_ref;
%
chie_mem = chie_loop;
chii_mem = chii_loop;
D_mem    = D_loop;
V_mem    = Vware_ref;
%
te_loop  = profil.tep;
ti_loop  = profil.tip;
ne_loop  = profil.nep;
%
te_mem  = profil.tep;
ti_mem  = profil.tip;
ne_mem  = profil.nep;
% anti aliasing
f = 0.1;
% mask de fusion
mask  = ones(size(profil.tep));
%mask(:,end - 4)  = 0.66;
%mask(:,end - 3) = 0.33;
mask(:,end - 2)  = 0;
mask(:,end - 1)  = 0;
mask(:,end)  = 0;
ind02 = 5;
switch prediction_domain
case 'Restricted'
    for k=1:length(indice_inv)
      ind = ceil(max(ind02,indice_inv(k)));
      mask(k,1:ind) = 0;
    %  mask(k,ind - 1) = 0.33;
    %  mask(k,ind) = 0.66;
    end
    disp('Call in restricted mode');
end
mask = double(mask);
% boucle de calcul
qr    = profil.qjli;
sr    = NaN .* qr;
for k=1:length(profil.temps)
  sr(k,:)    = pdederive(profil.rmx(k,:),qr(k,:),0,2,2,1) .* profil.rmx(k,:) ./ qr(k,:);
end
switch prediction_domain
case 'Clamp'
  qr(qr < 1)   = 1;
  sr(sr < 0.1) = 0.1;
  disp('Call in clamp mode');
end
RLFS = profil.Raxe .* (1 + profil.epsi);
geo.a  = (profil.Raxe(:,end) .* profil.epsi(:,end)) * ones(size(profil.xli));
geo.r0 = profil.Raxe(:,end) * ones(size(profil.xli));
geo.b0 = (profil.fdia(:,end) ./ profil.Raxe(:,end)) * ones(size(profil.xli));
rhomax = profil.rmx(:,end) * ones(size(profil.xli));
Amain = meff * ones(size(profil.xli));

indedge = length(profil.xli) - 2;
nbloop_max = 501;
recuit = 0;
err_mem = Inf;
err_stab = 1/eps;
for k=1:nbloop_max
  if (k >  (nbloop_max / 3))  || (err_stab > err_mem);
    if recuit == 0
      f = 0.3;
    end
    f = f .* 0.95;
    recuit = 1;
  else
     err_mem = err_stab;
  end
  % bar 
  switch mode
  case 'full'
    waitbar(0.5 + 0.5 .* k ./ nbloop_max,hwaitbar,'Please wait...');
  end  
  % variales pour l'appel des reseaux
  tid  = pdederive(RLFS,ti_loop,0,2,2,1) .* (1 + recuit .* f .* rand(size(RLFS)));
  ted  = pdederive(RLFS,te_loop,0,2,2,1) .* (1 + recuit .* f .* rand(size(RLFS)));
  rlti = -geo.r0 .* tid./ti_loop .* (1 + recuit .* f .* rand(size(RLFS)));;
  tite = ti_loop ./ te_loop .* (1 + recuit .* f .* rand(size(RLFS)));;

  % calcul des coefficients de transport
  switch model_qlk_version
  case '2017'
      [chie0_loc,chii0_loc] = qlkANN4Dkin(te_loop(:),ti_loop(:),ted(:),tid(:),qr(:),sr(:),geo.r0(:),geo.a(:),geo.b0(:),Amain(:),data.net);
      chii_loop = reshape(chii0_loc,size(te_loop));
      chie_loop = reshape(chie0_loc,size(te_loop));
  otherwise
      chii_loop = reshape(qlkANNk(rlti(:),qr(:),tite(:),sr(:),data.netief.max_qlk,data.netief.min_qlk,data.netief.net,data.netief.prepros),size(te_loop));
      chie_loop = reshape(qlkANNk(rlti(:),qr(:),tite(:),sr(:),data.neteef.max_qlk,data.neteef.min_qlk,data.neteef.net,data.neteef.prepros),size(ti_loop));
  end
  D_loop    = reshape(qlkANNk(rlti(:),qr(:),tite(:),sr(:),data.netdfe.max_qlk,data.netdfe.min_qlk,data.netdfe.net,data.netdfe.prepros),size(ne_loop));
  Vt        = reshape(qlkANNk(rlti(:),qr(:),tite(:),sr(:),data.netvte.max_qlk,data.netvte.min_qlk,data.netvte.net,data.netvte.prepros),size(ne_loop));
  Vc        = reshape(qlkANNk(rlti(:),qr(:),tite(:),sr(:),data.netvce.max_qlk,data.netvce.min_qlk,data.netvce.net,data.netvce.prepros),size(ne_loop));

  % turbulent
  chifac = (te_loop .* phys.e) .^ 1.5 .* (sqrt(meff .* phys.ua) * ones(size(profil.xli))) ./ (phys.e .^ 2 .* geo.b0 .^ 2 .* geo.a);

  % The NN was (mistakingly) fit to unnormalized D and V. We must normalize back to GB using the
  % Te, B0, and a values used when constructing the NN database, and then unnormalize again using
  % the simulation values

  % NN training set values: Te = 8 keV , a = 1 m, B0 = 3 T, Amain = 2, n = 5e19 [m^-3], R/Ln=2
  chifacNN = (8e3.*phys.e).^1.5.*sqrt(2.*phys.ua)./(phys.e.^2.*3.^2.*1);

  switch model_qlk_version
  case '2017'
    % already in SI
  otherwise
      chii_loop=(chifac .* chii_loop); %convert to SI
      chie_loop=(chifac .* chie_loop); %convert to SI
  end
  D_loop = D_loop .* chifac./chifacNN; %convert from SI to GB with NN units, then to SI with simulation units
  V_loop=(Vt + Vc) .*chifac./chifacNN .* 3./geo.r0;
  
  % filter
  filter = find((chie_loop(:) <=0 ) | (chii_loop(:) <=0 ) | (D_loop(:) <=0) |  ...
          ~isfinite(chie_loop(:)) | ~isfinite(chii_loop(:)) | ~isfinite(D_loop(:)) | ~isfinite(V_loop(:)));

  chii_loop(filter)=0;
  chie_loop(filter)=0;
  D_loop(filter)=0;
  V_loop(filter)=0;
  
%    % limits
%    chiD_max        = 30;                       % max. value of chi (for convergence)
%    chi_min         = 0.1;                      % min. value of chi (for convergence)
%    D_min           = 0.3;                      % min. value of D (for convergence)
%    V_min           = 0.1;                      % min. value of V (for convergence). Positive for inward pinch
%    V_max           = 300;
  chie_loop(chie_loop <  che_min) = che_min;
  chii_loop(chii_loop <  chi_min) = chi_min;
  V_loop(D_loop < D_min) = -V_min;
  D_loop(D_loop < D_min) = D_min;
  chie_loop(chie_loop > chiD_max) = chiD_max;
  chii_loop(chii_loop > chiD_max) = chiD_max;
  V_loop(V_loop > V_max)          = V_max;
  V_loop(V_loop < -V_max)         = -V_max;
  D_loop(D_loop > chiD_max)       = chiD_max;

  % neoclassical
  % neoclassical part (use se same on ion and electron)
  profil_loc = profil;
  profil_loc.tep = te_loop;
  profil_loc.tip = ti_loop;
  chii_neo  = neo_Hinton(profil_loc,option,Amain);
  chie_neo  = chii_neo .* profil.grho2 ./ profil.grho;
  VWare_neo = Vware_ref; % not updated in this first version

  % somme 
  chii_loop = chii_loop + chii_neo;
  V_loop    = V_loop + VWare_neo;
  % space filter
  for lk = 1:length(profil.temps)
    chie_loop(lk,:) = sgolayfilt(chie_loop(lk,:),1,3);
    chii_loop(lk,:) = sgolayfilt(chii_loop(lk,:),1,3);
    D_loop(lk,:)    = sgolayfilt(D_loop(lk,:),1,3);
    V_loop(lk,:)    = sgolayfilt(V_loop(lk,:),1,3);
  end
  
  % masquage 
  chie_loop = (1 - mask) .* chie_ref + mask .* chie_loop .* factor_time;
  chii_loop = (1 - mask) .* chii_ref + mask .* chii_loop .* factor_time;
  D_loop    = (1 - mask) .* D_ref    + mask .* D_loop .* factor_time;
  V_loop    = (1 - mask) .* V_ref    + mask .* V_loop .* factor_time;  
  
  % amorti
  chie_loop = (1 - f) .* chie_mem + f .* chie_loop;
  chii_loop = (1 - f) .* chii_mem + f .* chii_loop;
  D_loop    = (1 - f) .* D_mem    + f .* D_loop;
  V_loop    = (1 - f) .* V_mem    + f .* V_loop;  
  
  % calcul des profils 
  warning off
  % with fixed sources including Qei
  ted1_loop    = - (Qe - Qei) ./ (chie_loop .* ne_loop .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
  tid1_loop    = - (Qi + Qei) ./ (chii_loop .* ne_loop .* ae .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
  % flux equation not output (convention Qualikiz !)
  ge_loop      = ge_sts - V_loop .* ne_loop .* profil.grho ./ profil.grho2;
  ned1_loop    = - ge_loop ./ D_loop .* rhomax;
  warning on

  % maximum gradient
  tdmax =  pi .* max(max(abs(ted1(:))),max(abs(tid1(:))));
  ndmax =  pi .* max(abs(ned1(:)));
  ted1_loop = max(-tdmax,min(tdmax,ted1_loop));
  ned1_loop = max(-ndmax,min(ndmax,ned1_loop));
  tid1_loop = max(-tdmax,min(tdmax,tid1_loop));
    
 
  % security
  ted1_loop(~isfinite(ted1_loop)) =ted1(~isfinite(ted1_loop));
  tid1_loop(~isfinite(tid1_loop)) =tid1(~isfinite(tid1_loop));
  ned1_loop(~isfinite(ned1_loop)) =ned1(~isfinite(ned1_loop));
  
  % central part
  ted1_loop = mask .* ted1_loop + (1 - mask) .* ted1;
  tid1_loop = mask .* tid1_loop + (1 - mask) .* tid1;
  ned1_loop = mask .* ned1_loop + (1 - mask) .* ned1;
  
  % symetry
  ted1_loop(:,1) = 0;
  tid1_loop(:,1) = 0;
  ned1_loop(:,1) = 0;
 
  % integration
  te_loop = cumtrapz(profil.xli(end:-1:1),ted1_loop(:,end:-1:1),2);
  te_loop = te_loop(:,end:-1:1);
  ti_loop = cumtrapz(profil.xli(end:-1:1),tid1_loop(:,end:-1:1),2);
  ti_loop = ti_loop(:,end:-1:1);
  ne_loop = cumtrapz(profil.xli(end:-1:1),ned1_loop(:,end:-1:1),2);
  ne_loop = ne_loop(:,end:-1:1);

  % raccord
  te_loop = te_loop + (profil.tep(:,indedge) - te_loop(:,indedge)) * ones(size(profil.xli));
  te_loop(:,indedge:end) = profil.tep(:,indedge:end);
  ti_loop = ti_loop + (profil.tip(:,indedge) - ti_loop(:,indedge)) * ones(size(profil.xli));
  ti_loop(:,indedge:end) = profil.tip(:,indedge:end);
  ne_loop = ne_loop + (profil.nep(:,indedge) - ne_loop(:,indedge)) * ones(size(profil.xli));
  ne_loop(:,indedge:end) = profil.nep(:,indedge:end);
  
  % valeur minimum
  te_loop = max(te_loop,profil.tep(:,end) * ones(size(profil.xli)));
  ti_loop = max(ti_loop,profil.tip(:,end) * ones(size(profil.xli)));
  ne_loop = max(ne_loop,profil.nep(:,end) * ones(size(profil.xli)));
  
  % amorti
  te_loop = (1 - f) .* te_mem + f .* te_loop;
  ti_loop = (1 - f) .* ti_mem + f .* ti_loop;
  switch model_qlk_version
  case '2017'
    ne_loop = ne_mem;
  otherwise
    ne_loop = (1 - f) .* ne_mem + f .* ne_loop;
  end
  % new Qei value
  nHp_loop  = nHp  .*  ne_loop ./ max(1,profil.nep);
  nDp_loop  = nDp  .*  ne_loop ./ max(1,profil.nep);
  nTp_loop  = nTp  .*  ne_loop ./ max(1,profil.nep);
  nhep_loop = nhep .*  ne_loop ./ max(1,profil.nep);
  nzp_loop  = nzp  .*  ne_loop ./ max(1,profil.nep);
  nwp_loop  = nwp  .*  ne_loop ./ max(1,profil.nep);
  % terme source
  switch  extended_qei
      case 'on'
          qeib0 = equipartition_full(te_loop,ti_loop,ne_loop,nHp_loop,nDp_loop,nTp_loop,nhep_loop,nzp_loop, ...
              option.rimp .* nzp,(1 - Sn_fraction) .* nwp_loop, Sn_fraction .* nwp_loop,option.zimp,option.zmax);
      otherwise
          
          warning off
          lnei          =  14.9 - 0.5.*log(ne_loop ./ 1e20) + log(te_loop ./ 1e3);
          warning on
          ind = find(~isfinite(lnei) | (lnei <10));
          if ~isempty(ind)
              lnei(ind) = 10 .* ones(1,length(ind));
          end
          taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
              ((phys.e .* max(min(profil.tep(:)),te_loop)) .^ (3/2) ./ ne_loop ./ lnei);
          if Sn_fraction > 0
              qeib0    = 3 .* phys.me ./ phys.mp ./ taues .* (nHp_loop + nDp_loop ./ 2 + nTp_loop ./ 3 + nhep_loop  +  ...
                  (option.zimp./ (7/3)) .* nzp_loop + option.rimp .*  (option.zmax./ (7/3)) .* nzp_loop +  ...
                  (1 - Sn_fraction) .*  z0wavez(te_loop) ./ 183.84 .* nwp_loop +  Sn_fraction .*  z0snavez(te_loop) ./ 118.71 .* nwp_loop);
              
          else
              qeib0    = 3 .* phys.me ./ phys.mp ./ taues .* (nHp_loop + nDp_loop ./ 2 + nTp_loop ./ 3 + nhep_loop  +  ...
                  (option.zimp./ (7/3)) .* nzp_loop + option.rimp .*  (option.zmax./ (7/3)) .* nzp_loop + z0wavez(te_loop) ./ 183.84 .* nwp_loop);
          end
  end
  qeib     = qeib0 .* phys.e .* (te_loop - ti_loop);
  Qei = cumtrapz(profil.xli,qeib .* profil.vpr,2);

  % calcul de la stabilite
  err_stab = std(chie_loop(:) - chie_mem(:)) + std(chii_loop(:) - chii_mem(:)) + std(D_loop(:) - D_mem(:)) + std(V_loop(:) - V_mem(:)) + ...
             std(te_loop(:) - te_mem(:)) ./ 1e3 + std(ti_loop(:) - ti_mem(:)) ./ 1e3 + std(ne_loop(:) - ne_mem(:)) ./ 1e18;
  if ~isfinite(err_stab)
    keyboard
  end
  % condition de sortie
  if err_stab < max(1e-3,option.tol0d)
      break
  elseif recuit
      fprintf('err = %g @ Temp = %g\n',err_stab,f);
  else 
       fprintf('err = %g\n',err_stab); 
  end
  % recopy dans mem
  chie_mem = chie_loop;
  chii_mem = chii_loop;
  D_mem    = D_loop;
  V_mem    = V_loop;
  te_mem   = te_loop;
  ti_mem   = ti_loop;
  ne_mem   = ne_loop;
  
   
%    figure(21);
%    subplot(2,2,1)
%    plot(profil.xli,ted1_loop,'r',profil.xli,ted1,'b');
%    subplot(2,2,2)
%    plot(profil.xli,tid1_loop,'r',profil.xli,tid1,'b');
%    subplot(2,2,3)
%    plot(profil.xli,ned1_loop,'r',profil.xli,ned1,'b');
%    subplot(2,2,4)
%    plot(profil.xli,Qei,'r',profil.xli,profil.qei,'b');
%    drawnow

end
% fin de la boucle
fprintf('converged in %d loops\n',k);
% remove wait bar
switch mode
case 'full'
  delete(hwaitbar);
end
mask_err = mask;
mask(mask == 0) = NaN;
mask(mask > 0.5) = 0;

% output flux
flux_qe = - chie .* profil.nep .* ted1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
flux_qe_err = - chie_err .* profil.nep .* ted1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
flux_qi = - chii .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
flux_qi_err = - chii_err .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
flux_qineo = - chii_neo .* profil.nip .* tid1 .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
% flux equation not output (convention Qualikiz !)
flux_ne        = - D .* ned1 ./ rhomax + V .* profil.nep .* profil.grho ./ profil.grho2;
flux_ne_err    = abs(D_err .* ned1) ./ rhomax + abs(V_err .* profil.nep .* profil.grho ./ profil.grho2);

% energy content
wth_nn = (3/2) .* trapz(profil.xli,phys.e .* ne_loop .* (te_loop + ti_loop .* profil.nip ./ max(1,profil.nep)) .* profil.vpr,2);

% export data for external data mecanism of METIS
switch mode
case 'full'
  QLK_nn_data.time = profil.temps;
  QLK_nn_data.x    = profil.xli;
  QLK_nn_data.TE   = te_loop;
  QLK_nn_data.TI   = ti_loop;
  QLK_nn_data.NE   = ne_loop;
otherwise
    if evalin('base','exist(''QLK_nn_data'')')
	 err_press_mem = evalin('base','QLK_nn_data.err');
    else
         err_press_mem = Inf;
    end
    QLK_nn_data.time = profil.temps;
    QLK_nn_data.x    = profil.xli;
    QLK_nn_data.err  = sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2));
    [err_press,ind_best] = min(sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2)));
    QLK_nn_data.err  = err_press;
    vt = ones(size(profil.temps));
    %disp([err_press ,err_press_mem])
    if err_press > err_press_mem
      te_exp = 0.1 .* te_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.9 .* profil.tep(ind_best,:) .* profil.nep(ind_best,:);
      ti_exp = 0.1 .* ti_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.9 .* profil.tip(ind_best,:) .* profil.nep(ind_best,:);
      ne_exp = 0.1 .* ne_loop(ind_best,:) + 0.9 .* profil.nep(ind_best,:);
    else
      te_exp = 0.3 .* te_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.7 .* profil.tep(ind_best,:) .* profil.nep(ind_best,:);
      ti_exp = 0.3 .* ti_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.7 .* profil.tip(ind_best,:) .* profil.nep(ind_best,:);
      ne_exp = 0.3 .* ne_loop(ind_best,:) + 0.7 .* profil.nep(ind_best,:);
    end
    QLK_nn_data.TE   = vt * (te_exp ./ ne_exp);
    QLK_nn_data.TI   = vt * (ti_exp ./ ne_exp);
    QLK_nn_data.NE   = vt * ne_exp;
    QLK_nn_data.N_bar = vt * trapz(profil.xli,ne_exp);
end
zassignin('base','QLK_nn_data',QLK_nn_data);

switch mode
case 'full'

      hz =findobj(0,'type','figure','tag','flux_net_kin');
      if isempty(hz)
		hz=figure('tag','flux_net_kin','name','Flux comparison (ITG part only)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1])

      subplot(2,2,1)
      zplotprof(gca,profil.temps,profil.xli,profil.qe,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_qe,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_qe -  flux_qe_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_qe +  flux_qe_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_qe (W/m^{-3})');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,2)
      zplotprof(gca,profil.temps,profil.xli,profil.qi,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_qi,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_qi -  flux_qi_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_qi -  flux_qi_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_qi (W/m^{-3})');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,3)
      zplotprof(gca,profil.temps,profil.xli,profil.ge,'color','b');
      zplotprof(gca,profil.temps,profil.xli,flux_ne,'color','r');
      zplotprof(gca,profil.temps,profil.xli,(flux_ne - flux_ne_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(flux_ne + flux_ne_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('flux_ne (e^-/m^{-3}/s)');
      legend('METIS','QLKNN-4Dkin');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      hz =findobj(0,'type','figure','tag','Chi_net_kin');
      if isempty(hz)
		hz=figure('tag','Chi_net_kin','name','Transport coefficients comparison (ITG part only)');
      else
		figure(hz);
      end
      clf
      set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	      'defaultlinelinewidth',3,'color',[1 1 1])

      subplot(2,2,1)
      zplotprof(gca,profil.temps,profil.xli,chie_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,chie,'color','r');
      zplotprof(gca,profil.temps,profil.xli,chie_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(chie-chie_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(chie+chie_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('Chi_e (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,2)
      zplotprof(gca,profil.temps,profil.xli,chii_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,chii,'color','r');
      zplotprof(gca,profil.temps,profil.xli,chii_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,chii_neo,'color','g');
      zplotprof(gca,profil.temps,profil.xli,(chii-chii_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(chii+chii_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('Chi_i (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction','Neo');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,3)
      zplotprof(gca,profil.temps,profil.xli,D_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,D,'color','r');
      zplotprof(gca,profil.temps,profil.xli,D_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(D-D_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(D+D_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('D (m^2/s)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon

      subplot(2,2,4)
      VsD_err = abs(V_err ./ D) + abs(V.* D_err ./ D .^ 2);
      zplotprof(gca,profil.temps,profil.xli,V_ref ./ D_ref,'color','b');
      zplotprof(gca,profil.temps,profil.xli,V ./ D,'color','r');
      zplotprof(gca,profil.temps,profil.xli,V_loop ./ D_loop,'color','k');
      zplotprof(gca,profil.temps,profil.xli,(V ./ D - VsD_err),'color','m');
      zplotprof(gca,profil.temps,profil.xli,(V ./ D + VsD_err),'color','m');
      set(gca,'YScale','linear');
      ylabel('V / D (m)');
      legend('METIS','QLKNN-4Dkin with METIS profiles','QLKNN-4Dkin prediction');
      xlabel('r/a');
      %z0loglin(gca);
      zoom yon
end

hz =findobj(0,'type','figure','tag','Prof_net_kin');
if isempty(hz)
  	  hz=figure('tag','Prof_net_kin','name','Profiles prediction qlkAnn_kin_e (ITG part only, fixed sources but equipartition, no time dependance)');
else
  	  figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)
zplotprof(gca,profil.temps,profil.xli,profil.tep ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,te_loop ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('T_e (keV)');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.tip ./ 1e3,'color','b');
zplotprof(gca,profil.temps,profil.xli,ti_loop ./ 1e3,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('T_i (keV)');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,profil.nep ./ 1e19,'color','b');
zplotprof(gca,profil.temps,profil.xli,ne_loop ./ 1e19,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('n_e (10^{19} m^{-3})');
%legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,4)
zplotprof(gca,profil.temps,profil.xli,profil.qei ./ 1e6,'color','b');
zplotprof(gca,profil.temps,profil.xli, Qei ./ 1e6,'color','r');
zplotprof(gca,profil.temps,profil.xli,mask,'color','g');
set(gca,'YScale','linear');
ylabel('Q_{e,i} (MW)');
legend('METIS','QLKNN-4Dkin','interval of computation with QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

switch mode
case 'full'
    hz =findobj(0,'type','figure','tag','Wth_net_kin');
    if isempty(hz)
	      hz=figure('tag','Wth_net_kin','name','Energy content prediction qlkAnn_kin_e (ITG part only, fixed sources but equipartition, no time dependance)');
    else
	      figure(hz);
    end
    clf
    set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',3,'color',[1 1 1])

    subplot(2,1,1)
    plot(profil.temps,wth ./ 1e6,'b',profil.temps,wth_nn ./ 1e6,'r');
    xlabel('time (s)');
    ylabel('W_{th} (MJ)');
    legend('METIS','QLKNN-4Dkin');
    zoom on
    subplot(2,1,2);
    semilogy(profil.temps,wth_nn ./ max(1,wth),'b',profil.temps,factor_time,'r',profil.temps,ones(size(profil.temps)),'g');
    xlabel('time (s)');
    legend('W_{th,QLKNN-4Dkin} / W_{th,METIS}','time factor');
    z0loglin(gca);
    zoom on
    joint_axes(hz,2);
end

% backup internal data for future use
tempf = tempname;
clear hz
save(tempf)
zassignin('base','internal_data_QLK',load(tempf)); 
delete(sprintf('%s.mat',tempf));

function te_out = compute_temp(xli,ted1_out,tep);

ve     = ones(size(xli));

te_out = fliplr(cumtrapz(fliplr(xli),fliplr(ted1_out),2));
te_out = te_out  + (tep(:,end - 1) - te_out(:,end-1)) * ve;
te_out(:,end) = tep(:,end);
te_out(:,1) = te_out(:,2);


%ref: Wesson
function chii_neo = neo_Hinton(profil,option,Amain)

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
    alpha  = nI_ZI2 ./ profil.nhep;
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
           
