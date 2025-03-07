% Attention sorties = definitions Qualikiz !
%[chie,chii,chii_neo,D,V] = z0qlkANN_kin_e_2018(post.z0dinput.option,post.zerod,post.profil0d); 
function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0qlkANN_kin_e_2018(option,zerod,profil)

if isfield(option,'Sn_fraction') && (option.Sn_fraction > 0)
    error('Sn is not taken into account by this function (option.Sn_fraction > 0)');
end
if isfield(option,'extended_qei')
    switch  option.extended_qei
        case 'on'
            error('complete collisional heat transfert is not taken into account by this function (option.extended_qei == ''on'')');
    end
end

% init output
chie_loop =[];
chii_loop =[];
chii_neo  =[];
D_loop    =[];
V_loop    =[];

% ask for parameters
first_call = 1;
if isfield(option,'seuil')
  first_call = 0;
else
  prompt={'shift_chi_e:','shift_chi_i:','shift_DV:','threshold:','neo_ion_in_D:','neo_ion_in_Chie:','factor_neo_chii:'};
  name='QLKNN parameters';
  numlines=1;
  defaultanswer={'0','0','0','0','0','0','1'};
 
  answer=inputdlg(prompt,name,numlines,defaultanswer);
  if isempty(answer)
    return
  end
  option.shift_chi_e     = str2num(answer{1});
  option.shift_chi_i     = str2num(answer{2});
  option.shift_DV        = str2num(answer{3});
  option.seuil           = str2num(answer{4});
  option.neo_ion_in_D    = str2num(answer{5});
  option.neo_ion_in_Chie = str2num(answer{6});
  option.factor_neo_chii = str2num(answer{7});
  first_call = 1;
end                         
prediction_domain = 'Extrapolated';

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
w    = interp1(zerod.temps,zerod.w,profil.temps,'linear','extrap');
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

% loading NN data
if isdeployed
  net    = load('kinetic_e_5D_ITG_20180109.mat');
else
  net    = load(fullfile(fileparts(which('qlkANN.m')),'neuralnet_repo','kinetic_e_5D_ITG_20180109'));  
end

data.net    = net.net;
data.net.constrainrlti = true;
data.net.constraintite = true;
data.net.constrainq = true;
data.net.constrains = true;
data.net.rltimin     = 2;
data.net.rltimax     = 12;
data.net.titemin     = 0.3;
data.net.titemax     = 3.0;
data.net.qmin        = 0.5;
data.net.qmax        = 5;
data.net.shearmin    = 0.1;
data.net.shearmax    = 10;
data.net.shift_chi_e = option.shift_chi_e;
data.net.shift_chi_i = option.shift_chi_i;
data.net.shift_DV    = option.shift_DV;
data.net.Dmin        = 0.001;
data.net.Vmin        = 0;
data.net.chiemin     = 0.001;
data.net.chiimin     = 0.001;
data.net.allowpospinch = false;
data.net.smooth_space  = true;

% display choosen parameters
if first_call == 1
    disp('=============================================================')
    fprintf('This tool use a prototype version of neural network interpolation of QUALIKIZ;\nplease consider results with waryness\n');
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
    fprintf('shift_chi_i =  %g\n',data.net.shift_chi_i);
    fprintf('shift_DV =  %g\n',data.net.shift_DV);
    fprintf('chiemin =  %g\n',data.net.chiemin);
    fprintf('chiimin =  %g\n',data.net.chiimin);
    fprintf('Dmin =  %g\n',data.net.Dmin);
    fprintf('Vmin =  %g\n',data.net.Vmin);
    fprintf('allowpospinch =  %g\n',data.net.allowpospinch);
    fprintf('smooth_space =  %g\n',data.net.smooth_space);
    disp('=============================================================')
    disp('');
    hdisc = warndlg('This tool use a prototype version of neural network interpolation of QUALIKIZ; please consider results with waryness','Disclamer');
else
    hdisc = [];
end

if first_call == 1
  hwaitbar = waitbar(0,'Please wait...');
end

% Initial transport coefficient to start the convergence
% neoclassical part (use se same on ion and electron)
% Amain
Amain = meff;
chii_neo = option.factor_neo_chii .* neo_Hinton(profil,option,Amain * ones(size(profil.xli)));

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
Vware = profil.ware;   % have the same sign than in CRONOS
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
Vware_ref = - profil.ware  .* profil.grho2 ./ profil.grho;
%figure(21);plot(profil.xli,profil.grho2 ./ profil.grho);drawnow

% rhomax definition (2D)
rhomax = profil.rmx(:,end) * ve;

% electromagnetic correction
% for benchmark
if isfield(profil,'chii_neo')
    factor_em = ones(size(profil.qjli));
else
    factor_em = min(1,wth ./ w) * ones(size(profil.xli));
end
%figure(21);plot(profil.temps,factor_em);drawnow

% recalcul des profils
% factor time dependant
if isappdata(0,'TE_EXP') && isappdata(0,'TI_EXP') 
    factor_time = ones(size(profil.qjli));
else
    % previous version
%      dwthdt      = z0dxdt(wth,profil.temps);
%      taue_demi   = 0.5 .* (cat(1,taue(1),taue(1:end-1)) + cat(1,taue(2:end),taue(end)));
%      wth_demi   = 0.5 .* (cat(1,wth(1),wth(1:end-1)) + cat(1,wth(2:end),wth(end)));
%      rap         = taue_demi .* dwthdt ./ wth_demi;
%      rap         = min(9,max(-0.9,rap));
%      %figure(21);clf;subplot(2,1,1);plot(profil.temps,rap);
%      factor_time = 1 ./ (1 + rap);
%      %figure(21);subplot(2,1,2);plot(profil.temps,factor_time);drawnow
%      factor_time = factor_time * ones(size(profil.xli));

    % new version
    dwthdt      = z0dxdt(wth,profil.temps);
    power       = max(cumtrapz(profil.xli,(profil.qe + profil.qi) .* profil.vpr,2),[],2);
    factor_time = sgolayfilt((power - dwthdt) ./ max(1,power),1,3);
    indbad      = find((factor_time < 0.1) |(factor_time > 10));
    indok       = find((factor_time >= 0.1) & (factor_time <= 10));
    if length(indok) > 2
      factor_time(indbad) = interp1(profil.temps(indok),factor_time(indok),profil.temps(indbad),'linear',1);
    else
      factor_time = ones(size(factor_time));
    end
    factor_time = factor_time * ones(size(profil.xli));
end
% heat flux are kept as in METIS: not self consistent
Qe    = (profil.qe + profil.qei) .* factor_time;
Qi    = (profil.qi - profil.qei) .* factor_time;
Qei   = profil.qei .* factor_time;
% 
ae    = profil.nip ./ max(1e13,profil.nep);

% time independant particles flux
%  ge =   1 ./ max(eps,profli.vpr_tor .* profli.grho2) .*  ...
%         cumtrapz(profli.xli, profli.vpr .* stot - z0dxdt(nep .* profli.vpr,zs.temps) +   ...
%        ((z0dxdt(profli.rmx(:,end),zs.temps) ./ profli.rmx(:,end)) * profli.xli) .* ...
%         pdederive(profli.xli,nep .* profli.vpr,0,2,2,1),2);
%  

geoft =   1 ./ max(eps,profil.vpr_tor .* profil.grho2) .*  ...
       cumtrapz(profil.xli, - z0dxdt(profil.nep .* profil.vpr,profil.temps) +   ...
      ((z0dxdt(profil.rmx(:,end),profil.temps) ./ profil.rmx(:,end)) * profil.xli) .* ...
       pdederive(profil.xli,profil.nep .* profil.vpr,0,2,2,1),2);
geoft(:,1) = 0;
% for benchmark
if isfield(profil,'chii_neo')
    ge_sts = profil.ge;
else
    ge_sts = profil.ge - geoft;
end
ge_sts = ge_sts .* profil.grho2 ./ profil.grho;
%figure(22);plot(profil.xli,ge_sts,'r',profil.xli,profil.ge - geoft,'b.');drawnow
gi_sts = profil.nip ./ profil.nep .* ge_sts;

% boucle de convergence

% limits
chiD_max        = 30;                       % max. value of chi (for convergence)
chi_min         = 0.1;                      % min. value of chi (for convergence)
che_min         = 0.1;                      % min. value of chi (for convergence)
D_min           = 0.01;                      % min. value of D (for convergence)
V_min           = 0;                      % min. value of V (for convergence). Positive for inward pinch
V_max           = 300;
if isfield(data.net,'chiemin');
 che_min = data.net.chiemin;
end
if isfield(data.net,'chiimin');
 chi_min = data.net.chiimin;
end
if isfield(data.net,'Dmin');
 D_min = data.net.Dmin;
end
if isfield(data.net,'Vmin');
 V_min = data.net.Vmin;
end

%filter = find((chie_ref(:) <=0 ) | (chii_ref(:) <=0 ) | (D_ref(:) <=0) |  ...
%	~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
filter = find(~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
chii_ref(filter)=0;
chie_ref(filter)=0;
D_ref(filter)=0;
V_ref(filter)=0;

chie_ref(chie_ref <  che_min) = che_min;
chii_ref(chii_ref <  chi_min) = chi_min;
V_ref(V_ref > -V_min) = -V_min;
D_ref(D_ref < D_min) = D_min;
chie_ref(chie_ref > chiD_max) = chiD_max;
chii_ref(chii_ref > chiD_max) = chiD_max;
V_ref(V_ref > V_max)          = V_max;
V_ref(V_ref < -V_max)         = -V_max;
D_ref(D_ref > chiD_max)       = chiD_max;

% initialisation
chie_loop = ones(size(profil.tep));
chii_loop = 3 .* ones(size(profil.tep));
D_loop    = ones(size(profil.tep)) / 3;
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
%mask(:,end - 2)  = 0;
mask(:,end - 1)  = 0;
mask(:,end)  = 0;
indedge = length(profil.xli) - 1;
%indedge = length(profil.xli) - 2;
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

% Local gradient (for maximum allows value)
ted1  = pdederive(profil.xli,profil.tep,0,2,2,1);
tid1  = pdederive(profil.xli,profil.tip,0,2,2,1);
ned1  = pdederive(profil.xli,profil.nep,0,2,2,1);
tdmax =  pi .* max(max(abs(ted1(:))),max(abs(tid1(:))));
ndmax =  pi .* max(abs(ned1(:)));

nbloop_max = 501;
recuit = 1e-6;
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
  if first_call == 1
     waitbar(k ./ nbloop_max,hwaitbar,'Please wait...');
  end
  % variales pour l'appel des reseaux
  tid  = factor_em .* pdederive(RLFS,ti_loop,0,2,2,1) .* (1 + recuit .* f .* rand(size(RLFS)));
  ted  = pdederive(RLFS,te_loop,0,2,2,1) .* (1 + recuit .* f .* rand(size(RLFS)));
  ned  = pdederive(RLFS,ne_loop,0,2,2,1) .* (1 + recuit .* f .* rand(size(RLFS)));
  rlti = -geo.r0 .* tid./ti_loop .* (1 + recuit .* f .* rand(size(RLFS)));;
  tite = ti_loop ./ te_loop .* (1 + recuit .* f .* rand(size(RLFS)));;

  % calcul des coefficients de transport
  % qlkANN4Dkin(te,ti,tepa,tipa,q,shear,R0,a0,B0,Amain,pp,transflag)
  [chie0_loc,chii0_loc] = qlkANN4Dkin(te_loop(:),ti_loop(:),ted(:),tid(:),qr(:),sr(:),geo.r0(:),geo.a(:),geo.b0(:),Amain(:),data.net,1);
  chii_loop = reshape(chii0_loc,size(te_loop));
  chie_loop = reshape(chie0_loc,size(te_loop));
    
  [D_loc,V_loc] = qlkANN4Dkin(te_loop(:),ti_loop(:),ted(:),tid(:),qr(:),sr(:),geo.r0(:),geo.a(:),geo.b0(:),Amain(:),data.net,2);
  D_loop = reshape(D_loc,size(te_loop));
  V_loop = reshape(V_loc,size(te_loop));    
   
   
%    figure(21);
%    subplot(2,2,1)
%    plot(profil.xli,chie_loop);
%    subplot(2,2,2)
%    plot(profil.xli,chii_loop);
%    subplot(2,2,3)
%    plot(profil.xli,D_loop);
%    subplot(2,2,4)
%    plot(profil.xli,V_loop);
%    drawnow
  
  
  % boucle sur les shift
  if (option.seuil > 0) 
      seuil = option.seuil;
      
      % flux normalized
      flux_qe = abs(chie_loop .* ted ./ te_loop); 
      flux_qi = abs(chii_loop .* ae .* tid ./ ti_loop); 
      flux_ne = abs(D_loop .* ned ./ ne_loop);

      % points � recalculer 
      indi = find((flux_qe(:) ./ flux_qi(:)) > seuil);
      % correction i
      if ~isempty(indi)
           chii_loop(indi) = min(chiD_max,abs(chie_loop(indi) .* ted(indi) ./ tid(indi) .* ti_loop(indi) ./ te_loop(indi) ./  ae(indi)) ./ seuil);
      end
      
      % flux normalized
      flux_qe = abs(chie_loop .* ted ./ te_loop); 
      flux_qi = abs(chii_loop .* ae .* tid ./ ti_loop); 
      flux_ne = abs(D_loop .* ned ./ ne_loop);

      % correction e
      inde = find((flux_qi(:) ./ flux_qe(:)) > seuil);
      if ~isempty(inde)
           chie_loop(inde) = min(chiD_max,abs(chii_loop(inde) .* tid(inde) ./ ted(inde) .* te_loop(inde) ./ ti_loop(inde) .*  ae(inde)) ./ seuil);
      end
      
      % flux normalized
      flux_qe = abs(chie_loop .* ted ./ te_loop); 
      flux_qi = abs(chii_loop .* ae .* tid ./ ti_loop); 
      flux_ne = abs(D_loop .* ned ./ ne_loop);

      % correction n
      indn = find((max(flux_qi(:),flux_qe(:)) ./ flux_ne(:)) > seuil);

      if ~isempty(indn)
           fluxn = max(flux_qi(indn),flux_qe(indn));
	   D_loop(indn) = min(chiD_max,abs(fluxn ./ ne_loop(indn) .* ned(indn) ./ seuil));
      end
      
      nbcor = length(indi) + length(inde) + length(indn);
      nbi = length(te_loop(:));
	  
   else
	  nbi       = length(te_loop(:));
          nbcor     = 0;
   end
   %disp([shift,length(inde) ./ nbi ,length(indi) ./ nbi ,length(indn) ./ nbi]);

%    figure(21)
%    %plot(ted(:) ./ te_loop(:),flux_qe(:),'.r',tid(:) ./ ti_loop(:),flux_qi(:),'.b',ned(:) ./ ne_loop(:),flux_ne(:),'.g')
%    subplot(2,2,1)
%    hist(flux_qe(:) ./ flux_qi(:),0:0.1:10)
%    set(gca,'xlim',[0,10]);
%    subplot(2,2,2)
%    hist(flux_qi(:) ./ flux_qe(:),0:0.1:10)
%    set(gca,'xlim',[0,10]);
%    subplot(2,2,3)
%    hist(flux_ne(:) ./ max(flux_qi(:),flux_qe(:)),0:0.1:10)
%    set(gca,'xlim',[0,10]);
%    subplot(2,2,4)
%    hist( min(flux_qi(:),flux_qe(:)) ./ flux_ne(:),0:0.1:10)
%    set(gca,'xlim',[0,10]);
%    drawnow
  
  
  % how to evaluate the shift ?
  % filter
  filter = ~isfinite(chie_loop(:)) | ~isfinite(chii_loop(:)) | ~isfinite(D_loop(:)) | ~isfinite(V_loop(:));

  if any(filter)
    disp('bad values in NN results');
  end

  chii_loop(filter) = chie_mem(filter);
  chie_loop(filter) = chii_mem(filter);
  D_loop(filter)    = V_mem(filter);
  V_loop(filter)    = D_mem(filter);
  
  test = (chie_loop > chiD_max) | (chii_loop > chiD_max) | (D_loop > chiD_max) | (abs(V_loop) > V_max);
  if any(test(:))
      overshoot = sum(double(test(:)))./ numel(test)
  end
  chie_loop(chie_loop > chiD_max) = chie_mem(chie_loop > chiD_max);
  chii_loop(chii_loop > chiD_max) = chie_mem(chii_loop > chiD_max);
  D_loop(D_loop > chiD_max)       = D_mem(D_loop > chiD_max);
  V_loop(V_loop > V_max)          = V_mem(V_loop > V_max);
  V_loop(V_loop < -V_max)         = V_mem(V_loop < -V_max);
  %test = ((chie_loop < che_min) |(chii_loop < chi_min) |(D_loop < D_min));
  %sum(double(test(:)))./ numel(test)
  
  chie_loop(chie_loop < che_min)  = che_min;
  chii_loop(chii_loop < chi_min)  = chi_min;
  V_loop(V_loop > -V_min)     = -V_min;
  D_loop(D_loop < D_min)     = D_min;
  
  % neoclassical
  % neoclassical part (use se same on ion and electron)
  profil_loc = profil;
  profil_loc.tep = te_loop;
  profil_loc.tip = ti_loop;
  switch option.gaz
      case 5
          iso = zerod.nhem./ zerod.n1m;
      case 11
          iso = zerod.nTm ./ zerod.n1m;
      otherwise
          iso = zeros(size(zerod.nem));
  end
  chii_neo  = option.factor_neo_chii .* neo_Hinton(profil_loc,option,Amain,iso);
  chie_neo  = chii_neo .* profil.grho2 ./ profil.grho;
  VWare_neo = Vware_ref; % not updated in this first version

  % somme 
  chii_loop = chii_loop + chii_neo;
  V_loop    = V_loop + VWare_neo;
  D_loop    = max(D_loop,option.neo_ion_in_D .* chii_neo);
  chie_loop = max(chie_loop,option.neo_ion_in_Chie .* chii_neo);
  
  % space filter
  if data.net.smooth_space
      for lk = 1:length(profil.temps)
	chie_loop(lk,:) = sgolayfilt(chie_loop(lk,:),1,3);
	chii_loop(lk,:) = sgolayfilt(chii_loop(lk,:),1,3);
	D_loop(lk,:)    = sgolayfilt(D_loop(lk,:),1,3);
	V_loop(lk,:)    = sgolayfilt(V_loop(lk,:),1,3);
      end
  end

  % masquage  (only evolved part of the profile is changed)
  chie_loop = (1 - mask) .* chie_ref + mask .* chie_loop;
  chii_loop = (1 - mask) .* chii_ref + mask .* chii_loop;
  D_loop    = (1 - mask) .* D_ref    + mask .* D_loop;
  V_loop    = (1 - mask) .* V_ref    + mask .* V_loop;  
  
  % amorti
  chie_loop = (1 - f) .* chie_mem + f .* chie_loop;
  chii_loop = (1 - f) .* chii_mem + f .* chii_loop;
  D_loop    = (1 - f) .* D_mem    + f .* D_loop;
  V_loop    = (1 - f) .* V_mem    + f .* V_loop;  
  
  % calcul des profils 
  warning off
  % missing flux due to particles (3/2 or 5/2 ?)
  % with fixed sources including Qei
  ted1_loop    = - (Qe - Qei - 5/2 .* phys.e .* te_loop .* ge_sts) ./ (chie_loop .* ne_loop .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
  tid1_loop    = - (Qi + Qei - 5/2 .* phys.e .* ti_loop .* gi_sts) ./ (chii_loop .* ne_loop .* ae .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
  % flux equation not output (convention Qualikiz !)
  ge_loop      = ge_sts - V_loop .* ne_loop;
  ned1_loop    = - ge_loop ./ D_loop .* rhomax;
  warning on

  % maximum gradient (to prevent overshoot)
  overgradient = (abs(ted1_loop) > tdmax) | (abs(tid1_loop) > tdmax) | (abs(ned1_loop) > ndmax);
  if any(overgradient)
    overgradient = sum(overgradient(:)) ./ numel(overgradient(:))
  end  
  ted1_loop = max(-tdmax,min(tdmax,ted1_loop));
  ned1_loop = max(-ndmax,min(ndmax,ned1_loop));
  tid1_loop = max(-tdmax,min(tdmax,tid1_loop));
  
  % security
  ted1_loop(~isfinite(ted1_loop)) =ted1(~isfinite(ted1_loop));
  tid1_loop(~isfinite(tid1_loop)) =tid1(~isfinite(tid1_loop));
  ned1_loop(~isfinite(ned1_loop)) =ned1(~isfinite(ned1_loop));
  
  % masquage  (only evolved part of the profile is changed)
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
  
  % valeur minimum (boundary condition)
  te_loop = max(te_loop,profil.tep(:,end) * ones(size(profil.xli)));
  ti_loop = max(ti_loop,profil.tip(:,end) * ones(size(profil.xli)));
  ne_loop = max(ne_loop,profil.nep(:,end) * ones(size(profil.xli)));
  
  %figure(21);plot(profil.xli,ne_loop);drawnow
  
  % amorti
  te_loop = (1 - f) .* te_mem + f .* te_loop;
  ti_loop = (1 - f) .* ti_mem + f .* ti_loop;
  ne_loop = (1 - f) .* ne_mem + f .* ne_loop;
  
  % to turn off density transport
  %  ne_loop = ne_mem;
  % to turn off heqt transport
  % te_loop = te_mem;
  % ti_loop = ti_mem;
  
  % new Qei value
  nHp_loop  = nHp  .*  ne_loop ./ max(1,profil.nep);
  nDp_loop  = nDp  .*  ne_loop ./ max(1,profil.nep);
  nTp_loop  = nTp  .*  ne_loop ./ max(1,profil.nep);
  nhep_loop = nhep .*  ne_loop ./ max(1,profil.nep);
  nzp_loop  = nzp  .*  ne_loop ./ max(1,profil.nep);
  nwp_loop  = nwp  .*  ne_loop ./ max(1,profil.nep);
  % terme source
  warning off
  lnei          =  14.9 - 0.5.*log(ne_loop ./ 1e20) + log(te_loop ./ 1e3);
  warning on
  ind = find(~isfinite(lnei) | (lnei <10));
  if ~isempty(ind)
	  lnei(ind) = 10 .* ones(1,length(ind));
  end
  taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
		  ((phys.e .* max(min(profil.tep(:)),te_loop)) .^ (3/2) ./ ne_loop ./ lnei);
  qeib0    = 3 .* phys.me ./ phys.mp ./ taues .* (nHp_loop + nDp_loop ./ 2 + nTp_loop ./ 3 + nhep_loop  +  ...
	    (option.zimp./ (7/3)) .* nzp_loop + option.rimp .*  (option.zmax./ (7/3)) .* nzp_loop + z0wavez(te_loop) ./ 183.84 .* nwp_loop);
  qeib     = qeib0 .* phys.e .* (te_loop - ti_loop);
  Qei = cumtrapz(profil.xli,qeib .* profil.vpr,2);

  % calcul de la stabilite
  err_stab = std(chie_loop(:) - chie_mem(:)) + std(chii_loop(:) - chii_mem(:)) + std(D_loop(:) - D_mem(:)) + std(V_loop(:) - V_mem(:)) + ...
             std(te_loop(:) - te_mem(:)) ./ 1e3 + std(ti_loop(:) - ti_mem(:)) ./ 1e3 + std(ne_loop(:) - ne_mem(:)) ./ 1e18;
  % condition de sortie
  if err_stab < max(1e-3,option.tol0d)
      break
  elseif recuit
      fprintf('err = %g @ Temp = %g (%g)\n',err_stab,f,nbcor ./ nbi);
  else 
      fprintf('err = %g (%g)\n',err_stab,nbcor ./ nbi); 
  end
  % recopy dans mem
  chie_mem = chie_loop;
  chii_mem = chii_loop;
  D_mem    = D_loop;
  V_mem    = V_loop;
  te_mem   = te_loop;
  ti_mem   = ti_loop;
  ne_mem   = ne_loop;
  
end
% fin de la boucle
fprintf('converged in %d loops\n',k);
% remove wait bar
if first_call == 1
   delete(hwaitbar);
end
mask_err = mask;
mask(mask == 0) = NaN;
mask(mask > 0.5) = 0;

% energy content
wth_nn = (3/2) .* trapz(profil.xli,phys.e .* ne_loop .* (te_loop + ti_loop .* profil.nip ./ max(1,profil.nep)) .* profil.vpr,2);

% output flux
flux_qe = - chie_loop .* ne_loop .* ted1_loop .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
ni_loop = ne_loop .* ae;
flux_qi = - chii_loop .* ni_loop .* tid1_loop .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
flux_qineo = - chii_neo .* ni_loop .* tid1_loop .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2; 
% flux equation not output (convention Qualikiz !)
flux_ne        = - D_loop .* ned1_loop ./ rhomax + V_loop .* ne_loop;



% export data for external data mecanism of METIS
if first_call == 1
  QLK_nn_data.time = profil.temps;
  QLK_nn_data.x    = profil.xli;
  QLK_nn_data.TE   = te_loop;
  QLK_nn_data.TI   = ti_loop;
  QLK_nn_data.NE   = ne_loop;
else
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


% this plot can be improved
if first_call == 1
    hz =findobj(0,'type','figure','tag','Wth_net_kin');
    if isempty(hz)
	      hz=figure('tag','Wth_net_kin','name','Energy content prediction qlkAnn_kin_e (ITG part only, fixed sources but equipartition, no time dependance)');
    else
	      figure(hz);
    end
    clf
    set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	    'defaultlinelinewidth',3,'color',[1 1 1])

    subplot(3,1,1)
    plot(profil.temps,wth ./ 1e6,'b',profil.temps,wth_nn ./ 1e6,'r');
    %xlabel('time (s)');
    ylabel('W_{th} (MJ)');
    legend('METIS','QLKNN-4Dkin');
    zoom on
    subplot(3,1,2);
    semilogy(profil.temps,wth_nn ./ max(1,wth),'b',profil.temps,ones(size(profil.temps)),'g');
    %xlabel('time (s)');
    legend('W_{th,QLKNN-4Dkin} / W_{th,METIS}');
    z0loglin(gca);
    zoom on
    subplot(3,1,3);
    semilogy(profil.temps,factor_time(:,1),'r',profil.temps,ones(size(profil.temps)),'g');
    xlabel('time (s)');
    legend('time factor');
    z0loglin(gca);
    zoom on
    joint_axes(hz,3);
end

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
zplotprof(gca,profil.temps,profil.xli,profil.qe  + profil.qei - 5/2 .* phys.e .* te_loop .* ge_sts,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_qe + Qei,'color','r');
set(gca,'YScale','linear');
ylabel('flux_qe (W/m^{-3})');
%      legend('METIS','QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
title('Without Qei flux and particles convected flux');
zoom yon

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.qi   - profil.qei - 5/2 .* phys.e .* ti_loop .* gi_sts,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_qi - Qei,'color','r');
set(gca,'YScale','linear');
ylabel('flux_qi (W/m^{-3})');
%legend('METIS','QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,ge_sts,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_ne,'color','r');
set(gca,'YScale','linear');
ylabel('flux_ne (e^-/m^{-3}/s)');
legend('METIS (input)','QLKNN-4Dkin (ouput/balance)');
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
zplotprof(gca,profil.temps,profil.xli,profil.xie .* profil.grho2 ./ profil.grho,'color','b');
zplotprof(gca,profil.temps,profil.xli,chie_loop,'color','r');
set(gca,'YScale','linear');
ylabel('Chi_e (m^2/s)');
%legend('METIS','QLKNN-4Dkin prediction');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.xii .* profil.grho2 ./ profil.grho,'color','b');
zplotprof(gca,profil.temps,profil.xli,chii_loop,'color','r');
zplotprof(gca,profil.temps,profil.xli,chii_neo,'color','g');
set(gca,'YScale','linear');
ylabel('Chi_i (m^2/s)');
%legend('METIS','QLKNN-4Dkin prediction','Neo');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,profil.dn,'color','b');
zplotprof(gca,profil.temps,profil.xli,D_loop,'color','r');
set(gca,'YScale','linear');
ylabel('D (m^2/s)');
%legend('METIS','QLKNN-4Dkin prediction');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,4)
zplotprof(gca,profil.temps,profil.xli,-profil.vn ./ profil.dn .* profil.grho2 ./ profil.grho,'color','b');
zplotprof(gca,profil.temps,profil.xli,V_loop ./ D_loop,'color','r');
zplotprof(gca,profil.temps,profil.xli,-VWare_neo ./ profil.dn .* profil.grho2 ./ profil.grho,'color','g');
set(gca,'YScale','linear');
ylabel('V / D (m)');
%legend('METIS','QLKNN-4Dkin prediction');
legend('METIS','QLKNN-4Dkin prediction','Neo');
xlabel('r/a');
%z0loglin(gca);
zoom yon


% backup internal data for future use
tempf = tempname;
clear hz
save(tempf)
zassignin('base','internal_data_QLK',load(tempf)); 
delete(sprintf('%s.mat',tempf));


if ishandle(hdisc)
  delete(hdisc);
end

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
        nI_ZI2 = nI_ZI2 + profil.n1p .* (iso * ones(1,size(profil.n1p,2)));
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
    
% for benchmark
if isfield(profil,'chii_neo')
  chii_neo = profil.chii_neo;
end

%  figure(21);
%  rmx   = profil.rmx ./ (max(profil.rmx,[],2) * ones(size(profil.xli))); 
%  plot(rmx',chii_neo','b');
%  hold on
%  evalin('base','plot(param.gene.x,data.neo.coef.ii ./ data.prof.ni,''r'');');
%  drawnow

