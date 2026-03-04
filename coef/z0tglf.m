% This function is the main function for the use of TGLF in METIS
% It computes as a postprocessing of METIS, using METIS sources and equilirbium Te, Ti and Ne.
% It solve time independant tranport equation (that can be written as an integral)
% Time variation is integrated true energy content via modification of effective sources.
% The convergence is based on simulated annealing algorithm.
% TE, TI and ne computed by this function can later be used as external data in METIS
%
% Warning outputs = TGLF definitions !
%
% test :
% info = declare_TGLF_op_parameters;
% tglfparams = info.valeur;
%[chie,chii,chii_neo,D,V] = z0qlkz_std(tglfparams,post.z0dinput.option,post.zerod,post.profil0d);
%
function [chie_loop,chii_loop,chii_neo,D_loop,V_loop] = z0tglf(tglfparams,option,zerod,profil,loop_in,klz_idx)

disp(loop_in)

% for disclamer display
persistent disclamer_on
persistent data_dir
persistent data 
persistent data_file

if isempty(data)
    data = struct();
end

disp(pwd)
ext_name = sprintf('external%d', klz_idx);


% init output
chie_loop =[];
chii_loop =[];
chii_neo  =[];
D_loop    =[];
V_loop    =[];

% first_call
if isempty(disclamer_on)
    disclamer_on = true;
else
    disclamer_on = false;
end

% for later 
if (nargin > 4) && (loop_in > 1) % loop_in = 2 in z0workingpoint_tglf.m
   first_call = 0;
else
   first_call = 1;
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

%rechantillonage en temps   resampling in time
meff  = interp1(zerod.temps,zerod.meff,profil.temps,'linear','extrap');
modeh = interp1(zerod.temps,zerod.modeh,profil.temps,'linear','extrap');
indice_inv  = interp1(zerod.temps,zerod.indice_inv,profil.temps,'linear','extrap');
wth  = interp1(zerod.temps,zerod.wth,profil.temps,'linear','extrap');		% thermal plasma energy
w    = interp1(zerod.temps,zerod.w,profil.temps,'linear','extrap');		% total plasma energy (joules)
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
n1m  = interp1(zerod.temps,zerod.n1m,profil.temps,'linear','extrap'); % volume­ averaged density of H + D + T ions (m-­3 )
nDm  = interp1(zerod.temps,zerod.nDm,profil.temps,'linear','extrap'); % volume­ averaged density of D ions (m^-­3 )
nTm  = interp1(zerod.temps,zerod.nTm,profil.temps,'linear','extrap'); % volume­ averaged density of T ions (m^-­3 )
nhep = profil.nhep;						      % density of helium (m^-­3 )
nHp  = profil.n1p .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);		% n1p: density of HDT ions (m -3 )
nDp  = profil.n1p .* ((nDm ./ n1m) * ve);
nTp  = profil.n1p .* ((nTm ./ n1m) * ve);
nzp  = profil.nzp;							% density of main impurity (m -3 )
if isfield(profil,'nwp')
    nwp  = profil.nwp;
else
    nwp  = 0 .* nzp;
end

% display chosen parameters
disclamer = 'METIS Te, Ti and ne computation using TGLF (https://gafusion.github.io/doc/tglf.html).';
if first_call == 1
    disp('=============================================================')
    fprintf(disclamer);
    disp('=============================================================')
    disp('TGLF standalone is called with following parameters:')
    fms = get(0,'FormatSpacing');
    format compact
    noms = fieldnames(tglfparams);
    for k= 1:length(noms)
        fprintf('%s = ',noms{k});disp(tglfparams.(noms{k}));
    end
    disp('=============================================================')
    format(fms);
    disp('');
    if disclamer_on
        hdisc = helpdlg(sprintf(disclamer),'Disclamer');
    else
        hdisc = [];
    end
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
chii_neo = tglfparams.factor_neo_chii .* neo_Hinton(profil,option,Amain * ones(size(profil.xli)));

% for Temperature recomputation
chie_ref = profil.xie .* profil.grho2 ./ profil.grho;
chii_ref = profil.xii .* profil.grho2 ./ profil.grho;

% from BgBs model
Dbgbs    = max(tglfparams.coef_min,profil.xie) .* max(tglfparams.coef_min,profil.xii) ./ max(2 .* tglfparams.coef_min,profil.xie + profil.xii) ./ profil.qjli;
Dbgbs(modeh ~= 0,end) = tglfparams.coef_min;
Dbgbs(modeh ~= 0,end - 1) = tglfparams.coef_min;
Dbgbs(modeh ~= 0,end - 2) = 0.5 .* Dbgbs(modeh ~= 0,end - 2) + 0.5 .* tglfparams.coef_min;

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


%filter = find((chie_ref(:) <=0 ) | (chii_ref(:) <=0 ) | (D_ref(:) <=0) |  ...
%	~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
filter = find(~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
chii_ref(filter)=0;
chie_ref(filter)=0;
D_ref(filter)=0;
V_ref(filter)=0;

chie_ref(chie_ref <  tglfparams.minChie) = tglfparams.minChie;
chii_ref(chii_ref <  tglfparams.minChii) = tglfparams.minChii;
V_ref(V_ref < tglfparams.minVe)          = tglfparams.minVe;
D_ref(D_ref < tglfparams.minDe)          = tglfparams.minDe;
chie_ref(chie_ref > tglfparams.maxChie)  = tglfparams.maxChie;
chii_ref(chii_ref > tglfparams.maxChii)  = tglfparams.maxChii;
V_ref(V_ref > tglfparams.maxVe)          = tglfparams.maxVe;
D_ref(D_ref > tglfparams.maxDe)          = tglfparams.maxDe;

RLFS = profil.Raxe .* (1 + profil.epsi);
geo.a  = (profil.Raxe(:,end) .* profil.epsi(:,end));% * ones(size(profil.xli));
%geo.r0 = profil.Raxe(:,end) * ones(size(profil.xli));
geo.R = profil.Raxe(:,end);% * ones(size(profil.xli));
geo.b0 = (profil.fdia(:,end) ./ profil.Raxe(:,end));% * ones(size(profil.xli));

% initialisation
chie_loop = chie_ref;
chii_loop = chii_ref;
D_loop    = D_ref;
V_loop    = V_ref;
disp(['chie_loop_initially=',mat2str(size(chie_loop))]);
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
f =tglfparams.anti_aliasing; %0.1;
% mask de fusion
mask  = ones(size(profil.tep));
%mask(:,end - 4)  = 0.66;
%mask(:,end - 3) = 0.33;
%mask(:,end - 2)  = 0;
% PM: Adapt radial domain with mask=0
[tmp indmask] = min(abs(profil.xli(:)-tglfparams.domain_rho));
simbound = length(profil.xli)-indmask;
%mask(:,end - 1)  = 0;
mask(:,end-simbound)  = 0;

disp(['========= BOUNDARY SET AT ' num2str(profil.xli(end-simbound)) '==========']) 

% PM: Modified edge boundary
indedge = length(profil.xli) - simbound;
%indedge = length(profil.xli) - 1;
%indedge = length(profil.xli) - 2;
ind02 = 5;
switch tglfparams.prediction_domain
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
disp('mask=');
disp(mask);
% boucle de calcul
qr    = profil.qjli;
switch tglfparams.prediction_domain
    case 'Clamped'
        qr(qr < 1)   = 1;
        disp('Call in clamped mode');
end

rhomax = profil.rmx(:,end) * ones(size(profil.xli));
Amain = meff * ones(size(profil.xli));

% Local gradient (for maximum allows value)
ted1  = pdederive(profil.xli,profil.tep,0,2,2,1);
tid1  = pdederive(profil.xli,profil.tip,0,2,2,1);
ned1  = pdederive(profil.xli,profil.nep,0,2,2,1);
tdmax = max(max(abs(ted1(:))),max(abs(tid1(:))));
ndmax = max(abs(ned1(:)));


nbloop_max = tglfparams.nbloop_max; %501
recuit = tglfparams.temperature_factor; %1e-6;
err_mem = Inf;
err_stab = 1/eps;
err_evolution = zeros(1, nbloop_max);  % added

path_to_METIS_coef = fullfile(fileparts(which('install_QLKNN10D')));
path_to_METIS_TGLF = [path_to_METIS_coef '/TGLF_TOOLS'];
addpath(path_to_METIS_TGLF);

if isempty(getappdata(0,'local'))
  output_path = fullfile(path_to_METIS_TGLF,'runs_tglf_data');
else
  output_path = fullfile(getappdata(0,'local'),'runs_tglf_data');  
end  

unix(sprintf(['mkdir -p ' output_path]));

if isempty(data_dir)
    ts = datestr(now,'yyyymmdd_HHMMSS');  % e.g. 20251118_231530123   HHMMSSFFF')
    run_folder_name = ['run_' ts]
    
    data_dir = fullfile(output_path, run_folder_name);
    unix(sprintf(['mkdir -p ' data_dir]));    

    fprintf('====> Plot directory %s\n', data_dir);
    
    % full path to the mat file for this run
    data_file = fullfile(data_dir, ['TGLF_coupling_data_' ts '.mat']);
    
    % create the file immediately with empty data struct
    save(data_file, 'data', '-v7.3');
    fprintf('====== Data will be saved to %s if save flag activated ======\n', data_file);
end

path_plots_internal = fullfile(data_dir, 'Figures_Internal_Convergence');
unix(sprintf(['mkdir -p ' path_plots_internal]));    

path_plots_external = fullfile(data_dir, 'Figures_External_Convergence');
unix(sprintf(['mkdir -p ' path_plots_external]));    

internal_convergence_time = 0;
% Internal loop for convergence 
for k=1:nbloop_max

    fprintf('>>>>>>>>>>>>> TGLF internal iteration %d / %d\n >>>>>>>>>', k, nbloop_max);
    tstart = tic;
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

    % variales pour l'appel des reseaux (avec recuit)
    % % variables for network calls (with annealing)
    profil_tmp=profil;
    if tglfparams.calc_heat_transport == 0
        profil_tmp.tep   = te_loop;
        profil_tmp.tip   = ti_loop;
    else
        profil_tmp.tip   = ti_loop .* (1 + recuit .* f .* rand(size(RLFS)));
        profil_tmp.tep   = te_loop .* (1 + recuit .* f .* rand(size(RLFS)));
    end

    if tglfparams.calc_part_transport == 0
        profil_tmp.nep   = ne_loop;        
    else
        profil_tmp.nep   = ne_loop .* (1 + recuit .* f .* rand(size(RLFS)));        
    end

    % Computes turbulent transport coefficients
    clear post_tmp;
    post_tmp.profil0d=profil_tmp;
    post_tmp.z0dinput.geo = geo;
    post_tmp.z0dinput.option = option;
    post_tmp.zerod = zerod;

   switch tglfparams.SWITCH_tglf
   case 'fast'
     disp('faster one is running......')
     chie_tmp=zeros(1,size(ne_loop,2));
     chii_tmp=zeros(1,size(ne_loop,2));
     D_tmp=zeros(1,size(ne_loop,2));
     V_tmp=zeros(1,size(ne_loop,2));     
   	for kk=3
       		[chie_tmp,chii_tmp,D_tmp,V_tmp] = run_tglf(tglfparams,post_tmp,kk);
    	end
    	chie_loop = repmat(chie_tmp,size(chie_loop,1),1);
   	chii_loop = repmat(chii_tmp,size(chie_loop,1),1);
  	D_loop = repmat(D_tmp,size(chie_loop,1),1);
   	V_loop = repmat(V_tmp,size(chie_loop,1),1);
    
    case 'long'
     	% This option is longer    
     	disp('Attention: This will take long time to finish!')
    	for kk=1:size(chie_loop,1)
      	[chie_loop(kk,:),chii_loop(kk,:),D_loop(kk,:),V_loop(kk,:)] = run_tglf(tglfparams,post_tmp,kk);
    	end
    end
   
    TGLF_data.Chi_e_TGLF = chie_loop; 
    TGLF_data.Chi_i_TGLF = chii_loop;
    TGLF_data.D_ne_TGLF  = D_loop;
    TGLF_data.V_ne_TGLF  = V_loop;   

    %
    if tglfparams.verbosity > 5
        figure(300);
        subplot(2,2,1)
        plot(profil.xli,chie_loop);
        subplot(2,2,2)
        plot(profil.xli,chii_loop);
        subplot(2,2,3)
        plot(profil.xli,D_loop);
        subplot(2,2,4)
        plot(profil.xli,V_loop);
        drawnow
    end

    % filter
    filter = ~isfinite(chie_loop(:)) | ~isfinite(chii_loop(:)) | ~isfinite(D_loop(:)) | ~isfinite(V_loop(:));
    
    if any(filter)
        disp('bad values in TGLF results');
    end

    chii_loop(filter) = chii_mem(filter);
    chie_loop(filter) = chie_mem(filter);
    D_loop(filter)    = D_mem(filter);
    V_loop(filter)    = V_mem(filter);

    
    test = (chie_loop > tglfparams.maxChie) | (chii_loop > tglfparams.maxChii) |  ...
        (D_loop > tglfparams.maxDe) | (V_loop > tglfparams.maxVe) | (V_loop < tglfparams.minVe);
    overshoot = sum(double(test(:)))./ numel(test);
    
    chie_loop(chie_loop > tglfparams.maxChie) = chie_mem(chie_loop > tglfparams.maxChie);
    chii_loop(chii_loop > tglfparams.maxChii) = chii_mem(chii_loop > tglfparams.maxChii);
    D_loop(D_loop > tglfparams.maxDe)         = D_mem(D_loop > tglfparams.maxDe);
    V_loop(V_loop > tglfparams.maxVe)         = V_mem(V_loop > tglfparams.maxVe);
    chie_loop(chie_loop < tglfparams.minChie) = tglfparams.minChie;
    chii_loop(chii_loop < tglfparams.minChii) = tglfparams.minChii;
    V_loop(V_loop < tglfparams.minVe)         = tglfparams.minVe;
    D_loop(D_loop < tglfparams.minDe)         = tglfparams.minDe;
    
    % Neoclassical transport
    profil_loc = profil;
    profil_loc.tep = te_loop;
    profil_loc.tip = ti_loop;
    profil_loc.nep = ne_loop;
    % change on other species
    profil_loc.nHp  = nHp  .*  ne_loop ./ max(1,profil.nep);
    profil_loc.nDp  = nDp  .*  ne_loop ./ max(1,profil.nep);
    profil_loc.nTp  = nTp  .*  ne_loop ./ max(1,profil.nep);
    profil_loc.nhep = nhep .*  ne_loop ./ max(1,profil.nep);
    profil_loc.nzp  = nzp  .*  ne_loop ./ max(1,profil.nep);
    profil_loc.nwp  = nwp  .*  ne_loop ./ max(1,profil.nep);
    profil_loc.n1p  = profil_loc.nHp + profil_loc.nDp + profil_loc.nTp;
    profil_loc.nip  = profil_loc.n1p + profil_loc.nhep + profil_loc.nzp + profil_loc.nwp;
    chii_neo  = tglfparams.factor_neo_chii .* neo_Hinton(profil_loc,option,Amain);
    chii_neo  = chii_neo .* profil.grho2 ./ profil.grho;
    VWare_neo = Vware_ref; % not updated
    
    % sum
    chii_loop = chii_loop + chii_neo;
    V_loop    = V_loop + VWare_neo;
    D_loop    = max(D_loop,tglfparams.neo_ion_in_D .* chii_neo);
    chie_loop = max(chie_loop,tglfparams.neo_ion_in_Chie .* chii_neo);
    
    % space filter
    if tglfparams.smooth_space
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
    
    
    switch tglfparams.SWITCH_cushion_diff_coef
    case 'yes'
 	chie_loop = (1 - f) .* chie_mem + f .* chie_loop;
    	chii_loop = (1 - f) .* chii_mem + f .* chii_loop;
   	D_loop    = (1 - f) .* D_mem    + f .* D_loop;
    	V_loop    = (1 - f) .* V_mem    + f .* V_loop;
    
    case 'no'
    	disp('no cushoning happened.')
    end

    %
    TGLF_data.Chi_e = chie_loop;
    TGLF_data.Chi_i = chii_loop;
    TGLF_data.D_ne  = D_loop;
    TGLF_data.V_ne  = V_loop;
    %
    
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
    overgradient = double(mask)  .* double((abs(ted1_loop) > (tdmax/sqrt(f))) | (abs(tid1_loop) > (tdmax/sqrt(f))) | (abs(ned1_loop) > (ndmax/sqrt(f))));
    overgradient = sum(overgradient(:)) ./ numel(overgradient(:));
    
    ted1_loop = max(-tdmax/sqrt(f),min(tdmax/sqrt(f),ted1_loop));
    ned1_loop = max(-ndmax/sqrt(f),min(ndmax/sqrt(f),ned1_loop));
    tid1_loop = max(-tdmax/sqrt(f),min(tdmax/sqrt(f),tid1_loop));
    
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
    
    % synchronistion with option in QLK-NN
    if tglfparams.calc_heat_transport == 0
	te_loop = te_mem;
	ti_loop = ti_mem;
    end
    if tglfparams.calc_part_transport == 0
	ne_loop = ne_mem;
    end
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
    err_evolution(k) = err_stab; % added 
    
   
   % save data and plots for every 10 iterations.
    if mod(k,5) == 0 || k == 1 || k == nbloop_max || err_stab < max(1e-3,option.tol0d)

        int_name = sprintf('internal%d', k);   % e.g. 'internal1', 'internal50'

        data.xli = profil.xli;
        data.(ext_name).(int_name).te_loop_tglf = te_loop;
        data.(ext_name).(int_name).ti_loop_tglf = ti_loop;
        data.(ext_name).(int_name).ne_loop_tglf = ne_loop;
        
        data.(ext_name).(int_name).te_metis = profil.tep;
        data.(ext_name).(int_name).ti_metis = profil.tip;
        data.(ext_name).(int_name).ne_metis = profil.nep;

        data.(ext_name).(int_name).chie_tglf_raw = TGLF_data.Chi_e_TGLF;
        data.(ext_name).(int_name).chii_tglf_raw = TGLF_data.Chi_i_TGLF;
        data.(ext_name).(int_name).D_tglf_raw = TGLF_data.D_ne_TGLF;
        data.(ext_name).(int_name).V_tglf_raw = TGLF_data.V_ne_TGLF;
     
        data.(ext_name).(int_name).chie_loop_tglf  = chie_loop;
        data.(ext_name).(int_name).chii_loop_tglf  = chii_loop;
        data.(ext_name).(int_name).D_loop_tglf     = D_loop;
        data.(ext_name).(int_name).V_loop_tglf     = V_loop;
          
        data.(ext_name).(int_name).chie_metis = profil.xie .* profil.grho2 ./ profil.grho;
        data.(ext_name).(int_name).chii_metis = profil.xii .* profil.grho2 ./ profil.grho;
        data.(ext_name).(int_name).D_metis = profil.dn;
        data.(ext_name).(int_name).V_metis = profil.vn;
        
        % convergence error for this k
        data.(ext_name).(int_name).error = err_stab;             

	if tglfparams.save_intermediate_data
          save(data_file, 'data', '-v7.3');
        end

        % Transport coefficients plots %%%%%%%%%%
        hz_transp =findobj(0,'type','figure','tag','Transp. Coef.');
        if isempty(hz_transp)
          hz_transp=figure('tag','Transp. Coef.','name',sprintf('Transport coefficients comparison between METIS & TGLF'));
        else
          figure(hz_transp);
        end
        clf
        set(hz_transp, 'Position', [100 100 1200 900]);
	
    
        % --- Chi_e    
        subplot(2,2,1)
        plot(profil.xli, profil.xie(1,:) .* profil.grho2(1,:) ./ profil.grho(1,:),'b', 'LineWidth', 2); hold on; % metis
        plot(profil.xli, chie_loop(1,:), 'r', 'LineWidth', 2);hold on;						 % trimmed from tglf
        %plot(profil.xli,TGLF_data.Chi_e(1,:), 'g', 'LineWidth', 2);hold on;
        plot(profil.xli, TGLF_data.Chi_e_TGLF(1,:), 'm', 'LineWidth', 2,'linestyle','-');					% tglf raw 
        xlabel('r/a'); ylabel('chi_e (m^2/s)', 'Interpreter', 'tex');
        %title(sprintf('\\chi_e @ k = %d', k));
        %legend('METIS', 'TGLF', 'Raw TGLF');
        grid on;
        zoom yon 


        % --- Chi_i
        subplot(2,2,2)
        plot(profil.xli, profil.xii(1,:) .* profil.grho2(1,:) ./profil.grho(1,:),'b', 'LineWidth', 2); hold on;
        plot(profil.xli, chii_loop(1,:), 'r', 'LineWidth', 2);hold on;
        plot(profil.xli, TGLF_data.Chi_i_TGLF(1,:), 'm', 'LineWidth', 2,'linestyle','-');hold on;
        plot(profil.xli, chii_neo(1,:), 'g', 'LineWidth', 2);
        xlabel('r/a'); ylabel('chi_i (m^2/s)');
        %title(sprintf('\\chi_e and \\chi_i @ k = %d', k));
        %legend('METIS', 'TGLF','Raw TGLF','Neo'), 'Location', 'Best');
        grid on;
        zoom yon 

        % --- D ---
        subplot(2,2,3)
        plot(profil.xli, profil.dn(1,:), 'b', 'LineWidth', 2,'Displayname','METIS'); hold on;
        plot(profil.xli, D_loop(1,:), 'r', 'LineWidth',2,'Displayname','TGLF');hold on;
        plot(profil.xli, TGLF_data.D_ne_TGLF(1,:), 'm', 'LineWidth',2,'Displayname','Raw TGLF','linestyle','-');
        xlabel('r/a'); ylabel('D (m^2/s)');
        %title(sprintf(' D @ k = %d', k));
        legend('location','best');
        grid on;
        zoom yon 

        % --- V ---
        subplot(2,2,4)
        plot(profil.xli, -profil.vn(1,:) .* profil.grho2(1,:) ./ profil.grho(1,:), 'b', 'LineWidth', 2,'Displayname','METIS'); hold on;
        plot(profil.xli, V_loop(1,:) , 'r', 'LineWidth', 2,'Displayname','TGLF'); hold on;
        plot(profil.xli, TGLF_data.V_ne_TGLF(1,:).* profil.grho2(1,:) ./ profil.grho(1,:), 'm', 'LineWidth', 2,'Displayname','Raw TGLF','linestyle','-'); hold on;
        plot(profil.xli, -VWare_neo(1,:) .* profil.grho2(1,:) ./ profil.grho(1,:), 'g', 'LineWidth', 2,'Displayname','Neo');
        xlabel('r/a'); ylabel('V (m/s)');
        %title(sprintf(' V @ k = %d', k));
        legend('location','best');
        grid on;  

	if tglfparams.save_all_figs
          saveas(hz_transp,fullfile(path_plots_internal,sprintf('trasnport_coeff_klz%02d_it%03d.png',klz_idx,k)));
	end        

        % ---------- Save kinetic profiles ----------
        hz_prof =findobj(0,'type','figure','tag','Profiles');
        if isempty(hz_prof)
          hz_prof=figure('tag','Profiles','name',sprintf('Profiles prediction with TGLF'));
        else
          figure(hz_prof);
        end
        clf
        set(hz_prof, 'Position', [100 100 1200 900]);
        
        % --- T_e profile ---
        subplot(2,2,1)
        plot(profil.xli, profil.tep(1,:) / 1e3, 'b', 'LineWidth',2); hold on;
        plot(profil.xli, te_loop(1,:) / 1e3, 'r', 'LineWidth', 2);
        xlabel('r/a'); ylabel('T_e (keV)');
        grid on;
        %title(sprintf('Electron Temperature @ k = %d', k));
        %legend('METIS', 'TGLF'); grid on;

        % --- T_i profile ---
        subplot(2,2,2)
        plot(profil.xli, profil.tip(1,:) / 1e3, 'b', 'LineWidth', 2); hold on;
        plot(profil.xli, ti_loop(1,:) / 1e3, 'r', 'LineWidth', 2);
        xlabel('r/a'); ylabel('T_i (keV)');
        grid on;
        %title(sprintf('Ion Temperature @ k = %d', k));
        %legend('METIS', 'TGLF'); grid on;

        % --- n_e profile ---
        subplot(2,2,3)
        plot(profil.xli, profil.nep(1,:) / 1e19, 'b', 'LineWidth', 2); hold on;
        plot(profil.xli, ne_loop(1,:) / 1e19, 'r', 'LineWidth', 2);
        xlabel('r/a'); ylabel('n_e (10^{19} m^{-3})');
        %title(sprintf('Electron Density @ k = %d', k));
        legend('METIS', 'TGLF'); 
        grid on; 

	if tglfparams.save_all_figs
          saveas(hz_prof,fullfile(path_plots_internal,sprintf('profiles_klz%02d_it%03d.png',klz_idx,k)));
        end
	
	% --- Error evolution so far ---
        hz_err =findobj(0,'type','figure','tag','Error');
        if isempty(hz_err)
          hz_err=figure('tag','Error','name',sprintf('Error evolution from internal transport loop.'));
          set(hz_err, 'Position', [100 100 1200 900]);	  
        else
          figure(hz_err);
        end
        hold on

        semilogy(1:k, err_evolution(1:k), 'ko-', 'LineWidth', 2); 		
        grid on;
        xlabel('Iteration'); ylabel('Error');
        title('Convergence Error');
	
	if tglfparams.save_all_figs	
          saveas(hz_err,fullfile(path_plots_internal,sprintf('error_klz%02d_it%03d.png',klz_idx,k)));
        end


    end

    if err_stab < max(1e-3,option.tol0d)
        break
    elseif recuit
        fprintf('error = %g @ temperature = %g (overshoot = %g %% , overgradient = %g %%)\n',err_stab,f,overshoot*100,overgradient*100);
    else
        fprintf('error = %g (overshoot = %g %% , overgradient = %g %%)\n',err_stab,overshoot*100,overgradient*100);
    end
    % recopy dans mem
    chie_mem = chie_loop;
    chii_mem = chii_loop;
    D_mem    = D_loop;
    V_mem    = V_loop;
    te_mem   = te_loop;
    ti_mem   = ti_loop;
    ne_mem   = ne_loop;
    
t_per_it = toc(tstart);
internal_convergence_time =  internal_convergence_time + t_per_it ;
%fprintf('internal convergence completed in %d iterations and %.2f seconds \n',internal_it_count,internal_convergence_time);

end
% fin de la boucle end of loop

internal_it_count = k;
%t_per_it = toc(tstart);
%internal_convergence_time =  internal_convergence_time + t_per_it ;
fprintf('internal convergence completed in %d iterations and %.2f seconds \n',internal_it_count,internal_convergence_time);


%figure; % I added this
%semilogy(1:k, err_evolution(1:k), 'o-','LineWidth', 2);
%xlabel('Iteration');
%ylabel('Convergence Error (err\_stab)');
%title('Evolution of Convergence Error in z0tglf');
%grid on;

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

% export data for external data mechanism of METIS
if first_call == 1
    TGLF_data.time = profil.temps;
    TGLF_data.x    = profil.xli;
    TGLF_data.TE   = te_loop;
    TGLF_data.TI   = ti_loop;
    TGLF_data.NE   = ne_loop; 
    TGLF_data.err  = Inf;
else
    if evalin('base','exist(''TGLF_data'')')
        err_press_mem = evalin('base','TGLF_data.err');
    else
        err_press_mem = Inf;
    end
    TGLF_data.time = profil.temps;
    TGLF_data.x    = profil.xli;
    TGLF_data.err  = sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2));
    fprintf('TGLF_data.err = %g',TGLF_data.err) % this error is for the external loop
    [err_press,ind_best] = min(sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2)));
    TGLF_data.err  = err_press;
    vt = ones(size(profil.temps));
    fprintf('TGLF_data.err = err_press = %g',TGLF_data.err)
    if err_press > err_press_mem
        te_exp = 0.1 .* te_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.9 .* profil.tep(ind_best,:) .* profil.nep(ind_best,:);
        ti_exp = 0.1 .* ti_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.9 .* profil.tip(ind_best,:) .* profil.nep(ind_best,:);
        ne_exp = 0.1 .* ne_loop(ind_best,:) + 0.9 .* profil.nep(ind_best,:);
    else
        te_exp = 0.3 .* te_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.7 .* profil.tep(ind_best,:) .* profil.nep(ind_best,:);
        ti_exp = 0.3 .* ti_loop(ind_best,:) .* ne_loop(ind_best,:) + 0.7 .* profil.tip(ind_best,:) .* profil.nep(ind_best,:);
        ne_exp = 0.3 .* ne_loop(ind_best,:) + 0.7 .* profil.nep(ind_best,:);
    end
    TGLF_data.TE   = vt * (te_exp ./ ne_exp);
    TGLF_data.TI   = vt * (ti_exp ./ ne_exp);
    TGLF_data.NE   = vt * ne_exp;
    TGLF_data.N_bar = vt * trapz(profil.xli,ne_exp);
end
zassignin('base','TGLF_data',TGLF_data);
disp(TGLF_data)


data.(ext_name).final.te_tglf = te_loop;   % what goes back to METIS
data.(ext_name).final.ti_tglf = ti_loop;
data.(ext_name).final.ne_tglf = ne_loop;

data.(ext_name).final.chie_tglf    = chie_loop;
data.(ext_name).final.chii_tglf    = chii_loop;
data.(ext_name).final.D_tglf       = D_loop;
data.(ext_name).final.V_tglf       = V_loop;

if first_call == 0
    data.(ext_name).final.Te_passed  = TGLF_data.TE;
    data.(ext_name).final.Ti_passed  = TGLF_data.TI;
    data.(ext_name).final.ne_passed  = TGLF_data.NE;
end


data.(ext_name).final.te_metis = profil.tep;
data.(ext_name).final.ti_metis = profil.tip;
data.(ext_name).final.ne_metis = profil.nep;

data.(ext_name).final.chie_metis = profil.xie .* profil.grho2 ./ profil.grho;
data.(ext_name).final.chii_metis = profil.xii .* profil.grho2 ./ profil.grho;
data.(ext_name).final.D_metis = profil.dn;
data.(ext_name).final.V_metis = profil.vn;

if tglfparams.save_intermediate_data
  save(data_file, 'data', '-v7.3');
  assignin('base','data',data);
end

%outer_it_suffix = sprintf('outer_it_klz%02d',klz_idx);

%%% net Energy Plots %%%%%%%%%%%%%%%%%%%%%%%%%5 
% this plot can be improved
if first_call == 1
    hz_energy =findobj(0,'type','figure','tag','Wth_net_kin');
    if isempty(hz_energy)
        hz_energy=figure('tag','Wth_net_kin','name',sprintf('Energy content prediction %s (fixed sources but equipartition, no time dependance)',tglfparams.channel_tag));
    else
        figure(hz_energy);
    end
    clf
    set(hz_energy,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
        'defaultlinelinewidth',3,'color',[1 1 1])
    
    subplot(3,1,1)
    plot(profil.temps,wth ./ 1e6,'b',profil.temps,wth_nn ./ 1e6,'r');
    %xlabel('time (s)');
    ylabel('W_{th} (MJ)');
    legend('METIS',tglfparams.channel_tag);
    zoom on
    subplot(3,1,2);
    semilogy(profil.temps,wth_nn ./ max(1,wth),'b',profil.temps,ones(size(profil.temps)),'g');
    %xlabel('time (s)');
    legend(sprintf('W_{th,%s} / W_{th,METIS}',tglfparams.channel_tag));
    z0loglin(gca);
    zoom on
    subplot(3,1,3);
    semilogy(profil.temps,factor_time(:,1),'r',profil.temps,ones(size(profil.temps)),'g');
    xlabel('time (s)');
    legend('time factor');
    z0loglin(gca);
    zoom on
    joint_axes(hz_energy,3);
    %saveas(hz_energy,fullfile(path_plots_external,sprintf('net_energy_klz%02d.png',klz_idx)));
end
%saveas(hz_energy,fullfile(path_plots_external,sprintf('net_energy_klz%02d.png',klz_idx)));
%saveas(hz_energy,fullfile(path_plots_external,sprintf('net_energy_klz%02d.fig',klz_idx)));
%close(hz_energy);

if ishandle(hdisc)
    delete(hdisc);
end

% backup internal data for future use
tempf = tempname;
clear hz_prof hz_err hz_energy hz_flux hz_transp ans hdisc hwaitbar
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
        alpha  = nI_ZI2 ./ profil.nhep ./ 4;
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
part2   = part2 .* (1 + (1.33 .* alpha .* (1 + 0.60 .* alpha)) ./ (1 + 1.79 .* alpha));

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

