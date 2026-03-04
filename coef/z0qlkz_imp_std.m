% This function is the main function for the use of QLK-NN 10D in METIS
% It computes as a postprocessing of METIS, using METIS sources and equilirbium Te, Ti and Ne.
% It solve time independant tranport equation (taht can be write as an integral)
% Time variation is integrated true energy content via modification of effective sources.
% The convergence in based on simulated annealing algorithm.
% TE, TI and ne computed by this function can later be used as external  data in METIS
%
%
%%%%%%%%%%%%%%%%% IMPURITY TRANSPORT VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Warning outputs = QUALIKIZ definitions !
%
% test :
% info = declare_QLKNN10D_parameters;
% qlkparams = info.valeur;
%[chie,chii,chii_neo,D,V] = z0qlkz_std(qlkparams,post.z0dinput.option,post.zerod,post.profil0d);
%
function varargout = z0qlkz_imp_std(qlkparams,option,zerod,profil,loop_in,twp,z0dinput)
clear QLK_nn_data
% for disclamer display
persistent disclamer_on

% init output
chie_loop =[];
chii_loop =[];
chii_neo  =[];
D_loop    =[];
V_loop    =[];
DD_loop   =[];
VV_loop   =[];
G_loop    =[];

DD_loop_n   =[];
VV_loop_n   =[];
G_loop_n    =[];

n_loop    =[];
% n_loop
% first_call
if isempty(disclamer_on)
    disclamer_on = true;
else
    disclamer_on = false;
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

% for later 
if (nargin > 4) && (loop_in > 1)
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

%rechantillonage en tempsefe_GB
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
nep  = profil.nep;
nHp  = profil.n1p .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);
nDp  = profil.n1p .* ((nDm ./ n1m) * ve);
nTp  = profil.n1p .* ((nTm ./ n1m) * ve);
nzp  = profil.nzp;
if isfield(profil,'nwp')
    nwp  = profil.nwp;
else
    nwp  = 0 .* nzp;
end
switch option.gaz
    case 5
        nhe3p = nhep;
        nBp   = zeros(size(nTp));
        nhep  = option.frhe0 .* nep;
    case 11
        nBp = nTp;
        nTp = zeros(size(nTp));
        nhe3p = zeros(size(nTp));
    otherwise
        nBp = zeros(size(nTp));
        nhe3p = zeros(size(nTp));
end

% improve precision
[A_el,Z_el,name_el] = chargemasse;
dd   = abs(Z_el - option.zimp);
mask = (dd == min(dd));
Aimp = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(Aimp)
    Aimp = 7/3 .* option.zimp;
end
dd   = abs(Z_el - option.zmax);
mask = (dd == min(dd));
Amax = sum(A_el .* mask) ./ max(1,sum(mask));
if ~isfinite(Amax)
    Amax = 7/3 .* option.zmax;
end



% display choosen parameters
disclamer = 'METIS Te, Ti and ne computation using QUALIKIZ (https://gitlab.com/qualikiz-group/QuaLiKiz/-/wikis/Physics-of-QuaLiKiz).\n';
if first_call == 1
    disp('=============================================================')
    fprintf(disclamer);
    disp('=============================================================')
    disp('QuaLiKiz standalone is called with following parameters:')
    fms = get(0,'FormatSpacing');
    format compact
    noms = fieldnames(qlkparams);
    for k= 1:length(noms)
        fprintf('%s = ',noms{k});disp(qlkparams.(noms{k}));
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

%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial transport coefficient to start the convergence
% neoclassical part (use se same on ion and electron)
% Amain
Amain = meff;
chii_neo = qlkparams.factor_neo_chii .* neo_Hinton(profil,option,Amain * ones(size(profil.xli)));

% for Temperature recomputation
chie_ref = profil.xie .* profil.grho2 ./ profil.grho;
chii_ref = profil.xii .* profil.grho2 ./ profil.grho;

% from BgBs model
Dbgbs    = max(qlkparams.coef_min,profil.xie) .* max(qlkparams.coef_min,profil.xii) ./ max(2 .* qlkparams.coef_min,profil.xie + profil.xii) ./ profil.qjli;
Dbgbs(modeh ~= 0,end) = qlkparams.coef_min;
Dbgbs(modeh ~= 0,end - 1) = qlkparams.coef_min;
Dbgbs(modeh ~= 0,end - 2) = 0.5 .* Dbgbs(modeh ~= 0,end - 2) + 0.5 .* qlkparams.coef_min;

n_loop = repmat(profil.nep, [1, 1, 6]);  %electrons are on the 6th positions
n_loop(:,:,1) = profil.nip;
n_loop(:,:,2) = profil.nzp;
n_loop(:,:,3) = profil.nzp.*option.rimp;
n_loop(:,:,4) = profil.nwp;
n_loop(:,:,5) = profil.nip.*option.cmin;
n_init = n_loop;

for j=1:6
    nd1(:,:,j) = pdederive(profil.xli,n_init(:,:,j),0,2,2,1);
end


ned1  = pdederive(profil.xli,profil.nep,0,2,2,1);
Vware = profil.ware;   % have the same sign than in CRONOS
ve = ones(size(profil.xli));

for j=1:6   % Not true
    warning off
    DD_ref(:,:,j) = - (Vware .* n_init(:,:,j) + profil.ge ) ./ nd1(:,:,j) .* (profil.rmx(:,end) * ve);
    warning on
    max(Dbgbs,DD_ref(:,:,j));

    VV_ref(:,:,j) = - (profil.ge + DD_ref(:,:,j) .* nd1(:,:,j) ./ (profil.rmx(:,end) * ve)) ./ n_init(:,:,j);
    VV_ref(:,:,j) = VV_ref(:,:,j) - Vware;
    VV_ref(:,:,j) = max(0, VV_ref(:,:,j));

    warning off
    DD_ref(:,:,j) = - ((VV_ref(:,:,j) + Vware) .* n_init(:,:,j) + profil.ge) ./ nd1(:,:,j) .* (profil.rmx(:,end) * ve);
    DD_ref(:,1,j) = DD_ref(:,2,j);
    warning on
end


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
% calcul final de D (initial)
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
prad  =  cumtrapz(profil.xli,profil.prad .* profil.vpr,2);
Qe    = (profil.qe + profil.qei + prad) .* factor_time;  %net electron heat profil
Qi    = (profil.qi - profil.qei) .* factor_time;
Qei   = profil.qei .* factor_time;
Qrad  = prad .* factor_time;
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

% Initialisation tableau d’erreur
err_stab_vec = [];
err_n_vec = [];
err_rad_W_vec = [];
err_rad_z1_vec = [];
err_rad_z2_vec = [];
err_T_vec = [];
err_chi_vec = [];
err_D_vec = [];
err_D_n_vec = [];

% palette, on retient juste la 1ère et la 6e couleur
colors = lines(6);
colIon  = colors(1,:);   % ions  : bleu-vert
colElec = colors(6,:);   % électrons : violet

% Atomic data
%atom=Lz_zave_HDT0_default;
atom_tmp = load('Lz_zave.mat');
atom=atom_tmp.tabmat;

if qlkparams.imp_trans == 3 
    [tmp itime] = min(abs(profil.temps-twp));
    facitparams.ION_REGIME = 1;
    facitparams.OPT_ASYM = 0;
    facitparams.LAB_GEO = 0;
    facitparams.MULTIS = 0;
    facitparams.PLOT = 0;
   
end

if qlkparams.imp_trans == 4
    [tmp itime] = min(abs(profil.temps-twp));
    neoparams.SIM_MODEL = 1;
    neoparams.EQUILIBRIUM_MODEL = 2;
    neoparams.COLLISION_MODEL = 4;
    neoparams.N_ENERGY = 6;
    neoparams.N_XI = 11;
    neoparams.N_THETA = 21;
    neoparams.DV = 2;

    
end

%%%%%%%%%%%%%%% END INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% CONVERGENCE LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%filter = find((chie_ref(:) <=0 ) | (chii_ref(:) <=0 ) | (D_ref(:) <=0) |  ...
%	~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
filter = find(~isfinite(chie_ref(:)) | ~isfinite(chii_ref(:)) | ~isfinite(D_ref(:)) | ~isfinite(V_ref(:)));
filter_bis = find(~isfinite(DD_ref(:)) | ~isfinite(VV_ref(:)));
chii_ref(filter)=0;
chie_ref(filter)=0;
D_ref(filter)=0;
V_ref(filter)=0;

DD_ref(filter_bis)=0;
VV_ref(filter_bis)=0;

chie_ref(chie_ref <  qlkparams.minChie) = qlkparams.minChie;
chii_ref(chii_ref <  qlkparams.minChii) = qlkparams.minChii;
V_ref(V_ref < qlkparams.minVe)          = qlkparams.minVe;
D_ref(D_ref < qlkparams.minDe)          = qlkparams.minDe;
chie_ref(chie_ref > qlkparams.maxChie)  = qlkparams.maxChie;
chii_ref(chii_ref > qlkparams.maxChii)  = qlkparams.maxChii;
V_ref(V_ref > qlkparams.maxVe)          = qlkparams.maxVe;
D_ref(D_ref > qlkparams.maxDe)          = qlkparams.maxDe;

VV_ref(VV_ref < qlkparams.minVe)          = qlkparams.minVe;
DD_ref(DD_ref < qlkparams.minDe)          = qlkparams.minDe;
VV_ref(VV_ref > qlkparams.maxVe)          = qlkparams.maxVe;
DD_ref(DD_ref > qlkparams.maxDe)          = qlkparams.maxDe;
% initialisation
chie_loop = ones(size(profil.tep));    %size(profil.tep) = (insize,21)
chii_loop = 3 .* ones(size(profil.tep));
D_loop    = ones(size(profil.tep)) / 3;
V_loop    = Vware_ref;

DD_loop = repmat(D_loop, [1, 1, 6]);
VV_loop = repmat(Vware_ref, [1, 1, 6]);

%
chie_mem = chie_loop;
chii_mem = chii_loop;
D_mem    = D_loop;
V_mem    = V_loop;

DD_mem = DD_loop;
VV_mem = VV_loop;

DD_loop_n = repmat(D_loop, [1, 1, 6]);
VV_loop_n = repmat(Vware_ref, [1, 1, 6]);

DD_mem_n = ones(size(DD_loop_n))*0.001;
VV_mem_n = ones(size(VV_loop_n))*0.001;


%
te_loop  = profil.tep;
ti_loop  = profil.tip;
ne_loop  = profil.nep;

% n_loop = repmat(ne_loop, [1, 1, 6]);  %electrons are on the 6th positions
% n_loop(:,:,1) = profil.nip;
% n_loop(:,:,2) = profil.nzp;
% n_loop(:,:,3) = profil.nzp.*option.rimp;
% n_loop(:,:,4) = profil.nwp;
% n_loop(:,:,5) = profil.nip.*option.cmin;
% n_init = n_loop;
%
te_mem  = profil.tep;
ti_mem  = profil.tip;
ne_mem  = ne_loop;

n_mem = n_loop;



% anti aliasing
f =qlkparams.anti_aliasing; %0.1;
% mask de fusion
mask  = ones(size(profil.tep));
%mask(:,end - 4)  = 0.66;
%mask(:,end - 3) = 0.33;
%mask(:,end - 2)  = 0;
% PM: Adapt radial domain with mask=0
[tmp indmask] = min(abs(profil.xli(:)-qlkparams.domain_rho));
simbound = length(profil.xli)-indmask;
%mask(:,end - 1)  = 0;
mask(:,end-simbound)  = 0;

disp(['========= BOUNDARY SET AT ' num2str(profil.xli(end-simbound)) '==========']) 

% PM: Modified edge boundary
indedge = length(profil.xli) - simbound;
%indedge = length(profil.xli) - 1;
%indedge = length(profil.xli) - 2;
ind02 = 5;
switch qlkparams.prediction_domain
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
% Loop de calcul
qr    = profil.qjli;
switch qlkparams.prediction_domain
    case 'Clamped'
        qr(qr < 1)   = 1;
        disp('Call in clamped mode');
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
% for j=1:6
%     nd1(:,:,j) = pdederive(profil.xli,n_init(:,:,j),0,2,2,1);DD_loop_n
% end
ned1  = nd1(:,:,end);
tdmax = max(max(abs(ted1(:))),max(abs(tid1(:))));
ndmax = max(abs(nd1(:)));


nbloop_max = qlkparams.nbloop_max; %501
recuit = qlkparams.temperature_factor; %1e-6;
err_mem = Inf;
err_stab = 1/eps;

err_mem_n   = Inf;
err_mem_rad_W = Inf;
err_mem_rad_z1 = Inf;
err_mem_rad_z2 = Inf;
err_mem_T   = Inf;
err_mem_chi   = Inf;
err_mem_D   = Inf;
err_mem_D_n   = Inf;

err_n = 1/eps;
err_rad_W = 1/eps;
err_rad_z1 = 1/eps;
err_rad_z2 = 1/eps;
err_T = 1/eps;
err_chi = 1/eps;
err_D = 1/eps;
err_D_n = 1/eps;

for k = 1:nbloop_max
    if (k > nbloop_max / 3) || ...
    (err_n   > err_mem_n) || ...
    (err_rad_W > err_mem_rad_W) || ...
    (err_rad_z1 > err_mem_rad_z1) || ...
    (err_rad_z2 > err_mem_rad_z2) || ...
    (err_T   > err_mem_T) || ...
    (err_chi   > err_mem_chi) || ...
    (err_D   > err_mem_D) || ...
    (err_D_n   > err_mem_D_n)


        if recuit == 0
            f = 0.3;
        end
        f = f * 0.95;
        recuit = 1;

    else
        err_mem_n   = err_n;
        err_mem_rad_W = err_rad_W;
        err_mem_rad_z1 = err_rad_z1;
        err_mem_rad_z2 = err_rad_z2;
        err_mem_T   = err_T;
        err_mem_chi   = err_chi;
        err_mem_D   = err_D;
        err_mem_D_n   = err_D_n;
    end

    if first_call == 1
        waitbar(k ./ nbloop_max,hwaitbar,'Please wait...');
    end

    % variales pour l'appel des reseaux (avec recuit)
    if qlkparams.calc_heat_transport == 0
        te   = te_loop;
        ti   = ti_loop;
    else
        ti   = ti_loop .* (1 + recuit .* f .* rand(size(RLFS)));
        te   = te_loop .* (1 + recuit .* f .* rand(size(RLFS)));
    end

    z1 = option.zimp;  %for this one we specify a Z and a profile
    z2 = option.zmax;  %for this one we specify a Z and take the profile of z1
    z3 = round((1 - Sn_fraction) .* z0wavez(te*1000) + Sn_fraction .* z0snavez(te*1000));  %Tungsten
    z1 = z1.*ones(size(te,1),size(te,2));  
    z2 = z2.*ones(size(te,1),size(te,2));

    if qlkparams.calc_part_transport == 0
        profil_var = profil;

        n1 = n_loop(:,:,2);
        profil_var.nzp = n1;

        nW = n_loop(:,:,4);
        profil_var.nwp = nW;  

        ne   = n_loop(:,:,end); 
        
        n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5));
        ni = n_loop(:,:,1);
        profil_var.nip = ni;
    else
        profil_var = profil;

        n1 = n_loop(:,:,2) .* (1 + recuit .* f .* rand(size(RLFS)));
        profil_var.nzp = n1;

        nW = n_loop(:,:,4) .* (1 + recuit .* f .* rand(size(RLFS)));
        profil_var.nwp = nW;  

        ne   = n_loop(:,:,end) .* (1 + recuit .* f .* rand(size(RLFS)));     
        
        n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5))./(1+option.cmin);  %Quasi-neutrality
        ni = n_loop(:,:,1);
        profil_var.nip = ni;   
    end

    
    % calcul des coefficients de transport
    % Zeff passed to QuaLiKiz should be low Z contributions only because in QLKNN trained with Carbon species (will affect dilution impact)

    if qlkparams.imp_trans == 4   %auto consistent QLKZ + NEO (facitparams,option,z0dinput, profil, tmp, itime,fac)
        yneo_out = neo_std(neoparams,option,z0dinput,profil_var,itime,tmp);
        G_n = yneo_out.G; DD_loop_n = yneo_out.D; VV_loop_n = yneo_out.V;
        DD_loop_n = repmat(DD_loop_n, 4, 1, 1);
        VV_loop_n = repmat(VV_loop_n, 4, 1, 1);
        G_n = repmat(G_n, 4, 1, 1);

        [chie_loop,chii_loop,DD_loop,VV_loop,G] = qlkz_std(profil.xli,te,ti,ne,profil_var.web,profil_var.zeff,qr,factor_em,meff,geo.r0,geo.a,geo.b0,Amain,qlkparams,profil_var,option,zerod);

    elseif qlkparams.imp_trans == 3   %auto consistent QLKZ + FACIT
        [G_n,DD_loop_n,VV_loop_n] = facit_std(facitparams,option,z0dinput, profil_var, tmp, itime,1);
        DD_loop_n(:,:,5) = 0;
        VV_loop_n(:,:,5) = 0;
        G_n(:,:,5) = 0;
        DD_loop_n = repmat(DD_loop_n, 4, 1, 1);
        VV_loop_n = repmat(VV_loop_n, 4, 1, 1);
        G_n = repmat(G_n, 4, 1, 1);

        [chie_loop,chii_loop,DD_loop,VV_loop,G] = qlkz_std(profil.xli,te,ti,ne,profil_var.web,profil_var.zeff,qr,factor_em,meff,geo.r0,geo.a,geo.b0,Amain,qlkparams,profil_var,option,zerod);

    elseif qlkparams.imp_trans == 2   %auto consistent QLKZ
        [chie_loop,chii_loop,DD_loop,VV_loop,G] = qlkz_std(profil.xli,te,ti,ne,profil_var.web,profil_var.zeff,qr,factor_em,meff,geo.r0,geo.a,geo.b0,Amain,qlkparams,profil_var,option,zerod);

    elseif qlkparams.imp_trans == 1   %raw QLKZ
        profil_var.nip  = n_init(:,:,1);
        profil_var.nzp  = n_init(:,:,2);
        profil_var.nwp  = n_init(:,:,4);
        ne              = n_init(:,:,end);
        [chie_loop,chii_loop,DD_raw,VV_raw,G_raw] = qlkz_std(profil.xli,te,ti,ne,profil_var.web,profil_var.zeff,qr,factor_em,meff,geo.r0,geo.a,geo.b0,Amain,qlkparams,profil_var,option,zerod);
    else 
    end
    %Densities
    ne_loop = n_loop(:,:,end);
    n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5));

    % ======= Figure des profils de densité =======
    if ~exist('hz', 'var') || ~isgraphics(hz)
        hz = figure('tag', 'Prof_n_species', ...
                    'name', 'Density profiles', ...
                    'color', [1 1 1]);
    else
        figure(hz);
    end
    
    clf(hz);  % Efface proprement la figure

    set(hz, 'defaultaxesfontsize', 18, ...
            'defaultaxesfontweight', 'bold', ...
            'defaultaxesfontname', 'times', ...
            'defaultlinelinewidth', 2);
    set(hz, 'Position', [100, 100, 1000, 700]);

    colors = lines(6);
    species_names = {'Main Ions', 'Imp 1', 'Imp 2', 'W', 'H', 'e'};

    for i = 1:6
        ax = subplot(3,2,i); cla(ax); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

        % METIS (trait plein)
        zplotprof(ax, profil.temps, profil.xli, n_init(:,:,i), ...
                'color', colors(i,:), 'linestyle', '-', 'linewidth', 2);

        % Loop (pointillé)
        zplotprof(ax, profil.temps, profil.xli, n_loop(:,:,i), ...
                'color', colors(i,:), 'linestyle', ':', 'linewidth', 2);

        xlabel(ax,'r / a');
        ylabel(ax,'n (m^{-3})');
        title(ax, species_names{i});
        set(ax, 'yscale','linear');
        zoom(ax,'yon');

        % -- Récupération des handles du 1er subplot pour la légende --
        if i == 1
            ch = findobj(ax, 'Type', 'line');
            h_loop  = ch(1);  % pointillé (n_loop)
            h_metis = ch(2);  % plein (n_init)
        end
    end

    % Légende commune
    hl = legend([h_metis, h_loop], ...
                {'METIS (n\_init)', 'Loop (n\_loop)'}, ...
                'Orientation', 'horizontal', ...
                'FontSize', 14, 'Box', 'off');
    set(hl, 'Units', 'normalized', ...
            'Position', [0.35 0.02 0.3 0.04]);  % ajuster si besoin

    drawnow;


    QLK_nn_data.Chi_e_QLK = chie_loop;
    QLK_nn_data.Chi_i_QLK = chii_loop;

    QLK_nn_data.D_n_QLK  = DD_loop;
    QLK_nn_data.V_n_QLK  = VV_loop;   
    QLK_nn_data.D_n_QLK_n  = DD_loop_n;
    QLK_nn_data.V_n_QLK_n  = VV_loop_n;  
    QLK_nn_data.D_ne_QLK  = D_loop;
    QLK_nn_data.V_ne_QLK  = V_loop;   

    %
    if qlkparams.verbosity > 5
        figure(21);
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
    filter_bis = ~isfinite(DD_loop(:)) | ~isfinite(VV_loop(:)) | ~isfinite(DD_loop_n(:)) | ~isfinite(VV_loop_n(:));
    if any(filter)
        disp('bad values in NN results');
    end
    if any(filter_bis)
        disp('bad values in NN results');
    end
    
% PM MODIF
%disp('WARNING: modifierd MEM LOOP')
%    chii_loop(filter) = chie_mem(filter);
%    chie_loop(filter) = chii_mem(filter);
%    D_loop(filter)    = V_mem(filter);
%    V_loop(filter)    = D_mem(filter);

    chii_loop(filter) = chii_mem(filter);
    chie_loop(filter) = chie_mem(filter);
    D_loop(filter)    = D_mem(filter);
    V_loop(filter)    = V_mem(filter);

    DD_loop(filter_bis)    = DD_mem(filter_bis);
    VV_loop(filter_bis)    = VV_mem(filter_bis);

    DD_loop_n(filter_bis)    = DD_mem_n(filter_bis);
    VV_loop_n(filter_bis)    = VV_mem_n(filter_bis);

    
    % test = (chie_loop > qlkparams.maxChie) | (chii_loop > qlkparams.maxChii) |  ...
    %     (D_loop(:,:,end) > qlkparams.maxDe) | (V_loop(:,:,end) > qlkparams.maxVe) | (V_loop(:,:,end) < qlkparams.minVe);
    % overshoot = sum(double(test(:)))./ numel(test);

    test = (chie_loop > qlkparams.maxChie) | (chii_loop > qlkparams.maxChii) |  ...
        (DD_loop > qlkparams.maxDe) | (VV_loop > qlkparams.maxVe) | (VV_loop < qlkparams.minVe);
    overshoot = sum(double(test(:)))./ numel(test);
    
    chie_loop(chie_loop > qlkparams.maxChie) = chie_mem(chie_loop > qlkparams.maxChie);
    chii_loop(chii_loop > qlkparams.maxChii) = chie_mem(chii_loop > qlkparams.maxChii);
    D_loop(D_loop > qlkparams.maxDe)         = D_mem(D_loop > qlkparams.maxDe);
    V_loop(V_loop > qlkparams.maxVe)         = V_mem(V_loop > qlkparams.maxVe);

    DD_loop(DD_loop > qlkparams.maxDe)         = DD_mem(DD_loop > qlkparams.maxDe);
    VV_loop(VV_loop > qlkparams.maxVe)         = VV_mem(VV_loop > qlkparams.maxVe);

    chie_loop(chie_loop < qlkparams.minChie) = qlkparams.minChie;
    chii_loop(chii_loop < qlkparams.minChii) = qlkparams.minChii;
    V_loop(V_loop < qlkparams.minVe)         = qlkparams.minVe;
    D_loop(D_loop < qlkparams.minDe)         = qlkparams.minDe;

    VV_loop(VV_loop < qlkparams.minVe)         = qlkparams.minVe;
    DD_loop(DD_loop < qlkparams.minDe)         = qlkparams.minDe;
    
    % Quasi-neutrality
    n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5));

    % neoclassical
    % neoclassical part (use the same on ion and electron)
    profil_loc = profil_var;
    profil_loc.tep = te_loop;
    profil_loc.tip = ti_loop;
    profil_loc.nep = n_loop(:,:,end);

    profil_loc.nHp    = n_loop(:,:,1) .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);
    profil_loc.nDp    = n_loop(:,:,1) .* ((nDm ./ n1m) * ve);
    profil_loc.nTp    = n_loop(:,:,1) .* ((nTm ./ n1m) * ve);
    profil_loc.nzp    = n_loop(:,:,2);
    profil_loc.nwp    = n_loop(:,:,4);
    switch option.gaz
        case 5
            profil_loc.nBp   = zeros(size(nTp));
            profil_loc.nhep  = option.frhe0 .* n_loop(:,:,end);
            profil_loc.nhe3p = profil_loc.nhep;
        case 11
            profil_loc.nBp = profil_loc.nTp;
            profil_loc.nTp = zeros(size(nTp));
            profil_loc.nhe3p = zeros(size(nTp));
        otherwise
            profil_loc.nBp = zeros(size(nTp));
            profil_loc.nhe3p = zeros(size(nTp));
    end
    
    % change on other species
    % profil_loc.nHp   = nHp  .*  ne_loop ./ max(1,profil.nep);
    % profil_loc.nDp   = n_loop(:,:,1);
    % profil_loc.nTp   = nTp  .*  ne_loop ./ max(1,profil.nep);
    % profil_loc.nBp   = nBp  .*  ne_loop ./ max(1,profil.nep);
    % profil_loc.nhep  = nhep .*  ne_loop ./ max(1,profil.nep);
    % profil_loc.nhe3p = nhe3p .*  ne_loop ./ max(1,profil.nep);
    % profil_loc.nzp   = n_loop(:,:,2);
    % profil_loc.nwp   = n_loop(:,:,1);
    profil_loc.n1p   = profil_loc.nHp + profil_loc.nDp + profil_loc.nTp;
    profil_loc.nip   = profil_loc.n1p; %+ profil_loc.nhep + profil_loc.nzp + profil_loc.nwp;
    chii_neo  = qlkparams.factor_neo_chii .* neo_Hinton(profil_loc,option,Amain);
    chii_neo  = chii_neo .* profil.grho2 ./ profil.grho;
    VWare_neo = Vware_ref; % not updated in this first version
    
    % sum
    chii_loop = chii_loop + chii_neo;
    V_loop    = V_loop + VWare_neo;
    VV_loop(:,:,end) = VV_loop(:,:,end) + VWare_neo;
    D_loop    = max(D_loop,qlkparams.neo_ion_in_D .* chii_neo);
    chie_loop = max(chie_loop,qlkparams.neo_ion_in_Chie .* chii_neo);
    
    % space filter
    if qlkparams.smooth_space
        for lk = 1:length(profil.temps)
            chie_loop(lk,:) = sgolayfilt(chie_loop(lk,:),1,3);
            chii_loop(lk,:) = sgolayfilt(chii_loop(lk,:),1,3);
            D_loop(lk,:)    = sgolayfilt(D_loop(lk,:),1,3);
            V_loop(lk,:)    = sgolayfilt(V_loop(lk,:),1,3);
            
            for p = 1:6
                DD_loop(lk,:,p)    = sgolayfilt(DD_loop(lk,:,p),1,3);
                VV_loop(lk,:,p)    = sgolayfilt(VV_loop(lk,:,p),1,3);

                DD_loop_n(lk,:,p)    = sgolayfilt(DD_loop_n(lk,:,p),1,3);
                VV_loop_n(lk,:,p)    = sgolayfilt(VV_loop_n(lk,:,p),1,3);
            end
        end
    end
    
    % masquage  (only evolved part of the profile is changed)
    chie_loop = (1 - mask) .* chie_ref + mask .* chie_loop;
    chii_loop = (1 - mask) .* chii_ref + mask .* chii_loop;
    D_loop    = (1 - mask) .* D_ref    + mask .* D_loop;
    V_loop    = (1 - mask) .* V_ref    + mask .* V_loop;
    % DD_loop = (1 - mask) .* repmat(D_ref, [1, 1, 6])  + mask .* DD_loop;
    % VV_loop = (1 - mask) .* repmat(V_ref, [1, 1, 6])  + mask .* VV_loop;

    DD_loop(:,:,end) = (1 - mask) .* DD_ref(:,:,end)  + mask .* DD_loop(:,:,end);
    VV_loop(:,:,end) = (1 - mask) .* VV_ref(:,:,end)  + mask .* VV_loop(:,:,end);
    

    % amorti
    chie_loop = (1 - f) .* chie_mem + f .* chie_loop;
    chii_loop = (1 - f) .* chii_mem + f .* chii_loop;
    D_loop    = (1 - f) .* D_mem    + f .* D_loop;
    V_loop    = (1 - f) .* V_mem    + f .* V_loop;

    DD_loop    = (1 - f) .* DD_mem    + f .* DD_loop;
    VV_loop    = (1 - f) .* VV_mem    + f .* VV_loop;


    DD_loop_n    = (1 - f) .* DD_mem_n    + f .* DD_loop_n;
    VV_loop_n    = (1 - f) .* VV_mem_n    + f .* VV_loop_n;

    %
    QLK_nn_data.Chi_e = chie_loop;
    QLK_nn_data.Chi_i = chii_loop;
    QLK_nn_data.D_ne  = D_loop;
    QLK_nn_data.V_ne  = V_loop;

    QLK_nn_data.D_n  = DD_loop;
    QLK_nn_data.V_n  = VV_loop;

    QLK_nn_data.D_n_n  = DD_loop_n;
    QLK_nn_data.V_n_n  = VV_loop_n;
    %

    % calcul des profils
    warning off
    % missing flux due to particles (3/2 or 5/2 ?)
    % with fixed sources including  (ion-e heat flux exchange due to collisions)
    ted1_loop    = - (Qe - Qei - Qrad - 5/2 .* phys.e .* te_loop .* ge_sts) ./ (chie_loop .* n_loop(:,:,end) .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
    % ted1_loop    = - (Qe - Qei - 5/2 .* phys.e .* te_loop .* ge_sts) ./ (chie_loop .* n_loop(:,:,end) .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);
    tid1_loop    = - (Qi + Qei - 5/2 .* phys.e .* ti_loop .* gi_sts) ./ (chii_loop .* n_loop(:,:,1) .* phys.e .* profil.vpr .* profil.grho ./ rhomax .^ 2);

    %Radiated power by impurities
    fields = fieldnames(atom);
    z1_n = fields{cell2mat(cellfun(@(f) atom.(f).Z == option.zimp, fields, 'UniformOutput', false))};
    z2_n = fields{cell2mat(cellfun(@(f) atom.(f).Z == option.zmax, fields, 'UniformOutput', false))};


    te_adas_tab = [atom.(z1_n).data(:,1), atom.(z2_n).data(:,1), atom.W.data(:,1)] * 1e3;  % [eV]
    lz_adas_tab = [atom.(z1_n).data(:,2), atom.(z2_n).data(:,2), atom.W.data(:,2)];
    Lz_tab = zeros(size(te_loop,1), size(te_loop,2), 3); 
    for i = 1:3
        Lz_tab(:,:,i) = interp1(te_adas_tab(:,i), lz_adas_tab(:,i), te_loop, 'linear', 'extrap');
    end

    Prad_loop = n_loop(:,:,end) .* n_loop(:,:,2:4) .* Lz_tab;  
    Prad = sum(Prad_loop, 3);
    Qrad_loop = zeros(size(Prad_loop));  % Pré-allocation
    for j = 1:3
        Qrad_loop(:,:,j) = cumtrapz(profil.xli, Prad_loop(:,:,j).* profil.vpr, 2);
    end
    Qrad = sum(Qrad_loop, 3);

    %Neoclassical fluxes
    if qlkparams.imp_trans == 4
        G_loop_n(:,:,(size(n_loop,3))) = - VV_loop_n(:,:,end) .* n_loop(:,:,end);
    else 
        G_loop_n(:,:,(size(n_loop,3))) = - Vware .* n_loop(:,:,end);
    end
    G_loop_n(:,:,2:(size(n_loop,3)-1)) = - VV_loop_n(:,:,2:(end-1)) .* n_loop(:,:,2:(end-1));
    %Turbulent fluxes
    ge_loop  = ge_sts - V_loop .* n_loop(:,:,end);
    G_loop(:,:,(size(n_loop,3))) = ge_sts - VV_loop(:,:,end) .* n_loop(:,:,end); 
    G_loop(:,:,2:(size(n_loop,3)-1)) = - VV_loop(:,:,2:(end-1)) .* n_loop(:,:,2:(end-1));


    %Gradients
    ned1_loop   = - ge_loop ./ D_loop .* rhomax;
    nd1_loop(:,:,2:(size(n_loop,3))) = - (G_loop(:,:,2:(size(n_loop,3))) + G_loop_n(:,:,2:(size(n_loop,3)))) ./ (DD_loop(:,:,2:(size(n_loop,3))) + DD_loop_n(:,:,2:(size(n_loop,3)))) .* rhomax;
    nd1_loop(:,:,1) = nd1_loop(:,:,end) - z1.*nd1_loop(:,:,2)-z2.*nd1_loop(:,:,3)-z3.*nd1_loop(:,:,4)-nd1_loop(:,:,5);
    warning on

    % maximum gradient (to prevent overshoot)
    % overgradient = double(mask)  .* double((abs(ted1_loop) > (tdmax/sqrt(f))) | (abs(tid1_loop) > (tdmax/sqrt(f))) | (abs(ned1_loop) > (ndmax/sqrt(f))));
    % overgradient = sum(overgradient(:)) ./ numel(overgradient(:));

    overgradient = double(mask)  .* double((abs(ted1_loop) > (tdmax/sqrt(f))) | (abs(tid1_loop) > (tdmax/sqrt(f))) | (abs(nd1_loop) > (ndmax/sqrt(f))));
    overgradient = sum(overgradient(:)) ./ numel(overgradient(:));
    
    ted1_loop = max(-tdmax/sqrt(f),min(tdmax/sqrt(f),ted1_loop));
    ned1_loop = max(-ndmax/sqrt(f),min(ndmax/sqrt(f),ned1_loop));
    tid1_loop = max(-tdmax/sqrt(f),min(tdmax/sqrt(f),tid1_loop));

    nd1_loop = max(-ndmax/sqrt(f),min(ndmax/sqrt(f),nd1_loop));
    
    % security
    ted1_loop(~isfinite(ted1_loop)) =ted1(~isfinite(ted1_loop));
    tid1_loop(~isfinite(tid1_loop)) =tid1(~isfinite(tid1_loop));
    ned1_loop(~isfinite(ned1_loop)) =ned1(~isfinite(ned1_loop));

    nd1_loop(~isfinite(nd1_loop)) =nd1(~isfinite(nd1_loop));
    
    % masquage  (only evolved part of the profile is changed)
    ted1_loop = mask .* ted1_loop + (1 - mask) .* ted1;
    tid1_loop = mask .* tid1_loop + (1 - mask) .* tid1;
    ned1_loop = mask .* ned1_loop + (1 - mask) .* ned1;

    nd1_loop = mask .* nd1_loop + (1 - mask) .* nd1; 
    
    % symetry
    ted1_loop(:,1) = 0;
    tid1_loop(:,1) = 0;
    ned1_loop(:,1) = 0;

    nd1_loop(:,1,:) = 0;
    
    % integration
    te_loop = cumtrapz(profil.xli(end:-1:1),ted1_loop(:,end:-1:1),2);
    te_loop = te_loop(:,end:-1:1);
    ti_loop = cumtrapz(profil.xli(end:-1:1),tid1_loop(:,end:-1:1),2);
    ti_loop = ti_loop(:,end:-1:1);
    ne_loop = cumtrapz(profil.xli(end:-1:1),ned1_loop(:,end:-1:1),2);
    ne_loop = ne_loop(:,end:-1:1);

    for j = 1:size(n_loop,3)
        n_loop(:,:,j) = cumtrapz(profil.xli(end:-1:1), nd1_loop(:,end:-1:1,j), 2);
        n_loop(:,:,j) = n_loop(:,end:-1:1,j); 
    end

    % Quasi-neutrality
    n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)+n_loop(:,:,5));
    % raccord
    te_loop = te_loop + (profil.tep(:,indedge) - te_loop(:,indedge)) * ones(size(profil.xli));
    te_loop(:,indedge:end) = profil.tep(:,indedge:end);
    ti_loop = ti_loop + (profil.tip(:,indedge) - ti_loop(:,indedge)) * ones(size(profil.xli));
    ti_loop(:,indedge:end) = profil.tip(:,indedge:end);
    ne_loop = ne_loop + (profil.nep(:,indedge) - ne_loop(:,indedge)) * ones(size(profil.xli));
    ne_loop(:,indedge:end) = profil.nep(:,indedge:end);
    for j = 1:size(n_loop, 3)
        n_loop(:,:,j) = n_loop(:,:,j) + (n_init(:,indedge,j) - n_loop(:,indedge,j)) * ones(size(profil.xli));
        n_loop(:,indedge:end,j) = n_init(:,indedge:end,j);
    end

    % Quasi-neutrality
    n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5));

    % valeur minimum (boundary condition)
    te_loop = max(te_loop,profil.tep(:,end) * ones(size(profil.xli)));
    ti_loop = max(ti_loop,profil.tip(:,end) * ones(size(profil.xli)));
    ne_loop = max(ne_loop,profil.nep(:,end) * ones(size(profil.xli)));
    
    n_loop = max(n_loop,n_init(:,end,:) * ones(size(profil.xli, 1)));
    %figure(21);plot(profil.xli,ne_loop);drawnow
    
    % amorti
    te_loop = (1 - f) .* te_mem + f .* te_loop;
    ti_loop = (1 - f) .* ti_mem + f .* ti_loop;
    ne_loop = (1 - f) .* ne_mem + f .* ne_loop;
    
    n_loop = (1 - f) .* n_mem + f .* n_loop;
    
    % Quasi-neutrality
    n_loop(:,:,1) = (n_loop(:,:,end) - z1.*n_loop(:,:,2)-z2.*n_loop(:,:,3)-z3.*n_loop(:,:,4)-n_loop(:,:,5));

    % synchronistion with option in QLK-NN
    if qlkparams.calc_heat_transport == 0
	te_loop = te_mem;
	ti_loop = ti_mem;
    end
    if qlkparams.calc_part_transport == 0
	ne_loop = ne_mem;
    n_loop = n_mem;
    end
    % to turn off density transport
    %  ne_loop = ne_mem;
    % to turn off heqt transport
    % te_loop = te_mem;
    % ti_loop = ti_mem;

    % new Qei value
    nHp_loop    = n_loop(:,:,1) .* ((max(0,n1m - nDm - nTm) ./ n1m) * ve);
    nDp_loop    = n_loop(:,:,1) .* ((nDm ./ n1m) * ve);
    nTp_loop    = n_loop(:,:,1) .* ((nTm ./ n1m) * ve);
    nzp_loop    = n_loop(:,:,2);
    nwp_loop    = n_loop(:,:,4);
    switch option.gaz
        case 5
            nBp_loop   = zeros(size(nTp));
            nhep_loop  = option.frhe0 .* n_loop(:,:,end);
            nhe3p_loop = nhep_loop;
        case 11
            nBp_loop = nTp_loop;
            nTp_loop = zeros(size(nTp));
            nhe3p_loop = zeros(size(nTp));
        otherwise
            nhep_loop  = option.frhe0 .* n_loop(:,:,end);
            nBp_loop = zeros(size(nTp));
            nhe3p_loop = zeros(size(nTp));
    end

    % terme source
    switch  extended_qei
        case 'on'
            qeib0 = equipartition_full(te_loop,ti_loop,n_loop(:,:,end),nHp_loop,nDp_loop,nTp_loop,nhep_loop,nzp_loop, ...
                option.rimp .* nzp,(1 - Sn_fraction) .* nwp_loop, Sn_fraction .* nwp_loop,option.zimp,option.zmax,nBp_loop,nhe3p_loop);
        otherwise
            
            warning off
            lnei          =  15.2 - 0.5.*log(n_loop(:,:,end) ./ 1e20) + log(te_loop ./ 1e3);
            warning on
            ind = find(~isfinite(lnei) | (lnei <10));
            if ~isempty(ind)
                lnei(ind) = 10 .* ones(1,length(ind));
            end
            taues  = (12 .* pi .^ (3/2) ./ sqrt(2)) .* (phys.epsi0 .^ 2 ./ phys.e .^ 4 .* sqrt(phys.me))  .* ...
                ((phys.e .* max(min(profil.tep(:)),te_loop)) .^ (3/2) ./ n_loop(:,:,end) ./ lnei);
            
            if Sn_fraction > 0
                qeib0    = 3 .* phys.me ./ phys.mp ./ taues .* (nHp_loop + nDp_loop ./ 2 + nTp_loop ./ 3 + nhep_loop + nhe3p_loop .* (4/3.02) +  ...
                    (option.zimp .^ 2 ./ Aimp) .* nzp_loop + option.rimp .*  (option.zmax .^ 2 ./ Amax) .* nzp_loop +  ...
                    (1 - Sn_fraction) .* z0wavez(te_loop) ./ 183.84 .* nwp_loop + ...
                    Sn_fraction .* z0snavez(te_loop) ./ 118.71 .* nwp_loop ) + 25/11 .* nBp_loop;
            else
                qeib0    = 3 .* phys.me ./ phys.mp ./ taues .* (nHp_loop + nDp_loop ./ 2 + nTp_loop ./ 3 + nhep_loop + nhe3p_loop .* (4/3.02) +  ...
                    (option.zimp .^ 2 ./ Aimp) .* nzp_loop + option.rimp .*  (option.zmax .^ 2 ./ Amax) .* nzp_loop + z0wavez(te_loop) ./ 183.84 .* nwp_loop) + 25/11 .* nBp_loop;
            end
    end
    qeib     = qeib0 .* phys.e .* (te_loop - ti_loop);
    Qei = cumtrapz(profil.xli,qeib .* profil.vpr,2);
    
    %Radiated power by impurities
    fields = fieldnames(atom);
    z1_n = fields{cell2mat(cellfun(@(f) atom.(f).Z == option.zimp, fields, 'UniformOutput', false))};
    z2_n = fields{cell2mat(cellfun(@(f) atom.(f).Z == option.zmax, fields, 'UniformOutput', false))};


    te_adas_tab = [atom.(z1_n).data(:,1), atom.(z2_n).data(:,1), atom.W.data(:,1)] * 1e3;  % [eV]
    lz_adas_tab = [atom.(z1_n).data(:,2), atom.(z2_n).data(:,2), atom.W.data(:,2)];
    Lz_tab = zeros(size(te_loop,1), size(te_loop,2), 3); 
    Lz_mem_tab = Lz_tab;
    Lz_init_tab = Lz_tab;
    for i = 1:3
        Lz_tab(:,:,i) = interp1(te_adas_tab(:,i), lz_adas_tab(:,i), te_loop, 'linear', 'extrap');
        Lz_mem_tab(:,:,i) = interp1(te_adas_tab(:,i), lz_adas_tab(:,i), te_mem, 'linear', 'extrap');
        Lz_init_tab(:,:,i) = interp1(te_adas_tab(:,i), lz_adas_tab(:,i), profil.tep, 'linear', 'extrap');
    end

    Prad_loop = n_loop(:,:,end) .* n_loop(:,:,2:4) .* Lz_tab;
    Prad_mem  = n_mem(:,:,end)  .* n_mem(:,:,2:4)  .* Lz_mem_tab;
    Prad_init = n_init(:,:,end) .* n_init(:,:,2:4) .* Lz_init_tab;
    % calcul de la stabilite de base
    err_rad_z1 = std((Prad_loop(:,:,1) - Prad_mem(:,:,1)) ./ mean(Prad_init(:,:,1), 'all'), 0, 'all');
    err_rad_z2 = std((Prad_loop(:,:,2) - Prad_mem(:,:,2)) ./ mean(Prad_init(:,:,2), 'all'), 0, 'all');
    err_rad_W = std((Prad_loop(:,:,3) - Prad_mem(:,:,3)) ./ mean(Prad_init(:,:,3), 'all'), 0, 'all');


    err_n = std((n_loop(:,:,1) - n_mem(:,:,1)) ./ mean(n_init(:,:,1),'all'), 0, 'all');
    err_T = std((te_loop(:) - te_mem(:)) ./ mean(profil.tep(:),'all'), 0, 'all') + std((ti_loop(:) - ti_mem(:)) ./ mean(profil.tip(:),'all'), 0, 'all');
    err_chi = std(chie_loop(:) - chie_mem(:)) +  std(chii_loop(:) - chii_mem(:));
    err_D = std((DD_loop(:,:,:) - DD_mem(:,:,:)), 0, 'all');
    err_D_n = std((DD_loop_n(:,:,2:(end-1)) - DD_mem_n(:,:,2:(end-1))), 0, 'all');
    err_stab = err_n + (err_rad_W + err_rad_z1 + err_rad_z2) + err_T + err_chi + err_D + err_D_n;
    
        
%    if err_stab < max(1e-3,option.tol0d)
    if [err_stab] < max(1e-3, option.tol0d)
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

    DD_mem = DD_loop;
    VV_mem = VV_loop;

    DD_mem_n = DD_loop_n;
    VV_mem_n = VV_loop_n;
    n_mem = n_loop;    

    err_stab_vec(end+1) = err_stab;
    err_n_vec(end+1) = err_n;
    err_rad_W_vec(end+1) = err_rad_W;
    err_rad_z1_vec(end+1) = err_rad_z1;
    err_rad_z2_vec(end+1) = err_rad_z2;
    err_T_vec(end+1) = err_T;
    err_chi_vec(end+1) = err_chi;
    err_D_vec(end+1) = err_D;
    err_D_n_vec(end+1) = err_D_n;


    % ======= Figure des erreurs de stabilité =======
    if ~exist('hf_errstab', 'var') || ~isgraphics(hf_errstab)
        hf_errstab = figure('tag','Fig_ErrStab', ...
                            'name','Stability error', ...
                            'color',[1 1 1]);
    else
        figure(hf_errstab);
    end
    
    clf(hf_errstab);  % Efface proprement la figure ciblée

    axes('yscale', 'log');  % Axe Y en échelle log
    hold on; box on; grid on;
    set(gca, 'fontsize', 14, 'fontweight', 'bold', 'fontname', 'times');
    xlabel('k');
    ylabel('err\_stab');
    title('Convergence of err\_stab');

    semilogy(1:k, err_n_vec, '-og', 'LineWidth', 2);
    semilogy(1:k, err_rad_W_vec, '-or', 'LineWidth', 2);
    semilogy(1:k, err_rad_z1_vec, '-oc', 'LineWidth', 2);
    semilogy(1:k, err_rad_z2_vec, '-om', 'LineWidth', 2);
    semilogy(1:k, err_T_vec, '-ob', 'LineWidth', 2);
    semilogy(1:k, err_D_vec, '--og', 'LineWidth', 2);
    semilogy(1:k, err_D_n_vec, ':og', 'LineWidth', 2);
    semilogy(1:k, err_chi_vec, '--ob', 'LineWidth', 2);
    semilogy(1:k, err_stab_vec, '-ok', 'LineWidth', 2);

    legend('\sigma_n', '\sigma_{rad,W}', '\sigma_{rad,z1}', '\sigma_{rad,z2}', ...
        '\sigma_T', '\sigma_D', '\sigma_{D_n}', '\sigma_{\chi}', '\sigma_{tot}');
    drawnow;

    
    % ======= Figure Te / Ti =========
    if ~exist('hf_TeTi', 'var') || ~isgraphics(hf_TeTi)
        hf_TeTi = figure('tag','Prof_TeTi', ...
                         'name','Profils Te / Ti – METIS vs loop', ...
                         'color',[1 1 1]);
    else
        figure(hf_TeTi);
    end    
    clf(hf_TeTi);  % Efface proprement cette figure

    set(hf_TeTi, 'defaultaxesfontsize', 18, ...
                'defaultaxesfontweight', 'bold', ...
                'defaultaxesfontname', 'times', ...
                'defaultlinelinewidth', 2);

    set(hf_TeTi, 'Position', [150 150 900 600]);

    % -- Subplot 1 : Ti
    ax1 = subplot(2,1,1); cla(ax1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
    zplotprof(ax1, profil.temps, profil.xli, profil.tip ./ 1e3, ...
            'color', colIon, 'linestyle','-', 'linewidth',2);
    zplotprof(ax1, profil.temps, profil.xli, ti_loop ./ 1e3, ...
            'color', colIon, 'linestyle',':', 'linewidth',2);
    ylabel(ax1, 'T_i (keV)');
    title(ax1, 'Ions temperature');
    legend(ax1, {'METIS','Loop'}, 'location','best', 'fontsize', 12);

    % -- Subplot 2 : Te
    ax2 = subplot(2,1,2); cla(ax2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
    zplotprof(ax2, profil.temps, profil.xli, profil.tep ./ 1e3, ...
            'color', colElec, 'linestyle','-', 'linewidth',2);
    zplotprof(ax2, profil.temps, profil.xli, te_loop ./ 1e3, ...
            'color', colElec, 'linestyle',':', 'linewidth',2);
    xlabel(ax2, 'r / a');
    ylabel(ax2, 'T_e (keV)');
    title(ax2, 'Electrons temperature');
    legend(ax2, {'METIS','Loop'}, 'location','best', 'fontsize', 12);

    drawnow;


    QLK_nn_data.N    = n_loop;
    
end
% fin de la Loop
fprintf('converged in %d loops\n',k);
% remove wait bar
if first_call == 1
    delete(hwaitbar);
end
mask_err = mask;
mask(mask == 0) = NaN;
mask(mask > 0.5) = 0;


% Computation of the dominant mode frequency
path_to_METIS_QLKZ = fullfile(fileparts(which('install_QLKNN10D')));
username = getenv('USER');
runname_qlk=['METIS_QLKZ_tmp_' username];
runpath_qlk = [path_to_METIS_QLKZ '/runs_qlkz/from_METIS/'];
tmp_dir = tempdir;
temp_QLKZ = struct_qlkz(runname_qlk,tmp_dir,0);

kperp = reshape(temp_QLKZ.kperp2, size(n_init(:,:,end),1), size(n_init(:,:,end),2), []);
ome = reshape(temp_QLKZ.ome_GB, size(n_init(:,:,end),1), size(n_init(:,:,end),2), qlkparams.numsols, []);  % ome : [4 × 21 × 3 × 10]
gam = reshape(temp_QLKZ.gam_GB, size(n_init(:,:,end),1), size(n_init(:,:,end),2), qlkparams.numsols, []);  % gam : [4 × 21 × 3 × 10]

gam_kperp = gam ./ reshape(kperp, size(n_init(:,:,end),1), size(n_init(:,:,end),2), 1, []);  % taille résultante : [4 × 21 × 3 × 10]

omega_max = zeros(size(n_init(:,:,end),1), size(n_init(:,:,end),2));  % Pour stocker les valeurs de ome au maximum

for i = 1:size(n_init(:,:,end),1)
    for j = 1:size(n_init(:,:,end),2)
        slice = squeeze(gam_kperp(i,j,:,:));  % [3 × 10]
        [max_val, idex] = max(slice(:));
        if max_val == 0
            omega_max(i,j) = NaN;
        else
            [mode_idx, kp_idx] = ind2sub([qlkparams.numsols, size(temp_QLKZ.kthetarhos)], idex);
            omega_max(i,j) = ome(i,j,mode_idx,kp_idx);
        end
    end
end

% PLOTS
if ismember(qlkparams.imp_trans, [2, 3, 4])
    G_loop(:,:,2:end) = - DD_loop(:,:,2:end) .* nd1_loop(:,:,2:end) ./ rhomax + VV_loop(:,:,2:end) .* n_loop(:,:,2:end);
    G_loop(:,:,1) = G_loop(:,:,end) - z1.*G_loop(:,:,2) - z2.*G_loop(:,:,3) - z3.*G_loop(:,:,4) - G_loop(:,:,5); 

    if ismember(qlkparams.imp_trans, [3, 4])
        G_loop_n(:,:,2:end) = - DD_loop_n(:,:,2:end) .* nd1_loop(:,:,2:end) ./ rhomax + VV_loop_n(:,:,2:end) .* n_loop(:,:,2:end);
        G_loop_n(:,:,1) = G_loop_n(:,:,end) - z1.*G_loop_n(:,:,2) - z2.*G_loop_n(:,:,3) - z3.*G_loop_n(:,:,4) - G_loop_n(:,:,5); 
    else
    end


    % === Figure -V/D ===
    if ~exist('hf_VD', 'var') || ~isgraphics(hf_VD)
        hf_VD = figure('tag','Fig_VD', 'name','-V/D Profile', 'color',[1 1 1]);
    else
        figure(hf_VD);
    end
    clf(hf_VD);
    hold on; box on; grid on;

    colors = lines(6);
    labels = {'$1$', '$2$', '$W$', '$H$'};
    for i = 1:4
        idx = i + 1;
        plot(profil.xli, -VV_loop(2,:,idx)./DD_loop(2,:,idx), '-', ...
            'Color', colors(idx,:), 'LineWidth', 1.5, 'Marker', 'o');
    end
    xlabel('r/a', 'FontSize', 14);
    ylabel('-V/D', 'FontSize', 14, 'Interpreter', 'latex');
    legend(labels, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);




    % === Figure Q_s vs Q_rad ===
    if ~exist('hf_Qrad', 'var') || ~isgraphics(hf_Qrad)
        hf_Qrad = figure('tag','Fig_Qrad', 'name','Q_s vs Q_rad', 'color',[1 1 1]);
    else
        figure(hf_Qrad);
    end
    clf(hf_Qrad);
    hold on; box on; grid on;

    plot(profil.xli, Qe(2,:), '-', 'LineWidth', 1.5, 'Marker', 'o');
    if qlkparams.imp_trans >= 3
        for i = 1:3
            plot(profil.xli, Qrad_loop(2,:,i), '--', 'LineWidth', 1.2);
        end
        plot(profil.xli, Qrad(2,:), '--k', 'LineWidth', 1.5);
        legend({'$Q_s$','$Q_{rad,2}$','$Q_{rad,3}$','$Q_{rad,4}$','$Q_{rad}$'}, ...
            'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
    else
        plot(profil.xli, Qrad(2,:), '-', 'LineWidth', 1.5, 'Marker', 's');
        legend({'$Q_s$','$Q_{rad}$'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
    end
    xlabel('r/a', 'FontSize', 14);
    ylabel('Q [W/m^3]', 'FontSize', 14);
    set(gca, 'YScale', 'log', 'FontSize', 12);



    % === Figure D_eff ===
    D_eff = - G_loop(:,:,1) ./ nd1_loop(:,:,1) .* rhomax;
    if ~exist('hf_Deff', 'var') || ~isgraphics(hf_Deff)
        hf_Deff = figure('tag','Fig_Deff', 'name','D_eff Profile', 'color',[1 1 1]);
    else
        figure(hf_Deff);
    end
    clf(hf_Deff);
    hold on; box on; grid on;

    plot(profil.xli, D_eff(2,:,1), '-', 'LineWidth', 1.5, 'Marker', 's');
    xlabel('r/a', 'FontSize', 14);
    ylabel('D (m^2/s)', 'FontSize', 14, 'Interpreter', 'latex');
    legend({'$D_{eff}$'}, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);


    varargout = {chie_loop, chii_loop, chii_neo, DD_loop, VV_loop, G_loop};

    if qlkparams.imp_trans >= 3
        DD_loop_n(:,:,1) = 0; DD_loop_n(:,:,end) = 0;
        VV_loop_n(:,:,1) = 0; VV_loop_n(:,:,end) = 0;
        G_loop_n(:,:,1)  = 0; G_loop_n(:,:,end)  = 0;
        varargout(7:9) = {DD_loop_n, VV_loop_n, G_loop_n};
    end

elseif qlkparams.imp_trans == 1

    % === FIGURE : -V/D (raw) ===
    figure; hold on; box on; grid on;
    colors = lines(6); labels = {'$1$', '$2$', '$W$', '$H$'};
    for i = 1:4
        idx = i + 1;
        plot(profil.xli, -VV_raw(2,:,idx)./DD_raw(2,:,idx), '-', ...
            'Color', colors(idx,:), 'LineWidth', 1.5, 'Marker', 'o');
    end
    xlabel('r/a'); ylabel('-V/D');
    legend(labels, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);

    varargout = {chie_loop, chii_loop, chii_neo, DD_raw, VV_raw, G_raw};

else
    varargout = {chie_loop, chii_loop, chii_neo, D_loop, V_loop};
end

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
    QLK_nn_data.NE   = n_loop(:,:,end); 
    QLK_nn_data.N    = n_loop;
    QLK_nn_data.err  = Inf;
else
    if evalin('base','exist(''QLK_nn_data'')')
        err_press_mem = evalin('base','QLK_nn_data.err');
    else
        err_press_mem = Inf;
    end
    QLK_nn_data.time = profil.temps;
    QLK_nn_data.x    = profil.xli;
    % QLK_nn_data.err  = sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2));
    % [err_press,ind_best] = min(sqrt(sum(mask_err .* ((profil.tep .* profil.nep + profil.tip .* profil.nep) - (te_loop .* ne_loop + ti_loop .* ne_loop)) .^ 2,2)));

    % QLK_nn_data.err  = sqrt(sum(mask_err .* ((profil.tep .* n_init(:,:,end) + profil.tip .* n_init(:,:,1)) - (te_loop .* n_loop(:,:,end) + ti_loop .* n_init(:,:,1))) .^ 2,2));
    % [err_press,ind_best] = min(sqrt(sum(mask_err .* ((profil.tep .* n_init(:,:,end) + profil.tip .* n_init(:,:,1)) - (te_loop .* n_loop(:,:,end) + ti_loop .* n_init(:,:,1))) .^ 2,2)));

    ion_sum_init = zeros(size(n_init(:,:,1)));
    ion_sum_loop = zeros(size(n_init(:,:,1)));
    for i = 1:5
        ion_sum_init = ion_sum_init + profil.tip .* n_init(:,:,i);
        ion_sum_loop = ion_sum_loop + ti_loop .* n_loop(:,:,i);
    end
    % électrons (déjà en "end")
    elec_init = profil.tep .* n_init(:,:,end);
    elec_loop = te_loop .* n_loop(:,:,end);
    % calcul de l’erreur
    QLK_nn_data.err = sqrt(sum(mask_err .* ((elec_init + ion_sum_init) - (elec_loop + ion_sum_loop)) .^ 2, 2));
    [err_press, ind_best] = min(QLK_nn_data.err);
    QLK_nn_data.err  = err_press;
    vt = ones(size(profil.temps));
    %disp([err_press ,err_press_mem])
    if err_press > err_press_mem
        te_exp = 0.1 .* te_loop(ind_best,:) .* n_loop(ind_best,:,end) + 0.9 .* profil.tep(ind_best,:) .* n_init(ind_best,:,end);
        ti_exp = 0.1 .* ti_loop(ind_best,:) .* n_loop(ind_best,:,1) + 0.9 .* profil.tip(ind_best,:) .* n_init(ind_best,:,1);
        ne_exp = 0.1 .* n_loop(ind_best,:,end) + 0.9 .* n_init(ind_best,:,end);
        n1_exp = 0.1 .* n_loop(ind_best,:,2) + 0.9 .* n_init(ind_best,:,2);
        nW_exp = 0.1 .* n_loop(ind_best,:,4) + 0.9 .* n_init(ind_best,:,4);
    else
        te_exp = 0.3 .* te_loop(ind_best,:) .* n_loop(ind_best,:,end) + 0.7 .* profil.tep(ind_best,:) .* n_init(ind_best,:,end);
        ti_exp = 0.3 .* ti_loop(ind_best,:) .* n_loop(ind_best,:,1) + 0.7 .* profil.tip(ind_best,:) .* n_init(ind_best,:,1);
        ne_exp = 0.3 .* n_loop(ind_best,:,end) + 0.7 .* n_init(ind_best,:,end);
        n1_exp = 0.3 .* n_loop(ind_best,:,2) + 0.7 .* n_init(ind_best,:,2);
        nW_exp = 0.3 .* n_loop(ind_best,:,4) + 0.7 .* n_init(ind_best,:,4);
    end
    QLK_nn_data.TE   = vt * (te_exp ./ ne_exp);
    QLK_nn_data.TI   = vt * (ti_exp ./ ne_exp);
    QLK_nn_data.NE   = vt * ne_exp;
    QLK_nn_data.N(:,:,2) = vt * n1_exp;
    QLK_nn_data.N(:,:,4) = vt * nW_exp;
    QLK_nn_data.N_bar = vt * trapz(profil.xli,ne_exp);
end
zassignin('base','QLK_nn_data',QLK_nn_data);

hz =findobj(0,'type','figure','tag','Prof_net_kin');
if isempty(hz)
    hz=figure('tag','Prof_net_kin','name',sprintf('Profiles prediction %s  (fixed sources but equipartition, no time dependance)',qlkparams.channel_tag));
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
if first_call == 0
   zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.TE ./ 1e3,'color','k','linestyle',':'); 
end
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
if first_call == 0
    zplotprof(gca,profil.temps,profil.xli, QLK_nn_data.TI ./ 1e3,'color','k','linestyle',':');
end
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
if first_call == 0
   zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.NE ./ 1e19,'color','k','linestyle',':');
end
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
if first_call == 0
    zplotprof(gca,profil.temps,profil.xli,NaN * QLK_nn_data.NE ./ 1e19,'color','k','linestyle',':');
%    legend('METIS',qlkparams.channel_tag,'interval of computation with QLKNN-4Dkin','External data');
else
%    legend('METIS',qlkparams.channel_tag,'interval of computation with QLKNN-4Dkin');
end
set(gca,'YScale','linear');
ylabel('Q_{e,i} (MW)');
xlabel('r/a');
%z0loglin(gca);
zoom yon


% this plot can be improved
if first_call == 1
    hz =findobj(0,'type','figure','tag','Wth_net_kin');
    if isempty(hz)
        hz=figure('tag','Wth_net_kin','name',sprintf('Energy content prediction %s (fixed sources but equipartition, no time dependance)',qlkparams.channel_tag));
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
    legend('METIS',qlkparams.channel_tag);
    zoom on
    subplot(3,1,2);
    semilogy(profil.temps,wth_nn ./ max(1,wth),'b',profil.temps,ones(size(profil.temps)),'g');
    %xlabel('time (s)');
    legend(sprintf('W_{th,%s} / W_{th,METIS}',qlkparams.channel_tag));
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

mask_plot = mask;
hz =findobj(0,'type','figure','tag','flux_net_kin');
if isempty(hz)
    hz=figure('tag','flux_net_kin','name',sprintf('Flux comparison between METIS & %s (fixed sources but equipartition, no time dependance)',qlkparams.channel_tag));
else
    figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
    'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)
zplotprof(gca,profil.temps,profil.xli,profil.qe  + profil.qei - 5/2 .* phys.e .* te_loop .* ge_sts + mask_plot,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_qe + Qei + mask_plot,'color','r');
set(gca,'YScale','linear');
ylabel('flux_qe (W/m^{-3})');
%      legend('METIS','QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
title('Without Qei flux and particles convected flux');
zoom yon

subplot(2,2,2)
zplotprof(gca,profil.temps,profil.xli,profil.qi   - profil.qei - 5/2 .* phys.e .* ti_loop .* gi_sts + mask_plot,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_qi - Qei + mask_plot,'color','r');
set(gca,'YScale','linear');
ylabel('flux_qi (W/m^{-3})');
%legend('METIS','QLKNN-4Dkin');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,ge_sts + mask_plot,'color','b');
zplotprof(gca,profil.temps,profil.xli,flux_ne + mask_plot,'color','r');
set(gca,'YScale','linear');
ylabel('flux_ne (e^-/m^{-3}/s)');
legend('METIS (input)',sprintf('%s (ouput/balance)',qlkparams.channel_tag));
xlabel('r/a');
%z0loglin(gca);
zoom yon

hz =findobj(0,'type','figure','tag','Chi_net_kin');
if isempty(hz)
    hz=figure('tag','Chi_net_kin','name',sprintf('Transport coefficients comparison between METIS & %s (fixed sources but equipartition, no time dependance)',qlkparams.channel_tag));
else
    figure(hz);
end
clf
set(hz,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
    'defaultlinelinewidth',3,'color',[1 1 1])

subplot(2,2,1)
zplotprof(gca,profil.temps,profil.xli,profil.xie .* profil.grho2 ./ profil.grho,'color','b');
zplotprof(gca,profil.temps,profil.xli,chie_loop,'color','r');
zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.Chi_e_QLK,'color','m'); 
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
zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.Chi_i_QLK,'color','m'); 
set(gca,'YScale','linear');
ylabel('Chi_i (m^2/s)');
%legend('METIS','QLKNN-4Dkin prediction','Neo');
xlabel('r/a');
%z0loglin(gca);
zoom yon

subplot(2,2,3)
zplotprof(gca,profil.temps,profil.xli,profil.dn,'color','b');
zplotprof(gca,profil.temps,profil.xli,D_loop,'color','r');
zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.D_ne_QLK,'color','m'); 
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
zplotprof(gca,profil.temps,profil.xli,QLK_nn_data.V_ne_QLK ./ QLK_nn_data.D_ne_QLK .* profil.grho2 ./ profil.grho,'color','m'); 
set(gca,'YScale','linear');
ylabel('V / D (m)');
%legend('METIS','QLKNN-4Dkin prediction');
legend('METIS',sprintf('%s prediction',qlkparams.channel_tag),'Neoclassic (Hinton)',sprintf('Raw %s',qlkparams.channel_tag));
xlabel('r/a');
%z0loglin(gca);
zoom yon

if ishandle(hdisc)
    delete(hdisc);
end

% backup internal data for future use
tempf = tempname;
clear hz ans hdisc hwaitbar
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
    case 5
        alpha  = nI_ZI2 ./ profil.n1p ./ (1 + 4 .* (profil.nhep + profil.nhe3p) ./ max(1,profil.nDp));
    case 11
        iso    = profil.nBp ./ max(1,profil.n1p);
        nI_ZI2 = nI_ZI2 + profil.n1p .* (1 + iso);
        alpha  = nI_ZI2 ./ profil.n1p ./ (1 + 25 .* iso);        
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

