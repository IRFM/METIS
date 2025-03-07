function z0dinput = zerod_init_st40_from_astra(astra_file,any_metis_model_file,parameters_filename)

% handle of input
if nargin == 0 || isempty(astra_file)
   error('syntax: z0dinput = zerod_init_st40_from_astra(astra_file,any_metis_model_file)'); 
end
if nargin < 2 || isempty(any_metis_model_file)
    % load model
    any_metis_model_file = fullfile(fileparts(which('metis')),'certification','metis','MAST_like.mat');
    metis_model = load(any_metis_model_file,'post');
    z0dinput = metis_model.post.z0dinput;
elseif ischar(any_metis_model_file)
    metis_model = load(any_metis_model_file,'post');
    z0dinput = metis_model.post.z0dinput;
elseif isstruct(any_metis_model_file)
    if isfield(any_metis_model_file,'post')
        z0dinput = any_metis_model_file.post.z0dinput;
    elseif isfield(any_metis_model_file,'z0dinput')
        z0dinput = any_metis_model_file.z0dinput;
    else 
        z0dinput = any_metis_model_file;
    end
else
    error('unable to use input "any_metis_model_file" as METIS data model');
end
if nargin < 3
    parameters_filename = [];
end

% load Astra data
astra_data = [];
load(astra_file)

% clean exp0d
if isfield(z0dinput,'exp0d')
    noms = fieldnames(z0dinput.exp0d);
    for k =1:length(noms)
        z0dinput.exp0d.(noms{k}) = [];
    end
end
%
langue                 =  'anglais';
% default 
% list of option change
option = z0dinput.option;   % use by default loaded option 
% structure de donnees vide pour METIS
if ischar(parameters_filename)
    if ~isempty(parameters_filename)
        option = z0doverwriteparam(parameters_filename,option);
    end
elseif isstruct(parameters_filename)
    option = parameters_filename;
else
    % ST40 tuning
    option.gaz          = 1;
    option.zmax         = 8;
    option.zimp         = 6;
    option.rimp         = 0.1;    % must be provided by Sycomore
    %option.frhe0        = 0;
    %option.tauhemul     = -3; % default value in Sycomore is 5
    option.neasser      = 1;
    option.Recycling    = 0.95;
    option.natural      = 0;
    option.fnbar_nat    = 1;
    %option.nea_factor   = 1;   % Sycomore model ?
    %option.ftaup        = 1;
    option.fn0a             = 1;
    option.fn0a_div         = 0.1;
    option.ane              = 11;
    %option.vane             = 1;
    %option.ne_shape         = 'Auto';
    option.ne_free          = 3;       % better desciption of density profile
    %option.neped_expo       = -0.7000; % not wellknown from present experiment
    %option.pix              = 0.7000;
    %option.piw              = 0;
    %option.pif              = 0;
    option.scaling          = 0;         % 12 is better but not yet implemanted in Sycomore (= 0)
    option.dilution         = 0;          % effect not taken into account in Sycomore
    option.tau_limitation   =  'Off';      % effect not taken into account in Sycomore but usefull for rampup and rampdown
    option.ploss_exp        = 'max_power'; % best choice from DEMO1 and DEMO2 studies for Ploss computation
    %option.fprad            = 0.3333;
    %option.HH_li            = 0;
    %option.l2hscaling       = 0;          % not optimal, but must be the same as in Sycomore
    option.l2hscaling        = 2;          %loi ITER LH02Zeff (Takizuka,PPCF, 2004)
    option.modeh            = 2;
    option.l2hmul           = 1;
    option.plhthr           = 'P_LCFS';    % used the real power crossing the LCFS
    %option.l2hslope         = 0;
    option.hysteresis       = 0;           % No hysteresis for the back transition to  L-mode
    option.fpped            = 1;
    option.hmore_pped       = 2;
    %option.fstiff           = 1;
    option.usepped_scl      = 2;         % use minmum pressure between scaling for predestal pressure and standard METIS rule Pped = K (W_H -W_L)
    %option.taurotmul        = 0;
    option.fintrinsic       = 0;         % machine dependent factor: intrinsic rotation must be tuned (used scaling in collisionality)
    %option.omega_shape      = 0;
    option.xiioxie          = -4.5;      % which are the right value: consistant computation (-4.5) or equal value ?
    %option.kishape          = 3;
    %option.ki_expo          = 2;
    %option.xieorkie         = 0;
    option.grad_ped         = 3;         % allows radiative collapse
    option.qdds             = 1;         % continuous sawtooth model on q =1 (initite frequency)
    %option.w1               = 0.5;
    %option.epsq             = 0.01;
    %option.ddsmode          = 0;
    option.kidds            = 3;         % increase factor of the transport in inversion radius
    %option.peeling          = 0;
    %option.dwow_elm         = 0;
    %option.tau_elm_factor   = 10;
    option.runaway          = 0;        % no ruaway for the first studies
    option.modeboot         = 1;        % default Sauter-Angioni formalation
    %option.bootmul          = 1;
    %option.ffit_ped         = 1;
    %option.fspot            = 0.15;
    %option.vloop            = 0;       % default option: Ip given
    %option.vref              = 0.1;
    %option.tswitch           = 1e6;
    %option.li               = 1;
    %option.breakdown        = -11.6867;
    %option.berror           = 0;        % for the first studies breadown phas eis no included
    %option.L_eddy           = 0;
    %option.R_eddy           = 0;
    %option.p_prefill        = 1.0000e-03;
    option.zeff             = 0;         % Zeff without He ashes is provided
    option.faccu            = 0;         % No impurities accumulation
    %option.heat_acc         = 0;
    %option.fne_acc          = 0;
    option.W_effect         = 0;          % No tunsten 
    option.density_model    = 'minconv';  % to have a feeling of transport coefficients need to obtain the design density profile
    option.frad             = 1;          % we trust the METIS radiative model
    option.matthews         = 1;          % line radiation computed with Matthews
    %option.z_prad          = 'zmax'
    option.gaunt            = 1;          % used tabulated Gaunt factor instead of 1.2 (higher accuracy)
    %option.noncoronal       = 1;          % non coronal correction for ramp up
    option.noncoronal       = 0;          % non coronal don't work for the moment)
    option.rw               = 0.7000;     % faction of cyclotron radiation reflected (must be provided by Sycomore)
    option.configuration    = 2;          % always in divertor mode
    option.lambda_scale     = 3;          % SOL width scaling
    option.factor_scale     = 1;
    %option.sol_lscale       = 0;
    option.eioniz           = 0;          % used tabulated law
    %option.de               = 0.5000;
    %option.alpha_e          = 0.8200;
    %option.fnesol           = 0;
    %option.sol_model        = '2_points'; % 2 points model actived
    option.sol_model        = 'scaling';  % for testing
    option.sol_rad          = 'decoupled';% convinient value for sol_rad
    option.lcx              = 1;          % connexion length
    option.fpower           = 0.6000;     % fraction of powe ron outer target
    option.fR_target        = 1;          % position of outer target if known
    option.fcond            = -1;         % with kinetic correction
    option.fmom             = 0;          % compute friction term on neutral
    option.mach_corr        = 1;          % allows Mach number > 1 at the target
    option.yield_model      = 'Javev';    % default Sputtering model
    option.ftwleak          = -1;         % DiVImp fit for leakage
    option.cw_factor        = 0;          % desactivate W source from target, reserved for future application with Sycomore
    option.cw_offset        = 0;       % tungsten concentration in core plasma   , must be provide by Sycomore
    option.fzmax_div        = 1;          % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
    %option.angle_ece        = 180;        % optimisation of current drive by using HFS resonnance
    %option.synergie          = 0;
    %option.sens              = 1;         % co current ECCD
    %option.eccdmul           = 1;
    option.angle_nbi         = 90; % co current
    option.rtang             = 0.42;
    option.zext              = 0;
    option.einj              = 25e3 * (0.7 + 0.5 * 0.1 + 0.33 * 0.2);    % energy of injected neutral, first ijector (eV)
    option.nbicdmul          = 1;
    %option.drs1 = 0         % tune later
    %option.dzs1 = 0
    option.e_shielding       = 'Honda-Sauter'; % backcurrent formula including collisionality effect.
    %option.cur_nbi_time      = 0;
    option.nb_nbi            = 2;          % two NBI modules
    option.angle_nbi2        = 90; % co current
    option.rtang2            = 0.42;
    option.zext2             = 0;
    option.einj2             = 55e3 * (0.65 + 0.5 * 0.25 + 0.33 * 0.1);    % energy of injected neutral, second ijector (eV)
    option.nbicdmul2         = 1;
    %option.drs2 = 0             % tune later
    %option.dzs2 = 0
    option.fast_ion_sbp      = 3;
    %option.lhmode            = 5;          % LHCD channel used for 2nd ECCD system
    %option.upshiftmode       = 'newmodel';
    %option.fupshift          = 1;
    %option.etalh             = 1;          % 2nd eccd system in co-current
    %option.npar0             = 2;
    %option.freqlh            = 3.7000;
    %option.wlh               = 0;
    %option.xlh               = 0;          % 2nd ECCD system position of maximum heat depostion
    %option.dlh               = 0.4000;     % 2nd ECCD system width of heat depostion profile
    %option.npar_neg          = 0;
    %option.angle_ece2        = 180;        % optimisation of current drive by using HFS resonnance
    %option.fwcd              = 0;
    %option.mino              = 'T';
    %option.cmin              = 1;
    %option.nphi              = 25;
    %option.freq              = 55.0357;
    %option.icrh_width        = 1;
    %option.fact_mino         = 0;
    option.sitb              = 0;          % no ITB, reserved for a future version
    %option.itb_sensitivity   = 1;
    %option.itb_slope_max     = 2;
    option.tae               = 0;          % no TAe fast alpha profile spreading
    option.alpha_channeling  = 0;          % no alpha channeling, needs external infomration to be tuned
    option.smhd              = 100;        % no MHD limit
    option.tmhd              = 0;          % no MHD limi
    %option.rip               = 0;
    option.signe             = 1;          % sign of Bt . J_phi: must be provided by the machine description
    option.laochange         = 1;          % always take into account toricity
    %option.morphing          = 5;
    option.mode_expo_inte    = 1;          % large time step integration
    option.cronos_regul      = 4;          % better for plot
    %option.impur_rot         = 'imp';
    %option.carnot            = 0.4200;    % better computation in Sycomore
    %option.mul_blanket       = 1.2000;    % better computation in Sycomore
    %option.aux               = 0.0500;    % better computation in Sycomore
    %option.effinj            = 0.7000;    % better computation in Sycomore
    option.available_flux    = 7;        % to be defined
    %option.machine           = 'ITER@6.2_m&15_MA'; % dynamically generated
    %option.shot              = 1;         % provided by Kepler
    %option.reference_parameters= ''       % not used
    %option.dwdt_method       = 'implicit';
    %option.nbmax             = 31;
    %option.tol0d= 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UAL writing control parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    option.init_output_cpo = 0;
    option.restart= '';
    option.coreprof= 1;
    option.coretransp= 1;
    option.coresource_lhcd= 0;
    option.coresource_eccd= 0;
    option.coresource_icrh= 0;
    option.coresource_nbicd= 1;
    option.coresource_fusion= 1;
    option.coreneutrals= 1;
    option.coresource_radiation= 1;
    option.coresource_cyclotron= 1;
    option.neoclassic= 1;
    option.coresource_full= 1;
    option.equilibrium= 1;
    option.grid_equi= 0;
    option.scenario_occurrence= '';
    option.coreprof_occurrence= '';
    option.coretransp_occurrence= '';
    option.coreneutrals_occurrence= '';
    option.neoclassic_occurrence= '';
    option.equilibrium_occurrence= '';
    option.coresources_occurrence= '';
    
    % default value for ITM
    option.COCOS  = 11;
end
z0dinput.option = option;
% parameter changing for shot to shot
z0dinput.cons.temps    = double(astra_data.TIME.data(:));
vt = ones(size(z0dinput.cons.temps));
signe_ip = sign(mean(double(astra_data.GLOBAL.IPL.data)));
z0dinput.cons.ip = abs(double(astra_data.GLOBAL.IPL.data(:)));
% default values
z0dinput.cons.flux    = -double(astra_data.GLOBAL.FBND.data(:))/2/pi;
z0dinput.cons.nbar    = mean(double(astra_data.PROFILES.ASTRA.NE.data),1)' * 1e19;
z0dinput.cons.picrh   = 0 .* vt;
z0dinput.cons.plh     = 0 .* vt;
% maybe theorder is the right one
z0dinput.cons.pnbi    = 0 * vt; % not availale in ASTRA file;
% set all in hydrogen
z0dinput.cons.ftnbi    = vt * (1 + sqrt(-1));
z0dinput.cons.iso    = vt;
%
z0dinput.cons.pecrh   = 0 .* vt;
z0dinput.cons.zeff    = double(astra_data.GLOBAL.ZEFF.data(:));
z0dinput.cons.xece    = 0 .* vt;
z0dinput.cons.hmore   = double(astra_data.GLOBAL.H98.data(:));
    
% plasma shape
R_LCFS = double(astra_data.P_BOUNDARY.RBND.data');
Z_LCFS = double(astra_data.P_BOUNDARY.ZBND.data');

z0dinput.geo.a        = (max(R_LCFS,[],2) - min(R_LCFS,[],2)) / 2;
z0dinput.geo.R        = (max(R_LCFS,[],2) + min(R_LCFS,[],2)) / 2;
z0dinput.geo.K        = (max(Z_LCFS,[],2) - min(Z_LCFS,[],2)) ./ (max(R_LCFS,[],2) - min(R_LCFS,[],2));
z0dinput.geo.z0       = (max(Z_LCFS,[],2) + min(Z_LCFS,[],2)) / 2;
z0dinput.geo.d        = 0      .* vt;
for k=1:length(z0dinput.geo.z0)
    indzmax = find(Z_LCFS(k,:) ==  max(Z_LCFS(k,:)),1);
    dup     =  (R_LCFS(k,indzmax) - z0dinput.geo.R(k)) / z0dinput.geo.a(k); 
    indzmin = find(Z_LCFS(k,:) ==  min(Z_LCFS(k,:)),1);
    dlow    =  (R_LCFS(k,indzmin) - z0dinput.geo.R(k)) / z0dinput.geo.a(k); 
    z0dinput.geo.d(k)        =  (dup + dlow) / 2;
end
z0dinput.geo.b0       = abs(double(astra_data.GLOBAL.BTVAC.data(:)) * 0.5 ./ z0dinput.geo.R);
z0dinput.option.orientation  = sign(double(mean(astra_data.GLOBAL.BTVAC.data))); 
z0dinput.option.signe = signe_ip * z0dinput.option.orientation ;
z0dinput.geo.vp       = [];
z0dinput.geo.sp       = [];
z0dinput.geo.sext     = [];

z0dinput = z0dsepageo(z0dinput, R_LCFS, Z_LCFS, z0dinput.cons.temps);

% relative to aminor
z0dinput.option.drs1 = sqrt(pi * (0.25/2) ^ 2) / 2 /  mean(z0dinput.geo.a);
z0dinput.option.dzs1 = sqrt(pi * (0.25/2) ^ 2) / 2 /  mean(z0dinput.geo.a);
z0dinput.option.drs2 = sqrt(pi * (0.15/2) ^ 2) / 2 /  mean(z0dinput.geo.a);
z0dinput.option.dzs2 = sqrt(pi * (0.15/2) ^ 2) / 2 /  mean(z0dinput.geo.a);

% experimental data
z0dinput.exp0d.temps  = z0dinput.cons.temps;
z0dinput.exp0d.li     = double(astra_data.GLOBAL.LI3.data(:));
z0dinput.option.li    = z0dinput.exp0d.li(1);
z0dinput.exp0d.wth    = double(astra_data.GLOBAL.WTH.data(:));
z0dinput.exp0d.ploss  = double(astra_data.GLOBAL.P_TOT_E.data(:))*1e6 + double(astra_data.GLOBAL.P_TOT_I.data(:))*1e6;
z0dinput.exp0d.vloop  = zdxdt(double(astra_data.GLOBAL.FBND.data(:)),z0dinput.cons.temps);
z0dinput.exp0d.ip     = double(astra_data.GLOBAL.IPL.data(:));
z0dinput.exp0d.nbar   = z0dinput.cons.nbar;
z0dinput.exp0d.ne0    = double(astra_data.PROFILES.ASTRA.NE.data(1,:))'* 1e19;
z0dinput.exp0d.te0    = double(astra_data.PROFILES.ASTRA.TE.data(1,:))'* 1e3;
z0dinput.exp0d.nem    = double(astra_data.GLOBAL.NEV.data(:))* 1e19;
z0dinput.exp0d.tem    = double(astra_data.GLOBAL.TEV.data(:))* 1e3;
z0dinput.exp0d.edgeflux  = - double(astra_data.GLOBAL.FBND.data(:)) / 2 /pi;
z0dinput.exp0d.betap  = double(astra_data.GLOBAL.BETP.data(:));
z0dinput.exp0d.zeff  = z0dinput.cons.zeff;
z0dinput.exp0d.iboot = double(astra_data.GLOBAL.I_BS.data(:)) * 1e6;
z0dinput.exp0d.icd   = double(astra_data.GLOBAL.I_NBI.data(:)) * 1e6 + double(astra_data.GLOBAL.I_RF.data(:)) * 1e6;
z0dinput.exp0d.inbicd   = double(astra_data.GLOBAL.I_NBI.data(:)) * 1e6 ;
z0dinput.exp0d.ini   = z0dinput.exp0d.icd + z0dinput.exp0d.iboot;
z0dinput.exp0d.iohm   = double(astra_data.GLOBAL.I_OH.data(:)) * 1e6 ;   
z0dinput.exp0d.nebord = double(astra_data.PROFILES.ASTRA.NE.data(end,:))' * 1e19;
z0dinput.exp0d.ni0    = double(astra_data.PROFILES.ASTRA.NI.data(1,:))' * 1e19;          
z0dinput.exp0d.nibord = double(astra_data.PROFILES.ASTRA.NI.data(end,:))' * 1e19;      
z0dinput.exp0d.pfus   = double(astra_data.GLOBAL.P_FUS_TOT.data(:)) * 1e6;    
z0dinput.exp0d.pion_nbi = double(astra_data.GLOBAL.P_NBI_I.data(:)) * 1e6; 
z0dinput.exp0d.pnbi     = double(astra_data.GLOBAL.P_NBI_E.data(:)) * 1e6 + double(astra_data.GLOBAL.P_NBI_I.data(:)) * 1e6;        
z0dinput.exp0d.pohm     = double(astra_data.GLOBAL.P_OH.data(:)) * 1e6;      
z0dinput.exp0d.q0       = double(astra_data.PROFILES.PSI_NORM.Q.data(1,:))';        
z0dinput.exp0d.q95      = double(astra_data.GLOBAL.Q95.data(:));     
z0dinput.exp0d.qa       = double(astra_data.GLOBAL.QWL.data(:));           
z0dinput.exp0d.qmin     = min(double(astra_data.PROFILES.PSI_NORM.Q.data),[],1)';         
z0dinput.exp0d.tebord   = double(astra_data.PROFILES.ASTRA.TE.data(end,:)) * 1e3;
z0dinput.exp0d.tem      = double(astra_data.GLOBAL.TEV.data(:)) * 1e3;
z0dinput.exp0d.tibord   = double(astra_data.PROFILES.ASTRA.TI.data(end,:))' * 1e3;
z0dinput.exp0d.tite     = double(astra_data.GLOBAL.TIV.data(:)) ./double(astra_data.GLOBAL.TEV.data(:));
z0dinput.exp0d.vp       = double(astra_data.GLOBAL.VOL.data(:));
z0dinput.machine     = 'ST40';
% try to find shot number from filename
z0dinput.shot = [];
if ischar(astra_file)
    astra_file = strrep(astra_file,'_',' ');
    for k=1:30
       [rep,astra_file] = strtok(strtrim(astra_file),' '); 
       try
           num = str2num(rep);
       catch          
           num = [];
       end
       if ~isempty(num)
          z0dinput.shot = num;
          break;
       end
       if isempty(astra_file)
           break;
       end
    end
end
% make shot number from comment
if isempty(z0dinput.shot)
    [s,t] = unix(sprintf('echo "%s" | md5sum',astra_data.HELP));
    if s == 0
        t = strtok(t,' ');
        t(t<=32) = [];
        z0dinput.shot   = hex2num(t);
    else
        z0dinput.shot        = NaN;
    end
end
z0dinput.mode_exp = 201; % for ST40 nexw entry
% copy duplicated info
z0dinput.option.shot    = z0dinput.shot;
z0dinput.option.machine = z0dinput.machine;

% complete exp0d structure
noms = fieldnames(z0dinput.exp0d);
for k=1:length(noms)
    if ~isfield(z0dinput.exp0d,noms{k}) || isempty(z0dinput.exp0d.(noms{k}))
        z0dinput.exp0d.(noms{k}) = NaN * vt;
    end
end


function yi = interp10d(x,y,xi,methode)

    if (size(x,1)>1) && (size(x,2)==1)
        indnok = 0;
        nb     = 100;
        while (~isempty(indnok)) && (nb >0)
            indnok = find(diff(x)<=0);
            if ~isempty(indnok)
                x(indnok) = [];
                y(indnok,:) = [];
        
            end
            nb = nb - 1;
        end
    end
    yi = interp1(x,y,xi,methode);



function data = decimate_data(t, data, t_out)

    % fast resample data using a running average filter

%    interp1(data_out.U_loop.time_axis.data(t_resample/2:t_resample:end-t_resample/2)/1e3,compass_resample_data(data_out.U_loop.data,t_resample),shotime,'linear');%

    % enforce column vector
    % TODO adjust automatically for row/column vectors
    data = reshape(data, [], 1);
    ww = floor(min(diff(t_out)) / min(diff(t)));
    % cut the last data points
    % TODO properly treat the last data points
    n_mod = mod(length(data), ww);
    if n_mod > 0
        last_data = mean(data(end-n_mod:end));
    else
        last_data = [];
    end
    data = data(1:end-n_mod);
    % reshape the data and calculate the mean
    data = reshape(data, ww, []);
    data = mean(data, 1);


