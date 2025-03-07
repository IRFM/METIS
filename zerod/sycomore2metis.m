% SYCOMORE2METIS : prototype of scenarion generator for reactor; starting point for the coupling to Sycomore
%-----------------------------------------------------------------------------------------------------------
% function Matlab 2012b: sycomore2metis.m -> sycomore2metis
%
% This function is the starting point for METIS coupling to SYCOMORE
% syntax: 
%
%   testing :       
%       [z0dinput,post,option,sepa_option,sycomore_data] = sycomore2metis;
%
%   computation :
%       [z0dinput,post,option,sepa_option,sycomore_data] = sycomore2metis(sepa_option,parameters_filename,sycomore_data,0);
%
% input:
%
%     sepa_option: data structure for LCFS description (same as HELIOS)
%
%     parameters_filename: name of METIS parametrs file or METIS paramters data structure (see metis4itm.m and zerod_param.m)
%
%     sycomore_data: Input data coming from SYCOMORE
%
%     testing: flag(0/1/2). If =1, fast computatiion only with optimisation; if = 2, fast computation without optimisation
%        
% structure details for sycomore_data (with data example for ITER):
% 
%      sycomore_data.ip = 15e6;                % flat top plasma current (A).
%      sycomore_data.q95 = 3;                  % flat top q95.
%      sycomore_data.rb0 = 5.3 .* 6.2;         % magnetic rigidity (flat top)
%      sycomore_data.available_flux = 300;     % Poloidal available flux (Wb)
%      sycomore_data.device = 'ITER';          % device name (used to generate file name and comments)
%      sycomore_data.scaling  = 0;             % code for energy plasma content scaling law (same as METIS)
%      sycomore_data.H_H  = 1;                 % enhancement factor for energy content on flat top
%      sycomore_data.f_Greenwald = 0.85;       % Greenwald density fraction on flat top
%      sycomore_data.ne_peak = 1.2;            % electron density profile peaking factor; if isempty or not define, used METIS model
%      sycomore_data.shot = 1;                 % shot number, will be used to write data in UAL
%      sycomore_data.run  = 1;                 % run  number, will be used to write data in UAL
%      sycomore_data.f_ni = 0;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
%      sycomore_data.rimp = 0.1;               % ratio  between berillium and argon in core plasma
%      sycomore_data.Recycling = 0.95;         % recycling  @ divertor 
%      sycomore_data.zeff      = 1.3;          % line averaged Zeff without He ashes contribution
%      sycomore_data.tau_He_o_tau_E_core = 3;  % ratio between core confinement time of He over energy confinement time.
%      sycomore_data.rw = 0.7;                 % cyclotron radiation reflection coefficient
%      sycomore_data.PNBI = 53e6;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
%      sycomore_data.PECRH = 0;                % EC power.
%      sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
%      sycomore_data.nbi_e_inj = 1e6;          % neutral beam energy injection.
%      sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
%      sycomore_data.flux_expansion = 3;       % divertor flux expansion 
%      sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
%      sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
%      sycomore_data.cW               = 1e-5;  % tungsten concentration in core plasma   , must be provided by Sycomore 
%      sycomore_data.fzmax_div        = 0;     % Argon concentration enhancement in divertor (%)  , must be provided by Sycomore
%      sycomore_data.f_DSOL             = 1;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
%      sycomore_data.f_ne_LCFS          = 1;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
%      sycomore_data.couplage         = 0.15;  % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
%      sycomore_data.duration         = 7200;  % estimated duration of the shot in s
%      sycomore_data.shine_through_limit = 1e6;  % limit of power lost in shine throught (and first orbit losses, W)
%
% optionnal fields 
%
%	sycomore_data.tokamak  = '';           % new UAL tokamak database name
%	sycomore_data.user     = '';           % new UAL user database name
%       sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
%	sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])
%
%       sycomore_data.path4metis_files = '';   % path for files save as a postprocessing of METIS (METIS data, figures, ...)  
%                                              % kept empty to remove files creation
%
%
% structure details for sepa_option (with data example for ITER):
%
%      sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
%      sepa_option.zxup      = 1.687;     % upper altitude X point (minor radius unit)
%      sepa_option.apup      = 0;         % upper separatrix angle (R,X)  (LFS, degrees)
%      sepa_option.amup      = 0;         % upper separatrix angle (-R,X) (HFS, degrees)
%      sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
%      sepa_option.za        = 0;         % altitude of the magnetic axis (m) [0.9]
%      sepa_option.a         = 2;         % minor radius (m) [2]
%      sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
%      sepa_option.zxdo      = 2.001;     % lower altitude X point (minor radius unit)
%      sepa_option.apdo      = 22.46;     % lower separatrix angle (R,X)  (LFS, degrees)
%      sepa_option.amdo      = 67.92;     % lower separatrix angle (-R,X)  (HFS, degrees)
%      sepa_option.b0        = 11.1;      % magnetic field at R0
%      sepa_option.delta     = 1.23;      % magnetic field at R0     
%
% output:
%
%      z0dinput:      standard metis input data structure (see METIS documentation), contains data used for the computation or sample data structure in testing mode
%      post:          standard metis output data structure (see METIS documentation) or sample data structure in testing mode
%      option:        METIS parameters as used during computation (see metis4itm.m and zerod_param.m) or sample data structure in testing mode
%      sepa_option:   LCFS desciption as used during computation or sample data structure in testing mode
%      sycomore_data: same as in input or sample data structure in testing mode
%      tsnapshot:     time slice selected for performances evaluation
%
% function writed by J-F Artaud
% CVS version (created 2014/07/15)
%-----------------------------------------------------------------------
%
function [z0dinput,post,option,sepa_option,sycomore_data,tsnapshot] = sycomore2metis(sepa_option,parameters_filename,sycomore_data,testing)

% gestion of input number and contents
if nargin < 4
  testing = 0;
end
if nargin < 3
     sycomore_data = []; 
end
if nargin < 2 
    parameters_filename = [];
end
if nargin < 1
      sepa_option = [];
end 

% structure de donnees vide pour METIS
if ischar(parameters_filename)
    z0dinput = zerod_init(-2);
    if ~isempty(parameters_filename)
	option = z0doverwriteparam(parameters_filename,z0dinput.option);
    else
        option = z0dinput.option;
    end
elseif isstruct(parameters_filename)
    option = parameters_filename;
else
    % default parametes (to obtain a up to date set).
    info                = metis4itm;
    option              = info.valeur;
    % ITER tuning
    option.gaz          = 3;
    option.zmax         = 4;
    option.zimp         = 18;
    option.rimp         = 0.001;    % must be provided by Sycomore    
    %option.frhe0        = 0;
    option.tauhemul     = -3; % default value in Sycomore is 5
    option.neasser      = 1; 
    option.Recycling    = 0.95;
    option.natural      = 1;
    option.fnbar_nat    = 1; 
    %option.nea_factor   = 1;   % Sycomore model ?
    %option.ftaup        = 1;
    option.fn0a             = 1;
    option.fn0a_div         = 0.1;
    option.ane              = 11;
    %option.vane             = 1;
    %option.ne_shape         = 'Auto';
    option.ne_free          = 3;       % better desciption of density profile  
    option.neped_expo       = -0.7000; % not wellknown from present experiment
    %option.pix              = 0.7000;
    %option.piw              = 0;     
    %option.pif              = 0;     
    option.scaling          = 12;         % 12 is better but not yet implemanted in Sycomore (= 0)      
    option.dilution         = 1;          % effect not taken into account in Sycomore      
    option.tau_limitation   =  'On';      % effect not taken into account in Sycomore but usefull for rampup and rampdown
    option.ploss_exp        = 'max_power'; % best choice from DEMO1 and DEMO2 studies for Ploss computation
    %option.fprad            = 0.3333;     
    %option.HH_li            = 0;          
    %option.l2hscaling       = 0;          % not optimal, but must be the same as in Sycomore          
    option.l2hscaling        = 2;          %loi ITER LH02Zeff (Takizuka,PPCF, 2004)
    %option.modeh            = 1;              
    %option.l2hmul           = 0;          
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
    option.vref              = 0.1;     
    option.tswitch           = 1e6; 
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
    option.W_effect         = 1;          % tunsten coming from divertor is taken  into account         
    option.density_model    = 'minconv';  % to have a feeling of transport coefficients need to obtain the design density profile  
    option.frad             = 1;          % we trust the METIS radiative model         
    option.matthews         = 0;          % line radiation computed using cooling rate          
    %option.z_prad= 'zmax'       
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
    option.cw_offset        = 1e-5;       % tungsten concentration in core plasma   , must be provide by Sycomore 
    option.fzmax_div        = 1;          % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore      
    option.angle_ece        = 180;        % optimisation of current drive by using HFS resonnance
   %option.synergie          = 0;          
    option.sens              = 1;         % co current ECCD      
   %option.eccdmul           = 1;          
   %option.angle_nbi         = 90;         
   %option.rtang             = 5.2952;    
   %option.zext              = 0;          
   option.einj              = 1000000;    % energy of injected neutral, first ijector (eV)     
   %option.nbicdmul          = 1;          
   option.e_shielding       = 'Honda-Sauter'; % backcurrent formula including collisionality effect.  
   %option.cur_nbi_time      = 0;          
   %option.nb_nbi            = 2;          % two NBI modules          
   %option.angle_nbi2        = 90;
   %option.rtang2            = 5.2952;
   %option.zext2             = 0;
   option.einj2             = 1000000;    % energy of injected neutral, second ijector (eV) 
   %option.nbicdmul2         = 1;
   option.lhmode            = 5;          % LHCD channel used for 2nd ECCD system
   %option.upshiftmode       = 'newmodel';
   %option.fupshift          = 1;
   option.etalh             = 1;          % 2nd eccd system in co-current
   %option.npar0             = 2;
   %option.freqlh            = 3.7000;
   %option.wlh               = 0;
   option.xlh               = 0;          % 2nd ECCD system position of maximum heat depostion
   option.dlh               = 0.4000;     % 2nd ECCD system width of heat depostion profile
   %option.npar_neg          = 0;
   option.angle_ece2        = 180;        % optimisation of current drive by using HFS resonnance
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
   option.cronos_regul      = 2;          % better for plot
   %option.impur_rot         = 'imp';
   %option.carnot            = 0.4200;    % better computation in Sycomore
   %option.mul_blanket       = 1.2000;    % better computation in Sycomore
   %option.aux               = 0.0500;    % better computation in Sycomore
   %option.effinj            = 0.7000;    % better computation in Sycomore
   option.available_flux    = 300;        % must be provided by Sycomore
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
    option.coresource_lhcd= 1;
    option.coresource_eccd= 1;
    option.coresource_icrh= 1;
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
    option.COCOS  = 13;


end

if isempty(sycomore_data)
    sycomore_data.ip = 15e6;                % flat top plasma current (A).
    sycomore_data.q95 = 3;                  % flat top q95.
    sycomore_data.rb0 = 5.3 .* 6.2;         % magnetic rigidity (flat top)
    sycomore_data.available_flux = 300;     % Poloidal available flux (Wb)
    sycomore_data.device = 'ITER';          % device name (used to generate file name and comments)
    sycomore_data.scaling  = 0;             % code for energy plasma content scaling law (same as METIS)
    sycomore_data.H_H  = 1;                 % enhancement factor for energy content on flat top
    sycomore_data.f_Greenwald = 0.85;       % Greenwald density fraction on flat top
    sycomore_data.ne_peak = 1.2;            % electron density profile peaking factor; if isempty or not define, used METIS model
    sycomore_data.shot = 1;                 % shot number, will be used to write data in UAL
    sycomore_data.run  = 1;                 % run  number, will be used to write data in UAL
    sycomore_data.f_ni = 0;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
    sycomore_data.rimp = 0.1;               % ratio  between berillium and argon in core plasma
    sycomore_data.Recycling = 0.95;         % recycling  @ divertor 
    sycomore_data.zeff      = 1.3;          % line averaged Zeff without He ashes contribution
    sycomore_data.tau_He_o_tau_E_core = 3;  % ratio between core confinement time of He over energy confinement time.
    sycomore_data.rw = 0.7;                 % cyclotron radiation reflection coefficient
    sycomore_data.PNBI = 53e6;              % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
    sycomore_data.PECRH = 53e6;             % ECRH
    sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
    sycomore_data.nbi_e_inj = 1e6;          % neutral beam energy injection.
    sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
    sycomore_data.flux_expansion = 3;       % divertor flux expansion 
    sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
    sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
    sycomore_data.cW               = 1e-5;  % tungsten concentration in core plasma   , must be provide by Sycomore 
    sycomore_data.fzmax_div        = 1;     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore      
    sycomore_data.f_DSOL             = 1;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
    sycomore_data.f_ne_LCFS          = 1;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
    sycomore_data.couplage         = 0.15   % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
    sycomore_data.duration         = 1000;  % estimated duration of the shot in s
    % use constant P/R for allowed shine througth (and first orbit losses)
    % take < 5% of 33 MW for ITER @ 6.2m -> 2e5 W/m
    sycomore_data.shine_through_limit =  2e5 .* 6.2;  % limit of power lost in shine through (W)
    [p,fname] = fileparts(tempname);
    mkdir(fname)
    %
    sycomore_data.path4metis_files = fname;   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       

end

% input data from sycomore
% LCFS data
if isempty(sepa_option)
    sepa_option.rxup      = 0.466;     % upper triangularity (minor radius unit)
    sepa_option.zxup      = 1.687;     % upper altitude X point (minor radius unit)
    sepa_option.apup      = 0;         % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 0;         % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = 6.2;       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;         % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = 2;         % minor radius (m) [2]
    sepa_option.rxdo      = 0.568;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = 2.001;     % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 22.46;     % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 67.92;     % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = 11.1;      % magnetic field at R0
    sepa_option.delta     = 1.23;      % magnetic field at R0
end

% open wrinting on standard output
if ~isempty(sycomore_data.path4metis_files);
      if ~isempty(sycomore_data.device)
	  tokamak = sycomore_data.device;
      else
	  tokamak = sprintf('DEMO-R%d-a%d-RBt%d-Ip%d',ceil(sepa_option.ra*100),ceil(sepa_option.a*100),ceil(sycomore_data.rb0),ceil(sycomore_data.ip/1e6));
      end
      root_name = fullfile(sycomore_data.path4metis_files,sprintf('sycomore2metis_%s@%d_run_%d',tokamak,sycomore_data.shot,sycomore_data.run))
      diary off
      diary(sprintf('%s_verbatim.txt',root_name));
      diary on
      fprintf('Starting METIS simumlation from SYCOMORE result for %s\n',tokamak);
else
      root_name = '';
      tokamak = sycomore_data.device;
end



% other information
sepa_option.nbp       = 201;                 % number of points for the separatrix (depends on equilibrium module) [201]
sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
sepa_option.filename  = '';
sepa_option.za        = -(sepa_option.zxup - sepa_option.zxdo) .* sepa_option.a;       % altitude of the magnetic axis (m) [0.9]

% from separatrix
[R,a,z0,K,d]   = sepamoment(sepa_option);
%  % major radius (m)
%  R = sepa_option.ra;
%  % vertical shift
%  z0 = sepa_option.za;
%  % minor radius
%  a = sepa_option.a;
%  % elongation
%  K = (max(sepa.Z)  -  min(sepa.Z)) ./ (max(sepa.R)  -  min(sepa.R));
%  % triangularity
%  d = (sepa_option.rxup + sepa_option.rxdo) ./ 2;

% end rampdown vertical position 
%zlow = sepa_option.za - K .* a + a / 2;
zlow = z0 - K .* a + a / 2;

% volume ratio to scale power
rap_power  = R .* a .^ 2 .* K ./ (6.2 .* 2 .^ 2 .* 1.844);

% magnetic rigidity
rb0 = sycomore_data.rb0;
b0  = rb0 ./ R;
% available poloidal flux from CS and PF coils (Wb)
available_flux = sycomore_data.available_flux; %800;
% flattop plasma current
ip   = sycomore_data.ip; %15e6;

% maximal available power
% LHCD : used as a second ECRH system in DEMO
%plh   = 0.0e6;
% ICRH
picrh = 0e6;

% will scale with volume on ITER design
%  option for ECRH only or mixed ECRH/NBI
% the maximum power will also scale with L2H threshold 
% ECRH flattop
if isfield(sycomore_data,'PECRH') && (sycomore_data.PECRH >= 0)
    pecrh1   = sycomore_data.PECRH ./ rap_power;
else
    pecrh1   = 0;
end
% ecrh power are lineraly interpolated durant ramp-up and ramp-down
% preheat  (from start to xpoint, limiter mode)
p0_ecrh1 = 0e6;
% ramp up (from xpoint formation to H mode transition)
p1_ecrh1 = 0e6;
% end rampup (from hmode transition to end of ramp-up)
p2_ecrh1 =  max(20e6, pecrh1/2);
% full power, up to ignition
pecrh_plus1 = pecrh1;

% the second ecrh is used for pre-ionisation and rampup
 % ECRH flattop
pecrh2   = 0e6;
% ecrh power are lineraly interpolated durant ramp-up and ramp-down
% preheat  (from start to xpoint, limiter mode)
p0_ecrh2 = 2e6;
% ramp up (from xpoint formation to H mode transition)
p1_ecrh2 =  max(8e6,pecrh1/3);
% end rampup (from hmode transition to end of ramp-up)
p2_ecrh2 = 0e6;
% full power, up to ignition
pecrh_plus2 = pecrh2;

% NBI : added ICRH power from ITER
pnbi1  =  (10e6 + 16.5e6);
pnbi2  =  (10e6 + 16.5e6);
onnbi2 = 1;

% breakdown + rampup information
% li    = 0.5;           % given as input parameter             
if isfield(sycomore_data,'CEjima') && (sycomore_data.CEjima > 0)
  % for use in METIS GUI
  CEjima    = sycomore_data.CEjima;
else
  CEjima    = 0.35;         % heated help rampup
end
frampdown = 0.75;


% heat source defautl parameters
xece  = 0.7;
xece_itb  = 0;

% physics model assumption
% Greenwald  fraction 
fnegr =  sycomore_data.f_Greenwald;
% conversion (10^20 m^-3)
nbar  = fnegr .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
 % enhancement factor for energy content on flat top
hmore = sycomore_data.H_H;

% plasma composition
iso   = 1;
ftnbi = 0;
zeff  = sycomore_data.zeff;

%%%%%%%%%%%%%%%%%%%%%%%%
% start data estimation
%%%%%%%%%%%%%%%%%%%%%%%%
% scaling fir dIp/dt
% we kept E// = 0.3V/m as for  ITER breakdown,
breakdown = 0.3 .* 2 .* pi .* 6.2 .^ 2 ./ R;
% constant
mu0  = (4*pi.*1e-7);
% we assume CEjima ITER = 0.45.
CEjima_iter  =  0.45;
dt_iter   = 80;
ip_iter   = 15e6;
% ip ramp up rate for ITER
dipdt_iter = ip_iter/dt_iter;
% ITER reference loop voltage
Vloop   = (CEjima_iter + (log(8.*6.2 ./ 2 ./ sqrt(1.84)) - 3/2)) .* dipdt_iter  .* mu0 .* 6.2;
% using the same loop voltage for DEMOS (L_iter * dIdt_iter = L_demo * dIdt_demo)
if isfield(sycomore_data,'rampup_dipdt_factor')
    % for use in METIS GUI
    dipdt   = sycomore_data.rampup_dipdt_factor .* (Vloop ./ mu0 ./ R) ./ (CEjima + (log(8 .* R ./ a ./ sqrt(K)) - 3/2));
else
    dipdt   = (Vloop ./ mu0 ./ R) ./ (CEjima + (log(8 .* R ./ a ./ sqrt(K)) - 3/2));
end
fprintf('dIp/dt_{rampup} = %g (kA/s)\n',dipdt/1e3);
% we have to estimate the ramp down current decrease rate 
% to be ajusted to have no flux consumption
% there is trade between flux consumption and l_i change
%dipdt_down = (2/3) .* (Vloop ./ mu0 ./ R) ./ (log(8 .* R ./ a ./ sqrt(K)) - 3/2);
% ITER ramp down is between 60 and 300 s for 15 MA , betwween 250 kA/s and 50 kA/s
% we rescale from ITER studies : -0.06 MA/s (DINA/CRONOS simulation).
Vloop_rampdown = 0.06e6 .* ((log(8.*6.2 ./ 2 ./ sqrt(1.84)) - 3/2) .* mu0 .* 6.2); 
if isfield(sycomore_data,'rampdown_dipdt_factor')
    % for use in METIS GUI
    dipdt_down = sycomore_data.rampdown_dipdt_factor .* Vloop_rampdown ./ (mu0 .* R .* (log(8 .* R ./ a ./ sqrt(K)) - 3/2));
else
    dipdt_down = Vloop_rampdown ./ (mu0 .* R .* (log(8 .* R ./ a ./ sqrt(K)) - 3/2));
end
fprintf('dIp/dt_{rampdown, end of plasma} = %g (kA/s)\n',dipdt_down/1e3);

% dynamical parameters
% flattop current
ip_flattop = ip;
% current at witch the full power is applied
ip_full_power = 0.8 .* ip;
% current after breakdown
fip_ini = 15 ./ 0.4; % Like for ITER
ip_ini = ip / fip_ini;
% plasam current for the back transition to L-mode (reduced magnetic energy by a factor 2)
%ip_dn = 0.8 .* ip;
ip_dn = ip ./ sqrt(2);
% density at application of full power
nbar_full_power = (2/3) .* ip_full_power ./ ip .* min(nbar,1e20 .* (ip_full_power / 1e6) ./ (pi.* a .^ 2));
% used of scaling on density for minimal L2H power threshold. reference: F. Ryter N.F. 2014 
nbar_l2h_min    = 0.7e19 .* (ip_full_power ./ 1e6) .^ 0.34 .* a .^ -0.95 .* b0 .^ 0.62 .* (R ./ a) .^ 0.4;
nbar_full_power = min(nbar_full_power,nbar_l2h_min);
% denstity at the start of flatop
nbar_flattop    = min(nbar,1e20 .* (ip / 1e6) ./ (pi.* a .^ 2));
% density after breakdown
fr_a_ini = 0.5; % fraction of minor radius after breakdown
a_ini = fr_a_ini .* a;
R_ini = R .* (1 - fr_a_ini .* a/R );
q95_ini = 5 .* a_ini .^ 2 .* b0 ./ (ip_ini./1e6) ./ R_ini .* (1.17 - 0.65 .* a_ini ./ R_ini) ./ (1 - (a_ini./R_ini) .^ 2 ) .^ 2;
negr_ini = 1;
nini  = min(negr_ini .* 1e20 .* (ip_ini / 1e6) ./ (pi.* a_ini .^ 2), ...
            1e20 .* (b0 ./ q95_ini ./ R_ini) .^ 0.6  .*  0.25);


% times of interrest
% start of the simulation
t_start = 1.5;
% rampup duration
dt_ramp_up   = ip./ dipdt;
% X point formation (as soon as possible, scale on  ITER scenario)
t_xpoint     = dt_ramp_up / 10;
% plasma current at  X-point formation
ip_xpoint    = ip_ini + t_xpoint .* dipdt;
% time for full power
t_full_power = dt_ramp_up .* ip_full_power ./ ip;
% start of the flattop
t_flattop    = dt_ramp_up;
% time for full density (just a estimation)
t_flattop_plus = t_flattop + 3 .* dt_ramp_up;


% density at x-point transition
q95_x = 5 .* a .^ 2 .* b0 ./ (ip_xpoint ./ 1e6) ./ R .* (1 + K .^ 2 .*  ...
       (1 + 2 .* d .^ 2 - 1.2 .* d .^ 3) ./ 2) .* (1.17 - 0.65 .* a ./ R) ./  ...
       (1 - (a./R) .^ 2 ) .^ 2;
n_xpoint  = max(negr_ini .* 1e20 .* (ip_xpoint ./ 1e6) ./ (pi.* a .^ 2), ...
            1e20 .* (b0 ./ q95_x ./ R) .^ 0.6  .*  0.25);
n_xpoint = max(nini,n_xpoint);
nbar_full_power = max(n_xpoint,nbar_full_power);

% data for rampdown
if isfield(sycomore_data,'rampup_dipdt_factor') && isfield(sycomore_data,'rampdown_dipdt_factor')
  dipdt_adj = dipdt ./ sycomore_data.rampup_dipdt_factor .* sycomore_data.rampdown_dipdt_factor;
else
  dipdt_adj = dipdt;
end
% simulation duration
% end time (wee look for a pulse duration greater than 2 hours)
% duration is increased to be sure to capture the time when all poloial flux is consummed
tend = t_flattop_plus + (3/2) .* sycomore_data.duration;
% time for back transition in L-mode
% there is trade between flux consumption and l_i change
tdn          = tend + (ip_flattop - ip_dn) ./ (dipdt_adj + dipdt_down) .* 2;
% time to back transition to limiter 
tdip         = tdn  +  (ip_dn - ip_xpoint) ./ dipdt_down;
% end of ramp-down simulation (as fast as possible => limter mode)
%tv0          = tdip + (ip_xpoint - ip_ini) ./ dipdt_down ./ 2;
tv0          = tdip + (ip_xpoint - ip_ini) ./ dipdt_adj;

disp('Key times:');
fprintf('Starting time of the simulation: %g s @ Ip = %g MA\n',t_start,ip_ini/1e6);
fprintf('X-point formation: %g s @ Ip = %g MA\n',t_xpoint,ip_xpoint/1e6);
fprintf('Full power time: %g s @ Ip = %g MA\n',t_full_power,ip_full_power/1e6);
fprintf('Start of flat-top: %g s @ Ip = %g MA\n',t_flattop,ip_flattop/1e6);
fprintf('End of full power assited phase: %g s @ Ip = %g MA\n',t_flattop_plus,ip_flattop/1e6);
fprintf('End of flat-top: %g s @ Ip = %g MA\n',tend,ip_flattop/1e6);
fprintf('H to L back transition: %g s @ Ip = %g MA\n',tdn,ip_dn/1e6);
fprintf('back transition to limiter: %g s @ Ip = %g MA\n',tdip,ip_xpoint/1e6);
fprintf('End of simulation/plasma  : %g s @ Ip = %g MA\n',tv0,ip_ini/1e6);

% Parameters for Heating sources
% NBI
rtang = 5.2952 ./ 6.2 .* R;
% parametres chauffage
zext1 = 0.376 ./ 2 ./ 1.84;
zext2 = 0.950 ./ 2 ./ 1.84;

% ICRH
rres  = 6.05 ./ 6.2 .* R;
bres  = b0 .* R ./ rres; 
ag    = 3;
zg    = 1;
harm  = 2;
freq  = harm .* bres .* (95.5e6 .* zg ./ ag) ./ (2 .* pi .* 1e6); % MHz

% LH
freqlh = 5; %GHz
wlh    =  0.5390 .* 3.7 ./ freqlh ./ 6.2 .* R;
npar0  = 1.9 ./ 1e40 .* nbar .^ 2;

% time slices for the computation
temps    = union(union(linspace(t_start,2 .* t_flattop,101), linspace(2 .* t_flattop,tend,301)),linspace(tend,tv0,71))';

% pre-parametrised METIS input data
z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,pecrh2,picrh,pecrh1,pnbi1+pnbi2,zeff,xece,hmore,iso,ftnbi,0);
% load standard reference parameters
z0dinput.option = option;

% local change in parameters
% during early ramp up, follows natural density scaling law
z0dinput.option.neasser = 1;
z0dinput.option.natural  = 1;
% pedestal confinement goes like global confinement
%z0dinput.option.fpped = sycomore_data.H_H; % now included in METIS
% switch off ip control during ramp down
z0dinput.option.vref = 0;
z0dinput.option.tswitch = tv0;
z0dinput.option.vloop = 0;
% decoupled model for radiation			 
z0dinput.option.sol_rad = 'decoupled';
% transition in mode H : prevent transition before X point formation (in limiter the threshold is 2 time that in deverted mode)
z0dinput.option.fpl2h_lim = 2;
% offset to prevent to early H mode transition after X point formation (in the error bar of scaling for L 2 H transition)
z0dinput.option.l2hmul = (p0_ecrh1 + p0_ecrh2 + ip_ini) ./ 1e6;

% NBI is always co-current
z0dinput.option.angle_nbi = 90;
% Rtangence is computed depending on machine size
z0dinput.option.rtang = rtang;
z0dinput.option.zext  = zext1; % given as an input parameter
z0dinput.option.einj  = sycomore_data.nbi_e_inj;          % neutral beam energy injection.

% always 2 NBI
z0dinput.option.nb_nbi = 2;
z0dinput.option.angle_nbi2 = 90;
z0dinput.option.rtang2 = rtang;
z0dinput.option.zext2  = zext2;  % given as an input parameter
z0dinput.option.einj2  = sycomore_data.nbi_e_inj;          % neutral beam energy injection.

% No LHCD in DEMO
% LHCD is used as a second ECRH system
z0dinput.option.lhmode = 5;
z0dinput.option.etalh  = 1;
z0dinput.option.xlh    = 0;     
z0dinput.option.dlh    = 0.4;      
z0dinput.option.angle_ece = z0dinput.option.angle_ece;    

% ECCD or ECRH to assist L-> H transition depending on inductive or non inductive shot
if isfield(sycomore_data,'PECRH') && (sycomore_data.PECRH >= 0)
    z0dinput.option.sens = 1;
elseif sycomore_data.f_ni > 0.9
    z0dinput.option.sens = 1;
else
    z0dinput.option.sens = 0;
end
% NO ICRH in DEMO
% ICRH is configured as for ITER
% if needed, given in input parameters
z0dinput.option.fwcd = 0;
z0dinput.option.mino = 'T';
z0dinput.option.cmin = 1;
z0dinput.option.nphi = 25;
z0dinput.option.freq =freq;
% depending of DEMO pulsed or stady state with or without ITB
% given in input in the set of option parameters

% plasma initialisation
z0dinput.option.breakdown  = -breakdown;

% random shot number
z0dinput.option.shot = sycomore_data.shot;
% decoration
z0dinput.option.machine = tokamak;
% for backward compatibility
z0dinput.machine = z0dinput.option.machine;
z0dinput.shot = z0dinput.option.shot;

% flux consumption
z0dinput.option.available_flux = sycomore_data.available_flux;


% other parameters comming from Sycomore
if isfield(sycomore_data,'ne_peak') && ~isempty(sycomore_data.ne_peak)
    z0dinput.option.ane = 4;
    z0dinput.option.vane =  sycomore_data.ne_peak;
    disp('Density profile peaking factor provides by SYCOMORE');
end
z0dinput.option.scaling     = sycomore_data.scaling;
z0dinput.option.rimp        = sycomore_data.rimp;
z0dinput.option.Recycling   = sycomore_data.Recycling;
% leave METIS compute taken into account Recycling
%z0dinput.option.tauhemul    = -3;  % will be added recycling effect
z0dinput.option.tauhemul    = sycomore_data.tau_He_o_tau_E_core;
z0dinput.option.rw          = sycomore_data.rw;
z0dinput.option.fR_target   = sycomore_data.fR_target;
z0dinput.option.signe       = sycomore_data.signe;     
z0dinput.option.cw_offset   = sycomore_data.cW;   
z0dinput.option.fzmax_div   = sycomore_data.fzmax_div; 
z0dinput.option.factor_scale = sycomore_data.f_DSOL;  % to be compatible with Sycomore; information must be directly transmited from Sycomore          
z0dinput.option.nea_factor   = sycomore_data.f_ne_LCFS;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore



% waveforms control nodes generation
tv      = [];     % time  nodes
ipv     = [];     % plasma currentnodes
nv      = [];     % line averaged density nodes
picrhv  = [];     % ICRH power nodes
plhv    = [];     % LHCD/ECRH2 power nodes
pecrhv  = [];     % ECRH power nodes
pnbi1v  = [];     % NBI1 power nodes
pnbi2v  = [];     % NBI2 power nodes
Rv      = [];     % major radius nodes
av      = [];     % minor radius nodes
Kv      = [];     % averaged elongation nodes
dv      = [];     % averaged triangularity nodes
z0v     = [];     % plasma vertical shift nodes
isov    = [];     % isotopic ration nT/nD nodes
xecev   = [];     % ECRH maximum depostion position nodes
hmore   = [];     % confinement enhancement nodes

% first time after breakdown
tv(end+1)      = 1.5;
ipv(end+1)     = ip_ini;
nv(end+1)      = nini;
picrhv(end+1)  = 0;
plhv(end+1)    = p0_ecrh2;
pecrhv(end+1)  = p0_ecrh1;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = R_ini ./ R;
av(end+1)      = a_ini ./ a;
Kv(end+1)      = 0;
dv(end+1)      = 0;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

% time for X point formation
tv(end+1)      = t_xpoint;
ipv(end+1)     = ip_xpoint;
nv(end+1)      = n_xpoint;
picrhv(end+1)  = 0;
plhv(end+1)    = p1_ecrh2;
pecrhv(end+1)  = p1_ecrh1;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

%  % time just before full power
%  tv(end+1)      = t_before;
%  ipv(end+1)     = ip_before;
%  nv(end+1)      = 3 .* nini;
%  picrhv(end+1)  = 1;
%  plhv(end+1)    = p2_ecrh2;
%  pecrhv(end+1)  = p2_ecrh1;
%  pnbi1v(end+1)  = 0;
%  pnbi2v(end+1)  = 0;
%  Rv(end+1)      = 1;
%  av(end+1)      = 1;
%  Kv(end+1)      = 1;
%  dv(end+1)      = 1; 
%  z0v(end+1)     = sepa_option.za;
%  isov(end+1)    = 1;
%  xecev(end+1)   = xece_itb;
%  hmore(end+1)   = 1;

% time for full power and H-mode transition
tv(end+1)      = t_full_power;
ipv(end+1)     = ip_full_power;
nv(end+1)      = nbar_full_power;
picrhv(end+1)  = 1;
plhv(end+1)    = p2_ecrh2;
pecrhv(end+1)  = p2_ecrh1;
pnbi1v(end+1)  = 1;
pnbi2v(end+1)  = 1;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = xece_itb;
hmore(end+1)   = sycomore_data.H_H;

% start of flattop
tv(end+1)      = t_flattop;
ipv(end+1)     = ip_flattop;
nv(end+1)      = nbar_flattop;
picrhv(end+1)  = 1;
plhv(end+1)    = pecrh_plus2;
pecrhv(end+1)  = pecrh_plus1;
pnbi1v(end+1)  = 1;
pnbi2v(end+1)  = 1;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = xece;
hmore(end+1)   = sycomore_data.H_H;

% full density and decrease of auxiliary heating
tv(end+1)      = t_flattop_plus;
ipv(end+1)     = ip_flattop;
nv(end+1)      = nbar;
picrhv(end+1)  = 1;
plhv(end+1)    = pecrh2;
pecrhv(end+1)  = pecrh1;
pnbi1v(end+1)  = 1;
pnbi2v(end+1)  = 1;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = xece;
hmore(end+1)   = sycomore_data.H_H;

% end of flattop
tv(end+1)      =  tend;
ipv(end+1)     =  ip_flattop;
nv(end+1)      =  nbar;
picrhv(end+1)  =  0;
plhv(end+1)    =  pecrh2;
pecrhv(end+1)  =  pecrh1;
pnbi1v(end+1)  =  1;
pnbi2v(end+1)  =  1;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  1;
xecev(end+1)   =  xece;
hmore(end+1)   =  sycomore_data.H_H;;

% shutdown NBI (and ICRH if present) 50%
% during start of rampdown
tv(end+1)      =  tend + (tdn - tend) ./ 3;
ipv(end+1)     =  ip_flattop - (ip_flattop - ip_dn) ./ 3;
nv(end+1)      =  nbar;
picrhv(end+1)  =  0;
plhv(end+1)    =  pecrh2;
pecrhv(end+1)  =  p1_ecrh1;
pnbi1v(end+1)  =  1;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  1;
xecev(end+1)   =  xece;
hmore(end+1)   =  1;

% shutdown NBI (and ICRH if present) 100%
% during start of rampdown
tv(end+1)      =  tend + 2 .* (tdn - tend) ./ 3;
ipv(end+1)     =  ip_flattop - 2 .* (ip_flattop - ip_dn) ./ 3;
nv(end+1)      =  nbar;
picrhv(end+1)  =  0;
plhv(end+1)    =  p1_ecrh2;
pecrhv(end+1)  =  p1_ecrh1;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  1;
xecev(end+1)   =  xece;
hmore(end+1)   =  1;

% back transition to L-mode
tv(end+1)      =  tdn;
ipv(end+1)     =  ip_dn;
nv(end+1)      =  nini;
picrhv(end+1)  =  0;
plhv(end+1)    =  0;
pecrhv(end+1)  =  0;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  1;
xecev(end+1)   =  0;
hmore(end+1)   =  1;

% back transition to limiter
tv(end+1)      =  tdip;
ipv(end+1)     =  ip_xpoint;
nv(end+1)      =  nini;
picrhv(end+1)  =  0;
plhv(end+1)    =  0;
pecrhv(end+1)  =  0;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  0.7;
xecev(end+1)   =  0;
hmore(end+1)   =  1;

% end of controlled ramp-down
tv(end+1)      =  tv0;
ipv(end+1)     =  ip_ini;
nv(end+1)      =  nini;
picrhv(end+1)  =  0;
plhv(end+1)    =  0;
pecrhv(end+1)  =  0;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1 - 0.7 .* sepa_option.a ./ sepa_option.ra .* sepa_option.rxdo;
av(end+1)      =  0.5;
Kv(end+1)      =  0;
dv(end+1)      =  0;
z0v(end+1)     =  zlow;
isov(end+1)    =  0.5;
xecev(end+1)   =  0;
hmore(end+1)   =  1;

% scales
fprintf('power scale to ITER: %g\n',rap_power);
picrhv  = picrhv  .* picrh .* rap_power;
plhv    = plhv    .* rap_power;
pecrhv  = pecrhv  .* rap_power;
if isfield(sycomore_data,'PECRH') && (sycomore_data.PECRH >= 0)
    pecrhv   = min(sycomore_data.PECRH,pecrhv)
end
if sycomore_data.PNBI >= 0
  pnbi1v  = pnbi1v  .* sycomore_data.PNBI ./ 2;
  pnbi2v  = pnbi2v  .* sycomore_data.PNBI ./ 2 ;
else
  pnbi1v  = pnbi1v  .* pnbi1 .* rap_power;
  pnbi2v  = pnbi2v  .* pnbi2 .* rap_power;
end
Rv      = Rv .* R;
av      = av .* a;
Kv      = Kv .* (K - 1) + 1;
dv      = dv .* d;

% time interpolation
z0dinput.cons.ip = pchip(tv,ipv,temps);
z0dinput.cons.picrh = zinterpnc(tv,picrhv,temps);
z0dinput.cons.pnbi = zinterpnc(tv,pnbi1v,temps) + sqrt(-1) .* zinterpnc(tv,pnbi2v,temps);
if onnbi2 == 1
 %ind_late = find(temps <= (t_full_power + t_xpoint));
 ind_late = find((temps >= (t_full_power - t_xpoint)) & (temps <= t_full_power));
 z0dinput.cons.pnbi(ind_late) = max(real(z0dinput.cons.pnbi(ind_late + 1)));
end
z0dinput.cons.plh = zinterpnc(tv,plhv,temps);
z0dinput.cons.pecrh = zinterpnc(tv,pecrhv,temps);
indft = find((temps > t_xpoint) & (temps < t_full_power));
z0dinput.cons.pecrh(indft) = pchip(tv,pecrhv,temps(indft));
z0dinput.cons.iso = pchip(tv,isov,temps);
if z0dinput.option.nb_nbi > 1
	z0dinput.cons.ftnbi = (1 + sqrt(-1)) .* ftnbi .* ones(size(z0dinput.cons.ftnbi));
else
	z0dinput.cons.ftnbi = ftnbi .* ones(size(z0dinput.cons.ftnbi));
end
z0dinput.geo.a   = pchip(tv,av,temps);
z0dinput.geo.R   = pchip(tv,Rv,temps);
z0dinput.geo.K   = pchip(tv,Kv,temps);
z0dinput.geo.d   = pchip(tv,dv,temps);
z0dinput.geo.z0  = pchip(tv,z0v,temps);
z0dinput.cons.nbar = pchip(tv,nv,temps);
z0dinput.cons.xece = pchip(tv,xecev,temps);
z0dinput.cons.hmore = zinterpnc(tv,hmore,temps);
% calcul de l'evolution de zeff (nimp = Ct, ne = K * ip) ; just to fill the field
% Zeff increase up to X point  formation and get it reference value on flatop due to seeding injection
if (z0dinput.option.W_effect > 0 ) || (z0dinput.option.zimp == 74) || (z0dinput.option.zmax == 74)
      % Zeff during early ramp up starting for pure DT > DT + impurities : with W we hope to have Zeff = 1.3
      z0dinput.cons.zeff =  1.2 - 0.2 .* max(0,(z0dinput.cons.temps - t_xpoint) ./ (z0dinput.cons.temps(1) - t_xpoint));
      % zeff for ramp down
      zeff_dens =  min(zeff - 1,0.2) .*  min(ip_xpoint  ./ z0dinput.cons.ip,max(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) + 1;
else
      % Zeff during early ramp up starting for pure DT > DT + impurities : with C  we hope to have Zeff = 2.3
      z0dinput.cons.zeff =  min(zeff,2.3) - min(zeff - 1,1.3) .* max(0,(z0dinput.cons.temps - t_xpoint) ./ (z0dinput.cons.temps(1) - t_xpoint));
      % Zeff for ramp down
      zeff_dens =  (zeff - 1) .*  min(ip_xpoint  ./ z0dinput.cons.ip,max(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) + 1;
end
% we have Zeff after early ramp up that depend on input power (impurities sources and seeding to protect divertor)
% we assume Q=25 drived by NBI
ffus  =  1 + 2 .* double((z0dinput.cons.temps > ((3/4) .* t_flattop + (1/4) .* t_flattop_plus)) & (z0dinput.cons.temps < tdn)) + 2 .*  (z0dinput.cons.nbar ./ max(z0dinput.cons.nbar)) .^ 2;
%pin   = z0dinput.cons.ip + z0dinput.cons.picrh + z0dinput.cons.plh + z0dinput.cons.pecrh + ffus .* real(z0dinput.cons.pnbi)  +  ffus .* imag(z0dinput.cons.pnbi);
pin   = z0dinput.cons.ip + ffus .* real(z0dinput.cons.pnbi)  +  ffus .* imag(z0dinput.cons.pnbi);
fzeff = min(1, pin ./ max(real(z0dinput.cons.pnbi)  + imag(z0dinput.cons.pnbi) + 1) ./ max(ffus)) .* (z0dinput.cons.temps > t_xpoint);
% seeding effect and ramp-down
%  fzeff = ones(size(z0dinput.cons.zeff));
%  fzeff(z0dinput.cons.temps < t_full_power) = 0;
fzeff = sgolayfilt(fzeff,1,3);
z0dinput.cons.zeff = z0dinput.cons.zeff .* (1 - fzeff) + zeff .* fzeff;
% Zeff in unchanged during ramp-down up to the limited phase; the fuelling is limited during this phase, it can be difficult to remove impurities
z0dinput.cons.zeff(z0dinput.cons.temps > tend) = z0dinput.cons.zeff(z0dinput.cons.temps == tend);
% Zeff evolution during limited final phase
fzeff = ones(size(z0dinput.cons.zeff));
fzeff(z0dinput.cons.temps > tdip) = zeff_dens(z0dinput.cons.temps > tdip) < zeff;
fzeff = sgolayfilt(fzeff,1,3);
z0dinput.cons.zeff = zeff_dens .* (1 - fzeff) + z0dinput.cons.zeff .* fzeff;
% limits for formula inside METIS 
z0dinput.cons.zeff = max(1,min(7, sgolayfilt(z0dinput.cons.zeff,1,3)));


% LCFS computation
%  if 0
%  sepa_option.ton       = t_xpoint;
%  sepa_option.toff      = tdip;
%  z0dinput = z0separatrix(z0dinput,sepa_option,0);
%  z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;
%  
%  % print separtrix parameters:
%  disp('Separatrix parameters (flat top phase):');
%  sepa_option
%  
%  % decrease of elongation during start of rampdown
%  indkreduc = find(z0dinput.cons.temps > tend);
%  % reference xpoint position
%  z_xpoint_ref = min(z0dinput.exp0d.Zsepa(indkreduc(1) - 1,:)) + z0dinput.geo.z0(indkreduc(1) - 1);
%  % minimal value of elongation after decrease
%  Kmin_reduc = (sepa_option.zxdo + sepa_option.zxup) ./ 2;
%  % try to keep q_95 = Constant during first part of ramp down and try to increase vertical stability
%  Kreduc    = sqrt(max(Kmin_reduc,(K .^ 2 + 1) .* z0dinput.cons.ip ./ ip_flattop - 1));
%  zxdo_mem = sepa_option.zxdo;
%  zxup_mem = sepa_option.zxup;
%  rxup_mem = sepa_option.rxup;
%  rxdo_mem = sepa_option.rxdo;
%  
%  
%  for k = indkreduc'
%    sepa_option.zxup = Kreduc(k) ./ K .* zxup_mem;
%    sepa_option.zxdo = Kreduc(k) ./ K .* zxdo_mem;
%    sepa_option.rxup = rxup_mem .* z0dinput.cons.ip(k) ./ ip_flattop;
%    sepa_option.rxdo = rxdo_mem .* z0dinput.cons.ip(k) ./ ip_flattop;
%    z0d1t = zerod_get1t(z0dinput,k);
%    z0d1t.geo.K = min(z0d1t.geo.K,Kreduc(k));
%    rep = z0separatrix(z0d1t,sepa_option,0);
%    noms = fieldnames(rep.geo);
%    for l=1:length(noms)
%        z0dinput.geo.(noms{l})(k) = rep.geo.(noms{l});
%    end
%    z0dinput.exp0d.Rsepa(k,:) = rep.exp0d.Rsepa(1,:);
%    % kept xpoint position fixed
%    z0dinput.exp0d.Zsepa(k,:) = rep.exp0d.Zsepa(1,:) - min(rep.exp0d.Zsepa(1,:)) - z0d1t.geo.z0 + z_xpoint_ref;
%    % regularisation of d
%    z0dinput.geo.d(k) = min(z0dinput.geo.d(k),z0dinput.geo.d(k - 1));
%  end
%  z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;
%  end



sepa_option.ton       = t_xpoint;
sepa_option.toff      = tdip;
z0dinput = z0separatrix(z0dinput,sepa_option,0);
z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;

% print separtrix parameters:
disp('Separatrix parameters (flat top phase):');
sepa_option

% decrease of elongation during start of rampdown
indkreduc = find(z0dinput.cons.temps > tend);
% reference xpoint position
z_xpoint_ref = min(z0dinput.exp0d.Zsepa(indkreduc(1) - 1,:)) + z0dinput.geo.z0(indkreduc(1) - 1);
% minimal value of elongation after decrease
% try to keep q_95 = Constant during first part of ramp down and try to increase vertical stability
fKreduc    = sqrt(max(1,(K .^ 2 + 1) .* z0dinput.cons.ip ./ ip_flattop - 1)) ./ K;

zxdo_mem = sepa_option.zxdo;
zxup_mem = sepa_option.zxup;
rxup_mem = sepa_option.rxup;
rxdo_mem = sepa_option.rxdo;
apup_mem = sepa_option.apup;
amup_mem = sepa_option.amup;
Kref = z0dinput.geo.K(indkreduc(1) - 1);
dref = z0dinput.geo.d(indkreduc(1) - 1);

for k = indkreduc'
    sepa_option.zxup = max(1,fKreduc(k) .* zxup_mem);
    sepa_option.zxdo = max(1,fKreduc(k) .* zxdo_mem);
    sepa_option.rxup = rxup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.rxdo = rxdo_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.apup = apup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.amup = amup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    z0d1t = zerod_get1t(z0dinput,k);
    z0d1t.geo.K = min(z0d1t.geo.K,max(1,Kref .* fKreduc(k)));
    z0d1t.geo.d = min(z0d1t.geo.d,dref .* z0dinput.cons.ip(k) ./ ip_flattop);
    rep = z0separatrix(z0d1t,sepa_option,0);
    z0dinput.geo.K(k) = rep.geo.K;
    z0dinput.geo.d(k) = rep.geo.d;
    z0dinput.exp0d.Rsepa(k,:) = rep.exp0d.Rsepa(1,:) - (max(rep.exp0d.Rsepa(1,:)) + min(rep.exp0d.Rsepa(1,:))) ./ 2 + z0d1t.geo.R;
    % kept xpoint position fixed
    z0dinput.exp0d.Zsepa(k,:) = rep.exp0d.Zsepa(1,:) - min(rep.exp0d.Zsepa(1,:)) - z0d1t.geo.z0 + z_xpoint_ref;
    % regularisation of d
    z0dinput.geo.d(k) = max(0,min(z0dinput.geo.d(k),z0dinput.geo.d(k - 1)));
    z0dinput.geo.K(k) = max(1,min(z0dinput.geo.K(k),z0dinput.geo.K(k - 1)));
end

%  % graph for separatrix 
%  h = findobj(0,'type','figure','tag','z0geosepa');
%  if isempty(h)
%      h=figure('tag','z0geosepa');
%  else
%      figure(h);
%  end
%  clf
%  plot(z0dinput.exp0d.Rsepa',(z0dinput.exp0d.Zsepa + z0dinput.geo.z0 * ones(1,size(z0dinput.exp0d.Zsepa,2)))');
%  axis('square')
%  axis('equal')

% limitation of Zeff
% evolution of Zeff following gaz density 
ulh = 0.25;
fh  = zeros(size(temps));
fh(temps >= t_full_power) = 1;
negr     = (z0dinput.cons.ip ./ 1e6) ./ z0dinput.geo.a .^ 2 ./ pi .* 1e20;
% d95  = d .* 0.95 .^ 2;         (de varie en x^ 2)
% formule a partir de k0 = (K+1)/2 et K(x) = k0 + (K-k0) x^ 4
K95 = z0dinput.geo.K .* (1 + 0.95 .^ 4 ) ./ 2 + 0.5 .* (1 - 0.95 .^4);
d95 = z0dinput.geo.d .* (0.95 .^ 2);
q95 = 5 .* z0dinput.geo.a .^ 2 .* z0dinput.geo.b0 ./ (ip./1e6) ./ z0dinput.geo.R .* (1 + K95 .^ 2 .*  ...
     (1 + 2 .* d95 .^ 2 - 1.2 .* d95 .^ 3) ./ 2) .* (1.17 - 0.65 .* z0dinput.geo.a ./ z0dinput.geo.R) ./  ...
      (1 - (z0dinput.geo.a./z0dinput.geo.R) .^ 2 ) .^ 2;
nbar_nat = min(negr,max(1e13, 1e20 .* (sycomore_data.rb0 ./ q95) .^ 0.6  .*  (ulh + (1 - ulh) .* fh) .* z0dinput.option.fnbar_nat));
nbar_star     = max(z0dinput.cons.nbar,nbar_nat);
zeff_max           =  min((zeff - 1) .* max(nbar_star)  ./ nbar_star + 1,max(z0dinput.option.zmax,z0dinput.option.zimp));
z0dinput.cons.zeff =  min(zeff_max,z0dinput.cons.zeff);
%figure(11);clf;plot(z0dinput.cons.temps,z0dinput.cons.zeff);drawnow

% tune parameter to have a feedback on NBI power when all current must be non inductive
% dont't work with fusion power as main additionnal heating
%  if sycomore_data.f_ni   == 1              
%    z0dinput.option.vloop     = 6;
%    z0dinput.option.vref      = 0.1;     
%    z0dinput.option.tswitch   = t_flattop_plus; 
%  end


% if only input structure z0dinput is required
if nargout == 1
  return
end

% loop on H mode transition
switch z0dinput.option.sol_model
case '2_points'
    z0dinput.option.sol_model    = 'scaling'
    m2points = 0;
otherwise 
    m2points = 1;
end
% first computation
indHtest   = find(z0dinput.cons.temps > t_flattop_plus,1);
ind_overshoot = find((z0dinput.cons.temps > t_flattop) & (z0dinput.cons.temps < t_flattop_plus));
ind_flattop   = find((z0dinput.cons.temps > t_flattop_plus) & (z0dinput.cons.temps < tend));
indnoHtest = find((z0dinput.cons.temps < t_flattop) & (z0dinput.cons.pecrh >  0) & (z0dinput.cons.pecrh <  max(z0dinput.cons.pecrh)));
indnodisrup= find((z0dinput.cons.temps > tend) & (z0dinput.cons.temps <= tdip));
indoverdrive = find((z0dinput.cons.temps < t_flattop) & (z0dinput.cons.plh >  0));
indrampup = find((z0dinput.cons.temps < t_flattop));
indmax_pecrh = find(z0dinput.cons.pecrh ==  max(z0dinput.cons.pecrh));
ind_push_hmode = find(z0dinput.cons.temps >= t_full_power,1):max(find(z0dinput.cons.pecrh ==  max(z0dinput.cons.pecrh)));
% flag set to one if nothing can be change due to other rules
nochange  = 0;
nochange2 = 0;
nochange3 = 0;
nochange4 = 0;
% memorize original NBI power
pnbi_mem = z0dinput.cons.pnbi;
%fpnbi_fact = 2;
fpnbi_fact = sqrt(2);
% loop for optimisation
for k=1:31
    fprintf('=====> start of optimization round %d (over maximum of 31) <=====\n',k);
    % attenuation de la variation de  lapuissance nbi
    if k > 21
	fpnbi_fact = fpnbi_fact ./ sqrt(2);
    end
    % call of METIS
    [zs,infovoid,profli] = zerodfast(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
    post.z0dinput = z0dinput;
    post.zerod    = zs;
    post.profil0d =profli;
    z0plotsc;drawnow;
    % if optimisation is not required
    if testing == 2
	  disp('No optimisation loop: testing mode required without optimisation');
	  break;
    end
    % index for rampdown disrution avoidance
    indnodisrup = find((z0dinput.cons.temps > tend) & (post.zerod.ip >= ip_xpoint) & (z0dinput.cons.temps <= tdip));
    if isempty(indnodisrup)
	indnodisrup = find(z0dinput.cons.temps > tend,1);
    end
    % overdrive testing
    test_overdrive = double(((zs.ini(indoverdrive) >= zs.ipar(indoverdrive)) | (zs.li(indoverdrive) < 0.5)) .* (k <= 11));
    % over shoot power
    test_overshoot = double((max(zs.pfus(ind_overshoot)) > mean(zs.pfus(ind_flattop))) .* (k <= 7));
    % start of optimisation part
    if (zs.modeh(indHtest) >= 1) && ~any(zs.modeh(indnoHtest) >= 1) && (all(zs.disrup(indnodisrup) == 0) || (nochange3 == 1)) &&  ...
       (~any(test_overdrive) || (nochange == 1)) && (~any(test_overshoot) || (nochange2 == 1)) && (all(zs.disrup(indrampup) == 0) || (nochange4 == 1))
      if m2points == 2
	  break;
      elseif m2points == 0
	  m2points = 1;
	  z0dinput.option.sol_model = '2_points'; 
 	  disp('start run with 2 points model');
      else
	  % one more run for NBI optimisation
	  m2points = 2;
      end
    end
    if (zs.modeh(indHtest)  < 1)
	disp('L2H transition missed: new try with increased ECRH power');
    end
    if any(zs.modeh(indnoHtest)  >= 1)
	disp('to early L2H transition: new try with adjusted ECRH power');
    end
    if any(zs.disrup(indnodisrup) > 0) && (nochange3 == 0)
 	disp('diruption during ramp-down: new try with adjusted ECRH power');
    end
    if any(zs.disrup(indrampup) > 0) && (nochange4 == 0)
 	disp('diruption during ramp-up: new try with adjusted ECRH2 (in LH) power');
    end
    if any(test_overdrive > 0) && (nochange == 0)
 	disp('overdrive during ramp-up: new try with adjusted ECRH2 (in LH) power');
    end
    if any(test_overshoot > 0) && (nochange2 == 0)
 	disp('fusion power overshoot: new try with adjusted ECRH power');
    end
    if (zs.modeh(indHtest) < 1)
      %indmax_pecrh
      delta_pecrh = post.zerod.plhthr(indHtest) - (post.zerod.plossl2h(indHtest) + 1e6 .*  post.z0dinput.option.l2hmul);
      if delta_pecrh <= 0
	    delta_pecrh = 0.2 .* max(z0dinput.cons.pecrh);
      end 
      %z0dinput.cons.pecrh(indmax_pecrh) = z0dinput.cons.pecrh(indmax_pecrh) .* (1 + delta_pecrh ./ max(1,max(z0dinput.cons.pecrh)));
      z0dinput.cons.pecrh(ind_push_hmode) = z0dinput.cons.pecrh(ind_push_hmode) .* (1 + delta_pecrh ./ max(1,max(z0dinput.cons.pecrh)));
    elseif any(test_overshoot)
      % don't really work : no known solution at this probl?me
      % overshoot of fusion power after transition to H mode_expo_inte
      % try to reduce the margin for L to H transition
      delta_pecrh = post.zerod.plhthr(indmax_pecrh) - 1.5 .* (post.zerod.plossl2h(indmax_pecrh) + 1e6 .*  post.z0dinput.option.l2hmul);
      z0dinput.cons.pecrh(indmax_pecrh) = z0dinput.cons.pecrh(indmax_pecrh) .* max(0.3,1 - max(0,min(delta_pecrh, ...
				    (max(zs.pfus(ind_overshoot)) - 1.2 .* mean(post.zerod.pfus(ind_flattop)))) ./ max(1,max(post.zerod.pin))));
      if all(z0dinput.cons.pecrh == post.z0dinput.cons.pecrh)
	  nochange2 = 1;
      else 
	  nochange2 = 0;
      end
    else
      nochange2 = 1;
    end 
%      if  (zs.modeh(indHtest) >= 1)
%  	% change tswich to opimize ramp down
%  	%post.z0dinput.cons.temps(max(find(post.zerod.modeh)))
%  	%z0dinput.option.tswitch = max(tend,post.z0dinput.cons.temps(max(find(post.zerod.modeh))));
%  	z0dinput.option.tswitch = max(tdn,post.z0dinput.cons.temps(max(find(post.zerod.modeh))));
%      end
    if any(zs.modeh(indnoHtest)  >= 1)
	ind_change = indnoHtest(zs.modeh(indnoHtest)  > 0);
	ind_change = min(ind_change):max(ind_change);
        %zs.modeh(ind_change)' 
	%delta_pecrh = zs.plhthr - (zs.plossl2h + 1e6 .*  z0dinput.option.l2hmul);
	%delta_pecrh =  sgolayfilt(zs.plhthr - (zs.plossl2h + 1e6 .*  z0dinput.option.l2hmul),1,5);
	%delta_pecrh =  sgolayfilt(max(0,zs.pin - (0.8 .* zs.plossl2h + 1e6 .*  z0dinput.option.l2hmul)),1,5);
	delta_pecrh =  max(0,post.zerod.pin - (0.8 .* post.zerod.plossl2h + 1e6 .*  post.z0dinput.option.l2hmul));
	%z0dinput.cons.pecrh(ind_change)'
	%z0dinput.cons.pecrh(ind_change) = max(1,z0dinput.cons.pecrh(ind_change) - 1.5 .* max(0,delta_pecrh(ind_change)) - 1e6 .* (zs.modeh(ind_change)>=1));
	z0dinput.cons.pecrh(ind_change) = max(1,z0dinput.cons.pecrh(ind_change) - 0.3 .* delta_pecrh(ind_change));
	%z0dinput.cons.pecrh(ind_change)'
	%sum(zs.modeh(indnoHtest))
    else
	 ind_change = [];
    end   
    % ramp down tuning to prevent disruption
    if any(zs.disrup(indnodisrup) > 0)
		%zs.disrup(indnodisrup)'
                pdisrup  = max(zs.prad + zs.pcyclo + zs.pbrem + zs.pioniz - (zs.pin + zs.dwdt .* (zs.dwdt >0)) , ...
                              ((0.01 .* zs.pin) - (zs.pel-min(0,zs.pei))));
                %pdisrup  = max(-zs.plhthr,pdisrup);
                pdisrup  = max(zs.ip ./ 10,pdisrup);
                %pdisrup(indnodisrup)'
                %z0dinput.cons.pecrh(indnodisrup)'
                perch_nodisrup_mem = z0dinput.cons.pecrh(indnodisrup);
		%z0dinput.cons.pecrh(indnodisrup) =   z0dinput.cons.pecrh(indnodisrup) + min(zs.ip(indnodisrup),(1 - 0.7 .* m2points) .* zs.disrup(indnodisrup) .* pdisrup(indnodisrup));
		z0dinput.cons.pecrh(indnodisrup) =   z0dinput.cons.pecrh(indnodisrup) + min(zs.ip(indnodisrup),0.3 .* zs.disrup(indnodisrup) .* pdisrup(indnodisrup));
                %z0dinput.cons.pecrh(indnodisrup)'
	        % maximum available power for rampdown
	        pecrh_max_rd = max(0,max(z0dinput.cons.pecrh(z0dinput.cons.temps < tend) + z0dinput.cons.plh(z0dinput.cons.temps < tend)) - max(z0dinput.cons.plh(indnodisrup)));
	        % reverse loop to spread the energie
		delta_pecrh = 0;
	        for rk = indnodisrup(end - 1):-1:indnodisrup(1)
		    z0dinput.cons.pecrh(rk) = z0dinput.cons.pecrh(rk) + delta_pecrh;
		    delta_pecrh = max(0,z0dinput.cons.pecrh(rk) - pecrh_max_rd);
		    z0dinput.cons.pecrh(rk) = min(pecrh_max_rd,z0dinput.cons.pecrh(rk));
	        end
	        %z0dinput.cons.pecrh(indnodisrup) = min(pecrh_max_rd,z0dinput.cons.pecrh(indnodisrup));
	        if all(perch_nodisrup_mem == z0dinput.cons.pecrh(indnodisrup))
		  nochange3 = 1;
		else
		      %delta_nodisrup = - perch_nodisrup_mem + z0dinput.cons.pecrh(indnodisrup)
		end
                %z0dinput.cons.pecrh(indnodisrup)'
    end
    if any(test_overdrive > 0)
          % temproral shift
          ind_overdrive = setdiff(find(test_overdrive),ind_change);
          ind_overdrive = max(1,min(ind_overdrive) - 1):max(ind_overdrive);
          if ~isempty(ind_overdrive)
              % kept minimum power to prevent disruption during early ramp-up
              % reduce plh (second ECRH launcher in his case) to prevent over current drive and to low li value.           
	      z0dinput.cons.plh(ind_overdrive) = max(p0_ecrh2 + z0dinput.cons.ip(ind_overdrive) ,0.7 .* z0dinput.cons.plh(ind_overdrive));
	      if all(z0dinput.cons.plh == post.z0dinput.cons.plh)
		  nochange = 1;
	      else
		  nochange = 0;
	      end
	  end
    else
	  nochange = 1;	
	  ind_overdrive = [];
    end
    if any(zs.disrup(indrampup) > 0)
	  ind_change_rampup = setdiff(intersect(indrampup,find(zs.disrup>0)),ind_overdrive);
	  pdisrup  = max(zs.prad + zs.pcyclo + zs.pbrem + zs.pioniz - (zs.pin + zs.dwdt .* (zs.dwdt >0)) , ...
			((0.01 .* zs.pin) - (zs.pel-min(0,zs.pei))));
	  %pdisrup  = max(-zs.plhthr,pdisrup);
	  pdisrup  = max(zs.ip ./ 10,pdisrup);
	  plh_nodisrup_mem = z0dinput.cons.plh(ind_change_rampup);
	  if min(ind_change_rampup) > 1
	        delta_plh_ru =  min(zs.ip(ind_change_rampup),0.3 .* zs.disrup(ind_change_rampup) .* pdisrup(ind_change_rampup));
		z0dinput.cons.plh(ind_change_rampup)    =   z0dinput.cons.plh(ind_change_rampup)     + delta_plh_ru;
		z0dinput.cons.plh(ind_change_rampup - 1) =  z0dinput.cons.plh(ind_change_rampup - 1) + delta_plh_ru;
	  else
		z0dinput.cons.plh(ind_change_rampup) =   z0dinput.cons.plh(ind_change_rampup) + min(zs.ip(ind_change_rampup),0.3 .* zs.disrup(ind_change_rampup) .* pdisrup(ind_change_rampup));
          end
          %z0dinput.cons.temps(ind_change_rampup) 
          %z0dinput.cons.plh(ind_change_rampup)
          %z0dinput.cons.plh(ind_change_rampup) - plh_nodisrup_mem
          if all(plh_nodisrup_mem == z0dinput.cons.plh(ind_change_rampup))
	     nochange4 = 1;	
          else
	     nochange4 = 0;	         
          end
    else
	  nochange4 = 1;	
    end
    % filering to remove peak
%      ind_filter = find(z0dinput.cons.temps < t_flattop);
%      z0dinput.cons.plh(ind_filter) = sgolayfilt(z0dinput.cons.plh(ind_filter),1,3);
%      z0dinput.cons.pecrh(ind_filter) = sgolayfilt(z0dinput.cons.pecrh(ind_filter),1,3);
    
    %figure(21);clf;plot(post.z0dinput.cons.temps,post.z0dinput.cons.pecrh,'b',z0dinput.cons.temps,z0dinput.cons.pecrh,'r',z0dinput.cons.temps,z0dinput.cons.plh,'g');drawnow
%      nochange   
%      z0dinput.cons.plh' - post.z0dinput.cons.plh'
%      nochange2   
%      z0dinput.cons.pecrh' - post.z0dinput.cons.pecrh'
%      figure(21);
%      clf;
%      subplot(2,1,1)
%      plot(post.z0dinput.cons.temps,post.z0dinput.cons.plh,'b',z0dinput.cons.temps,z0dinput.cons.plh,'r');
%      subplot(2,1,2)
%      plot(post.z0dinput.cons.temps,post.z0dinput.cons.pecrh,'b',z0dinput.cons.temps,z0dinput.cons.pecrh,'r',z0dinput.cons.temps,z0dinput.cons.plh,'g');
%      drawnow
    %keyboard
    % NBI optimisation
    disp('Optimisation of NBI geometry:')
    % X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
    options_optim = optimset('display','iter','Algorithm','active-set','MaxFunEvals',1e4);
    if sycomore_data.f_ni >= 1
      res = cat(2,post.z0dinput.option.rtang,post.z0dinput.option.rtang2,post.z0dinput.option.zext,post.z0dinput.option.zext2,1,1);
      lb  = [max(z0dinput.geo.R) - 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) -  0.8 .* max(z0dinput.geo.a),0,0,1./fpnbi_fact,1./fpnbi_fact];
      ub  = [max(z0dinput.geo.R) + 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) +  0.8 .* max(z0dinput.geo.a),0.7,0.7,fpnbi_fact,fpnbi_fact];
    else
      res = cat(2,post.z0dinput.option.rtang,post.z0dinput.option.rtang2,post.z0dinput.option.zext,post.z0dinput.option.zext2);
      lb  = [max(z0dinput.geo.R) - 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) -  0.8 .* max(z0dinput.geo.a),0,0];
      ub  = [max(z0dinput.geo.R) + 0.8 .* max(z0dinput.geo.a),max(z0dinput.geo.R) +  0.8 .* max(z0dinput.geo.a),0.7,0.7];
    end
    [res,objective_end,exit_flag] = fmincon(@(x) z0optimnbigeo(x,post,sycomore_data.shine_through_limit),res,[],[],[],[],lb,ub,'',options_optim);
    [objective,rtang,rtang2,zext,zext2,ialign,fini,soq,shine,qnot,qmin] = z0optimnbigeo(res,post,sycomore_data.shine_through_limit);
    if k > 11
	z0dinput.option.rtang  = 0.3 .* res(1) + 0.7 .* z0dinput.option.rtang;
	z0dinput.option.rtang2 = 0.3 .* res(2) + 0.7 .* z0dinput.option.rtang2;
	z0dinput.option.zext   = 0.3 .* res(3) + 0.7 .* z0dinput.option.zext;
	z0dinput.option.zext2  = 0.3 .* res(4) + 0.7 .* z0dinput.option.zext2;
    else
	z0dinput.option.rtang  = res(1);
	z0dinput.option.rtang2 = res(2);
	z0dinput.option.zext   = res(3);
	z0dinput.option.zext2  = res(4);   
    end
    fprintf('objective (to be minimized): %g\n',objective);
    fprintf('non inductive fraction: %g\n',fini);
    fprintf('current alignement: %g\n',ialign);
    fprintf('qnot = %g & qmin = %g\n',qnot,qmin);
    fprintf('shine through (and first orbit losses) & maximum allowed shine through (MW): %g <? %g\n',shine ./ 1e6 ,sycomore_data.shine_through_limit ./ 1e6 );
    fprintf('confinement criterium (<|s|/q>): %g\n',soq);
    fprintf('NBI1: Rtang = %g m & vertical shift = %g (normalized)\n',res(1),res(3));
    fprintf('NBI2: Rtang = %g m & vertical shift = %g (normalized)\n',res(2),res(4));
    if (sycomore_data.f_ni >= 1) && (zs.modeh(indHtest)  >= 1)
      if k > 11
	  z0dinput.cons.pnbi = (0.3 .* res(5) + 0.7) .* real(z0dinput.cons.pnbi) + sqrt(-1) .* (0.3 .* res(6) + 0.7) .* imag(z0dinput.cons.pnbi);
      else
	  z0dinput.cons.pnbi = res(5) .* real(z0dinput.cons.pnbi) + sqrt(-1) .* res(6) .* imag(z0dinput.cons.pnbi);
      end
      fprintf('NBI power multiplied by: NBI_1 = %g & NBI_2 = %g\n', ...
            max(real(z0dinput.cons.pnbi)) ./ max(real(pnbi_mem)), ...
            max(imag(z0dinput.cons.pnbi)) ./ max(imag(pnbi_mem)));	  
    end
    pause(3);
end

% faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value.
% will be implemented later using a optimizer 
% loop on non inductive current
% Ialign must included in the loop through zext and zext2
%  if sycomore_data.f_ni > 0
%    switch sycomore_data.f_ni  
%    case 1
%        % case fully non inductive, test on flux comsumption.
%    
%    otherwise
%        % case with ohmic residal current, test on flux comsumption.
%   
%    end
%  end

% print option data
disp('METIS parameters:')
z0dinput.option

% full computation
% turn off during testing
if ~ testing
  [zs,infovoid,profli] = zerod(z0dinput.option,z0dinput.cons,z0dinput.geo,z0dinput.exp0d);
end
post.z0dinput = z0dinput;
post.zerod    = zs;
post.profil0d =profli;


% calcul of peak power on divetor using result of 2 points model.
%      sycomore_data.flux_expansion = 3;       % divertor flux expansion 
%      sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
%      sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
post.zerod.peakdiv = compute_peak_div(post,sycomore_data.target_orientation,sycomore_data.S_factor,sycomore_data.flux_expansion);

% text information
tsnapshot = (1/4) .* t_flattop_plus  + (3/4) .* tend;
ind_fusmax = find(zs.temps >= tsnapshot,1); 
tsnapshot = zs.temps(ind_fusmax);
fprintf('snapshot point @ %g s\n',tsnapshot);
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = post.z0dinput.cons.picrh + post.z0dinput.cons.pecrh + real(post.z0dinput.cons.pnbi) + imag(post.z0dinput.cons.pnbi) + zs.pohm + zs.plh;
betanth = zs.wth .* (1.6.*pi./3) .* post.z0dinput.geo.a ./ zs.vp ./ post.z0dinput.geo.b0 ./ zs.ip;
fprintf('b0 = %g, R = %g, a = %g, K_area = %g, d_area = %g, epsi = %g\n', post.z0dinput.geo.b0(ind_fusmax),  post.z0dinput.geo.R(ind_fusmax), ...
        post.z0dinput.geo.a(ind_fusmax),post.z0dinput.geo.K(ind_fusmax),post.z0dinput.geo.d(ind_fusmax),post.z0dinput.geo.R(ind_fusmax)./post.z0dinput.geo.a(ind_fusmax));
fprintf('Ip = %g, qa = %g, q95 = %g, q0 = %g, qmin = %g, li = %g\n', ...
        zs.ip(ind_fusmax)./1e6,zs.qa(ind_fusmax),zs.q95(ind_fusmax), ...
        zs.q0(ind_fusmax),zs.qmin(ind_fusmax),zs.li(ind_fusmax));
fprintf('nbar/ngr = %g,<nD> = %g, <nT> = %g, <nHe> =%g, <nimp> =%g, <nW> =%g\n',zs.nbar(ind_fusmax) ./1e20./  (zs.ip(ind_fusmax) / 1e6) .*  ...
        (pi.* post.z0dinput.geo.a(ind_fusmax) .^ 2),zs.nDm(ind_fusmax)./1e19,zs.nTm(ind_fusmax)./1e19,zs.nhem(ind_fusmax)./1e19, ...
        zs.nimpm(ind_fusmax)./1e19,zs.nwm(ind_fusmax)./1e19);
fprintf('taue = %g, tauhe = %g , HH = %g, Wth = %g, W = %g, Wfast = %g\n', ...
        zs.taue(ind_fusmax),zs.tauhe(ind_fusmax),zs.taue(ind_fusmax)./zs.tauh(ind_fusmax).*zs.hitb(ind_fusmax),  ...
	zs.wth(ind_fusmax)./1e6,zs.w(ind_fusmax)./1e6,(zs.w(ind_fusmax) -  zs.wth(ind_fusmax)) ./1e6);
fprintf('Plh = %g, Picrh = %g, Pecrh = %g, Pnbi = %g, Pline = %g, Pbrem = %g, Pcyclo = %g\n',zs.plh(ind_fusmax)/1e6, ...
        zs.picrh(ind_fusmax)/1e6,zs.pecrh(ind_fusmax)/1e6,real(zs.pnbi(ind_fusmax)/1e6) + imag(zs.pnbi(ind_fusmax)/1e6),(zs.prad(ind_fusmax) + zs.pioniz(ind_fusmax) + zs.pradsol(ind_fusmax)) /1e6,zs.pbrem(ind_fusmax)/1e6, zs.pcyclo(ind_fusmax) ./1e6);
fprintf('<ne> = %g, ne0 = %g, ne0/<ne> = %g, <Te> = %g , Te0 = %g, Te0/<Te> = %g,  Ti/Te = %g\n',zs.nem(ind_fusmax)/1e19,zs.ne0(ind_fusmax)/1e19, ...
         zs.ne0(ind_fusmax) ./ zs.nem(ind_fusmax) ,zs.tem(ind_fusmax)/1e3,zs.te0(ind_fusmax)/1e3,zs.te0(ind_fusmax) ./ zs.tem(ind_fusmax),zs.tite(ind_fusmax));
fprintf('Palpha = %g, Pfus = %g, Q = %g, betan(total) = %g , betan(th) = %g, frad = %g  \n',zs.pfus(ind_fusmax)/1e6,zs.pfus(ind_fusmax).* rfan./1e6, ...
         rfan .* zs.pfus(ind_fusmax) ./ padd(ind_fusmax),zs.betan(ind_fusmax).*100,betanth(ind_fusmax).*100, ...
	 (zs.prad(ind_fusmax) + zs.pbrem(ind_fusmax) + zs.pcyclo(ind_fusmax)) ./ zs.pin(ind_fusmax));
fprintf('Iboot/Ipar = %g, Icd/Ipar = %g, Vloop = %g , Zeff(line) = %g, li =%g\n',zs.iboot(ind_fusmax) ./zs.ipar(ind_fusmax),zs.icd(ind_fusmax) ./ zs.ipar(ind_fusmax), ...
        zs.vloop(ind_fusmax),zs.zeff(ind_fusmax)./zs.zmszl(ind_fusmax),zs.li(ind_fusmax));
fprintf('Vol = %g, Section = %g, Sext = %g , Llcms = %g\n',zs.vp(ind_fusmax),zs.sp(ind_fusmax),zs.sext(ind_fusmax),zs.peri(ind_fusmax));
fprintf('frHe = %g, frImp = %g @ Z= %d  &  %g @ Z = %d, <Zeff> = %g, fr_W = %g \n', ...
         zs.nhem(ind_fusmax)./zs.nem(ind_fusmax),zs.nimpm(ind_fusmax)./zs.nem(ind_fusmax),z0dinput.option.zimp, ...
	 zs.nimpm(ind_fusmax)./zs.nem(ind_fusmax).* z0dinput.option.rimp,z0dinput.option.zmax, ...
	 (zs.n1m(ind_fusmax) + 4 .* zs.nhem(ind_fusmax) +zs.nimpm(ind_fusmax) .* z0dinput.option.zimp .^ 2 + ...
	 zs.nimpm(ind_fusmax) .* z0dinput.option.zmax  .^ 2 .* z0dinput.option.rimp) ./ zs.nem(ind_fusmax),  zs.nwm(ind_fusmax)./zs.nem(ind_fusmax));

% maximum power used for each heating sources:
fprintf('max(P_nbi)  = %g MW\n',(max(real(post.z0dinput.cons.pnbi)) + max(imag(post.z0dinput.cons.pnbi))) ./ 1e6)
fprintf('max(P_ecrh) = %g MW\n',max(post.z0dinput.cons.pecrh + post.z0dinput.cons.plh) ./ 1e6);


% flattop duration computation
post.flux_data = poloidal_flux(post,sycomore_data.couplage);
% adjustement of duration taking into account rampdwon.
% to be done

% ual output
% test if UAL is connected to METIS
if ~isappdata(0,'UAL_exist')
	setappdata(0,'UAL_exist',0);
	try
		rep = javaclasspath('-all');
		for l=1:length(rep)
			sl = rep{l}; 
			ind = findstr(sl,'ualjava.jar');
			if ~isempty(ind)
			    setappdata(0,'UAL_exist',1);			    
			end
		end
	end
end
if getappdata(0,'UAL_exist') == 1
    % preparing occurence structure
    if ~isempty(sycomore_data.tokamak) && ~isempty(sycomore_data.user) && ~isempty(sycomore_data.user) && ~isempty(sycomore_data.dataversion)
	occurrence.tokamak     = sycomore_data.tokamak;
	occurrence.user        = sycomore_data.user;
	occurrence.dataversion = sycomore_data.dataversion;
	occurrence.occurrence  = sycomore_data.occurrence;   
    else
      occurrence = [];
    end
    % writing data to UAL
    error_flag = metis4itm(sycomore_data.shot,sycomore_data.run,occurrence,post);
end

% metis file writing
metis_save(sprintf('%s_METIS_simulation',root_name),post)

% automatic graphs
if ~isempty(root_name)
    couplage = sycomore_data.couplage; 
    NOINTER  = 1;
    z0drapport_png(sprintf('%s_plot',root_name),1);
end
% end of print on standard ouput ?
% flush diary file
diary off 
% the reccord continue just the end of actor execution
diary on 

% computation of separatrix moment
function [R,a,z0,K,d] = sepamoment(option)
% calcul des moments
sepa  = z0dsepanew2(1,option);
% centre pour angle d'integration
rc = mean(sepa.R,2);
zc = mean(sepa.Z,2);
vc = ones(1,size(sepa.R,2));
uc = unwrap(angle((sepa.R-rc*vc) + sqrt(-1) .* (sepa.Z  -zc*vc)));
uc    = uc .* (uc >0) + (uc + 2*pi) .* (uc<= 0);
uc(:,1)   = uc(:,end) + 2 .* pi;
xu    = linspace(0,1,length(vc));
%dudx  = pdederive(xu,uc,2,2,2,1);
%dudx(:,1) = (dudx(:,1) +dudx(:,end)) ./ 2;
%dudx(:,end) = dudx(:,1);
dRdx  = pdederive(xu,sepa.R,2,2,2,1);
dZdx  = pdederive(xu,sepa.Z,2,2,2,1);
% calcul de R0 et Z0
maskrmax  = (sepa.R == (max(sepa.R,[],2) * vc));
% recalcul des parametres sur le vecteur final
rmin  = min(sepa.R,[],2);
rmax  = max(sepa.R,[],2);
a = 0.5 .* (rmax - rmin);
R = 0.5 .* (rmax + rmin);
zmin  = min(sepa.Z,[],2);
zmax  = max(sepa.Z,[],2);
z0   = (zmax + zmin + sum(sepa.Z .* maskrmax,2) ./ sum(maskrmax,2)) ./ 3;
K     = (abs(trapz(xu,sepa.Z .*  dRdx,2) ./ pi ./ a) + (zmax - zmin)) ./ 3 ./ a;
rzmax = R;
rzmin = R;
for k = 1:size(sepa.Z,1)
	rzmax(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmax(k))));
	rzmin(k) = sepa.R(k,min(find(sepa.Z(k,:) == zmin(k))));
end
uu   =  angle(rzmax - R + sqrt(-1) .* (zmax - z0));
ul   =  angle(rzmin - R + sqrt(-1) .* (zmin - z0));
tu   =  abs((acos((rzmax - R) ./ a) - acos(cos(uu))) ./ sin(uu));
tl   =  abs((acos((rzmin - R) ./ a) - acos(cos(ul))) ./ sin(ul));
tm   =  (tl + tu) ./ 2;
d    =   abs(rzmax + rzmin -  2 .* R) ./ 2 ./ a;
d    =  0.6 .* d + 0.4  .* sin(tm);



function flux_data_output_ = poloidal_flux(post,couplage)

% z0plotflux.m
NOPLUTFLUX_HELIOS = 1;
z0plotflux;

% create data structure
liste_of_fields_ = who;
for kkk_ = 1:length(liste_of_fields_)
    flux_data_output_.(liste_of_fields_{kkk_}) = eval(liste_of_fields_{kkk_});
end

function peakdiv = compute_peak_div(post,angle,S,flux_exp)
% reference  : Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER, 
% T. Eich et al, Nucl. Fusion 53 (2013) 093031 (7pp); doi:10.1088/0029-5515/53/9/093031
% plot des donnees du modele a 2 points
zs      = post.zerod;
option  = post.z0dinput.option; 
geo     = post.z0dinput.geo; 
cons    = post.z0dinput.cons; 
profli  = post.profil0d;

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
phys.pam3        =   (4.41e-4 .* phys.avo);    % conversion d'un nombre de particules en en Pa.m^3

% compatibilite
if ~isfield(option,'yield_model')
 	option.yield_model      = 'Javev';
end

% puissance conduite a la separatrice
pl        = max(zs.pin ./ 100,zs.pin - zs.prad - zs.pbrem - zs.pcyclo - zs.pioniz - zs.pradsol);
% fraction perdue en volume dans la sol par rayonnement:
%fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
switch option.sol_rad
case 'coupled'
	fesol   = max(0,min(1, (zs.pradsol + max(0,1 - option.fprad) .* zs.prad) ./ max(1,pl)));
        pradsol = zs.pradsol + max(0,1 - option.fprad) .* zs.prad;
otherwise
	fesol   = max(0,min(1, zs.pradsol ./ max(1,pl)));
        pradsol = zs.pradsol;
end

% these E. Tsitrone
lclim = pi .* geo.R .* zs.qa;
lcpol = pi .* geo.R;
lcx = sqrt(zs.peri .^ 2  + (pi .* geo.R .* option.lcx .* zs.qa) .^ 2);  
switch option.configuration
case 0
	lc = lcpol;
case 1
	lc = lclim;
case 2
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lcpol;
case 3
	lc  = zs.xpoint .* lcx + (~zs.xpoint) .* lclim;
otherwise
	lc  = lcx;
end
%lc  = zs.modeh .* lcx + (~zs.modeh) .* lclim;

if isfield(zs,'dsol')
	dsol = zs.dsol;
elseif option.sol_lscale  == 0
    dsol        = geo.a ./ 100;
elseif option.sol_lscale  > 0
    dsol        = geo.a .* option.sol_lscale;
else
    dsol        = - geo.R .* option.sol_lscale;
end


% flux //  (formula 5.64 et 5.75 Stangeby)
x      = profli.xli;
ve     =  ones(size(x));
Raxea  = interp1(profli.temps,profli.Raxe,zs.temps,'pchip','extrap');
Fa     = interp1(profli.temps,profli.fdia,zs.temps,'pchip','extrap');
psi    = interp1(profli.temps,profli.psi,zs.temps,'pchip','extrap');
rmx    = interp1(profli.temps,profli.rmx,zs.temps,'pchip','extrap');
zeffp  = interp1(profli.temps,profli.zeff,zs.temps,'pchip','extrap');
nzp    = interp1(profli.temps,profli.nzp,zs.temps,'pchip','extrap');
nwp    = interp1(profli.temps,profli.nwp,zs.temps,'pchip','extrap');
nhep   = interp1(profli.temps,profli.nhep,zs.temps,'pchip','extrap');
n1p    = interp1(profli.temps,profli.n1p,zs.temps,'pchip','extrap');
nep    = interp1(profli.temps,profli.nep,zs.temps,'pchip','extrap');
Qe     = interp1(profli.temps,profli.qe(:,end),zs.temps,'pchip','extrap');
Qi     = interp1(profli.temps,profli.qi(:,end),zs.temps,'pchip','extrap');
rext         = Raxea + geo.a * x;
btor         = Fa ./ rext;
grho         = abs((rmx(:,end) * ve)./ max(eps,pdederive(x,rext,0,2,2,1)));
grho(:,1)    = grho(:,2);
bpol         = -pdederive(x,psi,0,2,2,1) ./ rext .* grho ./ (rmx(:,end) * ve);
ut           = atan(abs(bpol(:,end) ./  btor(:,end)));
% flux //  (formula 5.64 et 5.75 Stangeby)
switch option.sol_model
case '2_points';
	% rien
otherwise
	warning('the 2 points model is not used in this simulation','2 points model');
	option.sol_model = '2_points';
end
[tebord,nelim,telim,qpl_target,err,nb,indbad,fmom,qpl_rad_div,qpl_neutral_div,qpl_tot,pl,zeff_div,gamma,mach_target] = ...
    z0convergence_2points_dic(option,post.z0dinput.cons,post.z0dinput.geo,post.zerod,post.profil0d);
Asol_para = 2 .* pi  .* Raxea(:,end) .* dsol .* sin(ut);
% 50 % du rayonnement du divertor retourne sur les plaques.
pref = (qpl_target  + 0.5 .* qpl_rad_div) .* Asol_para;

maskx = ones(size(zs.temps));
maskx(zs.xpoint == 0) = NaN;

% estimation parametrique du la puissance deposee
[sb,qpdep,fnorm] = z0div_power_dep(angle,S,flux_exp,dsol,pref,option.fR_target .* geo.R,max(geo.a ./ 2));
peakdiv = max(qpdep,[],2);
peakdiv(zs.xpoint == 0)  = NaN;

