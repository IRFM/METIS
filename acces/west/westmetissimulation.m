% WESTMETISSIMULATION : prototype of scenarion generator for WEST
%-----------------------------------------------------------------------------------------------------------
% function Matlab 2012b: westmetissimulation.m -> westmetissimulation
%
% This function prototype of scenarion generator for WEST
%
% syntax:
%
%   testing :
%       z0dinput = westmetissimulation([],'',[]);
%
%   parameters declaration:
%       option = westmetissimulation(1);
%
%   computation :
%       [z0dinput,tsnapshot] = westmetissimulation(sepa_option,parameters_filename,reference_simulation);
%
% input:
%
%     sepa_option: etheir scenario label, CREATE reference or data structure for LCFS description (see below)
%                  labels in {1,2,3,4.1,4.2,5.1,5.2}
%
%     parameters_filename: name of METIS parametrs file or METIS paramters data structure (see metis4itm.m and zerod_param.m)
%                          leave empty to use default parameters
%
%     reference_simulation: Input data for scenario
%
% structure details for reference_simulation (with data example for dummu WEST scenario):
%
%      reference_simulation.gas = 2;               
%      reference_simulation.ip = 3.5e6;
%      reference_simulation.rb0 =  2.28.* 2.93 ;    
%      reference_simulation.ip_first = 0.4e6; 
%      reference_simulation.f_Greenwald = 0;      
%      reference_simulation.density     = 4.2e19;        
%      reference_simulation.edge_density_factor     = 1;        
%      reference_simulation.H_H = 1.4;                
%      reference_simulation.ITB = 'on'; 
%      reference_simulation.shot = 1;                 
%      reference_simulation.run  = 1;                 
%      reference_simulation.PLHCD  = 7e6;              
%      reference_simulation.PICRH  = 12e6;              
%      reference_simulation.PECCD    = 6e5;              
%      reference_simulation.PBREAK   = 6e5;              
%      reference_simulation.PRAMPUP  = 3e4;              
%      reference_simulation.Recycling = 0.7;       
%      reference_simulation.radiation  = 'Lz';               % model for line radiative power in core plasam (Lz = colling rate or Matthews)
%      reference_simulation.Zeff  = 1.2;               
%      reference_simulation.Cw  = 1e-4;               
%      reference_simulation.SOL_model  = '2_points';         % SOL model: scaling or 2_points
%      reference_simulation.runaway    = 'on';         
%      reference_simulation.breakdown  = 'on';         
%      reference_simulation.duration         = 30;      
%      reference_simulation.f_dipdt_rampup   = 0.6e6;   
%      reference_simulation.f_dipdt_rampdown = 0.6e6;   
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
%      tsnapshot:     time slice selected for performances evaluation
%
% function writed by J-F Artaud
% CVS version (created 2014/07/15)
%-----------------------------------------------------------------------
%
function [z0dinput,tsnapshot] = westmetissimulation(sepa_option,parameters_filename,reference_simulation)


% parameters declaration
if nargin < 1
  sepa_option = [];
end
if (nargin <= 1) && ~isstruct(sepa_option)
    rep = dir(fullfile(fileparts(which('westmetissimulation.m')),'ref_equilibrium','*.mat'));
    for kl =1:length(rep)
      name_list{kl}  = rep(kl).name;
    end

    valeur.LCFS = name_list{1};     
    type.LCFS   = 'string';
    borne.LCFS  = name_list;  
    defaut.LCFS = name_list{1};
    info.LCFS   = 'Plasma shape: LCFS reference from CEDRES++';    

    rep = dir(fullfile(fileparts(which('westmetissimulation.m')),'configuration','*.txt'));
    conf_list{1} = 'default tuning of the scenario generator';
    for kl =1:length(rep)
      conf_list{kl+1}  = rep(kl).name;
    end

    valeur.configuration = conf_list{1};     
    type.configuration   = 'string';
    borne.configuration  = conf_list;  
    defaut.configuration = conf_list{1};
    info.configuration   = 'METIS internal models: initial model parameter set;\nif = none, use default tuning of the scenario generator.';    

    valeur.gas = 2;   % gas species as in METIS (1=H, 2=D, 3=DT & 4=He)      
    type.gas   = 'integer';
    borne.gas     = {1,2,4};  
    defaut.gas    = 2;
    info.gas      = 'Plasma composition: main gas species: 1 -> H, 2 -> D, 4 -> He';
    
    valeur.ip     = 0.8;   % plasma current (MA)  
    type.ip       = 'float';
    borne.ip      = [0.4,1];  
    defaut.ip     = 3.5;
    info.ip       = 'plasma current (MA)';

    valeur.b0     = 3.7;   
    type.b0       = 'float';
    borne.b0      = [1.8,4.2];  
    defaut.b0     = 3.7;
    info.b0       = 'vacuum magnetic field @ R_0 (T)';

    valeur.ip_first = 0.4;       
    type.ip_first   = 'real';
    borne.ip_first  = [0.05 0.4];  
    defaut.ip_first = 0.4;
    info.ip_first   = 'plasma current at time when PCS take control of ramp-up rate (MA)';
    mode.ip_first   = 'advanced';
    
    valeur.f_Greenwald     = 0;   
    type.f_Greenwald       = 'float';
    borne.f_Greenwald      = [0,1.4];  
    defaut.f_Greenwaldd    = 0;
    info.f_Greenwald       = 'Greenwald fraction. If set to 0, generator uses density parameter';
    
    valeur.density     = 4.2;   
    type.density       = 'float';
    borne.density      = [1,15];  
    defaut.density     = 4.2;
    info.density       = 'line averaged electron density during flattop (10^{19} m^{-3})';
    
    valeur.edge_density_factor     = 1;   
    type.edge_density_factor       = 'float';
    borne.edge_density_factor      = [0.1,10];  
    defaut.edge_density_factor     = 1;
    info.edge_density_factor       = 'edge density: multiplication factor applied to edge density scaling law:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar';
    mode.edge_density_factor   = 'advanced';
    
    valeur.H_H     = 0.8;   
    type.H_H       = 'float';
    borne.H_H      = [0.5,2];  
    defaut.H_H     = 0.8;
    info.H_H       = 'time confinement multiplication factor during H-mode phase';
    
    valeur.ITB = 'on';       
    type.ITB   = 'string';
    borne.ITB  = {'on','off'};  
    defaut.ITB = 'on';
    info.ITB   = 'allows or not ITB formation (on/off)';

    valeur.PICRH     = 0;   
    type.PICRH       = 'float';
    borne.PICRH      = [0 12];  
    defaut.PICRH     = 0;
    info.PICRH       = 'ICRH source: maximum power for PICRH during flattop (MW)';
    
    valeur.PLHCD     = 0;   
    type.PLHCD       = 'float';
    borne.PLHCD      = [0 7];  
    defaut.LHCD      = 0;
    info.PLHCD       = 'LHCD source: maximum power for LHCD during flattop (MW)';
    
    valeur.PECCD     = 0;   
    type.PECCD       = 'float';
    borne.PECCD      = [0 0.6];  
    defaut.PECCD     = 0;
    info.PECCD       = 'maximum power for ECRH/ECCD during flattop (MW):\nbaseline maximum power is 0 MW';
    mode.PECCD       = 'advanced';

    valeur.PBREAK     = 0;   
    type.PBREAK       = 'float';
    borne.PBREAK      = [0 0.6];  
    defaut.PBREAK     = 0;
    info.PBREAK       = 'Plasma initiation: injected ECRH / ECCD power for assisted breakdown (MW);\nnote that only a small fraction of this power is absorbed:\nbaseline maximum power is 0 MW';
    mode.PBREAK       = 'advanced';

    valeur.PRAMPUP     = 0;   
    type.PRAMPUP       = 'float';
    borne.PRAMPUP      = [0 0.6];  
    defaut.PRAMPUP     = 0;
    info.PRAMPUP       = 'ECRH / ECCD assisted ramp-up:  power during ramp-up (MW)';

    valeur.Recycling     = 0.7;   
    type.Recycling       = 'float';
    borne.Recycling      = [0,1-1e-4];  
    defaut.Recycling     = 0.7;
    info.Recycling       = 'Recycling coefficient';
    
    valeur.radiation = 'Lz';       
    type.radiation   = 'string';
    borne.radiation  = {'Lz','Matthews'};  
    defaut.radiation = 'Lz';
    info.radiation   = 'Line radiation model: model for line radiative power in core plasma (Lz = cooling rate or Matthews = Matthews scaling law)';
    mode.radiation   = 'advanced';

    valeur.Cw            = 1e-4;   
    type.Cw              = 'float';
    borne.Cw             = [0,1e-3];  
    defaut.Cw            = 1e-4;
    info.Cw              = 'Tungsten concentration with respect to electron density';
    mode.Cw              = 'advanced';
    
    valeur.Zeff            = 1.2;   
    type.Zeff              = 'float';
    borne.Zeff             = [1.1,6];  
    defaut.Zeff            = 1.2;
    info.Zeff              = 'Effective charge during flattop';
    mode.Zeff              = 'advanced';
    
    valeur.SOL_model = 'scaling';       
    type.SOL_model   = 'string';
    borne.SOL_model  =  {'scaling','2_points'};  
    defaut.SOL_model = 'scaling';
    info.SOL_model   = 'SOL model: scaling law or 2 points model';
    mode.SOL_model   = 'advanced';

    valeur.runaway = 'off';       
    type.runaway   = 'string';
    borne.runaway  = {'on','off'};  
    defaut.runaway = 'off';
    info.runaway   = 'allows or not runaway electrons in the discharge (on/off)';
    mode.runaway   = 'advanced';

    valeur.breakdown = 'off';       
    type.breakdown   = 'string';
    borne.breakdown  = {'on','off'};  
    defaut.breakdown = 'off';
    info.breakdown   = 'turn on or off the breakdown model for plasma initiation';
    mode.breakdown   = 'advanced';

    valeur.duration     = 30;   
    type.duration       = 'float';
    borne.duration      = [10,1000];  
    defaut.duration     = 10;
    info.duration       = 'shot duration: time of end of flat-top (s)';
    
    valeur.f_dipdt_rampup     = 0.6;   
    type.f_dipdt_rampup       = 'float';
    borne.f_dipdt_rampup      = [0.1,1.5];  
    defaut.f_dipdt_rampup     = 0.6;
    info.f_dipdt_rampup       = 'Ramp-up rate (MA/s)';
      
    valeur.f_dipdt_rampdown     = 1;   
    type.f_dipdt_rampdown       = 'float';
    borne.f_dipdt_rampdown      = [0.3,3];  
    defaut.f_dipdt_rampdown     = 1;
    info.f_dipdt_rampdown       = 'Ramp-down rate (MA/s)';
      
    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

    z0dinput.valeur     = valeur;
    z0dinput.type       = type;
    z0dinput.borne      = borne;
    z0dinput.defaut     = defaut;
    z0dinput.info       = info;
    z0dinput.interface  = interface;
    z0dinput.mode  	= mode;
    

    z0dinput.description = sprintf('WEST scenario generator; references:\nWEST Research Plan');   % description (une ligne) de la fonction

    z0dinput.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    z0dinput.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    z0dinput.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

    return

end
% gestion of input number and contents
if nargin < 3
    reference_simulation = '';
end
if nargin < 2
    parameters_filename = [];
end
if nargin < 1
    sepa_option = [];
end

if (nargin == 1) && isstruct(sepa_option)
      west_param = sepa_option;
      clear sepa_option;
      sepa_option = west_param.LCFS;
else
      west_param = [];
end

% set parameter file if not already provided
if isempty(parameters_filename) && isfield(west_param,'configuration') && ~isempty(west_param.configuration)

 switch west_param.configuration
 case {'none','default tuning of the scenario generator'}
  %rien
 otherwise
      parameters_filename = fullfile(fileparts(which('westmetissimulation.m')),'configuration',west_param.configuration);
 end

end
% structure of parameters for METIS
%z0dinput = zerod_init(-2);
%option = z0dinput.option;
info = metis4imas;
option = info.valeur;
option.gaz = 2;
option.neasser = 1;
option.Recycling = 0.7;
option.natural = 1;
option.ane = 11;
option.fn0a = 1;
option.fn0a_div = 0.1;
%
option.scaling = 0;
option.dilution = 1;
option.tau_limitation = 'On';
option.l2hscaling = 3;
option.pl2h_mass_charge = 1;
option.modeh = 1;
option.hysteresis = 0;
option.configuration = 2;
option.l2hslope =  0.5;
option.usepped_scl = 1;
option.taurotmul = 0;
option.fintrinsic = 0.2;
option.xiioxie = -4.5;
option.kishape = 0;
option.xieorkie = 0;
option.omega_shape = 0;
option.fstiff = 1;
option.ploss_exp = 'max_power';
option.xiioxie_ped = 0;
option.hmore_pped  = 2;
option.fpl2h_lim   = 2;
option.ki_expo     = 2;
option.plhthr      = 'P_LCFS';
option.grad_ped    = 3;
option.ode_pped    = 1;
option.adiabatic   = 1;

%
option.qdds = 1;
option.kidds = 3;
option.sitb = 2;
option.itb_sensitivity = 1;
option.itb_slope_max = 2;
option.smhd = 100;
option.tmhd = 0;
%
option.runaway = 5;
option.modeboot = 2;
%
option.li = 1;
option.breakdown = - 10;
option.berror = 1e-3;
option.L_eddy = 4.76e-4;
option.R_eddy = 4.1e-5;
option.C_eddy = 1;
option.B_eddy = 0.1;
option.I_eddy = 0;
option.p_prefill = 7e-03;
option.VV_volume = 58;
%
option.zeff = 0;
option.faccu = 0;
option.heat_acc = 0;
option.fne_acc = 0;
option.zmax = 8;
option.zimp = 8;
option.rimp = 1;
option.density_model ='minconv';
%
option.frad = 1;
option.matthews = -1;
option.fte_edge = 1;
option.gaunt = 1;
option.noncoronal = 2;
option.z_prad = 'Stangeby';
%
option.sol_lscale = 0;
option.eioniz     = 0;
option.fnesol     = 0;
option.sol_model  = '2_points';
option.lcx = 1;
option.fcond = -1;
option.fmom = 0;
option.lambda_scale = 3;
option.sol_rad = 'decoupled';
option.fzmax_div = -1;
option.W_effect = 1;
option.cw_factor = 0;
option.cw_offset = 5e-5;
option.factor_scale = 1;
option.fpower = 0.6000;
option.fR_target = 2.7 / 2.94;
option.mach_corr = 1;

%
option.angle_ece = 90;
option.synergie  = 1;
option.sens      = 1;
option.eccdmul   = 1;
%
option.angle_nbi = 90;
option.rtang     = 2.85;
option.zext      = 0.3;
option.einj      = 500000;
option.nbicdmul  = 1;
option.nb_nbi    = 2;
option.e_shielding = 'Honda-Sauter';
option.drs1        = 0;
option.dzs1        = 0;
%
option.angle_nbi2 = 0;
option.rtang2 = 2.85;
option.zext2  = 0.1;
option.einj2  = 85000;
option.nbicdmul2 = 1;
option.drs2   = 0;
option.dzs2   = 0;
%
option.lhmode = 0; 
option.etalh  = 0.8;  
option.wlh = 0.58;
option.xlh = 0.3;
option.dlh = 0.2;
option.angle_ece2 = 90;
option.npar0 = 2;
option.npar_neg = -4;

% used as third injector
option.fwcd = 0;
option.mino = 'H';
option.cmin = 0.15;
option.nphi = 30;
option.freq = 55.5;
option.icrh_width = 0.7;
option.fact_mino  = 0;
option.orbit_width  = 1;
option.icrh_model = 'Dumont-Vu';
%
option.equi_ppar = 3;
option.signe = 1;
option.cronos_regul= 2;
option.available_flux =  9.8  - (-7.82);  %WB
option.machine = 'WEST';
option.evolution = 0;
%

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

% default value for IMAS
option.COCOS  = 11;
option.signe = 1;
option.orientation = -1;

% overwriting
if ischar(parameters_filename)
    if ~isempty(parameters_filename)
        option = z0doverwriteparam(parameters_filename,option);
        option.reference_parameters  = '';
    end
elseif isstruct(parameters_filename)
    option = parameters_filename;
end

if isempty(reference_simulation)
    reference_simulation.gas = 2;               
    reference_simulation.ip = 1e6;
    reference_simulation.rb0 =  2.55 .* 3;    
    reference_simulation.ip_first = 0.4e6; 
    reference_simulation.f_Greenwald = 0.4;      
    reference_simulation.density     = 4e19;        
    reference_simulation.edge_density_factor     = 1;        
    reference_simulation.H_H = 0.8;                
    reference_simulation.ITB = 'on'; 
    reference_simulation.shot = 1;                 
    reference_simulation.run  = 1;                 
    reference_simulation.PLHCD  = 7e6;              
    reference_simulation.PICRH  = 12e6;              
    reference_simulation.PECCD    = 6e5;              
    reference_simulation.PBREAK   = 6e5;              
    reference_simulation.PRAMPUP  = 6e5;              
    reference_simulation.Recycling = 0.7;       
    reference_simulation.radiation  = 'Lz';               % model for line radiative power in core plasam (Lz = colling rate or Matthews)
    reference_simulation.Zeff       =  1.2;       
    reference_simulation.Cw         = 1e-4;       
    reference_simulation.SOL_model  = '2_points';         % SOL model: scaling or 2_points
    reference_simulation.runaway    = 'off';         
    reference_simulation.breakdown  = 'off';         
    reference_simulation.duration         = 30;      
    reference_simulation.f_dipdt_rampup   = 0.6e6;   
    reference_simulation.f_dipdt_rampdown = 0.6e6;   
end

% input data
external_lcfs = 0;
Kref = (1.687 + 2.001) / 2;
dref = (0.466 + 0.568) / 2;
fKup = 1.687 ./ Kref;
fKdo = 2.001 ./ Kref;
fdup = 0.46 ./ dref;
fddo = 0.568 ./ dref;
% LCFS data
if isempty(sepa_option)
    rep = dir(fullfile(fileparts(which('westmetissimulation.m')),'ref_equilibrium','*.mat'));
    sepa_option = rep(1).name;
end
if isstruct(sepa_option)
    % nothing to be done
    % other information
    sepa_option.nbp       = 201;                 % number of points for the separatrix (depends on equilibrium module) [201]
    sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
    sepa_option.filename  = '';
    %sepa_option.za        = -(sepa_option.zxup - sepa_option.zxdo) .* sepa_option.a;       % altitude of the magnetic axis (m) [0.9]
    z0_before_xpoint      = -(sepa_option.zxup - sepa_option.zxdo) .* sepa_option.a;       % altitude of the magnetic axis (m) [0.9]
    % from separatrix
    [R,a,z0,K,d]   = sepamoment(sepa_option);
else
      % loading LCFS
      lcfs_data = load(fullfile(fileparts(which('westmetissimulation.m')),'ref_equilibrium',sepa_option));
      KH = sort(unique(convhull(lcfs_data.eqgeometry.boundary.r,lcfs_data.eqgeometry.boundary.z)));
      if (length(KH) ~= length(lcfs_data.eqgeometry.boundary.r))
	    index_full = 1:length(lcfs_data.eqgeometry.boundary.r);
	    Rsepa = lcfs_data.eqgeometry.boundary.r(KH);
	    Zsepa = lcfs_data.eqgeometry.boundary.z(KH);
	    lcfs_data.eqgeometry.boundary.r = interp1(KH,Rsepa,index_full,'linear');
	    lcfs_data.eqgeometry.boundary.z = interp1(KH,Zsepa,index_full,'linear');
	    indbad_lcfs = find(~isfinite(lcfs_data.eqgeometry.boundary.r) | ~isfinite( lcfs_data.eqgeometry.boundary.z));
	    if ~isempty(indbad_lcfs)
		lcfs_data.eqgeometry.boundary.r(indbad_lcfs) = [];
		lcfs_data.eqgeometry.boundary.z(indbad_lcfs) = [];
	    end
      end
      if (lcfs_data.eqgeometry.boundary.r(1) ~= lcfs_data.eqgeometry.boundary.r(end)) || ...
         (lcfs_data.eqgeometry.boundary.z(1) ~= lcfs_data.eqgeometry.boundary.z(end))
         lcfs_data.eqgeometry.boundary.r(end+1) = lcfs_data.eqgeometry.boundary.r(1);
         lcfs_data.eqgeometry.boundary.z(end+1) = lcfs_data.eqgeometry.boundary.z(1);
      end
      control =sepa_moments(lcfs_data.eqgeometry.boundary.r,lcfs_data.eqgeometry.boundary.z);
      if ~isempty(control.x_point)
	xup = 0;
	xdown = 0;
	for k = 1:length(control.x_point)
	  if control.x_point{k}.z < 0
	    xdown = 1;
	  elseif control.x_point{k}.z > 0
	    xuup = 1;
	  end
	end
      end
      clear sepa_option;
      sepa_option.rxup      = control.du;     % upper triangularity (minor radius unit)
      sepa_option.zxup      = control.Ku;    % upper altitude X point (minor radius unit)
      if xup
	  sepa_option.apup      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	  sepa_option.amup      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
      else
	sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
      end
      sepa_option.ra        = control.R0;       % major radius R0 (m) [6.2]
      sepa_option.za        = control.z0_geo;       % altitude of the magnetic axis (m) [0.9]
      sepa_option.a         = control.a;         % minor radius (m) [2]
      sepa_option.rxdo      = control.dl;     % lower triangularity (minor radius unit)
      sepa_option.zxdo      = control.Kl;       % lower altitude X point (minor radius unit)
      if xdown
	sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
      else
 	sepa_option.apdo      = 0;   % lower separatrix angle (R,X)  (LFS, degrees)
	sepa_option.amdo      = 0;   % lower separatrix angle (-R,X)  (HFS, degrees)
      end
      sepa_option.delta     = 0.73;      % magnetic field at R0
      sepa_option.b0        = 11;
      sepa_option.nbp       = 201;       % number of points for the separatrix (depends on equilibrium module) [201]
      sepa_option.mode       = 'elliptical';       % number of points for the separatrix (depends on equilibrium module) [201]
      

      % create file with LCFS data
      sepa_option.filename = sprintf('%s.mat',tempname);
      zassignin('base','sepa_option',sepa_option);
      R = lcfs_data.eqgeometry.boundary.r;
      Z = lcfs_data.eqgeometry.boundary.z;
      save(sepa_option.filename,'R','Z');
      clear R Z
      
      % data for generator
      R  = control.R0;
      a  = control.a;
      z0 = control.z0_geo;
      K  = control.K;
      d  = min(control.dl,control.du);

end

if ~isempty(west_param)

    reference_simulation.gas = west_param.gas;               
    reference_simulation.ip  = west_param.ip .* 1e6;
    reference_simulation.rb0 = west_param.b0 .* sepa_option.ra;    
    reference_simulation.ip_first = west_param.ip_first .* 1e6; 
    reference_simulation.f_Greenwald = west_param.f_Greenwald;   
    reference_simulation.density     = west_param.density .* 1e19;        
    reference_simulation.edge_density_factor     = west_param.edge_density_factor;        
    reference_simulation.H_H = west_param.H_H;                
    reference_simulation.ITB = west_param.ITB; 
    try
	reference_simulation.shot = evalin('base','post.z0dinput.shot');  
    catch
	reference_simulation.shot = 1;
    end 
    try
	reference_simulation.run  = evalin('base','post.z0dinput.run') + 1;  
    catch
	reference_simulation.run  = 1; 
    end
    reference_simulation.PLHCD = west_param.PLHCD .* 1e6;            
    reference_simulation.PICRH = west_param.PICRH .* 1e6;              
    reference_simulation.PECCD    = west_param.PECCD .* 1e6;              
    reference_simulation.PBREAK   = west_param.PBREAK .* 1e6;              
    reference_simulation.PRAMPUP  = west_param.PRAMPUP .* 1e6;              
    reference_simulation.Recycling = west_param.Recycling;       
    reference_simulation.radiation  = west_param.radiation;               % model for line radiative power in core plasam (Lz = colling rate or Matthews)
    reference_simulation.Zeff  = west_param.Zeff;              
    reference_simulation.Cw  = west_param.Cw;               
    reference_simulation.SOL_model  = west_param.SOL_model;         % SOL model: scaling or 2_points
    reference_simulation.runaway    = west_param.runaway;         
    reference_simulation.breakdown    = west_param.breakdown;         
    reference_simulation.duration   = west_param.duration;     
    reference_simulation.f_dipdt_rampup     = west_param.f_dipdt_rampup   .* 1e6;   
    reference_simulation.f_dipdt_rampdown   = west_param.f_dipdt_rampdown .* 1e6;   
end


% open wrinting on standard output
tokamak = 'WEST';
fprintf('Starting METIS simumlation for %s\n',tokamak);
root_name = '';


% end rampdown vertical position
%zlow = sepa_option.za - K .* a + a / 2;
zlow = z0 - K .* a + a / 2;
z0_before_xpoint = 0;
% magnetic rigidity
rb0 = reference_simulation.rb0;
b0  = rb0 ./ R;
% flattop plasma current
ip   = reference_simulation.ip; %15e6;

% will scale with volume on ITER design
%  option for ECRH only or mixed ECRH/NBI
% the maximum power will also scale with L2H threshold
% ECRH flattop
pecrh1   = reference_simulation.PECCD;
% ecrh power are lineraly interpolated durant ramp-up and ramp-down
% preheat  (from start to xpoint, limiter mode)
p0_ecrh1 = 0e6;
% ramp up (from xpoint formation to H mode transition)
p1_ecrh1 = 0e6;
% end rampup (from hmode transition to end of ramp-up)
p2_ecrh1 = reference_simulation.PECCD;
% full power, up to ignition
% rampdown
p3_ecrh1 = reference_simulation.PECCD;
%
pecrh2     = reference_simulation.PBREAK; 
p0_ecrh2   = reference_simulation.PBREAK;

% NBI :
plhcd  =  reference_simulation.PLHCD;
picrh  =  reference_simulation.PICRH;

% heat source default parameters
xece  = 0;

% physics model assumption
if reference_simulation.f_Greenwald > 0
    % Greenwald  fraction
    fnegr =  reference_simulation.f_Greenwald;
    % conversion (10^20 m^-3)
    nbar  = fnegr .* 1e20 .* (ip / 1e6) ./ (pi.* a .^ 2);
else
    nbar  = reference_simulation.density;
end
% enhancement factor for energy content on flat top
hmore = reference_simulation.H_H;

% plasma composition
iso   = 0;
ftnbi = 0;
switch reference_simulation.gas
    case 4
        zeff  = max(2.1,reference_simulation.Zeff);
    otherwise
        zeff  = max(1.1,reference_simulation.Zeff);
end
%%%%%%%%%%%%%%%%%%%%%%%%
% start data estimation
%%%%%%%%%%%%%%%%%%%%%%%%
% constant
mu0  = (4*pi.*1e-7);

% from reference:
% H Urano, Fusion engineering and design  100 (2015)  345-356
dipdt_ref   = reference_simulation.f_dipdt_rampup; % A/s
dipdt       = reference_simulation.f_dipdt_rampup;
dipdt_down  = reference_simulation.f_dipdt_rampdown;
fprintf('dIp/dt_{rampup} = %g (kA/s)\n',dipdt/1e3);
fprintf('dIp/dt_{rampdown, end of plasma} = %g (kA/s)\n',dipdt_down/1e3);

% dynamical parameters
% current after breakdown @ 60 ms
ip_ini = 40e3;
ip_ms   = max(ip_ini,reference_simulation.ip_first);
% ip @ x-point formation
ip_xpoint = (ip + reference_simulation.ip_first) ./ 2;
% flattop current
ip_flattop = ip;
% plasam current for the back transition to L-mode (reduced magnetic energy by a factor 2)
ip_dn = ip ./ sqrt(2);

% used of scaling on density for minimal L2H power threshold. reference: F. Ryter N.F. 2014
nbar_l2h_min    = min(nbar,0.7e19 .* (ip ./ 1e6) .^ 0.34 .* a .^ -0.95 .* b0 .^ 0.62 .* (R ./ a) .^ 0.4);
% denstity at the start of flatop
nbar_flattop    = min(0.5 .* (nbar + nbar_l2h_min),1e20 .* (ip / 1e6) ./ (pi.* a .^ 2));
% density after breakdown @ 60ms
R_ini = 2.1;
a_ini = R_ini - 1.84;
q95_ini = 5 .* a_ini .^ 2 .* b0 ./ (ip_ini./1e6) ./ R_ini .* (1.17 - 0.65 .* a_ini ./ R_ini) ./ (1 - (a_ini./R_ini) .^ 2 ) .^ 2;
negr_ini = 1;
nini  = min(negr_ini .* 1e20 .* (ip_ini / 1e6) ./ (pi.* a_ini .^ 2), ...
    1e20 .* (b0 ./ q95_ini ./ R_ini) .^ 0.6  .*  0.25);


% times of interrest
% start of the simulation
t_start = 0;
% time for breakdown
t_break = 60e-3;
% time for switch resistor and end of ECRH breakdown
t_ms    =  0.16;
q95_ms = 5 .* a_ini .^ 2 .* b0 ./ (ip_ms./1e6) ./ R_ini .* (1.17 - 0.65 .* a_ini ./ R_ini) ./ (1 - (a_ini./R_ini) .^ 2 ) .^ 2;
nms  = min(negr_ini .* 1e20 .* (ip_ms / 1e6) ./ (pi.* a_ini .^ 2), ...
    1e20 .* (b0 ./ q95_ms ./ R_ini) .^ 0.6  .*  0.25);
nms = max(nini,nms);

t_xpoint        = t_ms + (ip_xpoint - ip_ms) ./ reference_simulation.f_dipdt_rampup;
t_before_xpoint = t_xpoint - 0.1 /2;
t_after_xpoint  = t_xpoint + 0.1 /2;
% start of the flattop
t_flattop    = t_xpoint + (ip - ip_xpoint) ./ reference_simulation.f_dipdt_rampup;
% time for full density (just a estimation)
t_flattop_plus = t_flattop + 1;
        


% density at x-point transition
q95_x = 5 .* a .^ 2 .* b0 ./ (ip_xpoint ./ 1e6) ./ R .* (1 + K .^ 2 .*  ...
    (1 + 2 .* d .^ 2 - 1.2 .* d .^ 3) ./ 2) .* (1.17 - 0.65 .* a ./ R) ./  ...
    (1 - (a./R) .^ 2 ) .^ 2;
n_xpoint  = max(negr_ini .* 1e20 .* (ip_xpoint ./ 1e6) ./ (pi.* a .^ 2), ...
    1e20 .* (b0 ./ q95_x ./ R) .^ 0.6  .*  0.25);
n_xpoint = min(n_xpoint,nbar_l2h_min);
n_xpoint = max(nms,n_xpoint);
nbar_flattop    = max(n_xpoint,nbar_flattop);

% data for rampdown
% simulation duration
% end time (wee look for a pulse duration greater than 2 hours)
% duration is increased to be sure to capture the time when all poloial flux is consummed
tend = reference_simulation.duration;
% time for back transition in L-mode
% there is trade between flux consumption and l_i change
tdn          = tend + (ip_flattop - ip_dn) ./ dipdt;
% time to back transition to limiter
tdip         = tdn  +  (ip_dn - ip_xpoint) ./ dipdt_down;
% end of ramp-down simulation (as fast as possible => limter mode)
tv0          = tdip + (ip_xpoint - ip_ini) ./ (dipdt + dipdt_down);

disp('Key times:');
fprintf('Starting time of the simulation: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_start,1e3/1e6,1e17/1e19);
fprintf('Beakdown time : %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_break,ip_ini/1e6,nini./1e19);
fprintf('End resitor network use: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_ms,ip_ms/1e6,nms./1e19);
fprintf('X-point formation: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_xpoint,ip_xpoint/1e6,n_xpoint./1e19);
fprintf('Start of flat-top: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_flattop,ip_flattop/1e6,nbar_flattop./1e19);
fprintf('Full power: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',t_flattop_plus,ip_flattop/1e6,nbar./1e19);
fprintf('End of flat-top: %g s @ Ip = %g MA & n_bar = %g 10^{19} m^{-3}\n',tend,ip_flattop/1e6,nbar./1e19);
fprintf('H to L back transition: %g s @ Ip = %g MA & n_bar > %g 10^{19} m^{-3}\n',tdn,ip_dn/1e6,nini./1e19);
fprintf('back transition to limiter: %g s @ Ip = %g MA & n_bar > %g 10^{19} m^{-3}\n',tdip,ip_xpoint/1e6,nini./1e19);
fprintf('End of simulation/plasma  : %g s @ Ip = %g MA & n_bar > %g 10^{19} m^{-3}\n',tv0,ip_ini/1e6,1e17./1e19);

% ICRH
rres  = R  - 0.2 .* a;
bres  = b0 .* R ./ rres;
ag    = 1;
zg    = 1;
harm  = 1;
freq  = harm .* bres .* (95.5e6 .* zg ./ ag) ./ (2 .* pi .* 1e6); % MHz
if freq < ((52 + 55.5)/2)
  freq = 52;
else
  freq = 55.5;
end

% time slices for the computation
temps = t_start:1e-3:t_break;
temps = union(temps,t_break:5e-3:t_ms);
temps = union(temps,t_ms:1e-2:t_after_xpoint);
temps = union(temps,t_after_xpoint:0.1:t_flattop_plus);
temps = union(temps,t_flattop_plus:0.5:tend);
temps = union(temps,tend:0.1:tv0);
if temps(end) < tv0
    temps = union(temps,tv0);
end
temps = temps';
% pre-parametrised METIS input data
z0dinput = zerod_scalaire(temps,b0,R,a,K,d,ip,nbar,plhcd,picrh,pecrh1,0,zeff,xece,hmore,iso,ftnbi,0);
% load standard reference parameters
z0dinput.option = option;
% first wall
z0dinput.option.first_wall = fullfile(fileparts(which('westmetissimulation.m')),'west_wall.mat');

%% switch off ip control during ramp down
%%z0dinput.option.vref = 0;
%%z0dinput.option.tswitch = tv0;
z0dinput.option.vloop = 0;
% offset to prevent to early H mode transition after X point formation (in the error bar of scaling for L 2 H transition)
z0dinput.option.l2hmul = 2 + (p0_ecrh1 + p0_ecrh2 + ip_ini) ./ 1e6;
% ICRH is configured to be used as NBI like source
z0dinput.option.freq =freq;
% random shot number
z0dinput.option.shot = reference_simulation.shot;
% decoration
z0dinput.option.machine = tokamak;
% for backward compatibility
z0dinput.machine = z0dinput.option.machine;
z0dinput.shot = z0dinput.option.shot;

% main gas
z0dinput.option.gaz = reference_simulation.gas;
% edge density
z0dinput.option.nea_factor = reference_simulation.edge_density_factor;        
% recycling
z0dinput.option.Recycling = reference_simulation.Recycling;        
% model for line radiative power in core plasam (Lz = colling rate or Matthews)
switch reference_simulation.radiation 
case 'Matthews'
  z0dinput.option.matthews  =  1;
otherwise
  z0dinput.option.matthews  =  -1;
end
z0dinput.option.cw_offset = reference_simulation.Cw;        


 % SOL model: scaling or 2_points
z0dinput.option.sol_model = reference_simulation.SOL_model;        

% with or without ITB
switch upper(reference_simulation.ITB)
case 'ON'
  z0dinput.option.sitb = 2;
  z0dinput.option.hmore_pped  = 0;
otherwise
  z0dinput.option.sitb = 0;
  z0dinput.option.hmore_pped  = 2;
end
% runaway model
switch upper(reference_simulation.runaway)
case 'ON'
  z0dinput.option.runaway = 5;
otherwise
  z0dinput.option.runaway = 0;
end
% breakdown model
switch upper(reference_simulation.breakdown)
case 'ON'
  z0dinput.option.berror = 1e-4;
otherwise
  z0dinput.option.berror = 0;
end

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

% first time
tv(end+1)      = t_start;
switch upper(reference_simulation.breakdown)
case 'ON'
  ipv(end+1)     = 1;
otherwise
  ipv(end+1)     = 1e3;
end
nv(end+1)      = 1e17;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PBREAK;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = R_ini ./ R;
av(end+1)      = a_ini ./ a;
Kv(end+1)      = 0;
dv(end+1)      = 0;
z0v(end+1)     = 0;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

% first time after breakdown
tv(end+1)      = t_break;
ipv(end+1)     = ip_ini;
nv(end+1)      = nini;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PBREAK;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = R_ini ./ R;
av(end+1)      = a_ini ./ a;
Kv(end+1)      = 0;
dv(end+1)      = 0;
z0v(end+1)     = 0;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

% MS time
tv(end+1)      = t_ms;
ipv(end+1)     = ip_ms;
nv(end+1)      = nms;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PBREAK;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = R_ini ./ R;
av(end+1)      = a_ini ./ a;
Kv(end+1)      = 0;
dv(end+1)      = 0;
z0v(end+1)     = (z0 + z0_before_xpoint) ./ 2;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

% time before X point formation
tv(end+1)      = t_before_xpoint;
ipv(end+1)     = ip_xpoint;
nv(end+1)      = n_xpoint;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PRAMPUP;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = (z0 + z0_before_xpoint) ./ 2;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;


% time for X point formation
tv(end+1)      = t_xpoint;
ipv(end+1)     = ip_xpoint;
nv(end+1)      = n_xpoint;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PRAMPUP;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = 3/4 .* z0 + z0_before_xpoint .* 1./ 4;
isov(end+1)    = 1;
xecev(end+1)   = 0;
hmore(end+1)   = 1;

% time after X point formation
tv(end+1)      = t_after_xpoint;
ipv(end+1)     = ip_xpoint;
nv(end+1)      = n_xpoint;
picrhv(end+1)  = 0;
plhv(end+1)    = 0;
pecrhv(end+1)  = reference_simulation.PRAMPUP;
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

%  % time for full power and H-mode transition
%  tv(end+1)      = t_full_power;
%  ipv(end+1)     = ip_full_power;
%  nv(end+1)      = nbar_full_power;
%  picrhv(end+1)  = 1;
%  plhv(end+1)    = 0;
%  pecrhv(end+1)  = p2_ecrh1;
%  pnbi1v(end+1)  = 0;
%  pnbi2v(end+1)  = 1;
%  Rv(end+1)      = 1;
%  av(end+1)      = 1;
%  Kv(end+1)      = 1;
%  dv(end+1)      = 1;
%  z0v(end+1)     = z0;
%  isov(end+1)    = 1;
%  xecev(end+1)   = xece;
%  hmore(end+1)   = reference_simulation.H_H;

% start of flattop
tv(end+1)      = t_flattop;
ipv(end+1)     = ip_flattop;
nv(end+1)      = nbar_flattop;
picrhv(end+1)  = 0;
plhv(end+1)    = 1;
pecrhv(end+1)  = pecrh1;
pnbi1v(end+1)  = 0;
pnbi2v(end+1)  = 0;
Rv(end+1)      = 1;
av(end+1)      = 1;
Kv(end+1)      = 1;
dv(end+1)      = 1;
z0v(end+1)     = z0;
isov(end+1)    = 1;
xecev(end+1)   = xece;
hmore(end+1)   = reference_simulation.H_H;

% full density 
tv(end+1)      = t_flattop_plus;
ipv(end+1)     = ip_flattop;
nv(end+1)      = nbar;
picrhv(end+1)  = 1;
plhv(end+1)    = 1;
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
hmore(end+1)   = reference_simulation.H_H;

% end of flattop
tv(end+1)      =  tend;
ipv(end+1)     =  ip_flattop;
nv(end+1)      =  nbar;
picrhv(end+1)  =  1;
plhv(end+1)    =  1;
pecrhv(end+1)  =  pecrh1;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  1;
av(end+1)      =  1;
Kv(end+1)      =  1;
dv(end+1)      =  1;
z0v(end+1)     =  z0;
isov(end+1)    =  1;
xecev(end+1)   =  xece;
hmore(end+1)   =  reference_simulation.H_H;

% shutdown NBI (and ICRH if present) 50%
% during start of rampdown
tv(end+1)      =  tend + (tdn - tend) ./ 3;
ipv(end+1)     =  ip_flattop - (ip_flattop - ip_dn) ./ 3;
nv(end+1)      =  nbar_flattop;
picrhv(end+1)  =  0;
plhv(end+1)    =  0.5;
pecrhv(end+1)  =  p3_ecrh1;
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

% shutdown NBI
% during start of rampdown
tv(end+1)      =  tend + 2 .* (tdn - tend) ./ 3;
ipv(end+1)     =  ip_flattop - 2 .* (ip_flattop - ip_dn) ./ 3;
nv(end+1)      =  nbar_flattop;
picrhv(end+1)  =  0;
plhv(end+1)    =  0;
pecrhv(end+1)  =  p3_ecrh1;
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
nv(end+1)      =  1e17;
picrhv(end+1)  =  0;
plhv(end+1)    =  0;
pecrhv(end+1)  =  0;
pnbi1v(end+1)  =  0;
pnbi2v(end+1)  =  0;
Rv(end+1)      =  R_ini/R;
av(end+1)      =  0.3;
Kv(end+1)      =  0;
dv(end+1)      =  0;
z0v(end+1)     =  zlow;
isov(end+1)    =  0.5;
xecev(end+1)   =  0;
hmore(end+1)   =  1;

% scales
picrhv  = picrhv  .* picrh;
plhv    = plhv    .* plhcd;
pnbi2v  = pnbi2v  .* 0;
Rv      = Rv .* R;
av      = av .* a;
Kv      = Kv .* (K - 1) + 1;
dv      = dv .* d;

% time interpolation
z0dinput.cons.ip = pchip(tv,ipv,temps);
z0dinput.cons.picrh = zinterpnc(tv,picrhv,temps);
%z0dinput.cons.pnbi = zinterpnc(tv,pnbi1v,temps) + sqrt(-1) .* zinterpnc(tv,pnbi2v,temps);
z0dinput.cons.plh = zinterpnc(tv,plhv,temps);
% change for break down 500 ms ECRH
z0dinput.cons.pecrh = zinterpnc(tv,pecrhv,temps);
z0dinput.cons.pecrh(temps < 0.16) = reference_simulation.PBREAK;
%indft = find((temps > t_xpoint) & (temps < t_full_power));
%z0dinput.cons.pecrh(indft) = pchip(tv,pecrhv,temps(indft));
z0dinput.cons.iso(:) = 0;
z0dinput.geo.a   = pchip(tv,av,temps);
z0dinput.geo.R   = pchip(tv,Rv,temps);
z0dinput.geo.K   = pchip(tv,Kv,temps);
z0dinput.geo.d   = pchip(tv,dv,temps);
z0dinput.geo.z0  = pchip(tv,z0v,temps);
z0dinput.cons.nbar = pchip(tv,nv,temps);
z0dinput.cons.xece = pchip(tv,xecev,temps);
z0dinput.cons.hmore = zinterpnc(tv,hmore,temps);
z0dinput.cons.zeff =  min(zeff,2.3) - min(zeff - 1,1.3) .* max(0,(z0dinput.cons.temps - t_xpoint) ./ (z0dinput.cons.temps(1) - t_xpoint));
zeff_dens =  (zeff - 1) .*  min(ip_xpoint  ./ z0dinput.cons.ip,max(z0dinput.cons.nbar) ./ z0dinput.cons.nbar) + 1;

% we have Zeff after early ramp up that depend on input power (impurities sources and seeding to protect divertor)
pin   = z0dinput.cons.ip + z0dinput.cons.picrh  + z0dinput.cons.plh + z0dinput.cons.pecrh;
fzeff = min(1, pin ./ max(pin)) .* (z0dinput.cons.temps > t_xpoint);
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
sepa_option.ton       = t_xpoint;
sepa_option.toff      = tdip;
z0dinput = z0separatrix(z0dinput,sepa_option,0);
z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;

% print separtrix parameters:
disp('Separatrix parameters (flat top phase):');
sepa_option

% wall information
rwall = lcfs_data.wall.r(:);
zwall = lcfs_data.wall.z(:);

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
a_mem    = sepa_option.a;
ra_mem   = sepa_option.ra;
Kref = z0dinput.geo.K(indkreduc(1) - 1);
dref = z0dinput.geo.d(indkreduc(1) - 1);

for k = indkreduc'
    sepa_option.a = min(z0dinput.geo.a(k),max(0.1 .* a_mem,fKreduc(k) .* a_mem));
    sepa_option.ra = max(min(rwall) + 1e-2,max(min(z0dinput.geo.R(indkreduc)),ra_mem + (sepa_option.a - a_mem)));
    sepa_option.zxup = max(1,fKreduc(k) .* zxup_mem);
    sepa_option.zxdo = max(1,fKreduc(k) .* zxdo_mem);
    sepa_option.rxup = rxup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.rxdo = rxdo_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.apup = apup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    sepa_option.amup = amup_mem .* z0dinput.cons.ip(k) ./ ip_flattop .* fKreduc(k);
    z0d1t = zerod_get1t(z0dinput,k);
    z0d1t.geo.K = min(z0d1t.geo.K,max(1,Kref .* fKreduc(k)));
    z0d1t.geo.d = min(z0d1t.geo.d,dref .* z0dinput.cons.ip(k) ./ ip_flattop);
    z0d1t.geo.a = sepa_option.a;
    z0d1t.geo.R = sepa_option.ra;
    rep = z0separatrix(z0d1t,sepa_option,0);
    z0dinput.geo.K(k) = rep.geo.K;
    z0dinput.geo.d(k) = rep.geo.d;
    z0dinput.geo.R(k) = rep.geo.R;
    z0dinput.geo.a(k) = rep.geo.a;
    z0dinput.exp0d.Rsepa(k,:) = rep.exp0d.Rsepa(1,:) - (max(rep.exp0d.Rsepa(1,:)) + min(rep.exp0d.Rsepa(1,:))) ./ 2 + z0d1t.geo.R;
    % kept xpoint position fixed
    z0dinput.exp0d.Zsepa(k,:) = rep.exp0d.Zsepa(1,:) - min(rep.exp0d.Zsepa(1,:)) - z0d1t.geo.z0 + z_xpoint_ref;
    % wall collision detection
    mask = zinout(rwall,zwall-z0dinput.geo.z0(k),z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:));
    if any(mask == 0)
	dz = 1e-3;
	nb = sum(mask == 0);
	nb_mem = nb;
	iter_loop = 1000;
	while (nb ~= 0) && (iter_loop > 0)
	    mask = zinout(rwall,zwall-z0dinput.geo.z0(k),z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:));
	    nb = sum(mask == 0);
	    if nb > nb_mem
		break;
	    end
	    nb_mem = nb;
	    iter_loop = iter_loop - 1;
	    z0dinput.geo.z0(k) = z0dinput.geo.z0(k) + dz;
	end
	z0dinput.geo.z0(k) = z0dinput.geo.z0(k) + 1e-2; % don't touch
	%disp(nb)
	%figure(21);clf
	%plot(rwall,zwall,z0dinput.exp0d.Rsepa(k,:),z0dinput.exp0d.Zsepa(k,:) + z0dinput.geo.z0(k));
	%drawnow
    end
    % time regularisation 
    z0dinput.geo.d(k) = max(0,min(z0dinput.geo.d(k),z0dinput.geo.d(k - 1)));
    z0dinput.geo.K(k) = max(1,min(z0dinput.geo.K(k),z0dinput.geo.K(k - 1)));
    z0dinput.geo.a(k) = min(z0dinput.geo.a(k),z0dinput.geo.a(k - 1));
    z0dinput.geo.R(k) = min(z0dinput.geo.R(k),z0dinput.geo.R(k - 1));
end
z0dinput.geo.b0 = rb0 ./ z0dinput.geo.R;


% limitation of Zeff
% evolution of Zeff following gaz density
ulh = 0.25;
fh  = zeros(size(temps));
fh(temps >= t_flattop) = 1;
negr     = (z0dinput.cons.ip ./ 1e6) ./ z0dinput.geo.a .^ 2 ./ pi .* 1e20;
K95 = z0dinput.geo.K .* (1 + 0.95 .^ 4 ) ./ 2 + 0.5 .* (1 - 0.95 .^4);
d95 = z0dinput.geo.d .* (0.95 .^ 2);
q95 = 5 .* z0dinput.geo.a .^ 2 .* z0dinput.geo.b0 ./ (ip./1e6) ./ z0dinput.geo.R .* (1 + K95 .^ 2 .*  ...
    (1 + 2 .* d95 .^ 2 - 1.2 .* d95 .^ 3) ./ 2) .* (1.17 - 0.65 .* z0dinput.geo.a ./ z0dinput.geo.R) ./  ...
    (1 - (z0dinput.geo.a./z0dinput.geo.R) .^ 2 ) .^ 2;
nbar_nat = min(negr,max(1e13, 1e20 .* (reference_simulation.rb0 ./ q95) .^ 0.6  .*  (ulh + (1 - ulh) .* fh) .* z0dinput.option.fnbar_nat));
nbar_star     = max(z0dinput.cons.nbar,nbar_nat);
zeff_max           =  min((zeff - 1) .* max(nbar_star)  ./ nbar_star + 1,max(z0dinput.option.zmax,z0dinput.option.zimp));
z0dinput.cons.zeff =  min(zeff_max,z0dinput.cons.zeff);
%figure(11);clf;plot(z0dinput.cons.temps,z0dinput.cons.zeff);drawnow
% print option data
disp('METIS parameters:')
z0dinput.option
tsnapshot = (1/4) .* t_flattop_plus  + (3/4) .* tend;

% maximum power used for each heating sources:
fprintf('max(P_LHCD)  = %g MW\n',max(z0dinput.cons.picrh) ./ 1e6);
fprintf('max(P_ICRH) = %g MW\n',max(z0dinput.cons.plh) ./ 1e6);

% set METIS main windows title
txt = 'Metis : Fast tokamak simulator';
txt = sprintf('%s (%s@%d)',txt,z0dinput.machine,z0dinput.shot);
setappdata(0,'METIS_INTERFACE_TITLE',txt);
if isappdata(0,'METIS_FILENAME');
    rmappdata(0,'METIS_FILENAME');
end

function control =sepa_moments(Rsepa,Zsepa)

% calcul des moments
% calcul de R0 et Z0
maskrmax  = (Rsepa == max(Rsepa,[],2));
% recalcul des parametres sur le vecteur final
rmin  = min(Rsepa,[],2);
rmax  = max(Rsepa,[],2);
a = 0.5 .* (rmax - rmin);
R0 = 0.5 .* (rmax + rmin);
control.R0  = R0;
control.a   = a;
zmin  = min(Zsepa,[],2);
zmax  = max(Zsepa,[],2);
control.K    = (zmax - zmin) ./ 2 ./ a;
rzmax = Rsepa(min(find(Zsepa == zmax)));
rzmin = Rsepa(min(find(Zsepa == zmin)));
z0    = Zsepa(min(find(Rsepa == rmax)));
% IMAS definition
control.z0_geo= (zmax + zmin) ./ 2;
control.z0    = z0;
control.Ku    = (zmax - z0)  ./ a;
control.Kl    = (z0 - zmin)  ./ a;
control.d     = abs(rzmax + rzmin -  2 .* R0) ./ 2 ./ a;
control.du    = abs(rzmax - R0) ./ a;
control.dl    = abs(rzmin - R0) ./ a;

% Xpoint detection
indh     = find(Zsepa == zmax,1);
indh     = cat(2,indh - 3, indh - 2,indh - 1,indh,indh + 1,indh + 2,indh + 3); 
indh     = mod(indh-1,size(Zsepa,2))+1;
rh       = Rsepa(indh);
zh       = Zsepa(indh);
ph       = polyfit(rh,zh,2);
eh       = sqrt(mean((zh - polyval(ph,rh)).^ 2)) ./ (max(rh) - min(rh));
indl     = find(Zsepa == zmin,1);
indl     = cat(2,indl - 3, indl - 2,indl - 1,indl,indl + 1,indl + 2,indl + 3);
indl     = mod(indl-1,size(Zsepa,2))+1;
rl       = Rsepa(indl);
zl       = Zsepa(indl);
pl       = polyfit(rl,zl,2);
el       = sqrt(mean((zl - polyval(pl,rl)).^ 2)) ./ (max(rl) - min(rl));
control.x_point = {};
if el > 2e-2
  indlz  = find(zl == min(zl),1);
  control.x_point{end + 1}.r = rl(indlz);
  control.x_point{end}.z = zl(indlz);
end
if eh > 2e-2
  indhz  = find(zh == max(zh),1);
  control.x_point{end + 1}.r = rh(indhz);
  control.x_point{end}.z = zh(indhz);
end

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
    
