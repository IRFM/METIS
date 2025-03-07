function z0dinput = reactormetissimulation(reactor_option)

% parameters declaration
if nargin < 1
   reactor_option = [];
end
if (nargin <= 1) && ~isstruct(reactor_option)


    valeur.sepa_type = 'SN';       
    type.sepa_type   = 'string';
    borne.sepa_type  = {'SN','DN'};  
    defaut.sepa_type = 'SN';
    info.sepa_type   = 'type of LCFS: SN = single null & DN = double null';

    valeur.gas    = 3;   % gas species as in METIS (1=H, 2=D, 3=DT & 4=He)      
    type.gas      = 'integer';
    borne.gas     = {1,2,3,4};  
    defaut.gas    = 3;
    info.gas      = 'main gas species as in METIS: 1 -> H, 2 -> D, 3-> D/T & 4 -> He';

    valeur.R     = 6.2;    
    type.R       = 'float';
    borne.R      = [1,15];  
    defaut.R     = 6.2;
    info.R       = 'major radius (m)';

    valeur.a     = 2;    
    type.a       = 'float';
    borne.a      = [.1,4];  
    defaut.a     = 2;
    info.a       = 'minor radius (m)';

    valeur.K     = 1.844;    
    type.K       = 'float';
    borne.K      = [1,3];  
    defaut.K     = 1.844;
    info.K       = 'LCFS elongation';

    valeur.d     = 0.52;    
    type.d       = 'float';
    borne.d      = [-0.5,1];  
    defaut.d     = 0.52;
    info.d       = 'LCFS averaged up and down triangularity';

    valeur.ip     = 15;   % plasma current (MA)  
    type.ip       = 'float';
    borne.ip      = [0.5,100];  
    defaut.ip     = 15;
    info.ip      = 'plasma current (MA)';

    valeur.b0     = 5.3;  
    type.b0       = 'float';
    borne.b0      = [0.5,15];  
    defaut.b0     = 5.3;
    info.b0       = 'vacuum magnetic field @ R (T)';
    
    valeur.delta_int     = 1.23;    
    type.delta_int       = 'float';
    borne.delta_int      = [0.7,3];  
    defaut.delta_int     = 1.23;
    info.delta_int       = 'minimum distance between TF conductor and plasam (m)';
    
    valeur.f_Greenwald     = 0.85;   
    type.f_Greenwald       = 'float';
    borne.f_Greenwald      = [0.3,1.5];  
    defaut.f_Greenwald     = 0.85;
    info.f_Greenwaldd      = 'electron density Greenwald fraction';
    
    valeur.edge_density_factor     = 1;   
    type.edge_density_factor       = 'float';
    borne.edge_density_factor      = [0.1,10];  
    defaut.edge_density_factor     = 1;
    info.edge_density_factor       = 'factor applied to edge density scaling law';
    
    valeur.zeff     = 1.3;   
    type.zeff       = 'float';
    borne.zeff      = [1.1,7];  
    defaut.zeff     = 1.3;
    info.zeff       = 'flattop line averaged Zeff without He ashes contribution';
    
    valeur.H_H     = 1;   
    type.H_H       = 'float';
    borne.H_H      = [0.5,2];  
    defaut.H_H     = 1;
    info.H_H       = 'enhancement for the selected scaling law during H-mode phase';
    
    valeur.tau_He_o_tau_E_core     = 5;   
    type.tau_He_o_tau_E_core       = 'float';
    borne.tau_He_o_tau_E_core      = [1,10];  
    defaut.tau_He_o_tau_E_core     = 5;
    info.tau_He_o_tau_E_core       = 'ratio between core confinement time of He over energy confinement time';
    
    valeur.rw     = 0.7;   
    type.rw       = 'float';
    borne.rw      = [0,1];  
    defaut.rw     = 0.7;
    info.rw       = 'cyclotron radiation reflection coefficient';
    
    valeur.c_W     = 1e-5;   
    type.c_W       = 'float';
    borne.c_W      = [1e-10,1e-3];  
    defaut.c_W     = 1e-5;
    info.c_W       = 'tungsten concentration in core plasma';
        
    valeur.P_NBI     = 53;   
    type.P_NBI       = 'float';
    borne.P_NBI      = [5 500];  
    defaut.P_NBI     = 53;
    info.P_NBI       = 'maximum power for NBI during flattop (MW)';
    
    valeur.E_NBI     = 1;   
    type.E_NBI       = 'float';
    borne.E_NBI      = [0.2 2];  
    defaut.E_NBI     = 1;
    info.E_NBI       = 'neutral beam energy injection (MV)';
    
    valeur.Recycling     = 0.97;   
    type.Recycling       = 'float';
    borne.Recycling      = [0,1-1e-4];  
    defaut.Recycling     = 0.97;
    info.Recycling       = 'Recycling coefficient';
    
    valeur.available_flux     = 180;   
    type.available_flux       = 'float';
    borne.available_flux      = [3,3000];  
    defaut.available_flux     = 180;
    info.available_flux       = 'Poloidal available flux (Wb)';
       
    valeur.duration     = 400;   
    type.duration       = 'float';
    borne.duration      = [30,100000];  
    defaut.duration     = 400;
    info.duration       = 'shot duration including ramp-up and flat-top but without ramp-down (s)';
    
    valeur.CEjima     = 0;   
    type.CEjima       = 'float';
    borne.CEjima      = [0,1];  
    defaut.CEjima     = 0;
    info.CEjima       = 'Resistive flux comsumption factor during ramp-up;\nif = 0 use SYCOMORE setting (reserved to expert physicist)';
    
    valeur.rampup_dipdt_factor     = 1;   
    type.rampup_dipdt_factor       = 'float';
    borne.rampup_dipdt_factor      = [0.1,10];  
    defaut.rampup_dipdt_factor     = 1;
    info.rampup_dipdt_factor       = 'Factor applied to rampup rate current decrease (reserved to expert physicist)';
            
    valeur.rampdown_dipdt_factor     = 1;   
    type.rampdown_dipdt_factor       = 'float';
    borne.rampdown_dipdt_factor      = [0.1,10];  
    defaut.rampdown_dipdt_factor     = 1;
    info.rampdown_dipdt_factor       = 'Factor applied to rampdown rate current decrease (reserved to expert physicist)';
           
    valeur.device = 'TEST';       
    type.device   = 'string';
    borne.device  =  '';  
    defaut.device = 'TEST';
    info.device   = 'name of the machine';

    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

    z0dinput.valeur     = valeur;
    z0dinput.type       = type;
    z0dinput.borne      = borne;
    z0dinput.defaut     = defaut;
    z0dinput.info       = info;
    z0dinput.interface  = interface;

    z0dinput.description = 'METIS reactor scenario generator';   % description (une ligne) de la fonction

    z0dinput.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    z0dinput.gui      ='';                             % nom de l'interface graphique specifique si elle existe
    z0dinput.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

    return

end
% gestion of input number and contents
% geometry
delta_int = reactor_option.delta_int;
R         = reactor_option.R;
a         = reactor_option.a;
d         = reactor_option.d;
K         = reactor_option.K;
b0        = reactor_option.b0;
delta_int = reactor_option.delta_int;
% mapping the data
sycomore_data.ip = reactor_option.ip .* 1e6;                % flat top plasma current (A).
sycomore_data.q95 = 5 .* a .^ 2 .* b0 ./ (sycomore_data.ip ./ 1e6) ./ R .* (1 + K .^ 2 .*  ...
                   (1 + 2 .* d .^ 2 - 1.2 .* d .^ 3) ./ 2) .* (1.17 - 0.65 .* a ./ R) ./  ...
                   (1 - (a./R) .^ 2 ) .^ 2;               % flat top safety factor
sycomore_data.rb0 = b0 .* R;         % magnetic rigidity (flat top)
sycomore_data.available_flux = reactor_option.available_flux;     % Poloidal available flux (Wb)
sycomore_data.device = sprintf('%s-R%d-a%d-RBt%d-Ip%d',reactor_option.device,ceil(R*100),ceil(a*100),ceil(sycomore_data.rb0),ceil(sycomore_data.ip/1e6)); % device name (used to generate file name and comments)
sycomore_data.scaling  = 12;             % code for energy plasma content scaling law (same as METIS)
sycomore_data.H_H  = reactor_option.H_H;                 % enhancement factor for energy content on flat top
sycomore_data.f_Greenwald = reactor_option.f_Greenwald;       % Greenwald density fraction on flat top
sycomore_data.shot = ceil(rand*1e4);                 % shot number, will be used to write data in UAL
sycomore_data.f_ni = 0;                 % faction of non inductive current, if f_ni > 0, auxiliary power will be adjusted to have at least this value. 
sycomore_data.rimp = 0.1;               % ratio  between berillium and argon in core plasma
sycomore_data.Recycling = reactor_option.Recycling;         % recycling  @ divertor 
sycomore_data.zeff      = reactor_option.zeff ;          % line averaged Zeff without He ashes contribution
sycomore_data.tau_He_o_tau_E_core = - reactor_option.tau_He_o_tau_E_core;  % ratio between core confinement time of He over energy confinement time.
sycomore_data.rw = reactor_option.rw;                 % cyclotron radiation reflection coefficient
sycomore_data.PNBI = reactor_option.P_NBI .* 1e6;        % NBI power used in Sycomore simulation (W, for ITER like default  case contains also ICRH power, if PNBI <  0, scale on ITER power).
sycomore_data.fR_target = 1;            % position of outer target in unit of R0.
sycomore_data.nbi_e_inj = reactor_option.E_NBI .* 1e6;          % neutral beam energy injection.
sycomore_data.signe = 1;                % sign of Bt . J_phi: must be provided by the machine description
sycomore_data.flux_expansion = 3;       % divertor flux expansion 
sycomore_data.target_orientation = 90;  % poloidal orientation of outer divertor target  (degrees)
sycomore_data.S_factor = 3e-3;          % divertor flux spreading factor (m) 
sycomore_data.cW               = reactor_option.c_W;  % tungsten concentration in core plasma   , must be provide by Sycomore 
% There no enhancement of impurities concentration id divertor. have been checked  in soldiv module of Sycomore.
sycomore_data.fzmax_div        = 3;     % Argon concentration enhancement in divertor (%)  , must be provide by Sycomore
sycomore_data.f_DSOL           = 1;   % multiplicator applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling).
sycomore_data.f_ne_LCFS        = reactor_option.edge_density_factor;   % factor applied to edge scaling law for density:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar'electron density at LCFS (m^-3)  , must be provide by Sycomore
sycomore_data.couplage         = 0.15;   % coupling coefficient given the ratio between poloidal flux consumption in the plasam and in the central selenoid (CS) during breakdown to take into account dissipated flux in passive structure.
sycomore_data.duration         = reactor_option.duration;  % estimated duration of the shot in s
% expert keys
sycomore_data.CEjima                 = reactor_option.CEjima;
sycomore_data.rampup_dipdt_factor    = reactor_option.rampup_dipdt_factor;
sycomore_data.rampdown_dipdt_factor  = reactor_option.rampdown_dipdt_factor;
% use constant P/R for allowed shine througth (and first orbit losses)
% take < 5% of 33 MW for ITER @ 6.2m -> 2e5 W/m
%  2%
sycomore_data.shine_through_limit =  2e6;  % limit of power lost in shine through (W)

%  optionnal fields 
sycomore_data.tokamak  = '';           % new UAL tokamak database name
sycomore_data.user     = '';           % new UAL user database name
sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])
%
sycomore_data.path4metis_files = '';   % path for files save as a postprocessing of METIS (METIS data, figures, ...); if left empty, then no wrinting       
 
%  structure details for sepa_option (with data example for ITER):
switch reactor_option.sepa_type
case 'SN'
    Kref = (1.687 + 2.001) / 2;
    dref = (0.466 + 0.568) / 2;
    sepa_option.rxup      = 0.466 .* d ./ dref;     % upper triangularity (minor radius unit)
    sepa_option.zxup      = 1.687 .* K / Kref;    % upper altitude X point (minor radius unit)
    sepa_option.apup      = 0;       % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 0;       % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = R;       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = a;         % minor radius (m) [2]
    sepa_option.rxdo      = 0.568 .* d ./ dref;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = 2.001 .* K / Kref;       % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = 14.2;      % magnetic field at R0
    sepa_option.delta     = delta_int;      % magnetic field at R0
otherwise
    sepa_option.rxup      = d;     % upper triangularity (minor radius unit)
    sepa_option.zxup      = K;    % upper altitude X point (minor radius unit)
    sepa_option.apup      = 22.46;       % upper separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amup      = 67.92;       % upper separatrix angle (-R,X) (HFS, degrees)
    sepa_option.ra        = R;       % major radius R0 (m) [6.2]
    sepa_option.za        = 0;       % altitude of the magnetic axis (m) [0.9]
    sepa_option.a         = a;         % minor radius (m) [2]
    sepa_option.rxdo      = d;     % lower triangularity (minor radius unit)
    sepa_option.zxdo      = K;       % lower altitude X point (minor radius unit)
    sepa_option.apdo      = 22.46;   % lower separatrix angle (R,X)  (LFS, degrees)
    sepa_option.amdo      = 67.92;   % lower separatrix angle (-R,X)  (HFS, degrees)
    sepa_option.b0        = 14.2;      % magnetic field at R0
    sepa_option.delta     = delta_int;      % magnetic field at R0
end
% metis call
parameters_filename ='';
computation_mode = 2;
%  optionnal fields 
sycomore_data.tokamak  = '';           % new UAL tokamak database name
sycomore_data.user     = '';           % new UAL user database name
sycomore_data.dataversion = '';        % selected a differente version of data (not the last one)
sycomore_data.occurrence = '';         % input cpo scenario occurrence in Kepler (default = [])

z0dinput = zerod_init(-2);
option = z0dinput.option;
option.neasser = 1;
option.natural = 1;
option.ane = 11;
option.fn0a = 1;
option.fn0a_div = 0.1;
option.dilution = 1;
option.tau_limitation = 'On';
option.l2hscaling = 3;
option.pl2h_mass_charge = 1;
option.modeh = 1;
option.hysteresis = 0;
option.configuration = 2;
option.l2hslope =  0.5;
option.usepped_scl = 4;
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
option.sitb = 0;
option.itb_sensitivity = 1;
option.itb_slope_max = 2;
option.smhd = 100;
option.tmhd = 0;
%
option.runaway = 0;
option.modeboot = 2;
%
option.zeff = 0;
option.faccu = 0;
option.heat_acc = 0;
option.fne_acc = 0;
option.zmax = 4;
option.zimp = 4;
option.rimp = 0.1;
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
option.sol_model  = 'scaling';
option.lcx = 1;
option.fcond = -1;
option.fmom = 0;
option.lambda_scale = 3;
option.sol_rad = 'decoupled';
option.W_effect = 1;
option.factor_scale = 1;
option.fpower = 0.6000;
option.mach_corr = 1;
option.cw_factor = 0;
% reactivate fast ions ionisation of neutral beam
option.fast_ion_sbp   = 1;


z0dinput = sycomore2metis(sepa_option,option,sycomore_data,computation_mode);

z0dinput.option.gaz = reactor_option.gas;
switch reactor_option.gas
case 3
  % nothing
otherwise
 z0dinput.cons.iso(:) = 0; 
end

% add possible missing field on profile description dependning on function
% creating the simulation
if ~isfield(z0dinput,'profinfo')
    z0dinput.profinfo      = z0dprofinfo;
end

% set METIS main windows title
txt = 'Metis : Fast tokamak simulator';
txt = sprintf('%s (%s@%d)',txt,z0dinput.machine,z0dinput.shot);
setappdata(0,'METIS_INTERFACE_TITLE',txt);
if isappdata(0,'METIS_FILENAME')
    rmappdata(0,'METIS_FILENAME');
end