%Function to define input TGLF parameters which can be changed by the user
function param = declare_TGLF_op_parameters

% METIS implementation parameters
valeur.numproc  	= maxNumCompThreads;
type.numproc    	= 'real';
borne.numproc   	= [0,100];
defaut.numproc  	= maxNumCompThreads;
info.numproc    	= 'Number of processors to be used';
  
valeur.save_all_figs	= 0;
type.save_all_figs	= 'real';
borne.save_all_figs	= {0,1};
defaut.save_all_figs	= 0;
info.save_all_figs	= 'Save intermediate figures for internal and external loops. Used for debugging.';

valeur.save_intermediate_data	= 0;
type.save_intermediate_data	= 'real';
borne.save_intermediate_data	= {0,1};
defaut.save_intermediate_data	= 0;
info.save_intermediate_data	= 'Save intermediate data, including transport coefficients and kinetic profiles. Used for debugging.';
  
valeur.smooth_space	= 1;
type.smooth_space	= 'real';
borne.smooth_space	= {0,1};
defaut.smooth_space	= 1;
info.smooth_space	= 'if = 1, applied a filter on space direction to smooth transport coeffecient';

valeur.anti_aliasing		= 0.1;
type.anti_aliasing		= 'real';
borne.anti_aliasing		= [0,1];
defaut.anti_aliasing		= 0.1;
info.anti_aliasing		= 'anti aliasing factor';
%label.anti_aliasing
  
valeur.temperature_factor	= 1e-6;
type.temperature_factor		= 'real';
borne.temperature_factor	= [0,1];
defaut.temperature_factor	= 1e-6;
info.temperature_factor		= 'temperature factor for simulated annealing';
%label.temperature_factor
  
valeur.nbloop_max	= 501;
type.nbloop_max		= 'real';
borne.nbloop_max	= [0,1];
defaut.nbloop_max	= 501;
info.nbloop_max		= 'maximum number of loop';
%label.nbloop_max

valeur.prediction_domain	= 'Full';
type.prediction_domain		= 'string';
%borne.prediction_domain		= {'Full','Restricted','Clamped'};
borne.prediction_domain		= {'Full','Restricted'};
defaut.prediction_domain	= 'Full';
info.prediction_domain		= 'if = Full, use QLKZ to predict up to magnetic axis otherwise in restricted mode, predict Te, Ti  and Ne only up to x= 0.25 or sawtooth inversion radius if this one is mode external';

valeur.domain_rho		= 1;
type.domain_rho	    	        = 'real';
borne.domain_rho		= [0,1];
defaut.domain_rho		= 1;
info.domain_rho		        = 'Normalised radial coordinate for the boundary condition. Must be used with prediction domain=Full';

valeur.factor_neo_chii		= 1;
type.factor_neo_chii		= 'real';
borne.factor_neo_chii		= [0,100];
defaut.factor_neo_chii		= 1;
info.factor_neo_chii		= 'factor applied ot Chi_neo (computed with Hinton formula) that is added to turbulent coefficient Chi_ion';
%label.factor_neo_chii

valeur.neo_ion_in_D		= 0;
type.neo_ion_in_D		= 'real';
borne.neo_ion_in_D		= [0,1];
defaut.neo_ion_in_D		= 0;
info.neo_ion_in_D		= 'fraction of Chi_ion_neo limiting the minimal value of De';
%label.neo_ion_in_D

valeur.neo_ion_in_Chie		= 0;
type.neo_ion_in_Chie		= 'real';
borne.neo_ion_in_Chie		= [0,1];
defaut.neo_ion_in_Chie		= 0;
info.neo_ion_in_Chie		= 'fraction of Chi_ion_neo limiting the minimal value Chi_e';
%label.neo_ion_in_Chie
   
valeur.minChie	        = 0.05;
type.minChie		= 'real';
borne.minChie		= [0,1];
defaut.minChie		= 0.05;
info.minChie		= 'minimum value of Chi_e';

valeur.minChii	        = 0.05;
type.minChii		= 'real';
borne.minChii		= [0,1];
defaut.minChii		= 0.05;
info.minChii		= 'minimum value of Chi_i';

valeur.minDe	        = 0.05;
type.minDe		= 'real';
borne.minDe		= [0,1];
defaut.minDe		= 0.05;
info.minDe		= 'minimum value of D_e';

valeur.minVe	        = -100;
type.minVe		= 'real';
borne.minVe		= [-1000,0];
defaut.minVe		= -100;
info.minVe		= 'minimum value of V_e';

valeur.maxChie	        = 100;
type.maxChie		= 'real';
borne.maxChie		= [1,1000];
defaut.maxChie		= 100;
info.maxChie		= 'maximum value of Chi_e';

valeur.maxChii	        = 100;
type.maxChii		= 'real';
borne.maxChii		= [1,1000];
defaut.maxChii		= 100;
info.maxChii		= 'maximum value of Chi_i';

valeur.maxDe	        = 100;
type.maxDe		= 'real';
borne.maxDe		= [1,1000];
defaut.maxDe		= 100;
info.maxDe		= 'maximum value of D_e';

valeur.maxVe	        = 100;
type.maxVe		= 'real';
borne.maxVe		= [1,1000];
defaut.maxVe		= 100;
info.maxVe		= 'maximum value of V_e';


valeur.coef_min		= 1e-3;
type.coef_min		= 'real';
borne.coef_min		= [0,10];
defaut.coef_min		= 1e-3;
info.coef_min		= 'minimum value of diffusion coefficient what ever is other tunning for convergence of the solver';
%label.coef_min

valeur.calc_heat_transport	= 1;
type.calc_heat_transport	= 'real';
borne.calc_heat_transport	= {0,1};
defaut.calc_heat_transport	= 1;
info.calc_heat_transport	= 'Calculate q_e and q_i';
%label.calc_heat_transport

valeur.calc_part_transport	= 0;
type.calc_part_transport	= 'real';
borne.calc_part_transport	= {0,1};
defaut.calc_part_transport	= 0;
info.calc_part_transport	= 'Calculate Gamma';
%label.calc_part_transport

valeur.verbosity	= 0;
type.verbosity		= 'real';
borne.verbosity		= {0,1,2,3,4,5,10};
defaut.verbosity	= 0;
info.verbosity		= 'Level of verbosity (debugging parameter)';
%label.verbosity

valeur.channel_tag	= 'TGLF';
type.channel_tag	= 'string';
borne.channel_tag	= {'TGLF'};
defaut.channel_tag	= 'TGLF';
info.channel_tag	= 'Tag in plots for this model';
%label.channel_tag

%valeur.SWITCH_MODEL = 1;
%type.SWITCH_MODEL   = 'real';
%borne.SWITCH_MODEL  = {0,1}; % {} for the list , [] for the range 
%defaut.SWITCH_MODEL = 1;
%info.SWITCH_MODEL   = 'How to run TGLF:\n 0: with 3 species, ions, electrons, light impurity,\n 1 all METIS species (6)';

valeur.SWITCH_tglf = 'fast';
type.SWITCH_tglf  = 'string';
borne.SWITCH_tglf  = {'fast', 'long'}; % {} for the list , [] for the range 
defaut.SWITCH_tglf = 'fast';
info.SWITCH_tglf   = 'How to run TGLF:\n fast: with just one time; long: with four time slices';

valeur.SWITCH_cushion_diff_coef = 'yes';
type.SWITCH_cushion_diff_coef  = 'string';
borne.SWITCH_cushion_diff_coef = {'yes', 'no'}; % {} for the list , [] for the range 
defaut.SWITCH_cushion_diff_coef = 'yes';
info.SWITCH_cushion_diff_coef   = 'yes: The diffusion coefficients are modified using (1 - f) .* chie_mem + f .* chie_loop before calculating the temperature. For testing purposes.';

valeur.FILTER = 2;
type.FILTER   = 'real';
borne.FILTER  = {0,1,2,3};
defaut.FILTER = 2;
info.FILTER  = 'Sets threshold for frequency/drift frequnecy to filter out non-driftwave instabilities';

%valeur.NS = 6;
%type.NS   = 'real';
%borne.NS  = [1,7];
%defaut.NS = 4;
%info.NS   = 'Number of species including both electrons and ions.';

valeur.USE_TRANSPORT_MODEL = 'true';
type.USE_TRANSPORT_MODEL   = 'string';
borne.USE_TRANSPORT_MODEL  = {'true','false'}; 
defaut.USE_TRANSPORT_MODEL = 'true';
info.USE_TRANSPORT_MODEL   = 'Transport Model';

valeur.GEOMETRY_FLAG = 1;
type.GEOMETRY_FLAG   = 'real';
borne.GEOMETRY_FLAG  = {0,1,2,3}; % {} for the list , [] for the range 
defaut.GEOMETRY_FLAG = 1;
info.GEOMETRY_FLAG   = 'Geometry type. 0: s-alpha, 1: Miller, 2: Fourier, 3: Elite';

valeur.USE_BPER = 'false';
type.USE_BPER   = 'string';
borne.USE_BPER  = {'true','false'}; 
defaut.USE_BPER = 'false';
info.USE_BPER   = 'Include transverse magnetic fluctuations';

valeur.USE_BPAR = 'false';
type.USE_BPAR   = 'string';
borne.USE_BPAR  = {'true','false'}; 
defaut.USE_BPAR = 'false';
info.USE_BPAR   = 'Include compressional magnetic fluctuations';

valeur.USE_BISECTION = 'true';
type.USE_BISECTION  = 'string';
borne.USE_BISECTION  = {'true','false'};
defaut.USE_BISECTION = 'true';
info.USE_BISECTION   = 'Use bisection search method to find width that maximizes the growth rate.';

valeur.USE_MHD_RULE  = 'true';
type.USE_MHD_RULE   = 'string';
borne.USE_MHD_RULE  = {'true','false'};
defaut.USE_MHD_RULE= 'true';
info.USE_MHD_RULE   = 'Ignore pressure gradient contribution to curvature drift. Recommeded:  Set false for high beta.';

valeur.USE_INBOARD_DETRAPPED = 'false';
type.USE_INBOARD_DETRAPPED   = 'string';
borne.USE_INBOARD_DETRAPPED  = {'true','false'};
defaut.USE_INBOARD_DETRAPPED = 'false';
info.USE_INBOARD_DETRAPPED   = 'Set trapped fraction to zero if eigenmode is inward balloonong.';

valeur.SAT_RULE = 2;
type.SAT_RULE   = 'real';
borne.SAT_RULE  = {0,1,2,3};
defaut.SAT_RULE = 0;
info.SAT_RULE   = '0: finds zonal flow shear at each ky, 1: dominant saturation mechanism (ZF mixing rate or drift-wave growth rate) and includes ky coupling, 2: builds on 1 with refined geometric effects, improved TEM physics and simplified spectral shift, 3: builds on 2, captures anti-gyrobohm scaling of fluxes, treats saturation of ITG and TEM differently, has quasi linear model approximations to align quasi linear weights with NLGK simulations.';

valeur.KYGRID_MODEL = 4;
type.KYGRID_MODEL   = 'real';
borne.KYGRID_MODEL  = {0,1,4};
defaut.KYGRID_MODEL = 4;
info.KYGRID_MODEL   = '0: user defined ky grid with nky modes up to maimum ky, equal spaced with kymin = ky/nky, 1: standard ky spectrum often used for sat0 and sat1, where kymin=0.1/rho_ion, 4: additional low-ky modes preferred for sat2, where kymin = 0.05*grad_r0/rho_ion ';

valeur.XNU_MODEL  = 3;
type.XNU_MODEL   = 'real';
borne.XNU_MODEL  = {2,3};
defaut.XNU_MODEL = 3;
info.XNU_MODEL   = ' Collision model, 2: defaut preset for sat0 and sat1, 3: defaut preset for sat2.';

valeur.VPAR_MODEL  = 0;
type.VPAR_MODEL   = 'real';
borne.VPAR_MODEL  = {0};
defaut.VPAR_MODEL = 0;
info.VPAR_MODEL   = '0: low-Mach number limit.';

valeur.SIGN_BT  = -1;
type.SIGN_BT   = 'real';
borne.SIGN_BT  = {-1,1};
defaut.SIGN_BT = -1;
info.SIGN_BT   = 'Sign of B_t w.r.t counter clockwise toroidal direction from top.';

valeur.SIGN_IT  = -1;
type.SIGN_IT    = 'real';
borne.SIGN_IT   = {-1,1};
defaut.SIGN_IT  = -1;
info.SIGN_IT    = 'Sign of I_t w.r.t counter clockwise toroidal direction from top.';

valeur.KY = 0.3;
type.KY   = 'real';
borne.KY  = [0,1];
defaut.KY = 0.3;
info.KY   = 'k_theta rho_s,unit for single mode call to TGLF or user-defined ky grid.';

valeur.NEW_EIKONAL  = 'true';
type.NEW_EIKONAL   = 'string';
borne.NEW_EIKONAL  = {'true','false'};
defaut.NEW_EIKONAL= 'true';
info.NEW_EIKONAL   = 'Compute the eikonal';

valeur.VEXB_SHEAR   = 0.0;
type.VEXB_SHEAR    = 'real';
borne.VEXB_SHEAR   = [0,1];
defaut.VEXB_SHEAR = 0.0;
info.VEXB_SHEAR    = 'Normalised toroidal ExB velocity Doppler shift gradient common to all species. Not coupled to METIS data yet.';

valeur.DEBYE = 0.0;
type.DEBYE   = 'real';
borne.DEBYE  = [0,1];
defaut.DEBYE = 0.0;
info.DEBYE   = 'Debye Length/Gyroradius';

valeur.IFLUX = 'true';
type.IFLUX    = 'string';
borne.IFLUX   = ['true,false'];
defaut.IFLUX  = true;
info.IFLUX    = 'Compute quasi linear weights and mode amplitudes.';

valeur.IBRANCH = -1;
type.IBRANCH    = 'real';
borne.IBRANCH   = {-1,0};
defaut.IBRANCH  = -1;
info.IBRANCH    = '0: Find two most unstable modes one for each sign of frequency, electron drift direction (1), ion drift direction (2), -1:sort the unstable modes by growthrate in rank order.';

valeur.NMODES = 2;
type.NMODES   = 'real';
borne.NMODES  = {2};
defaut.NMODES = 2;
info.NMODES   = 'Number of modes to store for tglf_ibranch_in = -1';

valeur.NBASIS_MAX = 6;
type.NBASIS_MAX   = 'real';
borne.NBASIS_MAX  = {4,5,6};
defaut.NBASIS_MAX = 6;
info.NBASIS_MAX   = 'Maximum number of parallel basis functions';

valeur.NBASIS_MIN = 6;
type.NBASIS_MIN   = 'real';
borne.NBASIS_MIN  = {2,4,6};
defaut.NBASIS_MIN = 6;
info.NBASIS_MIN   = 'Manimum number of parallel basis functions';

valeur.NXGRID = 16;
type.NXGRID   = 'real';
borne.NXGRID  = {16};
defaut.NXGRID = 16;
info.NXGRID   = 'Number of nodes in Gauss-Hermite quadrature';

valeur.NKY = 19;
type.NKY   = 'real';
borne.NKY  = {12,13,14,15,16,17,18,19,20};
defaut.NKY = 19;
info.NKY   = 'Number of poloidal modes in the high k-spectrum of TGLF, i.e. number of modes with logarithmic spacing to cover roughly 1<ky<24';

valeur.ADIABATIC_ELEC = 'false';
type.ADIABATIC_ELEC   = 'string';
borne.ADIABATIC_ELEC  = {'true','false'};
defaut.ADIABATIC_ELEC = 'false';
info.ADIABATIC_ELEC   = 'use adiabatic electrons. They are assumed to be perfectly Maxwellian and remain so  '; % Look at it later


valeur.ALPHA_P = 1.0;
type.ALPHA_P   = 'real';
borne.ALPHA_P  = {1.0};
defaut.ALPHA_P = 1.0;
info.ALPHA_P   = 'Multiplies parallel velocity shear for all species.';


valeur.ALPHA_MACH = 0.0;
type.ALPHA_MACH   = 'real';
borne.ALPHA_MACH  = {0.0};
defaut.ALPHA_MACH = 0.0;
info.ALPHA_MACH   = 'Multiplies parallel velocity for all species.';

valeur.ALPHA_E = 1.0;
type.ALPHA_E  = 'real';
borne.ALPHA_E  = {1.0};
defaut.ALPHA_E = 1.0;
info.ALPHA_E   = 'Multiplies ExB velocity shear for spectral shift model.';

valeur.ALPHA_QUENCH = 0.0;
type.ALPHA_QUENCH  = 'real';
borne.ALPHA_QUENCH  = {0.0, 1.0};
defaut.ALPHA_QUENCH = 0.0;
info.ALPHA_QUENCH   = '0: use new spectral shift model, 1: use quench rule.';

valeur.XNU_FACTOR = 1.0;
type.XNU_FACTOR   = 'real';
borne.XNU_FACTOR  = [0,2];
defaut.XNU_FACTOR =1.0;
info.XNU_FACTOR   = 'Multiplies the trapped/passing boundary electron-ion collision terms';

valeur.DEBYE_FACTOR = 1.0;
type.DEBYE_FACTOR   = 'real';
borne.DEBYE_FACTOR  = [0,2];
defaut.DEBYE_FACTOR = 1.0;
info.DEBYE_FACTOR  = 'Multiplies Debye Length/Gyroradius';

valeur.ETG_FACTOR = 1.25;
type.ETG_FACTOR   = 'real';
borne.ETG_FACTOR  = [0,2];
defaut.ETG_FACTOR = 1.25;
info.ETG_FACTOR  = 'Exponent for ETG saturation rule.';

valeur.WRITE_WAVEFUNCTION_FLAG = 0;
type.WRITE_WAVEFUNCTION_FLAG   = 'real';
borne.WRITE_WAVEFUNCTION_FLAG  = {0,1};
defaut.WRITE_WAVEFUNCTION_FLAG = 0;
info.WRITE_WAVEFUNCTION_FLAG  = 'Self-explanatory';

valeur.UNITS = 'GYRO';
type.UNITS   = 'string';
borne.UNITS  = {'GYRO','CGYRO'};
defaut.UNITS = 'GYRO';
info.UNITS  = 'GYRO: sat0, CGYRO: sat2';

interface.ts = '';    % name of the interfacing function with TS data 
interface.jet = '';   % name of the interfacing function with JET data

param.valeur     = valeur;
param.type       = type;
param.borne      = borne;
param.defaut     = defaut;
param.info       = info;
param.interface  = interface;
param.description = 'GUI for TGLF';   % description (une ligne) de la fonction
param.help     = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
param.gui      ='';                             % nom de l'interface graphique specifique si elle existe
param.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
