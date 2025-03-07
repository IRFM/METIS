% default METIS paramaters for WEST
function option = map_option_west(option)

option.gaz = 2;
option.neasser = 0;
option.Recycling = 0.7;
option.natural = 0;
option.ane = 11;
option.fn0a = 1;
option.fn0a_div = 0.1;
%
option.scaling = 12;
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
option.hollow      = 1;
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
option.berror = 0; % no breakdown simulation by default
option.L_eddy = 4.76e-4;
option.R_eddy = 4.1e-5;
option.C_eddy = 1;
option.B_eddy = 0.1;
option.I_eddy = 0;
option.p_prefill = 0.7e-03;
option.VV_volume = 58;
%
option.zeff = 0;
option.faccu = 0;
option.heat_acc = 0;
option.fne_acc = 0;
option.zmax = 8;
option.zimp = 7;
option.rimp = 1;
option.density_model ='minconv';
%
option.frad = 1;
option.matthews = -1;
option.fte_edge = 1;
option.gaunt = 1;
option.noncoronal = 0;
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
option.fzmax_div = -1;
option.W_effect = 1;
option.cw_factor = 0;
option.cw_offset = 5e-5;
option.factor_scale = 1;
option.fpower = 0.6000;
option.fR_target = 2.7 / 2.94;
option.mach_corr = 1;
option.Sq = -1; 
option.Recycling_plate = 1-eps;
option.imp_div = 'auto';
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
option.lhmode = 0;  % used as ECCD system for breakdown
option.etalh  = 0.8;  % perpendicular
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
option.cronos_regul= 3;
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

