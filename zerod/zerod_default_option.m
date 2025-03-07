% compatibilite entre version concernant les parametres
function [option,cons] = zerod_default_option(option,cons)

% compatibilite entre version
option.transitoire = 1;
if ~isfield(option,'vane')
	option.vane = 1;
end
if ~isfield(option,'tauhemul')
	option.tauhemul = 0;
end
if ~isfield(option,'l2hscaling')
	option.l2hscaling = 0;
end
if ~isfield(option,'frad')
	option.frad = 1;
end
if ~isfield(option,'fpped')
	option.fpped = 1;
end
if ~isfield(option,'wlh')
	option.wlh = 0;
end
if ~isfield(option,'npar0')
	option.npar0 = 2;
end
if ~isfield(option,'tae')
	option.tae = 0;
end
if ~isfield(option,'freqlh')
	option.freqlh = 3.7;
end
if ~isfield(option,'runaway')
	option.runaway = 1;
end
if ~isfield(option,'limode')
	option.limode = 1;
end
if ~isfield(option,'xiioxie')
	option.xiioxie = 2;
end
if ~isfield(option,'kishape')
	option.kishape = 3;
end
if ~isfield(option,'qdds')
	option.qdds = 0;
end
if ~isfield(option,'kidds')
	option.kidds = 1;
end
if ~isfield(option,'nbmax')
	option.nbmax = 31;
end
if ~isfield(option,'plotconv')
	option.plotconv = 0;
end
if ~isfield(option,'modeboot')
	option.modeboot = 0;
end
if ~isfield(option,'zimp')
	option.zimp = fix(max(3,option.zmax - 1));
end
if ~isfield(option,'rimp')
	option.rimp = 0.1;
end
if ~isfield(option,'df')
	option.df = 0;
end
if isfield(option,'li')
 	cons.li = option.li .* ones(size(cons.temps));
end
if ~isfield(cons,'xece') && isfield(option,'xece')
 	cons.xece = option.xece .* ones(size(cons.temps));
end
if isfield(cons,'xece') && ~isfield(option,'xece')
 	option.xece = mean(cons.xece);
end
if ~isfield(cons,'flux')
	cons.flux = zeros(size(cons.temps));
end 
if ~isfield(option,'HH_li')
   option.HH_li  = 0;
end
if ~isfield(option,'morphing')
   option.morphing  = 9.5;
end
if ~isfield(option,'mode_expo_inte')
   option.mode_expo_inte  = 0;
end
if ~isfield(option,'configuration')
   switch option.modeh
   case 0
   	option.configuration  = 1;
   case 1
   	option.configuration  = 3;
   otherwise
   	option.configuration  = 4;
   end
end
if ~isfield(option,'evolution')
   option.evolution  = 0;
end
if ~isfield(option,'init_geo')
   option.init_geo  = 0;
end
if ~isfield(option,'laochange')
   option.laochange  = 1;
end
if ~isfield(option,'taurotmul')
   option.taurotmul  = 3;
end
if ~isfield(option,'breakdown')
   option.breakdown  = 0.03;
end
if ~isfield(option,'nbicdmul')
   option.nbicdmul = 1;
end
if ~isfield(option,'dilution')
   option.dilution = 1;
end
if ~isfield(option,'tau_limitation')
   option.tau_limitation = 'Off';
end
if ~isfield(option,'eccdmul')
   option.eccdmul  = 1;
end
if ~isfield(option,'Recycling')
   option.Recycling = 0;
end
if ~isfield(option,'natural')
   option.natural = 0;
end
if ~isfield(option,'cronos_regul')
   option.cronos_regul = 0;
end
if ~isfield(option,'impur_rot')
   option.impur_rot = 'imp';
end
if ~isfield(option,'sol_lscale')
   option.sol_lscale = 0;
end
if ~isfield(option,'dwdt_method')
   option.dwdt_method = 'explicit';
end
if ~isfield(option,'tswitch')
   option.tswitch = Inf;
end
if ~isfield(option,'eioniz')
   option.eioniz = 25;
end
if ~isfield(option,'sol_model')
   option.sol_model    = 'scaling';
end
if ~isfield(option,'fcond')
   option.fcond   = 1;
end
if ~isfield(option,'fmom')
   option.fmom    = 1;
end
if ~isfield(option,'lcx')
    option.lcx    = 7;
end
if ~isfield(option,'fnesol')
    option.fnesol    = 0;
end
%  if ~isfield(option,'target_recycle')
%      option.target_recycle    = 0;
%  end
if ~isfield(option,'dwow_elm')
    option.dwow_elm = 0;
end
if ~isfield(option,'tau_elm_factor')
    option.tau_elm_factor = 1;
end
if ~isfield(option,'bootmul')
    option.bootmul = 1;
end
if ~isfield(option,'ffit_ped')
    option.ffit_ped = 1;
end
if ~isfield(option,'fintrinsic')
	option.fintrinsic = 1;
end
if ~isfield(option,'upshiftmode')
	option.upshiftmode = 'linear';
end
if ~isfield(option,'fupshift')
	option.fupshift = 1;
end
if ~isfield(option,'ddsmode')
    option.ddsmode = 0;
end
if ~isfield(option,'w1')
    option.w1 = 0.5 ;
end
if ~isfield(option,'epsq')
    option.epsq = 1e-2;
end
if ~isfield(option,'nb_nbi')
    option.nb_nbi = 1;
end
if ~isfield(option,'einj2')
    option.einj2 = 1;
end
if ~isfield(option,'nea_factor')
       option.nea_factor = 1;
end
if ~isfield(option,'fnbar_nat')
  option.fnbar_nat = 1;
end
if ~isfield(option,'density_model')
  option.density_model = 'default';
end
if ~isfield(option,'npar_neg')
  option.npar_neg = 0;
end
if ~isfield(option,'ne_shape')
  option.ne_shape = 'Auto';
end
if ~isfield(option,'itb_sensitivity')
  option.itb_sensitivity = 1;
end
if ~isfield(option,'itb_slope_max')
  option.itb_slope_max = 2;
end
if ~isfield(option,'faccu')
  option.faccu = 0.5;
end
if ~isfield(option,'xieorkie')
  option.xieorkie = 0;
end
if ~isfield(option,'usepped_scl')
  option.usepped_scl = 0;
end
if ~isfield(option,'omega_shape')
  option.omega_shape = 0;
end
if ~isfield(option,'peeling')
  option.peeling = 0;
end
if ~isfield(option,'e_shielding')
  option.e_shielding    = 'Lin-Liu';
end
if ~isfield(option,'z_prad')
  option.z_prad    = 'zmax';
end
if ~isfield(option,'W_effect')
  option.W_effect  = 0;
end
if ~isfield(option,'sol_rad')
  option.sol_rad  = 'coupled';
end
if ~isfield(option,'lambda_scale')
  option.lambda_scale  = 0;
end
if ~isfield(option,'heat_acc')
  option.heat_acc  = 0;
end
if ~isfield(option,'heat_rot')
  option.heat_rot  = 0;
end
if ~isfield(option,'gradient_Wacc')
  option.gradient_Wacc  = 0;
end
if ~isfield(option,'fzmax_div')
  option.fzmax_div  = 0;
end

if ~isfield(option,'mul_blanket')
  option.mul_blanket = 1;
end

if ~isfield(option,'yield_model')
 	option.yield_model      = 'Javev';
end

if ~isfield(option,'berror')
	option.berror    = 0;
end
if ~isfield(option,'L_eddy')
	option.L_eddy    = 0;
end
if ~isfield(option,'R_eddy')
	option.R_eddy    = 0;
end
if ~isfield(option,'ki_expo')
	option.ki_expo    = 2;
end
if ~isfield(option,'p_prefill')
	option.p_prefill    = 2;
end
if ~isfield(option,'fstiff')
	option.fstiff    = option.fpped;
end
if ~isfield(option,'cw_factor')
	option.cw_factor    = 1;
end
if ~isfield(option,'cw_offset')
	option.cw_offset    = 0;
end
if ~isfield(option,'ftwleak')
	option.ftwleak    = -0.5;
end
if ~isfield(option,'plhthr')
	option.plhthr    = 'pel+pion';
end
if ~isfield(option,'ploss_exp')
	option.ploss_exp ='with_prad';
end
if ~isfield(option,'factor_scale')
	option.factor_scale = 1;
end
if ~isfield(option,'fpower')
	option.fpower = 1;
end
if ~isfield(option,'fR_target')
	option.fR_target = 1;
end
if ~isfield(option,'de')
	option.de = 0.5;
end
if ~isfield(option,'mach_corr')
	option.mach_corr = 0;
end
if ~isfield(option,'fne_acc')
	option.fne_acc = 0;
end
if ~isfield(option,'grad_ped')
	option.grad_ped = 0;
end
if ~isfield(option,'angle_ece2')
	option.angle_ece2    = 90;
end
%  % internal option
%  if ~isfield(option,'notimederivative')
%  	option.notimederivative = 0;
%  end
if ~isfield(option,'alpha_e')
	option.alpha_e  = 0.82;
end
if ~isfield(option,'hysteresis')
	option.hysteresis = 1;
end
if ~isfield(option,'gaunt')
	option.gaunt = 0;
end
if ~isfield(option,'fspot')
	option.fspot = 0.15;
end
if ~isfield(option,'noncoronal')
	option.noncoronal = 0;
end
if ~isfield(option,'cur_nbi_time')
	option.cur_nbi_time = 0;
end
if ~isfield(option,'icrh_width')
	option.icrh_width = 1;
end
if ~isfield(option,'fact_mino')
	option.fact_mino = 0;
end
if ~isfield(option,'alpha_channeling')
	option.alpha_channeling = 0;
end
if ~isfield(option,'ftaup')
	option.ftaup = 1;
end
if ~isfield(option,'fn0a')
	option.fn0a = 1;
end
if ~isfield(option,'fn0a_div')
	option.fn0a_div = option.fn0a;
end
if ~isfield(option,'neped_expo')
    option.neped_expo = - 0.7;
end
if ~isfield(option,'ne_free')
    option.ne_free = 3;
end
if ~isfield(option,'fpl2h_lim')
	option.fpl2h_lim = 1;
end
if ~isfield(option,'ton_modeh')
	option.ton_modeh = Inf;
end
if ~isfield(option,'toff_modeh')
	option.toff_modeh = Inf;
end
if ~isfield(option,'equi_ppar')
	option.equi_ppar = 0;
end
if ~isfield(option,'xiioxie_ped')
	option.xiioxie_ped = 0;
end
if ~isfield(option,'HH_gas_puff')
	option.HH_gas_puff = 0;
end
if ~isfield(option,'hmore_pped')
	option.hmore_pped = 0;
end
if ~isfield(option,'pl2h_mass_charge')
	option.pl2h_mass_charge= 0;
end
if ~isfield(option,'fte_edge')
	option.fte_edge = 1;
end
if ~isfield(option,'drs1')
        option.drs1 = 0;
end
if ~isfield(option,'dzs1')
        option.dzs1 = 0;
end    
if ~isfield(option,'drs2')
        option.drs2 = 0;
end
if ~isfield(option,'dzs2')
        option.dzs2 = 0;
end    
if ~isfield(option,'C_eddy')
        option.C_eddy = 1;
end    
if ~isfield(option,'B_eddy')
        option.B_eddy = 1;
end    
if ~isfield(option,'I_eddy')
        option.I_eddy = 1;
end    
if ~isfield(option,'PSI_eddy')
        option.PSI_eddy = 1;
end    
if ~isfield(option,'carbonblow')
        option.carbonblow = 0;
end    
if ~isfield(option,'temp_vac')
        option.temp_vac = 300;
end    
if ~isfield(option,'VV_volume')
        option.VV_volume = 0;
end    
if ~isfield(option,'initiation_only')
        option.initiation_only = 0;
end    
if ~isfield(option,'ode_pped')
        option.ode_pped = 0;
end    
if ~isfield(option,'coef_shape')
        option.coef_shape = 'bgb';
end    
if ~isfield(option,'s1crit')
        option.s1crit = 0;
end    
if ~isfield(option,'rotation_scl')
    if option.fintrinsic < 0
        option.rotation_scl = 'Barnes & Parra';
    else
        option.rotation_scl = 'Rice';
    end
end
if ~isfield(option,'eta_gas_puff')
        option.eta_gas_puff = 1;
end    
if ~isfield(option,'adiabatic')
        option.adiabatic = 0;
end    
if ~isfield(option,'icrh_model')
        option.icrh_model = 'PION_fit-Stix';
end    
if ~isfield(option,'ifast_icrh')
        option.ifast_icrh = 0;
end    
if ~isfield(option,'fabs_fw')
        option.fabs_fw = 0;
end    
if ~isfield(option,'fMC_loss')
        option.fMC_loss = 0;
end    
if ~isfield(option,'orbit_width')
        option.orbit_width = 0;
end    
if ~isfield(option,'mode_vtheta')
        option.mode_vtheta    = 'Neoclassical V_pol';
end    
if ~isfield(option,'rot_jr_loss')
        option.rot_jr_loss    = 'off';
end    
if ~isfield(option,'width_ecrh')
        option.width_ecrh    =  0;
end    
if ~isfield(option,'cw_ecrh')
        option.cw_ecrh    =  0;
end    
if ~isfield(option,'cw_icrh')
        option.cw_icrh    =  0;
end    
if ~isfield(option,'cw_lhcd')
        option.cw_lhcd    =  0;
end    
if ~isfield(option,'cw_nbi1')
        option.cw_nbi1    =  0;
end    
if ~isfield(option,'cw_nbi2')
        option.cw_nbi2    =  0;
end    
if ~isfield(option,'isotope_stiff')
        option.isotope_stiff = 0;
end
if ~isfield(option,'refined_ptot')
        option.refined_ptot = 0;
end
if ~isfield(option,'imp_div')
        option.imp_div = 'fixed';
end
if ~isfield(option,'nea_model')
        option.nea_model = 'Mahdavi';
end
if ~isfield(option,'fast_ion_sbp')
        option.fast_ion_sbp = 0;
end
if ~isfield(option,'shinethrough')
        option.shinethrough = 0;
end
if ~isfield(option,'disrup')
        option.disrup = 0;
end
if ~isfield(option,'orientation')
        option.orientation = 1;
end
if ~isfield(option,'acc_col')
        option.acc_col = 'off';
end
if ~isfield(option,'solid_rotation')
        option.solid_rotation = 1;
end
if ~isfield(option,'itb_density')
        option.itb_density = 1;
end
if ~isfield(option,'residence_time')
        option.residence_time = 0;
end
if ~isfield(option,'hollow')
        option.hollow = 0;
end
if ~isfield(option,'Sq')
        option.Sq = 0;
end
if ~isfield(option,'Recycling_target')
        option.Recycling_target = 0;
end
if ~isfield(option,'detach')
        option.detach = 0;
end
if ~isfield(option,'short_dt')
        option.short_dt = 'off';
end
if ~isfield(option,'exp_shape')
        option.exp_shape = 0;
end
if ~isfield(option,'MC_onoff')    
    option.MC_onoff  = 'on';
end
if ~isfield(option,'moments_mode')    
    option.moments_mode  = 1;
end
if ~isfield(option,'q0_dds_trig')    
    option.q0_dds_trig  = 0;
end
if ~isfield(option,'betap1crit')    
    option.betap1crit = 0;
end
if ~isfield(option,'HH_delta')    
    option.HH_delta = 0;
end
if ~isfield(option,'protect_sepa_z0')    
    option.protect_sepa_z0 = 'none';
end
if ~isfield(option,'fpolarized')    
    option.fpolarized = 0;
end
if ~isfield(option,'Kappa_xpoint')
    option.Kappa_xpoint = 0;
end
if ~isfield(option,'delta_xpoint')
    option.delta_xpoint = 0;
end
if ~isfield(option,'ane_factor')
    option.ane_factor = 1;
end
if ~isfield(option,'force_spitzer')
    option.force_spitzer = 0;
end
if ~isfield(option,'neutral_friction')
    option.neutral_friction = 0;
end
if ~isfield(option,'f_eta_turb')
    option.f_eta_turb = 0;
end
if ~isfield(option,'R_LFS_xpoint')
    option.R_LFS_xpoint = 0;
end
if ~isfield(option,'R_HFS_xpoint')
    option.R_HFS_xpoint = 0;
end
if ~isfield(option,'cx_ion')
    option.cx_ion = 0;
end
if ~isfield(option,'te_edge_fixed')
	option.te_edge_fixed = 0;
end
if ~isfield(option,'forced_H_NBI')
	option.forced_H_NBI = 0;
end
if ~isfield(option,'Sn_fraction')
    option.Sn_fraction = 0;
end
if ~isfield(option,'te_max')
    option.te_max = 1e5;
end
if ~isfield(option,'extended_qei')
    option.extended_qei = 'off';
end
if ~isfield(option,'cor_rel_spitzer')
    option.cor_rel_spitzer = 'off';
end
if ~isfield(option,'cor_rel_brem')
    option.cor_rel_brem = 'standart';
end
if ~isfield(option,'natural_nD_o_nH')
    option.natural_nD_o_nH    = 0.000115;
end
if ~isfield(option,'collapse')
    option.collapse    = 0;
end




