function zs = zerod_param(void)

    valeur.gaz    = 3;
    type.gaz      = 'integer';                    
    borne.gaz     = {1,2,3,4,5,11};  
    defaut.gaz    = 3;
    info.gaz      = '1 -> H, 2 -> D , 3 -> D-T, 4 -> He, 5 -> D-He3 & 11 -> p-B11';
    label.gaz     = 'gas';
    section.gaz   = 'Composition';
 
    valeur.frhe0  = 0;
    type.frhe0    = 'real';
    borne.frhe0   = [0,0.5];
    defaut.frhe0  = 0;
    info.frhe0    = 'ratio of residual helium, other than alpha from fusion DT';
    section.frhe0  = 'Composition';
    mode.frhe0    = 'advanced';
    
    valeur.tauhemul   = 0;
    type.tauhemul     = 'real';
    borne.tauhemul    = [-20,20];
    defaut.tauhemul   = 0;
    info.tauhemul  = 'helium confinement time: tau_He_star = tauhemul * taue; if = 0 , use the scaling law;\nif < 0 take into account recycling of neutral in the divertor :\n tau_He_star = tauhemul * taue + Recycling / (1 - Recycling) * tau_p';
    section.tauhemul  = 'Composition';

    valeur.neasser  = 0;
    type.neasser    = 'integer';                    
    borne.neasser   = {0,1,2};  
    defaut.neasser  = 0;                
    info.neasser    = 'ODE for density evolution: if = 0: zs.nbar = cons.nbar (no ODE solved);\nif = 1: density controlled by gas puff with reference cons.nbar using electron density confinment times;\n if = 2, as 1 and density limited to prevent disruption';
    section.neasser = 'Density';
    mode.neasser     = 'advanced';

    valeur.Recycling    = 0;
    type.Recycling      = 'real';
    borne.Recycling     = [0,0.999];
    defaut.Recycling    = 0;
    info.Recycling      = 'Global recycling coefficient used with neasser > 0; also used for recycling coefficient at divertor target in two points model, if Recycling_target parameter is equal to 0';
    section.Recycling   = 'Density';

    valeur.eta_gas_puff    = 1;
    type.eta_gas_puff      = 'real';
    borne.eta_gas_puff     = [0,1];
    defaut.eta_gas_puff    = 1;
    info.eta_gas_puff      = 'Gas puff fuelling efficiency (allows to use directly measurements from gas injection)';
    section.eta_gas_puff   = 'Density';
    mode.eta_gas_puff      = 'advanced';

    valeur.natural    = 0;
    type.natural      = 'integer';
    borne.natural     = {0,1,2};
    defaut.natural    = 0;
    info.natural      = 'natural density: if = 1, impose density to be higher than natural density;\nif = 2, impose natural density instead of nbar reference';
    section.natural   = 'Density';
    mode.natural      = 'advanced';

    valeur.fnbar_nat   = 1;
    type.fnbar_nat     = 'real';
    borne.fnbar_nat    = [0.1,10];
    defaut.fnbar_nat   = 1;
    info.fnbar_nat     = 'multiplication factor applied to natural density scaling law';
    section.fnbar_nat  = 'Density';
    mode.fnbar_nat     = 'advanced';

    valeur.nea_factor    = 1;
    type.nea_factor      = 'real';
    borne.nea_factor     = [-0.7 10];
    defaut.nea_factor    = 1;
    info.nea_factor      = 'edge density: multiplication factor applied to edge density scaling law:\nif > 0, ne_edge =  nea_factor * LCFS_denstity_scaling_law;\nif < 0,  ne_edge =  abs(nea_factor) * n_bar';
    section.nea_factor   = 'Density';

    valeur.nea_model    = 'Mahdavi';
    type.nea_model      = 'string';
    borne.nea_model     = {'Mahdavi','Eich','fixed ratio'};
    defaut.nea_model    = 'Mahdavi';
    info.nea_model      = 'LCFS density model for H-mode diverted plasma:\n Mahdavi model (default) or\nEich model based on critical density limit due to ballooning mode;\nfixed ratio: min(n_Greenwald,n_bar)/3';
    section.nea_model   = 'Density';
    mode.nea_model      = 'advanced';

    valeur.ftaup    = 1;
    type.ftaup      = 'real';
    borne.ftaup     = [0.1 10];
    defaut.ftaup    = 1;
    info.ftaup      = 'factor multiplying particle confinement time (taup)';
    section.ftaup   = 'Density';
    mode.ftaup      = 'advanced';

    valeur.fn0a    = 1;
    type.fn0a      = 'real';
    borne.fn0a     = [0.1 10];
    defaut.fn0a    = 1;
    info.fn0a      = 'cold neutral source: factor multiplying the core plasma source of neutral for limited plasma (allows to choose the fraction that goes directly in core plasma and in SOL)';
    section.fn0a   = 'Density';
    mode.fn0a      = 'advanced';

    valeur.fn0a_div    = 0.1;
    type.fn0a_div      = 'real';
    borne.fn0a_div     = [0 10];
    defaut.fn0a_div    = 0.1;
    info.fn0a_div      = 'cold neutral source: factor multiplying the core plasma source of neutral for diverted plasma (allows to choose the fraction that goes directly in core plasma and in SOL)';
    section.fn0a_div   = 'Density';
    mode.fn0a_div      = 'advanced';


    valeur.ane    = 0;
    type.ane      = 'integer';
    borne.ane     = {0,1,2,3,4,5,10,11,12};
    defaut.ane    = 0;
    info.ane       = 'density  peaking: (central density / volume average density) choice:\n0 -> f(nsat / nbar) in L-mode, where nsat is the saturation density (LOC/SOC) and f(n_Gr / nbar) in H-mode;\n1 -> flat profile;\n2 -> peaking factor function of li;\n3 -> peaking factor function of collisionality scaling law;\n4 -> fixed value given by parameter vane;\n5 -> proportional to Ti;\n10 -> C. Angioni formula 5 NF 2007 (depends on Greenwald fraction, NBI fuelling and R0);\n11 -> Angioni 2007 formula 5 in H mode and new fit scaling law from Lmode data base for L-mode;\n12 -> SPARC scaling similar to Angioni 2007 (J. Plasma Phys. (2020), vol. 86, 865860502)  in H mode and new fit scaling law from Lmode data base for L-mode';
    section.ane = 'Density';
    section.ane = 'Density';

    valeur.vane    = 1;
    type.vane      = 'real';
    borne.vane     = [0.5,5];
    defaut.vane    = 1;
    info.vane      = 'if ane = 4, value of density peaking factor';
    section.vane   = 'Density';
    
    valeur.ne_shape    = 'Auto';
    type.ne_shape      = 'string';
    borne.ne_shape     = {'Auto','Lmode','Hmode'};
    defaut.ne_shape    = 'Auto';
    info.ne_shape      = 'method used to compute density profile shape: \nAuto -> depending of L or H mode; \nHmode -> force Hmode shape also in Lmode; \nLmode -> force Lmode shape also in Hmode';
    section.ne_shape   = 'Density';
    mode.ne_shape      = 'advanced';
    
    valeur.ne_free    = 3;
    type.ne_free      = 'integer';
    borne.ne_free     = {0,3,4};
    defaut.ne_free    = 3;
    info.ne_free      = 'number of parameters to define density profile in H mode (with pedestal, but depending on value of parameter ne_shape):\nif = 3, central density, edge density and peaking factor;\nif = 4,central density, edge density, peaking factor and pedestal density;\nif=0, TEP model: formule 6.5 in Phd of Alexey Zabolotskiy p 121, up to the top of the pedestal';
    section.ne_free   = 'Density';
    
    valeur.neped_expo    = - 0.7;
    type.neped_expo      = 'real';
    borne.neped_expo     = [-2.5,0];
    defaut.neped_expo    = - 0.7;
    info.neped_expo      = 'density at pedestal top: exponent in scaling law for pedestal density:\nne_ped = ne_a * (ne_a / ne_Gr)^neped_expo';
    section.neped_expo   = 'Density';
    
    valeur.pix   = 0.7;
    type.pix      = 'real';
    borne.pix     = [0,1];
    defaut.pix    = 0.7;
    info.pix      = 'position of maximum of pellet deposition profile';
    section.pix   = 'Pellet';
    
    valeur.piw   = 0;
    type.piw      = 'real';
    borne.piw     = [0,3];
    defaut.piw    = 0;
    info.piw      = 'width of pellet deposition profile (gaussian); if = 0, NGS model is used to compute the shape of deposition';
    section.piw   = 'Pellet';
    
    valeur.pif    = 0;
    type.pif      = 'real';
    borne.pif     = [0,1];
    defaut.pif    = 0;
    info.pif      = 'fraction of fuelling due to pellet injection; if = 1, automatic detection of pellet injection (detection of peaks in nbar waveform)';
    section.pif   = 'Pellet';

    valeur.fpolarized    = 0;
    type.fpolarized      = 'real';
    borne.fpolarized     = [0,1];
    defaut.fpolarized    = 0;
    info.fpolarized      = 'effective fraction of material in pellet that is polarized and then enhance fusion reactivity by a factor 1.5 (reference: L. Baylor N.F. 2023 https://doi.org/10.1088/1741-4326/acc3ae)';
    section.fpolarized   = 'Pellet';

    valeur.scaling    = 0;
    type.scaling      = 'integer';
    borne.scaling     = {0,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21};
    defaut.scaling    = 0;
    info.scaling      = sprintf('%s','choice of scaling law for energy confinement time:\n0 = ITERL-96P(th) + ITERH-98P(y,2);\n1 = Ohmic (to be used for startup phase);\n', ...
                                     '2 = ITPA 2 terms;\n3 = DS03 (no beta dependence);\n 4 = adjusted to match experimental value of Wdia;\n5 = scaling ITER EIV,Std;\n6 = Ohmic scaling in Tokamak Wesson;\n', ...
                                     '7 = ITERH-98P(y,2)/2 L-mode + ITERH-98P(y,2) H-mode;\n8 = user defined scaling (as matlab function);\n9 =  J. Garcia PRL J_pol = 0 (for hybrid scenario);\n', ...
                                     '10 = Sauter & Martin in H-mode + quasi analytical in L-mode;\n11 = Sauter & Martin in H-mode + Elbeze EIV 2005 L-mode;\n', ...
                                     '12 = same as 11 with limitation on beta_N for burning plasma with high radiative fraction;\n13 = reserved;\n', ...
                                     '14 = Robust scaling H1 and L1 (A.Murari, NF 57,2017,120617); must be used with ploss_exp= ''no_prad'';\n',...
                                     '15 = Robust scaling H1 and L1 (A.Murari, NF 57,2017,120617) with limitation on beta_N; must be used with ploss_exp= ''no_prad'';\n', ...
                                     '16 = ITER89-P + ITERH-98P(y,2);\n17 = ITER89-P in Lmode and ITER89-P + Cordey tau_{pedestal} in H mode;\n', ...
                                     '18 = ITER89-P in L-mode and Petty 2008 gyroBohm scaling in H-mode\n', ...
                                     '19 = ITER89-P in L-mode and  ITER89-P  + W_ped scaling incorporating experimental data and prediction from MHD code for ITER and DEMO\n', ...
                                     '20 = must not be used (preprint version of ITPA20);\n21= ITPA20 scaling for ITER 2020 with triangularity dependance');
    section.scaling   = 'Confinement & Transport';

    valeur.dilution      = 1;
    type.dilution        = 'integer';
    borne.dilution       = {0,1};
    defaut.dilution      = 1;
    info.dilution        = 'dilution effect on energy confinement:\nif = 0, no dilution effect;\nif = 1, take into account the ions density dilution on plasma energy contents (W_new = (1 + ni/ne) * W_scaling)';
    section.dilution     = 'Confinement & Transport';
    mode.dilution        = 'advanced';
 
    valeur.tau_limitation    = 'On';
    type.tau_limitation      = 'string';
    borne.tau_limitation     = {'On','Off','Saturate'};
    defaut.tau_limitation    = 'On';
    info.tau_limitation      = 'SOC / LOC transition:\nif = Off, no limitation of confinement time;\nif = On, the confinement time is the minimum of confinement time given by Neo-Alcator scaling and L-mode scaling law (useful for ramp up or at low density);\nif = Saturate, the confinement time is the minimum of confinement time given by SOC / LOC transition point and L-mode scaling law';
    section.tau_limitation   = 'Confinement & Transport';
    mode.tau_limitation      = 'advanced';

    valeur.ploss_exp      = 'with_prad';
    type.ploss_exp        = 'string';
    borne.ploss_exp       = {'with_prad','no_prad','max_power','max(pel)+max(pion)'};
    defaut.ploss_exp      = 'with_prad';
    info.ploss_exp        = 'method used to compute ploss :\nif = with_prad, radiative power is subtracted from input power (pin), including brem, cyclo and fraction (fprad) of line radiation;\n if = no_prad, ploss = pin;\n if = max_power, ploss is the maximum of volume integrated total source power (Qe + Qi) as a function of radial position;\nif = max(pel)+max(pion), ploss is the maximum of volume integrated electron source power Qe  + the maximum of volume integrated ion source power Qi';
    section.ploss_exp     = 'Confinement & Transport';

    valeur.fprad      = 1/3;
    type.fprad        = 'real';
    borne.fprad       = [0,1];
    defaut.fprad      = 1/3;
    info.fprad        = 'fraction of line radiative power (core plasma) substracted from input power to compute Ploss (Ploss = Pin - fprad * Pline)';
    section.fprad     = 'Confinement & Transport';
   
    valeur.HH_delta  = 0;
    type.HH_delta    = 'real';
    borne.HH_delta   = [-1,1];
    defaut.HH_delta  = 0;
    info.HH_delta    = 'if # 0, add triangularity effect on energy confinement: W = ((1+delta)/(1+delta0)) ^ HH_delta * W_scaling,\nwhere delta0 is the neutral triangularity (delta0 = 0.35) and recommended value for HH_delta = -0.35\n(https://doi.org/10.1088/0029-5515/39/11y/321);\nremark: delta = z0dinput.geo.d';
    section.HH_delta = 'Confinement & Transport';
    mode.HH_delta    = 'advanced';

    valeur.HH_li     = 0;
    type.HH_li       = 'real';
    borne.HH_li      = [0,3];
    defaut.HH_li     = 0;
    info.HH_li       = 'if > 0, li variation effect on plasma energy content W = (li/HH_li) ^ (2/3) * W_scaling';
    section.HH_li    = 'Confinement & Transport';
    mode.HH_li       = 'advanced';

    valeur.HH_gas_puff     = 0;
    type.HH_gas_puff       = 'real';
    borne.HH_gas_puff      = [0,100];
    defaut.HH_gas_puff     = 0;
    info.HH_gas_puff       = 'if > 0, energy confinement reduction with gas puff fuelling: W = (1 -tanh(HH_gas_puff .* P_ioniz ./ P_in)) * W_scaling';
    section.HH_gas_puff    = 'Confinement & Transport';
    mode.HH_gas_puff       = 'advanced';

    valeur.l2hscaling    = 0;
    type.l2hscaling      = 'real';
    borne.l2hscaling     = {-10,-5,-3,-2,-1.5,-1,0,1,2,3,4,5,6,10,28,30,31};
    defaut.l2hscaling    = 0;
    info.l2hscaling      = sprintf('%s\n','L to H power threshold scaling law:', ...'
                                              'if = 0  -> LH99(1);', ...
                                              'if = 1  -> LH2002;', ...
                                              'if = 2  -> LH2002 + Zeff;', ...
                                              'if = 3  -> YR Martin 2008;', ...
                                              'if = 4  -> NLM-7 Murari 2012;', ...
                                              'if = 5  -> NLM-11 Murari 2012;', ...
                                              'if = 6  -> Jpol change of sign in edge region (E. R. Solano rule);', ...
                                              'if = 10 -> Multimachine scaling law from Murari 2013 (BUEMS);', ...
                                              'if = 28 -> Low density case - Ryter et al, NF 54 (2014) 083003, equation 4', ...
                                              'if = 30 -> Fit of metalic tokamaks database for horizontal targets (E . Delabie et al, 2025 ?);',...
                                              'if = 31 -> Fit of metalic tokamaks database for vertical targets and corner configuration (E . Delabie et al, 2025 ?);',...
                                              'if < 0, criterion based on plasma rotation (abs(l2hscaling) = value of Gamma_ExB / Gamma_ITG for transition)');
    section.l2hscaling   = 'H mode transition';

    valeur.modeh    = 1;
    type.modeh      = 'integer';
    borne.modeh     = {0,1,2};  
    defaut.modeh    = 1;                
    info.modeh      = 'L-Mode to H-mode allowed transition:\n0 -> force L-mode;\n1-> L-Mode to H-mode transition allowed;\n2 -> force H-mode';
    section.modeh   = 'H mode transition';

    valeur.l2hmul    = 0;
    type.l2hmul      = 'real';                    
    borne.l2hmul     = [-10,100];  
    defaut.l2hmul    = 0;
    info.l2hmul      = 'offset added to the threshold power for the transition L-> H (MW)';
    label.l2hmul     = 'L-H offset';
    section.l2hmul   = 'H mode transition';
    mode.l2hmul      = 'advanced';
    
    valeur.plhthr    = 'pel+pion';
    type.plhthr      = 'string';                    
    borne.plhthr     = {'pel+pion','2*pion','P_LCFS','P_LCFS_dwdt'};  
    defaut.plhthr    = 'pel+pion';
    info.plhthr      = 'Power compared to scaling for L to H transition, either:\nPloss as defined for scaling law  (pel + pion - fraction of prad);\ntwice ion  power (2*pion);\npower conducted to LCFS without dWdt term(P_LCFS);\npower conducted to LCFS with dWdt term(P_LCFS_dwdt)';	
    section.plhthr   = 'H mode transition';
    mode.plhthr      = 'advanced';

    valeur.l2hslope    = 0;
    type.l2hslope      = 'real';                    
    borne.l2hslope     = [-1,1] ;  
    defaut.l2hslope    = 0;
    info.l2hslope = 'slope of linear transition between tau_L and tau_H controlled by difference between conducted and threshold power;\n if = 0, on / off transition;\nif < 0 , additionally decrease of confinement when density is close to Greenwald limit';	
    section.l2hslope   = 'H mode transition';
    mode.l2hslope      = 'advanced';

    valeur.fpl2h_lim    = 1;
    type.fpl2h_lim      = 'real';                    
    borne.fpl2h_lim     = [0.5,10] ;  
    defaut.fpl2h_lim    = 1;
    info.fpl2h_lim      = 'factor applied to L-> H scaling power threshold in limiter configuration';	
    section.fpl2h_lim   = 'H mode transition';
    mode.fpl2h_lim      = 'advanced';

    valeur.pl2h_mass_charge   = 0;
    type.pl2h_mass_charge     = 'integer';                    
    borne.pl2h_mass_charge    = {0,1} ;  
    defaut.pl2h_mass_charge   = 0;
    info.pl2h_mass_charge     = 'L to H transition:\nif = 0, scaling law as defined by l2hscaling;\nif = 1, adds dependences on mass and charge of main ion (R. Behn et all, PPCF 2015)';	
    section.pl2h_mass_charge  = 'H mode transition';
    mode.pl2h_mass_charge     = 'advanced';

    valeur.hysteresis    = 1;
    type.hysteresis      = 'real';                    
    borne.hysteresis     = [0,1];  
    defaut.hysteresis    = 1;
    info.hysteresis      = 'Control of histeresis for the back transition H-> L mode:\nP_{H->L} = hysteresis * P_{in} + (1 - hysteresis) * P_{lh,thr}\n';	
    section.hysteresis   = 'H mode transition';
    mode.hysteresis      = 'advanced';

    valeur.ton_modeh    = Inf;
    type.ton_modeh      = 'real';                    
    borne.ton_modeh     = [-Inf,Inf];  
    defaut.ton_modeh    = Inf;
    info.ton_modeh      = sprintf('if it is defined, start time for the H mode phase: transition to mode H cannot start earlier than ton_modeh;\n(undefined value = Inf; mode H state between ton_modeh and toff_modeh is controlled by the parameter modeh)');	
    section.ton_modeh   = 'H mode transition';

    valeur.toff_modeh    = Inf;
    type.toff_modeh      = 'real';                    
    borne.toff_modeh     = [-Inf,Inf];  
    defaut.toff_modeh    = Inf;
    info.toff_modeh      = sprintf('if it is defined, end time for the H mode phase: mode H  phase cannot continue latter than toff_modeh\n(undefined value = Inf; mode H state between ton_modeh and toff_modeh is controlled by the parameter modeh)');	
    section.toff_modeh   = 'H mode transition';

    valeur.fpped    = 1;
    type.fpped      = 'real';
    borne.fpped     = [-100,10];
    defaut.fpped    = 1;
    info.fpped      = 'if > 0, pedestal pressure multiplier (pressure deduced from scaling law);\n if = 0, switch off limitation  of pedestal pressure due to MHD and experimental limit;\n if < 0, stiff model is activated (in this case, pedestal pressure multiplier is abs(fpped))';
    section.fpped   = 'Confinement & Transport';

    valeur.hmore_pped    = 0;
    type.hmore_pped      = 'integer';
    borne.hmore_pped     = {0,1,2,3,4,5,6};
    defaut.hmore_pped    = 0;
    info.hmore_pped      = sprintf('%s','if = 0, no effect;\nif = 1, pedestal pressure is multiplied by the H factor waveform: pped_use = H_factor * pped_predicted;\n', ...
                                         'if = 2, pedestal pressure is multiplied by the H factor waveform when H factor is less than 1: pped_use = min(1,H_factor) * pped_predicted;\n', ...
                                         'if = 3, pedestal pressure is multiplied by the H factor waveform when H factor is greater than 1: pped_use = max(1,H_factor) * pped_predicted;\n', ...
                                         'if = 4, as option 1, but maximum pedestal pressure is not affected by fpped;\nif = 5, as option 2, but maximum pedestal pressure is not affected by fpped;\n',...
                                         'if = 6, as option 0, but maximum pedestal pressure is not affected by fpped');
    section.hmore_pped   = 'Confinement & Transport';
    mode.hmore_pped       = 'advanced';

    valeur.ode_pped    = 0;
    type.ode_pped      = 'integer';
    borne.ode_pped     = {0,1};
    defaut.ode_pped    = 0;
    info.ode_pped      = 'if = 0, no effect;\nif = 1, pressure at the top of pedestal evolves in time with pedestal confinement time:\nd(3/2*Vp*Pped)/dt = - (3/2*Vp*Pped)/ tau_ped + Power_LCFS with tau_ped = W_ped_steady_state / Power_LCFS';
    section.ode_pped   = 'Confinement & Transport';
    mode.ode_pped       = 'advanced';

    valeur.fstiff    = 1;
    type.fstiff      = 'real';
    borne.fstiff     = [-10,10];
    defaut.fstiff    = 1;
    info.fstiff      = 'stiff transport model: when stiff model is activated (fpped < 0),abs(fstiff) gives the temperature gradient in the core in eV per electron banana width; default value = 1;\n if fstiff<0, use ftrap*rho_banana+(1-ftrap)*rho_larmor instead of rho_banana;\n if "alpha" or "alpha+neo" shape is selected, then it becomes the multiplication factor of pressure gradient provided by the s_alpha limit formula.';
    section.fstiff   = 'Confinement & Transport';
    mode.fstiff      = 'advanced';

    valeur.isotope_stiff    = 1;
    type.isotope_stiff      = 'real';
    borne.isotope_stiff     = [-2,2];
    defaut.isotope_stiff    = 1;
    info.isotope_stiff      = 'if stiff transport is activated and fstiff < 0, isotope_stiff is the exponent of masse dependance on ITG part (fraction of trapped particules):\n confinement time is proportional to (meff/2)^isotope_stiff.';
    section.isotope_stiff   = 'Confinement & Transport';
    mode.isotope_stiff      = 'advanced';

    valeur.usepped_scl    = 0;
    type.usepped_scl      = 'integer';
    borne.usepped_scl     = {0,1,2,3,4,5};
    defaut.usepped_scl    = 0;
    info.usepped_scl      = 'if = 0, pedestal energy content is the difference between H-mode and L-mode energy content or half of this difference (see scaling options);\n if = 1, tau_ped scaling law (ITPA McDonald) is used to compute Pped (does not work with stiff model);\nif = 2, minimum between standard rule (Pped = K (W_Hmode - W_Lmode) and scaling law prediction is used;\nif = 3, P_ped scaling incorporating experimental data and prediction from MHD code for ITER and DEMO is used\nif = 4, minimum between standard rule (Pped = K (W_Hmode - W_Lmode) and scaling incorporating prediction from MHD code for ITER and DEMO is used\nif = 5, use model from reference F.D. Halpern et al, PoP 15 (2008) p 062505';
    section.usepped_scl   = 'Confinement & Transport';

    valeur.taurotmul   = 0;
    type.taurotmul     = 'real';
    borne.taurotmul    = [0,10];
    defaut.taurotmul   = 0;
    info.taurotmul     = 'toroidal rotation momentum confinement time: multiplication factor applied to energy confinement time to obtain toroidal rotation momentum confinement time (tau_rotation = taurotmul * taue);\n if = 0, toroidal rotation momentum confinement time = ion heat confinement time';
    section.taurotmul  = 'Rotation';
    mode.taurotmul     = 'advanced';
   
    valeur.rotation_scl    = 'Rice';
    type.rotation_scl      = 'string';
    borne.rotation_scl     = {'Rice','Barnes & Parra','deGrassie'};
    defaut.rotation_scl    = 'Rice';
    info.rotation_scl      = 'intrinsic rotation scaling: switch between differents scaling for spontaneous (intrinsic) rotation';
    section.rotation_scl   = 'Rotation';
    mode.rotation_scl      = 'advanced';
  
    valeur.fintrinsic   = 1;
    type.fintrinsic     = 'real';
    borne.fintrinsic    = [0,10];
    defaut.fintrinsic   = 1;
    info.fintrinsic     = 'multiplication factor applied to intrinsic rotation scaling.The constant in the scaling in not determined;\n Suggested value with Rice scaling is between 0.15-0.25;\n Suggested value with deGrassie scaling is about 1;\n if = 0,a factor depending on collisionality is used, computed from  J. C. Hillesheim et al, arxiv:1407.2121v1 physics.plasma-ph 8 Jul 2014';
    section.fintrinsic  = 'Rotation';
    mode.fintrinsic     = 'advanced';
   
    valeur.omega_shape  = 0;
    type.omega_shape     = 'integer';
    borne.omega_shape    = {0,1,2,3,4,5,6};
    defaut.omega_shape   = 0;
    info.omega_shape     = 'shape of rotation profile: if = 0, proportional to Ti;\nif = 1, given by scaling law from H. Weisen et al paper (NF 52 2012 p 042001);\nif = 2, proportionnal to Pion;\nif = 3, proportional to Ptot;\nif = 4, proportional to n_ion;\nif = 5, proportional local value of deGrassie scaling;\nif = 6, interpolation between Ti shape and deGrassie shape depending of relative weight of NBI source on rotation compared to intrinsic rotation';
    section.omega_shape  = 'Rotation';
    mode.omega_shape     = 'advanced';
   
    valeur.solid_rotation   = 1;
    type.solid_rotation     = 'integer';
    borne.solid_rotation    = {-2,-1,0,1};
    defaut.solid_rotation   = 1;
    info.solid_rotation     = 'In radial electric field computation:\n if = 1, assume plasma solid rotation; i.e. neglect term in poloidal rotation;\nif = 0, take into account term  in poloidal rotation (main ion contribution from Kim model);\n if = -1, take into account term in poloidal rotation (simplified formulation with k_neo * grad(T) / eB formulation for each species);\n if = -2,take into account term in poloidal rotation (simplified formulation with k_neo * grad(T) / eB formulation for each species) and additionnal term for poloidal rotation induced by NBI';
    section.solid_rotation  = 'Rotation';
    mode.solid_rotation     = 'advanced';
    
    valeur.xiioxie    = 2;
    type.xiioxie      = 'real';
    borne.xiioxie     = [-10,10];
    defaut.xiioxie    = 2;
    info.xiioxie      = 'ratio  Xii over Xie (more precisely: (ni*Xii) / (ne*Xie));\n if = 0 -> compute from ITG / TEM stability diagram; \n if > 0, assigned value of xii // xie; \n if < 0 , xii // xie is computed using the critical gradient model in which the stiffness parameter is given by xiioxie;\n in this case the theoretical value is 4.5';
    section.xiioxie   = 'Confinement & Transport';

    valeur.xiioxie_ped    = 0;
    type.xiioxie_ped      = 'real';
    borne.xiioxie_ped    = [0,10];
    defaut.xiioxie_ped    = 0;
    info.xiioxie_ped      = 'ratio Xii over Xie for edge an pedestal (more precisely: (ni*Xii_ped) / (ne*Xie_ped));\n if = 0,the value for the core plasma is used (xiioxie)';
    section.xiioxie_ped   = 'Confinement & Transport';
    mode.xiioxie_ped     = 'advanced';

    valeur.kishape    = 3;
    type.kishape      = 'real';
    borne.kishape     = [-10,10];
    defaut.kishape    = 3;
    info.kishape      = 'radial shape of heat tranport cofficient:\nif > 0, Kappa = C  ( 1 + kishape * x ^ ki_expo);\nif = 0, Kappa = C * model_based_shape (see coef_shape);\nif < 0 , Kappa =  q^abs(kishape)';
    section.kishape   = 'Confinement & Transport';

    valeur.ki_expo    = 2;
    type.ki_expo      = 'real';
    borne.ki_expo     = [-3,20];
    defaut.ki_expo    = 2;
    info.ki_expo      = 'if > 0, exponent of Kappa shape fonction;\nif < 0, the shape is Kappa = C*(x-2 / 3*x^2+abs(ki_expo) / x / 30+kishape*x^20);\n in this case ki_expo controls the centre and kishape controls the edge;\nthis gives a rather linear temperature profile in the gradient zone';
    section.ki_expo   = 'Confinement & Transport';

    valeur.xieorkie   = 1;
    type.xieorkie     = 'integer';
    borne.xieorkie    = {0,1};
    defaut.xieorkie   = 1;
    info.xieorkie     = 'if = 0, radial shape of heat tranport cofficient is Kappa;\nif = 1, instead of given Kappa shape, Chi shape is fixed (Kappa ~ Ne * Chi)';
    section.xieorkie  = 'Confinement & Transport';
    mode.xieorkie    = 'advanced';

    valeur.coef_shape   = 'bgb';
    type.coef_shape     = 'integer';
    borne.coef_shape    = {'bgb','cdbm','bgb+neo','cdbm+neo','stiff','stiff+neo','stiff_limited','stiff_limited+neo','alpha','alpha+neo'};
    defaut.coef_shape   = 'bgb';
    info.coef_shape     = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s','Shape of transport coefficient (Kappa) when kishape = 0:', ...
                                                   '1/ Bohm-gyroBohm = Bohm-gyro Bohm shape,',... 
                                                   '2/ CDBM model = CDBM model shape,', ...
                                                   '3/ Stiff model associated to Stiff confinment,', ...
                                                   '4/ Alpha is associated with model where transport is provided by limitation of pressure gradient due to ballooning limit from s-alpah diagram,', ...
                                                   '5/ Stiff_limited correspond to the stiff model limited by the alpha limit;', ...
                                                   'Remarks:', ...
                                                   '1/ in stiff and aplha case, parameter xieorkie is ignored', ...
                                                   '2/ if +neo, then neoclasical transport coeficient for ions, computed with Hinton model, is added turbulent tranport model (both to electron and ion channel).');
    section.coef_shape  = 'Confinement & Transport';

    valeur.exp_shape   = 0;
    type.exp_shape     = 'integer';
    borne.exp_shape    = {0,1};
    defaut.exp_shape   = 0;
    info.exp_shape     = 'if = 1, when external Te and Ti are provided, keep only the shape and rescale Te and Ti to get W_th computed in  METIS';
    section.exp_shape  = 'Confinement & Transport';
    mode.exp_shape    = 'advanced';

    valeur.grad_ped   = 3;
    type.grad_ped     = 'integer';
    borne.grad_ped    = {0,1,2,3};
    defaut.grad_ped   = 3;
    info.grad_ped     = 'method to compute the pressure gradient at the top of pedestal:\nif = 0, use the method that does not take into account the discontinuity in gradient;\nif = 1, take into account the discontinuity in gradient;\n if = 2, improved integration method for equation is used;\n if = 3, improved integration method for equation is used and no limitation in flux;';
    section.grad_ped  = 'Confinement & Transport';
    mode.grad_ped    = 'advanced';

    valeur.adiabatic   = 0;
    type.adiabatic     = 'integer';
    borne.adiabatic    = {0,1};
    defaut.adiabatic   = 0;
    info.adiabatic     = 'if = 0, turn off adiabatic compression term in the ODE for energy content evolution;\n if = 1, turn on adiabatic compression term in the ODE for energy content evolution; this term is generally negligible but not during breakdown, fast build-up of plasma volume or disruption;\n Warning:this term can add some numerical noise in simulation';
    section.adiabatic  = 'Confinement & Transport';
    mode.adiabatic     = 'advanced';

    valeur.hollow   = 0;
    type.hollow     = 'integer';
    borne.hollow    = {0,1};
    defaut.hollow   = 0;
    info.hollow     = 'if = 0, does not allow hollow temperature profiles;\n if = 1, allows hollow temperature profiles;\n hollow profiles may exist in presence of heavy impurities (i.e W or Mo) accumulation in core plasma.';
    section.hollow  = 'Confinement & Transport';
    mode.hollow     = 'advanced';

    valeur.disrup    = 0;
    type.disrup      = 'integer';
    borne.disrup     = {0,1,2,3};
    defaut.disrup    = 0;
    info.disrup      = 'Radiative limit induced disruption:\nif = 0, no effect on energy confinement time when radiative power exceed input power;\nif = 1, reduction of confinement time at time slices where radiative power exceed input power;\nif = 2, reduction of confinement time since first time where radiative power exceed input power;\nmust be used with parameter ploss_exp set to max_power or to max(pel)+max(pion);\nif = 3, reduce confiment accordingly to fraction of volume where enegy flux become negative.';
    section.disrup   = 'Confinement & Transport';
    mode.disrup      = 'advanced';

    valeur.collapse    = 0;
    type.collapse      = 'integer';
    borne.collapse     = {0,1};
    defaut.collapse    = 0;
    info.collapse      = 'Radiative colapse:\nif = 1, allows progressive radiative collapse to squeeze radialy the plasma pressure.';
    section.collapse   = 'Confinement & Transport';
    mode.collapse      = 'advanced';

    valeur.te_max    = 1e5;
    type.te_max     = 'real';
    borne.te_max     = [1e4,1e6];
    defaut.te_max    = 1e5;
    info.te_max     = sprintf('%s\n%s\n%s\n%s\n%s','Maximum allowed internal electron temperature (and ion temperature) in METIS solver (in eV).', ...
                              'The default value is 1e5 eV corresponding to the maximum tabulated temperature for thermal cross section and radiative cooling rate.', ...
                              'The upper limit allowed by this parameter is set just below the energy threshold for pair prodcution.', ...
                              'This limit should be increase for some aneutronic fusion reaction and accordingly models for cross section, radiative cooling rate,', ...
                              'relativist bremsstrahlung and enhanced collisional heat exchange between eletrons and ions should be selected');                      
    section.te_max  = 'Confinement & Transport';
    mode.te_max     = 'advanced';
    
    valeur.extended_qei    = 'off';
    type.extended_qei      = 'string';
    borne.extended_qei     = {'on','off'};
    defaut.extended_qei    = 'off';
    info.extended_qei       = sprintf('%s\n%s\n%s\n%s\n%s','if = on, the collisional heat exchange term between eletrons and ions is computed', ...
                                      'with the formula including relastivistic correction and large T_i/T_e effect.', ...
                                      'reference: Modification of classical Spitzer ionelectron energy transfer rate', ...
                                      'for large ratios of ion to electron temperatures, T.H. Rider and P.J. Catto,', ...
                                      'Phys. Plasmas 2, 18731885 (1995), https://doi.org/10.1063/1.871274');                      
    section.extended_qei  = 'Confinement & Transport';
    mode.extended_qei    = 'advanced';

    valeur.qdds    = 0;
    type.qdds      = 'real';
    borne.qdds     = [-2,3];
    defaut.qdds    = 0;
    info.qdds       = 'Sawtooth model:\nif > 0,time averaged effect with clamping of safety factor at q_st value inside q_st radius;\nif = 0, no effect;\nif < 0, time resolved sawteeth, triggered when q0 <= abs(q_st)';
    label.qdds     = 'q_st';
    type.qdds      = 'real';
    section.qdds   = 'MHD & ITB';

    valeur.ddsmode         = 0;   % mode 0= Porcelli , 1 = Kadomtsev, 2 = Kadomtsev incomplet
    type.ddsmode           = 'integer';                % type entier
    borne.ddsmode          = {0,1,2,3};                  % valeurs possible 
    defaut.ddsmode         = 0;                        % valeurs par defaut
    info.ddsmode           = 'Sawthooth reconnection:\n0 = simple clamping;\n1 = Porcelli;\n2 = Kadomtsev;\n3 = partial Kadomtsev';
    section.ddsmode        = 'MHD & ITB';
    label.ddsmode          = 'q_st_mode';   % mode 0= Porcelli , 1 = Kadomtsev, 2 = Kadomtsev incomplet
    mode.ddsmode           = 'advanced';

    valeur.q0_dds_trig    = 0;
    type.q0_dds_trig      = 'real';
    borne.q0_dds_trig     = [0,2];
    defaut.q0_dds_trig    = 0;
    info.q0_dds_trig      = 'value of q0 triggering a sawtooth independently of others conditions (with qdds < 0)';
    section.q0_dds_trig   = 'MHD & ITB';
    mode.q0_dds_trig      = 'advanced';

    valeur.s1crit    = 0;
    type.s1crit      = 'real';
    borne.s1crit     = [0,2];
    defaut.s1crit    = 0;
    info.s1crit      = 'critical shear for sawtooth triggering: if = 0, use criterium on q_st,\notherwise trigger a ST when magnetic shear @ q=1 is above s1crit;\nq_st must be < 0 to activate this mechanism (Porcelli PPCF 1996)';
    section.s1crit   = 'MHD & ITB';
    mode.s1crit      = 'advanced';

    valeur.betap1crit    = 0;
    type.betap1crit      = 'real';
    borne.betap1crit     = [0,2];
    defaut.betap1crit    = 0;
    info.betap1crit      = 'critical betap for sawtooth triggering: if = 0, use criterium on q_st,\notherwise trigger a ST when magnetic betap1 @ q=1 is above betap1crit;\nq_st must be < 0 to activate this mechanism (Jardin PoP 2020)';
    section.betap1crit   = 'MHD & ITB';
    mode.betap1crit      = 'advanced';

    valeur.w1         = 0.5;   % fraction de la zone reconnectee (en fraction de rmix)
    type.w1           = 'float';                % type reel
    borne.w1          = [0.1,1];               % valeurs possible 
    defaut.w1         = 0.1;                      % valeurs par defaut
    info.w1           = 'Sawtooth reconnection: width of null magnetic shear zone in unit of psistar(q=1) [0.1,1] {0.5} (Porcelli PPCF 1996)';
    section.w1        = 'MHD & ITB';
    mode.w1           = 'advanced';
	
    valeur.epsq         = 1e-2;   % pente de q dans la zone q = 1
    type.epsq           = 'float';                % type reel
    borne.epsq          = [0,0.1];               % valeurs possible 
    defaut.epsq         = 1e-2;                      % valeurs par defaut
    info.epsq           = 'Sawtooth reconnection: slope of q inside mixing radius (Porcelli PPCF 1996)';
    section.epsq        = 'MHD & ITB';
    mode.epsq           = 'advanced';

    valeur.kidds    = 1;
    type.kidds      = 'real';
    borne.kidds     = [1,1e3];
    defaut.kidds    = 1;
    info.kidds      = 'Sawtooth model: transport multiplicator inside q <= q_st flux surface';
    label.kidds     = 'Chi_st';
    section.kidds   = 'MHD & ITB';

    valeur.peeling    = 0;
    type.peeling      = 'integer';
    borne.peeling     = {0,1,2};
    defaut.peeling    = 0;
    info.peeling      = 'ELMs type:\nif = 0, ballooning limit only;\nif = 1, ballooning and peeling limits;\nif = 2, peeling limit only (use as peeling limit the proxy <j>_top_pedestal > (Ip / Splasma))';
    section.peeling   = 'MHD & ITB';
    mode.peeling      = 'advanced';

    valeur.dwow_elm    = 0;
    type.dwow_elm      = 'real';
    borne.dwow_elm     = [-10,1];
    defaut.dwow_elm    = 0;
    info.dwow_elm      = 'ELMs model: if > 0, fraction of pedestal energy losses during one ELM (crash  is triggered when pped exceeds ppedmax value);\nif = 0,  no ELM (default mode);\n if = 1, use MHD limit for threshold and energy scaling for energy content after crash;\n if < 0, use scaling law f(nustar) * abs(dwow_elm)';
    section.dwow_elm   = 'MHD & ITB';
    mode.dwow_elm      = 'advanced';

    valeur.tau_elm_factor    = 10;
    type.tau_elm_factor      = 'real';
    borne.tau_elm_factor     = [1,1000];
    defaut.tau_elm_factor    = 10;
    info.tau_elm_factor      = 'ELMs model: factor between turbulent tranport and neoclassical transport in pedestal;\ni.e. multiplicator of energy confinement time used to obtain effective confinement time in H-mode between ELMs crashes \n(taue_etb = tau_elm_factor * taue)';
    section.tau_elm_factor   = 'MHD & ITB';
    mode.tau_elm_factor      = 'advanced';

    valeur.runaway    = 0;
    type.runaway      = 'integer';
    borne.runaway     = {0,1,2,3,4,5};
    defaut.runaway    = 0;
    info.runaway   = 'runaway current model:\nif = 0, no runaway current;\nif = 1, runaways appear if electric field is above critical field limit;\nif = 2, model for runaway including LH waves effect;\nif = 3, include model for LH and delay for runaway acceleration;\nif = 4, as 3 profiles (Te, Ne, Zeff) instead of averaged quantities;\nif = 5, as 4 with additional collisions with neutrals (dedicated to breakdown studies)';
    section.runaway   = 'Current diffusion & Equilibrium';
    mode.runaway      = 'advanced';

    valeur.modeboot    = 5;
    type.modeboot      = 'integer';
    borne.modeboot     = {0,1,2,3,4,5,6};
    defaut.modeboot    = 5;
    label.modeboot     = 'model';
    info.modeboot      = 'bootstrap current model:\nif = 0, scaling law G. T. HOANG;\nif  = 1, Sauter formula;\nif = 2, Sauter formula + asymmetric current;\nif  = 3, Hager & Chang modified Sauter formula;\nif = 4, Hager & Chang modified Sauter formula + asymmetric current;\nif = 5, NEO fit (A; Redl et al, PoP 2021);\nif = 6, NEO fit (A; Redl et al, PoP 2021) + asymmetric current';
    section.modeboot   = 'Bootstrap';

    valeur.bootmul    = 1;
    type.bootmul      = 'real';
    borne.bootmul     = [0,10];
    defaut.bootmul    = 1;
    info.bootmul      = 'multiplication factor applied to bootstrap current';
    section.bootmul   = 'Bootstrap';
    mode.bootmul      = 'advanced';

    valeur.ffit_ped    = 1;
    type.ffit_ped      = 'real';
    borne.ffit_ped     = [0,10];
    defaut.ffit_ped    = 1;
    info.ffit_ped      = 'bootstrap current model: multiplication factor applied to gradients computation in pedestal\nif = 0 use standard 3 points derivative';
    section.ffit_ped   = 'Bootstrap';
    mode.ffit_ped      = 'advanced';

    valeur.fspot    = 0.15;
    type.fspot      = 'real';
    borne.fspot     = [0,1];
    defaut.fspot    = 0.15;
    info.fspot      = 'multiplication factor applied to bootstrap like current due to fast alpha particles (0.05 - 0.15; must be computed with the help of SPOT or other MC codes)';
    section.fspot   = 'Bootstrap';
    mode.fspot      = 'advanced';
    
    valeur.force_spitzer    = 0;
    type.force_spitzer      = 'integer';
    borne.force_spitzer     = {0,1};
    defaut.force_spitzer    = 0;
    info.force_spitzer      = 'if = 1, replace neoclassical resistivity by Spitzer formula where Z = Zeff (use equation 18 of Sauter PoP 1999)';
    section.force_spitzer   = 'Bootstrap';
    mode.force_spitzer      = 'advanced';
    
    valeur.neutral_friction   = 0;
    type.neutral_friction     = 'real';
    borne.neutral_friction    = [0,10];
    defaut.neutral_friction   = 0;
    info.neutral_friction     = 'if > 0, add neutral friction effect on resistivity (V. A. Belyakov et al , PhysCon 2003 Saint Petersburg Russia (IEEE));\n neutral_friction is a mutiplicator factor applied to the formula;\nThis effect is already taken into account if breakdown model is switch on:\nin this case if  neutral_friction = 0, the factor 1 is used.';
    section.neutral_friction  = 'Bootstrap';
    mode.neutral_friction     = 'advanced';
    
    valeur.f_eta_turb         = 0;
    type.f_eta_turb           = 'real';
    borne.f_eta_turb          = [-2.7,2.7]; 
    defaut.f_eta_turb         = 0;
    info.f_eta_turb           = 'if ~= 0, adding (if >0) or sustracting (if < 0) turbulent resistivity computed from reference: L.Colas 1993 Nucl. Fus. 33 156.\nf_eta_turb is a multiplicative factor and D_tild = electron heat diffusivity;\n It is generally assumed that turbulence decrease resistivity in a tokamak.';
    section.f_eta_turb        = 'Bootstrap';
    mode.f_eta_turb           = 'advanced';
  
    valeur.vloop    = 0;
    type.vloop      = 'integer';
    borne.vloop     = {0,1,2,3,4,5,6};  
    defaut.vloop    = 0;
    info.vloop      = 'Current diffusion, boundary condition:\nif = 0, Ip  waveform given;\nif = 1, vloop = 0 as reference, Ip free;\nif = 2, vloop = vref, Ip free;\nif = 3, PLH is computed to follow Ip waveform @ vloop = 0;\n if = 4, poloidal edge flux given by flux waveform\n;if = 5, hybrid boundary condition for the coupling with FREEBIE;\nif = 6, P_NBI is computed to follow Ip waveform @ vloop = 0';
    section.vloop   = 'Current diffusion & Equilibrium';
    mode.vloop      = 'advanced';

    valeur.vref    = 0;
    type.vref      = 'real';
    borne.vref     = [-3,3];
    defaut.vref    = 0;
    info.vref      = 'Current diffusion, boundary condition: reference value  for Vloop;\nSwitch from Ip reference to poloidal flux reference when Vloop <= vref';
    section.vref   = 'Current diffusion & Equilibrium';
    mode.vref      = 'advanced';

    valeur.tswitch    = Inf;
    type.tswitch      = 'real';
    borne.tswitch     = [0,Inf];
    defaut.tswitch    = Inf;
    info.tswitch      = 'Current diffusion, boundary condition: time at which switching from Ip reference to  poloidal flux reference (or vref) is allowed (s)\n (used only with option.vloop > 0, inactive in evolution mode)';
    section.tswitch   = 'Current diffusion & Equilibrium';
    mode.tswitch      = 'advanced';

    valeur.li    = 1;
    type.li      = 'real';
    borne.li     = [0.1,10];
    defaut.li    = 1;
    label.li     = 'initial li';
    info.li      = 'value of normalised internal inductance at the first time (used to compute initial current profile)';
    section.li   = 'Breakdown and burn-through';

    valeur.breakdown    = 0.03;
    type.breakdown      = 'real';
    borne.breakdown     = [-100,1];
    defaut.breakdown    = 0.03;
    label.breakdown     = 'initial Vloop';
    info.breakdown      = 'Value of Vloop at the first time step:\nif > 0, electric field at the breakdown time normalised to the Dreicer electric field;\nif < 0, initial Vloop in Volt per turn';
    section.breakdown   = 'Breakdown and burn-through';

    valeur.berror    = 0;
    type.berror      = 'real';
    borne.berror     = [0,10];
    defaut.berror    = 0;
    info.berror      = 'residual magnetic field at breakdown time (T);\nif = 0, computation taue_breakdown is desactivated';
    section.berror   = 'Breakdown and burn-through';
    mode.berror      = 'advanced';

    valeur.L_eddy    = 0;
    type.L_eddy      = 'real';
    borne.L_eddy     = [0,1];
    defaut.L_eddy    = 0;
    info.L_eddy      = 'characteristic inductance value of passive structure that can carry current during breakdown (H);\nif = 0, use typical value taking into account major radius of the plasma';
    section.L_eddy   = 'Breakdown and burn-through';
    mode.L_eddy      = 'advanced';

    valeur.R_eddy    = 0;
    type.R_eddy      = 'real';
    borne.R_eddy     = [0,10];
    defaut.R_eddy    = 0;
    info.R_eddy      = 'characteristic resistance value of passive structure that can carry current during breakdown (Ohm);\nif = 0, use typical value taking into account major radius of the plasma';
    section.R_eddy   = 'Breakdown and burn-through';
    mode.R_eddy      = 'advanced';

    valeur.C_eddy    = 1;
    type.C_eddy      = 'real';
    borne.C_eddy     = [0,1];
    defaut.C_eddy    = 1;
    info.C_eddy      = 'fraction of eddy current that is removed from plasma current waveform (I_p = I_p_ref - C_eddy * I_eddy)';
    section.C_eddy   = 'Breakdown and burn-through';
    mode.C_eddy      = 'advanced';

    valeur.B_eddy    = 1;
    type.B_eddy      = 'real';
    borne.B_eddy     = [0,10];
    defaut.B_eddy    = 1;
    info.B_eddy      = 'multiplication factor of error magnetic field created by eddy current (B_RorZ = B_eddy * mu0 * I_eddy / R / pi)';
    section.B_eddy   = 'Breakdown and burn-through';
    mode.B_eddy      = 'advanced';

    valeur.PSI_eddy    = 1;
    type.PSI_eddy      = 'real';
    borne.PSI_eddy     = [0,2];
    defaut.PSI_eddy    = 1; 
    info.PSI_eddy      = 'multiplication factor of pertubation due to eddy current on poloidal flux waveform (Flux = Flux_ref - PSI_eddy * I_eddy)';
    section.PSI_eddy   = 'Breakdown and burn-through';
    mode.PSI_eddy      = 'advanced';

    valeur.I_eddy    = 1;
    type.I_eddy      = 'real';
    borne.I_eddy     = {-1,0,1};
    defaut.I_eddy    = 1;
    info.I_eddy      = 'Initial eddy current:\nif = 0, set inital eddy current to 0;\nif = 1, take the maximum between breakdown voltage divided by R_eddy, and initial plasma current waveform as  inital eddy current;\nif = -1, take minus the maximum between breakdown voltage divided by R_eddy, and initial plasma current waveform as  inital eddy current';
    section.I_eddy   = 'Breakdown and burn-through';
    mode.I_eddy      = 'advanced';

    valeur.p_prefill   = 1e-3;
    type.p_prefill     = 'real';
    borne.p_prefill    = [0,1000];
    defaut.p_prefill   = 1e-3;
    info.p_prefill     = 'prefill pressure (Pa, 1 Torr = 133.3 Pa)';
    section.p_prefill  = 'Breakdown and burn-through';
    mode.p_prefill     = 'advanced';
   
    valeur.temp_vac   = 300;
    type.temp_vac     = 'real';
    borne.temp_vac    = [1.8,3000];
    defaut.temp_vac   = 300;
    info.temp_vac     = 'temperature of vacuum vessel and cold neutrals (K)';
    section.temp_vac  = 'Breakdown and burn-through';
    mode.temp_vac     = 'advanced';

    valeur.VV_volume   = 0;
    type.VV_volume     = 'real';
    borne.VV_volume    = [0.1,10000];
    defaut.VV_volume   = 0;
    info.VV_volume     = 'vacuum vessel volume; if = 0, the maximum plasma volume defined by LCFS is used (m^3)';
    section.VV_volume  = 'Breakdown and burn-through';
    mode.VV_volume     = 'advanced';
  
    valeur.initiation_only   = 0;
    type.initiation_only     = 'integer';
    borne.initiation_only    = {0,1};
    defaut.initiation_only   = 0;
    info.initiation_only     = 'if = 1, stay in initiation mode for the whole simulation';
    section.initiation_only  = 'Breakdown and burn-through';
    mode.initiation_only     = 'advanced';
  
    valeur.zeff    = 0;
    type.zeff      = 'integer';
    borne.zeff     = {0,1,2,3,4,5,6,7,8,9};
    defaut.zeff    = 0;
    info.zeff      = sprintf('Zeff: 0-> reference & flat,\n1 -> average given + profile effect,\n 2,3 & 4  -> scaling + profile effect for wall/divertor in C,Be or W;\n 5 -> Tore Supra scaling,\n 6-> Matthews scaling if not used for radiative power,\n 7 -> universal scaling law (J. G. Cordey rapport JET-P(85) 28),\n 8 -> universal scaling law +  profile effect,\n 9 -> I. Erofeev scaling expressed in Greenwald fraction (for rampup)');
    section.zeff   = 'Composition';
    mode.zeff      = 'advanced';
    
    valeur.faccu    = 0.5;
    type.faccu      = 'real';
    borne.faccu     = [-1,10];
    defaut.faccu    = 0.5;
    info.faccu      = 'factor of accumulation for heavy impurity in the core plasma: 1) with W_effect = 0: \n works with zeff key values that allow profile effect (1,2,3,4 & 8);\n n_heavy = faccu * n_accumulation(x) + (1 - faccu) * r .* n_e(x);\n if = 0, no accumulation;\n if > 0, use d neoclassical simplified formula depending on density peaking and temperature peaking;\n if < 0, use d neoclassical simplified formula depending only on density peaking;\n 2) with W_effect = 1:\na) if acc_col is off, factor in exponential exp(acc_inte * faccu);\nb) if acc_col is on, amplitude of turbulent transport added to neoclassical part with D_W = faccu * Chi_ion (thermal);\nif Sn_fraction > 0, a tin fraction will be add to the tungsten..';
    section.faccu   = 'Composition';
    mode.faccu      = 'advanced';

    valeur.acc_col     = 'off';
    type.acc_col      = 'string';
    borne.acc_col     = {'off','on'};
    defaut.acc_col    = 'off';
    info.acc_col      = 'turn on or off collisionality dependance of factor in neoclassical formulation for tungsten accumulation or tunsten and tin if Sn_fraction > 0.';
    section.acc_col   = 'Composition';
    mode.acc_col      = 'advanced';

    valeur.heat_acc    = 0;
    type.heat_acc      = 'real';
    borne.heat_acc     = [-10,10];
    defaut.heat_acc    = 0;
    info.heat_acc      = 'factor for the plasma heating decontamination (works with W_effect = 1 only): if > 0, electron heating tends to reduce impurity accumulation and ion heating tends to increase impurity accumulation in core plasma; reverse sign adds inverse effect.';
    section.heat_acc   = 'Composition';
    mode.heat_acc      = 'advanced';

    valeur.rot_acc    = 0;
    type.rot_acc      = 'real';
    borne.rot_acc     = [-10,10];
    defaut.rot_acc    = 0;
    info.rot_acc      = 'factor for the plasma rotation decontamination (works with W_effect = 1 only): if > 0, rotation tends to reduce impurity accumulation;\notherwise if < 0,rotation tends to increase impurity accumulation;\n if = 0, no effect of rotation on impurity accumulation.';
    section.rot_acc   = 'Composition';
    mode.rot_acc      = 'advanced';

    valeur.fne_acc   = 0;
    type.fne_acc     = 'real';
    borne.fne_acc    = [-2,10];
    defaut.fne_acc   = 0;
    info.fne_acc     = 'with W_effect = 1: exponent applied to the normalized electron density for the shape of tungsten density;\nn_W(x) = (C(SOL,divertor)+ cw_offset * ne_edge)  * (n_e(x)/n_edge)^fne_acc *  exp(acc_inte * faccu);\nif Sn_fraction > 0, a tin fraction will be add to the tungsten.';
    section.fne_acc   = 'Composition';
    mode.fne_acc      = 'advanced';

    valeur.gradient_Wacc    = 0;
    type.gradient_Wacc      = 'integer';
    borne.gradient_Wacc     = [0,10];
    defaut.gradient_Wacc    = 0;
    info.gradient_Wacc      = 'factor applied to gradients computation (2 points formula) in pedestal; used in  W accumulation formula (value must be the same as ffit_ped)\nif = 0 use standard 3 points derivative';
    section.gradient_Wacc   = 'Composition';
    mode.gradient_Wacc      = 'advanced';

    valeur.W_effect    = 0;
    type.W_effect      = 'integer';
    borne.W_effect     = {0,1};
    defaut.W_effect    = 0;
    info.W_effect      = 'if = 1, take into account specific effect of Tungsten (as variation of the ionisation state in plasma).\nif used, SOL parameters must be adjusted (at least cw_factor and cw_offset);\nif Sn_fraction > 0, a tin fraction will be add to the tungsten.';
    section.W_effect   = 'Composition';

    valeur.Sn_fraction    = 0;
    type.Sn_fraction      = 'real';
    borne.Sn_fraction     = [0,1];
    defaut.Sn_fraction    = 0;
    info.Sn_fraction      = 'The purpose of this parameter is to allow study of liquid Sn divertor: if > 0, fraction of Sn in heavy impurities compostion: n_Sn = Sn_fraction * n_heavy; n_W = (1 -Sn_fractio) * n_heavy and n_heavy is identified to zerod.nwm and profile.nwp';
    section.Sn_fraction   = 'Composition';
    mode.Sn_fraction      = 'advanced';

    valeur.zmax    = 8;
    type.zmax      = 'integer';
    borne.zmax     = [3,100];
    defaut.zmax    = 8;
    info.zmax      = 'charge number of main impurity responsible for radiated power (example: O, Ar, Ne, Xe);\nif 0 chosen, the default value of 3 will be imposed;\nZmax is used to constrained the max Zeff value possible (plasma of only Zmax ions will have Zeff = Zmax)';
    section.zmax   = 'Composition';

    valeur.zimp    = 6;
    type.zimp      = 'integer';
    borne.zimp     = [3,100];
    defaut.zimp    = 6;
    info.zimp      = 'charge number of light impurity (example: C, Be);\nif 0 chosen, the default value of 3 will be imposed';
    section.zimp   = 'Composition';

    valeur.rimp     = 0.1;
    type.rimp      = 'real';
    borne.rimp     = [0,1];
    defaut.rimp    = 0.1;
    info.rimp      = 'ratio between density of main impurity (n_zmax) and density of light impurity (n_zimp);\nhas to be between 0 and 1';
    section.rimp   = 'Composition';
    
    valeur.natural_nD_o_nH    = 0.000115;
    type.natural_nD_o_nH      = 'real';
    borne.natural_nD_o_nH     = [0,1];
    defaut.natural_nD_o_nH    = 0.000115;
    info.natural_nD_o_nH      = 'ratio between density deuterium and hydrogen for proton boron fusion (option.gaz ==11); remark: the default value in the natural concentration on Earth.';
    section.natural_nD_o_nH   = 'Composition';
    mode.natural_nD_o_nH      = 'advanced';
 
    valeur.density_model    = 'minconv';
    type.density_model      = 'string';
    borne.density_model     = {'default','control','curvature','minconv'};
    defaut.density_model    = 'minconv';
    info.density_model      = 'model used to identify density transport coefficients (post processing; have no impact on density profile)';
    section.density_model   = 'Density';
    mode.density_model      = 'advanced';
    
    valeur.ane_factor    = 1 ;
    type.ane_factor      = 'real';
    borne.ane_factor     = [0.1,10];
    defaut.ane_factor    = 1;
    info.ane_factor      = 'multiplication factor applied to density peaking prediction for ane >= 10 (used to modulate scaling law prediction)';
    section.ane_factor   = 'Density';
    mode.ane_factor      = 'advanced';

    valeur.frad    = 1 ;
    type.frad      = 'real';
    borne.frad     = [0.1,10];
    defaut.frad    = 1;
    info.frad      = 'multiplication factor for line radiative power';
    section.frad   = 'Radiation';

    valeur.fte_edge    = 1 ;
    type.fte_edge      = 'real';
    borne.fte_edge     = [0.1,10];
    defaut.fte_edge    = 1;
    info.fte_edge      = 'multiplication factor applied to LCFS temperature (for studies of egde radiation and poloidal flux comsumption)';
    section.fte_edge   = 'Radiation';
    mode.fte_edge      = 'advanced';

    valeur.te_edge_fixed    = 0 ;
    type.te_edge_fixed      = 'real';
    borne.te_edge_fixed     = [0,1000];
    defaut.te_edge_fixed    = 0;
    info.te_edge_fixed      = 'if > 0, force the value of LCFS electron temperature to te_edge_fixed (eV) whatever is other settings (if Te profile is provided as a external data, LCFS electron temperature is read from the profile. A upper limiting rule is also applied to keep max(Te(x)) > Te_a )';
    section.te_edge_fixed   = 'Radiation';
    mode.te_edge_fixed      = 'advanced';

   valeur.matthews    = 1;
    type.matthews      = 'integer';
    borne.matthews     = {-1,0,1,2};
    defaut.matthews    = 1;
    label.matthews     = 'line radiation model';
    info.matthews      = 'line radiation model:\nif = 1, line radiative power (Pline) is normalised to Matthews scaling law;\nif = 2, Pline is normalised to Rapp scaling law;\nif = 0, Pline is given by cooling rate formulation (coronal equilibrium) for the core plasma. Difference between Matthews scaling and cooling rate formulation is used for SOL + Divertor radiative power;\nif = -1, same as 0 but use improved computation for radiative mantle';
    section.matthews  = 'Radiation';

    valeur.z_prad    = 'zmax';
    type.z_prad      = 'string';
    borne.z_prad     = {'zmax','zimp','Stangeby'};
    defaut.z_prad    = 'zmax';
    info.z_prad      = 'Charge number used in scaling formula: zmax (original version), zimp or Z_Stangeby = formula defined in Stangeby book';
    section.z_prad   = 'Radiation';
    mode.z_prad      = 'advanced';

    valeur.gaunt    = 0;
    type.gaunt      = 'integer';
    borne.gaunt     = {0,1};
    defaut.gaunt    = 0;
    info.gaunt      = 'Gaunt factor:\nif = 0, use fixed Gaunt factor (1.1027);\nif = 1, compute Gaunt factor depending on Te and main ion Z for bremsstrahlung';
    section.gaunt   = 'Radiation';
    mode.gaunt      = 'advanced';

    valeur.cor_rel_brem    = 'standart';
    type.cor_rel_brem      = 'string';
    borne.cor_rel_brem     = {'standart','improved'};
    defaut.cor_rel_brem    = 'standart';
    info.cor_rel_brem      = 'if = improved, use improved formula for Bremsstrahlung relativistic correction';
    section.cor_rel_brem   = 'Radiation';
    mode.cor_rel_brem      = 'advanced';

    valeur.noncoronal    = 0;
    type.noncoronal      = 'integer';
    borne.noncoronal     = {-1,0,1,2,3,4,5};
    defaut.noncoronal    = 0;
    info.noncoronal      = 'non coronal radiative power:\nif = -1, use non coronal equilibrium (only available for He, Li, Be, C, O, N, Ne, Ar);\nif = 0, use coronal equilibrium;\nif = 1, use a simplified model for non coronal radiative power computation;\nif > 1 , use a simplified model for non coronal radiative power computation in divertor only;\nif = 2,  full effect;\n if = 3, half effect;\n if = 4, 0.2 * full effect;\nif = 5, 0.1 * full effect';
    section.noncoronal   = 'Radiation';
    mode.noncoronal      = 'advanced';

    valeur.rw    = 0.7;
    type.rw      = 'real';                    
    borne.rw     = [-1,1];
    defaut.rw    = 0.7;
    info.rw      = 'wall reflection coefficient of the cyclotron radiation; if < 0, use LATF model instead of Albajar scaling';
    section.rw   = 'Radiation';

    valeur.configuration    = 3;
    type.configuration      = 'integer';
    borne.configuration     = {0,1,2,3,4};  
    defaut.configuration    = 3;                
    info.configuration      = 'PFC configuration:\nif = 0,poloidal limiter in L and H mode;\nif = 1, toroidal limiter in L and H mode;\nif = 2, poloidal limiter or divertor depending on LCFS shape (X point automatic detection);\nif = 3, toroidal limiter  or divertor depending on LCFS shape (X point automatic detection);\nif = 4, divertor in both L and H mode';	
    section.configuration   = 'SOL';

    valeur.lambda_scale    = 3;
    type.lambda_scale      = 'real';                    
    borne.lambda_scale     = {0,1,2,3,4};
    defaut.lambda_scale    = 3;
    info.lambda_scale      = 'scaling used for SOL width:\nif = 0, uses fraction of a or R (see sol_lscale);\nif = 1, uses Goldston model (gives similar result to Eich scaling);\nif = 2,  uses Goldston model in H-mode diverted plasma, otherwise uses given width;\nif = 3, uses Goldston model in H-mode diverted plasma, otherwise uses Halpern scaling law\nif = 4, uses D. Brunner et al 2018 Nucl. Fusion 58 094002, equation 4 on volume averaged plasma pressure';
    section.lambda_scale   = 'SOL';

    valeur.factor_scale    = 1;
    type.factor_scale      = 'real';                    
    borne.factor_scale     = [0.1,100];
    defaut.factor_scale    = 1;
    info.factor_scale      = 'multiplication factor applied to the SOL width when it is defined by a scaling law, for H mode only (Dsol_Hmode = factor_scale * Goldston scaling)';
    section.factor_scale   = 'SOL';
    mode.factor_scale      = 'advanced';

    valeur.sol_lscale    = 0;
    type.sol_lscale      = 'real';                    
    borne.sol_lscale     = [-0.1,0.1];
    defaut.sol_lscale    = 0;
    info.sol_lscale      = 'scaling parameter of folding length in the SOL:\nif = 0, a / 100 (previous metis version);\nif > 0, scale as a * sol_lscale (typical value 0.01);\nif < 0, scale as R0 * sol_lscale (typical value 0.003)';
    section.sol_lscale   = 'SOL';

    valeur.eioniz    = 25;
    type.eioniz      = 'real';                    
    borne.eioniz     = [0,1000];
    defaut.eioniz    = 25;
    info.eioniz      = 'average ionisation energy per atom;\nif = 0, tabulated value dependent on Te and ne is used';
    section.eioniz   = 'SOL';
    mode.eioniz      = 'advanced';
    
    valeur.cx_ion    = 0;
    type.cx_ion      = 'real';                    
    borne.cx_ion     = {0,1};
    defaut.cx_ion    = 0;
    info.cx_ion      = 'if = 1, take into account power loss by ions due to charge exchange with neutrals';
    section.cx_ion   = 'SOL';
    mode.cx_ion      = 'advanced';
    

    valeur.de	        = 0.5;   
    type.de		= 'real';     
    borne.de	        = [0.1,10];      
    defaut.de	        = 0.5; 
    info.de		= 'two points model: secondary electron emission coefficient';
    section.de		= 'SOL';
    mode.de             = 'advanced';

    valeur.alpha_e	= 0.82;   
    type.alpha_e	= 'real';     
    borne.alpha_e	= [0.15,3];      
    defaut.alpha_e	= 0.82; 
    info.alpha_e        = 'two points model: multiplication factor applied to parallel heat flux in formula for the kinetic effect correction to fluid formulation';
    section.alpha_e     = 'SOL';
    mode.alpha_e        = 'advanced';

    valeur.fnesol    = 0;
    type.fnesol      = 'real';                    
    borne.fnesol     = [0,1];
    defaut.fnesol    = 0;
    info.fnesol      = 'interpolation factor to compute average density in SOL along magnetic field line:\n(used for SOL radiative power computation)\nne_sol = (1 -fnesol) * ne_LCFS + fnesol * ne_target';
    section.fnesol   = 'SOL';
    mode.fnesol      = 'advanced';

    valeur.sol_model    = 'scaling';
    type.sol_model      = 'string';
    borne.sol_model     = {'scaling','2_points'};
    defaut.sol_model    = 'scaling';
    info.sol_model      = 'model used to compute egde (LCFS) values';
    section.sol_model   = 'SOL';

    valeur.sol_rad    = 'coupled';
    type.sol_rad      = 'string';                    
    borne.sol_rad     = {'coupled','decoupled'};
    defaut.sol_rad    = 'coupled';
    info.sol_rad      = 'SOL radiative power: if = coupled, the radiative power in SOL and core plasma are coupled (Cooling rate is used to compute core plasma radiative power and difference between scaling law and core radiative power gives SOL radiative power);\nif = decoupled, the radiative powers in SOL and divertor are separately computed (must be used with 2-points model)';
    section.sol_rad   = 'SOL';

    valeur.lcx    = 1;
    type.lcx      = 'real';                    
    borne.lcx     = [1,20];
    defaut.lcx    = 1;
    info.lcx      = 'connexion length multiplication factor to take into account LCFS safety factor divergence due to X point (connexion length ~ lcx * R * q_eff for diverted plasma)';
    section.lcx   = 'SOL';
    mode.lcx      = 'advanced';

    valeur.fpower    = 0.6;
    type.fpower      = 'real';                    
    borne.fpower     = [0.1,1];
    defaut.fpower    = 0.6;
    info.fpower      = 'two points model: fraction of power crossing the LCFS that goes to the target (only one target is considered in two points model, generally the outer one)';
    section.fpower   = 'SOL';
    mode.fpower      = 'advanced';

    valeur.fR_target    = 1;
    type.fR_target      = 'real';                    
    borne.fR_target     = [-3,3];
    defaut.fR_target    = 1;
    info.fR_target      = 'two points model: multiplication factor that gives the target radius as R_target = abs(fR_target) * R0;\nif < 0, the radial position of the target is taken into account in two point model (model from T.W. Petrie et al, Nuc. Fus. 53 (2013) 113024)';
    section.fR_target   = 'SOL';
    mode.fR_target      = 'advanced';

    valeur.Recycling_target    = 0;
    type.Recycling_target      = 'real';
    borne.Recycling_target     = [0,1];
    defaut.Recycling_target    = 0;
    info.Recycling_target      = 'Recycling coefficient at divertor target; if = 0, the recycling coefficient at divertor target is the global recycling coefficient';
    section.Recycling_target   = 'SOL';

    valeur.detach    = 0;
    type.detach      = 'real';
    borne.detach     = [0,10];
    defaut.detach    = 0;
    info.detach      = 'if > 0 switch on additional term for detachement from Apiwat/Pegourie model. In this case the value of parameter detach gives the exponent of the model (original model use the value of 1)';
    section.detach   = 'SOL';

    valeur.fcond    = 1;
    type.fcond      = 'real';                    
    borne.fcond     = [-1,1];
    defaut.fcond    = 1;
    info.fcond      = 'two points model: fraction of conductive // power versus convective+conductive;\nif > 0, without kinetic correction;\nif < 0, with kinetic correction';
    section.fcond   = 'SOL';
    mode.fcond      = 'advanced';

    valeur.fmom    = 1;
    type.fmom      = 'real';                    
    borne.fmom     = [0,1];
    defaut.fmom    = 1;
    info.fmom      = 'two points model: fmom * (1 + Tiu/Teu);\nfmom = factor to take into account momentum loss by ions (friction on neutrals, ...);\nTiu = LCFS Ti; Teu = LCFS Te;\nif fmom = 0, uses the simple model based on the balance between charge exchange and ionisation rate at the target (used for detachment studies)';
    section.fmom   = 'SOL';
    mode.fmom      = 'advanced';
    
    valeur.Sq    = 0;
    type.Sq      = 'real';                    
    borne.Sq     = [-1,100];
    defaut.Sq    = 0;
    info.Sq      = 'two points model: if >= 0, parallel heat flux spreading due to turbulence (mm);\nif = -1, use scalings from A. Scarabosio paper (Journal of Nuclear Materials 463 (2015) 49-54)';
    section.Sq   = 'SOL';
    mode.Sq      = 'advanced';

    valeur.mach_corr    = 0;
    type.mach_corr     = 'integer';                    
    borne.mach_corr     = {0,1};
    defaut.mach_corr    = 0;
    info.mach_corr     = 'two points model: supersonic flow;\nif = 0, assume Mach number = 1 at the target;\nif = 1, take into account possible Mach number > 1 near the target';
    section.mach_corr  = 'SOL';
    mode.mach_corr      = 'advanced';

    valeur.yield_model      = 'Javev';
    type.yield_model        = 'string';
    borne.yield_model       = {'fit','Javev','Matsunami'};
    defaut.yield_model      = 'Javev';
    info.yield_model        = 'Tungsten density - model used to compute sputtering yields:\nif = fit, use fit of DIVIMP simulation;\nif = Javev, use model from Janev paper;\nif = Matsunami, use model from Matsunami report';
    section.yield_model    = 'SOL';
    mode.yield_model      = 'advanced';

    valeur.ftwleak      = -0.5;
    type.ftwleak        = 'real';                    
    borne.ftwleak       = [-1,1];
    defaut.ftwleak      = 0.5;
    info.ftwleak        = 'Tungsten density - interpolation factor used to compute W leakage from target to LCFS: Tleak=(abs(ftwleak)*Zw(Te_LCFS)+(1-abs(ftwleak))*Zw(Te_plate))*sqrt(ne*S);\nif > 0, fw_leak = exp(-(Tleak /Tplate)^2);\n if < 0, fw_leak = fit_DIVIMP(Tleak /Tplate).;\n if Sn_fraction > 0, applied to W + Sn using the averaged charge weighted by Sn_fraction';
    section.ftwleak     = 'SOL';
    mode.ftwleak      = 'advanced';

    valeur.cw_factor    = 1;
    type.cw_factor      = 'real';                    
    borne.cw_factor     = [-10,10];
    defaut.cw_factor    = 1;
    info.cw_factor      = 'Tungsten density: abs(cw_factor) = multiplication factor applied to the source of tungsten,\nused to compute tungsten density profile in core plasma;\nif > 0, no redeposition of tungsten on target;\nif < 0 uses a simple prompt redeposition model on target;\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_factor   = 'SOL';

    valeur.cw_offset    = 0;
    type.cw_offset      = 'real';                    
    borne.cw_offset     = [0,0.01];
    defaut.cw_offset    = 0;
    info.cw_offset      = 'Tungsten density: tungsten concentration due to sources other than divertor targets and auxiliary heating system,\n expressed in concentration relative to electron density (<nw_residual> / <n_e>);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_offset   = 'SOL';

    valeur.cw_ecrh    = 0;
    type.cw_ecrh      = 'real';                    
    borne.cw_ecrh     = [-0.01,0.01];
    defaut.cw_ecrh    = 0;
    info.cw_ecrh      = 'Tungsten density: tungsten concentration due to sources linked with ECRH, expressed in concentration relative to electron density (<nw_ecrh> / <n_e> per MW of ECRH);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_ecrh   = 'SOL';
    mode.cw_ecrh      = 'advanced';

    valeur.cw_icrh    = 0;
    type.cw_icrh      = 'real';                    
    borne.cw_icrh     = [-0.01,0.01];
    defaut.cw_icrh    = 0;
    info.cw_icrh      = 'Tungsten density: tungsten concentration due to sources linked with ICRH, expressed in concentration relative to electron density (<nw_ecrh> / <n_e> per MW of ICRH);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_icrh   = 'SOL';
    mode.cw_icrh      = 'advanced';

    valeur.cw_lhcd    = 0;
    type.cw_lhcd      = 'real';                    
    borne.cw_lhcd     = [-0.01,0.01];
    defaut.cw_lhcd    = 0;
    info.cw_lhcd      = 'Tungsten density: tungsten concentration due to sources linked with LHCD, expressed in concentration relative to electron density (<nw_ecrh> / <n_e> per MW of LHCD);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_lhcd   = 'SOL';
    mode.cw_lhcd      = 'advanced';

    valeur.cw_nbi1    = 0;
    type.cw_nbi1      = 'real';                    
    borne.cw_nbi1     = [-0.01,0.01];
    defaut.cw_nbi1    = 0;
    info.cw_nbi1      = 'Tungsten density: tungsten concentration due to sources linked with first NBI, expressed in concentration relative to electron density (<nw_ecrh> / <n_e> per MW of NBI);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_nbi1   = 'SOL';
    mode.cw_nbi1      = 'advanced';

    valeur.cw_nbi2    = 0;
    type.cw_nbi2      = 'real';                    
    borne.cw_nbi2     = [-0.01,0.01];
    defaut.cw_nbi2    = 0;
    info.cw_nbi2      = 'Tungsten density: tungsten concentration due to sources linked with second NBI, expressed in concentration relative to electron density (<nw_ecrh> / <n_e> per MW of NBI);\n if Sn_fraction > 0, applied to W + Sn.';
    section.cw_nbi2   = 'SOL';
    mode.cw_nbi2      = 'advanced';

    valeur.imp_div    = 'fixed';
    type.imp_div     = 'string';                    
    borne.imp_div     = {'fixed','auto'};
    defaut.imp_div    = 'fixed';
    info.imp_div      = 'two points model:if = ''auto'', compute radiative impurity (zmax) concentration from core concentration and leakage from divertor;\n otherwise use parameter fzmax_div to fixe divertor concentration enrichment in radiative impurity';
    section.imp_div   = 'SOL';
    mode.imp_div      = 'advanced';

    valeur.fzmax_div    = 0;
    type.fzmax_div     = 'real';                    
    borne.fzmax_div     = [-100,100];
    defaut.fzmax_div    = 0;
    info.fzmax_div      = 'two points model: fraction (in percent) of radiative impurity, associated to charge zmax, in divertor;\ncontrols the radiative fraction in divertor;\nif > 0, the minimal temperature in the radiative power integral along field line is 0, as in original model(R. Clark et al, Journal of Nuclear Materials 1995);\n if < 0, the fraction is abs(fzmax_div) and the minimal temperature in integral is te_lim';
    section.fzmax_div   = 'SOL';
    mode.fzmax_div      = 'advanced';

    valeur.carbonblow    = 0;
    type.carbonblow     = 'real';                    
    borne.carbonblow     = [-3,3];
    defaut.carbonblow    = 0;
    info.carbonblow      = 'two points model - carbon density enhancement in divertor due to sputtering:\nif> 0, the density of carbon in divertor is increased proportionally to physical sputtering * carbonblow;\nif <0, the density of carbon in divertor is increased proportionally to total sputtering (physical + chemical + ...) * abs(carbonblow)';
    section.carbonblow   = 'SOL';
    mode.carbonblow      = 'advanced';

    valeur.residence_time    = 0;
    type.residence_time     = 'real';                    
    borne.residence_time     = [0,1];
    defaut.residence_time    = 0;
    info.residence_time      = 'Impurity residence time in SOL and divertor (s). This parameter is used to compute radiative power in divertor with two points model out of ionisation equilibrium;\n if = 0, the impurity residence time is then set to the connexion length divided by sound speed averaged on mahgnetic line';
    section.residence_time   = 'SOL';
    mode.residence_time      = 'advanced';

    valeur.angle_ece    = 90;
    type.angle_ece      = 'real';                    
    borne.angle_ece     = {0,90,180};  
    defaut.angle_ece    = 90;                
    info.angle_ece      = 'ECCD efficiency - trapped particles effect: ECRH effectif deposition poloidal location: 0 -> LFS, 90 -> top, 180 -> HFS.';
    section.angle_ece   = 'ECRH/ECCD';
    label.angle_ece     = 'poloidal location';

    valeur.synergie    = 1;
    type.synergie      = 'real';                    
    borne.synergie     = [0,10];  
    defaut.synergie    = 1;                
    info.synergie      = 'ECCD efficiency - synergy effect:ECCD effciency multiplication factor due to synergy with LHCD;\nif = 1, no  synergy;\nif = 0, computed from overlap of current sources due to ECCD and LHCD';
    section.synergie   = 'ECRH/ECCD';
    mode.synergie      = 'advanced';

    valeur.sens    = 0;
    type.sens      = 'integer';                    
    borne.sens     = {-1,0,1};  
    defaut.sens    = 0;   
    label.sens     = 'CD direction';
    info.sens      = 'ECCD current drive orientation:\nif = -1,  counter-current;\nif = 0, normal injection (no current);\nif = 1, co-current (same as Ip)';
    section.sens   = 'ECRH/ECCD';

    valeur.eccdmul    = 1;
    type.eccdmul     = 'real';                    
    borne.eccdmul     = [0.1,10];  
    defaut.eccdmul    = 1;                
    info.eccdmul      = 'ECCD efficiency: multiplication factor applied to ECRH current drive efficiency';	
    section.eccdmul   = 'ECRH/ECCD';

    valeur.width_ecrh    = 0;
    type.width_ecrh      = 'real';                    
    borne.width_ecrh     = [-10,1];  
    defaut.width_ecrh    = 0;                
    info.width_ecrh      = 'ECCD deposition width: if = 0, use internal metis formula to compute ECRH/ECCD width;\nif > 0, width_ecrh is the width of Gaussian deposition;\nif < 0, the internal metis formula value is mutiplied by width_ecrh to obtain the width of ECRH/ECCD deposition';	
    section.width_ecrh   = 'ECRH/ECCD';

    valeur.angle_nbi    = 90;
    type.angle_nbi      = 'real';
    borne.angle_nbi     = [-90,90];
    defaut.angle_nbi    = 90;
    info.angle_nbi = 'first NB injector: NBI beam angle  (<0 -> counter-current , 0 = normal, >0 -> co-current);\n used to compute the fraction of power injected perpendiculary (cos(angle_nbi)) and not perpendiculary (sin(angle_nbi))';
    section.angle_nbi   = 'NBI/NBICD';

    valeur.rtang    = 0;
    type.rtang      = 'real';
    borne.rtang     = [0,100];
    defaut.rtang    = 0;
    info.rtang      = 'first NB injector: tangency radius of neutral beam;\nif = 0, use angle_nbi to compute Rtang (m)';
    section.rtang   = 'NBI/NBICD';

    valeur.zext    = 0;
    type.zext      = 'real';
    borne.zext     = [0,0.5];
    defaut.zext    = 0;
    info.zext      = 'first NB injector: vertical shift at the center of the plasma of the neutral beam trajectory (normalized radius)';
    section.zext   = 'NBI/NBICD';
    
    valeur.drs1    = 0;
    type.drs1      = 'real';
    borne.drs1     = [0,0.25];
    defaut.drs1    = 0;
    info.drs1      = 'first NB injector: half width of neutral beam in toroidal direction expressed in units of minor radius (computed as variation of tangential radius);\nif = 0, use default METIS value (1/6)';
    section.drs1   = 'NBI/NBICD';
    mode.drs1      = 'advanced';
    
    valeur.dzs1    = 0;
    type.dzs1      = 'real';
    borne.dzs1     = [0,0.25];
    defaut.dzs1    = 0;
    info.dzs1      = 'first NB injector: normalised half width of neutral beam in vertical direction;\nif = 0, use default METIS value (0.05)';
    section.dzs1   = 'NBI/NBICD';
    mode.dzs1      = 'advanced';
    
    valeur.einj    = 1e6;
    type.einj      = 'real';                    
    borne.einj     = [1e4,1e7];  
    defaut.einj    = 1e6;                
    info.einj      = 'first NB injector: NBI beam energy (eV)';
    section.einj   = 'NBI/NBICD';

    valeur.nbicdmul    = 1;
    type.nbicdmul      = 'real';                    
    borne.nbicdmul     = [0.1,10];  
    defaut.nbicdmul    = 1;                
    info.nbicdmul      = 'first NB injector: multiplication factor applied to NBI current drive efficiency';	
    section.nbicdmul   = 'NBI/NBICD';

    valeur.e_shielding    = 'Honda-NEO';
    type.e_shielding      = 'string';                    
    borne.e_shielding     = {'Lin-Liu','Honda-Sauter','Honda-NEO'};  
    defaut.e_shielding    = 'Honda-NEO';                
    info.e_shielding      = 'current shielding factor model :\nif = Lin-Liu, use Y.R Lin-Liu & F. L. Hilton  model (Physics of Plasmas 4 (11) 1997);\nif = Honda-Sauter, use model with collisionality depedence (M. Honda et al NF 52 (2012) p 023021);\nif = Honda-NEO, same as Honda-Sauter but with L31 fitted on NOE code (A. Redl et al, Phys. Plasmas 28, 022502 (2021); https://doi.org/10.1063/5.0012664)';	
    section.e_shielding   = 'NBI/NBICD';
    mode.e_shielding      = 'advanced';

    valeur.fast_ion_sbp    = 3;
    type.fast_ion_sbp     = 'integer';                    
    borne.fast_ion_sbp     = {0,1,2,3,4};  
    defaut.fast_ion_sbp    = 3;                
    info.fast_ion_sbp      = 'Take or not into account increment of the stopping cross section due to fats ions (K. Okano et al, ECA vol 25A (2001) p 809):\nif = 0, no increment & Janev 1989 cross section (without impurities effect);\nif = 1, take into account increment & Janev 1989 cross section (with impurities effect);\n if = 2, no increment & Janev 1989 cross section (with impurities effect);\nif = 3, take into account increment & Suzuki 1998 cross section (with impurities effect);\n if = 4, no increment & Suziki 1998 cross section (with impurities effect)';	
    section.fast_ion_sbp   = 'NBI/NBICD';
    mode.fast_ion_sbp      = 'advanced';

    valeur.shinethrough    =  0;
    type.shinethrough     = 'integer';                    
    borne.shinethrough     = {0,1,2};  
    defaut.shinethrough    = 0;                
    info.shinethrough      = 'For testing; allows to turn off first orbit losses and/or shinethrough:\nif = 0, both first orbit losses and shinethrough are taking into account;\nif = 1, only shinethrough is taking into account (first orbit losses are discarded);\nif = 2, no losses are taking into account (both first orbit losses and shinethrough are discarded)';	
    section.shinethrough   = 'NBI/NBICD';
    mode.shinethrough      = 'advanced';

    valeur.cur_nbi_time    = 0;
    type.cur_nbi_time      = 'real';                    
    borne.cur_nbi_time     = [0,1];  
    defaut.cur_nbi_time    = 0;                
    info.cur_nbi_time      = 'Choose temporal evolution of I_NBICD:\nif = 0, same behaviour as PNBI waveform;\nif = 1; same behaviour as thermal PNBI_th;\nintermediate value gives behaviour proportional to cur_nbi_time * PNBI_th + (1-cur_nbi_time) * PNBI (reference);\nPossibly inacurate in evolution mode (inside Simulink/Kepler): must be set to 0 in this case';
    section.cur_nbi_time   = 'NBI/NBICD';
    mode.cur_nbi_time      = 'advanced';
    
    valeur.forced_H_NBI    = 0;
    type.forced_H_NBI      = 'integer';                    
    borne.forced_H_NBI     = {0,1};  
    defaut.forced_H_NBI    = 0;                
    info.forced_H_NBI      = 'if == 1, force the composition of neutral beams to be hydrogen what ever is set in option.gaz or reference ftnbi';
    section.forced_H_NBI   = 'NBI/NBICD';
    mode.forced_H_NBI      = 'advanced';

    valeur.nb_nbi    = 1;
    type.nb_nbi      = 'integer';                    
    borne.nb_nbi     = {1,2};  
    defaut.nb_nbi    = 1;                
    info.nb_nbi      = 'number of NBI injectors used in METIS';
    section.nb_nbi  = 'NBI/NBICD@2';

    valeur.angle_nbi2    = 90;
    type.angle_nbi2      = 'real';
    borne.angle_nbi2     = [-90,90];
    defaut.angle_nbi2    = 90;
    info.angle_nbi2      = 'second NB injector: NBI beam angle  (<0 -> counter-current , 0 = normal, >0 -> co-current);\n used to compute the fraction of power injected perpendiculary (cos(angle_nbi)) and not perpendiculary (sin(angle_nbi))';
    section.angle_nbi2   = 'NBI/NBICD@2';

    valeur.rtang2    = 0;
    type.rtang2      = 'real';
    borne.rtang2     = [0,100];
    defaut.rtang2    = 0;
    info.rtang2      = 'second NB injector: tangency radius of neutral beam;\nif = 0, use angle_nbi to compute Rtang (m)';
    section.rtang2   = 'NBI/NBICD@2';

    valeur.zext2    = 0;
    type.zext2      = 'real';
    borne.zext2     = [0,0.5];
    defaut.zext2    = 0;
    info.zext2      = 'second NB injector: vertical shift at the center of the plasma of the neutral beam trajectory (normalized radius)';
    section.zext2   = 'NBI/NBICD@2';
    
    valeur.drs2    = 0;
    type.drs2      = 'real';
    borne.drs2     = [0,0.25];
    defaut.drs2    = 0;
    info.drs2      = 'second NB injector: half width of neutral beam in toroidal direction expressed in units of minor radius (computed as variation of tangential radius);\nif = 0, use default METIS value (1/6)';
    section.drs2   = 'NBI/NBICD@2';
    mode.drs2      = 'advanced';
   
    valeur.dzs2    = 0;
    type.dzs2      = 'real';
    borne.dzs2     = [0,0.25];
    defaut.dzs2    = 0;
    info.dzs2      = 'second NB injector: normalised half width of neutral beam in vertical direction;\nif = 0, use default METIS value (0.05)';
    section.dzs2   = 'NBI/NBICD@2';
    mode.dzs2      = 'advanced';
    
    valeur.einj2    = 1e6;
    type.einj2      = 'real';                    
    borne.einj2     = [1e4,1e7];  
    defaut.einj2    = 1e6;                
    info.einj2      = 'second NB injector: NBI beam energy (eV)';
    section.einj2   = 'NBI/NBICD@2';

    valeur.nbicdmul2    = 1;
    type.nbicdmul2      = 'real';                    
    borne.nbicdmul2     = [0.1,10];  
    defaut.nbicdmul2    = 1;                
    info.nbicdmul2      = 'second NB injector: multiplication factor applied to NBI current drive efficiency';
    section.nbicdmul2   = 'NBI/NBICD@2';

    valeur.lhmode    = 2;
    type.lhmode      = 'real';
    borne.lhmode     = {0,1,2,3,4,5};
    defaut.lhmode    = 2;
    info.lhmode      = 'LHCD efficiency model:\nif = 0, Fisch like law  when wlh is defined, otherwise ITER basis scaling (Constant * <Te>);\nif = 1, adjusted to fit of chosen value of Vloop (vref);\nif = 2,  fixed value (etalh);\nif = 3,  Goniche scaling law;\nif = 4, simulTS scaling (for Tore Supra only);\nif = 5, this is use to describe a second ECCD system instead of LHCD system';
    section.lhmode   = 'LHCD';

    valeur.upshiftmode    = 'newmodel retune';
    type.upshiftmode      = 'string';                    
    borne.upshiftmode     = {'newmodel','newmodel + tail','newmodel retune','linear','1/q','Bpol','x^2','sqrt(x)','null','step@edge'};  
    defaut.upshiftmode    = 'newmodel retune';     
    info.upshiftmode      = 'parallel refractive index upshift model:\nif = "newmodel", then use the new formulation taking into account ALOHA/C3PO/LUKE  2012/2013 results;\nif = "newmodel + tail", same as "newemodel" with tail model effect;\nif = linear, upshift increases linearly from napr0 at the edge to fupshift*npar0 at plasma center;\nif = 1/q, npar proportional to 1/q;\n = Bpol, npar proportional to Bpol;\n = x^2,  npar proportional to x^2;\nif=sqrt(x), npar proportional to sqrt(x);\nif = null, no upshift;\nif = step@edge, npar becomes npar0 + fupshift';	
    section.upshiftmode   = 'LHCD';
    mode.upshiftmode      = 'advanced';

    valeur.fupshift        = 1;
    type.fupshift          = 'real';
    borne.fupshift     	   = [0,10];
    defaut.fupshift        = 1;
    info.fupshift          = 'parallel refractive index upshift model: parameter for upshift model.\nWhen upshiftmode = ''newmodel'', then factor applied to kinetic resonance position: n_par_Landau = fupshift * 6.5 / sqrt(Te).\nIn this case, for backward compatibility, if fupshift=0, then fupshift is reset internally to 1';
    section.fupshift       = 'LHCD';
    mode.fupshift      = 'advanced';

    valeur.etalh    = 0.8;
    type.etalh      = 'real';                    
    borne.etalh     = [-3e19,3e19];  
    defaut.etalh    = 0.8;                
    info.etalh      = 'LHCD efficiency or directivity:\nif lhmode = 2, value of normalised LHCD efficiency (etaLH in A/Wm^2);\notherwise launcher directivity defined as the fraction of total LH power in the co-current peak;\nif lhmode = 5,   multiplication factor applied to ECRH current drive efficiency for the second EC system (gives also sign of current source: if > 0, co-current and if < 0 counter-current)';	
    section.etalh   = 'LHCD';

    valeur.npar0    = 2;
    type.npar0      = 'real';
    borne.npar0     = [1,10];
    defaut.npar0    = 2;
    info.npar0      = 'launched parallel refractive index of LH at antenna';
    section.npar0   = 'LHCD';
    mode.npar0      = 'advanced';

    valeur.freqlh    = 3.7;
    type.freqlh      = 'real';
    borne.freqlh     = [1,10];
    defaut.freqlh    = 3.7;
    info.freqlh      = 'Lower Hybrid frequency (GHz)';
    section.freqlh   = 'LHCD';
    mode.freqlh      = 'advanced';

    valeur.wlh      = 0;
    type.wlh        = 'real';
    borne.wlh       = [0,3];
    defaut.wlh      = 0;
    info.wlh        = 'LH power deposition model:\nif = 0, LH source profile is computed with xlh et dlh;\notherwise wlh is the width of LH antenna active part (m);\nin this case, the source shape is computed with a simple model;\n if lhmode = 5, must be set to 0';
    section.wlh     = 'LHCD';
    mode.wlh        = 'advanced';
 
    valeur.xlh      = 0;
    type.xlh        = 'real';
    borne.xlh       = [0,0.8];
    defaut.xlh      = 0;
    info.xlh        = 'LH power deposition model: position of the maximum of LH current profile.\n(or position of the maximum of the second EC deposition profile, if lhmode = 5)';	
    section.xlh     = 'LHCD';

    valeur.dlh      = 0;
    type.dlh        = 'real';
    borne.dlh       = [0,0.7];
    defaut.dlh      = 0;
    info.dlh        = 'LH power deposition model: width of the LH current profile;\n(for Tore Supra, if = 0, use of profile from hard x-ray diagnostic; or width of second ECRH deposition profile if  lhmode = 5)';
    section.dlh     = 'LHCD';

    valeur.npar_neg    = 0;
    type.npar_neg      = 'real';
    borne.npar_neg     = [-10,0];
    defaut.npar_neg    = 0;
    info.npar_neg      = 'LH power deposition model:parallel refractive index of negative peak in the spectrum at the launcher; if = 0, used npar_neg = -npar0';
    section.npar_neg   = 'LHCD';
    mode.npar_neg      = 'advanced';

    valeur.angle_ece2    = 90;
    type.angle_ece2      = 'real';                    
    borne.angle_ece2     = {0,90,180};  
    defaut.angle_ece2    = 90;                
    info.angle_ece2      = 'second ECCD system efficiency - trapped particles effect: ECRH effectif deposition poloidal location: 0 -> LFS, 90 -> top, 180 -> HFS';
    section.angle_ece2   = 'LHCD';
    label.angle_ece2     = 'angle_ecrh2';
    mode.angle_ece2      = 'advanced';

    valeur.fwcd    = 0;
    type.fwcd      = 'integer';                    
    borne.fwcd     = {-1,0,1,2};
    defaut.fwcd    = 0;                
    info.fwcd      = 'ICRH heating scheme:\nif = -1, counter current FWCD;\nif = 0, minority ion heating;\nif =1, FWCD mode;\nif = 2, FW mode';
    section.fwcd   = 'ICRH/FW/FWCD';

    valeur.icrh_model    = 'PION_fit-Stix';
    type.icrh_model      = 'string';                    
    borne.icrh_model     = {'PION_fit-Stix','Dumont-Vu'};  
    defaut.icrh_model    = 'PION_fit-Stix';                
    info.icrh_model      = 'Selection of the ICRH model:\nif = PION_fit-Stix, deposition width fitted from PION +  Stix formulation for ion distribution function;\nif = Dumont-Vu, new model using resonance width and convergence on quaisilinear diffusion coefficient';	
    section.icrh_model   = 'ICRH/FW/FWCD';

    valeur.mino    = 'H';
    type.mino      = 'string';                    
    borne.mino     = {'H','He3','He4','T','D','B11'};  
    defaut.mino    = 'H';                
    info.mino      = 'ICRH heating scheme: minority species for ICRH scheme';
    section.mino   = 'ICRH/FW/FWCD';

    valeur.cmin    = 0.1;
    type.cmin      = 'real';                    
    borne.cmin     = [0,1];  
    defaut.cmin    = 0;                
    info.cmin      = 'ICRH heating scheme: fraction of the first minority species (nX/nD or nX/n_main if no deuterium in plasma discharge);\nif = 0, no suprathermal ions in plasma, direct heating of electrons and ions';	
    section.cmin   = 'ICRH/FW/FWCD';

    valeur.nphi    = 25;
    type.nphi      = 'real';                    
    borne.nphi     = [1,100];
    defaut.nphi    = 25;                
    info.nphi      = 'ICRH heating scheme: main toroidal wave number at ICRH launcher\n(n_phi ~ (2*pi*freq*R*n_par) / c with assumption K_phi ~ K_par)';	
    section.nphi   = 'ICRH/FW/FWCD';

    valeur.freq    = 57;
    type.freq      = 'real';                    
    borne.freq     = [1,1000];  
    defaut.freq    = 57;                
    info.freq      = 'ICRH heating scheme: ICRH frequency in MHz';
    section.freq   = 'ICRH/FW/FWCD';

    valeur.icrh_width   = 1;
    type.icrh_width     = 'real';                    
    borne.icrh_width    = [0.1,10];  
    defaut.icrh_width    = 1;                
    info.icrh_width      = 'ICRH heating scheme: multiplication factor for ICRH power deposition profile width (for PION_fit-Stix model)\nor vertical extention of ICRH power deposition in units of normalised radius (for Dumont-Vu model)';	
    section.icrh_width   = 'ICRH/FW/FWCD';

    valeur.fact_mino   = 0;
    type.fact_mino     = 'real';                    
    borne.fact_mino    = [0,10];  
    defaut.fact_mino    = 0;                
    info.fact_mino      = 'ICRH heating scheme (only used in PION_fit-Stix model): multiplication factor applied to formula giving fraction of minority ions accelerated by ICRH wave;\nif = 0, default value depending on species is used';	
    section.fact_mino   = 'ICRH/FW/FWCD';
    mode.fact_mino      = 'advanced';

    valeur.orbit_width   = 0;
    type.orbit_width     = 'real';                    
    borne.orbit_width    = {0,1};  
    defaut.orbit_width    = 0;                
    info.orbit_width      = 'Orbit width effect - with model Dumont-Vu:\nif = 0, no orbit width effect;\nif = 1, broaden the fast ion source profile according to analytical formula';	
    section.orbit_width   = 'ICRH/FW/FWCD';
    mode.orbit_width      = 'advanced';

    valeur.ifast_icrh  = 0;
    type.ifast_icrh    = 'real';                    
    borne.ifast_icrh   = [0,1];  
    defaut.ifast_icrh  = 0;                
    info.ifast_icrh    = 'Fast ion induced current: with model Dumont-Vu, multiplication factor applied to current due to fast minority ions; must be computed with the help of SPOT';	
    section.ifast_icrh = 'ICRH/FW/FWCD';
    mode.ifast_icrh      = 'advanced';
   
    valeur.fabs_fw  = 0;
    type.fabs_fw    = 'real';                    
    borne.fabs_fw   = [-1,1];  
    defaut.fabs_fw  = 0;                
    info.fabs_fw    = 'Direct electron power absorption in minority scheme: with model Dumont-Vu, abs(fabs_fw) = fraction of power absorbed directly by electrons;\nif > 0, fast wave current drive, otherwise if < 0, fast wave electron heating';	
    section.fabs_fw = 'ICRH/FW/FWCD';   
    mode.fabs_fw    = 'advanced';

    valeur.MC_onoff  = 'on';
    type.MC_onoff   = 'string';                    
    borne.MC_onoff  = {'on','off'};  
    defaut.MC_onoff = 'on';                
    info.MC_onoff   = 'Mode conversion: compute (on) or do not compute (off) fraction of input power lost due to mode conversion (IBW)';	
    section.MC_onoff = 'ICRH/FW/FWCD';   
    mode.MC_onoff    = 'advanced';

    valeur.fMC_loss  = 0;
    type.fMC_loss    = 'real';                    
    borne.fMC_loss   = [0,1];  
    defaut.fMC_loss  = 0;                
    info.fMC_loss    = 'Mode conversion: with model Dumont-Vu, fraction of power coming from mode conversion (IBW) that is lost (otherwise, this power is assume to heat the electron);\nUse the model from: L.G. Eriksson and T. Hellsten, Physica Scripta vol. 52, 70-79, 1995.';	
    section.fMC_loss = 'ICRH/FW/FWCD';   
    mode.fMC_loss    = 'advanced';

    valeur.sitb    = 0;
    type.sitb      = 'integer';
    borne.sitb     = {0,1,2,3};
    defaut.sitb    = 0;
    info.sitb      = 'ITB control:\nif = 0, no ITB;\nif = 1, allow ITB with null or negative magnetic shear;if = 2, same  as 1  + rotation effect on ITB;\n if = 3, same as 2 + MHD rational q effect';	
    section.sitb   = 'MHD & ITB';
    mode.sitb      = 'advanced';

    valeur.itb_sensitivity    = 1;
    type.itb_sensitivity      = 'real';
    borne.itb_sensitivity     = [0,10];
    defaut.itb_sensitivity    = 1;
    info.itb_sensitivity      = 'ITB control: sensitivity for the creation of the barrier [1]';
    section.itb_sensitivity  = 'MHD & ITB';
    mode.itb_sensitivity      = 'advanced';

    valeur.itb_slope_max    = 2;
    type.itb_slope_max      = 'real';
    borne.itb_slope_max     = [0,10];
    defaut.itb_slope_max    = 2;
    info.itb_slope_max      = 'ITB control: controls the maximum pressure gradient inside barrier [2];\nlarger value gives lager gradient';
    section.itb_slope_max   = 'MHD & ITB';
    mode.itb_slope_max      = 'advanced';

    valeur.itb_density        = 1;
    type.itb_density          = 'real';
    borne.itb_density         = [0,10];
    defaut.itb_density        = 1;
    info.itb_density          = 'ITB control: sensitivity for the barrier on density with NBI [1]';
    section.itb_density       = 'MHD & ITB';
    mode.itb_density          = 'advanced';
    
    valeur.tae     = 0;
    type.tae      = 'integer';
    borne.tae     = {0,1};
    defaut.tae    = 0;
    info.tae      = 'TAE control: if = 0 no TAE;\nif = 1, take into account TAE effect in alpha fusion power losses (decrease fast alpha pressure gradient)';
    section.tae   = 'MHD & ITB';
    mode.tae      = 'advanced';

    valeur.alpha_channeling    = 0;
    type.alpha_channeling      = 'real';
    borne.alpha_channeling     = [0,1];
    defaut.alpha_channeling    = 0;
    info.alpha_channeling      = 'factor of enhancement of ion heating due to alpha channeling:\nif = 0, no effect;\n if = 1, all power goes to ions';
    section.alpha_channeling   = 'MHD & ITB';
    mode.alpha_channeling      = 'advanced';

    valeur.smhd    = 100;
    type.smhd      = 'real';                    
    borne.smhd     = [-10,100];
    defaut.smhd    = 0;                
    info.smhd      = 'threshold for ideal no wall limit:\n if > 0, decreases confinement time when beta_N exceeds limit smhd (in percent);\nif = 100, no MHD limit;\nif = 0, threshold at 4*li;\nif < 0, threshold at abs(smhd)*li';
    section.smhd   = 'MHD & ITB';
    mode.smhd      = 'advanced';

    valeur.tmhd    = 0;
    type.tmhd      = 'real';
    borne.tmhd     = [0,inf];
    defaut.tmhd    = 0;
    info.tmhd      = 'threshold for ideal no wall limit: first time ideal no wall limit allowed (s)';	
    section.tmhd   = 'MHD & ITB';
    mode.tmhd      = 'advanced';

    valeur.rip    = 0;
    type.rip      = 'integer';                    
    borne.rip     = {0,1};  
    defaut.rip    = 0;                
    label.rip     = 'ripple';                
    info.rip      = 'Ripple (Tore Supra only):\nif = 0, no ripple effect;\nif = 1 , take in account the ripple in Tore Supra';
    section.rip   = 'Axisymmetry';
    mode.rip      = 'advanced';

    valeur.signe    = 1;
    type.signe      = 'integer';                    
    borne.signe     = {1,-1};  
    defaut.signe    = 1;                
    label.signe     = 'sign';                
    info.signe      = 'sign of toroidal field projected on plasma current direction';	
    section.signe   = 'Current diffusion & Equilibrium';

    valeur.laochange    = 1;
    type.laochange      = 'integer';
    borne.laochange     = {0,1};
    defaut.laochange    = 1;
    info.laochange      = 'Current diffusion coordinate:\nif = 0,current diffusion equation is solved using Lao coordinate (r/a);\nif = 1, current diffusion equation is solved using flux coordinate (rho)';
    section.laochange   = 'Current diffusion & Equilibrium';
    mode.laochange      = 'advanced';

    valeur.moments_mode  = 1;
    type.moments_mode    = 'integer';
    borne.moments_mode   = {0,1};
    defaut.moments_mode  = 1;
    info.moments_mode    = 'Choice of formula for moments (elongation K(x) and triangularity d(x) profiles) computation in equlibrium:\nif =0, leading order (K is constant and d is lineary decresing from LCFS to magnetic axis);\nif =1, solves ODE for K(x) and d(x) (small inverse aspect ratio approximation)';
    section.moments_mode = 'Current diffusion & Equilibrium';
    mode.moments_mode    = 'advanced';
    
    valeur.morphing    = 5;
    type.morphing      = 'real';
    borne.morphing     = [0,15];
    defaut.morphing    = 5;
    info.morphing      = '2D equilibrium shape: exponent of morphing curve for the matching of LCFS\n(when LCFS is given by points, not used if LCFS is defined only by moments)';    
    section.morphing   = 'Current diffusion & Equilibrium';
    mode.morphing      = 'advanced';
   
    valeur.mode_expo_inte   = 1;
    type.mode_expo_inte     = 'integer';
    borne.mode_expo_inte    = {0,1};
    defaut.mode_expo_inte   = 1;
    info.mode_expo_inte     = 'Current diffusion, numerical method:\nif = 0,poloidal flux diffusion equation is solved using Crank-Nicholson method;\nif = 1, poloidal flux diffusion equation is solved using exponential integrator method';
    section.mode_expo_inte  = 'Current diffusion & Equilibrium';
    mode.mode_expo_inte     = 'advanced';
    
    valeur.protect_sepa_z0   = 'none';
    type.protect_sepa_z0     = 'string';
    borne.protect_sepa_z0    = {'none','minmax','rmax'};
    defaut.protect_sepa_z0   = 'none';
    info.protect_sepa_z0     = 'Allowing to vertically recenter LCFS given by point. LCFS must be centered in METIS; the vertical position is given by geo.z0. Possible choices are:\nnone = no correction;\nmimax = z0 iscomputed as (max(Z)+min(Z))/2;\nrmax = z0 is the value of Z where R is maximum on LCFS';
    section.protect_sepa_z0  = 'Current diffusion & Equilibrium';
    mode.protect_sepa_z0     = 'advanced';
    
    valeur.short_dt   = 'off';
    type.short_dt     = 'integer';
    borne.short_dt    = {'on','full','off'};
    defaut.short_dt   = 'off';
    info.short_dt     = 'Anti aliasing filter: if = on, switch on antialiasig filter to reduce numerical noise with short time step, applied to elected fields;\nif  = full, switch on antialiasig filter to reduce numerical noise with short time step, applied on every fields;\nif = off, no filtering is applied';
    section.short_dt  = 'Current diffusion & Equilibrium';
    mode.short_dt     = 'advanced';
   

    valeur.cronos_regul    = 0;
    type.cronos_regul      = 'integer';
    borne.cronos_regul     = {0,1,2,3,4,5};
    defaut.cronos_regul    = 0;
    label.cronos_regul     = 'q(0) regularisation';
    info.cronos_regul      = 'method used to compute q(0); q(0) is constrained to be above 0.5:\nif = 0, use simple extrapolation (metis default);\nif = 1, use analytic formula (cronos default);\nif = 2, average of two previous formulas;\nif = 3, remove limitation to be above 0.5 and extrapolate q(0) imposing null second derivative on magnetic axis (for use with time resolved sawtooth model);\nif = 4, minimal modification using physical assumptions;\nif = 5, change in poloidal flux near magnetic axis in order to use all formulas whitout extrapolation';
    section.cronos_regul   = 'Current diffusion & Equilibrium';

    valeur.equi_ppar    = 0;
    type.equi_ppar      = 'integer';
    borne.equi_ppar     = {0,1,2,3};
    defaut.equi_ppar    = 0;
    info.equi_ppar      = 'equilibrium solver: input pressure for Grad Shafranov equation:\nif = 0, perpendicular or isotropic pressure;\nif = 1, parallel pressure;\nif = 2, total pressure;\nif = 3, parallel pressure + rotational energy';
    section.equi_ppar   = 'Current diffusion & Equilibrium';
    mode.equi_ppar      = 'advanced';

    valeur.refined_ptot    = 0;
    type.refined_ptot      = 'integer';
    borne.refined_ptot     = {0,1};
    defaut.refined_ptot    = 0;
    info.refined_ptot      = 'method for computing suprathermal pressure profile:\nif = 0, suprathermal pressure profile is proportional to thermal profile pressure;\nif = 1, suprathermal pressure profile is the sum of profiles computed from suprathermal energy content (0D) taking the shape of source deposition (1D).';
    section.refined_ptot   = 'Current diffusion & Equilibrium';
    mode.refined_ptot      = 'advanced';

    valeur.cor_rel_spitzer    = 'off';
    type.cor_rel_spitzer      = 'string';
    borne.cor_rel_spitzer     = {'on','off'};
    defaut.cor_rel_spitzer    = 'off';
    info.cor_rel_spitzer      = 'if = on, applied relativistic correction to Spitzer resistivity, and by extension to the neoclassic resistivity, even if the computation for neoclassic transport is not available';
    section.cor_rel_spitzer   = 'Current diffusion & Equilibrium';
    mode.cor_rel_spitzer      = 'advanced';
    
    valeur.Kappa_xpoint   = 0;
    type.Kappa_xpoint     = 'float';
    borne.Kappa_xpoint    = [0,2];
    defaut.Kappa_xpoint   = 0;
    info.Kappa_xpoint     = 'if LCFS is defined by moments, and option.configuration = 2 or 3, value of elongation above which x-point is set if triangularity if sufficient (see  delta_xpoint) and\nplasma width is sufficient (see R_HFS_xpoint and R_LFS_xpoint);\n if = 0, this is not taken into account';
    section.Kappa_xpoint  = 'Current diffusion & Equilibrium';
    mode.Kappa_xpoint     = 'advanced';

    valeur.delta_xpoint   = 0;
    type.delta_xpoint     = 'float';
    borne.delta_xpoint    = [-1,1];
    defaut.delta_xpoint   = 0;
    info.delta_xpoint     = 'if LCFS is defined by moments, and option.configuration = 2 or 3, value of triangularity above (id >0) or under (if < 0) which x-point is set if elongation if sufficient (see  Kappa_xpoint) and\nplasma width is sufficient (see R_HFS_xpoint and R_LFS_xpoint);\n if = 0, this is not taken into account';
    section.delta_xpoint  = 'Current diffusion & Equilibrium';
    mode.delta_xpoint     = 'advanced';

    valeur.R_HFS_xpoint   = 0;
    type.R_HFS_xpoint     = 'float';
    borne.R_HFS_xpoint    = [0,Inf];
    defaut.R_HFS_xpoint   = 0;
    info.R_HFS_xpoint     = 'if LCFS is defined by moments, and option.configuration = 2 or 3, value of minimum LCFS radius above which x-point is set if condition on LFS radius is met and\n if shapping is sufficient (see Kappa_xpoint and delta_xpoint);\n if = 0, this is not taken into account';
    section.R_HFS_xpoint  = 'Current diffusion & Equilibrium';
    mode.R_HFS_xpoint     = 'advanced';

    valeur.R_LFS_xpoint   = 0;
    type.R_LFS_xpoint     = 'float';
    borne.R_LFS_xpoint    = [0,Inf];
    defaut.R_LFS_xpoint   = 0;
    info.R_LFS_xpoint     = 'if LCFS is defined by moments, and option.configuration = 2 or 3, value of maximum LCFS radius under which x-point is set if condition on HFS radius is met and\n if shapping is sufficient (see Kappa_xpoint and delta_xpoint);\n if = 0, this is not taken into account';
    section.R_LFS_xpoint  = 'Current diffusion & Equilibrium';
    mode.R_LFS_xpoint     = 'advanced';

    valeur.impur_rot    = 'imp';
    type.impur_rot      = 'string';
    borne.impur_rot     = {'imp','max'};
    defaut.impur_rot    = 'imp';
    info.impur_rot      = 'toroidal rotation, choice of impurity for rotation data output: \n imp = light impurity (charge zimp);\n max = impurity for radiation (charge zmax)';
    section.impur_rot   = 'Rotation';
    mode.impur_rot      = 'advanced';

    valeur.mode_vtheta    = 'Neoclassical V_pol';
    type.mode_vtheta      = 'string';
    borne.mode_vtheta     = {'Neoclassical V_pol','same v_tor'};
    defaut.mode_vtheta    = 'Neoclassical V_pol';
    info.mode_vtheta      = 'poloidal rotation - selection of the method used to compute poloidal rotation:\n''Neoclassical V_pol'' = use neoclassical formulation from ref: Y. B. Kim et all, Phys. Fluids. B 3  (8) 1991  p 2050-\n ''same v_tor'' = assume all species have the same toroidal rotation';
    section.mode_vtheta   = 'Rotation';
    mode.mode_vtheta      = 'advanced';

    valeur.rot_jr_loss    = 'off';
    type.rot_jr_loss      = 'string';
    borne.rot_jr_loss     = {'on','off'};
    defaut.rot_jr_loss    = 'off';
    info.rot_jr_loss      = 'Toroidal rotation:\nif = on, take into account forces due to backward current generated to compensate for for radial current induced by particles losses (fast and thermal);\nif = off, this effect is not included in the toroidal rotation computation';
    section.rot_jr_loss   = 'Rotation';
    mode.rot_jr_loss      = 'advanced';

    valeur.carnot    = 0.42;
    type.carnot      = 'real';                    
    borne.carnot     = [0,1];  
    defaut.carnot    = 0.42; 
    label.carnot     = 'thermodynamic efficiency';
    info.carnot      = 'Power plant: thermal power to electricity power conversion efficiency';	
    section.carnot   = 'Miscellaneous';
    mode.carnot      = 'advanced';

    valeur.mul_blanket  = 1.2;
    type.mul_blanket    = 'real';                    
    borne.mul_blanket   = [0.5,10];  
    defaut.mul_blanket  = 1.2;                
    info.mul_blanket    = 'Power plant: fusion power multiplication factor due to neutron multiplication in breeding blanket';	
    section.mul_blanket = 'Miscellaneous';
    mode.mul_blanket    = 'advanced';

    valeur.aux       = 0.05;
    type.aux         = 'real';                    
    borne.aux        = [0,1];  
    defaut.aux       = 0.05;                
    info.aux         = 'Power plant: fraction of electric power used by auxiliary systems other than additional heating sources';	
    section.aux      = 'Miscellaneous';
    mode.aux         = 'advanced';
  
    valeur.effinj    = 0.7;
    type.effinj      = 'real';                    
    borne.effinj     = [0,1];  
    defaut.effinj    = 0.7;                
    info.effinj      = 'Power plant: conversion efficiency of additional heating sources';	
    section.effinj   = 'Miscellaneous';
    mode.effinj      = 'advanced';
   
    valeur.available_flux    = Inf;
    type.available_flux      = 'real';                    
    borne.available_flux     = [0,Inf];  
    defaut.available_flux    = Inf;                
    info.available_flux      = 'available poloidal flux provided by central solenoid and poloidal field coils (Wb)';	
    section.available_flux   = 'Miscellaneous';

    valeur.machine    = '';
    type.machine      = 'string';                    
    borne.machine     = '';  
    defaut.machine    = '';                
    label.machine     = 'machine name';                
    info.machine      = 'name of the Device';	
    section.machine   = 'Miscellaneous';

    valeur.shot       = fix(rand(1) .* 1e6 -1);
    type.shot         = 'integer';                    
    borne.shot        = [0,1e6-1];  
    defaut.shot       = 0;                
    info.shot         = 'shot number or simulation identification number';	
    section.shot      = 'Miscellaneous';

    valeur.orientation       = 1;
    type.orientation         = 'integer';                    
    borne.orientation        = {-1,1};  
    defaut.orientation       = 1;                
    info.orientation         = 'orientation of toroidal magnetic field: follows ITER convention (for ITER the value is -1, i.e. clockwise for tokamak see from above)';	
    section.orientation      = 'Miscellaneous';
    mode.orientation         = 'advanced';

    valeur.reference_parameters 	= '';
    type.reference_parameters           = 'string';                    
    borne.reference_parameters          = '';  
    defaut.reference_parameters         = '';                
    info.reference_parameters           = 'file name that contains Parameters that will overwrite user defined Parameters (if empty, it is not used)';	
    section.reference_parameters        = 'Miscellaneous';
 
    valeur.first_wall 	      = '';
    type.first_wall           = 'string';                    
    borne.first_wall          = '';  
    defaut.first_wall         = '';                
    info.first_wall           = 'matfile name that contains (R,Z) points describing poloidal section of the first_wall.\nif empty, it is not used\nVariable names should be R and Z. R and Z must be vectors of same length';	
    section.first_wall        = 'Miscellaneous';

    valeur.dwdt_method    = 'implicit';
    type.dwdt_method      = 'string';
    % old n'est pas accessible, juste une compatibilite pour une version intermediare
    borne.dwdt_method     = {'implicit','explicit','none','v4.2','working_point'};
    defaut.dwdt_method    = 'implicit';
    info.dwdt_method      = 'Numerical method for evolution mode: method used to compute the time derivative of energy (dW / dt):\nif =  implicit, same method as in full shot simulation;\nif = explicit, explicit polynomial filtered method;\nif = none, dW / dt is set to 0;\nif = v4.2, explicit method used in previous versions of METIS;\n if = working_point, for Operation mode computation, method to find the steady state solution (all d/dt = 0 except dpsi/dt)';
    section.dwdt_method   = 'Convergence';
    mode.dwdt_method      = 'advanced';

    valeur.nbmax    = 31;
    type.nbmax      = 'integer';
    borne.nbmax     = [31 1000];
    defaut.nbmax    = 31;
    info.nbmax      = 'Numerical method: maximum number of convergence loops';
    section.nbmax   = 'Convergence';
    mode.nbmax      = 'advanced';

    valeur.tol0d    = 0;
    type.tol0d      = 'integer';
    borne.tol0d     = [0,0.1];
    defaut.tol0d    = 0;
    info.tol0d      = 'Numerical method: tolerance on relative error to stop the convergence loop;\nif = 0, use default values : 1e-2 for fast mode, 1e-3 for full run and evolution mode';
    section.tol0d   = 'Convergence';
    mode.tol0d      = 'advanced';

    % wrap tooltip
    noms = fieldnames(info);
    for k=1:length(noms)
          info.(noms{k}) = sprintf(info.(noms{k}));
    end

    interface.ts = '';                    % nom de la fonction d'interfacage avec les donnees TS
    interface.jet = '';                   % nom de la fonction d'interfacage avec les donnees Jet

    zs.valeur=valeur;
    zs.type=type;
    zs.borne=borne;
    zs.defaut=defaut;
    zs.info=info;
    zs.mode  = mode;
    zs.interface=interface;
    zs.multicol = 1;                      % force sur une colonne
    zs.label    = label;
    zs.section    = section;

    zs.description = 'calcul 0D du scenario';   % description (une ligne) de la fonction

    zs.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
    zs.gui  ='';                             % nom de l'interface graphique specifique si elle existe
    zs.controle = '';                        % nom de la fonction de controle des valeurs si elle existe

