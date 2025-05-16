% Experimental profile data for METIS 
% The presence of this data overhide all user choice made with METIS interface 
% For auxilliary heating and current drive sources, only the shape of the source is taken into account
% the total power is drive by the power reference
% the current source is propotionnal to the power reference time current drive efficiency deduce from the experimental data
% the current drive efficiency are computed using int(j*S')/int(p*V') and then use to compute the current source.
% Optionally, the current drive efficiency can be modulated by the line averaged density if reference value is provided.
%
% 1- Electron density 
% 	NE_EXP.temps =  time vector [N,1]
% 	NE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	NE_EXP.ne    =  electron density [N,M] in m^-3
%	setappdata(0,'NE_EXP',NE_EXP)
%
% 2- Electron temperature 
% 	TE_EXP.temps =  time vector [N,1]
% 	TE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	TE_EXP.te    =  electron temperature [N,M] in eV
%	setappdata(0,'TE_EXP',TE_EXP)
%
%    optionnal setting for Sawteeth:
%       TE_EXP.ton_trig  = start time for sawtooth trigger
%       TE_EXP.toff_trig = stop time for sawtooth trigger
%       TE_EXP.dte_trig  = trigger criterium
%       if (te0(k+1) - te0(k)) * 2 / (te0(k+1) + te0(k)) < dte_trig and q0 < q0_trigger then a sawtooth is triggered.
%       without dte_trig field in TE_EXP, METIS evolved as usual for sawtooth.
%
% 3- ion temperature
% 	TI_EXP.temps =  time vector [N,1]
% 	TI_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	TI_EXP.ti    =  ion temperature [N,M] in eV
%	setappdata(0,'TI_EXP',TI_EXP)
%
% 4- LHCD shape (scale on 0D data)
% 	LH_SHAPE_EXP.temps =  time vector [N,1]
% 	LH_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	LH_SHAPE_EXP.jlh   =  shape of current source[N,M] , if possible in A/m^2
% 	LH_SHAPE_EXP.plh   =  shape of power deposition  source[N,M]  , if possible in W/m^3
%
%       optional field :
%	 	LH_SHAPE_EXP.nbar   =  line averaged electron density [N,1] , if possible in m^-3 (unit must be compatible with jlh unit)
%
%	setappdata(0,'LH_SHAPE_EXP',LH_SHAPE_EXP)
%
% 5- ECCD shape (scale on 0D data)
% 	ECCD_SHAPE_EXP.temps   =  time vector [N,1]
% 	ECCD_SHAPE_EXP.x       =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	ECCD_SHAPE_EXP.jeccd   =  shape of current source[N,M] , if possible in A/m^2
% 	ECCD_SHAPE_EXP.peccd   =  shape of power deposition  source[N,M]  , if possible in W/m^3
%
%       optional field :
%	 	ECCD_SHAPE_EXP.nbar   =  line averaged electron density [N,1] , if possible in m^-3 (unit must be compatible with jlh unit)
%
%	setappdata(0,'ECCD_SHAPE_EXP',ECCD_SHAPE_EXP)
%
% 6- ICRH shape (scale on 0D data)
% 	ICRH_SHAPE_EXP.temps   =  time vector [N,1]
% 	ICRH_SHAPE_EXP.x       =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	ICRH_SHAPE_EXP.jfwcd   =  current source[N,M] , if possible in A/m^2
% 	ICRH_SHAPE_EXP.pel     =  shape of power deposition source on electron due to collision with fast ions, 
% 	                          same scale as pion [N,M]  , if possible in W/m^3
% 	ICRH_SHAPE_EXP.pfw     =  shape of power deposition source on electron due to fast wave absorption,
%                                 same scale as pion [N,M]  , if possible in W/m^3
% 	ICRH_SHAPE_EXP.pion    =  shape of power deposition source on electron, same scale as pel [N,M],
% 	                          if possible in W/m^3
%
%       optional field :
%	 	ICRH_SHAPE_EXP.nbar   =  line averaged electron density [N,1] , if possible in m^-3 (unit must be compatible with jlh unit)
%
%	setappdata(0,'ICRH_SHAPE_EXP',ICRH_SHAPE_EXP)
%
% 7- NBICD shape (scale on 0D data)
% 	NBICD_SHAPE_EXP.temps   =  time vector [N,1]
% 	NBICD_SHAPE_EXP.x       =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	NBICD_SHAPE_EXP.jnbicd  =  shape of current source[N,M] , if possible in A/m^2
% 	NBICD_SHAPE_EXP.pel     =  shape of power deposition source on electron, same scale as pion [N,M], if possible in W/m^3
% 	NBICD_SHAPE_EXP.pion    =  shape of power deposition source on electron, same scale as pel [N,M], if possible in W/m^3
%
%       optional field :
%	 	NBICD_SHAPE_EXP.nbar   =  line averaged electron density [N,1] , if possible in m^-3 (unit must be compatible with jlh unit)
%
%	setappdata(0,'NBICD_SHAPE_EXP',NBICD_SHAPE_EXP)
%
% 8- NBICD2 shape (scale on 0D data)
% 	NBICD2_SHAPE_EXP.temps   =  time vector [N,1]
% 	NBICD2_SHAPE_EXP.x       =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	NBICD2_SHAPE_EXP.jnbicd  =  shape of current source[N,M] , if possible in A/m^2
% 	NBICD2_SHAPE_EXP.pel     =  shape of power deposition source on electron, same scale as pion [N,M], if possible in W/m^3
% 	NBICD2_SHAPE_EXP.pion    =  shape of power deposition source on electron, same scale as pel [N,M], if possible in W/m^3
%
%       optional field :
%	 	NBICD2_SHAPE_EXP.nbar   =  line averaged electron density [N,1] , if possible in m^-3 (unit must be compatible with jlh unit)
%
%	setappdata(0,'NBICD2_SHAPE_EXP',NBICD2_SHAPE_EXP)
%
% 9- RUNAWAY  (scale on 0D data)
% 	RUNAWAY_EXP.temps  =  time vector [N,1]
% 	RUNAWAY_EXP.x      =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	RUNAWAY_EXP.jrun   =  shape of current source[N,M], if possible in A/m^2
% 	RUNAWAY_EXP.irun   =  total runaway current [N,1],  in A
% 	RUNAWAY_EXP.iohm   =  total ohmic current [N,1], in A
% 	RUNAWAY_EXP.ip     =  plasma current (projected on toroidal direction) [N,1], in A
%	setappdata(0,'RUNAWAY_EXP',RUNAWAY_EXP)
%
%
% 10- Line radiation
% 	PLINE_EXP.temps =  time vector [N,1]
% 	PLINE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	PLINE_EXP.prad  =  density of line radiative power [N,M] in W/m^3
%	setappdata(0,'PLINE_EXP',PLINE_EXP)
%
% 11- Toroidal rotation 
% 	VTOR_SHAPE_EXP.temps =  time vector [N,1]
% 	VTOR_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	VTOR_SHAPE_EXP.omega =  profile shape of solid like plasma toroidal rotation [N,M] if possible in turn/s
%	setappdata(0,'VTOR_SHAPE_EXP',VTOR_SHAPE_EXP)
%
% 12- Effective charge profile shape (the line averaged value is given by the reference) 
% 	ZEFF_SHAPE_EXP.temps =  time vector [N,1]
% 	ZEFF_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	ZEFF_SHAPE_EXP.zeff  =  profile shape of effective charge without He ashes coming from fusion reaction and without W contribution
%	setappdata(0,'ZEFF_SHAPE_EXP',ZEFF_SHAPE_EXP)
%
% 13- Shape density of tungsten (the concentration of tungsten at LCFS in
%     METIS is provided by METIS model, if option.Sn_fraction > 0, 
%     shape of the mixed tungsten and tin profile density)
% 	NW_SHAPE_EXP.temps =  time vector [N,1]
% 	NW_SHAPE_EXP.x     =  space vector [1,M], with M >= 3 , Lao coordinate r/a
% 	NW_SHAPE_EXP.nwp   =  profile shape of density of tungsten 
%	setappdata(0,'NW_SHAPE_EXP',NW_SHAPE_EXP)
%
% 14- Confinement times: any of the following quantities can be provided to
%     METIS as an external data. Fields that are not provided or are empty
%     is just skipped.
%
%       TAUE_EXP.temps                 = time vector [N,1]
%       TAUE_EXP.tau_bohm              = Bohm confinement time (s), vector [N,1]
%       TAUE_EXP.tauthl                = L-mode confinement time (s), vector [N,1]
%       TAUE_EXP.tauh                  = H-mode confinement time (s), vector [N,1] 
%       TAUE_EXP.tauhe_l               = Helium ashes L-mode confinement time (s), vector [N,1]
%       TAUE_EXP.tauhe_h               = Helium ashes H-mode confinement time (s), vector [N,1]
%       TAUE_EXP.fraction_closed_lines = fraction of closed line in plasma poloidal cross section (0 to 1) [N,1]
%	    setappdata(0,'TAUE_EXP',TAUE_EXP);
%
% 15- source of cold neutral at the LCFS (for the coupling with DYON)
%
%       NEUTRAL_EXP.temps      = time vector [N,1]
%       NEUTRAL_EXP.n0a        = number of neutral inputed at the LCFS by the recycling and the gas puff [N,1]
%	    setappdata(0,'NEUTRAL_EXP',NEUTRAL_EXP);
%
%





 


