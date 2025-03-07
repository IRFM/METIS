%
% convert data from metis to cpo scenario (mapping)
%
function scenario = mapscenario_imas(z0dstruct,data_zerod,profil0d,texte_diary,error_flag,scenario,run,occurrence,sigma_B0_eff,sigma_bvac_r)

% pour les donnees independantes du temps qui sont dependantes du temps dans le cpo
vt = ones(size(data_zerod.temps));


% calcul du fuelling
fuelling_data = metis_gaz_imas(z0dstruct,data_zerod,profil0d);
if ~isfield(scenario,'datainfo')
	scenario.datainfo = datainfo_empty;      
end
if isempty(scenario.datainfo.dataprovider) 
	switch z0dstruct.z0dinput.mode_exp
	case -2
		scenario.datainfo.dataprovider = 'METIS evolution';
	case -1
		scenario.datainfo.dataprovider = 'METIS simulation from scratch';	
	case 0
		scenario.datainfo.dataprovider = 'CRONOS';
	case 1
		scenario.datainfo.dataprovider = 'Tore Supra';	
	case 2 
		scenario.datainfo.dataprovider = 'JET';		
	case 11
		scenario.datainfo.dataprovider = 'Tore Supra (preparation)';	
	otherwise
		scenario.datainfo.dataprovider = getenv('USER');	
	end
end
if isempty(scenario.datainfo.putdate)
	scenario.datainfo.putdate      = sprintf('%f',clock2julday);
end
if isempty(scenario.datainfo.source)
	scenario.datainfo.source       = 'METIS';
end
if isempty(scenario.datainfo.comment)
	scenario.datainfo.comment      = 'ITM implementation of METIS';
end
if isempty(scenario.datainfo.isref)
	scenario.datainfo.isref        = 0;
end
%
if isempty(scenario.datainfo.whatref.user)
	scenario.datainfo.whatref.user = scenario.datainfo.dataprovider;
end
if isempty(scenario.datainfo.whatref.machine)
	scenario.datainfo.whatref.machine =  z0dstruct.z0dinput.machine;
end
if isempty(scenario.datainfo.whatref.shot)
	scenario.datainfo.whatref.shot    =  real(z0dstruct.z0dinput.shot(1));
end
if isempty(scenario.datainfo.whatref.run)
	scenario.datainfo.whatref.run    =  real(run);
end
if isempty(scenario.datainfo.whatref.occurrence)
	scenario.datainfo.whatref.occurrence    =  fix(str2num(occurrence));
end
scenario.datainfo.cocos = z0dstruct.z0dinput.option.COCOS;
%
scenario.centre.te0.value    = data_zerod.te0;
scenario.centre.ti0.value    = interp1_imas(profil0d.temps,profil0d.tip(:,1),data_zerod.temps,'pchip','extrap');
scenario.centre.ne0.value    = data_zerod.ne0;
scenario.centre.ni0.value    = data_zerod.ni0;
scenario.centre.shift0.value = data_zerod.d0;
scenario.centre.psi0.value   =  interp1_imas(profil0d.temps,profil0d.psi(:,1),data_zerod.temps,'pchip','extrap');
scenario.centre.phi0.value   = zeros(size(data_zerod.temps));
scenario.centre.q0.value     = data_zerod.q0;
scenario.centre.Rmag.value   = interp1_imas(profil0d.temps,profil0d.Raxe(:,1),data_zerod.temps,'pchip','extrap');
if length(scenario.centre.Rmag.value) == 1
	scenario.centre.Zmag.value   = z0dstruct.z0dinput.geo.z0(end);
else
	scenario.centre.Zmag.value   = z0dstruct.z0dinput.geo.z0;
end
scenario.centre.vtor_0.value = interp1_imas(profil0d.temps,profil0d.vtor(:,1),data_zerod.temps,'pchip','extrap');

%
scenario.configs.config.value = 0 * vt;
switch z0dstruct.z0dinput.option.configuration
case 0
	scenario.configs.config.value = 2  * vt;
case 1
	scenario.configs.config.value = 4  * vt;
case 2
	scenario.configs.config.value = 2 + 5 .* data_zerod.modeh;
case 3
	scenario.configs.config.value =  4 + 3 .* data_zerod.modeh;
case 4
	scenario.configs.config.value = 7  * vt;
end
scenario.configs.ecrh_freq.value 		= [];
scenario.configs.ecrh_loc.value 		= mean(z0dstruct.z0dinput.cons.xece) ;
scenario.configs.ecrh_mode.value 		= []; 
scenario.configs.ecrh_tor_ang.value 	= (z0dstruct.z0dinput.option.sens .* (pi / 8) ) * vt; 
scenario.configs.ecrh_pol_ang.value 	= (z0dstruct.z0dinput.option.angle_ece ./180 .* pi) * vt;
scenario.configs.ecrh_harm.value 		= [];
scenario.configs.enbi.value 		= z0dstruct.z0dinput.option.einj * vt ;
scenario.configs.r_nbi.value 		= z0dstruct.z0dinput.option.rtang * vt ;
scenario.configs.grad_b_drift.value 	= 1;
scenario.configs.icrh_freq.value 		= z0dstruct.z0dinput.option.freq *vt ;
scenario.configs.icrh_phase.value         = [];

if z0dstruct.z0dinput.option.fwcd == 2	
	scenario.configs.icrh_scheme = 'FW';	
elseif z0dstruct.z0dinput.option.fwcd == -1
	scenario.configs.icrh_scheme = 'FW_CCD';	
elseif z0dstruct.z0dinput.option.fwcd == 1
	scenario.configs.icrh_scheme = 'FW_CD';	
else
	scenario.configs.icrh_scheme = sprintf('%s_min_%d',z0dstruct.z0dinput.option.mino,ceil(mean(data_zerod.harm)));	
end

scenario.configs.LH_freq.value 		= (z0dstruct.z0dinput.option.freqlh .* 1e9) * vt ;
scenario.configs.LH_