% mise a jour automatique du  parametre nbeq_mode
function zmaj_nbeq_mode(nom_mode,ydata)

disp(nom_mode)

if any(ydata > 1)
	nbeq_mode = evalin('base','param.gene.nbeq_mode');
	switch nom_mode
		
	case {'data.mode.pe','data.mode.pion','data.mode.nel'}
		zassignin('base','param.gene.nbeq_mode', max(nbeq_mode,4));
		
	case 'data.mode.rot'
		zassignin('base','param.gene.nbeq_mode', max(nbeq_mode,5));
		
	case {'data.mode.fluce','data.mode.flucion'}
		zassignin('base','param.gene.nbeq_mode', max(nbeq_mode,7));
		
	end
end
