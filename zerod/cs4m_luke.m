% cette fonction prend les donnees de cronos pour en faire des donnees externes de METIS
% elle peut etre utilisee pour reprendre des donnees experimentales
function  cs4m_luke(choice)



if nargin == 0
	choice ='';
end

% donnees pour le changement de coordonnees
try
	time   = evalin('base','post.profil0d.temps');
	time0d = evalin('base','post.zerod.temps');
catch
	warning('CS4M: No data for METIS in matlab workspace')
        return
end
if isempty(choice)
	hf = findobj(0,'type','figure','tag','METIS2LUKE_interface');
	if ~isempty(hf)
		figure(hf)
	else
		sections_liste = {'LUKE computation', ...
				'Use LUKE data in METIS for LHCD & ECCD', ...
				'Use LUKE data in METIS for LHCD', ...
				'Use LUKE data in METIS for ECCD', ...
				'Clear LUKE data in METIS'};
		commande_liste = {'cs4m_luke(''LUKE'');', ...
				'cs4m_luke(''LH&ECCD'');', ...
				'cs4m_luke(''LHCD'');', ...
				'cs4m_luke(''ECCD'');', ...
				'cs4m_luke(''CLEAR'');'};
		hf = zuimenuliste('METIS2LUKE',sections_liste,commande_liste);
		set(hf,'tag','METIS2LUKE_interface');
		drawnow
        end
        return
end

switch choice 
case 'LUKE'
	evalin('base','z0plotsc;');
        title('Select one time slice, please')
        drawnow
        [t,void] = ginput(1);
        close(gcf)
        if ~isempty(t) && isfinite(t)
		ind = find(time >= t,1);
                fprintf('computing LUKE @ t = %g s (id_zerod = %d, id_profile = %d)\n',time(ind),find(time0d >= t,1),ind);
                evalin('base',sprintf('post.lukeinmetis = metis2luke(post,%d);',ind));
                evalin('base',sprintf('post.lukeinmetis.computation_time_slice = %d;',ind));
                
                % verification
                evalin('base','test_lukeinmetis');
	end
case 'LH&ECCD'
	warning('Not yet available ... coming soon !');
	
case 'LHCD'

	lukeinmetis = evalin('base','post.lukeinmetis');
	% 4- LHCD shape (scale on 0D data)
	LH_SHAPE_EXP.x     =  evalin('base','post.profil0d.xli');
	LH_SHAPE_EXP.jlh   =  ones(length(time),1) * (lukeinmetis.j_TOT - lukeinmetis.j_OHM);
	LH_SHAPE_EXP.plh   =  ones(length(time),1) * (lukeinmetis.P_TOT - lukeinmetis.P_OHM);
	LH_SHAPE_EXP.temps =  time;
	setappdata(0,'LH_SHAPE_EXP',LH_SHAPE_EXP);
        % reglage des options de metis
        % passage a efficacite fixee
        evalin('base','z0dinput.option.lhmode = 2;');
        evalin('base','z0dinput.option.dlh = 0.2;');
        evalin('base','z0dinput.option.wlh = 0;');
        % calcul de l'efficacite
	zs = evalin('base','post.zerod');
	geo = evalin('base','post.z0dinput.geo');
        ind0d = find(time0d >= time(lukeinmetis.computation_time_slice),1);
        xlh          =  max(1,zs.plh(ind0d)) ./ zs.nbar(ind0d) ./ geo.R(ind0d);
        % Luke avec la synergie E//
        etalh_luke   =  (lukeinmetis.I_TOT - lukeinmetis.I_OHM) ./ xlh;
        % separation pour METIS
        etalh_1 = zs.etalh1(ind0d);
	etalh_luke_0 = etalh_luke - etalh_1;
        fprintf('etalh_luke = %g, etalh_nosynergy = %g, etalh_hot = %g  (A/W/m^2)\n',etalh_luke,etalh_luke_0,etalh_1);
        evalin('base',sprintf('z0dinput.option.etalh = %g;',etalh_luke_0));

case 'ECCD'

	lukeinmetis = evalin('post.lukeinmetis');
	% 4- LHCD shape (scale on 0D data)
	ECCD_SHAPE_EXP.x     =  evalin('base','post.profil0d.xli');
	ECCD_SHAPE_EXP.jeccd =  ones(length(time),1) * (lukeinmetis.j_TOT - lukeinmetis.j_OHM);
	ECCD_SHAPE_EXP.peccd =  ones(length(time),1) * (lukeinmetis.P_TOT - lukeinmetis.P_OHM);
	ECCD_SHAPE_EXP.temps =  time;
	setappdata(0,'ECCD_SHAPE_EXP',ECCD_SHAPE_EXP);
        % reglage des options de metis (uniquement l'efficacite pour le moment)
        ind0d = find(time0d >= time(lukeinmetis.computation_time_slice),1);
	zs = evalin('base','post.zerod');
        etaeccdmul = (lukeinmetis.I_TOT - lukeinmetis.I_OHM) ./ max(1,zs.ieccd(ind0d));
        fprintf('eccdmul = %g\n',etaeccdmul);
        evalin('base',sprintf('z0dinput.option.eccdmul = %g;',etaeccdmul));

case 'CLEAR'
	noms = fieldnames(getappdata(0));
	for k = 1:length(noms)
		if findstr(noms{k},'_EXP')
			rmappdata(0,noms{k});
			fprintf('removing external data in METIS for %s\n',strtok(noms{k},'_'));
		end
	end
otherwise
	% cas cancel
	return
end






