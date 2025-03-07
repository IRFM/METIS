% gestion des donnees externe selon contexte
function  cs4m(choice,nointer)

% test entree
if nargin == 0
  choice = '';
end

liste = {};
% selection du contexte
contexte_QLK = 0;
try
	time = evalin('base','QLK_nn_data.time');
	liste{end+1}  = 'Qualikiz NN';       
	contexte_QLK = 1;
end

contexte_cronos = 0;
try
	time = evalin('base','data.gene.temps');
	liste{end+1}  ='CRONOS';
        contexte_cronos = 1;
end

contexte_PROFILE_MAKER = 0;
try
	time = evalin('base','PROFILE_MAKER_data.time');
	liste{end+1}  = 'Profile Maker';       
	contexte_PROFILE_MAKER = 1;
end

contexte_IMAS = 0;
if isappdata(0,'IMAS_EXIST') && ~isempty(strcmp(getappdata(0,'IMAS_EXIST'),'YES'))
	liste{end+1}  = 'IMAS';       
	contexte_IMAS = 1;
end
if (length(liste) > 1) && (nargin < 2)
        rep = menu('Select a source for external data:',liste);
	if isempty(rep)
		return
	else
		item = liste{rep};
	end
elseif ~isempty(liste)
	item = liste{1};
else 
        warndlg('No external data source available','External data for METIS');
        try
	    type external_data_rule_for_METIS	    
	end
	return
end


switch item
case 'CRONOS'
		cs4m_cronos(choice);
case 'Qualikiz NN'
		cs4m_QLK(choice);
case 'Profile Maker'
		cs4m_PROFILE_MAKER(choice);
case 'IMAS'
		cs4m_imas(choice);
otherwise

	if contexte_cronos
		cs4m_cronos(choice);
	elseif contexte_QLK
		cs4m_QLK(choice);
	elseif contexte_PROFILE_MAKER
		cs4m_PROFILE_MAKER(choice);
	end
end
% donnees en provenance des predictions de QLK_nn
function  cs4m_QLK(choice)

itemliste = {'Ne','Te','Ti','Cancel'};
header    = sprintf('Copy of QLK_nn data in METIS external data structure\nWhich data or source must be copied ?');

while(isempty(strmatch(choice,itemliste)) | isempty(choice))
	Button = menu(header,itemliste);
	
	if Button == 0
		return
	end
	choice = itemliste{Button};
end
% lecture des donnees
data = evalin('base','QLK_nn_data');
time = data.time;
x    = data.x;

switch choice 
case 'Ne'
	% 1- Electron density 
	NE_EXP.x     =  x;
	NE_EXP.ne    =  data.NE;
	NE_EXP.temps =  time;
	setappdata(0,'NE_EXP',NE_EXP)
	
case 'Te'
	% 1- Electron temperature 
	TE_EXP.x     =  x;
	TE_EXP.te    =  data.TE;
	TE_EXP.temps =  time;
	setappdata(0,'TE_EXP',TE_EXP)

case 'Ti'
	% 1- Ion temperature
	TI_EXP.x     =  x;
	TI_EXP.ti    =  data.TI;
	TI_EXP.temps =  time;
	setappdata(0,'TI_EXP',TI_EXP)
	
otherwise
	% cas cancel
	return
end



% donnees en provenance des predictions de Porifle Maker
function  cs4m_PROFILE_MAKER(choice)

error('to be developped')


% cette fonction prend les donnees de cronos pour en faire des donnees externes de METIS
% elle peut etre utilisee pour reprendre des donnees experimentales
function  cs4m_cronos(choice)

if nargin == 0
	choice ='';
end

% donnees pour le changement de coordonnees
try
	time = evalin('base','data.gene.temps');
catch
	error('CS4M: No data for CRONOS in matlab workspace')
end
x  = evalin('base','param.gene.x');
ra = evalin('base','data.equi.a');
ra = ra ./ (ra(:,end) * ones(size(x)));

itemliste = {'Ne','Te','Ti','LHCD','ECRH','ICRH','NBI','NBI2','PLINE','RUNAWAY','Cancel'};
header    = sprintf('Copy of CRONOS data in METIS external data structure\nWhich data or source must be copied ?');

while(isempty(strmatch(choice,itemliste)) | isempty(choice))
	Button = menu(header,itemliste);
	
	if Button == 0
		return
	end
	choice = itemliste{Button};
end

switch choice 
case 'Ne'
	ne = evalin('base','data.prof.ne');
	% 1- Electron density 
	NE_EXP.x     =  x;
	indok = find(all(isfinite(ne),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ Ne : no valid data');
	end
	NE_EXP.ne    =  zeros(length(indok),length(x));
	NE_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		NE_EXP.ne(k,:)    =  interp1(ra(ind,:),ne(ind,:),x,'pchip','extrap');
		NE_EXP.temps(k)   =  time(ind);
	end
	setappdata(0,'NE_EXP',NE_EXP)
	
case 'Te'
	te = evalin('base','data.prof.te');
	% 1- Electron temperature 
	TE_EXP.x     =  x;
	indok = find(all(isfinite(te),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ Te : no valid data');
	end
	TE_EXP.te    =  zeros(length(indok),length(x));
	TE_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		TE_EXP.te(k,:)    =  interp1(ra(ind,:),te(ind,:),x,'pchip','extrap');
		TE_EXP.temps(k)   =  time(ind);
	end
	setappdata(0,'TE_EXP',TE_EXP)

case 'Ti'
	ti = evalin('base','data.prof.ti');
	% 1- Ion temperature
	TI_EXP.x     =  x;
	indok = find(all(isfinite(ti),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ Te : no valid data');
	end
	TI_EXP.ti    =  zeros(length(indok),length(x));
	TI_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		TI_EXP.ti(k,:)    =  interp1(ra(ind,:),ti(ind,:),x,'pchip','extrap');
		TI_EXP.temps(k)   =  time(ind);
	end
	setappdata(0,'TI_EXP',TI_EXP)

case 'LHCD'
	pel = evalin('base','data.source.hyb.el');
	js  = evalin('base','data.source.hyb.j');
	% 4- LHCD shape (scale on 0D data)
	LH_SHAPE_EXP.x     =  x;
	indok = find(all(isfinite(pel),2) & all(isfinite(js),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ LHCD : no valid data')
	end
	LH_SHAPE_EXP.jlh   =  zeros(length(indok),length(x));
	LH_SHAPE_EXP.plh   =  zeros(length(indok),length(x));
	LH_SHAPE_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		LH_SHAPE_EXP.plh(k,:)    =  interp1(ra(ind,:),pel(ind,:),x,'pchip','extrap');
		LH_SHAPE_EXP.jlh(k,:)    =  interp1(ra(ind,:),js(ind,:),x,'pchip','extrap');
		LH_SHAPE_EXP.temps(k)    =  time(ind);
	end
	setappdata(0,'LH_SHAPE_EXP',LH_SHAPE_EXP)
case 'ECRH'
	pel = evalin('base','data.source.fce.el');
	js  = evalin('base','data.source.fce.j');
	% 4- LHCD shape (scale on 0D data)
	ECCD_SHAPE_EXP.x     =  x;
	indok = find(all(isfinite(pel),2) & all(isfinite(js),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ ECCD : no valid data')
	end
	ECCD_SHAPE_EXP.jeccd   =  zeros(length(indok),length(x));
	ECCD_SHAPE_EXP.peccd   =  zeros(length(indok),length(x));
	ECCD_SHAPE_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		ECCD_SHAPE_EXP.peccd(k,:)    =  interp1(ra(ind,:),pel(ind,:),x,'pchip','extrap');
		ECCD_SHAPE_EXP.jeccd(k,:)    =  interp1(ra(ind,:),js(ind,:),x,'pchip','extrap');
		ECCD_SHAPE_EXP.temps(k)    =  time(ind);
	end
	setappdata(0,'ECCD_SHAPE_EXP',ECCD_SHAPE_EXP)
case 'ICRH'
	pel  = evalin('base','data.source.fci.el');
	pion = evalin('base','data.source.fci.ion');
	js   = evalin('base','data.source.fci.j');
	% 4- LHCD shape (scale on 0D data)
	ICRH_SHAPE_EXP.x     =  x;
	indok = find(all(isfinite(pel),2) & all(isfinite(pion),2) & all(isfinite(js),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ ICRH : no valid data')
	end
	ICRH_SHAPE_EXP.jfwcd   =  zeros(length(indok),length(x));
	ICRH_SHAPE_EXP.pel   =  zeros(length(indok),length(x));
	ICRH_SHAPE_EXP.pion   =  zeros(length(indok),length(x));
	ICRH_SHAPE_EXP.pfw   =  zeros(length(indok),length(x));
	ICRH_SHAPE_EXP.temps =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		ICRH_SHAPE_EXP.pel(k,:)    =  interp1(ra(ind,:),pel(ind,:),x,'pchip','extrap');
		ICRH_SHAPE_EXP.pion(k,:)   =  interp1(ra(ind,:),pion(ind,:),x,'pchip','extrap');
		ICRH_SHAPE_EXP.jfwcd(k,:)  =  interp1(ra(ind,:),js(ind,:),x,'pchip','extrap');
		ICRH_SHAPE_EXP.temps(k)    =  time(ind);
	end
	setappdata(0,'ICRH_SHAPE_EXP',ICRH_SHAPE_EXP)
case 'NBI'
	pel  = evalin('base','data.source.idn.el');
	pion = evalin('base','data.source.idn.ion');
	js   = evalin('base','data.source.idn.j');
	% 4- LHCD shape (scale on 0D data)
	NBICD_SHAPE_EXP.x     =  x;
	indok = find(all(isfinite(pel),2) & all(isfinite(pion),2) & all(isfinite(js),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ ICRH : no valid data')
	end
	NBICD_SHAPE_EXP.jnbicd  =  zeros(length(indok),length(x));
	NBICD_SHAPE_EXP.pel     =  zeros(length(indok),length(x));
	NBICD_SHAPE_EXP.pion    =  zeros(length(indok),length(x));
	NBICD_SHAPE_EXP.temps   =  zeros(length(indok),1);
	for k=1:length(indok)
		ind = indok(k);
		NBICD_SHAPE_EXP.pel(k,:)     =  interp1(ra(ind,:),pel(ind,:),x,'pchip','extrap');
		NBICD_SHAPE_EXP.pion(k,:)    =  interp1(ra(ind,:),pion(ind,:),x,'pchip','extrap');
		NBICD_SHAPE_EXP.jnbicd(k,:)  =  interp1(ra(ind,:),js(ind,:),x,'pchip','extrap');
		NBICD_SHAPE_EXP.temps(k)     =  time(ind);
	end
	setappdata(0,'NBICD_SHAPE_EXP',NBICD_SHAPE_EXP)
case 'NBI2'
      error('CS4M @ NBI2 : not yet implemented')

case 'PLINE'
	prad  = evalin('base','data.source.prad');
	% 4- LHCD shape (scale on 0D data)
	PLINE_EXP.x     =  x;
	indok = find(all(isfinite(prad),2) & all(isfinite(ra),2));
	if isempty(indok)
		error('CS4M @ PLINE : no valid data')
	end
	PLINE_EXP.prad          =  zeros(length(indok),length(x));
	for k=1:length(indok)
		ind = indok(k);
		PLINE_EXP.prad(k,:)     =  interp1(ra(ind,:),prad(ind,:),x,'pchip','extrap');
		PLINE_EXP.temps(k)     =  time(ind);
	end
	setappdata(0,'PLINE_EXP',PLINE_EXP)

case 'RUNAWAY'
      error('CS4M @ RUNAWAY : not yet implemented')

otherwise
	% cas cancel
	return
end






