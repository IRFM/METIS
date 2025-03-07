% ZUIDEDIT_PARAM_INTERVALLE_FCT  gestion des callbacks du formulaire des intervalles
%--------------------------------------------------------------
% fichier zuiedit_param_intervalle_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	des intervalles , parametres  g��aux
%	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_intervalle_fct(action)
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 24/09/2001.
% 
% liste des modifications : 
%  * 24/09/2001 -> correction petits bugs (ajouter, nom valide, affichage)
%
%--------------------------------------------------------------
function zuiedit_param_intervalle_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback de zuiedit_param_intervalle_fct: ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('intervalle');
% information pour l'assistant
zuicr(hfig,action) ;
if ~strcmp(action,'init')
	hoc = getfield(h,action)  ;
else
	hoc =[];
end
set(h.commentaire,'visible','off') ;

% selon ation
switch lower(action)

case 'list_intervalle'
	xd   = getappdata(h.list_intervalle,'xd') ;
	xf   = getappdata(h.list_intervalle,'xf') ;
	libel = getappdata(h.list_intervalle,'libel') ;
      	ind  = min(get(h.list_intervalle,'value')) ; 

	xdeb = xd(ind) ;
	xfin = xf(ind) ;
 	libel = libel(ind,:) ;
      	set(h.valeur_xdeb,'string',num2str(xdeb)) ;
      	set(h.valeur_xfin,'string',num2str(xfin)) ;
	set(h.nom_intervalle,'string',libel) ;

	% trace de l'intervalle
	set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
	    'linewidth',3,'visible','on','linestyle','-') ;
		
case 'valeur_xdeb'
	ind   = get(h.list_intervalle,'value') ;

	% on ne modifie pas les intervalles de temps et calcul
	if ind>2
		xdeb = zuidata(h.valeur_xdeb) ;
		xfin = zuidata(h.valeur_xfin) ;
		
		% on memorise la valeur precedente
		xd = getappdata(h.list_intervalle,'xd') ;
		xdprec = xd(ind) ;

		xlim = evalin('base','data.gene.temps') ;
		x1 = min(xlim) ;
		x2 = max(xlim) ;
		if xdeb < x1
			xdeb = x1 ;
		elseif xdeb >= x2
			xdeb = x1 ;
		end
	
		if xdeb >= xfin
			xdeb = xdprec ;
		end

		zuidata(h.valeur_xdeb,xdeb) ;

		% mise a jour de la liste
		xd(ind) = xdeb ;
		setappdata(h.list_intervalle,'xd',xd) ;
		zuilisteintrvl('intervalle') ;

		% trace de l'intervalle
		set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
	 	   'linewidth',3,'visible','on','linestyle','-') ;
	else
		xd = getappdata(h.list_intervalle,'xd') ;
		zuidata(h.valeur_xdeb,xd(ind)) ;	

		ss = sprintf(' time and computation interval cannot be modified') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		%pause(1) ;
	end

case 'select_xdeb'
	ind   = get(h.list_intervalle,'value') ;

	% on ne modifie pas les intervalles de temps et calcul
	if ind>2
		[xdeb,y] = ginput(1) ;
		xfin = zuidata(h.valeur_xfin) ;

		% on memorise la valeur precedente
		xd = getappdata(h.list_intervalle,'xd') ;
		xdprec = xd(ind) ;

		xlim = evalin('base','data.gene.temps') ;
		x1 = min(xlim) ;
		x2 = max(xlim) ;
		if xdeb < x1
			xdeb = x1 ;
		elseif xdeb >= x2
			xdeb = x1 ;
		end

		if xdeb >= xfin
			xdeb = xdprec ;
		end

		zuidata(h.valeur_xdeb,xdeb) ;

		% mise a jour de la liste
		xd(ind) = xdeb ;
		setappdata(h.list_intervalle,'xd',xd) ;
		zuilisteintrvl('intervalle') ;

		% trace de l'intervalle
		set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
		    'linewidth',3,'visible','on','linestyle','-') ;
	else
		ss = sprintf(' time and computation interval cannot be modified') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		%pause(1) ;
	end
	zuireset(h.select_xdeb) ;

case 'valeur_xfin'
	ind   = get(h.list_intervalle,'value') ;

	% on ne modifie pas les intervalles de temps et calcul
	if ind>2
		xfin = zuidata(h.valeur_xfin) ;
		xdeb = zuidata(h.valeur_xdeb) ;

		% on memorise la valeur precedente
		xf = getappdata(h.list_intervalle,'xf') ;
		xfprec = xf(ind) ;

		xlim = evalin('base','data.gene.temps') ;
		x1 = min(xlim) ;
		x2 = max(xlim) ;
		if xfin < x1
			xfin = x1 ;
		elseif xfin > x2
			xfin = x2 ;
		end
		
		if xfin <= xdeb
			xfin = xfprec ;
		end

		zuidata(h.valeur_xfin,xfin) ;

		% mise a jour de la liste
		xf(ind) = xfin ;
		setappdata(h.list_intervalle,'xf',xf) ;
		zuilisteintrvl('intervalle') ;

		% trace de l'intervalle
		set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
	 	   'linewidth',3,'visible','on','linestyle','-') ;
	else
		xf = getappdata(h.list_intervalle,'xf') ;
		zuidata(h.valeur_xfin,xf(ind)) ;

		ss = sprintf(' time and computation interval cannot be modified') ;
		%pause(1) ;
		set(h.commentaire,'visible','off') ;
	end

case 'select_xfin'
	ind   = get(h.list_intervalle,'value') ;

	% on ne modifie pas les intervalles de temps et calcul
	if ind>2
		[xfin,y] = ginput(1) ;
		xdeb = zuidata(h.valeur_xdeb) ;

		% on memorise la valeur precedente
		xf = getappdata(h.list_intervalle,'xf') ;
		xfprec = xf(ind) ;

		xlim = evalin('base','data.gene.temps') ;
		x1 = min(xlim) ;
		x2 = max(xlim) ;
		if xfin < x1
			xfin = x1 ;
		elseif xfin > x2
			xfin = x2 ;
		end

		if xfin <= xdeb
			xfin = xfprec ;
		end

		zuidata(h.valeur_xfin,xfin) ;
		zuireset(h.select_xfin) ;

		% mise a jour de la liste
		xf(ind) = xfin ;
		setappdata(h.list_intervalle,'xf',xf) ;
		zuilisteintrvl('intervalle') ;

		% trace de l'intervalle
		set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
	 	   'linewidth',3,'visible','on','linestyle','-') ;
	else
		ss = sprintf(' time and computation interval cannot be modified') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		%pause(1) ;
	end
	zuireset(h.select_xfin) ;

case 'radio_enchaine'
	ind   = get(h.list_intervalle,'value') ;

	% on ne supprime pas les intervalles de temps et calcul
	if ind>2
		xlim = evalin('base','data.gene.temps') ;
		xlim = max(xlim) ;
		xdeb = zuidata(h.valeur_xfin) ;
		zuidata(h.valeur_xdeb,xdeb) ;
		zuidata(h.valeur_xfin,xlim) ;

		% mise a jour de la liste
		xd = getappdata(h.list_intervalle,'xd') ;
		xf = getappdata(h.list_intervalle,'xf') ;
		xd(ind) = xdeb ;
		xf(ind) = xlim ;
		setappdata(h.list_intervalle,'xd',xd) ;
		setappdata(h.list_intervalle,'xf',xf) ;
		zuilisteintrvl('intervalle') ;

		% trace de l'intervalle
		set(h.trait,'xdata',[xdeb xlim],'ydata',[0 0], ...
		    'linewidth',3,'visible','on','linestyle','-') ;
	else
		ss = sprintf(' time and computation interval cannot be modified') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		%pause(1) ;
	end
	zuireset(h.radio_enchaine) ;

case {'nom_intervalle','ajouter'}
	% on recupere la valeur des champs
	var  = zuidata(h.nom_intervalle) ;
	var  = deblank(var) ;
	% securite nom valide
	ind  = find(var <= sprintf(' '));
	if ~isempty(ind)
		var(ind) = char('~'.* ones(1,length(ind)));
		var      = char(var);
	end
	if ~test_valide(var)
		ss = sprintf(' invalid name :  MatLab name is needed') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		zuidata(h.nom_intervalle,char('a' + min(length(getappdata(h.list_intervalle,'xd'))-1,25)));
		return
   end

	xdeb = zuidata(h.valeur_xdeb) ;
	xfin = zuidata(h.valeur_xfin) ;

	xd    = getappdata(h.list_intervalle,'xd') ;
	xf    = getappdata(h.list_intervalle,'xf') ;
	libel = getappdata(h.list_intervalle,'libel') ;

	% on regarde si le nouvel intervalle existe deja
	%istr = min(strmatch(var,libel)) ;
	istr = min(strmatch(var,libel,'exact')) ;
	nom  = deblank(libel(istr,:)) ;
	%if ~isempty(istr) & (strcmp(nom,var))
	if ~isempty(istr) 
		% on s'y positionne
		%xd(istr) = xdeb ;
		%xf(istr) = xfin ;
		xdeb      = xd(istr);
		xfin      = xf(istr) ;
	else
		% on le rajoute
		xd    = [xd xdeb] ;
		xf    = [xf xfin] ;
		libel = strvcat(libel,var) ;
		[istr,l] =size(libel) ;
	end
	if isempty(istr)
	% on le rajoute
		xd    = [xd xdeb] ;
		xf    = [xf xfin] ;
		libel = strvcat(libel,var) ;
		[istr,l] =size(libel) ;
	else
	% on s'y positionne
		xd(istr) = xdeb ;
		xf(istr) = xfin ;
	end
	
	setappdata(h.list_intervalle,'xd',xd) ;
	setappdata(h.list_intervalle,'xf',xf) ;
	setappdata(h.list_intervalle,'libel',libel) ;

	% mise a jour de la liste
	zuilisteintrvl('intervalle') ;
	set(h.list_intervalle,'value',istr) ;

	% trace de l'intervalle
	set(h.trait,'xdata',[xdeb xfin],'ydata',[0 0], ...
	    'linewidth',3,'visible','on','linestyle','-') ;


	zuireset(h.ajouter) ;
	
case 'supprimer'
   	ind   = get(h.list_intervalle,'value') ;

	% on ne supprime pas les intervalles de temps et calcul
	if ind>2
		xd    = getappdata(h.list_intervalle,'xd') ;
		xf    = getappdata(h.list_intervalle,'xf') ;
		libel = getappdata(h.list_intervalle,'libel') ;
	
		xd(ind) = [] ;
		xf(ind) = [] ;
		libel(ind,:)=[] ;
	
		setappdata(h.list_intervalle,'xd',xd) ;
		setappdata(h.list_intervalle,'xf',xf) ;
		setappdata(h.list_intervalle,'libel',libel) ;

		set(h.list_intervalle,'value',ind-1) ;

		zuidata(h.nom_intervalle,libel(ind-1,:)) ;
		zuidata(h.valeur_xdeb,xd(ind-1)) ;
		zuidata(h.valeur_xfin,xf(ind-1)) ;

		% mise a jour de la liste
		zuilisteintrvl('intervalle') ;
	else
		disp(' time and computation interval cannot be modified')
		ss = sprintf(' time and computation interval cannot be modifiedl') ;
		set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		%pause(1) ;
	end

	zuireset(h.supprimer) ;

case {'cons_1','cons_2','cons_3','cons_4','cons_5'}
	numstr = action(end) ;
	value  = get(getfield(h,action),'value') ;
	string = get(getfield(h,action),'string') ;
	tab    = get(getfield(h,action),'userdata') ;
	info   = tab{value} ;
	color  = get(getfield(h,action),'BackgroundColor') ;
	if ~isempty(info{3})
		xr = evalin('base',info{1}) ;
		yr = evalin('base',info{2}) ;
		if ~isnan(yr)
			% Chan : PB avec data.cons.flux vecteur 1x389
			% if size(yr,2) >1 
			if size(yr,2) >1 & size(yr,1)~=1
				yr = sum(yr,2) ;
			end
			% choix du multiplicateur
			yy    = abs(yr);
			ind   = find(yy > 0);
			if strcmp(info{3},'o')
				ind = find(yr >0);
				set(getfield(h,strcat('line_',numstr)),'xdata',xr(ind),'ydata',yr(ind), ...
				    'linestyle','none','visible','on','marker','o');
			else
				set(getfield(h,strcat('line_',numstr)),'xdata',xr,'ydata',yr, ...
				    'linestyle',info{3},'visible','on','marker','none','color',color);
			end	
		else
			ss = sprintf(' Pas de valeurs pour consigne : %s ',string(value,:)) ;
			set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
			%pause(1) ;
		end	
	else
		set(getfield(h,strcat('line_',numstr)),'visible','off');
	end
	
case {'asser_1','asser_2','asser_3','asser_4','asser_5'}
	numstr = action(end) ;
	value  = get(getfield(h,action),'value') ;
	string = get(getfield(h,action),'string') ;
	tab    = get(getfield(h,action),'userdata') ;
	info   = tab{value} ;
	color  = get(getfield(h,action),'BackgroundColor') ;
	if ~isempty(info{3})
		xr = evalin('base',info{1}) ;
		yr = evalin('base',info{2}) ;
		if ~isnan(yr)
			if size(yr,2) >1
				yr = sum(yr,2) ;
			end
			% choix du multiplicateur
			yy    = abs(yr) ;
			ind   = find(yy > 0) ;

			if strcmp(info{3},'o')
				ind = find(yr >0);
				set(getfield(h,strcat('line_',numstr)),'xdata',xr(ind),'ydata',yr(ind), ...
				    'linestyle','none','visible','on','marker','o') ;
			else
				set(getfield(h,strcat('line_',numstr)),'xdata',xr,'ydata',yr, ...
				    'linestyle',info{3},'visible','on','marker','none','color',color) ;
			end
		else
			ss = sprintf(' no feedback value : %s ',string(value,:)) ;
			set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
			%pause(1) ;
		end	
	else
		set(getfield(h,strcat('line_',numstr)),'visible','off') ;
	end

case {'annulation','close'}
	zuicloseone(hfig) ;

case {'init','raz'}
	zuiformvisible(hfig) ;
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;

	int = evalin('base','param.intervalle') ;
	name_int = fieldnames(int) ;
	xd=[] ;
	xf=[] ;
	name=[] ;
	for i=1:length(name_int)
		val = getfield(int,name_int{i}) ;
 		xd = [xd val(1)] ;
		xf = [xf val(2)] ;
		name = strvcat(name,name_int{i}) ;
	end
	setappdata(h.list_intervalle,'libel',name) ;
	setappdata(h.list_intervalle,'xd',xd) ;
	setappdata(h.list_intervalle,'xf',xf) ;
	set(h.nom_intervalle,'string',name(1,:)) ;

	% mise a jour de la liste
	zuilisteintrvl('intervalle') ;

	% on nettoie le trace de l'intervalle
	set(h.trait,'xdata',[0 0],'ydata',[0 0], ...
	    'linewidth',3,'visible','on','linestyle','-') ;

	zuireset(h.raz);
	
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig) ;

	% mise �jour de la structure param.intervalle
	libel = getappdata(h.list_intervalle,'libel') ;
	xd    = getappdata(h.list_intervalle,'xd') ;
	xf    = getappdata(h.list_intervalle,'xf') ;
        libelc = libel(1,:);
        libelc(libelc <= ' ') = [];
	int = struct(libelc,[xd(1) xf(1)]) ;
	[l,c]=size(libel) ;
	for i=2:l
                libelc = libel(i,:);
                libelc(libelc <= ' ') = [];
		int = setfield(int,libelc,[xd(i) xf(i)]) ;
	end
	zassignin('base','param.intervalle',int) ;

	zuisavenonok ;
	zuireset(h.validation) ;
	zuicloseone(hfig) ;
	
otherwise
	warning('action not taking into account')
	   
end

function cr = test_valide(var)

cr = 1;
eval(sprintf('%s = pi;',var),'cr =0;');
