% ZUIDEDIT_MODE_FCT   gestion des callbacks du formulaire de mode 
%--------------------------------------------------------------
% fichier zuiedit_mode_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de mode 
%
% syntaxe :
%	zuiedit_mode_fct(action)
%
% entrees :
%  action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 11/09/2003.
% 
% liste des modifications : 
% * 11/12/2002 : interface en anglais
% * 11/09/2003 -> corrrection bug "toutes les" (infini)
%--------------------------------------------------------------
function zuiedit_mode_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback de zuiedit_mode_fct: ')
% disp(action)

% recupere le handle de la fenetre concernee
%[hfig,h] = zuiformhandle('mode') ;
% recupere le handle de la fenetre concernee
hfig = gcf;
h    = getappdata(hfig,'zhandle');
htag = get(hfig,'tag');
% information pour l'assistant
zuicr(hfig,action) ;

% On efface la ligne de commentaire
set(h.commentaire,'visible','off') ;

% Abscisses
tt = evalin('base','data.gene.temps') ;

% selon ation
switch lower(action)

case 'list_intervalle'
	% on selectionne un intervalle
	xd   = getappdata(h.list_intervalle,'xd') ;
	xf   = getappdata(h.list_intervalle,'xf') ;
	libel = getappdata(h.list_intervalle,'libel') ;
      	ind  = min(get(h.list_intervalle,'value')) ; 
	tdeb = xd(ind) ;
	tfin = xf(ind) ;

	% recherche des indices kmin et kmax 
	dt = abs(tt-tdeb) ;
	kmin = min(find(dt==min(dt))) ;
	dt = abs(tt-tfin) ;
	kmax = min(find(dt==min(dt))) ;

      	set(h.valeur_tmin,'string',tdeb) ;
      	set(h.valeur_tmax,'string',tfin) ;
      	set(h.valeur_kmin,'string',kmin) ;
      	set(h.valeur_kmax,'string',kmax) ;

	% trace de l'intervalle
	ydata = get(h.trait_intrvl,'ydata') ;
	set(h.trait_intrvl,'xdata',[tdeb tfin],'ydata',ydata, ...
	    'linewidth',3,'visible','on','linestyle','-') ;

case 'radio_modif'
	% Modifications des intervalles
	hout = zuiedit_param_intervalle ;

        %drawnow
        %while(ishandle(hout))
        %     pause(1)
        %end
        zuicache(hfig) ;
        zwaitfor(hout) ;
        zuivisible(hfig) ;
        
        %set(hout,'windowstyle','modal')
        %uiwait(hout)	

	% mise a jour de la liste des intervalles
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

	zuilisteintrvl(htag) ;
	%zuilisteintrvl('mode') ;

	zuireset(h.radio_modif) ;

case 'radio_import'
	hout = zuiedit_mode_import(hfig) ;

        zuicache(hfig) ;
        zwaitfor(hout) ;
        zuivisible(hfig) ;

	% disp('valeurs recuperées')
	y = getappdata(hfig,'y_import') ;
	% pause

	if ~isempty(y)	
		% tester si ces nouvelles valeurs sont compatibles avec
		% le mode en cours d'edition

		val = getappdata(hfig,'valeur') ;
		i = 1 ;
		cmd = strcat('ind = find(y~=val{',num2str(i),'}') ;
		for i=2,length(val)
			cmd = strcat(cmd,' & y~=val{',num2str(i),'}') ;
		end
		cmd = strcat(cmd,') ;') ;
		eval(cmd) ;

		if ~isempty(ind)	
			ss = sprintf(' unvalid value') ;
			set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
		else
			xdata = evalin('base','data.gene.temps') ;
		
			% on trace le nouveau graphique
			mode = getappdata(hfig,'ydata') ;
			mode{end+1} = y ;

			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;
		end
	end

	zuireset(h.radio_import) ;

case 'radio_xgrid'
	val = zuidata(h.radio_xgrid)
	if val==1
		set(gca,'xgrid','on') ;
	elseif val==0
		set(gca,'Xgrid','off') ;
	end

case 'pop_1'
	zuidata(h.text_com,' ') ;
	zuienable(h.edit_1) ;
	zuicache(h.text_com) ;

	quand = get(h.pop_1,'value') ;
	switch quand
		case 1
			% a l'indice
		case 2
			% au temps
		case 3
			% tous les temps
			set(h.text_com,'visible','on') ;
			zuidata(h.text_com,'pas') ;
		case 4
			% toutes les 
			set(h.text_com,'visible','on') ;
			zuidata(h.text_com,'s') ;
		case 5
			% sur tout l'intervalle 
			zuidisable(h.edit_1) ;
		otherwise
	end

case 'edit_1'
	%Temps ou indice de changement de valeur
	val = zuidata(h.edit_1) ;
	quand = get(h.pop_1,'value') ;
	switch quand
		case 1
		% a l'indice
			% prendre un integer
			val = round(val) ;

     			kmin = str2num(get(h.valeur_kmin,'string')) ;
      			kmax = str2num(get(h.valeur_kmax,'string')) ;
			if val > kmax
				val = kmax ;
				ss = sprintf(' index out of range [%g %g] ',kmin,kmax) ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end
			if val < kmin
				val = kmin ;
				ss = sprintf(' index out of range [%g %g] ',kmin,kmax) ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end
			%zuidata(h.edit_1,val) ;

		case 2
		% au temps
			tmin = str2num(get(h.valeur_tmin,'string')) ;
			tmax = str2num(get(h.valeur_tmax,'string')) ;
			if val > tmax
				val = tmax ;
				ss = sprintf(' time out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
			end
			if val < tmin
				val = tmin ;
				ss = sprintf(' time out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
			end
			zuidata(h.edit_1,val) ;
				
		case 3
		% tous les x pas de temps
		case 4
		% toutes les x secondes
		case 5
		% sur tout l'intervalle 
		otherwise
	end	

case 'pop_2'
	ind = get(h.pop_2,'value') ;
	val = get(h.pop_2,'userdata') ;
	val = val{ind} ;
	s = sprintf(' valeur => %g \n ',val) ; disp(s)

case 'retour'
	xdata = evalin('base','data.gene.temps') ;

	zoom out
	% effacer le graphique
	y = ones(1,length(xdata)) * NaN ;
	zuiedit_mode_plot(xdata,y) ;

	% on recupere les modes 
	mode = getappdata(hfig,'ydata') ;

	% on 'efface' le dernier
	if length(mode)~=1
		mode(end) = [] ;
	end
	% on memorise les modes
	setappdata(hfig,'ydata',mode) ;
	
	y = mode{end} ;
	zuiedit_mode_plot(xdata,y) ;
	
	zuireset(h.retour) ;

case 'go'
	xdata = evalin('base','data.gene.temps') ;
	y = [];

	mode = getappdata(hfig,'ydata') ;

	value = get(h.pop_2,'value') ;
	val   = get(h.pop_2,'userdata') ;
	val   = val{value} ;
	quand = get(h.pop_1,'value') ;

	switch quand
		case 1
			% a l'indice
			ind= zuidata(h.edit_1) ;
			if isempty(ind)| ind==' '
				ss = sprintf(' index unknown') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				zuireset(h.go) ;
				return
			end
			
			% prendre un integer
			ind = round(ind) ;

     			kmin = str2num(get(h.valeur_kmin,'string')) ;
      			kmax = str2num(get(h.valeur_kmax,'string')) ;
			if ind > kmax
				ind = kmax ;
				ss = sprintf(' index out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end
			if ind < kmin
				val = kmin ;
				ss = sprintf(' index out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end

			% on efface le graphique précédent
			y = ones(1,length(xdata)) * NaN ;
			zuiedit_mode_plot(xdata,y) ;

			% on trace le nouveau graphique
			zuidata(h.edit_1,ind) ;

			y      = mode{end} ;
			y(ind) = val ;
			mode{end+1} = y ;
			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;

		case 2
			% au temps
			tps = zuidata(h.edit_1) ;

			tmin = str2num(get(h.valeur_tmin,'string')) ;
			tmax = str2num(get(h.valeur_tmax,'string')) ;
			if tps > tmax
				tps = tmax ;
				ss = sprintf(' time out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end
			if tps < tmin
				tps = tmin ;
				ss = sprintf(' time out of range') ;
				set(h.commentaire,'string',ss,'visible','on','foregroundcolor','red') ;
				return
			end

			% on efface le graphique précédent
			y = ones(1,length(xdata)) * NaN ;
			zuiedit_mode_plot(xdata,y) ;

			% on trace le nouveau graphique
			zuidata(h.edit_1,tps) ;

			d   = abs(xdata-tps) ;
			ind = find(d==min(d)) ;
			y   = mode{end} ;
			y(ind) = val .* ones(1,length(ind)) ;
			mode{end+1} = y ;
			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;

		case 3
			% on efface le graphique précédent
			y = ones(1,length(xdata)) * NaN ;
			zuiedit_mode_plot(xdata,y) ;

			% tous les x pas
			dk = zuidata(h.edit_1) ;
			% prendre un integer
			dk = round(dk) ;

    			kmin = str2num(get(h.valeur_kmin,'string')) ;
      			kmax = str2num(get(h.valeur_kmax,'string')) ;

			ind = kmin:dk:kmax ;		

			% on trace le nouveau graphique
			y = mode{end} ;
			y(ind) = val .* ones(1,length(ind)) ;
			mode{end+1} = y ;
			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;

		case 4
			% on efface le graphique précédent
			y = ones(1,length(xdata)) * NaN ;
			zuiedit_mode_plot(xdata,y) ;

			% toutes les x secondes

    			kmin = str2num(get(h.valeur_kmin,'string')) ;
      			kmax = str2num(get(h.valeur_kmax,'string')) ;
			ind = kmin:kmax ;	
		
			tmin = str2num(get(h.valeur_tmin,'string')) ;
			tmax = str2num(get(h.valeur_tmax,'string')) ;
			dt  = str2num(get(h.edit_1,'string')) ;
			u   = tmin:dt:tmax ;
                        t   = xdata(ind) ;
                        
                        ind = round(interp1(t,ind,u,'nearest')) ;

                        % probleme d'arrondi en debut d'intervalle
                        if isempty(ind)
                           disp('inconsitance intervalle ')
                           return
                        end
                        if ~isfinite(ind(1))
                           ind(1) = kmin;
                        end
                        ind(~isfinite(ind)) = [];
			% on trace le nouveau graphique
			y = mode{end} ;    
			y(ind) = val .* ones(1,length(ind)) ;
			mode{end+1} = y ;
			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;

		case 5
			% on efface le graphique précédent
			y = ones(1,length(xdata)) * NaN ;
			zuiedit_mode_plot(xdata,y) ;
			
    			kmin = str2num(get(h.valeur_kmin,'string')) ;
      			kmax = str2num(get(h.valeur_kmax,'string')) ;
			ind = kmin:kmax ;	
		
			% on trace le nouveau graphique
			y = mode{end} ;    
			y(ind) = val .* ones(1,length(ind)) ;
			mode{end+1} = y ;
			setappdata(hfig,'ydata',mode) ;

			zuiedit_mode_plot(xdata,y) ;

		
		otherwise
	end

	% On efface la ligne de commentaire
	set(h.commentaire,'visible','off') ;

	zuireset(h.go) ;

case {'annulation','close'}
	zuicloseone(hfig) ;	
	
case {'init','raz'}
	zuiformvisible(hfig) ;
	zuiformreset(hfig) ;
	zuiuploadform(hfig);
	
	zuilisteintrvl(htag) ;
	%zuilisteintrvl('mode') ;
	mode = getappdata(hfig,'ydata') ;
	ydata = mode{1} ;
	zuiedit_mode_plot(tt,ydata) ;

	zuireset(h.raz) ;
	
case {'validation','close'}
	% on recupere les modes
	mode = getappdata(hfig,'ydata') ;
	ydata = mode{end} ;
	nom_mode = get(h.nom_mode,'string') ;
	zassignin('base',nom_mode,ydata) ;
	
	% mise a jour automatique du  parametre nbeq_mode
	zmaj_nbeq_mode(nom_mode,ydata)
	
	zuiformcache(hfig) ;
	zuidownloadform(hfig) ;

	zuisavenonok ;
	zuireset(h.validation) ;
	
otherwise
	warning('action not taken into account')
	   
end

