% ZUIEDITCONS_ACTION    callback de l'editeur de consignes
%---------------------------------------------------------------------
% fichier zuieditcons_action.m ->  zuieditcons_action
%
% fonction Matlab 5 :
%	Cette fonction gere les "callback" de l'editeur de consignes
%
% syntaxe  :
%	zuieditcons_action(action);
%    
% entree :
%	action = chaine de caratere, nom de l'action
%
% sorties : aucune
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 15/09/2003.
% 
% liste des modifications : 
% 
%  * 28/08/2001 -> ajout du mode d'edition de basetemps
%  * 13/09/2001 -> correction bug mode 'degres'
%  * 21/09/2001 -> correction bug : ajout de l'action liee au multiplicateur
%  * 26/10/2001 -> correction bug si valeur initiale contient des NaN
%  * 22/11/2001 -> ajout: plot de la reference selon contexte
%  * 29/03/2002 -> ajout des codes de retour pour l'IdN
%  * 17/10/2002 -> ajout des codes de retour pour te et ti
%  * 15/09/22003 -> ajout du mode nhnd et ntnd
%  * 15/09/22003 -> ajout des code retour pour fce
%
%--------------------------------------------------------------
%
function zuieditcons_action(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end
% disp('callback cons: ')
% disp(action)

% recupere le handle de la fenetre concernee
hfig = gcf;
h    = getappdata(hfig,'zhandle');
% information pour l'assistant
zuicr(hfig,action);

% position du curseur
pp = get(h.axes_plot,'currentpoint');
xc = pp(1,1);
yc = pp(1,2);


% mode escalier
modessc = length(getappdata(hfig,'x_mem')) <= 2;

if ~verLessThan('matlab','8.6') || ispc
    % change action to souris if close to point
    % recherche du point selectionne
    switch action
    case 'list_plot'
	xd   = get(h.cons,'xdata');
	yd   = get(h.cons,'ydata');
	xlim = get(h.axes_plot,'xlim');
	dx   = abs(xlim(1) - xlim(2));
	ylim = get(h.axes_plot,'ylim');
	dy   = abs(ylim(1) - ylim(2));
	e    = min(((xd -xc) ./ dx).^2 + ((yd -yc) ./ dy) .^2);
	if e < 1e-3
	  action ='souris';
	end
    end
end

% selon ation
switch lower(action)
	
case {'init','raz'}
	zuiformvisible(hfig);
	
case {'annulation','close'}
	var_y       = getappdata(hfig,'var_y');
	switch var_y
	case 'data.prof.te'
		evalin('base','zrecalcultemppres(''el'',''clear'');');
	case 'data.prof.pe'
		evalin('base','zrecalcultemppres(''el'',''clear'');');
	case 'data.prof.ti'
		evalin('base','zrecalcultemppres(''ion'',''clear'');');
	case 'data.prof.pion'
		evalin('base','zrecalcultemppres(''ion'',''clear'');');
	end	
	zuiformcache(hfig);
	zuicloseone(hfig);
	
case 'valeur_multi'
	% mise a jour multiplicateur
	multi = zuidata(h.valeur_multi);
	if ~isempty(multi)
		if isfinite(multi)
			setappdata(hfig,'multi_mem',multi);
		else
			multi       = getappdata(hfig,'multi_mem');
			zuidata(h.valeur_multi,multi);
		end
	else
		multi       = getappdata(hfig,'multi_mem');
		zuidata(h.valeur_multi,multi);
	end
	
case 'validation'
	zuiformcache(hfig);
	%recupperation des donnees
	var_x       = getappdata(hfig,'var_x');
	var_y       = getappdata(hfig,'var_y');
	xx          = getappdata(hfig,'x_mem');
	canal       = getappdata(hfig,'canal');
	multi       = getappdata(hfig,'multi_mem');
	code_retour = getappdata(hfig,'code_retour');
	var_prop    = getappdata(hfig,'var_prop');
	xd          = get(h.cons,'xdata');
	yd          = get(h.cons,'ydata').*multi;
	modulation  = (get(h.modulation,'value') > 1);
	
	% calcul de la consigne sur la base temps de depart
	%spline_flag = get(h.spline,'value');
	spline_flag = zuidata(h.spline);
	if spline_flag == 0
		if modessc == 0
			yy = interp1(xd,yd,xx,'linear');
		else
			% mode edition de base-temps
			yy = yd;
			xx = xd;
		end
	elseif spline_flag == 1
		yy = zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx);
	else
		yy = zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx);
	end
   
   	% assignation des variables dans le workspace
	if ~isempty(var_x)
		zassignin('base',var_x,xx(:));
	end
	if ~isempty(var_y)
		% test de code_retour
		code = code_retour;
		yyold   = evalin('base',strcat(var_y,'(:,',int2str(canal),')'));
		indold  = find(~isfinite(yyold));
		if ~isempty(indold)
			yyold(indold) = zeros(1,length(indold));
		end
		switch code
		case 'LCFS'
                        disp('Remove LCFS given by points');
                        zclear0dsepa;
			code_retour = '';
		case 'abs'
			yy = yy .* exp(sqrt(-1) .* angle(yyold));
			code_retour = '';
                case 'pfce'
                        [ypuiss,ytoro,ypolo] = zdecodefce(yyold); 
                        yy = zcodefce(yy+eps,ytoro,ypolo); 
			code_retour = '';
                case 'toro'
                        [ypuiss,ytoro,ypolo] = zdecodefce(yyold); 
                        yy = zcodefce(ypuiss+eps,yy,ypolo); 
			code_retour = '';
                case 'polo'
                        [ypuiss,ytoro,ypolo] = zdecodefce(yyold); 
                        yy = zcodefce(ypuiss+eps,ytoro,yy); 
			code_retour = '';
 		case 'angle'
			yy = abs(yyold) .* exp(sqrt(-1) .* yy);
			code_retour = '';
		case 'degres'
			yy = yy / 180 * pi;
			yy = abs(yyold) .* exp(sqrt(-1) .* yy);
			code_retour = '';
		case 'idn'
			yy = yy + sqrt(-1) .* imag(yyold);
			code_retour = '';
		case 'idn2'
			yy = real(yyold) + sqrt(-1) .* yy;
			code_retour = '';
		case 'didndt'
			yy = real(yyold) + sqrt(-1) .* yy;
			code_retour = '';
		case 'nbar'
			yy = yy + sqrt(-1) .* imag(yyold);
			code_retour = '';
		case 'gaspuff'
			yy = real(yyold) + sqrt(-1) .* yy;
			code_retour = '';
		case 'nhnd'
			yy = yy + sqrt(-1) .* imag(yyold);
			code_retour = '';
		case 'nhe3nd'
			yy = yy + sqrt(-1) .* imag(yyold);
			code_retour = '';
		case 'ntnd'
			yy = real(yyold) + sqrt(-1) .* yy;
			code_retour = '';
		otherwise
			% rien
		end
		
		% cas de la consigne et du profil
		if isempty(var_prop)
			if modessc ~= 0
				zassignin('base',var_y,yy(:));
			else
				zassignin('base',strcat(var_y,'(:,',int2str(canal),')'),yy(:));
			end
		else
			% recupration de la modulation
			tm     = evalin('base',var_prop);
			% suppression des NaN et inf
			ind    = find(~isfinite(tm));
			if ~isempty(ind)
				tm(ind) = zeros(1,length(ind));
			end
			% sommation du  module sur l'espace/voies
			tm     = sum(abs(tm),2);
			
			% normalisation du max a 1
			tmax   = max(tm);
			if tmax == 0
				tm = ones(size(tm));
			else
				tm = tm ./ tmax;
			end
			
			% selon le cas constant ou module
			if modulation == 1
	            		yy = tm * yy(:)';		
			else
	            		yy = ones(size(tm)) * yy(:)';		
	      end
	            
	     	 	% assignation 
	     		zassignin('base',var_y,yy);
 
		end
	end
	
	% commande supplementaire
	if isempty(code_retour)
		zuisavenonok;
	elseif strcmp(code_retour,'complexe')
		zassignin('base',var_y,yy)
		zuicloseone(hfig)
		evalin('base','zuiedit_profcmplx_fct(''complexe'')') ;
		return
	elseif strcmp(code_retour,'retouche')
		zassignin('base',var_y,yy)
		zuicloseone(hfig)
		evalin('base','zuijet_retouche_fct(''retouche'')') ;
		return
	else
		evalin('base',code_retour,'disp(''Erreur dans code_retour'')');
   end
	
	switch var_y
	case 'data.prof.te'
		evalin('base','zrecalcultemppres(''el'');');
	case 'data.prof.pe'
		evalin('base','zrecalcultemppres(''el'');');
	case 'data.prof.ti'
		evalin('base','zrecalcultemppres(''ion'');');
	case 'data.prof.pion'
		evalin('base','zrecalcultemppres(''ion'');');
	end	

	
	% il n'y a pas d'autre variable ici
	zuicloseone(hfig);

case 'souris'
      	% recherche du point selectionne
	xd   = get(h.cons,'xdata');
	yd   = get(h.cons,'ydata');
	e    = (xd -xc).^2 + (yd -yc) .^2;
	ind  = min(find(e == min(e)));
	if isempty(ind)
	    	ind =1;
	end
	% marquage
	set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
	set(h.liste_valeurs,'value',ind);
	% valeur champ x et y
	set(h.valeur_X,'string',num2str(xd(ind)));
	set(h.valeur_Y,'string',num2str(yd(ind)));
	
	% activation de la pousuite
	set(hfig,'WindowButtonMotionFcn', 'zuieditcons_action(''deplacement'')', ...
	         'WindowButtonupFcn', 'zuieditcons_action(''point_fixe'')', ...
	         'pointer','crosshair');
		 
case 'deplacement'
     	xd   = get(h.cons,'xdata');
     	yd   = get(h.cons,'ydata');
    	ind  = get(h.marque,'userdata');
     	if ind == 1 
        	xmin = xd(1);
	     	xmax = xd(1);
     	elseif ind == length(xd)
        	xmin = xd(end);
	     	xmax = xd(end);
     	else
        	xmin = xd(ind-1)+eps;
	     	xmax = xd(ind+1)-eps;
     	end
     	xd(ind)  = max(xmin,min(xc,xmax));
     	hh=get(h.axes_plot,'title');
     	if (xd(ind) ~= xc)
     	   	set(hh,'string','Depassement des limites en X')
     	elseif ~isempty(get(hh,'string'))
     	   	set(hh,'string','')
     	end	
	yd(ind)  = yc;
	%set(h.marque,'xdata',xd(ind),'ydata',yd(ind));
	%set(h.cons,'xdata',xd,'ydata',yd);
	% champ x et y
	set(h.valeur_X,'string',num2str(xd(ind)));
	set(h.valeur_Y,'string',num2str(yd(ind)));
	%drawnow

case 'point_fixe'
     	% arret de la poursuite
	set(hfig,'WindowButtonMotionFcn', '','WindowButtonupFcn', '','pointer','arrow');
     	hh=get(h.axes_plot,'title');
     	set(hh,'string','')
     	%retracer de la courbe continue
     	xd   = get(h.cons,'xdata');
     	yd   = get(h.cons,'ydata');
     	ind  = get(h.marque,'userdata');
     	if ind == 1 
        	xmin = xd(1);
	     	xmax = xd(1);
     	elseif ind == length(xd)
        	xmin = xd(end);
	     	xmax = xd(end);
     	else
		xmin = xd(ind-1)+eps;
		xmax = xd(ind+1)-eps;
     	end
     	xd(ind)  = max(xmin,min(xc,xmax));
     	yd(ind)  = yc;
     	set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
     	set(h.cons,'xdata',xd,'ydata',yd);
     	%spline_flag = get(h.spline,'value');
	spline_flag = zuidata(h.spline);
     	if spline_flag == 0
     		if modessc == 0
     			set(h.interp,'xdata',xd,'ydata',yd);
     		else
     			[xesc,yesc] = zuiesc(xd,yd);
     			set(h.interp,'xdata',xesc,'ydata',yesc);
     		end
     	elseif spline_flag == 1
     		xx = linspace(min(xd),max(xd),1024);
     		set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
	else
     		xx = linspace(min(xd),max(xd),1024);
     		set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
     	end

      	% mise a jour de la liste
      	ind = 1:length(xd);
      	dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[ind;xd;yd]);
      	indice  = get(h.marque,'userdata');
      	set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',indice);
      	% champ x et y
      	set(h.valeur_X,'string',num2str(xd(indice)));
      	set(h.valeur_Y,'string',num2str(yd(indice)));

case 'liste_valeurs'
	xd   = get(h.cons,'xdata');
      	yd   = get(h.cons,'ydata');
      	ind  = min(get(h.liste_valeurs,'value'));      
      	set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
      	set(h.valeur_X,'string',num2str(xd(ind)));
      	set(h.valeur_Y,'string',num2str(yd(ind)));

case 'spline'
      	xd   = get(h.cons,'xdata');
      	yd   = get(h.cons,'ydata');
	spline_flag = zuidata(h.spline);
      	%spline_flag = get(h.spline,'value');
      	if spline_flag == 0
      		if modessc == 0
      			set(h.interp,'xdata',xd,'ydata',yd);
      		else
      			[xesc,yesc] = zuiesc(xd,yd);
      			set(h.interp,'xdata',xesc,'ydata',yesc);
      		end
		%set(h.interp,'xdata',xd,'ydata',yd);
      	elseif spline_flag == 1
           	xx=linspace(min(xd),max(xd),1024);
           	set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      	else
           	xx=linspace(min(xd),max(xd),1024);
           	set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      	end


case {'valeur_x', 'valeur_y'}

      	% lecture des valeurs
      	xc = str2num(get(h.valeur_X,'string'));
      	yc = str2num(get(h.valeur_Y,'string'));
      
      	% securite
      	if ~isempty(xc) & ~isempty(yc)
      		if isfinite(xc) & isfinite(yc)
      			% prise en compte de la modif ( avec bornes sur x)
      			xd   = get(h.cons,'xdata');
      			yd   = get(h.cons,'ydata');
      			ind  = get(h.marque,'userdata');
      			if ind == 1 
      				xmin = xd(1);
      				xmax = xd(1);
      			elseif ind == length(xd)
      				xmin = xd(end);
      				xmax = xd(end);
      			else
      				xmin = xd(ind-1)+eps;
      				xmax = xd(ind+1)-eps;
      			end
      			xd(ind)  = max(xmin,min(xc,xmax));
      			yd(ind)  = yc;
      			set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
      			set(h.cons,'xdata',xd,'ydata',yd);
      			% mise a jour de la liste
      			indice = 1:length(xd);
      			dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[indice;xd;yd]);
      			set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',ind);
      			% courbe interpolee
      			%spline_flag = get(h.spline,'value');
			spline_flag = zuidata(h.spline);
      			if spline_flag == 0
           		     		if modessc == 0
           		     			set(h.interp,'xdata',xd,'ydata',yd);
           		     		else
           		     			[xesc,yesc] = zuiesc(xd,yd);
           		     			set(h.interp,'xdata',xesc,'ydata',yesc);
           		     		end
           		     		%set(h.interp,'xdata',xd,'ydata',yd);
      			elseif spline_flag == 1
      				xx=linspace(min(xd),max(xd),1024);
      				set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      			else
      				xx=linspace(min(xd),max(xd),1024);
      				set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      			end
      		end
      	end

case {'ajouter','ajouter_souris'}

  	xd   = get(h.cons,'xdata');
   	yd   = get(h.cons,'ydata');
   	ind  = get(h.liste_valeurs,'value');      
   	% il doit y avoir plusieurs points
   	if ~isempty(ind)
        	% on complete si selection simple
		if length(ind) == 1
	   		if ind == length(xd)
	      			ind = [ind-1,ind];
	   		else
	      			ind = [ind,ind+1];
	   		end
        	end
		% calcul des nouveaux x
		xn = 0.5 .* (xd(ind(1:(end-1)))+xd(ind(2:end)));
		% nouveaux y
		%spline_flag = get(h.spline,'value');
		spline_flag = zuidata(h.spline);
      		if spline_flag == 0
            		yn = interp1(xd,yd,xn,'linear');
      		elseif spline_flag == 1
            		yn = zspline(xd,yd,xn);
      		else
            		yn = zpchip(xd,yd,xn);
     		end
		% creation des nouveaux vecteurs
		[xd,indice]= sort(cat(2,xd,xn));
		yy         = cat(2,yd,yn);
		yd         = yy(indice);
	
		% on met a jour le plot ...
		ind = min(ind);
		set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
     		set(h.cons,'xdata',xd,'ydata',yd);
        	% mise a jour de la liste
        	indice = 1:length(xd);
        	dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[indice;xd;yd]);
        	set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',ind);
		% courbe interpolee
        	%spline_flag = get(h.spline,'value');
		spline_flag = zuidata(h.spline);
      		if spline_flag == 0
      			if modessc == 0
      				set(h.interp,'xdata',xd,'ydata',yd);
      			else
      				[xesc,yesc] = zuiesc(xd,yd);
      				set(h.interp,'xdata',xesc,'ydata',yesc);
      			end
      			%set(h.interp,'xdata',xd,'ydata',yd);
      		elseif spline_flag == 1
           		xx=linspace(min(xd),max(xd),1024);
           		set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
		else
           		xx=linspace(min(xd),max(xd),1024);
           		set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
		end
      		set(h.valeur_X,'string',num2str(xd(ind)));
      		set(h.valeur_Y,'string',num2str(yd(ind)));

   	end
   	set(h.ajouter,'value',0);
   	
case 'ajouter_axes'

	xd   = get(h.cons,'xdata');
  	yd   = get(h.cons,'ydata');
	% creation des nouveaux vecteurs
	[xd,indice]= sort(cat(2,xd,xc));
	yy         = cat(2,yd,yc);
	yd         = yy(indice);
		
	% recherhce du nouvel indice 
	e    = (xd -xc).^2 + (yd -yc) .^2;
	ind  = min(find(e == min(e)));
	if isempty(ind)
		ind =1;
	end
	
% indice nom mise a jour ? 	
	% on met a jour le plot ...
	ind = min(ind);
	set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
     	set(h.cons,'xdata',xd,'ydata',yd);
        % mise a jour de la liste
        indice = 1:length(xd);
        dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[indice;xd;yd]);
        set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',ind);
	% courbe interpolee
        %spline_flag = get(h.spline,'value');
	spline_flag = zuidata(h.spline);
      	if spline_flag == 0
      		if modessc == 0
      			set(h.interp,'xdata',xd,'ydata',yd);
      		else
      			[xesc,yesc] = zuiesc(xd,yd);
      			set(h.interp,'xdata',xesc,'ydata',yesc);
      		end
      		%set(h.interp,'xdata',xd,'ydata',yd);
      	elseif spline_flag == 1
           	xx=linspace(min(xd),max(xd),1024);
           	set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      	else
           	xx=linspace(min(xd),max(xd),1024);
           	set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      	end
      	set(h.valeur_X,'string',num2str(xd(ind)));
      	set(h.valeur_Y,'string',num2str(yd(ind)));


case {'supprimer','supprimer_souris'}
   	xd   = get(h.cons,'xdata');
   	yd   = get(h.cons,'ydata');
   	ind  = get(h.liste_valeurs,'value');      
   	% il doit y avoir au moins une valeur
   	% et en laisser 2
   	if ~isempty(ind) & ((length(ind) +2) <= length(xd))&(min(ind)>1)&(max(ind)<length(xd))
        	xd(ind) =[];
		yd(ind)=[];
   	
		% on met a jour le plot ...
		ind = max(min(ind)-1,1);
		set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
     		set(h.cons,'xdata',xd,'ydata',yd);
        	% mise a jour de la liste
        	indice = 1:length(xd);
        	dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[indice;xd;yd]);
        	set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',ind);
		% courbe interpolee
        	%spline_flag = get(h.spline,'value');
		spline_flag = zuidata(h.spline);
      		if spline_flag == 0
      			if modessc == 0
      				set(h.interp,'xdata',xd,'ydata',yd);
      			else
      				[xesc,yesc] = zuiesc(xd,yd);
      				set(h.interp,'xdata',xesc,'ydata',yesc);
      			end
      			%set(h.interp,'xdata',xd,'ydata',yd);
      		elseif spline_flag == 1
           		xx=linspace(min(xd),max(xd),1024);
           		set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      		else
           		xx=linspace(min(xd),max(xd),1024);
           		set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
      		end
      		set(h.valeur_X,'string',num2str(xd(ind)));
      		set(h.valeur_Y,'string',num2str(yd(ind)));

   	end
   	set(h.supprimer,'value',0);

case 'list_plot'
	if strcmp(get(hfig,'SelectionType'),'alt')
		% on fait rien (c'est utilise par le menu contextuel)
	elseif strcmp(get(hfig,'SelectionType'),'open')
		% zoom out 
		set(h.axes_plot,'xlimmode','auto')
		set(h.axes_plot,'ylimmode','auto')
	elseif strcmp(get(hfig,'SelectionType'),'normal')
		% memorisation de la position du curseur
		xc_mem = xc;
		yc_mem = yc;
		res = rbbox;
		if (res(3) ~=0) & (res(4) ~=0)
			% nouvelle position du curseur
			pp = get(h.axes_plot,'currentpoint');
			xc = pp(1,1);
			yc = pp(1,2);
			set(h.axes_plot,'xlim',[min(xc,xc_mem),max(xc,xc_mem)]);
			set(h.axes_plot,'ylim',[min(yc,yc_mem),max(yc,yc_mem)]);
		end
	end
	
case 'grid'
	c = get(h.grid,'value');
	switch c
	case 1
		set(h.axes_plot,'xgrid','on','ygrid','off');
	case 2
		set(h.axes_plot,'xgrid','off','ygrid','on');
	case 3
		set(h.axes_plot,'xgrid','on','ygrid','on');
	case 4
		set(h.axes_plot,'xgrid','off','ygrid','off');
	end
	
case {'ref_1','ref_2','ref_3'}
	numstr = action(end);
	value  = get(getfield(h,action),'value');
	tab    = get(getfield(h,action),'userdata');
	info   = tab{value};
    if ~isempty(info{3})
        xr = evalin('base',info{1});
        yr = evalin('base',info{2});
        if size(yr,2) >1
            if ~isempty(findstr(info{2},'cons'))
                yr = sum(yr,2);
            else
                % petite svd
                [ur,kr,vr] = svd(yr,0) ;
                if all(size(xr) >1)
                    % rien
                elseif size(xr,1) == 1
                    yr = vr(:,1)';
                else
                    yr = ur(:,1);
                end
            end
        end
        if length(info)> 3
            switch info{4 }
                case 'real'
                    yr = real(yr);
                case 'imag'
                    yr = imag(yr);
                case 'norm'
                    yr = norm(yr);
            end
        end
        % choix du multiplicateur
		yy    = abs(yr);
		ind   = find(yy > eps);
        multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
        if isempty(multi)
            multi = 1;
        elseif ~isfinite(multi)
            multi = 1;
        end        
        
        if strcmp(info{3},'o')
			ind = find(yr >0);
			set(getfield(h,strcat('line_',numstr)),'xdata',xr(ind),'ydata',yr(ind)./multi, ...
			    'linestyle','none','visible','on','marker','o');
		else
			set(getfield(h,strcat('line_',numstr)),'xdata',xr,'ydata',yr./multi, ...
			    'linestyle',info{3},'visible','on','marker','none');
		end	
	else
		set(getfield(h,strcat('line_',numstr)),'visible','off');
	end
case 'node_1'	

	numstr = action(end);
	value  = get(getfield(h,action),'value');
	tab    = get(getfield(h,action),'userdata');
	info   = tab{value};
	if ~isempty(info{3})
		xr = evalin('base',info{1});
		yr = evalin('base',info{2});
		if size(yr,2) >1
		    if ~isempty(findstr(info{2},'cons'))
			yr = sum(yr,2);
	            else
		        % petite svd
			[ur,kr,vr] = svd(yr,0) ;
			if all(size(xr) >1)
			    % rien
			elseif size(xr,1) == 1
			    yr = vr(:,1)';
			else
			    yr = ur(:,1);
			end
		    end
		end
		% extraction des noeuds	
                if any(imag(yr))
		    yr = abs(yr);
		end
		nodes = zuieditcons_nodes(xr,yr,'',zuidata(h.spline) - 1);
		xd   = get(h.cons,'xdata');
		yd   = get(h.cons,'ydata');
		ind  = get(h.liste_valeurs,'value');      
                if ~isempty(nodes) && ~isempty(ind)
		     % boucle d'ajout des noeuds
                     for k=1:length(nodes)
			  xn = nodes(k);
                          % le point ne doit pas deja exister
                          if ~any(xn == xd)
			      % nouveaux y
			      spline_flag = zuidata(h.spline);
			      if spline_flag == 0
				      yn = interp1(xd,yd,xn,'linear');
			      elseif spline_flag == 1
				      yn = zspline(xd,yd,xn);
			      else
				      yn = zpchip(xd,yd,xn);
			      end
			      % creation des nouveaux vecteurs
			      [xd,indice]= sort(cat(2,xd,xn));
			      yy         = cat(2,yd,yn);
			      yd         = yy(indice);
	                      ind        = find(xd == xn,1);
			      % on met a jour le plot ...
			      set(h.marque,'xdata',xd(ind),'ydata',yd(ind),'userdata',ind);
			      set(h.cons,'xdata',xd,'ydata',yd);
			      % mise a jour de la liste
			      indice = 1:length(xd);
			      dd=sprintf('#%-3d  [%-6.5g, %-6.5g\t]|',[indice;xd;yd]);
			      set(h.liste_valeurs,'string',dd,'min',1,'max',length(xd),'value',ind);
			      % courbe interpolee
			      %spline_flag = get(h.spline,'value');
			      spline_flag = zuidata(h.spline);
			      if spline_flag == 0
				      if modessc == 0
					      set(h.interp,'xdata',xd,'ydata',yd);
				      else
					      [xesc,yesc] = zuiesc(xd,yd);
					      set(h.interp,'xdata',xesc,'ydata',yesc);
				      end
				      %set(h.interp,'xdata',xd,'ydata',yd);
			      elseif spline_flag == 1
				      xx=linspace(min(xd),max(xd),1024);
				      set(h.interp,'xdata',xx,'ydata',zspline(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
			      else
				      xx=linspace(min(xd),max(xd),1024);
				      set(h.interp,'xdata',xx,'ydata',zpchip(cat(2,xd(1)-mean(diff(xd)),xd),cat(2,yd(1),yd),xx));
			      end
			      set(h.valeur_X,'string',num2str(xd(ind)));
			      set(h.valeur_Y,'string',num2str(yd(ind)));

			 end
		     end
		end
	end

otherwise
      warning('ation non prise en compte')
end
%drawnow


% fonction plot escalier
function  [xesc,yesc] = zuiesc(xd,yd)

nx = length(xd);
indx = ceil(1:0.5:nx);
indy = fix(1:0.5:nx);

xesc = xd(indx);
yesc = yd(indy);
