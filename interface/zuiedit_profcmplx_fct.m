% ZUIEDIT_PROFCMPLX_FCT  	fonction de gestion des callback associes a
%			chaque uicontrol de l'editeur de profils complexes
%--------------------------------------------------------------
% fichier zuiedit_profcmplx_fct.m
% 
% fonction Matlab 5 :
%	fonction de gestion des callback associes a 
%	chaque uicontrol de l'editeur de profils complexes
%	cette fonction est utilisee par zuiedit_profcmplx.m
%
% syntaxe
%	 zuiprfcomplx_fct(action)
%
% entrees :
%  	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 6119
% version  1.7 , du 08/10/2001.
% 
% liste des modifications :
%  * 08/10/2001 -> correction effet plot, imagesc -> remplacement par
%                  un line et memorisation des proprietes (le plot reset les 
%                  proprietes des axes)
%  * 11/10/2001 -> ajout du trace de la donn� initiale avant recomposition
%  * 17/10/2002 -> ajout de l'edition de te,ti, pe et pion
%
%--------------------------------------------------------------

function zuiprfcomplx_fct(action)

if nargin ==0
	action = ' ' ;
end

% recupere le handle de la fenetre concernee
hfig = gcf ;
h    = getappdata(hfig,'zhandle') ;
htag = get(hfig,'tag') ;

% si la fenetre existe deja, on l'active
if ishandle(hfig)
	zuiformvisible(hfig)
end

% information pour l'assistant
zuicr(hfig,action) ;

%
info    = zinfo;

% recuperation des valeurs pas defaut
nom_var = getappdata(hfig,'nom_var') ;
yref    = getappdata(hfig,'yref') ;
k1def   = getappdata(hfig,'k1def') ;
k2def   = getappdata(hfig,'k2def') ;
k3def   = getappdata(hfig,'k3def') ;

v1 = evalin('base',['complexe.' nom_var '.v1']) ;
v2 = evalin('base',['complexe.' nom_var '.v2']) ;
v3 = evalin('base',['complexe.' nom_var '.v3']) ;
u1 = evalin('base',['complexe.' nom_var '.u1']) ;
u2 = evalin('base',['complexe.' nom_var '.u2']) ;
u3 = evalin('base',['complexe.' nom_var '.u3']) ;
k1 = evalin('base',['complexe.' nom_var '.k1']) ;
k2 = evalin('base',['complexe.' nom_var '.k2']) ;
k3 = evalin('base',['complexe.' nom_var '.k3']);

t = evalin('base','data.gene.temps') ;
x = evalin('base','param.gene.x') ;

% selon ation
switch lower(action)

% edit_k1
case 'edit_k1'
	k1l = str2num(get(h.edit_k1,'String')) ;
	% si k1l est vide; isempty(k1l) =1
 	if isempty(k1l)
 		set(h.edit_k1,'String',k1def) ;
		k1l = k2def ;
 	end
	if ~isfinite(k1l)
		set(h.edit_k1,'String',k1def) ;
		k1l = k2def ;
	end
	% pour eviter k1 ='1 0'
	if length(k1l) > 1
		set(h.edit_k1,'String',k1def) ;
		k1l = k2def ;
	end
	k1 = k1l ;


% edit_k2
case 'edit_k2'
	k2l = str2num(get(h.edit_k2,'String')) ;
	% si k2l est vide; isempty(k2l) =1
 	if isempty(k2l)
 		set(h.edit_k2,'String',k2def) ;
		k2l = k2def ;
 	end
	if ~isfinite(k2l)
		set(h.edit_k2,'String',k2def) ;
		k2l = k2def ;
	end
	if length(k2l) > 1
		set(h.edit_k2,'String',k2def) ;
		k2l = k2def ;
	end
	k2 = k2l ;


% edit_k3
case 'edit_k3'
	k3l = str2num(get(h.edit_k3,'String')) ;
	% si k3l est vide; isempty(k3l) =1
 	if isempty(k3l)
 		set(h.edit_k3,'String',k3def) ;
		k3l = k3def ;
 	end
	if ~isfinite(k3l)
		set(h.edit_k3,'String',k3def) ;
		k3l = k3def ;
	end
	if length(k3l) > 1
		set(h.edit_k3,'String',k3def) ;
		k3l = k3def ;
	end
	k3 = k3l ;


% edit_tps
case 'edit_tps'
	tt = str2num(get(h.edit_tps,'String')) ;

	% test:  tmin <= t <= tmax
	if isempty(tt)
		set(h.edit_tps,'string',num2str(get(h.slider_tps,'value'))) ;
		return
	elseif ~isfinite(tt)
		set(h.edit_tps,'string',num2str(get(h.slider_tps,'value'))) ;
		return
	end

	if get(h.slider_tps,'Min') > tt
		tt = get(h.slider_tps,'Min') ;
		set(h.edit_tps,'string',num2str(tt)) ;
	end
	if get(h.slider_tps,'max') < tt
		tt = get(h.slider_tps,'Max') ;
		set(h.edit_tps,'string',num2str(tt)) ;
	end
	set(h.slider_tps,'Value',tt) ;



% slider
case 'slider_tps'
	tt = get(h.slider_tps,'Value') ;
	set(h.edit_tps,'String',num2str(tt)) ;


% 'clique' souris sur les profils
case 'edit_prof'
%	disp('appel de zuieditcons')
	val = get(h.radio_zoom,'value') ;
	if val==1
		return
	end
	hplot = gca ;
	tag   = get(hplot,'tag') ;
	disp(tag)

	% variables d'entr� de zuieditcons en mode profil
	nom = tag ;
	switch tag
		case 'v1'
			nom  = ['f1 de ' nom_var] ;
			aide = ['dependance spatiale du mode 1 de ' nom_var] ;
			texte_y = 'f1' ;
		case 'v2'
			nom  = ['f2 de ' nom_var] ;
			aide = ['dependance spatiale du mode 2 de ' nom_var] ;
			texte_y = 'f2' ;
		case 'v3'
			nom  = ['f3 de ' nom_var] ;
			aide = ['dependance spatiale du mode 3 de ' nom_var] ;
			texte_y = 'f3' ;
	otherwise
			aide = ['???'] ;
	end
	x       = evalin('base','param.gene.x') ;
	y = evalin('base',['complexe.' nom_var '.'.' tag]) ;
	texte_x = 'x' ;
	var_x   = 'void';
	var_y   = ['complexe.' nom_var '.'.' tag] ;
	canal   = 1 ;
 	code_retour = 'complexe' ;
	liste_ref   = {'Pe','Pion','Ne','Psi','Jmoy','\alphae','Zeff', ...
	               'f_1(x)','f_2(x)','f_3(x)','vide'} ;
	var_ref   = {{'param.gene.x','data.prof.pe','-'},...
                     {'param.gene.x','data.prof.pion','-'},...
                     {'param.gene.x','data.prof.ne','-'},...
                     {'param.gene.x','data.prof.psi','-'},...
                     {'param.gene.x','data.prof.jmoy','-'},...
                     {'param.gene.x','data.prof.ae','-'},...
                     {'param.gene.x','data.prof.zeff','-'},...
		     {'param.gene.x',['complexe.' nom_var '.v1'],':'}, ...
	             {'param.gene.x',['complexe.' nom_var '.v2'],':'}, ...
	             {'param.gene.x',['complexe.' nom_var '.v3'],':'}, ...
	             {'[]','[]',''}};
	%texte_prop = 'Ip' ;
	%var_prop   = ones(size(y)) ;
	hout_edit  = zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                         code_retour,liste_ref,var_ref) ;
	return


% 'clique' souris sur les consignes
case 'edit_cons'
%	disp('appel de zuieditcons')
	val = get(h.radio_zoom,'value') ;
	if val==1
		return
	end
	hplot = gca ;
	tag   = get(hplot,'tag') ;
	disp(tag)

	% variables d'entree de l'editeur de consignes zuieditcons
	switch tag
		case 'u1'
			nom  = ['g1 de ' nom_var] ;
			aide = ['dependance temporelle du mode 1 de ' nom_var] ;
			texte_y = 'g1' ;
		case 'u2'
			nom  = ['g2 de ' nom_var] ;
			aide = ['dependance temporelle du mode 2 de ' nom_var] ;
			texte_y = 'g2' ;
		case 'u3'
			nom  = ['g3 de ' nom_var] ;
			aide = ['dependance temporelle du mode 3 de ' nom_var] ;
			texte_y = 'g3' ;
		otherwise
			aide = ['???'] ;
	end
	x       = evalin('base','data.gene.temps') ;
	y       = evalin('base',['complexe.' nom_var '.'.' tag]) ;
	texte_x = 'temps' ;
	var_x   = 'void';
	var_y   = ['complexe.' nom_var '.'.' tag] ;
	canal   = 1 ;
	code_retour = 'complexe' ;
	liste_ref = ['     Ip    |Vloop|Flux|FCI|FCE|LH|IDN|Gaz|Pompage|Zeffm|Glacon|', ...
	               'g_1(t)|g_2(t)|g_3(t)|vide'] ;
	var_ref   = {{'data.gene.temps','data.cons.ip',':'}, ...
	             {'data.gene.temps','data.cons.vloop',':'}, ...
	             {'data.gene.temps','data.cons.flux',':'}, ...
	             {'data.gene.temps','abs(data.cons.fci)',':'}, ...
	             {'data.gene.temps','abs(data.cons.fce)',':'}, ...
	             {'data.gene.temps','abs(data.cons.hyb)',':'}, ...
	             {'data.gene.temps','abs(data.cons.idn)',':'}, ...
	             {'data.gene.temps','data.cons.c',':'}, ...
	             {'data.gene.temps','data.cons.pomp',':'}, ...
	             {'data.gene.temps','data.cons.zeffm',':'}, ...
	             {'data.gene.temps','data.cons.glacon','o'}, ...
		     {'data.gene.temps',['complexe.' nom_var '.u1'],':'}, ...
	             {'data.gene.temps',['complexe.' nom_var '.u2'],':'}, ...
	             {'data.gene.temps',['complexe.' nom_var '.u3'],':'}, ...
	             {'[]','[]',''}};
	texte_prop = '' ;
	var_prop   = '' ;
	hout_edit = zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	return


case 'radio_zoom'
	val = get(h.radio_zoom,'value') ;
	if val==1
		zoom on
	else
		zoom off
	end
	return

case 'complexe'
%	disp('retour de zuieditcons')
	% trace de v1,v2,v3,u1,u2,u3
	linech(h.axes_plotv1,x,v1);
	linech(h.axes_plotv2,x,v2);
	linech(h.axes_plotv3,x,v3);
	linech(h.axes_plotu1,t,u1);
	linech(h.axes_plotu2,t,u2);
	linech(h.axes_plotu3,t,u3);
	set(h.axes_plotu1,'Xlim',[t(1) t(end)]) ;
	set(h.axes_plotu2,'Xlim',[t(1) t(end)]) ;
	set(h.axes_plotu3,'Xlim',[t(1) t(end)]) ;


% Boutons Annulation
case {'btn_annul','close'}
if ishandle(hfig)
   switch nom_var
   case 'data.prof.te'
   	evalin('base','zrecalcultemppres(''el'');');
   case 'data.prof.pe'
   	evalin('base','zrecalcultemppres(''el'');');
   case 'data.prof.ti'
   	evalin('base','zrecalcultemppres(''ion'');');
   case 'data.prof.pion'
   	evalin('base','zrecalcultemppres(''ion'');');
   end
	zuiformcache(hfig) ;
	zuireset(h.btn_annul) ;
	zuicloseone(hfig) ;
	return
end

% Boutons Validation
case {'btn_valid'}
	if ishandle(hfig)
		% reconstitution de la variable
		 plus1 = get(h.radio1_plus,'value') ;
		 plus2 = get(h.radio2_plus,'value') ;
		 y =         u1 * k1 * v1' + ...
    		     plus1 * u2 * k2 * v2' + ...
    		     plus2 * u3 * k3 * v3' ;
		 zassignin('base',nom_var,y) ;
		 zuisavenonok;
		switch nom_var
		case 'data.prof.te'
			evalin('base','zrecalcultemppres(''el'');');
		case 'data.prof.pe'
			evalin('base','zrecalcultemppres(''el'');');
		case 'data.prof.ti'
			evalin('base','zrecalcultemppres(''ion'');');
		case 'data.prof.pion'
			evalin('base','zrecalcultemppres(''ion'');');
		end
		zuireset(h.btn_valid)
		zuicloseone(hfig)
		return
	end
case 'radio1_plus'
	% rien a faire
case 'radio2_plus'
	% rien a faire
otherwise
	warning('action non prise en compte')
end
%drawnow

% recherche de l'indice
tt = str2num(get(h.edit_tps,'String')) ;
dt = abs(t - tt) ;
ind = min(find(dt==min(dt))) ;

% calcul de y
plus1 = get(h.radio1_plus,'value') ;
plus2 = get(h.radio2_plus,'value') ;
if plus1==1
	if plus2==1
		titre = 'u1*k1*v1+*u2*k2*v2+*u3*k3*v3' ;
	else
		titre = 'u1*k1*v1+*u2*k2*v2' ;
	end
else
	if plus2==1
		titre = 'u1*k1*v1+*u3*k3*v3' ;
	else
		titre = 'u1*k1*v1' ;
	end
end

y =         u1 * k1 * v1' + ...
    plus1 * u2 * k2 * v2' + ...
    plus2 * u3 * k3 * v3' ;

% image de y
ha = h.axes_img ;
cbmem = get(ha,'ButtonDownFcn');
tgmem = get(ha,'tag');
set(ha,'ActivePositionProperty','outerposition');
axes(ha) ;
imagesc(x,t,y)
xlabel('x')
ylabel('t')
set(gca,'ydir','normal')
colormap('default')
colorbar
title('y')
set(ha,'tag',tgmem,'ButtonDownFcn',cbmem);

% trace de y
linech(h.axes_ploty,x,yref(ind,:),'yref') ;
linech(h.axes_ploty,x,y(ind,:),'y') ;
hx=get(h.axes_ploty,'xlabel') ;
set(hx,'string','x');
ht=get(h.axes_ploty,'Title') ;
set(ht,'String',titre) ;


% memorisation des nouvelles valeurs
zassignin('base',['complexe.' nom_var '.v1'],v1) ;
zassignin('base',['complexe.' nom_var '.v2'],v2) ;
zassignin('base',['complexe.' nom_var '.v3'],v3) ;
zassignin('base',['complexe.' nom_var '.u1'],u1) ;
zassignin('base',['complexe.' nom_var '.u2'],u2) ;
zassignin('base',['complexe.' nom_var '.u3'],u3) ;
zassignin('base',['complexe.' nom_var '.k1'],k1) ;
zassignin('base',['complexe.' nom_var '.k2'],k2) ;
zassignin('base',['complexe.' nom_var '.k3'],k3) ;


function linech(h,x,y,t)

if nargin > 3
   hh = min(findobj(h,'type','line','tag',t));
else
   hh = min(findobj(h,'type','line'));
end
if ~isempty(hh)
   set(hh,'xdata',x,'ydata',y);
end

